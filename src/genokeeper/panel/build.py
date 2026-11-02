"""The genome-wide claims of a ClinVar GRCh38 release (`build_claims`): every record in
scope, its allele normalised against the lab's GRCh38 FASTA, and every other record counted
by the reason it is not a claim.

Contigs follow the FASTA (`chr1`, `chrX`, `chrM`), not ClinVar (`1`, `MT`). Gene, CLNSIG,
review status, rsID, condition and its IDs come from the record itself.
"""

import re
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from genokeeper.annotate.vep import fasta_contigs
from genokeeper.evidence.clinvar import clnsig_class
from genokeeper.genelists import LISTED, GeneList
from genokeeper.tools import run_tool
from genokeeper.variants.key import make_vid, ucsc_contig

RISK_TERMS = frozenset(
    {
        "Established_risk_allele",
        "Likely_risk_allele",
        "Uncertain_risk_allele",
        "risk_factor",
    }
)
# The primary assembly (chr1..chr22, chrX, chrY, chrM in the FASTA); ClinVar's patch and alt
# contigs (e.g. NT_187361.1) are not claims.
PRIMARY_CONTIGS = frozenset([*map(str, range(1, 23)), "X", "Y", "MT"])

# A CLNSIG value joins terms with "/" (one classification spanning terms) and "|" (several);
# a comma belongs to a term ("Pathogenic,_low_penetrance") and never separates.
_TERM_SEP = re.compile(r"[/|]")


@dataclass(frozen=True)
class _Record:
    variation_id: int
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    info: dict[str, str]


def _text(value: Any) -> str:
    """An INFO value as ClinVar wrote it: pysam splits Number=. values on commas, so rejoin."""
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(str(part) for part in value)
    return str(value)


def _read(rec: pysam.VariantRecord) -> _Record:
    if rec.id is None or rec.ref is None:
        raise ValueError(f"ClinVar record without ID or REF at {rec.chrom}:{rec.pos}")
    info = {
        key: _text(rec.info.get(key))
        for key in ("CLNSIG", "CLNREVSTAT", "GENEINFO", "RS", "CLNDN", "CLNDISDB", "MC")
    }
    return _Record(
        variation_id=int(rec.id),
        chrom=rec.chrom,
        pos=rec.pos,
        ref=rec.ref,
        alts=tuple(rec.alts or ()),
        info=info,
    )


def _rsid(rs: str) -> str | None:
    return f"rs{rs.split(',')[0]}" if rs else None


def _symbols(geneinfo: str) -> list[str]:
    return [part.split(":")[0] for part in geneinfo.split("|") if part]


def _terms(clnsig: str) -> set[str]:
    return {term.strip("_ ") for term in _TERM_SEP.split(clnsig) if term}


def _gene(geneinfo: str, genes: set[str]) -> str:
    symbols = _symbols(geneinfo)
    listed = [symbol for symbol in symbols if symbol in genes]
    return (listed or symbols or [""])[0]


def is_primary(chrom: str) -> bool:
    """True for 1-22, X, Y and MT, in ClinVar (`1`, `MT`) or UCSC (`chr1`, `chrM`) style."""
    name = chrom.removeprefix("chr")
    return ("MT" if name == "M" else name) in PRIMARY_CONTIGS


# Where every record of the GRCh38 release that is no claim went.
CLAIM_DROP_REASONS = (
    "benign",
    "vus_outside_lists",
    "not_in_scope",
    "non_primary_contig",
    "contig_not_in_reference",
    "multiallelic",
    "no_alt",
    "symbolic_allele",
    "ref_disagrees_with_reference",
    "structural_variant",
)
# The usual cut between an indel and a structural variant. A ClinVar deletion or insertion
# written with explicit bases and longer than this is dropped as a symbolic one is: SNV and
# small-indel calling cannot genotype it, so the gVCF's reference blocks over it would read
# "not a carrier" and a no-call inside it "unknown". `carrier.NOT_COVERED` names
# the class for every gene. The cut is a convention, not a measured limit of DeepVariant.
STRUCTURAL_VARIANT_BASES = 50
# The classes evaluated in every gene (track design, 2a "Claims everywhere"): slice 1a's
# pathogenic, conflicting and risk classes, an unmapped CLNSIG (shown with its alert), and
# drug_response. A VUS is a claim only in a listed gene, where it is a watch item.
_EVERY_GENE = frozenset({"pathogenic", "conflicting", "risk", "drug_response"})
_BASES = re.compile(r"[ACGTN]+")


@dataclass(frozen=True)
class Claim:
    """One ClinVar GRCh38 record in scope: its claim, and its allele normalised against the
    lab's GRCh38 FASTA (`chrom` in `chr` names), which `vid` keys."""

    variation_id: int
    gene: str
    clnsig: str
    clnrevstat: str
    disease_ids: str  # CLNDISDB as ClinVar writes it
    condition: str
    chrom: str
    pos: int
    ref: str
    alt: str
    vid: str
    rsid: str | None
    # MC as ClinVar writes it (`SO:0001587|nonsense,...`): slice 1e's truncations.
    consequence: str = ""


@dataclass(frozen=True)
class ClaimSet:
    """The claims of one ClinVar GRCh38 release, its records counted, and where every
    record that is no claim went (`CLAIM_DROP_REASONS`)."""

    claims: list[Claim]
    input_rows: int
    dropped: dict[str, int]


def _claim_scope(record: _Record, genes: set[str]) -> str | None:
    """None for a record in scope, else the reason it is counted and not kept."""
    clnsig = record.info["CLNSIG"]
    if not clnsig:
        return "not_in_scope"  # no germline classification: no claim
    if _terms(clnsig) & RISK_TERMS:
        return None
    cls = clnsig_class(clnsig)
    if cls is None or cls in _EVERY_GENE:
        return None
    if cls == "vus":
        listed = any(s in genes for s in _symbols(record.info["GENEINFO"]))
        return None if listed else "vus_outside_lists"
    return "benign" if cls == "benign" else "not_in_scope"


def normalise_records(
    alleles: list[tuple[str, int, str, str]], *, fasta: Path
) -> dict[int, tuple[str, int, str, str]]:
    """Each `(chrom, pos, ref, alt)` (chrom in the FASTA's names) left-aligned and trimmed by
    `bcftools norm -f fasta -c x -m -both`, keyed by its index in `alleles`; an allele whose
    REF the FASTA does not have is left out.

    `-c x` drops such a record where `-c e` would stop the whole release on it (a ClinVar
    REF on an N, or one record of a release ClinVar wrote against another patch). The drop
    cannot silently empty the claims: its anchor comes from the md5 the bootstrap checks the
    FASTA against (NCBI's md5checksums.txt), which a wrong or truncated reference cannot
    move, because the bootstrap refuses a FASTA whose md5 differs; and every record left
    out is counted (`ref_disagrees_with_reference`)."""
    contigs = fasta_contigs(fasta)
    order = {name: i for i, name in enumerate(contigs)}
    lines = sorted(
        (order[chrom], pos, f"{chrom}\t{pos}\t{i}\t{ref}\t{alt}\t.\t.\t.\n")
        for i, (chrom, pos, ref, alt) in enumerate(alleles)
    )
    header = "##fileformat=VCFv4.2\n" + "".join(
        f"##contig=<ID={name},length={length}>\n" for name, length in contigs.items()
    )
    found: dict[int, tuple[str, int, str, str]] = {}
    with tempfile.TemporaryDirectory(prefix="genokeeper-claims-") as work:
        raw = Path(work) / "claims.vcf"
        raw.write_text(
            header
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join(line for _, _, line in lines)
        )
        out = Path(work) / "normalised.vcf"
        run_tool(
            [
                "bcftools", "norm", "--no-version", "-f", str(fasta), "-c", "x",
                "-m", "-both", "-Ov", "-o", str(out), str(raw),
            ]
        )  # fmt: skip
        with out.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                chrom, pos, number, ref, alt = line.split("\t")[:5]
                found[int(number)] = (chrom, int(pos), ref.upper(), alt.upper())
    return found


def build_claims(
    clinvar38: Path, lists: Mapping[str, GeneList], *, fasta: Path
) -> ClaimSet:
    """The claims of the ClinVar GRCh38 release `clinvar38` in every gene, normalised
    against `fasta` (the lab's GRCh38 no-alt FASTA); every other record counted by reason.

    Gene, CLNSIG, review status, rsID, condition and its IDs come from the record itself:
    the GRCh38 release alone, no GRCh37 join. The gene is the first listed symbol of
    GENEINFO, else its first symbol."""
    genes = {gene for name in LISTED for gene in lists[name].genes}
    contigs = fasta_contigs(fasta)
    dropped: dict[str, int] = dict.fromkeys(CLAIM_DROP_REASONS, 0)
    candidates: list[tuple[_Record, str]] = []
    input_rows = 0
    with pysam.VariantFile(str(clinvar38)) as vcf:
        for rec in vcf:
            input_rows += 1
            record = _read(rec)
            reason = _claim_scope(record, genes)
            chrom = ucsc_contig(record.chrom)
            alt = record.alts[0] if len(record.alts) == 1 else "."
            if reason is None and not is_primary(record.chrom):
                reason = "non_primary_contig"
            if reason is None and len(record.alts) > 1:
                reason = "multiallelic"
            if reason is None and alt == ".":
                reason = "no_alt"
            if reason is None and not (
                _BASES.fullmatch(record.ref.upper()) and _BASES.fullmatch(alt.upper())
            ):
                reason = "symbolic_allele"
            if reason is None and chrom not in contigs:
                reason = "contig_not_in_reference"
            if reason is not None:
                dropped[reason] += 1
            else:
                candidates.append((record, chrom))
    normalised = normalise_records(
        [(chrom, r.pos, r.ref.upper(), r.alts[0].upper()) for r, chrom in candidates],
        fasta=fasta,
    )
    claims: list[Claim] = []
    for i, (record, _) in enumerate(candidates):
        allele = normalised.get(i)
        if allele is None:
            dropped["ref_disagrees_with_reference"] += 1
            continue
        chrom, pos, ref, alt = allele
        if max(len(ref), len(alt)) - 1 > STRUCTURAL_VARIANT_BASES:
            dropped["structural_variant"] += 1
            continue
        claims.append(
            Claim(
                variation_id=record.variation_id,
                gene=_gene(record.info["GENEINFO"], genes),
                clnsig=record.info["CLNSIG"],
                clnrevstat=record.info["CLNREVSTAT"],
                disease_ids=record.info["CLNDISDB"],
                condition=record.info["CLNDN"],
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                vid=make_vid(chrom, pos, ref, alt),
                rsid=_rsid(record.info["RS"]),
                consequence=record.info["MC"],
            )
        )
    order = {name: i for i, name in enumerate(contigs)}
    claims.sort(key=lambda c: (order[c.chrom], c.pos, c.ref, c.alt, c.variation_id))
    return ClaimSet(claims=claims, input_rows=input_rows, dropped=dropped)


def check_clinvar38(clinvar38: Path) -> None:
    """Raise unless `clinvar38` is a ClinVar GRCh38 release the run's reader takes: every
    record reads as `build_claims` reads it (a numeric ID, a REF), and at least one carries
    a classification that is a claim in every gene. The bootstrap runs it on a download
    before filing it, so a release this reader would stop on never becomes the newest."""
    in_scope = 0
    with pysam.VariantFile(str(clinvar38)) as vcf:
        for rec in vcf:
            in_scope += _claim_scope(_read(rec), set()) is None
    if not in_scope:
        raise ValueError(
            f"{clinvar38}: no record with a germline classification in scope"
        )


def clinvar_vids(clinvar38: Path, *, fasta: Path) -> set[str]:
    """The vids of the alleles of every record of the release ClinVar classifies germline
    (any class but `not_in_scope`, whatever the record's shape): each ALT, in explicit
    bases, on a primary contig the FASTA holds, normalised against `fasta`. Slice 1e's VEP
    does not annotate them again; an allele of a record without a germline classification,
    or of no record, it does."""
    contigs = fasta_contigs(fasta)
    alleles: list[tuple[str, int, str, str]] = []
    with pysam.VariantFile(str(clinvar38)) as vcf:
        for rec in vcf:
            chrom = ucsc_contig(rec.chrom)
            if rec.ref is None or not is_primary(rec.chrom) or chrom not in contigs:
                continue
            # With no listed gene a VUS scopes `vus_outside_lists`: classified either way.
            if _claim_scope(_read(rec), set()) == "not_in_scope":
                continue
            ref = rec.ref.upper()
            for alt in rec.alts or ():
                if _BASES.fullmatch(ref) and _BASES.fullmatch(alt.upper()):
                    alleles.append((chrom, rec.pos, ref, alt.upper()))
    return {make_vid(*a) for a in normalise_records(alleles, fasta=fasta).values()}
