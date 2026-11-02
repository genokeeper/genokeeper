"""Single-SNP research associations from the GWAS Catalog (track design, slice 1d,
"Research associations").

The download (`DOWNLOAD`, the catalog's associations with their ontology annotations,
fetched whole by the bootstrap's step_gwas into `<lab>/gwas/<its build date>/`) holds one
TSV row per reported association; its one member's date in the zip names the release. A
row is kept (`load_associations`) when:

- its p-value is at most 5e-8, genome-wide significance;
- it names one SNP by its rsID and the strongest risk allele is one base of that SNP
  (`rs1801282-C`): a haplotype, an interaction, an unknown risk allele (`?`) or an indel is
  no single-SNP claim;
- the catalog maps it to GRCh38 (`CHR_ID`, `CHR_POS`);
- its mapped trait, the one text of the row the store keeps, speaks of no dose
  (`pharmcat.DOSING`; `methadone dose measurement` and some 70 rows like it in the release
  of 2026-10-03). This control, against dosing text in the store: its anchor comes from
  that match on the stored text itself, which a catalog row cannot slip past because no
  other column of it is stored as text.

Of the rows with one rsID, risk allele and set of traits (their EFO IDs in any order; the
mapped trait's text for a row the catalog maps to none), the one with the smallest p-value
stands, its study and PubMed ID with it; the rest are counted as
`duplicate_association` (the controller's ruling R12). The effect is the catalog's
`OR or BETA` with `95% CI (TEXT)`: a CI text ending in `increase` or `decrease` marks a
beta, which the catalog gives unsigned with its direction for the risk allele, so a
decrease is stored negative; any other value is an odds ratio. There is no absolute risk.

The catalog gives no ALT, and its risk allele may sit on either strand. So each SNP is
located (`locate`) on the HGDP+1kGP panel's GRCh38 variants (`prs.pgsc.panel_pvar`: every
one, chrX included; ruling R22), whose single-base alleles give its REF, checked against
the FASTA, and its ALTs. The risk allele must be one of them on the forward strand; one
whose complement is also among them (a palindromic A/T or C/G SNP, or a site with both) is
skipped, as its strand cannot be told from the alleles. A risk allele that is an ALT is
located at it; one that is the REF at the panel's first ALT, carried as the ploidy less
that ALT's dosage (`genotype.sites.risk_copies`; ruling R26): one on a male's chrX outside
the PARs; beside another allele its copies are not countable, and `sources.gwas` drops it
(`other_allele_at_site`). Every row that is no statement is counted by why
(`GWAS_DROP_REASONS`).
This control, against a risk allele read on the wrong strand (its carrier graded a
non-carrier): its anchor comes from the panel's alleles, REF checked against the GRCh38
FASTA, which a catalog row cannot move because its risk allele is only compared with them,
never written into the allele the call is read at.
"""

import csv
import io
import re
import zipfile
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path

import pysam

from genokeeper.genotype.gvcf import Allele
from genokeeper.pgx.pharmcat import DOSING
from genokeeper.prs.pgsc import read_pvar
from genokeeper.variants.key import ucsc_contig

DOWNLOAD = "gwas-catalog-associations_ontology-annotated-full.zip"
# Genome-wide significance.
MAX_P = Decimal("5E-8")
# The columns read; no other is.
_SNPS, _RISK, _CHR, _POS = "SNPS", "STRONGEST SNP-RISK ALLELE", "CHR_ID", "CHR_POS"
_EFFECT, _CI, _P = "OR or BETA", "95% CI (TEXT)", "P-VALUE"
_PUBMED, _TRAIT, _URI, _STUDY, _GENE = (
    "PUBMEDID",
    "MAPPED_TRAIT",
    "MAPPED_TRAIT_URI",
    "STUDY ACCESSION",
    "MAPPED_GENE",
)
_COLUMNS = (
    *(_SNPS, _RISK, _CHR, _POS, _EFFECT, _CI, _P, _PUBMED, _TRAIT, _URI, _STUDY, _GENE),
)
# Where every row that is no statement went.
GWAS_DROP_REASONS = (
    "p_value_above_5e-8",
    "dosing_trait",
    "not_one_snp",
    "risk_allele_unknown",
    "no_grch38_position",
    "duplicate_association",
    "not_in_panel",
    "reference_differs",
    "strand_ambiguous",
    "risk_allele_mismatch",
    "other_allele_at_site",
)
_RSID = re.compile(r"rs\d+")
_RISK_ALLELE = re.compile(r"(rs\d+)-(\S+)")
_CHROMS = frozenset([*map(str, range(1, 23)), "X", "Y", "MT"])
_DIRECTION = re.compile(r"\b(increase|decrease)\s*$")
_NUMBER = r"-?\d+(?:\.\d+)?(?:[eE]-?\d+)?"
_BOUNDS = re.compile(rf"\[\s*({_NUMBER})\s*-\s*({_NUMBER})\s*\]")
_BASES = frozenset("ACGT")
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
# One gene symbol: the catalog joins two flanking genes with ` - ` and overlapping ones
# with `, `.
_GENE_SYMBOL = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


@dataclass(frozen=True)
class Association:
    """One kept association: its lead SNP (rsID, GRCh38 position), risk allele, mapped
    trait and its EFO IDs, the gene the catalog maps the SNP into (None when it maps it
    between or into several), effect (`OR` or `beta`, its value and CI; None where the
    catalog gives none), p-value as the catalog writes it, PubMed ID and study."""

    rsid: str
    chrom: str
    pos: int
    risk_allele: str
    trait: str
    gene: str | None
    efo_ids: tuple[str, ...]
    measure: str | None
    value: float | None
    ci: tuple[float, float] | None
    p_value: str
    pubmed_id: str
    study: str


@dataclass(frozen=True)
class Associations:
    """The associations kept, how many rows the download holds, where the others went,
    and the release (the member's build date, YYYYMMDD), the statements' source version."""

    associations: list[Association]
    input_rows: int
    dropped: dict[str, int]
    release: str


def _p_value(text: str) -> Decimal | None:
    try:
        return Decimal(text.strip())
    except InvalidOperation:
        return None


def _effect(
    value: str, ci_text: str
) -> tuple[str | None, float | None, tuple[float, float] | None]:
    """The measure, value and CI of a row's `OR or BETA` and `95% CI (TEXT)`."""
    try:
        number = float(value)
    except ValueError:
        return None, None, None
    bounds = _BOUNDS.search(ci_text)
    ci = None if bounds is None else (float(bounds[1]), float(bounds[2]))
    direction = _DIRECTION.search(ci_text)
    if direction is None:
        return ("OR", number, ci) if number > 0 else (None, None, None)
    sign = -1.0 if direction[1] == "decrease" else 1.0
    if ci is not None:
        low, high = sorted(sign * abs(b) for b in ci)
        ci = (low, high)
    return "beta", sign * abs(number), ci


def _reason(row: Mapping[str, str]) -> str | None:
    """Why a row is no single-SNP genome-wide association on GRCh38, or None."""
    p = _p_value(row[_P])
    if p is None or not p.is_finite() or p > MAX_P:
        return "p_value_above_5e-8"
    rsid = row[_SNPS].strip()
    risk = _RISK_ALLELE.fullmatch(row[_RISK].strip())
    if not _RSID.fullmatch(rsid) or risk is None or risk[1] != rsid:
        return "not_one_snp"
    if risk[2] not in ("A", "C", "G", "T"):
        return "risk_allele_unknown"
    if row[_CHR].strip() not in _CHROMS or not row[_POS].strip().isdigit():
        return "no_grch38_position"
    return "dosing_trait" if DOSING.search(row[_TRAIT]) else None


def load_associations(path: Path) -> Associations:
    """The download at `path` read as the module docstring says. A file that is not the
    download (not a zip of one TSV with the columns read) is refused (`ValueError`)."""
    dropped: dict[str, int] = dict.fromkeys(GWAS_DROP_REASONS, 0)
    kept: dict[tuple[str, str, tuple[str, ...]], tuple[Decimal, Association]] = {}
    rows = 0
    try:
        with zipfile.ZipFile(path) as z:
            (member,) = z.infolist()
            with io.TextIOWrapper(z.open(member), encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
                missing = [c for c in _COLUMNS if c not in (reader.fieldnames or [])]
                if missing:
                    raise ValueError(f"no column {', '.join(missing)}")
                for row in reader:
                    rows += 1
                    reason = _reason(row)
                    if reason is not None:
                        dropped[reason] += 1
                        continue
                    association = _association(row)
                    key = (
                        association.rsid,
                        association.risk_allele,
                        tuple(sorted(association.efo_ids)) or (association.trait,),
                    )
                    p = Decimal(association.p_value)
                    if key in kept:
                        dropped["duplicate_association"] += 1
                        if p >= kept[key][0]:
                            continue
                    kept[key] = (p, association)
    except (zipfile.BadZipFile, UnicodeDecodeError, ValueError, csv.Error) as e:
        raise ValueError(
            f"{path}: not the GWAS Catalog associations download ({e})"
        ) from e
    year, month, day = member.date_time[:3]
    return Associations(
        [a for _, a in kept.values()], rows, dropped, f"{year:04d}{month:02d}{day:02d}"
    )


def _association(row: Mapping[str, str]) -> Association:
    risk = _RISK_ALLELE.fullmatch(row[_RISK].strip())
    assert risk is not None  # `_reason` passed it
    efo_ids = tuple(
        dict.fromkeys(
            uri.strip().rsplit("/", 1)[-1]
            for uri in row[_URI].split(",")
            if uri.strip()
        )
    )
    measure, value, ci = _effect(row[_EFFECT].strip(), row[_CI])
    gene = row[_GENE].strip()
    return Association(
        rsid=risk[1],
        chrom=ucsc_contig(row[_CHR].strip()),
        pos=int(row[_POS]),
        risk_allele=risk[2],
        trait=row[_TRAIT].strip(),
        gene=gene if _GENE_SYMBOL.fullmatch(gene) else None,
        efo_ids=efo_ids,
        measure=measure,
        value=value,
        ci=ci,
        p_value=row[_P].strip(),
        pubmed_id=row[_PUBMED].strip(),
        study=row[_STUDY].strip(),
    )


def panel_alleles(
    pvar: Path, positions: Collection[tuple[str, int]]
) -> dict[tuple[str, int], list[tuple[str, str]]]:
    """The panel's single-base (REF, ALT) pairs at each of `positions` (`(chrom, pos)`, the
    FASTA's contig names) in the plink2 variant file `pvar`, one per ALT; a position it
    lacks, or has only indels at, has none."""
    wanted = set(positions)
    found: dict[tuple[str, int], list[tuple[str, str]]] = {}
    for chrom, pos, ref, alt in read_pvar(pvar):
        pair = (ref.upper(), alt.upper())
        if (chrom, pos) in wanted and all(b in _BASES for b in pair):
            found.setdefault((chrom, pos), []).append(pair)
    return found


def locate(
    associations: Sequence[Association],
    panel: Mapping[tuple[str, int], Sequence[tuple[str, str]]],
    *,
    fasta: Path,
) -> tuple[list[tuple[Association, Allele]], dict[str, int]]:
    """Each association at its SNP's panel allele (`panel_alleles`), REF checked against
    `fasta`, and where the others went (see the module docstring)."""
    located: list[tuple[Association, Allele]] = []
    dropped: dict[str, int] = dict.fromkeys(GWAS_DROP_REASONS, 0)
    with pysam.FastaFile(str(fasta)) as reference:
        contigs = set(reference.references)
        for a in associations:
            pairs = panel.get((a.chrom, a.pos), ())
            refs = {ref for ref, _ in pairs}
            site = refs | {alt for _, alt in pairs}
            base = (
                reference.fetch(a.chrom, a.pos - 1, a.pos).upper()
                if a.chrom in contigs
                else ""
            )
            if not pairs:
                reason = "not_in_panel"
            elif refs != {base}:
                reason = "reference_differs"
            elif a.risk_allele not in site:
                reason = "risk_allele_mismatch"
            elif _COMPLEMENT[a.risk_allele] in site:
                reason = "strand_ambiguous"
            else:
                alt = a.risk_allele if a.risk_allele != base else pairs[0][1]
                located.append((a, (a.chrom, a.pos, base, alt)))
                continue
            dropped[reason] += 1
    return located, dropped
