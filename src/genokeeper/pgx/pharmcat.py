"""PharmCAT 3.4.0 on our own calls (track design, slice 1b).

PharmCAT's input is a GRCh38 VCF written from our gVCF calls (`genotype.gvcf.calls_at`) at
every position of `pharmcat_positions_3.4.0.vcf.bgz` (a release asset, `positions_path`),
with PharmCAT's own REF and ALT spellings, reference and no-call rows explicit, each
position's genotype by `genotype.sites.site_genotype`: a call that is unknown, graded low or
not graded, or beside an allele PharmCAT does not list, is a no-call, so PharmCAT reports it
missing and the gene call's reliability drops with it rather than resting a diplotype on a
weak call (the slice's ruling); on a male's chrX outside the PARs (G6PD) the genotype is one
allele. Slice 3f adds two inputs the calibration earns: UGT1A1*28's call (`UGT1A1_28`) used
at its low grade when the table trusts it, and read-backed phase sets (`GT:PS`) when the table
trusts phase; without either, the VCF is written as before.

PharmCAT's own `pharmcat` reads the VCF as written: a position the VCF lacks or calls `./.`
is a no-call. `--missing-to-ref` (`-0`) belongs to its VCF preprocessor, which is never run
(the VCF is already in PharmCAT's representation at its positions), and is never passed:
`pharmcat_args` builds the whole command. Its anchor comes from this module's argument list,
which neither the VCF nor an outside-call file can extend, because they are passed as paths.

`run_pharmcat` runs it in the pinned image (`containers.IMAGES["pharmcat"]`) with
`--network none`, and keeps the HTML report beside the JSON one for the user. `parse_report`
reads the JSON: per gene its diplotype, phenotype, activity score, call source and missing
positions, and per drug annotation of that gene its drug, source, classification and
implications. It never reads `drugRecommendation` or `dosingInformation`, and drops an
implication that speaks of a dose (`DOSING`): no dosing text leaves the report. Its anchor
comes from the keys this module reads, which a report's text cannot add to, because every
value is taken from a named key and the recommendation keys are never named.
"""

import gzip
import json
import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genokeeper.config import Grading, load_grading
from genokeeper.containers import run_container
from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.reliability import on_haploid_contig
from genokeeper.genotype.sites import site_genotype
from genokeeper.store.db import Level

PHARMCAT_VERSION = "3.4.0"
POSITIONS_NAME = f"pharmcat_positions_{PHARMCAT_VERSION}.vcf.bgz"
# The sample name written into PharmCAT's input.
SAMPLE = "SAMPLE"
# Where the container sees its inputs and writes its report.
CONTAINER_IN = "/pharmcat-in"
CONTAINER_OUTSIDE = "/pharmcat-outside"
CONTAINER_OUT = "/pharmcat-out"
# PharmCAT names its outputs `<base>.report.json`, `<base>.report.html` and so on.
BASE_NAME = "pharmcat"
# The genes PharmCAT 3.4.0 calls: those its positions file defines (INFO PX of the pinned
# file, read 2026-10-20) and those it takes only as outside calls (HLA-A, HLA-B, MT-RNR1).
# PharmGKB's annotations in these genes are PharmCAT's to answer where its report gives the
# gene a drug (`evidence.pharmgkb.left_by_pharmcat`).
PHARMCAT_GENES = frozenset(
    {
        *("ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6"),
        *("CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "F2", "F5", "G6PD", "IFNL3", "NAT2"),
        *("NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"),
        *("HLA-A", "HLA-B", "MT-RNR1"),
    }
)


def positions_path(lab: Path) -> Path:
    """Where the bootstrap files PharmCAT's positions (`step_pharmcat`)."""
    return lab / "pharmcat" / PHARMCAT_VERSION / POSITIONS_NAME


# UGT1A1*28 (rs3064744, the TA repeat): its three ALTs as PharmCAT 3.4.0's positions file
# spells them at chr2:233760233 (read from the pinned file, 2026-10-30). `pharmcat_positions`
# refuses a file that spells them otherwise. Its call is used at a low grade only when the
# calibration trusts it (`write_pharmcat_vcf`'s `trusted`).
UGT1A1_28: tuple[Allele, ...] = (
    ("chr2", 233760233, "CAT", "C"),
    ("chr2", 233760233, "CAT", "CATAT"),
    ("chr2", 233760233, "CAT", "CATATAT"),
)


def _rows(path: Path) -> Iterator[list[str]]:
    """The positions file's data rows, split into columns."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def pharmcat_positions(path: Path) -> list[Allele]:
    """Every `(chrom, pos, ref, alt)` PharmCAT defines, one per ALT, in the file's order.

    The file is bgzipped (gzip-compatible). A row without an ALT is refused: PharmCAT's
    positions always spell the alleles they define. So is a row at UGT1A1*28's position
    whose alleles are not `UGT1A1_28`'s."""
    alleles: list[Allele] = []
    for chrom, pos, _id, ref, alt in (row[:5] for row in _rows(path)):
        if alt in ("", "."):
            raise ValueError(f"{path}: {chrom}:{pos} has no ALT")
        row = [(chrom, int(pos), ref.upper(), a.upper()) for a in alt.split(",")]
        if (chrom, int(pos)) == UGT1A1_28[0][:2] and tuple(row) != UGT1A1_28:
            raise ValueError(
                f"{path}: rs3064744 at {chrom}:{pos} is {ref}>{alt}, not UGT1A1*28's"
                f" {UGT1A1_28[0][2]}>{','.join(a[3] for a in UGT1A1_28)}"
            )
        alleles.extend(row)
    return alleles


def gene_regions(path: Path, *, margin: int = 1000) -> list[tuple[str, str, int, int]]:
    """Each gene the positions file's `INFO/PX` names (the gene a position defines):
    `(gene, chrom, start, end)` from its first to its last position, ± `margin`, in the
    file's order of first appearance. A position without PX (a position of interest)
    defines no gene."""
    spans: dict[str, tuple[str, int, int]] = {}
    for row in _rows(path):
        chrom, pos = row[0], int(row[1])
        info = dict(field.split("=", 1) for field in row[7].split(";") if "=" in field)
        for gene in filter(None, info.get("PX", "").split(",")):
            seen = spans.get(gene)
            if seen is None:
                spans[gene] = (chrom, pos, pos)
            elif seen[0] != chrom:
                raise ValueError(f"{path}: {gene} on {seen[0]} and {chrom}")
            else:
                spans[gene] = (chrom, min(seen[1], pos), max(seen[2], pos))
    return [
        (gene, chrom, max(1, start - margin), end + margin)
        for gene, (chrom, start, end) in spans.items()
    ]


def write_pharmcat_vcf(
    calls: Mapping[Allele, Called],
    positions: Sequence[Allele],
    out: Path,
    *,
    male: bool,
    levels: Mapping[Allele, Level],
    grading: Grading | None = None,
    phases: Mapping[Allele, tuple[int, int]] | None = None,
    trusted: frozenset[Allele] = frozenset(),
) -> Path:
    """Write PharmCAT's input to `out`: one row per position of `positions` (as
    `pharmcat_positions` gives them), with the genotype `calls` give its ALTs, each call
    graded by `levels`. `grading` gives the PARs (default: the project's grading.toml).

    With `phases` (allele -> its phase set and the haplotype its ALT is on,
    `phasing.phase_sites`), every row is `GT:PS`: a heterozygous biallelic row whose allele
    has a phase set is written `0|1` or `1|0` with it, every other row with `.` for PS. A
    `trusted` allele's call is used at a low grade (the calibration measured it); it is
    still a no-call when its dosage is unknown or it is not graded.

    `out` holds the user's genotypes: the caller puts it under the run's work directory."""
    grading = grading if grading is not None else load_grading()
    sites: dict[tuple[str, int, str], list[str]] = {}
    for chrom, pos, ref, alt in positions:
        sites.setdefault((chrom, pos, ref), []).append(alt)
    contigs = list(dict.fromkeys(chrom for chrom, _, _ in sites))
    lines = [
        "##fileformat=VCFv4.2",
        # `genomelab`, the name before the rename: this text is in the step's key.
        f"##source=genomelab (our calls at PharmCAT {PHARMCAT_VERSION}'s positions)",
        *(f"##contig=<ID={chrom}>" for chrom in contigs),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        *(
            ['##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">']
            if phases is not None
            else []
        ),
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
    ]
    for (chrom, pos, ref), alts in sites.items():
        haploid = on_haploid_contig(chrom, pos, grading=grading, male=male)
        keys = [(chrom, pos, ref, a) for a in alts]
        # `site_genotype` writes a low call as a no-call; a trusted one counts as graded.
        graded = [
            Level.MEDIUM
            if k in trusted and levels.get(k) == Level.LOW
            else levels.get(k)
            for k in keys
        ]
        gt = site_genotype([calls.get(k) for k in keys], graded, haploid=haploid)
        fmt, value = "GT", gt
        if phases is not None:
            phase = phases.get(keys[0]) if len(keys) == 1 and gt == "0/1" else None
            fmt = "GT:PS"
            if phase is None:
                value = f"{gt}:."
            else:
                value = f"{'0|1' if phase[1] == 1 else '1|0'}:{phase[0]}"
        lines.append(
            f"{chrom}\t{pos}\t.\t{ref}\t{','.join(alts)}\t.\tPASS\t.\t{fmt}\t{value}"
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def pharmcat_args(vcf: Path, outside: Path | None, out_dir: Path) -> list[str]:
    """PharmCAT's command inside its container: `vcf` and `outside` (the `-po` file) by
    name under their mounts (`CONTAINER_IN`, `CONTAINER_OUTSIDE`), the JSON and HTML
    reports into `CONTAINER_OUT` (`out_dir` on the host) as `BASE_NAME.report.*`."""
    del out_dir  # mounted at CONTAINER_OUT
    args = [
        "pharmcat",
        "-vcf",
        f"{CONTAINER_IN}/{vcf.name}",
        "-o",
        CONTAINER_OUT,
        "-bf",
        BASE_NAME,
        "-reporterJson",
        "-reporterHtml",
    ]
    if outside is not None:
        args += ["-po", f"{CONTAINER_OUTSIDE}/{outside.name}"]
    return args


# --- the run and its report ----------------------------------------------------------------

# The core gate's label for PharmCAT's run.
LABEL = "dna:pgx"
REPORT_JSON = f"{BASE_NAME}.report.json"
# A dose, a dosing or dosage, or an amount: implications that speak of one are dropped.
DOSING = re.compile(
    r"\b(?:dos(?:e|es|ed|ing|age|ages)|mg|mcg|µg|units?/kg)\b", re.IGNORECASE
)
# The report's annotation sources, as the store names them.
_SOURCES = {
    "CPIC_GUIDELINE": "CPIC",
    "DPWG_GUIDELINE": "DPWG",
    "FDA_LABEL": "FDA label",
    "FDA_ASSOC": "FDA association",
}
# The allele PharmCAT names in a gene it has no call for: `Unknown/Unknown`, or `Unknown`
# for a one-allele gene (a male's G6PD, as its run showed).
_UNKNOWN = "Unknown"


class PharmcatFailed(RuntimeError):
    """PharmCAT ran but wrote no JSON report."""


def run_pharmcat(vcf: Path, outside: Path | None, out_dir: Path) -> Path:
    """Run PharmCAT offline on `vcf` with the outside calls `outside` (or none), behind the
    core gate as `LABEL`; its JSON report's path. The HTML report stays in `out_dir`, the
    run's report directory, for the user. Raises `PharmcatFailed` without a JSON report,
    and what `run_container` raises."""
    out_dir.mkdir(parents=True, exist_ok=True)
    mounts = [(vcf.parent, CONTAINER_IN, True)]
    if outside is not None:
        mounts.append((outside.parent, CONTAINER_OUTSIDE, True))
    mounts.append((out_dir, CONTAINER_OUT, False))
    run_container(
        "pharmcat", pharmcat_args(vcf, outside, out_dir), mounts=mounts, label=LABEL
    )
    report = out_dir / REPORT_JSON
    if not report.is_file():
        raise PharmcatFailed(f"PharmCAT wrote no {REPORT_JSON} in {out_dir}")
    return report


@dataclass(frozen=True)
class DrugRec:
    """A drug annotation of a gene's result: the drug, the source (`CPIC`, `DPWG`, `FDA
    label`, `FDA association`), its classification (CPIC's recommendation strength, e.g.
    `Strong`; `Unspecified` and `No recommendation` as PharmCAT writes them) and its
    implications, none of which speaks of a dose."""

    drug: str
    source: str
    classification: str | None
    implications: tuple[str, ...]


@dataclass(frozen=True)
class GeneResult:
    """One gene of PharmCAT's report: the diplotype (`or` between alternatives; None when
    PharmCAT has no call), the phenotype (a diplotype's phenotypes joined by `; `), the
    activity score as PharmCAT writes it (`2.0`, `≥3.0`), the call source (`MATCHER`,
    `OUTSIDE`, `NONE`), the positions PharmCAT found missing (`chr:pos`), and its drugs."""

    gene: str
    diplotype: str | None
    phenotype: str | None
    activity_score: str | None
    call_source: str
    missing_positions: tuple[str, ...]
    drugs: tuple[DrugRec, ...]


def _joined(parts: Sequence[str], separator: str) -> str | None:
    unique = list(dict.fromkeys(p for p in parts if p))
    return separator.join(unique) if unique else None


def _drugs(report: dict[str, Any]) -> dict[str, list[DrugRec]]:
    """Each gene's drug annotations, in the report's order; annotations of one drug, source
    and classification (populations) are one, their implications joined."""
    found: dict[str, dict[tuple[str, str, str | None], list[str]]] = {}
    for drugs in report.get("drugs", {}).values():
        for drug in drugs.values():
            source = str(_SOURCES.get(drug["source"], drug["source"]))
            for guideline in drug.get("guidelines", []):
                for annotation in guideline.get("annotations", []):
                    genes = {
                        d["gene"]
                        for genotype in annotation.get("genotypes", [])
                        for d in genotype.get("diplotypes", [])
                    }
                    key = (drug["name"], source, annotation.get("classification"))
                    kept = [
                        i
                        for i in annotation.get("implications", [])
                        if not DOSING.search(i)
                    ]
                    for gene in sorted(genes):
                        found.setdefault(gene, {}).setdefault(key, []).extend(kept)
    return {
        gene: [
            DrugRec(name, source, classification, tuple(dict.fromkeys(implications)))
            for (name, source, classification), implications in recs.items()
        ]
        for gene, recs in found.items()
    }


def parse_report(report: Path) -> list[GeneResult]:
    """The genes of PharmCAT's JSON report, in its order: `genes`, then the genes it calls
    without any drug annotation (`unannotatedGeneCalls`; F2 and F5 on 3.4.0), which have no
    drugs. A phenotype `n/a` is none."""
    data = json.loads(report.read_text(encoding="utf-8"))
    drugs = _drugs(data)
    results: list[GeneResult] = []
    unannotated = [(e["geneSymbol"], e) for e in data.get("unannotatedGeneCalls", [])]
    for gene, entry in [*data["genes"].items(), *unannotated]:
        diplotypes = [
            d
            for d in entry.get("sourceDiplotypes", [])
            if set(str(d.get("label")).split("/")) != {_UNKNOWN}
        ]
        chrom = entry.get("chr")
        results.append(
            GeneResult(
                gene=gene,
                diplotype=_joined([d["label"] for d in diplotypes], " or "),
                phenotype=_joined(
                    [
                        _joined(
                            [p for p in d.get("phenotypes", []) if p != "n/a"], "; "
                        )
                        or ""
                        for d in diplotypes
                    ],
                    " or ",
                ),
                activity_score=_joined(
                    [d.get("activityScore") or "" for d in diplotypes], " or "
                ),
                call_source=entry.get("callSource", "NONE"),
                missing_positions=tuple(
                    f"{v.get('chromosome', chrom)}:{v['position']}"
                    for v in entry.get("variants", [])
                    if not v.get("call")
                ),
                drugs=tuple(drugs.get(gene, [])),
            )
        )
    return results
