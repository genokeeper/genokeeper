"""The sample's sex and its sex-chromosome dosage: stage 2 measures them from the CRAM (and,
for the chrX allele balance, the run's VCF); the run takes `--sex` to override.

`measure_ratios`: the depth at MAPQ >= [sex] mapq over chrX's mappable non-PAR region
(`xnonpar`) and over chrY's mappable X-degenerate region (`xdr`), each divided by the depth
over the mappable autosomes (all three from `genokeeper.regions.SEX_REGIONS`): about 0.5 per
copy against two autosomal copies. `copies` turns a ratio into a copy count within
`[dosage] tolerance`, else None ("between").

The sex used for calling (`sex_from_ratios`) is the X copy count's: one X is "male" (chrX
called haploid), two "female", whatever the Y ratio, so a man with mosaic loss of Y keeps a
haploid X; an X ratio between one and two copies refuses and asks for `--sex`.

`x_het_share` reads the chrX heterozygous share from allele depths, not genotypes (a male
run calls chrX haploid, so it has no heterozygous genotype): among the SNV records inside
`xnonpar` (PASS or RefCall, DP >= [reliability] min_dp) whose ALT fraction is at least
`het_low`, the share whose fraction is at most `het_high`. `karyotype` names the pattern
the dosage suggests; every finding but the mosaic loss of Y is a screening result that
needs a karyotype.
"""

import shutil
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pysam

from genokeeper.config import DosageThresholds, Grading, grch38_fasta, load_paths
from genokeeper.qc import mosdepth
from genokeeper.regions import SEX_REGIONS, Interval, giab38_dir, merge
from genokeeper.store.db import Kind, Topic, Validity

type Sex = Literal["male", "female"]

SCREENING = "screening result, needs a karyotype"


class SexRefused(RuntimeError):
    """The X depth ratio is neither one copy nor two; the run needs `--sex`."""


@dataclass(frozen=True)
class Dosage:
    """The X and Y depth ratios over the autosomes, the chrX heterozygous share (None
    below `[dosage] min_x_sites`) and the supported-site count it was read from."""

    x_ratio: float
    y_ratio: float
    x_het_share: float | None
    x_sites: int


def measure_ratios(
    cram: Path,
    *,
    grading: Grading,
    fasta: Path | None = None,
    lab: Path | None = None,
    threads: int = 4,
) -> tuple[float, float]:
    """The `xnonpar` and `xdr` mean depths over the `autosomes` one, at `[sex] mapq`
    (`fasta` and the region BED default to the lab's)."""
    lab = lab if lab is not None else load_paths().lab
    fasta = fasta if fasta is not None else grch38_fasta(lab)
    regions = giab38_dir(lab) / SEX_REGIONS
    if not regions.is_file():
        raise FileNotFoundError(f"{regions}: run scripts/bootstrap.sh")
    work = Path(tempfile.mkdtemp(prefix=".sex-", dir=cram.parent))
    try:
        depth = mosdepth(
            cram, fasta, regions, work / "sex", mapq=grading.sex.mapq, threads=threads
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    autosomes = depth["autosomes"][0]
    if autosomes <= 0:
        raise SexRefused(f"{cram}: no autosomal depth at MAPQ >= {grading.sex.mapq}")
    if "xnonpar" not in depth:
        raise FileNotFoundError(
            f"{regions} has no xnonpar region: run scripts/bootstrap.sh"
        )
    x = depth["xnonpar"][0] / autosomes
    y = depth.get("xdr", (0.0, 0.0))[0] / autosomes
    return x, y


def copies(ratio: float, tolerance: float) -> int | None:
    """`n` when the ratio is within `tolerance` of n/2, else None ("between")."""
    n = round(ratio * 2)
    return n if abs(ratio - n / 2) <= tolerance else None


def sex_from_ratios(x_ratio: float, y_ratio: float, t: DosageThresholds) -> Sex:
    """One X copy is "male", two "female", whatever the Y ratio; any other X ratio
    refuses (`SexRefused`) naming both ratios and `--sex`."""
    x = copies(x_ratio, t.tolerance)
    if x == 1:
        return "male"
    if x == 2:
        return "female"
    pattern = "between one and two copies" if x is None else f"{x} copies"
    raise SexRefused(
        f"cannot tell the X ploidy for calling: chrX non-PAR depth is {x_ratio:.3f} of the"
        f" autosomal ({pattern}), chrY X-degenerate {y_ratio:.3f}; give --sex male or"
        " --sex female, and the sex-chromosome dosage statement will report the pattern"
    )


def infer_sex(
    cram: Path,
    *,
    grading: Grading,
    fasta: Path | None = None,
    lab: Path | None = None,
    threads: int = 4,
) -> Sex:
    """The sex `cram` shows: its X copy count (`fasta` and the region BED default to the
    lab's)."""
    x, y = measure_ratios(cram, grading=grading, fasta=fasta, lab=lab, threads=threads)
    return sex_from_ratios(x, y, grading.dosage)


def _inside(intervals: Sequence[Interval]) -> dict[str, list[tuple[int, int]]]:
    found: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in merge(intervals):
        found.setdefault(chrom, []).append((start, end))
    return found


def x_het_share(
    vcf: Path, *, xnonpar: Sequence[Interval], grading: Grading
) -> tuple[float | None, int]:
    """The share of supported chrX SNV records that are allele-balanced, and the
    supported-site count; the share is None below `[dosage] min_x_sites`. `vcf` must be
    indexed."""
    t = grading.dosage
    sites = balanced = 0
    with pysam.VariantFile(str(vcf)) as calls:
        sample = next(iter(calls.header.samples))
        for chrom, spans in _inside(xnonpar).items():
            for start, end in spans:
                for record in calls.fetch(chrom, start, end):
                    if not start <= record.start < end:
                        continue  # a record overlapping the span from before it
                    fraction = _alt_fraction(record, sample, grading.min_dp)
                    if fraction is None or fraction < t.het_low:
                        continue
                    sites += 1
                    balanced += fraction <= t.het_high
    share = balanced / sites if sites and sites >= t.min_x_sites else None
    return share, sites


def _alt_fraction(
    record: pysam.VariantRecord, sample: str, min_dp: int
) -> float | None:
    """The ALT allele fraction of a PASS or RefCall SNV at DP >= `min_dp`, else None."""
    alleles = record.alleles or ()
    if len(alleles) < 2 or any(len(a) != 1 for a in alleles):
        return None
    if list(record.filter.keys()) not in (["PASS"], ["RefCall"]):
        return None
    call = record.samples[sample]
    dp, ad = call.get("DP"), call.get("AD")
    if dp is None or dp < min_dp or not ad or any(n is None for n in ad):
        return None
    total = sum(ad)
    return sum(ad[1:]) / total if total else None


@dataclass(frozen=True)
class Karyotype:
    """The sex-chromosome pattern a dosage suggests, as its statement grades it."""

    label: str
    kind: Kind
    topic: Topic
    validity: Validity
    reason: str


def _count(x: int, y: int) -> str:
    return f"{44 + x + y},{'X' * x}{'Y' * y}"


def karyotype(d: Dosage, t: DosageThresholds) -> Karyotype:
    """The pattern, in order: one X with a Y ratio above `y_absent_max` and below the
    one-copy floor (however low) is possible mosaic loss of Y (research, limited); one X
    with no Y is 45,X; a ratio between copy counts is a possible mosaic pattern; one X and
    one Y with a chrX heterozygous share of `x_het_min` or more a possible XY/XXY mosaic;
    46,XY and 46,XX are non-findings; any other count is a finding. A pattern worded
    "possible" is a doubtful one, at validity limited (the design's "Where each kind of
    result goes"); a whole count is moderate. Every finding but the mosaic loss of Y is a
    screening result that needs a karyotype."""
    x, y = copies(d.x_ratio, t.tolerance), copies(d.y_ratio, t.tolerance)
    floor = 0.5 - t.tolerance
    risk, moderate = Topic.RISKS, Validity.MODERATE
    if x == 1 and t.y_absent_max < d.y_ratio < floor:
        lost = round(100 * (1 - d.y_ratio / 0.5))
        return Karyotype(
            "possible mosaic loss of Y",
            Kind.FINDING,
            Topic.RESEARCH,
            Validity.LIMITED,
            f"estimated share of cells without Y: {lost}% (information; not actionable)",
        )
    ratios = f"X ratio {d.x_ratio:.2f}, Y ratio {d.y_ratio:.2f}"
    if x == 1 and d.y_ratio <= t.y_absent_max:
        return Karyotype(
            "45,X", Kind.FINDING, risk, moderate, f"{SCREENING}: one X, no Y ({ratios})"
        )
    if x is None or y is None:
        return Karyotype(
            "possible mosaic sex-chromosome pattern",
            Kind.FINDING,
            risk,
            Validity.LIMITED,
            f"{SCREENING}: a ratio between copy counts ({ratios})",
        )
    het = d.x_het_share
    if (x, y) == (1, 1) and het is not None and het >= t.x_het_min:
        return Karyotype(
            "possible 46,XY/47,XXY mosaic",
            Kind.FINDING,
            risk,
            Validity.LIMITED,
            f"{SCREENING}: one X by depth, but {het:.0%} of {d.x_sites:,} chrX sites"
            " are allele-balanced",
        )
    if (x, y) == (1, 1):
        return Karyotype(
            "46,XY", Kind.NON_FINDING, risk, moderate, "XY, no sign of an extra X or Y"
        )
    if (x, y) == (2, 0):
        return Karyotype(
            "46,XX",
            Kind.NON_FINDING,
            risk,
            moderate,
            "XX, no sign of a missing or extra X",
        )
    return Karyotype(
        _count(x, y),
        Kind.FINDING,
        risk,
        moderate,
        f"{SCREENING}: another sex-chromosome count ({ratios})",
    )
