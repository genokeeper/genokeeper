"""Genotype reliability: the chance a call is wrong, its display level, flags and warnings.

A call from our own GRCh38 gVCF (`grade_called`, reference blocks among them) carries
`p_wrong`, the genotype error rate measured on HG003 sequenced on the DNBSEQ-T7 for the
call's cell of the calibration table (`genokeeper.calibrate.table`), with its Jeffreys 95%
interval (`p_wrong_ci`). Outside HG003's benchmark regions (chrX and chrY among them: the
truth covers chr1-22), or in a cell nothing was measured in, there is no measured rate: the
call is warned `OUTSIDE_CALIBRATED` and its display level comes from DeepVariant's own
estimate, 10^(-GQ/10), as before stage 2b. Since slice 3c, HG002's cells measured where
HG003's truth does not reach (`ExtraCells`, `calibrate.hg002_cells`) rate such a call first,
by ploidy. Context flags are the GIAB stratifications the site falls in; every flag blocks
the high level. `guard_allele_fraction` warns a heterozygote at a low allele fraction and,
in the clonal genes, caps it at low.
"""

import gzip
import re
from bisect import bisect_right
from collections import Counter
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from genokeeper.calibrate.table import CalibrationTable, CellKey, CellRate, gq_bin
from genokeeper.config import Grading
from genokeeper.genotype.gvcf import INDEL_BLOCK_WARNING, OTHER_ALLELE_WARNING, Called
from genokeeper.store.db import Level, RowType
from genokeeper.tools import run_tool
from genokeeper.variants.key import ucsc_contig

OUTSIDE_CALIBRATED = "outside calibrated regions"
HAPLOID_HET = "het on a haploid contig"
DEPTH_WARNING = "depth outside expected range"
AF_LOW = "allele fraction low"
AF_CLONAL = (
    "allele fraction low: may be clonal, not inherited; confirm in another tissue"
)

HAPLOID_CONTIGS = frozenset({"chrX", "chrY"})  # in a male, outside the PARs
# The calibration table's row type for each gVCF row type that holds a genotype.
_CELL_ROW_TYPES = {
    RowType.PASS: "variant",
    RowType.REFCALL: "RefCall",
    RowType.REF_BLOCK: "ref_block",
}

_GT_SEP = re.compile(r"[/|]")


@dataclass(frozen=True)
class Reliability:
    """How far a call can be trusted: the number (measured on HG003, with its Jeffreys 95%
    interval), the level it shows at, flags and warnings."""

    p_wrong: float | None
    level: Level
    flags: tuple[str, ...]
    warnings: tuple[str, ...]
    p_wrong_ci: tuple[float, float] | None = None


# Merged, sorted 0-based half-open interval starts and ends of one contig.
type _Intervals = tuple[list[int], list[int]]


@dataclass(frozen=True)
class StrataIndex:
    """GIAB stratum membership: flag -> contig -> merged, sorted intervals."""

    intervals: dict[str, dict[str, _Intervals]]

    @classmethod
    def from_beds(cls, beds: dict[str, Path]) -> StrataIndex:
        """Read each flag's BED (plain, gzip or BGZF); contigs are normalised to `chr` names."""
        return cls({flag: _read_bed(path) for flag, path in beds.items()})

    def flags_at(self, chrom: str, pos: int) -> tuple[str, ...]:
        """The flags whose intervals contain the 1-based position `pos`, sorted by name."""
        contig, offset = ucsc_contig(chrom), pos - 1
        flags: list[str] = []
        for flag, contigs in self.intervals.items():
            starts, ends = contigs.get(contig, ([], []))
            i = bisect_right(starts, offset) - 1
            if i >= 0 and offset < ends[i]:
                flags.append(flag)
        return tuple(sorted(flags))


@dataclass(frozen=True)
class ExtraCells:
    """A part of HG002's cells (`calibrate.hg002_cells`), measured where HG003's truth does
    not reach: its table, the regions it was measured in, and whether it rates haploid
    calls (a male's chrX and chrY outside the PARs) or diploid ones. Defined here, where
    `grade_called` reads it, so `hg002_cells` imports this module and never the reverse."""

    table: CalibrationTable
    regions: StrataIndex
    haploid: bool


def _lines(path: Path) -> Iterator[str]:
    """The lines of a plain, gzip or BGZF text file, one at a time, without line ends."""
    with path.open("rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":  # gzip; BGZF is a series of gzip members
        with gzip.open(path, "rt") as f:
            yield from (line.rstrip("\r\n") for line in f)
    else:
        with path.open() as f:
            yield from (line.rstrip("\r\n") for line in f)


def _read_bed(path: Path) -> dict[str, _Intervals]:
    raw: dict[str, list[tuple[int, int]]] = {}
    for number, line in enumerate(_lines(path), start=1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        try:
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        except IndexError, ValueError:
            raise ValueError(f"{path}:{number}: not a BED line: {line!r}") from None
        if not 0 <= start <= end:
            raise ValueError(f"{path}:{number}: bad interval {start}-{end}")
        raw.setdefault(ucsc_contig(chrom), []).append((start, end))
    merged: dict[str, _Intervals] = {}
    for contig, intervals in raw.items():
        starts: list[int] = []
        ends: list[int] = []
        for start, end in sorted(intervals):
            if starts and start <= ends[-1]:
                ends[-1] = max(ends[-1], end)
            else:
                starts.append(start)
                ends.append(end)
        merged[contig] = (starts, ends)
    return merged


def depth_medians(vcf: Path) -> dict[str, float]:
    """The median DP of the PASS rows of `vcf`, per chromosome, keyed by its `chr` name."""
    out = run_tool(
        ["bcftools", "query", "-i", 'FILTER="PASS"', "-f", "%CHROM\t[%DP]\n", str(vcf)]
    )
    counts: dict[str, Counter[int]] = {}
    for line in out.splitlines():
        chrom, dp = line.split("\t")
        if dp != ".":
            counts.setdefault(ucsc_contig(chrom), Counter())[int(dp)] += 1
    return {chrom: _median(hist) for chrom, hist in counts.items()}


def _median(hist: Counter[int]) -> float:
    """The median of a histogram value -> count, without expanding it."""

    def at(rank: int) -> int:  # the value at 0-based `rank` in sorted order
        seen = 0
        for value in sorted(hist):
            seen += hist[value]
            if rank < seen:
                return value
        raise IndexError(rank)

    n = hist.total()
    return (at((n - 1) // 2) + at(n // 2)) / 2


def _indices(gt: str) -> list[str]:
    return _GT_SEP.split(gt)


def _in_par(chrom: str, pos: int, grading: Grading) -> bool:
    return any(start <= pos <= end for start, end in grading.par.get(chrom, []))


def depth_usable(
    chrom: str, dp: int, *, grading: Grading, medians: Mapping[str, float]
) -> bool:
    """Whether `dp` reads at a site on `chrom` are in the expected range: at least
    `min_dp`, and at most `max_dp_factor` times the chromosome's median where one is known.
    Outside it a call is warned (`DEPTH_WARNING`), and nothing is inferred from its reads."""
    median = medians.get(ucsc_contig(chrom))
    return dp >= grading.min_dp and (
        median is None or dp <= grading.max_dp_factor * median
    )


def on_haploid_contig(chrom: str, pos: int, *, grading: Grading, male: bool) -> bool:
    """Whether GRCh38 `chrom:pos` has one copy: a male's chrX or chrY outside the PARs."""
    contig = ucsc_contig(chrom)
    return male and contig in HAPLOID_CONTIGS and not _in_par(contig, pos, grading)


def _is_het(gt: str) -> bool:
    indices = _indices(gt)
    return "." not in indices and len(set(indices)) > 1


def grade_called(
    called: Called,
    chrom: str,
    pos: int,
    *,
    strata: StrataIndex,
    medians: Mapping[str, float],
    grading: Grading,
    indel: bool,
    haploid: bool,
    table: CalibrationTable,
    benchmark: StrataIndex,
    extra: Sequence[ExtraCells] = (),
) -> Reliability:
    """Grade a call from our own gVCF (`genotype.gvcf`) at GRCh38 `chrom:pos`, where its
    flags come from; `indel` says the claim's allele is one, `haploid` that the locus has
    one copy. `benchmark` holds HG003's benchmark regions, where `table` was measured;
    outside them, the first of `extra` (HG002's cells) whose ploidy is the locus's and
    whose regions hold the position rates the call instead.

    `p_wrong` is the genotype error `table` measured for the call's cell: its row type, the
    claim's allele class, its region (difficult where a flag applies) and its GQ bin, with
    the cell's interval as `p_wrong_ci`. Outside the benchmark, or in a cell nothing was
    measured in, `p_wrong` is None, the call is warned `OUTSIDE_CALIBRATED`, and its level
    is judged on DeepVariant's own 10^(-GQ/10) instead.

    The level, on that rate: low with no GQ, an unknown genotype, a rate above `p_low`, or
    a heterozygote on a haploid locus; high when the rate is at most `p_high`, no flag
    applies and nothing caps it; medium otherwise. An indel's reference inferred from a
    reference block (`INDEL_BLOCK_WARNING`) is capped at medium only while unmeasured: a
    block's indel cell counts the truth indels inside blocks. One inferred from another
    allele called at its position is always capped: the cells score called records, so
    none measures an indel missed beside another call."""
    flags = strata.flags_at(chrom, pos)
    warnings = list(called.warnings)
    low = called.dosage == "unknown"
    measured = (
        None
        if low
        else _measured(
            called, chrom, pos, flags, indel, table, benchmark, extra, haploid
        )
    )
    if measured is not None:
        p_wrong: float | None = measured.rate
        rate = measured.rate
    else:
        p_wrong = None
        rate = None if called.gq is None else 10 ** (-called.gq / 10)
        if not low and called.row_type in _CELL_ROW_TYPES:
            warnings.append(OUTSIDE_CALIBRATED)
    capped = False
    if indel and not low:
        if called.row_type == RowType.REF_BLOCK:
            capped = measured is None
            warnings.append(INDEL_BLOCK_WARNING)
        capped = capped or OTHER_ALLELE_WARNING in called.warnings
    if haploid and called.gt is not None and _is_het(called.gt):
        low = True
        warnings.append(HAPLOID_HET)
    if called.dp is not None and not depth_usable(
        chrom, called.dp, grading=grading, medians=medians
    ):
        warnings.append(DEPTH_WARNING)
    if rate is None or rate > grading.p_low or low:
        level = Level.LOW
    elif rate <= grading.p_high and not flags and not capped:
        level = Level.HIGH
    else:
        level = Level.MEDIUM
    return Reliability(
        p_wrong=p_wrong,
        level=level,
        flags=flags,
        warnings=tuple(warnings),
        p_wrong_ci=None if measured is None else (measured.lo, measured.hi),
    )


def guard_allele_fraction(
    called: Called, r: Reliability, *, clonal: bool, haploid: bool, bound: float
) -> Reliability:
    """The allele-fraction guard (slice 3c): a heterozygous call (dosage 1, not on a
    haploid locus) whose VAF is below `bound` is warned `AF_LOW`, or, where `clonal` (in
    one of the clonal genes' spans, `genotype.clonal`), `AF_CLONAL` with its level capped
    at low: a somatic clone in the blood shows as a heterozygote at a low fraction. `vaf`
    is DeepVariant's own `VAF`, AD[alt] / DP, the quantity `[af_guard] min_het_vaf` is
    measured on. Any other call is returned unchanged."""
    if called.dosage != "1" or haploid or called.vaf is None or called.vaf >= bound:
        return r
    if not clonal:
        return replace(r, warnings=(*r.warnings, AF_LOW))
    return replace(r, level=Level.LOW, warnings=(*r.warnings, AF_CLONAL))


def _measured(
    called: Called,
    chrom: str,
    pos: int,
    flags: tuple[str, ...],
    indel: bool,
    table: CalibrationTable,
    benchmark: StrataIndex,
    extra: Sequence[ExtraCells],
    haploid: bool,
) -> CellRate | None:
    """The measured genotype error of the call's cell, or None: a row with no genotype or
    GQ, a site neither in HG003's benchmark nor in an `extra` part of the locus's ploidy,
    or a cell nothing was measured in. HG003's table is tried first, then each part in
    turn."""
    row_type = _CELL_ROW_TYPES.get(called.row_type)
    if row_type is None or called.gq is None:
        return None
    key = CellKey(
        row_type,
        "indel" if indel else "snv",
        "difficult" if flags else "easy",
        gq_bin(called.gq),
    )
    tables = [table] if benchmark.flags_at(chrom, pos) else []
    tables += [
        part.table
        for part in extra
        if part.haploid == haploid and part.regions.flags_at(chrom, pos)
    ]
    for measured in tables:
        try:
            return measured.rate_for(key)
        except KeyError:
            continue
    return None
