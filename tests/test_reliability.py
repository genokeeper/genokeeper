"""Genotype reliability's shared parts: GIAB strata, per-chromosome depth medians, the
depth range, and haploid loci; and `p_wrong`, the error rate measured on HG003 for the
call's cell. `grade_called` on the gVCF call layer's own calls is tested with that layer
(tests/test_gvcf.py)."""

import gzip
import subprocess
import sys
from pathlib import Path

import pysam
import pytest
from fixtures import gvcf_mini as gm

from genokeeper.calibrate.table import (
    GENOTYPE,
    VERSION_KEYS,
    CalibrationTable,
    CellKey,
    Counts,
)
from genokeeper.config import load_grading
from genokeeper.genotype import reliability
from genokeeper.genotype.gvcf import INDEL_BLOCK_WARNING, OTHER_ALLELE_WARNING, Called
from genokeeper.genotype.reliability import (
    AF_CLONAL,
    AF_LOW,
    OUTSIDE_CALIBRATED,
    ExtraCells,
    Reliability,
    StrataIndex,
    depth_medians,
    depth_usable,
    grade_called,
    guard_allele_fraction,
    on_haploid_contig,
)
from genokeeper.store.db import Level, RowType

GRADING = load_grading()
MEDIANS = {"chr1": 30.0, "chrX": 30.0}


def _write_bed(path: Path, rows: list[tuple[str, int, int]], *, bgzip: bool) -> Path:
    text = "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in rows)
    if not bgzip:
        path.write_text(text)
        return path
    plain = path.with_suffix("")
    plain.write_text(text)
    pysam.tabix_compress(str(plain), str(path), force=True)
    plain.unlink()
    return path


def test_strata_reads_bgzip_and_plain_bed_and_merges(tmp_path: Path) -> None:
    gz = _write_bed(
        tmp_path / "a.bed.gz",
        [("chr1", 99, 300), ("chr1", 150, 200), ("chr2", 0, 10)],
        bgzip=True,
    )
    plain = _write_bed(tmp_path / "b.bed", [("1", 249, 260)], bgzip=False)
    strata = StrataIndex.from_beds({"segdup": gz, "mhc": plain})
    # BED is 0-based half-open: [99, 300) covers 1-based 100..300.
    assert strata.flags_at("chr1", 99) == ()
    assert strata.flags_at("chr1", 100) == ("segdup",)
    # Past the nested [150, 200): only a merged index still finds the outer interval.
    assert strata.flags_at("chr1", 250) == ("mhc", "segdup")
    assert strata.flags_at("chr1", 300) == ("segdup",)
    assert strata.flags_at("chr1", 301) == ()
    assert strata.flags_at("chr2", 10) == ("segdup",)
    assert strata.flags_at("chr3", 5) == ()


def test_strata_reads_plain_gzip_and_skips_headers(tmp_path: Path) -> None:
    path = tmp_path / "c.bed.gz"
    with gzip.open(path, "wt") as f:
        f.write("track name=x\n#comment\nchr1\t9\t20\textra\n")
    strata = StrataIndex.from_beds({"kir": path})
    assert strata.flags_at("chr1", 10) == ("kir",)


def test_strata_refuses_a_malformed_line(tmp_path: Path) -> None:
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t20\t10\n")
    with pytest.raises(ValueError, match="bad.bed:1"):
        StrataIndex.from_beds({"vdj": bed})


@pytest.mark.tools
def test_depth_medians_from_pass_rows(tmp_path: Path) -> None:
    def at(pos: int, depth: int, filt: str = "PASS", gt: str = "0/1") -> str:
        half = depth // 2
        ad = [half, depth - half] if gt == "0/1" else [depth, 0]
        return gm.variant("chr1", pos, 1, ["T"], gt=gt, gq=50, ad=ad, filt=filt)

    records = [
        at(100, 10),
        at(200, 20),
        at(300, 30),
        at(400, 40),
        # RefCalls, NoCalls and reference blocks are not PASS rows: never in the median.
        at(500, 100, "RefCall", "0/0"),
        at(550, 110, "RefCall", "0/0"),
        at(600, 0, "NoCall", "./."),
        gm.block("chr1", 700, 800, min_dp=200),
    ]
    vcf = gm.write_gvcf(tmp_path / "dv.g.vcf.gz", records)
    assert depth_medians(vcf) == {"chr1": 25.0}


@pytest.mark.parametrize(
    ("chrom", "pos", "male", "haploid"),
    [
        ("chrX", 10_000_000, True, True),
        ("X", 10_000_000, True, True),
        ("chrY", 10_000_000, True, True),
        ("chrX", 100_000, True, False),  # GRCh38 PAR1
        ("chrX", 155_800_000, True, False),  # GRCh38 PAR2
        ("chrY", 56_900_000, True, False),  # GRCh38 PAR2
        ("chrX", 155_000_000, True, True),  # GRCh37's PAR2, not GRCh38's
        ("chrX", 10_000_000, False, False),
        ("chr1", 10_000_000, True, False),
    ],
)
def test_on_haploid_contig(chrom: str, pos: int, male: bool, haploid: bool) -> None:
    assert on_haploid_contig(chrom, pos, grading=GRADING, male=male) is haploid


@pytest.mark.parametrize(
    ("dp", "usable"), [(9, False), (10, True), (60, True), (61, False)]
)
def test_depth_usable_is_the_depth_warnings_range(dp: int, usable: bool) -> None:
    assert depth_usable("chr1", dp, grading=GRADING, medians=MEDIANS) is usable
    assert depth_usable("1", dp, grading=GRADING, medians=MEDIANS) is usable
    # Without a median only the lower bound holds.
    assert depth_usable("chr2", dp, grading=GRADING, medians=MEDIANS) is (dp >= 10)


# --- measured p_wrong (stage 2b) ---------------------------------------------------------

VERSIONS = {key: f"{key}-value" for key in VERSION_KEYS}


def _benchmark(tmp_path: Path) -> StrataIndex:
    """HG003's benchmark regions in these tests: chr1 from 1-based 1,001 to 9,000."""
    bed = tmp_path / "benchmark.bed"
    bed.write_text("chr1\t1000\t9000\n")
    return StrataIndex.from_beds({"benchmark": bed})


def _table(cells: dict[CellKey, tuple[int, int]]) -> CalibrationTable:
    """A table of genotype errors over units per cell, each cell a pool of its own."""
    counts = {
        (GENOTYPE, key): Counts(errors, total, 30.0 * total)
        for key, (errors, total) in cells.items()
    }
    return CalibrationTable.from_counts(counts, VERSIONS)


def _called(
    row_type: RowType, gq: int, *, gt: str = "0/1", warnings: tuple[str, ...] = ()
) -> Called:
    return Called(
        gt=gt,
        dosage="0" if gt == "0/0" else "1",
        row_type=row_type,
        gq=gq,
        dp=30,
        ad=None,
        vaf=None,
        warnings=warnings,
    )


def _grade(
    called: Called,
    pos: int,
    table: CalibrationTable,
    benchmark: StrataIndex,
    *,
    chrom: str = "chr1",
    strata: StrataIndex | None = None,
    indel: bool = False,
) -> Reliability:
    return grade_called(
        called,
        chrom,
        pos,
        strata=strata if strata is not None else StrataIndex.from_beds({}),
        medians=MEDIANS,
        grading=GRADING,
        indel=indel,
        haploid=False,
        table=table,
        benchmark=benchmark,
    )


def test_p_wrong_is_cell_rate(tmp_path: Path) -> None:
    """`p_wrong` is the measured rate of the call's cell (row type x allele class x region
    x GQ bin) with its interval; the display level's cut-offs apply to it, not to GQ."""
    benchmark = _benchmark(tmp_path)
    segdup = tmp_path / "segdup.bed"
    segdup.write_text("chr1\t4000\t5000\n")
    strata = StrataIndex.from_beds({"segdup": segdup})
    variant_cell = CellKey("variant", "snv", "easy", 40)
    refcall_cell = CellKey("RefCall", "snv", "easy", 30)
    block_cell = CellKey("ref_block", "snv", "easy", 50)
    difficult_cell = CellKey("variant", "snv", "difficult", 40)
    table = _table(
        {
            variant_cell: (50, 10_000),  # 5e-3: medium, where GQ 45 says 3e-5
            refcall_cell: (1, 1_000_000),  # 1e-6: high, where GQ 35 says 3e-4
            block_cell: (3, 1_000_000),
            difficult_cell: (200, 10_000),  # 2e-2: low
        }
    )
    cases = [
        (RowType.PASS, 45, "0/1", 2000, variant_cell, Level.MEDIUM),
        (RowType.REFCALL, 35, "0/0", 2000, refcall_cell, Level.HIGH),
        (RowType.REF_BLOCK, 50, "0/0", 2000, block_cell, Level.HIGH),
        (RowType.PASS, 45, "0/1", 4500, difficult_cell, Level.LOW),
    ]
    for row_type, gq, gt, pos, cell, level in cases:
        graded = _grade(
            _called(row_type, gq, gt=gt), pos, table, benchmark, strata=strata
        )
        rate = table.rate_for(cell)
        assert (graded.p_wrong, graded.p_wrong_ci) == (rate.rate, (rate.lo, rate.hi))
        assert graded.level == level, (row_type, pos)
        assert OUTSIDE_CALIBRATED not in graded.warnings
    # The caller's own estimate is gone, its warning with it.
    assert not hasattr(reliability, "CALLER_WARNING")


def test_outside_benchmark_warned(tmp_path: Path) -> None:
    """Outside HG003's benchmark (chrX and chrY among it: the truth covers chr1-22), or in
    a cell nothing was measured in, a call has no measured rate: `p_wrong` is None, it is
    warned, and its display level comes from its GQ as before."""
    benchmark = _benchmark(tmp_path)
    table = _table({CellKey("variant", "snv", "easy", 40): (50, 10_000)})
    for chrom, pos in (("chr1", 9500), ("chrX", 2000), ("chrY", 2000)):
        graded = _grade(_called(RowType.PASS, 45), pos, table, benchmark, chrom=chrom)
        assert (graded.p_wrong, graded.p_wrong_ci) == (None, None)
        assert OUTSIDE_CALIBRATED in graded.warnings
        assert graded.level == Level.HIGH  # GQ 45: 3e-5
    # Inside the benchmark, a cell the table lacks: a RefCall at GQ 25.
    unmeasured = _grade(_called(RowType.REFCALL, 25, gt="0/0"), 2000, table, benchmark)
    assert (unmeasured.p_wrong, unmeasured.p_wrong_ci) == (None, None)
    assert OUTSIDE_CALIBRATED in unmeasured.warnings
    assert unmeasured.level == Level.MEDIUM  # GQ 25: 3e-3
    # A no-call has no rate and is not "outside": its genotype is unknown.
    nocall = Called(
        gt="./.",
        dosage="unknown",
        row_type=RowType.NOCALL,
        gq=14,
        dp=16,
        ad=None,
        vaf=None,
        warnings=(),
    )
    graded = _grade(nocall, 2000, table, benchmark)
    assert (graded.p_wrong, graded.level) == (None, Level.LOW)
    assert OUTSIDE_CALIBRATED not in graded.warnings


def test_indel_caps_lifted_when_measured(tmp_path: Path) -> None:
    """An indel's reference inferred from a block is capped at medium only while
    unmeasured: its cell counts the indels blocks miss. One inferred from another allele
    called at its position stays capped: no cell measures an indel missed beside a call."""
    benchmark = _benchmark(tmp_path)
    table = _table(
        {
            CellKey("ref_block", "indel", "easy", 50): (1, 1_000_000),
            CellKey("variant", "indel", "easy", 40): (1, 1_000_000),
        }
    )
    in_block = _called(RowType.REF_BLOCK, 50, gt="0/0")
    measured = _grade(in_block, 2000, table, benchmark, indel=True)
    assert measured.level == Level.HIGH
    assert INDEL_BLOCK_WARNING in measured.warnings
    unmeasured = _grade(in_block, 9500, table, benchmark, indel=True)
    assert unmeasured.level == Level.MEDIUM
    assert {INDEL_BLOCK_WARNING, OUTSIDE_CALIBRATED} <= set(unmeasured.warnings)
    # The call layer's shape: another allele's record, the claim's dosage 0.
    other = Called(
        gt="0/1", dosage="0", row_type=RowType.PASS, gq=45, dp=30, ad=None, vaf=None,
        warnings=(OTHER_ALLELE_WARNING,),
    )  # fmt: skip
    for pos in (2000, 9500):
        graded = _grade(other, pos, table, benchmark, indel=True)
        assert graded.level == Level.MEDIUM
        assert OTHER_ALLELE_WARNING in graded.warnings


# --- the allele-fraction guard (slice 3c) ------------------------------------------------


def test_guard_clonal_cap_and_plain_warning() -> None:
    """A heterozygous call whose VAF is below the bound is warned; in a clonal gene its
    level is capped at low, elsewhere it keeps its level. The bound is strict, and a
    haploid het, a two-copy call or a call without a VAF is left alone."""
    high = Reliability(p_wrong=1e-4, level=Level.HIGH, flags=(), warnings=())

    def het(vaf: float | None, *, dosage: str = "1") -> Called:
        return Called(
            gt="0/1" if dosage == "1" else "1/1", dosage=dosage, row_type=RowType.PASS,
            gq=50, dp=40, ad=None, vaf=vaf, warnings=(),
        )  # fmt: skip

    def guard(called: Called, *, clonal: bool, haploid: bool = False) -> Reliability:
        return guard_allele_fraction(
            called, high, clonal=clonal, haploid=haploid, bound=0.30
        )

    clonal = guard(het(0.25), clonal=True)
    assert (clonal.level, clonal.warnings) == (Level.LOW, (AF_CLONAL,))
    assert clonal.p_wrong == high.p_wrong
    plain = guard(het(0.25), clonal=False)
    assert (plain.level, plain.warnings) == (Level.HIGH, (AF_LOW,))
    assert guard(het(0.30), clonal=True) == high  # strict <
    assert guard(het(0.25), clonal=True, haploid=True) == high
    assert guard(het(0.9, dosage="2"), clonal=True) == high
    assert guard(het(None), clonal=True) == high
    assert AF_CLONAL.startswith(AF_LOW)


# --- HG002's cells (slice 3c) -------------------------------------------------------------


def test_hg003_wins_and_haploid_cells_need_haploid(tmp_path: Path) -> None:
    """Inside HG003's benchmark HG003's cell rates the call even where an extra part's
    regions hold it; outside, an extra part rates only calls whose ploidy matches it."""
    benchmark = _benchmark(tmp_path)  # chr1 1,001-9,000
    cell = CellKey("variant", "snv", "easy", 40)
    hg003 = _table({cell: (50, 10_000)})
    regions = tmp_path / "hg002.bed"
    regions.write_text("chr1\t0\t20000\nchrX\t0\t5000\n")
    index = StrataIndex.from_beds({"hg002": regions})
    haploid_part = ExtraCells(
        table=_table({cell: (3, 10_000)}), regions=index, haploid=True
    )
    diploid_part = ExtraCells(
        table=_table({cell: (7, 10_000)}), regions=index, haploid=False
    )

    def grade(
        called: Called,
        chrom: str,
        pos: int,
        *,
        haploid: bool,
        extra: tuple[ExtraCells, ...],
    ) -> Reliability:
        return grade_called(
            called, chrom, pos, strata=StrataIndex.from_beds({}), medians=MEDIANS,
            grading=GRADING, indel=False, haploid=haploid, table=hg003,
            benchmark=benchmark, extra=extra,
        )  # fmt: skip

    both = (haploid_part, diploid_part)
    inside = grade(_called(RowType.PASS, 45), "chr1", 2000, haploid=False, extra=both)
    assert inside.p_wrong == 50 / 10_000
    outside = grade(_called(RowType.PASS, 45), "chr1", 9500, haploid=False, extra=both)
    assert outside.p_wrong == 7 / 10_000 and OUTSIDE_CALIBRATED not in outside.warnings
    hemizygous = _called(RowType.PASS, 45, gt="1")
    unrated = grade(hemizygous, "chrX", 2000, haploid=False, extra=(haploid_part,))
    assert unrated.p_wrong is None and OUTSIDE_CALIBRATED in unrated.warnings
    rated = grade(hemizygous, "chrX", 2000, haploid=True, extra=(haploid_part,))
    assert rated.p_wrong == 3 / 10_000
    assert rated.p_wrong_ci == (
        haploid_part.table.rate_for(cell).lo,
        haploid_part.table.rate_for(cell).hi,
    )
    assert rated.warnings == ()
    # Without extras, as before slice 3c.
    plain = grade(hemizygous, "chrX", 2000, haploid=True, extra=())
    assert OUTSIDE_CALIBRATED in plain.warnings


def test_reliability_imports_without_a_cycle() -> None:
    """`hg002_cells` imports `reliability` (for `ExtraCells`), never the reverse: either
    import order works in a fresh interpreter."""
    first, second = (
        "genokeeper.genotype.reliability",
        "genokeeper.calibrate.hg002_cells",
    )
    for a, b in ((first, second), (second, first)):
        done = subprocess.run(
            [sys.executable, "-c", f"import {a}; import {b}"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert done.returncode == 0, done.stderr
