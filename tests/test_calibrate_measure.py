"""The calibration's measurement: HG003's gVCF scored against GIAB's truth inside the
benchmark regions, error per GQ bin, row type, allele class and region, and the table that
keeps it with the versions of the pipeline it measured.

The scoring tests run the real `rtg vcfeval` (installed by the bootstrap) on the mini
reference of `gvcf_mini`; nothing here is committed as a genomic file."""

import itertools
import json
import math
import random
import stat
import subprocess
from collections.abc import Sequence
from pathlib import Path

import pysam
import pytest
from fixtures.gvcf_mini import (
    CONTIGS,
    block,
    other_base,
    ref,
    variant,
    write_fasta,
    write_gvcf,
)

from genokeeper import tools
from genokeeper.calibrate import rtg_path
from genokeeper.calibrate.build import MAX_OUTSIDE, SPLITS, held_out, merged, report
from genokeeper.calibrate.measure import Cells, measure
from genokeeper.calibrate.table import (
    ALLELE,
    GENOTYPE,
    GQ_BINS,
    VERSION_KEYS,
    CalibrationTable,
    CellKey,
    Counts,
    beta_binomial_interval,
    betainc,
    dispersion,
    jeffreys,
)
from genokeeper.config import load_paths
from genokeeper.genotype.reliability import StrataIndex

VERSIONS = {key: f"{key}-value" for key in VERSION_KEYS}
TRUTH_SAMPLE = "HG003"


@pytest.fixture(scope="module")
def reference(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    """The mini FASTA and its rtg SDF."""
    folder = tmp_path_factory.mktemp("reference")
    fasta = write_fasta(folder / "mini.fa")
    sdf = folder / "mini.sdf"
    subprocess.run(
        [str(rtg_path(load_paths().lab)), "format", "-o", str(sdf), str(fasta)],
        check=True,
        capture_output=True,
    )
    return fasta, sdf


def _truth(path: Path, records: Sequence[tuple[str, int, str, str, str]]) -> Path:
    """GIAB-shaped truth: (contig, pos, ref, alt, GT) rows, bgzipped and indexed."""
    header = (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CONTIGS.items())
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{TRUTH_SAMPLE}\n"
    )
    rows = sorted(records, key=lambda r: (list(CONTIGS).index(r[0]), r[1]))
    plain = path.with_suffix("")
    plain.write_text(
        header
        + "".join(
            f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t{gt}\n" for c, p, r, a, gt in rows
        )
    )
    return Path(pysam.tabix_index(str(plain), preset="vcf", force=True))


def _bed(path: Path, intervals: Sequence[tuple[str, int, int]]) -> Path:
    path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals))
    return path


def _snv(contig: str, pos: int, gt: str, k: int = 0) -> tuple[str, int, str, str, str]:
    return (contig, pos, ref(contig, pos), other_base(ref(contig, pos), k), gt)


def _call(contig: str, pos: int, *, gt: str, gq: int, k: int = 0) -> str:
    return variant(
        contig, pos, 1, [other_base(ref(contig, pos), k)], gt=gt, gq=gq, ad=[10, 10]
    )


def _measure(
    tmp_path: Path,
    reference: tuple[Path, Path],
    *,
    gvcf: list[str],
    truth: Sequence[tuple[str, int, str, str, str]],
    bench: Sequence[tuple[str, int, int]],
    difficult: Sequence[tuple[str, int, int]] = (),
) -> CalibrationTable:
    _, sdf = reference
    strata = StrataIndex.from_beds(
        {"difficult": _bed(tmp_path / "difficult.bed", difficult)}
    )
    return measure(
        write_gvcf(tmp_path / "calls.g.vcf.gz", gvcf),
        _truth(tmp_path / "truth.vcf.gz", truth),
        _bed(tmp_path / "bench.bed", bench),
        sdf=sdf,
        strata=strata,
        versions=VERSIONS,
    )


def _cell(
    table: CalibrationTable,
    kind: str,
    row_type: str,
    allele_class: str,
    region: str,
    gq: int,
) -> tuple[int, int, float | None]:
    """A cell's own errors, units and mean depth."""
    found = table.rate_for(CellKey(row_type, allele_class, region, gq), errors=kind)
    return found.errors, found.total, found.depth


# One truth variant matched by one call: vcfeval has something to score in every test.
ANCHOR_TRUTH = _snv("chr1", 3500, "0/1")
ANCHOR_CALL = _call("chr1", 3500, gt="0/1", gq=45)
ANCHOR_BENCH = ("chr1", 3400, 3600)


@pytest.mark.tools
def test_block_denominator_is_positions_in_benchmark(
    tmp_path: Path, reference: tuple[Path, Path]
) -> None:
    table = _measure(
        tmp_path,
        reference,
        gvcf=[
            block("chr1", 1, 1000, gq=50, min_dp=20),
            block("chr1", 1001, 1200, gq=12, min_dp=8),
            block("chr1", 2201, 2300, gq=1, min_dp=0),  # no coverage: never a call
            block("chr1", 2301, 2400, gq=0, gt="./."),  # a no-call block
            ANCHOR_CALL,
        ],
        truth=[ANCHOR_TRUTH],
        bench=[
            ("chr1", 0, 300),
            ("chr1", 500, 600),
            ("chr1", 1100, 1150),
            ("chr1", 2150, 2400),
            ANCHOR_BENCH,
        ],
        difficult=[("chr1", 250, 280)],
    )
    for kind in (GENOTYPE, ALLELE):
        for allele_class in ("snv", "indel"):
            cell = (kind, "ref_block", allele_class)
            # 0-300 and 500-600 of the GQ 50 block, 30 of them difficult.
            assert _cell(table, *cell, "easy", 50) == (0, 270 + 100, 20.0)
            assert _cell(table, *cell, "difficult", 50) == (0, 30, 20.0)
            # 1100-1150 of the GQ 12 block.
            assert _cell(table, *cell, "easy", 10) == (0, 50, 8.0)
    assert sorted(
        key.gq_bin
        for key in table.measured(errors=GENOTYPE)
        if key.row_type == "ref_block"
    ) == [10, 10, 50, 50, 50, 50]


@pytest.mark.tools
def test_truth_variant_in_block_is_error(
    tmp_path: Path, reference: tuple[Path, Path]
) -> None:
    deletion = ("chr1", 600, ref("chr1", 600, 3), ref("chr1", 600), "0/1")
    table = _measure(
        tmp_path,
        reference,
        gvcf=[
            block("chr1", 1, 600, gq=20, min_dp=20),
            block("chr1", 601, 1000, gq=50, min_dp=20),
            ANCHOR_CALL,
        ],
        truth=[
            _snv("chr1", 400, "0/1"),
            # Its anchor, 600, ends the GQ 20 block; its first deleted base, 601, is the
            # GQ 50 block's error.
            deletion,
            _snv("chr1", 1500, "1/1"),  # outside every block: not a block's error
            ANCHOR_TRUTH,
        ],
        bench=[("chr1", 0, 2000), ANCHOR_BENCH],
    )
    for kind in (GENOTYPE, ALLELE):
        snv20, snv50, indel20, indel50 = (
            _cell(table, kind, "ref_block", allele_class, "easy", gq)[:2]
            for allele_class in ("snv", "indel")
            for gq in (20, 50)
        )
        assert (snv20, snv50) == ((1, 600), (0, 400))
        assert (indel20, indel50) == ((0, 600), (1, 400))
        snv = table.rate_for(CellKey("ref_block", "snv", "easy", 20), errors=kind)
        assert snv.rate == pytest.approx(1 / 600)
        # The indel rate rises with GQ here: the two bins pool.
        indel = table.rate_for(CellKey("ref_block", "indel", "easy", 50), errors=kind)
        assert indel.rate == pytest.approx(1 / 1000)


@pytest.mark.tools
def test_a_multi_allelic_truth_site_is_one_error(
    tmp_path: Path, reference: tuple[Path, Path]
) -> None:
    # Two truth alleles at one position: one unit, so one error, not two.
    base = ref("chr1", 400)
    two_alts = ("chr1", 400, base, f"{other_base(base)},{other_base(base, 1)}", "1|2")
    table = _measure(
        tmp_path,
        reference,
        gvcf=[block("chr1", 1, 1000, gq=20, min_dp=20), ANCHOR_CALL],
        truth=[two_alts, ANCHOR_TRUTH],
        bench=[("chr1", 0, 1000), ANCHOR_BENCH],
    )
    for kind in (GENOTYPE, ALLELE):
        assert _cell(table, kind, "ref_block", "snv", "easy", 20)[:2] == (1, 1000)


@pytest.mark.tools
def test_truth_variant_at_refcall_is_error(
    tmp_path: Path, reference: tuple[Path, Path]
) -> None:
    table = _measure(
        tmp_path,
        reference,
        gvcf=[
            _call("chr1", 700, gt="0/0", gq=25).replace("\tPASS\t", "\tRefCall\t"),
            _call("chr1", 800, gt="0/0", gq=27).replace("\tPASS\t", "\tRefCall\t"),
            _call("chr1", 900, gt="0/0", gq=5).replace("\tPASS\t", "\tRefCall\t"),
            ANCHOR_CALL,
        ],
        truth=[_snv("chr1", 700, "0/1"), ANCHOR_TRUTH],
        bench=[("chr1", 0, 2000), ANCHOR_BENCH],
    )
    for kind in (GENOTYPE, ALLELE):
        at = table.rate_for(CellKey("RefCall", "snv", "easy", 20), errors=kind)
        assert (at.errors, at.total) == (1, 2)
        low = table.rate_for(CellKey("RefCall", "snv", "easy", 0), errors=kind)
        assert (low.errors, low.total) == (0, 1)
    # A RefCall is a reference call: vcfeval skips it, so it is no variant row.
    assert all(
        key.row_type != "variant" or key.gq_bin == 40 for key in table.measured()
    )


@pytest.mark.tools
def test_genotype_and_allele_errors_apart(
    tmp_path: Path, reference: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    log = tmp_path / "gate.log"
    gate = tmp_path / "coreslot-gate"
    gate.write_text(f'#!/usr/bin/env bash\necho "$*" >> "{log}"\nshift 2\nexec "$@"\n')
    gate.chmod(gate.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setattr(tools, "GATE", gate)
    deletion = ("chr1", 1300, ref("chr1", 1300, 2), ref("chr1", 1300), "0/1")
    table = _measure(
        tmp_path,
        reference,
        gvcf=[
            _call("chr1", 1000, gt="0/1", gq=35),  # right
            _call("chr1", 1100, gt="0/1", gq=33),  # right allele, wrong zygosity
            _call("chr1", 1200, gt="0/1", gq=31),  # no variant in the truth
            variant("chr1", 1300, 2, [ref("chr1", 1300)], gt="0/1", gq=38, ad=[9, 9]),
            _call("chr1", 1400, gt="0/1", gq=31, k=1).replace(
                "\tPASS\t", "\tNoCall\t"
            ),  # not PASS: not scored
            ANCHOR_CALL,
        ],
        truth=[
            _snv("chr1", 1000, "0/1"),
            _snv("chr1", 1100, "1/1"),
            deletion,
            ANCHOR_TRUTH,
        ],
        bench=[("chr1", 900, 1500), ANCHOR_BENCH],
    )
    snv = CellKey("variant", "snv", "easy", 30)
    genotype = table.rate_for(snv, errors=GENOTYPE)
    allele = table.rate_for(snv, errors=ALLELE)
    assert (genotype.errors, genotype.total) == (2, 3)
    assert (allele.errors, allele.total) == (1, 3)
    assert genotype.depth == 20.0  # DP of the three records
    indel = table.rate_for(CellKey("variant", "indel", "easy", 30), errors=GENOTYPE)
    assert (indel.errors, indel.total) == (0, 1)
    assert table.rate_for(snv) == genotype  # p_wrong is the genotype error
    # F1 over every genotype-exact call: SNPs 2 right of 4 called, 2 of 3 in the truth.
    assert table.f1 == pytest.approx({"snp": 4 / 7, "indel": 1.0})
    # vcfeval twice behind the core gate: genotype-exact, then zygosity ignored.
    gated = log.read_text().splitlines()
    assert all(line.startswith("dna:calibrate -- ") for line in gated)
    runs = [line for line in gated if " vcfeval " in line]
    assert ["--squash-ploidy" in run for run in runs] == [False, True]
    assert all("--output-mode annotate" in run for run in runs)


def _counts(rows: dict[int, tuple[int, int]]) -> dict[tuple[str, CellKey], Counts]:
    return {
        (GENOTYPE, CellKey("variant", "snv", "easy", gq)): Counts(e, n, 30.0 * n)
        for gq, (e, n) in rows.items()
    }


def test_rate_monotone_in_gq() -> None:
    # GQ 20 measures worse than GQ 10: the two pool (pool-adjacent-violators).
    rows = {0: (50, 100), 10: (10, 100), 20: (20, 100), 30: (1, 100), 50: (0, 1000)}
    table = CalibrationTable.from_counts(_counts(rows), VERSIONS)
    rates = [
        table.rate_for(CellKey("variant", "snv", "easy", gq)).rate
        for gq in sorted(rows)
    ]
    assert rates == pytest.approx([0.5, 0.15, 0.15, 0.01, 0.0])
    assert all(a >= b for a, b in itertools.pairwise(rates))
    # The cells keep their own counts; the rate and its interval are the pooled run's.
    pooled = table.rate_for(CellKey("variant", "snv", "easy", 20))
    assert (pooled.errors, pooled.total) == (20, 100)
    assert (pooled.lo, pooled.hi) == pytest.approx(jeffreys(30, 200))
    assert pooled.depth == 30.0
    # A cell nobody measured has no rate.
    with pytest.raises(KeyError):
        table.rate_for(CellKey("variant", "snv", "easy", 40))
    assert GQ_BINS == (0, 10, 20, 30, 40, 50)


def test_interval_brackets_rate() -> None:
    for errors, total in ((0, 10), (3, 10), (10, 10), (5, 1000), (1, 2_000_000_000)):
        lo, hi = jeffreys(errors, total)
        assert 0 <= lo <= errors / total <= hi <= 1, (errors, total)
    assert jeffreys(0, 10)[0] == 0.0 and jeffreys(10, 10)[1] == 1.0
    # Large counts: the normal approximation, an independent check of the quantiles.
    errors, total = 500, 10_000
    p = errors / total
    half = 1.959964 * math.sqrt(p * (1 - p) / total)
    assert jeffreys(errors, total) == pytest.approx((p - half, p + half), abs=2e-4)
    # The incomplete beta against its closed forms.
    for x in (0.001, 0.3, 0.9):
        assert betainc(1.0, 4.0, x) == pytest.approx(1 - (1 - x) ** 4, rel=1e-10)
        assert betainc(3.0, 1.0, x) == pytest.approx(x**3, rel=1e-10)
        assert betainc(0.5, 0.5, x) == pytest.approx(
            2 / math.pi * math.asin(math.sqrt(x)), rel=1e-10
        )
    # Each bound is the quantile it claims: 2.5% of Beta(e + 1/2, n - e + 1/2) below lo.
    lo, hi = jeffreys(3, 40)
    assert betainc(3.5, 37.5, lo) == pytest.approx(0.025, abs=1e-9)
    assert betainc(3.5, 37.5, hi) == pytest.approx(0.975, abs=1e-9)


def test_table_round_trips_with_versions(tmp_path: Path) -> None:
    counts = _counts({0: (5, 50), 50: (1, 5000)})
    counts[(ALLELE, CellKey("ref_block", "indel", "difficult", 50))] = Counts(
        3, 10**9, 2.5e10
    )
    table = CalibrationTable.from_counts(counts, VERSIONS, f1={"snp": 0.9971})
    folder = tmp_path / "calibration" / "table"
    table.write(folder)
    back = CalibrationTable.read(folder)
    assert back == table
    assert (back.versions, back.f1) == (VERSIONS, {"snp": 0.9971})
    manifest = json.loads((folder / "manifest.json").read_text())
    assert manifest["versions"] == VERSIONS
    # Written again over itself: replaced whole.
    table.write(folder)
    assert CalibrationTable.read(folder) == table
    # Each of the seven versions is required.
    for key in VERSION_KEYS:
        with pytest.raises(ValueError, match=key):
            CalibrationTable.from_counts(
                counts, {k: v for k, v in VERSIONS.items() if k != key}
            )
    # A cell edited after writing is refused.
    cells = folder / "cells.tsv"
    cells.write_text(cells.read_text().replace("\t5\t50\t", "\t4\t50\t"))
    with pytest.raises(ValueError, match="sha256"):
        CalibrationTable.read(folder)


def test_report_splits_odd_even_and_chr20_22() -> None:
    """The held-out checks: a table built on the odd chromosomes is checked on the even
    ones, and one built on chr1-19 on chr20-22. A held-out cell is outside when its errors
    fall outside what its training cell predicts for its units.

    Four cells, each in a group of its own (no pooling between them), place every
    chromosome: one measured on the even chromosomes only, one on chr20-22 only, one
    whose rate differs between odd and even chromosomes, and one whose rate differs
    between chr1-19 and chr20-22. Each split's training half sees the other split's rate
    jump between its chromosomes and widens that cell's interval by it. A split that took
    the wrong chromosomes, or swapped its halves, counts other cells."""
    autosomes = [f"chr{n}" for n in range(1, 23)]
    even_only = CellKey("variant", "snv", "easy", 50)
    tail_only = CellKey("RefCall", "snv", "easy", 50)
    parity = CellKey("ref_block", "indel", "difficult", 30)
    tail = CellKey("variant", "indel", "difficult", 30)
    by_chrom: dict[str, dict[tuple[str, CellKey], Counts]] = {}
    for chrom in autosomes:
        n = int(chrom.removeprefix("chr"))
        cells: dict[tuple[str, CellKey], Counts] = {}
        if n % 2 == 0:
            cells[(GENOTYPE, even_only)] = Counts(10, 10**6, 3.0e7)
        if n >= 20:
            cells[(GENOTYPE, tail_only)] = Counts(10, 10**6, 3.0e7)
        cells[(GENOTYPE, parity)] = Counts(0 if n % 2 else 1000, 10**6, 3.0e7)
        cells[(GENOTYPE, tail)] = Counts(1000 if n >= 20 else 0, 10**6, 3.0e7)
        by_chrom[chrom] = cells
    odd_even = held_out(
        by_chrom,
        train=[c for c in autosomes if int(c[3:]) % 2],
        test=[c for c in autosomes if not int(c[3:]) % 2],
    )
    assert (odd_even.checked, odd_even.outside, odd_even.unmeasured) == (
        3,
        [parity],
        [even_only],
    )
    apart = held_out(by_chrom, train=autosomes[:19], test=autosomes[19:])
    assert (apart.checked, apart.outside, apart.unmeasured) == (
        3,
        [tail],
        [tail_only],
    )
    # The report runs both splits on the chromosomes the spec names.
    merged: dict[tuple[str, CellKey], Counts] = {}
    for cells in by_chrom.values():
        for key, found in cells.items():
            merged.setdefault(key, Counts()).merge(found)
    f1 = {"snp": 0.996, "indel": 0.98}
    text = report(CalibrationTable.from_counts(merged, VERSIONS, f1=f1), by_chrom)
    # Task 21's stop rule on F1, read off the report.
    assert "\nSNP F1 0.9960 (the part stops below 0.995)\n" in text
    assert "\nINDEL F1 0.9800 (the part stops below 0.99); below the limit\n" in text
    assert (
        "held out, odd chromosomes -> even chromosomes: 3 cells checked, 1 outside"
        " their interval (33.3%), 1 not measured in the training half"
    ) in text
    assert (
        "held out, chr1-19 -> chr20-22: 3 cells checked, 1 outside their interval"
        " (33.3%), 1 not measured in the training half"
    ) in text
    # The table itself, one row per measured cell, and the versions measured on.
    assert text.count("\nvariant\tsnv\teasy\t50\t") == 1
    assert all(f"{key}: {value}" in text for key, value in VERSIONS.items())


# GRCh38's autosome lengths in Mb, and units per Mb and base error rate per row type.
_MB = (248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90)
_MB += (83, 80, 59, 64, 47, 51)
_ROWS = {"variant": (30, 0.05), "RefCall": (5, 0.1), "ref_block": (300, 0.005)}


def _simulated(
    seed: int, *, even_factor: float = 1.0, spread: float = 0.0
) -> dict[str, Cells]:
    """Each chromosome draws each of the 72 cells' errors from one true rate (falling
    with GQ), its units in proportion to its length; even chromosomes' rates times
    `even_factor`, and each difficult cell's rate on each chromosome times a lognormal
    factor of sd `spread` (its own mix of difficult sequence)."""
    draw = random.Random(seed)
    by_chrom: dict[str, Cells] = {}
    for n, mb in enumerate(_MB, start=1):
        cells: Cells = {}
        for (row_type, (per_mb, base)), indel, hard, gq in itertools.product(
            _ROWS.items(), (0, 1), (0, 1), GQ_BINS
        ):
            rate = base * (1 + indel) * (1 + 2 * hard) * 0.6 ** (gq / 10)
            rate *= 1 if n % 2 else even_factor
            rate = min(rate * math.exp(spread * draw.gauss() if hard else 0.0), 0.9)
            region = ("easy", "difficult")[hard]
            key = CellKey(row_type, ("snv", "indel")[indel], region, gq)
            units = mb * per_mb
            cells[(GENOTYPE, key)] = Counts(draw.binomialvariate(units, rate), units)
        by_chrom[f"chr{n}"] = cells
    return by_chrom


def _shares(by_chrom: dict[str, Cells]) -> list[float]:
    """The share of held-out cells outside their interval, per split of the report."""
    found = [held_out(by_chrom, train=train, test=test) for _, train, _, test in SPLITS]
    return [len(f.outside) / f.checked for f in found]


def test_predictive_interval_exact_and_at_large_totals() -> None:
    def log_beta(x: float, y: float) -> float:
        return math.lgamma(x) + math.lgamma(y) - math.lgamma(x + y)

    # Small n: the quantiles of the pmf summed from its closed form.
    n, a, b = 30, 2.5, 7.5
    cdf = list(
        itertools.accumulate(
            math.comb(n, k) * math.exp(log_beta(k + a, n - k + b) - log_beta(a, b))
            for k in range(n + 1)
        )
    )
    expected = [next(k for k, c in enumerate(cdf) if c >= q) for q in (0.025, 0.975)]
    assert list(beta_binomial_interval(n, a, b)) == expected
    # 10^9 units against a training cell of 10^6 errors in 10^9: the mean +- 1.96 sd of
    # the beta-binomial, sd^2 = n p (1 - p) (a + b + n) / (a + b + 1).
    n, a, b = 10**9, 10**6 + 0.5, 10**9 - 10**6 + 0.5
    p = a / (a + b)
    sd = math.sqrt(n * p * (1 - p) * (a + b + n) / (a + b + 1))
    lo, hi = beta_binomial_interval(n, a, b)
    assert (lo, hi) == pytest.approx(
        (n * p - 1.96 * sd, n * p + 1.96 * sd), abs=sd / 100
    )
    # No training errors: the interval starts at 0.
    assert beta_binomial_interval(10**8, 0.5, 10**8 + 0.5)[0] == 0


def test_held_out_check_passes_a_calibrated_table_and_fails_a_miscalibrated_one() -> (
    None
):
    """Both halves drawn from the same rates: about 5% of held-out cells fall outside, so
    the 10% limit holds in the large majority of draws. Rates half as high again on the
    even chromosomes: far above the limit in every draw."""
    seeds = range(20)
    calibrated = [_shares(_simulated(seed)) for seed in seeds]
    under = sum(all(s <= MAX_OUTSIDE for s in shares) for shares in calibrated)
    mean = sum(shares[0] for shares in calibrated) / len(seeds)
    assert under >= 16, calibrated
    assert mean < 0.08, mean
    for seed in seeds[:5]:
        odd_even, _ = _shares(_simulated(seed, even_factor=1.5))
        assert odd_even > 0.5, (seed, odd_even)


def test_held_out_check_passes_rates_that_vary_between_chromosomes() -> None:
    """Difficult cells whose true rate differs between chromosomes (sd 30%): binomial
    noise alone puts far more than 10% of held-out cells outside; with the dispersion the
    training half shows, about 5% fall outside, as for a table right on average."""
    seeds = range(20)
    shares = [_shares(_simulated(seed, spread=0.3)) for seed in seeds]
    under = sum(all(s <= MAX_OUTSIDE for s in found) for found in shares)
    mean = sum(found[0] for found in shares) / len(seeds)
    assert under >= 16, shares
    assert mean < 0.08, mean


def test_interval_widens_by_dispersion_between_chromosomes() -> None:
    """A pool's interval is Jeffreys' at its counts divided by Pearson's X^2 / df across
    chromosomes; its rate, and so every call's p_wrong, stays the pooled rate. Binomial
    data, one chromosome, or no errors leave the interval as it was."""
    gq10 = CellKey("ref_block", "snv", "difficult", 10)
    gq20 = CellKey("ref_block", "snv", "difficult", 20)  # worse: pooled with gq10
    rows = [(100, 400), (300, 400), (100, 400), (300, 400)]  # errors per chromosome
    by_chrom = {
        f"chr{n}": {
            (GENOTYPE, gq10): Counts(e10, 10**5, 3.0e6),
            (GENOTYPE, gq20): Counts(e20, 10**5, 3.0e6),
        }
        for n, (e10, e20) in enumerate(rows, start=1)
    }
    counts = merged(by_chrom)
    plain = CalibrationTable.from_counts(counts, VERSIONS)
    spread = CalibrationTable.from_counts(counts, VERSIONS, by_chrom=by_chrom)
    p = 2400 / 8e5
    phi = sum((e10 + e20 - 2e5 * p) ** 2 / (2e5 * p * (1 - p)) for e10, e20 in rows)
    phi /= len(rows) - 1
    for key in (gq10, gq20):
        wide, narrow = spread.rate_for(key), plain.rate_for(key)
        assert wide.rate == narrow.rate == p
        assert (wide.lo, wide.hi) == pytest.approx(jeffreys(2400 / phi, 8e5 / phi))
        assert (narrow.lo, narrow.hi) == pytest.approx(jeffreys(2400, 8e5))
        assert wide.hi - wide.lo > 4 * (narrow.hi - narrow.lo)  # sqrt(phi) = 4.7
    # Counts in proportion, one chromosome, no errors: binomial (1).
    assert dispersion([(10, 1000), (20, 2000), (30, 3000)]) == 1.0
    assert dispersion([(800, 4 * 10**5)]) == dispersion([(0, 100), (0, 500)]) == 1.0
    one = CalibrationTable.from_counts(counts, VERSIONS, by_chrom={"chr1": counts})
    assert one == plain
    # Binomial draws: X^2 / df floors at 1, so most intervals are unchanged and the
    # widest is widened by chance alone, about sqrt of chi-square's upper tail over df.
    ratios = []
    for seed in range(5):
        drawn = _simulated(seed)
        before = CalibrationTable.from_counts(merged(drawn), VERSIONS)
        after = CalibrationTable.from_counts(merged(drawn), VERSIONS, by_chrom=drawn)
        for key in before.measured():
            old, new = before.rate_for(key), after.rate_for(key)
            assert new.rate == old.rate
            ratios.append((new.hi - new.lo) / (old.hi - old.lo))
    assert sorted(ratios)[len(ratios) // 2] == pytest.approx(1.0)
    assert 1.0 - 1e-9 <= min(ratios) and max(ratios) < 1.6
