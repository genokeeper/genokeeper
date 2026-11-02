"""Slice 3b's deletion calibration (`genokeeper.calibrate.deletions`): the user's option D
arithmetic on counts built by hand, the table written whole and refused when missing or made
with other versions, gnomAD-SV read through bcftools from hand-written VCFs, and the run's
deletion steps (`pipeline.deletion_steps`) with their containers faked."""

import gzip
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import pysam
import pytest
from fixtures.genome_mini import FakeDeletions, lab_paths
from fixtures.grch38_mini import write_reference

from genokeeper import pipeline
from genokeeper.calibrate import hg002, tables
from genokeeper.calibrate.build import build_deletions
from genokeeper.calibrate.classes import (
    BandTerms,
    Geometry,
    Likelihoods,
    Share,
    geometry,
    share,
)
from genokeeper.calibrate.deletions import (
    FEATURE_VALUES,
    OUTSIDE,
    OUTSIDE_BENCHMARK,
    UNSEEN,
    DeletionTable,
    class_of,
    deletion_versions,
    deletions_dir,
    expected_true_carried,
    measure,
    open_deletions,
    read_population,
    summary,
    write_table,
)
from genokeeper.calibrate.hg002 import LIBRARIES, hg002_dir
from genokeeper.calibrate.table import jeffreys
from genokeeper.calibrate.tables import NOT_CALIBRATED
from genokeeper.checksum import sha256_file
from genokeeper.config import (
    DeletionThresholds,
    InsertionThresholds,
    Paths,
    grch38_fasta,
    load_grading,
)
from genokeeper.deletions import callers, truvari
from genokeeper.deletions.callers import caller_filter
from genokeeper.deletions.exons import Target
from genokeeper.deletions.files import sv_files, truth_files
from genokeeper.deletions.screen import (
    COMMON,
    CallClass,
    Deletion,
    PopulationSv,
    Screen,
    population_spans,
    rare,
    screen,
)
from genokeeper.deletions.screen import key as class_key
from genokeeper.deletions.truvari import read_union
from genokeeper.pipeline import RunRefused
from genokeeper.regions import giab38_dir
from genokeeper.tools import ToolMissing

T = DeletionThresholds()
# 4,000 targets of 250 bp, 100 kb apart: a 2 kb call's reach is 4,000 x 2,250 bp.
G = Geometry(
    targets=tuple(("chr1", n * 100_000, n * 100_000 + 250) for n in range(1, 4001)),
    region_bp=2_000_000_000,
)
LONG = CallClass("long", "both", False, False)
PRIORS = {"long": 0.063, "short": 0.01}
type Calls = list[tuple[CallClass, int, bool]]


def _ones() -> Likelihoods:
    """Every feature value's likelihood 1, for false and true calls."""
    one = {
        f: {v: Share(1.0, 1.0, 1.0) for v in values}
        for f, values in FEATURE_VALUES.items()
    }
    return Likelihoods(false=one, true={f: dict(v) for f, v in one.items()})


def _long_band() -> tuple[Calls, list[tuple[int, bool]]]:
    """One library's half of the long band: 20 false and 330 true 2 kb calls; 225 of 250
    truth deletions found."""
    calls = [(LONG, 2000, False)] * 20 + [(LONG, 2000, True)] * 330
    truth = [(2000, True)] * 225 + [(2000, False)] * 25
    return calls, truth


def _table(band: str, terms: BandTerms, likelihoods: Likelihoods) -> DeletionTable:
    return DeletionTable(
        versions={}, bands={band: terms}, likelihoods={band: likelihoods}, measured={}
    )


def test_band_p_wrong_option_d() -> None:
    calls, truth = _long_band()
    bands, _ = measure(
        {"pcrfree": calls, "pcr": calls}, {"pcrfree": truth, "pcr": truth}, G, PRIORS, T
    )
    band = bands["long"]
    assert (band.calls, band.libraries, band.prior) == (700, 2, 0.063)
    assert band.h == pytest.approx(0.0045)
    assert band.false_calls == share(40, 700)
    assert band.recall == share(450, 500)
    # The smoothing is its own test: 40/700 for the arithmetic.
    terms = replace(band, false_calls=replace(band.false_calls, point=40 / 700))
    table = _table("long", terms, _ones())
    false, true = table.expected(LONG, "point")
    assert false == pytest.approx(0.09)
    assert true == pytest.approx(0.0567, rel=2e-3)
    found = table.p_wrong(LONG)
    assert found is not None
    rate, lo, hi = found
    assert rate == pytest.approx(0.6135, abs=1e-3)
    assert lo < rate < hi


def test_reach_is_the_union_of_the_windows() -> None:
    """Two overlapping targets give their union, not their sum; a call longer than the gap
    between two targets counts the bases their windows share once; a window clips at 0;
    another chromosome's targets add their own."""
    overlapping = Geometry((("chr1", 1000, 1300),), 10_000)
    assert overlapping.reach(0) == 300
    assert overlapping.reach(100) == 400
    apart = Geometry((("chr1", 1000, 1100), ("chr1", 1300, 1400)), 10_000)
    assert apart.reach(50) == 150 + 150
    # Past the 200 bp gap: [700, 1100) and [1000, 1400) share 100 bases, not summed (800).
    assert apart.reach(300) == 700
    assert apart.reach(500) == 900
    assert Geometry((("chr1", 100, 200),), 10_000).reach(1000) == 200
    two = Geometry((("chr1", 1000, 1100), ("chr2", 1000, 1100)), 10_000)
    assert two.reach(2000) == 2 * 1100


def test_geometry_merges_overlapping_targets() -> None:
    """`geometry` merges the genes' targets, so h counts two overlapping targets' bases
    once."""
    found = geometry(
        {
            "A": [Target("chr1", 1000, 1200, (1,))],
            "B": [Target("chr1", 1100, 1300, (2,))],
        },
        [("chr1", 0, 10_000)],
    )
    assert found == Geometry((("chr1", 1000, 1300),), 10_000)
    calls = [(LONG, 2000, False)]
    bands, _ = measure({"pcrfree": calls}, {"pcrfree": []}, found, PRIORS, T)
    assert bands["long"].h == pytest.approx(1300 / 10_000)


def test_h_at_most_one() -> None:
    """A call longer than the benchmark's bases: h is at most 1."""
    g = Geometry((("chr1", 50_000, 60_000), ("chr2", 90_000, 95_000)), 100_000)
    calls = [(LONG, 10**9, False)]
    bands, _ = measure({"pcrfree": calls}, {"pcrfree": []}, g, PRIORS, T)
    assert bands["long"].h == 1.0


def test_class_likelihood_moves_p_wrong() -> None:
    calls, truth = _long_band()
    bands, _ = measure(
        {"pcrfree": calls, "pcr": calls}, {"pcrfree": truth, "pcr": truth}, G, PRIORS, T
    )
    terms = replace(
        bands["long"], false_calls=replace(bands["long"].false_calls, point=40 / 700)
    )
    likelihoods = _ones()
    likelihoods.false["evidence"]["dysgu"] = Share(0.5, 0.5, 0.5)
    likelihoods.true["evidence"]["dysgu"] = Share(0.1, 0.1, 0.1)
    likelihoods.false["evidence"]["both+depth"] = Share(0.25, 0.25, 0.25)
    likelihoods.true["evidence"]["both+depth"] = Share(0.9, 0.9, 0.9)
    table = _table("long", terms, likelihoods)
    dysgu = CallClass("long", "dysgu", False, False)
    false, true = table.expected(dysgu, "point")
    assert false == pytest.approx(0.045)
    assert true == pytest.approx(0.00567, rel=2e-3)
    found = table.p_wrong(dysgu)
    assert found is not None and found[0] == pytest.approx(0.888, abs=1e-3)
    found = table.p_wrong(CallClass("long", "both+depth", False, False))
    assert found is not None and found[0] == pytest.approx(0.306, abs=1e-3)


def test_sparse_class_wide_interval() -> None:
    """40 false and 660 true calls: `delly` seen in 1 false and 1 true, `both` in 20 false
    and 330 true; the sparse class's interval is the wider, and wider than 0.5."""

    def calls(evidence: str, n: int, confirmed: bool) -> Calls:
        return [(CallClass("long", evidence, False, False), 2000, confirmed)] * n

    library = [
        *calls("delly", 1, False), *calls("both", 20, False), *calls("dysgu", 19, False),
        *calls("delly", 1, True), *calls("both", 330, True), *calls("dysgu", 329, True),
    ]  # fmt: skip
    _, truth = _long_band()
    bands, likelihoods = measure(
        {"pcrfree": library, "pcr": []}, {"pcrfree": truth, "pcr": truth}, G, PRIORS, T
    )
    assert likelihoods["long"].false["evidence"]["delly"] == share(1, 40)
    assert likelihoods["long"].true["evidence"]["both"] == share(330, 660)
    table = DeletionTable({}, bands, likelihoods, {})
    sparse = table.p_wrong(CallClass("long", "delly", False, False))
    dense = table.p_wrong(CallClass("long", "both", False, False))
    assert sparse is not None and dense is not None
    assert sparse[2] - sparse[1] > dense[2] - dense[1]
    assert sparse[2] - sparse[1] > 0.5


def test_share_smoothing() -> None:
    none = share(0, 40)
    assert none.point == pytest.approx(0.5 / 41) and none.point > 0
    assert (none.lo, none.hi) == jeffreys(0, 40)
    assert share(0, 0) == Share(0.5, 0.0, 1.0)


def test_band_without_calls_not_measured() -> None:
    calls, truth = _long_band()
    bands, likelihoods = measure({"pcrfree": calls}, {"pcrfree": truth}, G, PRIORS, T)
    assert bands["short"].calls == 0 and bands["short"].h is None
    table = DeletionTable({}, bands, likelihoods, {})
    for evidence in FEATURE_VALUES["evidence"]:
        for repeat in (True, False):
            for masked in (True, False):
                assert (
                    table.p_wrong(CallClass("short", evidence, repeat, masked)) is None
                )
    assert table.p_wrong(LONG) is not None


def test_recall_per_band() -> None:
    calls, _ = _long_band()
    truth = [(2000, True)] * 9 + [(2000, False)] + [(200, False)] * 5 + [(40, True)] * 3
    bands, _ = measure({"pcrfree": calls}, {"pcrfree": truth}, G, PRIORS, T)
    assert bands["long"].recall == share(9, 10)
    assert bands["short"].recall == share(0, 5)
    more = truth + [(300, False)] * 20
    again, _ = measure({"pcrfree": calls}, {"pcrfree": more}, G, PRIORS, T)
    assert again["long"].recall == bands["long"].recall
    assert again["short"].recall == share(0, 25)


def test_recall_counts_per_library() -> None:
    """`DeletionTable.recall`: per band, each library's truth deletions found and their
    number from its `summary`, never summed over the one genome's libraries; a band no
    library measured is left out."""
    calls, _ = _long_band()
    one = [(2000, True)] * 9 + [(2000, False)] + [(200, True)] * 3 + [(200, False)] * 4
    two = [(2000, True)] * 5 + [(2000, False)] * 5
    measured = {"pcrfree": summary(calls, one, T), "pcr": summary(calls, two, T)}
    table = DeletionTable({}, {}, {}, measured)
    assert table.recall() == {
        "long": {"pcr": (5, 10), "pcrfree": (9, 10)},
        "short": {"pcr": (0, 0), "pcrfree": (3, 7)},
    }
    assert DeletionTable(
        {}, {}, {}, {"pcr": {"long": measured["pcr"]["long"]}}
    ).recall() == {"long": {"pcr": (5, 10)}}
    assert DeletionTable({}, {}, {}, {}).recall() == {}


VERSIONS = {
    "delly": "genomelab/delly:2.7.0-a",
    "dysgu": "genomelab/dysgu:1.9.0-b",
    "truvari": "genomelab/truvari:5.5.0-c",
    "between_callers": "--pctseq 0",
    "fasta": "f" * 64,
}


def _measured_table(measured: object = None) -> DeletionTable:
    calls, truth = _long_band()
    bands, likelihoods = measure({"pcrfree": calls}, {"pcrfree": truth}, G, PRIORS, T)
    return DeletionTable(
        dict(VERSIONS), bands, likelihoods, {"pcrfree": measured or {"calls": 350}}
    )


def test_cautions_name_an_unseen_class_and_a_target_outside_the_benchmark() -> None:
    """A class no library's calibration call had, and a deletion over a target the
    benchmark leaves out, each carry their line; p_wrong itself stays as it is."""
    seen = CallClass("long", "both+depth", False, False)
    measured = {"pcrfree": {"long": {"classes": [class_key(seen)]}, "short": {}}}
    table = DeletionTable(
        {}, {}, {}, measured, extra={OUTSIDE_BENCHMARK: [["chr1", 1000, 1200]]}
    )
    inside = Deletion("chr1", 2000, 3000, "N", "0/1", {}, None)
    over = Deletion("chr1", 1100, 2100, "N", "0/1", {}, None)
    assert table.cautions(inside, seen) == []
    assert table.cautions(inside, LONG) == [UNSEEN]
    assert table.cautions(over, seen) == [OUTSIDE]
    assert (
        table.cautions(Deletion("chr2", 1100, 2100, "N", "0/1", {}, None), seen) == []
    )


def test_open_deletions_missing_and_stale(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    assert open_deletions(lab, VERSIONS) == (None, NOT_CALIBRATED)
    assert NOT_CALIBRATED.endswith("run `genokeeper calibrate`")
    _put(lab, _measured_table())
    table, reason = open_deletions(lab, VERSIONS)
    assert reason is None and table == _measured_table()
    for key in ("dysgu", "between_callers"):
        other = VERSIONS | {key: "another"}
        table, reason = open_deletions(lab, other)
        assert table is None and reason is not None
        assert reason == f"stale: {key} changed; run `genokeeper calibrate`"


def _put(lab: Path, table: DeletionTable) -> None:
    tables.write_table(
        lab, "deletions", lambda s: write_table(s, table), table.versions
    )


def test_write_table_atomic(tmp_path: Path) -> None:
    folder = deletions_dir(tmp_path / "lab")
    broken = _measured_table(object())
    with pytest.raises(TypeError):
        _put(tmp_path / "lab", broken)
    assert not folder.exists() and list(folder.parent.iterdir()) == []
    _put(tmp_path / "lab", _measured_table())
    with pytest.raises(TypeError):
        _put(tmp_path / "lab", broken)
    assert open_deletions(tmp_path / "lab", VERSIONS) == (_measured_table(), None)
    assert [p.name for p in folder.parent.iterdir()] == [folder.name]


# --- gnomAD-SV through bcftools -----------------------------------------------------------

type Record = tuple[str, int, int, str, str, int, float]


def _gnomad(
    path: Path, records: Sequence[Record], svlens: Mapping[str, int] | None = None
) -> Path:
    """A gnomAD-SV-shaped sites VCF of (id, pos, end, SVTYPE, FILTER, AN, AF) on chr20,
    with a group frequency beside the overall one, `SVLEN` from `svlens` by id (else `end -
    pos`); bgzipped and indexed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr20,length=64444167>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=UNRESOLVED,Description="Unresolved">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Alleles">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Frequency">',
        '##INFO=<ID=GRPMAX_AF,Number=A,Type=Float,Description="Group maximum">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for vid, pos, end, svtype, filt, an, af in sorted(records, key=lambda r: r[1]):
        svlen = (svlens or {}).get(vid, end - pos)
        info = f"END={end};SVTYPE={svtype};SVLEN={svlen};AN={an};AF={af};GRPMAX_AF=0.4"
        lines.append(f"chr20\t{pos}\t{vid}\tN\t<{svtype}>\t999\t{filt}\t{info}")
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return path


def _strata(lab: Path) -> None:
    giab = giab38_dir(lab)
    giab.mkdir(parents=True, exist_ok=True)
    for name in load_grading().strata.values():
        (giab / name).write_bytes(gzip.compress(b"chr20\t100\t200\n"))


def test_load_screen_population_filters(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    _strata(lab)
    _gnomad(
        sv_files(lab).gnomad,
        [
            # Starts before the called deletion [2000, 3000) and spans it.
            ("kept", 900, 5000, "DEL", "PASS", 2000, 0.02),
            ("an1999", 1500, 3500, "DEL", "PASS", 1999, 0.02),
            ("filtered", 1500, 3500, "DEL", "UNRESOLVED", 5000, 0.02),
            ("dup", 1500, 3500, "DUP", "PASS", 5000, 0.02),
            ("elsewhere", 50_000, 52_000, "DEL", "PASS", 5000, 0.02),
            ("after", 3200, 4000, "DEL", "PASS", 5000, 0.02),
        ],
    )
    called = [Deletion("chr20", 2000, 3000, "N", "0/1", {"dysgu": "1 QUAL 9"}, None)]
    s = pipeline.load_screen(lab, load_grading(), called, tmp_path)
    assert list(s.population) == [PopulationSv("kept", "chr20", 900, 5000, 0.02)]
    assert s.mask == {"chr20": [("chr20", 100, 200)]}
    assert s.repeats == {"chr20": [("chr20", 100, 200)]}
    assert "BRCA1" in s.targets


def test_regions_under_work(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The regions file goes under the folder the caller gives, not the system's temp."""
    seen: list[Path] = []

    def query(args: list[str]) -> str:
        seen.append(Path(args[args.index("-R") + 1]))
        return ""

    monkeypatch.setattr("genokeeper.calibrate.deletions.run_tool", query)
    read_population(tmp_path / "g.vcf.gz", [("chr20", 1, 2)], 2000, tmp_path)
    assert [p.parent.parent for p in seen] == [tmp_path]


def test_expected_true_carried(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    _gnomad(
        sv_files(lab).gnomad,
        [
            ("long1", 10_000, 12_000, "DEL", "PASS", 5000, 0.001),
            ("long2", 10_500, 13_000, "DEL", "PASS", 5000, 0.002),
            ("short", 11_000, 11_200, "DEL", "PASS", 5000, 0.004),
            ("common", 10_000, 12_500, "DEL", "PASS", 5000, 0.05),
            ("filtered", 10_000, 12_000, "DEL", "UNRESOLVED", 5000, 0.001),
            ("an1500", 10_000, 12_000, "DEL", "PASS", 1500, 0.001),
            ("away", 40_000, 42_000, "DEL", "PASS", 5000, 0.003),
        ],
    )
    s = Screen(
        targets={"GENEA": [Target("chr20", 11_100, 11_150, (1,))]},
        mask={},
        repeats={},
        population=[],
        thresholds=T,
    )
    grading = load_grading()
    assert expected_true_carried(lab, grading, s, "long") == pytest.approx(0.006)
    assert expected_true_carried(lab, grading, s, "short") == pytest.approx(0.008)


def test_read_population_svtype(tmp_path: Path) -> None:
    """Each type reads only its own records; an insertion its `SVLEN`, at its point."""
    vcf = _gnomad(
        tmp_path / "g.vcf.gz",
        [
            ("del", 1000, 3000, "DEL", "PASS", 5000, 0.02),
            ("dup", 1500, 4000, "DUP", "PASS", 5000, 0.03),
            ("ins", 2000, 2001, "INS", "PASS", 5000, 0.04),
        ],
        svlens={"ins": 700},
    )
    spans = [("chr20", 0, 10_000)]
    read = {
        svtype: read_population(vcf, spans, 2000, tmp_path, svtype=svtype)
        for svtype in ("DEL", "DUP", "INS")
    }
    assert read == {
        "DEL": [PopulationSv("del", "chr20", 1000, 3000, 0.02)],
        "DUP": [PopulationSv("dup", "chr20", 1500, 4000, 0.03, svtype="DUP")],
        "INS": [
            PopulationSv("ins", "chr20", 2000, 2001, 0.04, svtype="INS", svlen=700)
        ],
    }


def _ins_gnomad(path: Path, *, with_a: bool) -> Path:
    """A in the padded span (size ratio 0.93), B outside it, C inside (ratio 0.33), D a
    deletion over the call; all PASS, AN 3000, AF 0.05."""
    records: list[Record] = [
        ("A", 10_400, 10_401, "INS", "PASS", 3000, 0.05),
        ("B", 10_600, 10_601, "INS", "PASS", 3000, 0.05),
        ("C", 9_800, 9_801, "INS", "PASS", 3000, 0.05),
        ("D", 9_900, 10_200, "DEL", "PASS", 3000, 0.05),
    ]
    kept = [r for r in records if with_a or r[0] != "A"]
    return _gnomad(path, kept, svlens={"A": 280, "B": 300, "C": 100})


def test_insertion_common_by_window(tmp_path: Path) -> None:
    """Through the real reader: the call's span padded by `max_distance` reads A and C,
    never B (outside) or D (another type); A makes it common, C alone does not."""
    t = InsertionThresholds()
    call = Deletion(
        "chr20", 10_000, 10_001, "A", "0/1", {"dysgu": "1 QUAL 9"}, None,
        svtype="INS", svlen=300,
    )  # fmt: skip
    spans = population_spans([call], pad=t.max_distance)
    assert spans == [("chr20", 9_500, 10_501)]

    def screened(with_a: bool) -> tuple[list[PopulationSv], Screen]:
        vcf = _ins_gnomad(tmp_path / f"a{with_a}" / "g.vcf.gz", with_a=with_a)
        population = read_population(vcf, spans, 2000, tmp_path, svtype="INS")
        s = Screen(
            targets={"GENEA": [Target("chr20", 9_990, 10_010, (1,))]},
            mask={},
            repeats={},
            population=population,
            thresholds=t,
            callers=("dysgu", "insurveyor"),
        )
        return population, s

    population, s = screened(True)
    assert [p.id for p in population] == ["C", "A"]
    assert [p.svlen for p in population] == [100, 280]
    assert all(p.end == p.pos + 1 for p in population)
    assert screen(call, s) == ({}, COMMON)
    assert rare([call], population, t) == []
    population, s = screened(False)
    assert [p.id for p in population] == ["C"]
    assert screen(call, s) == ({"GENEA": (1,)}, None)
    assert rare([call], population, t) == [call]


def test_expected_true_insertions(tmp_path: Path) -> None:
    """An insertion's length is its `SVLEN`, its span its point: neither falls under the
    floor as `end - pos` would put it; a deletion over the target is not an insertion."""
    lab = tmp_path / "lab"
    _gnomad(
        sv_files(lab).gnomad,
        [
            ("short", 11_120, 11_121, "INS", "PASS", 5000, 0.001),
            ("long", 11_130, 11_131, "INS", "PASS", 5000, 0.001),
            ("away", 30_000, 30_001, "INS", "PASS", 5000, 0.004),
            ("del", 10_000, 12_000, "DEL", "PASS", 5000, 0.003),
        ],
        svlens={"short": 300, "long": 1500, "away": 300},
    )
    s = Screen(
        targets={"GENEA": [Target("chr20", 11_100, 11_150, (1,))]},
        mask={},
        repeats={},
        population=[],
        thresholds=T,
    )
    grading = load_grading()
    for band in ("short", "long"):
        found = expected_true_carried(lab, grading, s, band, svtype="INS")
        assert found == pytest.approx(0.002)
    assert expected_true_carried(lab, grading, s, "long") == pytest.approx(0.006)


def test_key_roundtrip_new_fields() -> None:
    every = [
        CallClass(band, evidence, repeat, masked, junction)
        for band in ("short", "long", "loss_1-5kb", "gain_large", "mei")
        for evidence in (
            "both",
            "both+depth",
            "dysgu",
            "delly",
            "dellycnv",
            "insurveyor",
        )
        for repeat in (False, True)
        for masked in (False, True)
        for junction in (False, True)
    ]
    assert all(class_of(class_key(c)) == c for c in every)
    assert len({class_key(c) for c in every}) == len(every)
    assert class_key(CallClass("long", "both", False, False)) == (
        "long:both:norepeat:unmasked:nojunction"
    )
    with pytest.raises(ValueError, match="not a class key"):
        class_of("long:both:norepeat:unmasked")
    assert FEATURE_VALUES["junction"] == ("true", "false")


# --- the run's deletion steps ---------------------------------------------------------------


def _steps(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fake: FakeDeletions,
    *,
    male: bool = True,
) -> tuple[pipeline.DeletionResults, dict[str, pipeline._Step]]:  # pyright: ignore[reportPrivateUsage]
    fasta = write_reference(tmp_path / "ref" / "g.fa")
    cram = tmp_path / "in" / "sample.cram"
    cram.parent.mkdir(exist_ok=True)
    cram.write_bytes(b"")
    for module in (callers, truvari):
        monkeypatch.setattr(module, "run_container", fake)
    return pipeline.deletion_steps(
        tmp_path / "work",
        cram=cram,
        cram_key="k" * 64,
        fasta=fasta,
        fasta_sha="f" * 64,
        grading=load_grading(),
        threads=2,
        male=male,
    )


def test_deletion_steps_failure_recorded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeDeletions()
    fake.fail.add("delly")
    results, steps = _steps(tmp_path, monkeypatch, fake)
    assert results.called is None
    assert results.failed == {"delly": "exit 1: delly crashed"}
    assert set(steps) == {"dysgu"}
    assert "truvari" not in fake.runs


@pytest.mark.parametrize("step", ["deletions_only", "match_callers"])
@pytest.mark.parametrize(
    "error",
    [ToolMissing("bcftools"), OSError(28, "No space left on device"), None],
    ids=["missing tool", "full disk", "wrote no file"],
)
def test_deletion_steps_failure_leaves_the_screen_unassessed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    step: str,
    error: Exception | None,
) -> None:
    """A missing tool, a full disk or a step that wrote no file, in a caller step or the
    union, is a failed step, not an aborted run (`None`: the step returns an empty folder,
    writing nothing)."""
    empty = tmp_path / "empty"

    def fail(*_: object, **__: object) -> Path:
        if error is not None:
            raise error
        empty.mkdir(exist_ok=True)
        return empty

    monkeypatch.setattr(pipeline, step, fail)
    results, _ = _steps(tmp_path, monkeypatch, FakeDeletions())
    assert results.called is None and results.failed


def test_deletion_steps_union(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fake = FakeDeletions()
    results, steps = _steps(tmp_path, monkeypatch, fake)
    assert results.failed == {}
    assert results.called is not None and len(results.called) == 57 + 43 + 30
    assert set(steps) == {"delly", "dysgu", "deletions_union"}
    # Each caller one container per primary chromosome of the mini FASTA (chr1, chrX, chrY).
    assert fake.runs == {"delly": 3, "dysgu": 3, "truvari": 1}
    male = {"chr1", "chrX", "chrY"}
    assert fake.shards == {"delly": male, "dysgu": male}
    again, cached = _steps(tmp_path, monkeypatch, fake)
    assert again == results and all(s.reused for s in cached.values())
    assert fake.runs == {"delly": 3, "dysgu": 3, "truvari": 1}


def test_deletion_steps_female_runs_no_chry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A female sample's callers run no chrY shard (DELLY fails on one with almost no
    reads); the union is read as for a male one."""
    fake = FakeDeletions()
    results, _ = _steps(tmp_path, monkeypatch, fake, male=False)
    assert results.failed == {} and results.called is not None
    assert fake.shards == {"delly": {"chr1", "chrX"}, "dysgu": {"chr1", "chrX"}}
    assert fake.runs == {"delly": 2, "dysgu": 2, "truvari": 1}


def test_deletion_steps_key_a_female_sample_apart(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A sex given rather than read from the CRAM cannot reuse the other sex's shards: a
    female sample's caller steps run again over a male one's cache, which stays reused."""
    fake = FakeDeletions()
    _, male = _steps(tmp_path, monkeypatch, fake)
    _, female = _steps(tmp_path, monkeypatch, fake, male=False)
    assert not any(s.reused for s in female.values())
    assert female["delly"].key != male["delly"].key
    _, again = _steps(tmp_path, monkeypatch, fake)
    assert all(s.reused for s in again.values())


def test_calibrate_deletions_refuses_a_missing_input(tmp_path: Path) -> None:
    paths = lab_paths(tmp_path)
    truth, _ = truth_files(paths.lab)
    with pytest.raises(RunRefused, match="run scripts/bootstrap.sh --with-hg002") as e:
        build_deletions(paths)
    assert str(truth) in str(e.value)
    assert not deletions_dir(paths.lab).exists()


def _hg002_lab(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, truth_lengths: Sequence[int]
) -> Paths:
    paths = lab_paths(tmp_path)
    lab = paths.lab
    _strata(lab)
    _gnomad(sv_files(lab).gnomad, [("away", 100, 3000, "DEL", "PASS", 5000, 0.5)])
    fasta = grch38_fasta(lab)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text(">chr20\nN\n")
    Path(f"{fasta}.fai").write_text("chr20\t64444167\t7\t60\t61\n")
    truth, bench = truth_files(lab)
    bench.parent.mkdir(parents=True, exist_ok=True)
    bench.write_text("chr20\t0\t64444167\n")
    records = [
        f"chr20\t{pos}\tt{pos}\tN\t<DEL>\t.\tPASS\tEND={pos + n};SVTYPE=DEL\tGT\t0/1"
        for pos, n in zip((5_000_000, 6_000_000, 7_000_000, 8_000_000), truth_lengths)
    ]
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
    )
    truth.write_bytes(gzip.compress((header + "\n".join(records) + "\n").encode()))
    for library in LIBRARIES:
        folder = hg002_dir(lab) / library
        folder.mkdir(parents=True)
        cram = folder / f"{library}.cram"
        cram.write_bytes(b"")
        record = hg002.Record(
            cram=str(cram),
            key=library * 8,
            metrics={},
            gvcf=None,
            reference_sha256=sha256_file(fasta),
            accepted=True,
            snp_f1=None,
            replaced=None,
        )
        hg002.write_record(lab, library, record)
    fake = FakeDeletions()
    for module in (callers, truvari):
        monkeypatch.setattr(module, "run_container", fake)
    return paths


def test_calibrate_deletions_writes_a_table_the_run_accepts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Both libraries through the deletion steps (faked), their rare calls of
    `min_length` or more against the truth (faked: every call confirmed), into a table
    `open_deletions` accepts under the lab's own versions."""
    paths = _hg002_lab(tmp_path, monkeypatch, (2000, 200, 40))
    lab = paths.lab
    table = build_deletions(paths, threads=2)
    (both,) = (hg002_dir(lab) / "pcrfree").glob("deletions_union/*/both")
    lengths = [d.length for d in read_union(both) if d.length >= T.min_length]
    long = sum(n >= T.long_from for n in lengths)
    assert 0 < long < len(lengths)
    assert table.bands["long"].calls == 2 * long
    assert table.bands["short"].calls == 2 * (len(lengths) - long)
    assert table.bands["long"].false_calls == share(0, 2 * long)
    assert table.bands["long"].recall == share(2, 2)
    assert table.bands["short"].recall == share(2, 2)
    assert set(table.measured) == set(LIBRARIES)
    # The benchmark is chr20 whole: every target elsewhere lies outside it.
    outside = {chrom for chrom, _, _ in table.extra[OUTSIDE_BENCHMARK]}
    assert "chr1" in outside and "chr20" not in outside
    versions = deletion_versions(lab, sha256_file(grch38_fasta(lab)), load_grading())
    assert open_deletions(lab, versions) == (table, None)


def test_versions_leave_out_the_wording_rule(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`likely_below` words a shown deletion and is never read by the table, so changing it
    leaves the versions alone; a threshold the table reads changes them."""
    lab = _hg002_lab(tmp_path, monkeypatch, (2000,)).lab
    grading, fasta = load_grading(), "f" * 64
    d = grading.deletions
    versions = deletion_versions(lab, fasta, grading)
    wording = replace(grading, deletions=replace(d, likely_below=0.4))
    assert deletion_versions(lab, fasta, wording) == versions
    assert "likely_below" not in versions["deletions"]
    assert (versions["delly_filter"], versions["dysgu_filter"]) == (
        caller_filter("delly", d.min_length),
        caller_filter("dysgu", d.min_length),
    )
    floor = replace(grading, deletions=replace(d, max_depth_ratio=0.6))
    assert deletion_versions(lab, fasta, floor)["deletions"] != versions["deletions"]


def test_calibrate_deletions_looks_up_only_truth_deletions_truvari_compares(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A truth deletion under the truth comparison's `--sizemin` (35) never reaches a
    result, so the population lookup is not asked about it."""
    paths = _hg002_lab(tmp_path, monkeypatch, (2000, 200, 40, 30))
    asked: list[Sequence[tuple[str, int, int]]] = []
    read = pipeline.read_population

    def spy(
        vcf: Path, spans: Sequence[tuple[str, int, int]], *rest: Any, **kw: Any
    ) -> list[PopulationSv]:
        asked.append(spans)
        return read(vcf, spans, *rest, **kw)

    monkeypatch.setattr(pipeline, "read_population", spy)
    build_deletions(paths, threads=2)
    assert len(asked) == len(LIBRARIES)
    for spans in asked:
        assert ("chr20", 7_000_000, 7_000_040) in spans
        assert not [s for s in spans if s[0] == "chr20" and s[1] == 8_000_000]
