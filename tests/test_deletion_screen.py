"""The deletion screen (`deletions.screen`): its three drop reasons and the call classes,
on targets, regions and deletions built by hand."""

from collections.abc import Mapping, Sequence

import pytest

from genokeeper.config import DeletionThresholds, InsertionThresholds
from genokeeper.deletions.exons import Exon, Target, targets
from genokeeper.deletions.screen import (
    COMMON,
    NO_DEPTH,
    NO_EXON,
    SAME_POSITION,
    SHORT,
    SHORT_DEPTH,
    CallClass,
    Deletion,
    PopulationSv,
    Screen,
    classify,
    is_common,
    key,
    population_match,
    population_spans,
    rare,
    screen,
    screen_all,
    warnings,
)
from genokeeper.regions import Interval

T = DeletionThresholds()
# BRCA1-like gene on chr17: exons 1..3 at 0-based [1000, 1100), [2000, 2100), [3000, 3100).
GENE_A = [
    Exon("GENEA", "TA.1", n, "chr17", s, s + 100)
    for n, s in enumerate((1000, 2000, 3000), 1)
]
# A second gene further on, exons 1..2 at [5000, 5100) and [6000, 6100).
GENE_B = [
    Exon("GENEB", "TB.1", n, "chr17", s, s + 100) for n, s in enumerate((5000, 6000), 1)
]
BOTH = {"dysgu": "1 QUAL 45.0", "delly": "DEL00001 QUAL 300"}


def _screen(
    *,
    mask: Mapping[str, Sequence[Interval]] | None = None,
    repeats: Mapping[str, Sequence[Interval]] | None = None,
    population: Sequence[PopulationSv] = (),
) -> Screen:
    return Screen(
        targets=targets({"GENEA": GENE_A, "GENEB": GENE_B}, T.exon_margin),
        mask=mask or {},
        repeats=repeats or {},
        population=population,
        thresholds=T,
    )


def _deletion(
    pos: int,
    end: int,
    *,
    callers: Mapping[str, str] = BOTH,
    depth_ratio: float | None = 0.4,
    chrom: str = "chr17",
) -> Deletion:
    return Deletion(chrom, pos, end, "A", "0/1", callers, depth_ratio)


def test_short_dropped() -> None:
    s = _screen()
    assert screen(_deletion(1020, 1069), s) == ({}, SHORT)
    fifty = _deletion(1020, 1070)
    assert fifty.length == 50
    assert screen(fifty, s) == ({"GENEA": (1,)}, None)
    assert classify(fifty, s).band == "short"
    assert classify(_deletion(1020, 2020), s).band == "long"


def test_no_exon_dropped() -> None:
    s = _screen()
    assert screen(_deletion(1300, 1800), s) == ({}, NO_EXON)
    # [1915, 1985) ends 15 bp before exon 2 at 2000; its target starts at 1980.
    assert screen(_deletion(1915, 1985), s) == ({"GENEA": (2,)}, None)


def test_masked_end_is_a_class() -> None:
    mask = {"chr17": [("chr17", 500, 901)]}
    s = _screen(mask=mask)
    first_in = _deletion(
        900, 1500
    )  # first deleted base, 0-based 900, is the segdup's last
    assert screen(first_in, s) == ({"GENEA": (1,)}, None)
    c = classify(first_in, s)
    assert c.masked is True
    assert "an end in a segmental duplication or low-mappability region" in warnings(c)
    last_in = _deletion(
        200, 501
    )  # the last deleted base, 0-based 500, is the segdup's first
    s2 = _screen(mask={"chr17": [("chr17", 500, 1100)]})
    assert classify(last_in, s2).masked is True
    outside = _deletion(400, 1050)
    # END 1050 sits on the next segdup's first base: the base after the deletion.
    assert (
        classify(outside, _screen(mask={"chr17": [("chr17", 1050, 1200)]})).masked
        is False
    )


def test_repeat_end_short_band_only() -> None:
    repeats = {"chr17": [("chr17", 1290, 1310)]}
    s = _screen(repeats=repeats)
    short = _deletion(1000, 1300)
    assert classify(short, s).repeat is True
    assert "an end in a tandem repeat: may be a repeat-length change" in warnings(
        classify(short, s)
    )
    long = _deletion(1299, 4299)
    assert long.length == 3000
    assert classify(long, s).repeat is False


def test_evidence_classes() -> None:
    s = _screen()
    span = (1000, 2500)
    delly = {"delly": "DEL00002 QUAL 300"}
    dysgu = {"dysgu": "7 QUAL 30.0"}
    cases = {
        "both+depth": _deletion(*span, depth_ratio=0.4),
        "both": _deletion(*span, depth_ratio=0.8),
        "delly+depth": _deletion(*span, callers=delly, depth_ratio=0.4),
        "delly": _deletion(*span, callers=delly, depth_ratio=0.8),
        "dysgu": _deletion(*span, callers=dysgu, depth_ratio=None),
    }
    got = {name: classify(d, s) for name, d in cases.items()}
    assert {name: c.evidence for name, c in got.items()} == {n: n for n in cases}
    assert classify(_deletion(*span, depth_ratio=None), s).evidence == "both"
    assert warnings(got["both+depth"]) == []
    assert warnings(got["both"]) == ["read depth does not support it"]
    assert warnings(got["delly+depth"]) == ["found by one caller (DELLY)"]
    assert warnings(got["delly"]) == [
        "found by one caller (DELLY)",
        "read depth does not support it",
    ]
    assert warnings(got["dysgu"]) == [
        "found by one caller (Dysgu)",
        "read depth not measured",
    ]
    every = [
        CallClass(band, evidence, repeat, masked)
        for band in ("short", "long")
        for evidence in cases
        for repeat in (False, True)
        for masked in (False, True)
    ]
    assert len({key(c) for c in every}) == len(every)


def test_short_band_depth_line_is_not_evidence_against() -> None:
    """Under `long_from` read depth rarely moves, so a short deletion without depth
    support is told so, not that depth fails to support it."""
    assert warnings(CallClass("short", "both", False, False)) == [SHORT_DEPTH]
    assert warnings(CallClass("short", "delly", False, False)) == [
        "found by one caller (DELLY)",
        SHORT_DEPTH,
    ]
    assert warnings(CallClass("long", "both", False, False)) == [NO_DEPTH]
    assert warnings(CallClass("short", "both+depth", False, False)) == []


def test_flank_only_overlap_named_per_gene() -> None:
    """GENEA's exon 2 is [2000, 2100), its target [1980, 2120): a deletion of [1915, 1985)
    touches the margin only, one of [1915, 2001) the exon's first base."""
    s = _screen()
    (flank,), _ = screen_all([_deletion(1915, 1985)], s)
    (exon,), _ = screen_all([_deletion(1915, 2001)], s)
    assert flank.flank == frozenset({"GENEA"})
    assert exon.flank == frozenset()
    # Over GENEA's exon 3 and only the margin before GENEB's exon 1 at 5000.
    (both,), _ = screen_all([_deletion(3050, 4985)], s)
    assert both.exons == {"GENEA": (3,), "GENEB": (1,)}
    assert both.flank == frozenset({"GENEB"})


def test_whole_gene_deletion_shown() -> None:
    assert screen(_deletion(500, 3500), _screen()) == ({"GENEA": (1, 2, 3)}, None)


def test_common_dropped_and_rare_shown() -> None:
    d = _deletion(1000, 2000)
    # Reciprocal overlap 600/1000 = 0.6.
    sixty = PopulationSv("gnomAD-SV_v3_DEL_chr17_1", "chr17", 1400, 2000, 0.02)
    assert screen(d, _screen(population=[sixty])) == ({}, COMMON)
    assert is_common(d, _screen(population=[sixty]))
    at_cut = PopulationSv("x", "chr17", 1400, 2000, T.max_gnomad_af)
    assert screen(d, _screen(population=[at_cut])) == ({"GENEA": (1, 2)}, None)
    forty = PopulationSv("y", "chr17", 1600, 2000, 0.5)
    assert screen(d, _screen(population=[forty])) == ({"GENEA": (1, 2)}, None)
    assert not is_common(d, _screen(population=[forty, at_cut]))


def test_population_index_finds_every_overlapping_record() -> None:
    """The screen's index reaches a record starting long before the call, and none of
    another chromosome or far away: one common test for the run and the calibration."""
    population = [
        PopulationSv("a", "chr17", 1400, 2000, 0.02),
        PopulationSv("b", "chr17", 1600, 2000, 0.5),
        PopulationSv("c", "chr2", 1000, 2000, 0.5),
        PopulationSv("d", "chr17", 9000, 30000, 0.5),
    ]
    s = _screen(population=population)
    kept = _deletion(5000, 6000)
    deletions = [
        _deletion(1000, 2000),
        _deletion(1700, 2000),
        _deletion(1000, 2000, chrom="chr2"),
        kept,
        _deletion(10000, 29000),
    ]
    assert rare(deletions, s) == [kept]
    match = population_match(_deletion(10000, 29000), s)
    assert match is not None and match.id == "d"


def test_two_genes() -> None:
    shown, reason = screen(_deletion(2500, 5500), _screen())
    assert reason is None
    assert list(shown) == ["GENEA", "GENEB"]
    assert shown == {"GENEA": (3,), "GENEB": (1,)}


def test_screen_all_keeps_best_at_one_position() -> None:
    s = _screen()
    dysgu_long = _deletion(
        1500, 4500, callers={"dysgu": "3 QUAL 20.0"}, depth_ratio=None
    )
    both_short = _deletion(1500, 2700)
    intronic = _deletion(1300, 1800)
    tiny = _deletion(1020, 1040)
    shown, dropped = screen_all([dysgu_long, both_short, intronic, tiny], s)
    assert [x.deletion for x in shown] == [both_short]
    assert shown[0].exons == {"GENEA": (2,)}
    assert shown[0].call_class == classify(both_short, s)
    assert shown[0].population is None
    assert dropped[SAME_POSITION] == 1
    assert dropped[NO_EXON] == 1
    assert dropped[SHORT] == 1
    assert len(shown) + sum(dropped.values()) == 4


# --- any SV type (slice 3d) -----------------------------------------------------------------

INS_T = InsertionThresholds()
EXON = {"GENEA": [Target("chr17", 995, 1010, (1,))]}


def _insertion(
    pos: int, svlen: int | None = 300, *, lower_bound: bool = False
) -> Deletion:
    return Deletion(
        "chr17", pos, pos + 1, "A", "0/1", {"dysgu": "1 QUAL 9"}, None,
        svtype="INS", svlen=svlen, lower_bound=lower_bound,
    )  # fmt: skip


def _ins_screen(population: Sequence[PopulationSv] = ()) -> Screen:
    return Screen(
        targets=EXON,
        mask={},
        repeats={},
        population=population,
        thresholds=INS_T,
        callers=("dysgu", "insurveyor"),
    )


def test_insertion_targets_point() -> None:
    """An insertion overlaps a target at its point `[pos, pos + 1)`; its length is its
    `svlen`, not `end - pos`."""
    inside = _insertion(1000)
    assert (inside.end, inside.length) == (1001, 300)
    assert screen(inside, _ins_screen()) == ({"GENEA": (1,)}, None)
    assert screen(_insertion(1100), _ins_screen()) == ({}, NO_EXON)
    assert Deletion("chr17", 1000, 1500, "A", "0/1", BOTH, None).length == 500


def test_insertion_lower_bound_by_position() -> None:
    """A lower-bound length matches gnomAD-SV by position and type alone (C20); an exact
    one needs the size ratio as well (400 / 6000 = 0.07 here)."""
    big = PopulationSv(
        "gnomad_ins", "chr17", 1300, 1301, 0.05, svtype="INS", svlen=6000
    )
    bound = _insertion(1000, 400, lower_bound=True)
    assert screen(bound, _ins_screen([big])) == ({}, COMMON)
    exact = _insertion(1000, 400)
    assert screen(exact, _ins_screen([big])) == ({"GENEA": (1,)}, None)
    # Beyond `max_distance` (500 bp) nothing matches, lower bound or not.
    far = PopulationSv("far", "chr17", 1501, 1502, 0.05, svtype="INS", svlen=400)
    assert screen(bound, _ins_screen([far])) == ({"GENEA": (1,)}, None)
    # A record of another type never matches.
    deletion = PopulationSv("del", "chr17", 900, 1100, 0.05)
    assert screen(bound, _ins_screen([deletion])) == ({"GENEA": (1,)}, None)
    assert rare([bound, exact], _ins_screen([big])) == [exact]


def test_insertion_spans_padded() -> None:
    call = _insertion(10_000)
    deletion = _deletion(1000, 2000)
    assert population_spans([call], pad=500) == [("chr17", 9500, 10501)]
    assert population_spans([deletion], pad=0) == [("chr17", 1000, 2000)]
    assert population_spans([_insertion(100)], pad=500) == [("chr17", 0, 601)]


def test_warnings_per_callers() -> None:
    """The depth lines are 3b's DELLY ratio: no other pair of callers says them."""
    pair = ("dysgu", "insurveyor")
    assert warnings(CallClass("short", "both", False, False), pair) == []
    assert warnings(CallClass("short", "insurveyor", False, False), pair) == [
        "found by one caller (INSurVeyor)"
    ]
    depth = ("dellycnv", "cnvpytor")
    assert warnings(CallClass("loss_1-5kb", "dellycnv", False, False), depth) == [
        "found by one caller (DELLY cnv)"
    ]
    assert warnings(CallClass("loss_1-5kb", "cnvpytor", False, True), depth) == [
        "found by one caller (CNVpytor)",
        "an end in a segmental duplication or low-mappability region",
    ]
    assert warnings(CallClass("long", "both", False, False)) == [NO_DEPTH]
    assert warnings(CallClass("long", "delly", False, False), ("dysgu", "delly")) == [
        "found by one caller (DELLY)",
        "read depth does not support it",
    ]


def test_classify_generic_callers() -> None:
    s = _ins_screen()
    cases = {
        "both": {"dysgu": "1 QUAL 9", "insurveyor": "A_INS_1 QUAL ."},
        "dysgu": {"dysgu": "1 QUAL 9"},
        "insurveyor": {"insurveyor": "A_INS_1 QUAL ."},
    }
    for evidence, callers in cases.items():
        d = Deletion(
            "chr17", 1000, 1001, "A", "0/1", callers, 0.1, svtype="INS", svlen=300
        )
        c = classify(d, s)
        assert c.evidence == evidence
        assert c.junction is False
    delly = Deletion("chr17", 1000, 1001, "A", "0/1", {"delly": "D QUAL 9"}, 0.1)
    with pytest.raises(ValueError, match="unknown callers"):
        classify(delly, s)
