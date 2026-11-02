"""The depth copy-number screen (`cnv.screen`, slice 3d): which DELLY cnv and CNVpytor events
are shown, under which genes and kinds, which are large events, named loci or 3b deletions,
and why the others are dropped. Every event is built by hand; no genome is read."""

import gzip
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import pytest

from genokeeper.calibrate.deletions import class_of
from genokeeper.cnv.screen import (
    CNV_BANDS,
    CNV_FEATURE_VALUES,
    JOINED,
    LARGE_FILTERED,
    MASKED_REF,
    OWNED,
    SHORT,
    CnvScreen,
    NamedCnv,
    classify_cnv,
    gain_kind,
    load_named,
    merge_depth,
    named_deletions,
    owned_regions,
    screen_cnv,
    segdup_exons,
)
from genokeeper.config import CnvThresholds, DeletionThresholds
from genokeeper.deletions.exons import OWNED_CYP2D, Target, load_spans
from genokeeper.deletions.screen import (
    COMMON,
    NO_EXON,
    SAME_POSITION,
    CallClass,
    Deletion,
    PopulationSv,
    Screen,
    key,
)
from genokeeper.evidence.clingen import DosageList
from genokeeper.regions import Interval

T = CnvThresholds()
DELLY = {"dellycnv": "CNV00001 CN 1"}
PYTOR = {"cnvpytor": "chr1:1-2 e1 1e-9"}
BOTH = DELLY | PYTOR
# GENEA: three target exons on chr1, its span their extent (`gene_spans`).
TARGETS = {
    "GENEA": [
        Target("chr1", 10_000, 10_500, (2,)),
        Target("chr1", 12_000, 12_500, (3,)),
        Target("chr1", 30_000, 30_400, (4,)),
    ]
}
SPANS = {"GENEA": ("chr1", 10_000, 30_400)}
NO_DOSAGE = DosageList("20261023", frozenset())


def _event(
    pos: int, end: int, callers: dict[str, str] = DELLY, svtype: str = "DEL", **kw: str
) -> Deletion:
    chrom = kw.get("chrom", "chr1")
    gt = "0/1"
    return Deletion(chrom, pos, end, "N", gt, callers, None, svtype=svtype)


def _screen(
    mask: Sequence[Interval] = (), population: Sequence[PopulationSv] = ()
) -> Screen:
    return Screen(
        targets=TARGETS,
        mask={"segdup": list(mask)},
        repeats={},
        population=list(population),
        thresholds=DeletionThresholds(),
        callers=("dellycnv", "cnvpytor"),
    )


def _run(
    events: Sequence[Deletion],
    *,
    deletions: Sequence[Deletion] = (),
    junctions: Sequence[Deletion] = (),
    screen: Screen | None = None,
    masked: Sequence[Interval] = (),
    owned: Sequence[Interval] = (),
    dosage: DosageList = NO_DOSAGE,
    regions: Sequence[tuple[str, int, int, str]] = (),
    recurrent: Sequence[tuple[str, int, int, str]] = (),
    named: Sequence[NamedCnv] = (),
) -> CnvScreen:
    result = screen_cnv(
        events,
        deletions,
        junctions,
        screen or _screen(),
        T,
        masked=masked,
        owned=owned,
        spans=SPANS,
        dosage=dosage,
        regions=regions,
        recurrent=recurrent,
        named=named,
    )
    # Every event is shown, large, named or dropped for one reason.
    out = len(result.shown) + len(result.large) + len(result.named)
    assert out + sum(result.dropped.values()) == len(events)
    return result


def _only(result: CnvScreen) -> str:
    """The one reason the screen dropped its one event."""
    (reason,) = [r for r, n in result.dropped.items() if n]
    assert result.dropped[reason] == 1
    return reason


def test_depth_loss_names_plus_clinical_exons_apart() -> None:
    """A depth loss over a Plus-Clinical-only target names its number apart from Select's,
    as a 3b deletion does."""
    plus = Target("chr1", 20_000, 20_200, (3,), plus_clinical=True)
    screen = replace(_screen(), targets={"GENEA": [*TARGETS["GENEA"], plus]})
    (shown,) = _run([_event(19_000, 22_000)], screen=screen).shown
    assert (shown.exons, shown.plus_clinical) == ({"GENEA": ()}, {"GENEA": (3,)})


def test_depth_loss_over_deletion_joins_it() -> None:
    """Review Focus 1: a depth loss over a shown 3b deletion (reciprocal overlap 0.81, its
    ends a few hundred bp off, so another exon set) is that deletion's: dropped `JOINED`,
    nothing shown, and `joined` names the depth caller for the 3b statement's reason."""
    deletion = Deletion("chr1", 10_000, 13_000, "A", "0/1", {"delly": "D1"}, 0.4)
    result = _run([_event(9_700, 13_400)], deletions=[deletion])
    assert _only(result) == JOINED
    assert result.shown == []
    assert result.joined == {("chr1", 10_000, 13_000): "DELLY cnv"}
    both = _run([_event(9_700, 13_400, BOTH)], deletions=[deletion])
    assert both.joined == {("chr1", 10_000, 13_000): "DELLY cnv and CNVpytor"}
    # A gain over the same deletion is not it.
    gain = _run([_event(9_700, 13_400, svtype="DUP")], deletions=[deletion])
    assert len(gain.shown) == 1 and gain.joined == {}


def test_masked_ref_dropped() -> None:
    """Half its length in 3c's masked intervals, from either caller: `MASKED_REF`, before
    the owned, common and exon rules; a third of it, with no end inside, is not."""
    masked = [("chr1", 11_200, 12_300)]  # 1,100 of 2,000 bp, neither end
    common = [PopulationSv("gnomad1", "chr1", 10_900, 13_100, 0.2)]
    for callers in (DELLY, PYTOR):
        result = _run(
            [_event(11_000, 13_000, callers)],
            masked=masked,
            owned=[("chr1", 10_000, 14_000)],
            screen=_screen(population=common),
        )
        assert _only(result) == MASKED_REF
    third = _run([_event(10_000, 13_000)], masked=[("chr1", 11_000, 12_000)])
    assert len(third.shown) == 1
    end = _run([_event(10_000, 13_000)], masked=[("chr1", 12_990, 13_500)])
    assert _only(end) == MASKED_REF


def test_owned_region_dropped(tmp_path: Path) -> None:
    """SMN1's span, a Parascopy model's main and partner copies and the CYP2D locus are
    owned: a loss over any is dropped `OWNED`."""
    model = tmp_path / "STRC.gz"
    with gzip.open(model, "wt") as f:
        f.write("main\tchr15:43558001-43650000\t*\tname=STRC\n")
        f.write(
            "duplication\tchr15:43558938-43580541\tchr15:43658492-43680099:+\tix=0\n"
        )
    spans = load_spans()
    owned = owned_regions(spans, [model])
    smn1 = spans["SMN1"]
    assert smn1 in owned and spans["SMN2"] in owned and OWNED_CYP2D in owned
    assert ("chr15", 43_558_000, 43_650_000) in owned
    assert ("chr15", 43_658_491, 43_680_099) in owned
    for event in (
        _event(smn1[1] + 1_000, smn1[1] + 5_000, chrom="chr5"),
        _event(43_660_000, 43_670_000, chrom="chr15"),
        _event(42_139_575, 42_144_483, chrom="chr22"),  # CYP2D7
    ):
        assert _only(_run([event], owned=owned)) == OWNED
    bad = tmp_path / "BAD.gz"
    with gzip.open(bad, "wt") as f:
        f.write("duplication\tchr15:1-2\tchr15:3-4:+\tix=0\n")
    with pytest.raises(ValueError, match="main"):
        owned_regions(spans, [bad])


def test_short_and_common() -> None:
    """Under 1 kb: `SHORT`. A gain matched by a common gnomAD-SV duplication: `COMMON`; a
    common deletion at the same place does not cut a gain."""
    assert _only(_run([_event(10_000, 10_999)])) == SHORT
    dup = PopulationSv("gnomad_dup", "chr1", 9_900, 12_600, 0.05, svtype="DUP")
    dele = PopulationSv("gnomad_del", "chr1", 9_900, 12_600, 0.05)
    gain = _event(10_000, 12_600, svtype="DUP")
    assert _only(_run([gain], screen=_screen(population=[dup]))) == COMMON
    kept = _run([gain], screen=_screen(population=[dele]))
    assert len(kept.shown) == 1


def test_merge_depth() -> None:
    """DELLY's and CNVpytor's losses at reciprocal overlap 0.7 are one event, DELLY's
    coordinates and both callers; a DELLY loss and a CNVpytor gain stay two; without
    CNVpytor DELLY's alone."""
    delly = _event(10_000, 20_000)
    pytor = _event(13_000, 21_000, PYTOR)  # 7,000 of 10,000
    (merged,) = merge_depth([delly], [pytor], min_overlap=0.5)
    assert (merged.pos, merged.end, dict(merged.callers)) == (10_000, 20_000, BOTH)
    gain = _event(13_000, 21_000, PYTOR, svtype="DUP")
    assert merge_depth([delly], [gain], min_overlap=0.5) == [delly, gain]
    assert merge_depth([delly], None, min_overlap=0.5) == [delly]
    far = _event(19_000, 30_000, PYTOR)
    assert merge_depth([delly], [far], min_overlap=0.5) == [delly, far]


def test_classify_cnv_bands() -> None:
    s = _screen()
    assert classify_cnv(_event(0, 4_999), s, T, junction=False).band == "loss_1-5kb"
    assert classify_cnv(_event(0, 5_000), s, T, junction=False).band == "loss_5kb+"
    large = classify_cnv(_event(0, 500_000, BOTH, "DUP"), s, T, junction=True)
    assert large == CallClass("gain_large", "both", False, False, True)
    assert not classify_cnv(_event(0, 2_000), s, T, junction=True).junction
    masked = _screen(mask=[("chr1", 0, 100)])
    assert classify_cnv(_event(50, 2_000, PYTOR), masked, T, junction=False) == (
        CallClass("loss_1-5kb", "cnvpytor", False, True, False)
    )


def test_cnv_class_keys_roundtrip() -> None:
    values = CNV_FEATURE_VALUES
    for band in CNV_BANDS:
        for evidence in values["evidence"]:
            for repeat in values["repeat"]:
                for masked in values["masked"]:
                    for junction in values["junction"]:
                        c = CallClass(
                            band,
                            evidence,
                            repeat == "true",
                            masked == "true",
                            junction == "true",
                        )
                        assert class_of(key(c)) == c
    assert set(CNV_BANDS) == {
        f"{d}_{s}" for d in ("loss", "gain") for s in ("1-5kb", "5kb+", "large")
    }


def test_gain_kind() -> None:
    span = ("chr1", 10_000, 30_400)
    inside = _event(11_000, 20_000, svtype="DUP")
    assert gain_kind(inside, span, junction=True) == "tandem"
    assert gain_kind(inside, span, junction=False) == "intragenic"
    assert gain_kind(_event(9_000, 31_000, svtype="DUP"), span, True) == "whole"
    assert gain_kind(_event(9_000, 20_000, svtype="DUP"), span, True) == "partial"


def test_whole_gene_gain_with_inner_junction() -> None:
    """Review Focus 2: a gain over the whole gene whose DUP junction lies inside it is a
    whole-gene gain, not a tandem duplication."""
    gain = _event(9_000, 31_000, BOTH, "DUP")
    junction = Deletion(
        "chr1", 10_200, 30_000, "N", "0/1", {"dysgu": "J1"}, None, "DUP"
    )
    (shown,) = _run([gain], junctions=[junction]).shown
    assert shown.kinds == {"GENEA": "whole"}
    assert shown.call_class.junction
    assert shown.exons == {"GENEA": (2, 3, 4)}
    tandem = _event(11_000, 20_000, BOTH, "DUP")
    inner = Deletion("chr1", 11_100, 19_900, "N", "0/1", {"dysgu": "J2"}, None, "DUP")
    (shown,) = _run([tandem], junctions=[inner]).shown
    assert shown.kinds == {"GENEA": "tandem"}
    (alone,) = _run([tandem]).shown
    assert alone.kinds == {"GENEA": "intragenic"}
    (loss,) = _run([_event(11_000, 20_000)], junctions=[inner]).shown
    assert loss.kinds == {"GENEA": "loss"} and not loss.call_class.junction


def test_large_event_filters() -> None:
    """Review Focus 6: a 600 kb loss over a ClinGen region is a large event only when both
    depth callers agree; DELLY alone: `LARGE_FILTERED`; common in gnomAD-SV: `COMMON`;
    inside 3c's masked intervals: `MASKED_REF`. 400 kb is no large event."""
    region = [("chr2", 1_000_000, 1_400_000, "2q region")]
    loss = _event(900_000, 1_500_000, BOTH, chrom="chr2")
    (large,) = _run([loss], regions=region).large
    assert large.regions == ("2q region",)
    assert large.call_class.band == "loss_large"
    alone = _event(900_000, 1_500_000, chrom="chr2")
    assert _only(_run([alone], regions=region)) == LARGE_FILTERED
    common = PopulationSv("gnomad_big", "chr2", 900_000, 1_500_000, 0.02)
    assert _only(_run([loss], regions=region, screen=_screen(population=[common]))) == (
        COMMON
    )
    masked = [("chr2", 900_000, 1_300_000)]
    assert _only(_run([loss], regions=region, masked=masked)) == MASKED_REF
    small = _event(1_000_000, 1_400_000, BOTH, chrom="chr2")
    assert _only(_run([small], regions=region)) == NO_EXON
    # A dosage-sensitive gene's span counts as a region does.
    dosage = DosageList(
        "20261023",
        frozenset({"GENEH"}),
        frozenset(),
        {"GENEH": ("chr2", 1_200_000, 1_210_000)},
    )
    (by_gene,) = _run([loss], dosage=dosage).large
    assert by_gene.regions == ("GENEH",)
    # Over target exons but no region or dosage-sensitive gene: never an ordinary event.
    over_exons = _event(5_000, 605_000, BOTH)
    assert _only(_run([over_exons])) == LARGE_FILTERED


def test_large_segdup_ends_need_recurrent_region() -> None:
    """A 2.5 Mb loss both callers make, both ends in 3b's mask (segmental duplications),
    matching a 22q11.2-like recurrent region at reciprocal overlap 0.9: a large event. A
    700 kb loss with one end in the mask over an HI-3 gene, matching no recurrent region
    (0.3 with the nearest): `LARGE_FILTERED`."""
    mask = [("chr22", 18_600_000, 18_700_000), ("chr22", 21_000_000, 21_200_000)]
    recurrent = [("chr22", 18_750_000, 21_000_000, "22q11.2 recurrent (DGS)")]
    regions = [("chr22", 18_750_000, 21_000_000, "22q11.2 recurrent (DGS)")]
    loss = _event(18_650_000, 21_150_000, BOTH, chrom="chr22")  # 2.25 / 2.5 Mb
    (large,) = _run(
        [loss], screen=_screen(mask=mask), regions=regions, recurrent=recurrent
    ).large
    assert large.call_class.masked
    dosage = DosageList(
        "20261023",
        frozenset({"GENEH"}),
        frozenset(),
        {"GENEH": ("chr22", 21_500_000, 21_510_000)},
    )
    other = _event(21_100_000, 21_800_000, BOTH, chrom="chr22")
    nearest = [("chr22", 20_900_000, 21_310_000, "nearby recurrent")]  # 210 / 700 kb
    result = _run([other], screen=_screen(mask=mask), dosage=dosage, recurrent=nearest)
    assert _only(result) == LARGE_FILTERED
    assert result.large == []


def test_named_loss_and_named_deletions() -> None:
    """A depth loss at reciprocal overlap 0.95 with GJB6-D13S1830, over no screen exon, is
    the named locus, not `NO_EXON`; a 3b deletion there is too (`named_deletions`); at 0.3,
    neither."""
    (row,) = load_named()
    assert (row.name, row.chrom, row.start, row.end, row.gene) == (
        "GJB6-D13S1830",
        "chr13",
        20_223_037,
        20_531_806,
        "GJB2",
    )
    loss = _event(20_210_000, 20_520_000, BOTH, chrom="chr13")
    (named, shown) = _run([loss], named=[row]).named[0]
    assert named == row and shown.deletion == loss
    assert shown.call_class.band == "loss_5kb+"
    far = _event(20_438_000, 20_748_000, BOTH, chrom="chr13")
    assert _only(_run([far], named=[row])) == NO_EXON
    deletion = Deletion(
        "chr13", 20_210_000, 20_520_000, "A", "0/1", {"delly": "D"}, 0.4
    )
    s = Screen(TARGETS, {}, {}, [], DeletionThresholds())
    ((named3b, shown3b),) = named_deletions([deletion], [row], s)
    assert named3b == row and shown3b.call_class == CallClass(
        "long", "delly+depth", False, False
    )
    far3b = Deletion("chr13", 20_438_000, 20_748_000, "A", "0/1", {"delly": "D"}, 0.4)
    assert named_deletions([far3b], [row], s) == []


def test_same_position_keeps_the_best() -> None:
    """Two shown losses at one position: the one with more callers is kept."""
    result = _run([_event(10_000, 12_600), _event(10_000, 12_400, BOTH)])
    (shown,) = result.shown
    assert shown.call_class.evidence == "both"
    assert result.dropped[SAME_POSITION] == 1


def test_segdup_exons() -> None:
    """A target wholly in the mask and outside every owned region is listed; one half
    outside the mask, or in an owned region, is not."""
    targets = {
        "GENEA": [
            Target("chr1", 100, 200, (1,)),
            Target("chr1", 300, 400, (2,)),
            Target("chr1", 500, 600, (3,)),
        ],
        "GENEB": [Target("chr2", 100, 200, (1,))],
    }
    mask = [("chr1", 50, 250), ("chr1", 350, 450), ("chr1", 450, 700), ("chr2", 0, 500)]
    owned = [("chr1", 550, 560), ("chr2", 0, 1_000)]
    assert segdup_exons(targets, mask, owned) == {"GENEA": (1,)}
    touching = [("chr1", 300, 350), ("chr1", 350, 400)]
    assert segdup_exons(targets, touching, []) == {"GENEA": (2,)}
