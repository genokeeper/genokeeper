"""The class model every structural table shares (`genokeeper.calibrate.classes`): one
measure over any feature schema, a geometry per band, the table written and read whole
through 3c's table module, and a slice-3b-format deletion table refused as stale."""

import json
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import pytest

from genokeeper.calibrate import tables
from genokeeper.calibrate.classes import (
    FORMAT,
    FORMAT_KEY,
    OUTSIDE,
    OUTSIDE_BENCHMARK,
    TABLE,
    UNSEEN,
    ClassTable,
    Geometry,
    Likelihoods,
    Share,
    measure,
    outside_benchmark,
    read_class_table,
    write_class_table,
)
from genokeeper.calibrate.deletions import open_deletions
from genokeeper.calibrate.tables import Current, open_table
from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import CallClass, Deletion
from genokeeper.deletions.screen import key as class_key


@dataclass(frozen=True)
class Klass:
    """A call class of a hand-made schema: a band and two features."""

    band: str
    evidence: str
    junction: bool = False


# 4,000 targets of 250 bp, 100 kb apart: a 2 kb call's reach is 4,000 x 2,250 bp.
G = Geometry(
    targets=tuple(("chr1", n * 100_000, n * 100_000 + 250) for n in range(1, 4001)),
    region_bp=2_000_000_000,
)
EVIDENCE = {"evidence": ("a", "b")}
WITH_JUNCTION = {**EVIDENCE, "junction": ("true", "false")}
A = Klass("long", "a")


def _band(c: Klass) -> str:
    return c.band


def _library() -> tuple[
    list[tuple[Klass, int, bool]], Sequence[tuple[str | None, bool]]
]:
    """3b's arithmetic case, one library's half: 20 false and 330 true 2 kb calls; 225 of
    250 truth events found."""
    calls = [(A, 2000, False)] * 20 + [(A, 2000, True)] * 330
    truth: Sequence[tuple[str | None, bool]] = [("long", True)] * 225 + [
        ("long", False)
    ] * 25
    return calls, truth


def _ones(features: dict[str, tuple[str, ...]]) -> Likelihoods:
    one = {
        f: {v: Share(1.0, 1.0, 1.0) for v in values} for f, values in features.items()
    }
    return Likelihoods(false=one, true={f: dict(v) for f, v in one.items()})


def test_measure_two_schemas() -> None:
    """The same calls under two schemas: the same band terms, and the junction's
    likelihood ratio (false 0.5 / true 0.1) moves `p_wrong` from 3b's 0.6135 to 0.888."""
    calls, truth = _library()
    found = {}
    for features in (EVIDENCE, WITH_JUNCTION):
        found[len(features)] = measure(
            {"pcrfree": calls, "pcr": calls},
            {"pcrfree": truth, "pcr": truth},
            {"long": G},
            {"long": 0.063},
            features,
            band=_band,
        )
    (one, lk_one), (two, lk_two) = found[1], found[2]
    assert one == two
    assert set(lk_one["long"].false) == {"evidence"}
    assert set(lk_two["long"].false) == {"evidence", "junction"}
    terms = replace(
        two["long"], false_calls=replace(two["long"].false_calls, point=40 / 700)
    )
    plain = ClassTable({}, {"long": terms}, {"long": _ones(EVIDENCE)}, {}, EVIDENCE)
    found_one = plain.p_wrong(Klass("long", "a", True))
    assert found_one is not None and found_one[0] == pytest.approx(0.6135, abs=1e-3)
    likelihoods = _ones(WITH_JUNCTION)
    likelihoods.false["junction"]["true"] = Share(0.5, 0.5, 0.5)
    likelihoods.true["junction"]["true"] = Share(0.1, 0.1, 0.1)
    table = ClassTable({}, {"long": terms}, {"long": likelihoods}, {}, WITH_JUNCTION)
    found_two = table.p_wrong(Klass("long", "a", True))
    assert found_two is not None and found_two[0] == pytest.approx(0.888, abs=1e-3)
    unmoved = table.p_wrong(Klass("long", "a", False))
    assert unmoved is not None and unmoved[0] == pytest.approx(0.6135, abs=1e-3)


def test_geometry_per_band() -> None:
    """Two bands, each with its own geometry: each band's `h` is its own reach over its
    own benchmark."""
    near = Geometry((("chr1", 100_000, 100_250),), 1_000_000)
    far = Geometry((("chr1", 100_000, 100_250),), 2_000_000)
    calls = [(Klass("near", "a"), 2000, False), (Klass("far", "a"), 2000, False)]
    bands, _ = measure(
        {"pcrfree": calls},
        {"pcrfree": []},
        {"near": near, "far": far},
        {"near": 0.1, "far": 0.1},
        EVIDENCE,
        band=_band,
    )
    assert bands["near"].h == pytest.approx(2250 / 1_000_000)
    assert bands["far"].h == pytest.approx(2250 / 2_000_000)


def test_truth_band_none_counts_nowhere() -> None:
    """A truth row whose band is None (3b: under `min_length`) is in no band's recall."""
    calls, _ = _library()
    truth: Sequence[tuple[str | None, bool]] = [("long", True)] * 9 + [
        (None, False)
    ] * 5
    bands, _ = measure(
        {"pcrfree": calls}, {"pcrfree": truth}, {"long": G}, {"long": 0.1},
        EVIDENCE, band=_band,
    )  # fmt: skip
    assert bands["long"].recall.point == pytest.approx(9.5 / 10)


def _table() -> ClassTable:
    calls, truth = _library()
    bands, likelihoods = measure(
        {"pcrfree": calls}, {"pcrfree": truth}, {"long": G}, {"long": 0.1},
        WITH_JUNCTION, band=_band,
    )  # fmt: skip
    return ClassTable(
        versions={"delly": "genomelab/delly:2.7.0-a", FORMAT_KEY: str(FORMAT)},
        bands=bands,
        likelihoods=likelihoods,
        measured={"pcrfree": {"long": {"calls": 350}}},
        features=WITH_JUNCTION,
        extra={"cnvpytor": "kept", "large": {"HG002": {"loss_large": 1}}},
    )


def test_cautions_read_each_bands_own_region() -> None:
    """Any class table names a class its calibration never saw, and a call over a target
    part outside the region its band was measured in: a loss band measured on the
    benchmark has such parts, a gain band measured genome-wide none."""
    loss = CallClass("loss_5kb+", "both", False, False)
    gain = CallClass("gain_5kb+", "both", False, False, True)
    measured = {
        "hg002": {
            b: {"classes": [class_key(c)]}
            for b, c in (("loss_5kb+", loss), ("gain_5kb+", gain))
        }
    }
    outside = {"loss_5kb+": [["chr1", 1000, 1200]], "gain_5kb+": []}
    table = ClassTable({}, {}, {}, measured, {}, extra={OUTSIDE_BENCHMARK: outside})
    over = Deletion("chr1", 1100, 9000, "N", "0/1", {}, None)
    away = Deletion("chr1", 2000, 9000, "N", "0/1", {}, None)
    assert table.cautions(over, loss) == [OUTSIDE]
    assert table.cautions(over, gain) == []
    assert table.cautions(away, loss) == []
    assert table.cautions(away, CallClass("loss_5kb+", "dellycnv", False, False)) == [
        UNSEEN
    ]


def test_outside_benchmark_per_band() -> None:
    """Each band's targets minus the region it was measured in, as JSON stores them."""
    targets = {"G": [Target("chr1", 100, 300, (1,)), Target("chr2", 100, 200, (2,))]}
    regions = {
        "loss": [("chr1", 0, 200)],
        "gain": [("chr1", 0, 10**6), ("chr2", 0, 10**6)],
    }
    assert outside_benchmark(targets, regions) == {
        "loss": [["chr1", 200, 300], ["chr2", 100, 200]],
        "gain": [],
    }


def test_class_table_roundtrip(tmp_path: Path) -> None:
    """`write_class_table` then `read_class_table` gives the table back, features (as
    tuples) and extra included, in 3c's folder with 3c's `versions.json`."""
    lab = tmp_path / "lab"
    table = _table()
    write_class_table(lab, "cnv", table, report="a report")
    found = open_table(lab, "cnv", versions=table.versions)
    assert isinstance(found, Current)
    assert read_class_table(found.folder) == table
    assert (found.folder / tables.REPORT).read_text() == "a report\n"
    document = json.loads((found.folder / TABLE).read_text())
    assert document["format"] == FORMAT == 2
    assert document["features"] == {
        "evidence": ["a", "b"],
        "junction": ["true", "false"],
    }


def _format1(lab: Path, versions: dict[str, str]) -> None:
    """A deletion table as slice 3b wrote it: format 1, no features, extra or junction."""
    one = {"evidence": {"both": [0.5, 0, 1]}, "repeat": {"true": [0.5, 0, 1]}}
    document = {
        "format": 1,
        "versions": versions,
        "bands": {
            "long": {
                "false_calls": [0.1, 0.05, 0.2], "calls": 10, "libraries": 2,
                "recall": [0.9, 0.8, 0.95], "h": 0.004, "prior": 0.063,
            }
        },
        "likelihoods": {"long": {"false": one, "true": one}},
        "measured": {},
    }  # fmt: skip

    def write(staging: Path) -> None:
        (staging / TABLE).write_text(json.dumps(document))

    tables.write_table(lab, "deletions", write, versions)


def test_format1_deletion_table_unreadable(tmp_path: Path) -> None:
    """A slice-3b deletion table (format 1, whose classes have no `junction`, which failed
    a run with `KeyError: 'junction'`) opens through the run's opener as stale with 3c's
    reason, because the run's versions name the table format; read directly it is refused,
    never read as a `ClassTable`."""
    lab = tmp_path / "lab"
    three_b = {"delly": "genomelab/delly:2.7.0-a", "fasta": "f" * 64}
    _format1(lab, three_b)
    run = three_b | {FORMAT_KEY: str(FORMAT)}
    assert open_deletions(lab, run) == (
        None,
        f"stale: {FORMAT_KEY} changed; run `genokeeper calibrate`",
    )
    with pytest.raises(ValueError, match="format 1"):
        read_class_table(tables.folder(lab, "deletions"))
    # Recorded under the run's versions but in format 1: unreadable, never a crash.
    _format1(lab, run)
    table, reason = open_deletions(lab, run)
    assert table is None and reason is not None
    assert reason.startswith("deletions table unreadable (")
    assert reason.endswith("; run `genokeeper calibrate`")
