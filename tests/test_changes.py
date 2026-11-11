"""`genokeeper changes` and the refresh log."""

import json
from collections.abc import Sequence
from dataclasses import replace
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest
from fixtures import snapshots_mini as sm
from fixtures.snapshots_mini import call, file, run_id, started, statement

from genokeeper import changes, check
from genokeeper.changes import NO_BASE, report
from genokeeper.config import DEMO_MARKER, Paths
from genokeeper.query import Filters
from genokeeper.refresh_log import (
    LOG,
    RefreshRecord,
    SourceOutcome,
    SourceStatus,
    append,
    default_base,
    read,
)

TOOLS = {"clinvar": "20261011", "analyser": "a" * 64}
INPUTS = {"clinvar38": file("cv1")}
D1, D2, D3 = run_id(1), run_id(2), run_id(3)
T0 = datetime(2026, 10, 1, 11, 0, tzinfo=UTC)


def vid(pos: int) -> str:
    return f"GRCh38:chr17:{pos}:C:T"


def at(sid: str, pos: int, **kw: Any) -> dict[str, Any]:
    return statement(sid, vid=vid(pos), **kw)


def write(
    reports: Path, rid: str, rows: list[dict[str, Any]], calls: Sequence[int] = ()
) -> None:
    sm.write(
        reports / "snapshots",
        rid,
        started(rid),
        TOOLS,
        INPUTS,
        rows,
        [call(vid(p)) for p in calls],
    )


def record(
    run: str | None,
    outcome: str = "complete",
    sources: tuple[SourceOutcome, ...] = (),
    error: str | None = None,
) -> RefreshRecord:
    return RefreshRecord(
        started=T0,
        ended=T0.replace(hour=12),
        outcome=outcome,  # pyright: ignore[reportArgumentType]
        sources=sources,
        base_run=None,
        pre_run=None,
        run=run,
        error=error,
    )


def changes_of(reports: Path, lab: Path, **kw: Any) -> dict[str, Any]:
    args: dict[str, Any] = {
        "base": None,
        "run": None,
        "f": Filters(),
        "every": True,
        "gene": None,
        "marker": None,
        "text": None,
    }
    return report(reports, lab, **{**args, **kw})


def test_single_snapshot_says_no_base(tmp_path: Path) -> None:
    write(tmp_path, D1, [at("clinvar:1", 1001)], [1001])
    r = changes_of(tmp_path, tmp_path / "lab")
    assert r["note"] == NO_BASE and r["changes"] == []
    assert (r["base"], r["run"]) == (None, D1)


def test_no_snapshot_at_all_says_no_base(tmp_path: Path) -> None:
    append(tmp_path / "lab", record(None, "failed", error="RunRefused"))
    r = changes_of(tmp_path, tmp_path / "lab")
    assert (r["base"], r["run"], r["refresh"], r["note"]) == (None, None, None, NO_BASE)


def test_default_base_is_the_previous_refresh_not_a_plain_run(tmp_path: Path) -> None:
    lab, root = tmp_path / "lab", tmp_path / "snapshots"
    for rid in (D1, D2, D3):
        write(tmp_path, rid, [at("clinvar:1", 1001)], [1001])
    assert default_base(lab, root, D3) == D2  # no refresh yet: the latest snapshot
    append(lab, record(D1))
    append(lab, record(D2, "failed", error="RunRefused"))  # D2 is no refresh's run
    append(lab, record(D3, "partial"))
    assert default_base(lab, root, D3) == D1
    assert default_base(lab, root, D1) is None
    assert changes_of(tmp_path, lab)["base"] == D1
    with pytest.raises(ValueError, match="not a usable snapshot"):
        default_base(lab, root, run_id(9))


def test_failed_source_is_listed_with_its_before(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    write(tmp_path, D1, [at("clinvar:1", 1001)], [1001])
    write(tmp_path, D2, [at("clinvar:1", 1001)], [1001])
    sources = (
        SourceOutcome("clinvar", "20261011", "20261018", SourceStatus.NEW),
        SourceOutcome(
            "cpic", "20261001", None, SourceStatus.FAILED, "CalledProcessError"
        ),
    )
    append(lab, record(D1))
    append(lab, record(D2, "partial", sources))
    r = changes_of(tmp_path, lab)
    assert r["failed"] == [{"name": "cpic", "before": "20261001"}]
    assert r["refresh"]["outcome"] == "partial" and r["base"] == D1


def test_research_associations_are_counted_unless_named(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    gwas = {
        "gene": "PPARG",
        "rsid": "rs1",
        "source": "gwas_catalog",
        "topic": "research",
    }
    write(
        tmp_path,
        D1,
        [at("gwas:rs1-T:E", 1001, validity="limited", **gwas), at("clinvar:2", 1002)],
        [1001, 1002],
    )
    write(
        tmp_path,
        D2,
        [at("gwas:rs1-T:E", 1001, **gwas), at("clinvar:2", 1002, reason="x")],
        [1001, 1002],
    )
    r = changes_of(tmp_path, lab)
    assert [c["id"] for c in r["changes"]] == ["clinvar:2"]
    assert r["research"] == {"level": 1}
    for name in ({"gene": "pparg"}, {"marker": "rs1"}, {"text": "PPARG"}):
        r = changes_of(tmp_path, lab, **name)
        assert [(c["id"], c["kind"]) for c in r["changes"]] == [
            ("gwas:rs1-T:E", "level")
        ], name
        assert r["research"] == {}
    # The first tier alone holds neither: a research association never enters or leaves.
    assert changes_of(tmp_path, lab, every=False, gene="PPARG")["changes"] == []


def test_level_3_never_enters_change_feed(tmp_path: Path) -> None:
    """Two snapshots whose only difference is a new level-3 annotation: no row, one
    research association counted; named, it is listed as new."""
    lab = tmp_path / "lab"
    research = {
        "gene": "RGS4",
        "rsid": "rs951439",
        "validity": "limited",
        "source": "pharmgkb_research",
        "topic": "research",
    }
    write(tmp_path, D1, [at("clinvar:1", 1001)], [1001])
    write(
        tmp_path,
        D2,
        [at("clinvar:1", 1001), at("pharmgkb_research:655384602", 1002, **research)],
        [1001, 1002],
    )
    r = changes_of(tmp_path, lab)
    assert (r["changes"], r["research"]) == ([], {"new": 1})
    r = changes_of(tmp_path, lab, gene="RGS4")
    assert [(c["id"], c["kind"]) for c in r["changes"]] == [
        ("pharmgkb_research:655384602", "new")
    ]


def test_research_scores_are_counted_unless_named(tmp_path: Path) -> None:
    """A new research score, like any research statement: no row, one counted; named by
    its condition or its id, listed as new."""
    lab = tmp_path / "lab"
    research = {
        "gene": None,
        "validity": "limited",
        "source": "pgsc_calc_research",
        "topic": "research",
        "claim": {"type": "pgs_research", "condition": "ADHD"},
    }
    sid = "pgsc_calc:research:PGS002746"
    write(tmp_path, D1, [at("clinvar:1", 1001)], [1001])
    write(tmp_path, D2, [at("clinvar:1", 1001), at(sid, 1002, **research)], [1001])
    r = changes_of(tmp_path, lab)
    assert (r["changes"], r["research"]) == ([], {"new": 1})
    for name in ({"text": "adhd"}, {"text": "PGS002746"}):
        r = changes_of(tmp_path, lab, **name)
        assert [(c["id"], c["kind"]) for c in r["changes"]] == [(sid, "new")], name


def test_publications_first_list_is_a_baseline(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    for rid in (D1, D2):
        write(tmp_path, rid, [at("clinvar:1", 1001)], [1001])
    sm.refreshed(tmp_path, lab, D1, {})
    sm.refreshed(
        tmp_path, lab, D2, {"clinvar:1": {"litvar": ["11"], "clinvar": ["12"]}}
    )
    assert changes_of(tmp_path, lab)["changes"] == []


def test_publications_change_from_new_pmids(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    for rid in (D1, D2):
        write(
            tmp_path, rid, [at("clinvar:1", 1001), at("clinvar:2", 1002)], [1001, 1002]
        )
    sm.refreshed(tmp_path, lab, D1, {"clinvar:1": {"litvar": ["1"], "clinvar": ["5"]}})
    sm.refreshed(
        tmp_path,
        lab,
        D2,
        {
            "clinvar:1": {"litvar": ["1", "2"], "clinvar": ["5", "6"]},
            "clinvar:2": {"litvar": ["7"], "clinvar": []},  # entering the set: baseline
        },
    )
    # Publications are in the second tier: the first is entered and left only.
    assert changes_of(tmp_path, lab, every=False)["changes"] == []
    (row,) = changes_of(tmp_path, lab)["changes"]
    assert (row["id"], row["kind"]) == ("clinvar:1", "publications")
    assert row["publications"] == {"new": ["2"], "attached_by_lab": ["6"]}


def test_first_tier_keeps_entered_and_left(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    write(
        tmp_path,
        D1,
        [
            at("clinvar:1", 1001, validity="limited"),
            at("clinvar:2", 1002),
            at("clinvar:3", 1003),
        ],
        [1001, 1002, 1003],
    )
    write(
        tmp_path,
        D2,
        [
            at("clinvar:1", 1001),
            at("clinvar:2", 1002, validity="conflicting"),
            at("clinvar:3", 1003, validity="moderate"),
        ],
        [1001, 1002, 1003],
    )
    first = changes_of(tmp_path, lab, every=False)["changes"]
    assert [(c["id"], c["kind"]) for c in first] == [
        ("clinvar:1", "entered"),
        ("clinvar:2", "left"),
    ]
    every = changes_of(tmp_path, lab)["changes"]
    assert [c["kind"] for c in every] == ["entered", "left", "level"]
    assert every[2]["dates"] == {
        "first_found": "2026-10-01",
        "last_changed": "2026-10-02",
        "assessment_changed": "2026-10-02",
        "on_or_before": True,
    }


def test_refresh_log_round_trip_and_has_no_free_text(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    first = record(
        D1,
        "partial",
        (
            SourceOutcome("clinvar", "20261011", "20261018", SourceStatus.NEW),
            SourceOutcome("cpic", "20261001", None, SourceStatus.FAILED, "http:404"),
        ),
    )
    second = record(None, "failed", error="RunRefused")
    append(lab, first)
    append(lab, second)
    assert read(lab) == [first, second]
    log = lab / "refresh" / LOG
    before = log.read_text()
    for bad in (
        record(None, "failed", error="QC refused: BRCA1 chr17:43045712 depth 3"),
        record(
            D1,
            sources=(SourceOutcome("cpic", "1", None, SourceStatus.FAILED, "a b"),),
        ),
    ):
        with pytest.raises(ValueError, match="class or a code"):
            append(lab, bad)
    assert log.read_text() == before
    keys = {key for line in before.splitlines() for key in json.loads(line)}
    assert keys == {
        "started",
        "ended",
        "outcome",
        "sources",
        "base_run",
        "pre_run",
        "run",
        "error",
    }


def test_the_check_reads_the_last_refresh_that_refreshed(tmp_path: Path) -> None:
    """A write cut part-way (power loss, a full disk) leaves half a line: the records
    before it still read. A refresh in which every source failed refreshed nothing."""
    lab, done = tmp_path / "lab", record(D1)
    failed = (SourceOutcome("clinvar", "20261011", None, SourceStatus.FAILED, "x"),)
    offline = replace(record(D1, "partial", failed), ended=T0.replace(day=20))
    for written in (done, offline):
        append(lab, written)
    with (lab / "refresh" / LOG).open("a") as log:
        log.write('{"started": "2026-10-')
    assert read(lab) == [done, offline]
    assert check.last_refresh(lab) == done.ended


def test_changes_command_prints_the_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from genokeeper.cli import main

    write(tmp_path, D1, [at("clinvar:1", 1001, validity="limited")], [1001])
    write(tmp_path, D2, [at("clinvar:1", 1001)], [1001])
    paths = Paths(
        repo=tmp_path, raw=tmp_path / "raw", reports=tmp_path, lab=tmp_path / "lab"
    )
    monkeypatch.setattr(changes, "load_paths", lambda: paths)
    assert main(["changes", "--json"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert [(c["id"], c["kind"]) for c in printed["changes"]] == [
        ("clinvar:1", "entered")
    ]
    assert main(["changes"]) == 0
    assert "entered" in capsys.readouterr().out
    assert main(["changes", "--base", "nope"]) == 2
    assert "not a usable snapshot" in capsys.readouterr().err


def write_check(lab: Path, ran_at: datetime, *sources: dict[str, Any]) -> None:
    (lab / "refresh").mkdir(parents=True, exist_ok=True)
    (lab / "refresh" / "check.json").write_text(
        json.dumps({"ran_at": ran_at.isoformat(), "notified": {}, "sources": sources})
    )


def found(name: str, status: str, upstream: str | None = "20261101") -> dict[str, Any]:
    return {
        "name": name,
        "kind": "refreshed",
        "upstream": upstream,
        "local": "20261001",
        "status": status,
        "error": None,
    }


@pytest.mark.parametrize(
    ("age", "stale"), [(0, False), (10, False), (35, False), (36, True), (400, True)]
)
def test_stale_check_is_reported(tmp_path: Path, age: int, stale: bool) -> None:
    """Over 35 days the check is silenced: `stale`, with its age in days (the tab computes
    its own at view time, from `ran_at`)."""
    lab, now = tmp_path / "lab", datetime(2026, 11, 20, 9, 0, tzinfo=UTC)
    ran = now - timedelta(days=age, hours=2)
    write_check(
        lab,
        ran,
        found("clinvar", "newer"),
        found("cpic", "same"),
        found("gwas", "failed", None),
    )
    block = changes_of(tmp_path, lab, now=now)["check"]
    assert block["ran_at"] == ran.isoformat()
    assert (block["days"], block["stale"]) == (age, stale)
    assert [n["name"] for n in block["newer"]] == ["clinvar"]
    assert block["newer"][0] == {
        "name": "clinvar",
        "kind": "refreshed",
        "local": "20261001",
        "upstream": "20261101",
    }
    assert block["failed"] == ["gwas"]  # could not check is never "up to date"


def test_missing_check_json_reads_never_ran(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    never = {"ran_at": None, "newer": [], "failed": [], "days": None, "stale": True}
    assert changes_of(tmp_path, lab)["check"] == never
    # A file that is not a check reads the same, not as a check with nothing newer.
    (lab / "refresh").mkdir(parents=True)
    (lab / "refresh" / "check.json").write_text('{"ran_at": 5')
    assert changes_of(tmp_path, lab)["check"] == never
    # So does a check written for no run: the note and the check are independent.
    write(tmp_path, D1, [at("clinvar:1", 1001)], [1001])
    assert changes_of(tmp_path, lab)["check"] == never


def test_changes_report_demo_lab(tmp_path: Path) -> None:
    """`changes --json` in a demo lab (it holds `DEMO_MARKER`) has no `check`; in any other
    lab without a check it reads "never ran"."""
    lab = tmp_path / "lab"
    lab.mkdir()
    (lab / DEMO_MARKER).write_text("{}")
    write(tmp_path, D1, [at("clinvar:1", 1001, validity="limited")], [1001])
    write(tmp_path, D2, [at("clinvar:1", 1001)], [1001])
    document = changes_of(tmp_path, lab)
    assert "check" not in document
    assert [c["id"] for c in document["changes"]] == ["clinvar:1"]
    (lab / DEMO_MARKER).unlink()
    never = {"ran_at": None, "newer": [], "failed": [], "days": None, "stale": True}
    assert changes_of(tmp_path, lab)["check"] == never


def test_check_from_the_future_is_not_negative(tmp_path: Path) -> None:
    lab, now = tmp_path / "lab", datetime(2026, 11, 20, tzinfo=UTC)
    write_check(lab, now + timedelta(days=3))
    block = changes_of(tmp_path, lab, now=now)["check"]
    assert (block["days"], block["stale"]) == (0, False)


def test_changes_text_names_reference_changed() -> None:
    row = {
        "kind": "changed",
        "gene": "CBS",
        "marker": {"rsid": None, "grch38": "chr21:1000"},
        "id": "expansionhunter:CBS",
        "before": {"validity": "moderate"},
        "after": {"validity": "moderate"},
        "causes": ["reference", "analyser"],
        "publications": {"new": [], "attached_by_lab": []},
    }
    document = {
        "base": D1,
        "run": D2,
        "refresh": None,
        "failed": [],
        "note": None,
        "changes": [row],
        "research": {},
    }
    text = changes._text(document)  # pyright: ignore[reportPrivateUsage]
    assert "causes: reference changed, analyser" in text
