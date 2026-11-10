"""`genokeeper refresh --plan` and the durations it forecasts from."""

from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

from genokeeper import refresh_log, refresh_plan, timings
from genokeeper.config import Paths
from genokeeper.refresh_log import RefreshRecord


def test_timings_group_steps_under_the_run_that_built_them(tmp_path: Path) -> None:
    timings.record(tmp_path, 3600.4, step="vep")
    timings.record(tmp_path, 60, run="r1")
    timings.record(tmp_path, 40, run="r2")
    (tmp_path / timings.FILE).open("a").write('{"step": "cut')
    assert timings.runs(tmp_path) == [("r1", 60, ["vep"]), ("r2", 40, [])]
    assert timings.steps(tmp_path) == {"vep": 3600}


def _paths(tmp_path: Path) -> Paths:
    return Paths(
        repo=tmp_path,
        raw=tmp_path / "raw",
        reports=tmp_path / "reports",
        lab=tmp_path / "lab",
    )


def test_plan_without_measurements_says_so(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(refresh_plan, "_missing", lambda _paths: ["missing: x"])
    lines = refresh_plan.plan(_paths(tmp_path))
    assert "  a refresh with a run: not measured yet" in lines
    assert "  a run that rebuilds no step: not measured yet" in lines
    assert "  calibration tables current: 0 of 9" in lines
    assert "  bootstrap: 1 items missing (`bash scripts/bootstrap.sh --check`)" in lines


def test_plan_forecasts_from_past_refreshes_and_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(refresh_plan, "_missing", lambda _paths: [])
    paths = _paths(tmp_path)
    start = datetime(2026, 11, 1, tzinfo=UTC)
    for day, minutes, made in (
        (1, 40, True),
        (2, 50, True),
        (3, 60, True),
        (4, 5, False),
    ):
        began = start + timedelta(days=day)
        # As `_refresh` writes it: the record names the latest run, made in it or not.
        made_at = began + timedelta(minutes=1)
        run = f"{made_at:%Y%m%dT%H%M%SZ}-aaaaaa" if made else "20261001T000000Z-bbbbbb"
        ended = began + timedelta(minutes=minutes)
        record = RefreshRecord(began, ended, "complete", (), "base", None, run)
        refresh_log.append(paths.lab, record)
    work = paths.reports / "work"
    timings.record(work, 3 * 3600, step="vep")
    timings.record(work, 4 * 3600, run="r0")
    timings.record(work, 30 * 60, run="r1")
    lines = refresh_plan.plan(paths)
    assert "  a refresh with nothing new: about 5 min (median of 1)" in lines
    assert "  a refresh with a run: about 50 min (median of 3)" in lines
    assert "  a run that rebuilds no step: about 30 min (median of 1)" in lines
    assert "    vep: 3.0 h at its last build" in lines
    assert "  bootstrap: nothing missing" in lines


def test_plan_steps_count_one_run_after_a_code_change(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No run yet: step 1 is the full run. A code change: step 1 runs, and step 3 waits for
    a new release, as `_refresh` does."""
    monkeypatch.setattr(refresh_plan, "_missing", lambda _paths: None)
    lines = refresh_plan.plan(_paths(tmp_path))
    assert (
        "  1. a run on the releases already on disk: yes, a full run from the FASTQ: there is no run yet"
        in lines
    )
    assert (
        "  3. a run on the new releases: only if a fetch brings a new release" in lines
    )
    assert (
        "  bootstrap: could not check (`bash scripts/bootstrap.sh --check` failed)"
        in lines
    )


def test_a_timing_that_cannot_be_written_fails_nothing(tmp_path: Path) -> None:
    (tmp_path / timings.FILE).mkdir()  # a folder where the file goes
    timings.record(tmp_path, 1, step="x")
    assert timings.read(tmp_path.parent / "absent") == []


def test_a_line_cut_part_way_loses_only_itself(tmp_path: Path) -> None:
    (tmp_path / timings.FILE).write_text('{"step": "cut", "sec')
    timings.record(tmp_path, 7, step="vep")
    assert timings.steps(tmp_path) == {"vep": 7}
