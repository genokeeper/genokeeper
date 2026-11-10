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
    assert any(line.startswith("  0. no run yet") for line in lines)


def test_plan_forecasts_from_past_refreshes_and_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(refresh_plan, "_missing", lambda _paths: [])
    paths = _paths(tmp_path)
    start = datetime(2026, 11, 1, tzinfo=UTC)
    for minutes, run in ((40, "r1"), (50, "r2"), (60, "r3"), (5, None)):
        refresh_log.append(
            paths.lab,
            RefreshRecord(
                start,
                start + timedelta(minutes=minutes),
                "complete",
                (),
                None,
                None,
                run,
            ),
        )
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
