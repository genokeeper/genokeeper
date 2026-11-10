"""`genokeeper refresh --plan`: what the next `genokeeper refresh` would do, how long it
would take on this machine, and how complete the lab is. It fetches nothing, runs nothing
and takes no lock; comparing the inputs with the last run may fill the lab's digest cache
(`checksum.sha256_cached`), which the first time hashes the reference.

The forecast rests on what this machine measured: the refresh log's past refreshes
(`refresh_log`) and the step and run durations in `<work>/timings.jsonl` (`timings`). A
quantity with no measurement says so, never a guess. Which steps a run rebuilds after the
code changed is not known before it runs, so that case gives the full measured cost as the
upper bound.
"""

import statistics
import subprocess
from collections.abc import Sequence
from datetime import UTC, datetime

from genokeeper import refresh_log, timings
from genokeeper.calibrate.tables import TABLES, stale_tables
from genokeeper.check import read_state, summary
from genokeeper.config import Paths
from genokeeper.identity import analyser_identity
from genokeeper.refresh import FOLDERS, REFRESHED, _dated, _latest, run_needed

# Steps a run reuses unless the reads, the reference or a tool change: what a code update
# may rebuild, longest first in the output.
HEAVY = 5


def _duration(seconds: float) -> str:
    minutes = round(seconds / 60)
    return f"{minutes} min" if minutes < 90 else f"{seconds / 3600:.1f} h"


def _median(values: Sequence[float]) -> str:
    return f"about {_duration(statistics.median(values))} (median of {len(values)})"


def plan(paths: Paths) -> list[str]:
    """The plan's lines, for the terminal."""
    lab, work = paths.lab, paths.reports / "work"
    latest = _latest(paths.reports / "snapshots")
    code_changed = latest is None or latest.tools.get("analyser") != analyser_identity()
    newer = {s["name"]: s["upstream"] for s in summary(read_state(lab))["newer"]}
    lines = ["The next `genokeeper refresh`, in order:"]
    pre_run = "no"
    if latest is None:
        pre_run = "yes, a full run from the FASTQ: there is no run yet"
    elif code_changed:
        pre_run = "yes, the code changed since the last run"
    lines.append(f"  1. a run on the releases already on disk: {pre_run}")
    lines.append(f"  2. fetch the newest release of {len(REFRESHED)} sources:")
    for name in REFRESHED:
        dated = _dated(lab / FOLDERS[name])
        here = dated[-1].name if dated else "none"
        upstream = f", newer upstream: {newer[name]}" if name in newer else ""
        lines.append(f"     {name}: on disk {here}{upstream}")
    # After a pre-run, which records the inputs on disk, only a new release makes one due.
    due = not code_changed and run_needed(lab, latest, analyser_identity())
    lines.append(
        "  3. a run on the new releases: "
        + (
            "yes, an input on disk changed since the last run"
            if due
            else "only if a fetch brings a new release"
        )
    )
    lines.append("  4. the publications, the dashboard, and old dated folders pruned")

    lines.append("")
    lines.append("Forecast on this machine:")
    records = [r for r in refresh_log.read(lab) if r.outcome != "failed"]
    without = [_seconds(r) for r in records if not _made_a_run(r)]
    with_run = [_seconds(r) for r in records if _made_a_run(r)]
    lines.append(
        "  a refresh with nothing new: "
        + (_median(without) if without else "not measured yet")
    )
    lines.append(
        "  a refresh with a run: "
        + (_median(with_run) if with_run else "not measured yet")
    )
    light = [seconds for _, seconds, built in timings.runs(work) if not built]
    lines.append(
        "  a run that rebuilds no step: "
        + (_median(light) if light else "not measured yet")
    )
    measured = timings.steps(work)
    if code_changed and latest is not None:
        total = sum(measured.values())
        lines.append(
            "  the code changed: a run rebuilds the steps whose version changed, up to "
            + (
                f"{_duration(total)} if all {len(measured)} measured steps do"
                if measured
                else "a full run, not measured yet"
            )
        )
    for step, seconds in sorted(measured.items(), key=lambda kv: -kv[1])[:HEAVY]:
        lines.append(f"    {step}: {_duration(seconds)} at its last build")

    lines.append("")
    lines.append("How complete the lab is:")
    stale = stale_tables(lab)
    lines.append(
        f"  calibration tables current: {len(TABLES) - len(stale)} of {len(TABLES)}"
    )
    if stale:
        lines.append(f"    not current: {', '.join(stale)} (`genokeeper calibrate`)")
    missing = _missing(paths)
    lines.append(
        "  bootstrap: "
        + (
            "could not check (`bash scripts/bootstrap.sh --check` failed)"
            if missing is None
            else f"{len(missing)} items missing (`bash scripts/bootstrap.sh --check`)"
            if missing
            else "nothing missing"
        )
    )
    lines.append(
        "  last run: "
        + (
            "none"
            if latest is None
            else f"{latest.run_id}, {latest.started_at:%Y-%m-%d}"
        )
    )
    if records:
        last = records[-1]
        lines.append(
            f"  last refresh: {last.started:%Y-%m-%d}, {last.outcome},"
            f" took {_duration(_seconds(last))}"
        )
    return lines


def _seconds(record: refresh_log.RefreshRecord) -> float:
    return (record.ended - record.started).total_seconds()


def _made_a_run(record: refresh_log.RefreshRecord) -> bool:
    """Whether the refresh made a run: a record names the latest run even when it made
    none, so the run counts only if its id (its start, `%Y%m%dT%H%M%SZ`) falls within it."""
    made = [r for r in (record.pre_run, record.run) if r is not None]
    for run_id in made:
        try:
            began = datetime.strptime(run_id[:16], "%Y%m%dT%H%M%SZ").replace(tzinfo=UTC)
        except ValueError:
            continue
        if record.started <= began <= record.ended:
            return True
    return False


def _missing(paths: Paths) -> list[str] | None:
    """The bootstrap's `missing:` lines (`--check` reads only), None when it failed."""
    try:
        done = subprocess.run(
            ["bash", str(paths.repo / "scripts" / "bootstrap.sh"), "--check"],
            capture_output=True,
            text=True,
            check=False,
            timeout=600,
        )
    except OSError, subprocess.TimeoutExpired:
        return None
    if done.returncode != 0:
        return None
    return [line for line in done.stdout.splitlines() if line.startswith("missing:")]
