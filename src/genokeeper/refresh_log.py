"""The refresh log: one JSON line per refresh in `<lab>/refresh/log.jsonl` (spec,
"Failure").

The lab is outside `reports/`, so a record holds source names, versions, run ids, outcomes,
timestamps and error classes only: never a gene, a variant or an error's text, which can
carry tool output (that goes to `reports/refresh/<started>.error.txt`). `append` refuses an
error that is not a class or a code.
"""

import json
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal, cast

from genokeeper.history import snapshots

LOG = "log.jsonl"
# An error class or code (`RunRefused`, `subprocess.CalledProcessError`, `http:404`): a
# word, never a sentence.
ERROR_CLASS = re.compile(r"[A-Za-z][A-Za-z0-9_.:-]{0,79}")
_OUTCOMES = ("complete", "partial", "failed")


class SourceStatus(StrEnum):
    """What a refresh found for one source."""

    NEW = "new"
    UNCHANGED = "unchanged"
    FAILED = "failed"


@dataclass(frozen=True)
class SourceOutcome:
    """One source in a refresh: its release before and after, and its error's class."""

    name: str
    before: str | None
    after: str | None
    status: SourceStatus
    error: str | None = None  # an error class, never text


@dataclass(frozen=True)
class RefreshRecord:
    """One refresh: when, how it ended, each source, and the runs it compared and made."""

    started: datetime
    ended: datetime
    outcome: Literal["complete", "partial", "failed"]
    sources: tuple[SourceOutcome, ...]
    base_run: str | None
    pre_run: str | None
    run: str | None
    error: str | None = None


def stamp(started: datetime) -> str:
    """A refresh's name in `reports/refresh/`: its start, UTC, to the second."""
    return f"{started.astimezone(UTC):%Y%m%dT%H%M%SZ}"


def publications_path(reports: Path, started: datetime) -> Path:
    """`reports/refresh/<started>/publications.json`, the lists of that refresh."""
    return reports / "refresh" / stamp(started) / "publications.json"


def as_json(record: RefreshRecord) -> dict[str, Any]:
    """A record as its log line's JSON object."""
    return {
        "started": record.started.isoformat(),
        "ended": record.ended.isoformat(),
        "outcome": record.outcome,
        "sources": [
            {
                "name": s.name,
                "before": s.before,
                "after": s.after,
                "status": str(s.status),
                "error": s.error,
            }
            for s in record.sources
        ],
        "base_run": record.base_run,
        "pre_run": record.pre_run,
        "run": record.run,
        "error": record.error,
    }


def _record(document: dict[str, Any]) -> RefreshRecord:
    if document["outcome"] not in _OUTCOMES:
        raise ValueError(f"unknown refresh outcome {document['outcome']!r}")
    return RefreshRecord(
        started=datetime.fromisoformat(document["started"]),
        ended=datetime.fromisoformat(document["ended"]),
        outcome=cast(Literal["complete", "partial", "failed"], document["outcome"]),
        sources=tuple(
            SourceOutcome(
                name=s["name"],
                before=s["before"],
                after=s["after"],
                status=SourceStatus(s["status"]),
                error=s["error"],
            )
            for s in document["sources"]
        ),
        base_run=document["base_run"],
        pre_run=document["pre_run"],
        run=document["run"],
        error=document["error"],
    )


def append(lab: Path, record: RefreshRecord) -> None:
    """Add `record` to `<lab>/refresh/log.jsonl`; refuse (`ValueError`) an error that is not
    a class or a code."""
    for error in (record.error, *(s.error for s in record.sources)):
        if error is not None and not ERROR_CLASS.fullmatch(error):
            raise ValueError(f"a refresh log error is a class or a code, not {error!r}")
    folder = lab / "refresh"
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = folder / LOG
    logged = path.read_text(encoding="utf-8") if path.is_file() else ""
    torn = "\n" if logged and not logged.endswith("\n") else ""  # end a cut line first
    with path.open("a", encoding="utf-8") as out:
        out.write(torn + json.dumps(as_json(record), sort_keys=True) + "\n")


def read(lab: Path) -> list[RefreshRecord]:
    """Every record of the log, oldest first; none without a log. Only this code writes
    the log, so a line that is not a record is a write cut part-way, and is skipped."""
    path = lab / "refresh" / LOG
    if not path.is_file():
        return []
    records: list[RefreshRecord] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            records.append(_record(json.loads(line)))
        except ValueError, KeyError, TypeError:
            continue
    return records


def default_base(lab: Path, root: Path, run: str) -> str | None:
    """The comparison base for the snapshot `run`: the base the latest complete or partial
    refresh of `run` recorded; for a run no refresh made, the run of the latest complete or
    partial refresh older than it, else the latest usable snapshot older than it, else
    None. A plain `genokeeper run` is never the base while a refresh is."""
    order = {s.run_id: i for i, s in enumerate(snapshots(root))}
    if run not in order:
        raise ValueError(f"{run} is not a usable snapshot")
    done = [r for r in read(lab) if r.outcome in ("complete", "partial")]
    recorded = [r.base_run for r in done if r.run == run and r.base_run in order]
    if recorded:
        return recorded[-1]
    older = [r.run for r in done if r.run in order and order[r.run] < order[run]]
    if older:
        return max(older, key=order.__getitem__)
    earlier = [rid for rid, i in order.items() if i < order[run]]
    return max(earlier, key=order.__getitem__) if earlier else None
