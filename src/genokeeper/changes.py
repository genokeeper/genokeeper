"""`genokeeper changes`: what changed between two run snapshots, trusted tier first.

The comparison runs from a base snapshot (by default the previous refresh's run,
`refresh_log.default_base`) to a run (by default the latest snapshot). The first tier is
what entered or left the trusted findings (what `query` keeps at the filters); `every` adds
the second tier, every other change, new publications included. Research associations are
counted per change kind, never listed, unless a gene, marker or text filter names them.
The sources the run's refresh could not fetch are listed first, each with the release its
rows date from (`before`): theirs were not checked, never "unchanged". `check` says when the
daily update check last ran, what it found newer and what it could not check, and is `stale`
when it never ran or not for more than `STALE_DAYS`; a demo lab has no check and no `check`.

Publications are kept per refresh (`refresh_log.publications_path`): a PMID in the list of
the run's refresh and not in that of the refresh before it (no later than the base) is new;
a statement's first list is its baseline.
"""

import argparse
import json
import sys
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper import check, history, refresh_log
from genokeeper.config import DEMO_MARKER, load_paths
from genokeeper.history import REFERENCE_CAUSE, Change, ChangeKind, Dates, State
from genokeeper.query import Filters, _marker_matches, _text_matches
from genokeeper.schedule import STALE_DAYS
from genokeeper.sources.research import is_research_layer
from genokeeper.store.db import Level, Validity

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

NO_BASE = "no earlier snapshot"
FIRST_TIER = frozenset({ChangeKind.ENTERED, ChangeKind.LEFT})
_KIND_ORDER = {kind: i for i, kind in enumerate(ChangeKind)}
# A side of a change: the compared fields and the source's release.
_SIDE = (*history.COMPARED, "source_version")


def _grch38(state: State) -> str | None:
    return state.vid.removeprefix("GRCh38:") if state.vid is not None else None


def _side(state: State | None) -> dict[str, Any] | None:
    return None if state is None else {f: getattr(state, f) for f in _SIDE}


def _claims(root: Path, *run_ids: str) -> dict[str, Any]:
    """Each statement's claim in the snapshots `run_ids` (a later one wins)."""
    claims: dict[str, Any] = {}
    with duckdb.connect() as con:
        for rid in run_ids:
            for sid, claim in con.execute(
                "SELECT statement_id, claim FROM read_parquet(?)",
                [str(root / rid / "statement.parquet")],
            ).fetchall():
                claims[sid] = json.loads(claim) if claim is not None else None
    return claims


def _named(
    state: State,
    claim: Any,
    gene: str | None,
    marker: str | None,
    text: str | None,
) -> bool:
    """Whether `state` passes the gene, marker and text filters as `query`'s do."""
    finding = {
        "id": state.statement_id,
        "gene": state.gene,
        "claim": claim,
        "marker": {"rsid": state.rsid, "grch38": _grch38(state)},
    }
    if gene is not None and (state.gene or "").lower() != gene.lower():
        return False
    if marker is not None and not _marker_matches(finding, marker):
        return False
    return text is None or _text_matches(finding, text)


def _lists(reports: Path, lab: Path, base: str, run: str) -> tuple[Path | None, ...]:
    """The two `publications.json` to compare: the latest one a refresh of `run` wrote,
    and the latest one an earlier refresh wrote whose run is no later than `base` (a
    refresh whose publications step failed is skipped, so its PMIDs are not lost)."""
    snaps = history.snapshots(reports / "snapshots")
    order = {s.run_id: i for i, s in enumerate(snaps)}
    lists = [
        (order[r.run or ""], path)
        for r in refresh_log.read(lab)
        if r.outcome != "failed" and r.run in order
        if (path := refresh_log.publications_path(reports, r.started)).is_file()
    ]
    mine = [k for k, (at, _) in enumerate(lists) if at == order[run]]
    if not mine:
        return None, None
    earlier = [path for at, path in lists[: mine[-1]] if at <= order[base]]
    return (earlier[-1] if earlier else None), lists[mine[-1]][1]


def _new_publications(
    reports: Path, lab: Path, base: str, run: str
) -> dict[str, dict[str, list[str]]]:
    """Per statement, the PMIDs the run's refresh found and the base's did not: LitVar's
    as `new`, ClinVar's as `attached_by_lab`. A statement the base lists none for is a
    baseline."""
    first, last = _lists(reports, lab, base, run)
    before = history.read_publications(first)
    found: dict[str, dict[str, list[str]]] = {}
    for sid, lists in history.read_publications(last).items():
        if sid not in before:
            continue
        added = {
            key: sorted(set(lists[name]) - set(before[sid][name]), key=int)
            for key, name in (("new", "litvar"), ("attached_by_lab", "clinvar"))
        }
        if any(added.values()):
            found[sid] = added
    return found


def compared(
    reports: Path, lab: Path | None, base: str, run: str
) -> tuple[dict[str, Change], dict[str, dict[str, list[str]]]]:
    """The changes from the snapshot `base` to `run` by id, a statement whose only change
    is new publications among them, and those publications (none without the `lab`
    whose log names the refreshes)."""
    root = reports / "snapshots"
    changes = {c.statement_id: c for c in history.diff(root, base, run)}
    publications = {} if lab is None else _new_publications(reports, lab, base, run)
    if missing := set(publications) - set(changes):
        before = history.states(root, base, missing)
        after = history.states(root, run, missing)
        for sid in missing:
            if sid in before and sid in after:
                changes[sid] = Change(sid, before[sid], after[sid], (), ())
    return changes, publications


def refresh_of(lab: Path, run: str | None) -> refresh_log.RefreshRecord | None:
    """The latest refresh whose run is `run`."""
    return next(
        (
            r
            for r in reversed(refresh_log.read(lab))
            if run is not None and r.run == run
        ),
        None,
    )


def failed(refresh: refresh_log.RefreshRecord | None) -> list[dict[str, Any]]:
    """The sources `refresh` could not fetch, each with the release its rows date from."""
    return [
        {"name": s.name, "before": s.before}
        for s in (refresh.sources if refresh is not None else ())
        if s.status is refresh_log.SourceStatus.FAILED
    ]


def row(
    sid: str,
    change: Change,
    kind: ChangeKind,
    dates: Mapping[str, Dates],
    publications: Mapping[str, dict[str, list[str]]],
) -> dict[str, Any]:
    """A change as `changes --json` lists it."""
    state = change.after or change.before
    assert state is not None
    found = dates.get(sid)
    return {
        "id": sid,
        "kind": str(kind),
        "gene": state.gene,
        "marker": {"rsid": state.rsid, "grch38": _grch38(state)},
        "topic": state.topic,
        "before": _side(change.before),
        "after": _side(change.after),
        "fields": list(change.fields),
        "causes": list(change.causes),
        "dates": None if found is None else found.to_json(),
        "publications": publications.get(sid, {"new": [], "attached_by_lab": []}),
    }


def order(listed: Mapping[str, Any]) -> tuple[int, str, str]:
    """The listing order of a row: by change kind, gene and id."""
    kind = _KIND_ORDER[ChangeKind(listed["kind"])]
    return kind, listed["gene"] or "", listed["id"]


def check_block(lab: Path, now: datetime) -> dict[str, Any]:
    """`check.summary` of `<lab>/refresh/check.json` and its age in whole days at `now`
    (None when it never ran); `stale` when it never ran or is over `STALE_DAYS` old."""
    block = check.summary(check.read_state(lab))
    ran = block["ran_at"]
    days = None if ran is None else max(0, (now - datetime.fromisoformat(ran)).days)
    return {**block, "days": days, "stale": days is None or days > STALE_DAYS}


def report(
    reports: Path,
    lab: Path,
    *,
    base: str | None,
    run: str | None,
    f: Filters,
    every: bool,
    gene: str | None,
    marker: str | None,
    text: str | None,
    now: datetime | None = None,
) -> dict[str, Any]:
    """The changes from `base` to `run` (None: the defaults above) as `changes --json`
    prints them. Raises `ValueError` for a base or run that is not a usable snapshot."""
    root = reports / "snapshots"
    snaps = history.snapshots(root)
    usable = {s.run_id for s in snaps}
    for given in (base, run):
        if given is not None and given not in usable:
            raise ValueError(f"{given} is not a usable snapshot under {root}")
    if run is None and snaps:
        run = snaps[-1].run_id
    if base is None and run is not None:
        base = refresh_log.default_base(lab, root, run)
    refresh = refresh_of(lab, run)
    # A demo lab (it holds `DEMO_MARKER`) has no update check to report.
    demo = (lab / DEMO_MARKER).is_file()
    document: dict[str, Any] = {
        "base": base,
        "run": run,
        "refresh": None if refresh is None else refresh_log.as_json(refresh),
        "failed": failed(refresh),
        **({} if demo else {"check": check_block(lab, now or datetime.now(UTC))}),
        "changes": [],
        "research": {},
        "note": None,
    }
    if base is None or run is None:
        document["note"] = NO_BASE
        return document

    changes, publications = compared(reports, lab, base, run)
    trusted_before = history.trusted_ids(root, base, f)
    trusted_after = history.trusted_ids(root, run, f)
    named = gene is not None or marker is not None or text is not None
    claims = _claims(root, base, run) if text is not None else {}
    dates = history.dates(root, run)
    research: dict[ChangeKind, int] = {}
    rows: list[dict[str, Any]] = []
    for sid, change in changes.items():
        state = change.after or change.before
        assert state is not None
        kind = history.kind_of(
            change,
            f,
            change.before is not None and change.before.statement_id in trusted_before,
            change.after is not None and sid in trusted_after,
        )
        if is_research_layer(state.source) and not named:
            research[kind] = research.get(kind, 0) + 1
            continue
        if named and not _named(state, claims.get(sid), gene, marker, text):
            continue
        if not every and kind not in FIRST_TIER:
            continue
        rows.append(row(sid, change, kind, dates, publications))
    document["changes"] = sorted(rows, key=order)
    document["research"] = {
        str(kind): research[kind] for kind in ChangeKind if kind in research
    }
    return document


# --- command line ---------------------------------------------------------------------------


def _cause_text(cause: str) -> str:
    return "reference changed" if cause == REFERENCE_CAUSE else cause


def _text(document: dict[str, Any]) -> str:
    lines = [f"changes from {document['base'] or '-'} to {document['run'] or '-'}"]
    refresh = document["refresh"]
    if refresh is not None:
        lines.append(f"refresh {refresh['started']}: {refresh['outcome']}")
    for source in document["failed"]:
        lines.append(
            f"NOT CHECKED: {source['name']} (its rows date from {source['before']})"
        )
    if document["note"]:
        lines.append(document["note"])
    for row in document["changes"]:
        marker = row["marker"]["rsid"] or row["marker"]["grch38"] or "(gene-level)"
        before, after = row["before"] or {}, row["after"] or {}
        lines.append(
            f"{row['kind']:<12} {row['gene'] or '(no gene)'} {marker}  {row['id']}"
            f"  validity {before.get('validity')} -> {after.get('validity')}"
            f"  causes: {', '.join(map(_cause_text, row['causes'])) or '-'}"
        )
        pubs = row["publications"]
        if pubs["new"] or pubs["attached_by_lab"]:
            lines.append(
                f"{'':12} new PMIDs: {', '.join(pubs['new']) or '-'};"
                f" attached by a lab: {', '.join(pubs['attached_by_lab']) or '-'}"
            )
    if document["research"]:
        counts = ", ".join(f"{k} {n}" for k, n in document["research"].items())
        lines.append(f"research associations (not listed): {counts}")
    return "\n".join(lines)


def _command(args: argparse.Namespace) -> int:
    paths = load_paths()
    try:
        document = report(
            paths.reports,
            paths.lab,
            base=args.base,
            run=args.run,
            f=Filters(
                min_reliability=args.min_reliability, min_validity=args.min_validity
            ),
            every=args.all,
            gene=args.gene,
            marker=args.marker,
            text=args.text,
        )
    except (ValueError, duckdb.Error) as e:
        print(f"genokeeper changes: {e}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(document, indent=2, ensure_ascii=False))
    else:
        print(_text(document))
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `changes` sub-command."""
    c = subparsers.add_parser(
        "changes", help="what changed between two run snapshots, trusted tier first"
    )
    c.add_argument("--base", metavar="RUN_ID", help="compare with this snapshot")
    c.add_argument("--run", metavar="RUN_ID", help="the snapshot (default: the latest)")
    c.add_argument(
        "--all", action="store_true", help="every change, not only the trusted tier"
    )
    c.add_argument(
        "--min-reliability", type=Level, choices=list(Level), default=Level.MEDIUM
    )
    c.add_argument(
        "--min-validity",
        type=Validity,
        choices=list(Validity),
        default=Validity.MODERATE,
    )
    c.add_argument("--gene", help="a gene symbol (case-insensitive)")
    c.add_argument("--marker", help="an rsID or a GRCh38 vid (chr1:123:A:G)")
    c.add_argument("--text", help="a substring of the gene, condition, id or a drug")
    c.add_argument("--json", action="store_true", help="print the JSON contract")
    c.set_defaults(handler=_command)
