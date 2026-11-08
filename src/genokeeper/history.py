"""What changed between run snapshots: each statement's state, its changes and their
causes, and its three dates (stage 4, "Computing a change").

Every run exports `<root>/<run_id>/{run,statement,call,gene_call,variant}.parquet`
(`store.db.export_snapshot`); this module reads only those files. A snapshot is usable when
each file is there with the columns the comparison reads; one an older schema wrote is
skipped, and a date that falls on the oldest usable snapshot is "on or before" it.

A statement is the same across snapshots when its id matches, or, for ClinVar, when one
`clinvar:` statement disappears while one other appears at the same vid in one pair of
consecutive snapshots (ClinVar merged the records): it carries on under the new id. A
change is a difference in a compared field (`COMPARED`), never in the source's version
alone. A non-finding (a claim at an allele the user does not carry, or one classified
benign) is read as absent: a new ClinVar release adds such claims by the thousand, none
about the user. A finding or watch item that becomes one is left or gone, its non-finding
state kept as the change's `after`.

DuckDB compares each pair of consecutive snapshots and returns only the statements that
differ: a snapshot holds about a million statements, a pair's changes a few hundred.

Trusted is what `query.query` keeps at the given filters, run over the snapshot itself;
a research association neither enters nor leaves the trusted tier.
"""

import hashlib
import json
import re
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, date, datetime
from enum import StrEnum
from importlib.resources import files
from pathlib import Path
from typing import Any

import duckdb

from genokeeper.query import OWNS_GENE_CALL, Filters, query
from genokeeper.sources.research import is_research_layer
from genokeeper.store.db import TABLES, Kind

COMPARED: tuple[str, ...] = (
    "kind",
    "validity",
    "source_level",
    "reason",
    "effect",
    "actionable",
    "topic",
    "dosage",
    "level",
    "result",
)
# The fields whose change is a change of the assessment (validity or reliability).
_ASSESSMENT = frozenset({"validity", "level"})
# The fields that are the user's call: only they moving, under the same analyser and
# inputs, is a change of the genotype.
_GENOTYPE = frozenset({"dosage", "level", "result"})


class ChangeKind(StrEnum):
    """What a change is for the user (spec, "The Changes tab")."""

    ENTERED = "entered"
    LEFT = "left"
    LEVEL = "level"
    NEW = "new"
    GONE = "gone"
    DETAILS = "details"
    PUBLICATIONS = "publications"


@dataclass(frozen=True)
class Snapshot:
    """A usable snapshot's run: its id, start (UTC), tools and inputs."""

    run_id: str
    started_at: datetime
    tools: dict[str, Any] = field(hash=False)
    inputs: dict[str, Any] = field(hash=False)


@dataclass(frozen=True)
class State:
    """A statement as one snapshot holds it: the compared fields, and what names it.
    `effect` is canonical JSON text; `level` a gene-level statement's gene call's, and
    `result` its gene call's diplotype and phenotype, else its claim's haplogroup or
    diplotype (APOE's); null for a statement at a variant."""

    statement_id: str
    gene: str | None
    vid: str | None
    rsid: str | None
    kind: str
    validity: str | None
    source_level: str | None
    reason: str | None
    effect: str | None
    actionable: bool | None
    topic: str | None
    dosage: str | None
    level: str | None
    result: str | None
    source: str | None
    source_version: str | None


@dataclass(frozen=True)
class Change:
    """A statement that differs between two snapshots (by its id in the later one): its
    states (None where absent), the compared fields that differ (none when it appeared or
    disappeared) and the causes."""

    statement_id: str
    before: State | None
    after: State | None
    fields: tuple[str, ...]
    causes: tuple[str, ...]


@dataclass(frozen=True)
class Dates:
    """A statement's first found (first a finding), last changed and assessment changed
    dates; `on_or_before` when it was already in the oldest usable snapshot, so a date that
    falls on that snapshot reads "on or before" it."""

    first_found: date | None
    last_changed: date | None
    assessment_changed: date | None
    on_or_before: bool

    def to_json(self) -> dict[str, Any]:
        """The dates as ISO strings (null where none), as the JSON contract prints them."""
        return {
            "first_found": _iso(self.first_found),
            "last_changed": _iso(self.last_changed),
            "assessment_changed": _iso(self.assessment_changed),
            "on_or_before": self.on_or_before,
        }


def _iso(day: date | None) -> str | None:
    return None if day is None else day.isoformat()


# The `run.inputs` a source's statements are graded on (a dotted path names a nested
# input). A source not listed here (`pgsc_calc`, `haplogrep3`, stage 3's) reads every
# input as a possible cause. `grading.toml` is part of the analyser's identity.
SOURCE_INPUTS: dict[str, tuple[str, ...]] = {
    "clinvar": ("clinvar38", "gene_validity", "gencc", "cgd", "faf"),
    "pharmgkb": ("pgx.pharmgkb", "clinvar38"),
    "pharmcat": ("pgx.cpic", "pgx.positions"),
    "gwas_catalog": ("gwas",),
    "pharmgkb_research": ("pgx.pharmgkb", "clinvar38"),
    "vep": (
        "clinvar38",
        "vep",
        "clingen_dosage",
        "gene_validity",
        "gencc",
        "cgd",
        "faf",
    ),
    # The panel, and what its statements are graded on or repeat (`sources.panelapp`).
    "panelapp": (
        "panelapp",
        "clinvar38",
        "vep",
        "clingen_dosage",
        "gene_validity",
        "gencc",
        "cgd",
        "faf",
    ),
}
# Inputs that are no cause: a run's cache records, whose `reused` flag differs between two
# runs on the same inputs; the reference FASTA, whose change is `REFERENCE_CAUSE`; and slice
# 3c's checks of the run itself (the identity verdict, the coverage and dosage steps'
# failures), which grade no statement; and the PanelApp panel, which only its own layer
# reads, so a moved panel is no cause of a source this table does not list.
_NOT_INPUTS = frozenset(
    {"steps", "fasta", "identity", "coverage", "dosage", "panelapp"}
)
# A change the reference moved (the masked GRCh38's md5 in `tools`): new calls, not news.
REFERENCE_CAUSE = "reference"
# The tools key naming an input's release when its record carries none.
_RELEASE_TOOL = {"clinvar38": "clinvar", "pgx.cpic": "cpic"}
_RELEASE_KEYS = ("release", "date", "created")
_PMID = re.compile(r"[0-9]+")

# The columns each snapshot file must have for the comparison.
_REQUIRED: dict[str, frozenset[str]] = {
    "run": frozenset({"run_id", "started_at", "inputs", "tools"}),
    "statement": frozenset(
        {
            "statement_id",
            "kind",
            "vid",
            "gene",
            "claim",
            "validity",
            "source_level",
            "reason",
            "effect",
            "actionable",
            "topic",
            "source",
            "source_version",
        }
    ),
    "call": frozenset({"vid", "dosage", "level"}),
    "gene_call": frozenset({"gene", "level", "diplotype", "phenotype"}),
    "variant": frozenset({"vid", "rsid"}),
}


def _parquet(root: Path, run_id: str, table: str) -> str:
    return str(root / run_id / f"{table}.parquet")


def _usable(con: duckdb.DuckDBPyConnection, folder: Path) -> bool:
    for table, needed in _REQUIRED.items():
        path = folder / f"{table}.parquet"
        if not path.is_file():
            return False
        try:
            columns = {
                str(row[0])
                for row in con.execute(
                    "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
                ).fetchall()
            }
        except duckdb.Error:
            return False
        if not needed <= columns:
            return False
    return True


def _json(value: Any) -> Any:
    return json.loads(value) if isinstance(value, str) else (value or {})


def snapshots(root: Path) -> list[Snapshot]:
    """The usable snapshots under `root`, oldest first."""
    if not root.is_dir():
        return []
    found: list[Snapshot] = []
    with duckdb.connect() as con:
        for folder in sorted(p for p in root.iterdir() if p.is_dir()):
            if not _usable(con, folder):
                continue
            row = con.execute(
                "SELECT run_id, started_at, inputs, tools FROM read_parquet(?) LIMIT 1",
                [str(folder / "run.parquet")],
            ).fetchone()
            if row is None or row[0] != folder.name:
                continue
            found.append(
                Snapshot(
                    run_id=row[0],
                    started_at=row[1].replace(tzinfo=UTC),
                    tools=_json(row[3]),
                    inputs=_json(row[2]),
                )
            )
    return sorted(found, key=lambda s: (s.started_at, s.run_id))


def local_day(snapshot: Snapshot) -> date:
    """The user's calendar day (this machine's time zone) on which `snapshot`'s run
    started: the day every date shown is."""
    return snapshot.started_at.astimezone().date()


def _canonical(effect: str | None) -> str | None:
    return None if effect is None else json.dumps(json.loads(effect), sort_keys=True)


def _literal(path: str) -> str:
    return "'" + path.replace("'", "''") + "'"


def _select(root: Path, run_id: str, *, owner_only: bool) -> str:
    """A query for `run_id`'s statements in `State`'s columns; `owner_only` leaves out the
    non-findings."""
    path = {t: _literal(_parquet(root, run_id, t)) for t in _REQUIRED if t != "run"}
    where = f" WHERE s.kind <> '{Kind.NON_FINDING}'" if owner_only else ""
    return (
        "SELECT s.statement_id, s.gene, s.vid, v.rsid, s.kind, s.validity,"
        " s.source_level, s.reason, CAST(s.effect AS VARCHAR) AS effect, s.actionable,"
        " s.topic, c.dosage, CASE WHEN s.vid IS NULL THEN g.level ELSE c.level END"
        " AS level, CASE WHEN s.vid IS NULL THEN COALESCE("
        "NULLIF(concat_ws(' ', g.diplotype, g.phenotype), ''),"
        " json_extract_string(s.claim, '$.haplogroup'),"
        " json_extract_string(s.claim, '$.diplotype')) END AS result,"
        " s.source, s.source_version"
        f" FROM read_parquet({path['statement']}) AS s"
        f" LEFT JOIN read_parquet({path['call']}) AS c ON c.vid = s.vid"
        f" LEFT JOIN read_parquet({path['variant']}) AS v ON v.vid = s.vid"
        f" LEFT JOIN read_parquet({path['gene_call']}) AS g"
        f" ON s.vid IS NULL AND {OWNS_GENE_CALL.format(s='s')}" + where
    )


def _state(values: Sequence[Any]) -> State:
    row = list(values)
    row[8] = _canonical(row[8])
    return State(*row)


def states(
    root: Path, run_id: str, ids: Collection[str] | None = None
) -> dict[str, State]:
    """Every statement of the snapshot `run_id` (or those of `ids`), by id."""
    sql, params = _select(root, run_id, owner_only=False), []
    if ids is not None:
        sql = f"SELECT * FROM ({sql}) WHERE statement_id IN (SELECT UNNEST(?))"
        params = [sorted(ids)]
    with duckdb.connect() as con:
        rows = con.execute(sql, params).fetchall()
    return {state.statement_id: state for state in map(_state, rows)}


def _fields(before: State, after: State) -> tuple[str, ...]:
    return tuple(f for f in COMPARED if getattr(before, f) != getattr(after, f))


_NAMES = tuple(State.__dataclass_fields__)
type _Pair = tuple[State | None, State | None]


def _load(con: duckdb.DuckDBPyConnection, root: Path, run_id: str) -> None:
    """Make `run_id`'s statements about the user the table `after`, and the previous
    `after` the table `before`."""
    con.execute("DROP TABLE IF EXISTS before")
    con.execute("ALTER TABLE IF EXISTS after RENAME TO before")
    con.execute(f"CREATE TEMP TABLE after AS {_select(root, run_id, owner_only=True)}")


def _pairs(con: duckdb.DuckDBPyConnection) -> list[_Pair]:
    """The statements that differ between the tables `before` and `after` (None where
    absent), a ClinVar record merged into another paired with it whether or not a field
    changed."""
    columns = ", ".join(f"{side}.{n}" for side in ("b", "a") for n in _NAMES)
    differs = " OR ".join(f"b.{f} IS DISTINCT FROM a.{f}" for f in COMPARED)
    rows = con.execute(
        f"SELECT {columns} FROM before AS b FULL OUTER JOIN after AS a"
        " ON a.statement_id = b.statement_id"
        f" WHERE b.statement_id IS NULL OR a.statement_id IS NULL OR {differs}"
    ).fetchall()
    width = len(_NAMES)
    pairs: list[_Pair] = []
    for row in rows:
        before = _state(row[:width]) if row[0] is not None else None
        after = _state(row[width:]) if row[width] is not None else None
        # The effect differs as text only (its keys in another order): no change.
        if before is None or after is None or _fields(before, after):
            pairs.append((before, after))
    return _merge_renames(pairs)


def _merge_renames(pairs: list[_Pair]) -> list[_Pair]:
    """Pair a `clinvar:` statement gone with the one other that appeared at its vid."""
    gone: dict[str, list[State]] = {}
    new: dict[str, list[State]] = {}
    for before, after in pairs:
        for state, into in ((before, gone), (after, new)):
            if (
                state is not None
                and (before is None or after is None)
                and state.statement_id.startswith("clinvar:")
                and state.vid
            ):
                into.setdefault(state.vid, []).append(state)
    merged: dict[str, State] = {
        old[0].statement_id: new[vid][0]
        for vid, old in gone.items()
        if len(old) == 1 and len(new.get(vid, ())) == 1
    }
    taken = {state.statement_id for state in merged.values()}
    found: list[_Pair] = []
    for before, after in pairs:
        if before is not None and after is None and before.statement_id in merged:
            found.append((before, merged[before.statement_id]))
        elif not (before is None and after is not None and after.statement_id in taken):
            found.append((before, after))
    return found


def _node(inputs: Mapping[str, Any], path: str) -> Any:
    node: Any = inputs
    for part in path.split("."):
        node = node.get(part) if isinstance(node, Mapping) else None
    return node


def _digest(node: Any) -> str:
    """A file input's sha256; a nested map's, the sha256 of its members' sorted
    `name:digest` lines; anything else, the sha256 of its canonical JSON."""
    if isinstance(node, Mapping):
        sha = node.get("sha256")
        if isinstance(sha, str):
            return sha
        lines = sorted(f"{name}:{_digest(member)}" for name, member in node.items())
        text = "\n".join(lines)
    else:
        text = json.dumps(node, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest()


def _release(snapshot: Snapshot, path: str) -> str:
    """An input's release in `snapshot`: its record's, else its tool version's, else the
    first 12 hex digits of its digest."""
    node = _node(snapshot.inputs, path)
    if isinstance(node, Mapping):
        for key in _RELEASE_KEYS:
            if isinstance(node.get(key), str):
                return node[key]
    tool = snapshot.tools.get(_RELEASE_TOOL.get(path, ""))
    if isinstance(tool, str):
        return tool
    return _digest(node)[:12]


def _source_causes(a: Snapshot, b: Snapshot, source: str | None) -> tuple[str, ...]:
    """`source:<input>@<release>` for each input `source` reads that differs from `a` to
    `b`; for a source `SOURCE_INPUTS` does not list, every input that differs."""
    listed = SOURCE_INPUTS.get(source or "")
    names: Iterable[str] = (
        listed
        if listed is not None
        else sorted((set(a.inputs) | set(b.inputs)) - _NOT_INPUTS)
    )
    return tuple(
        f"source:{name}@{_release(b, name)}"
        for name in names
        if _digest(_node(a.inputs, name)) != _digest(_node(b.inputs, name))
    )


def _pair_causes(
    a: Snapshot,
    b: Snapshot,
    before: State | None,
    after: State | None,
    by_source: dict[str | None, tuple[str, ...]],
) -> set[str]:
    """Why a statement changed from `a` to `b`: the analyser, the reference, a source
    input, or else, when only its call moved, the genotype."""
    causes: set[str] = set()
    if a.tools.get("analyser") is None or a.tools.get("analyser") != b.tools.get(
        "analyser"
    ):
        causes.add("analyser")
    if a.tools.get("reference_md5") != b.tools.get("reference_md5"):
        causes.add(REFERENCE_CAUSE)
    state = after or before
    source = state.source if state is not None else None
    if source not in by_source:
        by_source[source] = _source_causes(a, b, source)
    causes.update(by_source[source])
    if (
        not causes
        and before is not None
        and after is not None
        and set(_fields(before, after)) <= _GENOTYPE
    ):
        causes.add("genotype")
    return causes


def _chain(snaps: Sequence[Snapshot], base: str, run: str) -> Sequence[Snapshot]:
    order = {s.run_id: i for i, s in enumerate(snaps)}
    for rid in (base, run):
        if rid not in order:
            raise ValueError(f"{rid} is not a usable snapshot")
    if order[base] > order[run]:
        raise ValueError(f"the base {base} is newer than the run {run}")
    return snaps[order[base] : order[run] + 1]


@dataclass
class _Track:
    before: State | None
    after: State | None = None
    causes: set[str] = field(default_factory=set)


def diff(root: Path, base: str, run: str) -> list[Change]:
    """The statements that changed from the snapshot `base` to `run`, by id; each one's
    causes are the union over the consecutive usable snapshots between them in which it
    changed."""
    chain = _chain(snapshots(root), base, run)
    tracks: dict[str, _Track] = {}
    with duckdb.connect() as con:
        for k, snap in enumerate(chain):
            _load(con, root, snap.run_id)
            if k == 0:
                continue
            by_source: dict[str | None, tuple[str, ...]] = {}
            for before, after in _pairs(con):
                was = before or after
                assert was is not None
                track = tracks.pop(was.statement_id, None) or _Track(before=before)
                track.after = after
                if before is None or after is None or _fields(before, after):
                    track.causes |= _pair_causes(
                        chain[k - 1], snap, before, after, by_source
                    )
                tracks[(after or was).statement_id] = track
    # A finding or watch item that became a non-finding (reclassified benign) is gone from
    # the user's statements, but its new classification is the news: keep it as `after`.
    gone = {
        sid for sid, t in tracks.items() if t.before is not None and t.after is None
    }
    for sid, state in states(root, run, gone).items():
        tracks[sid].after = state
    changes: list[Change] = []
    for sid, track in tracks.items():
        before, after = track.before, track.after
        fields = _fields(before, after) if before and after else ()
        if (before is None and after is None) or (before and after and not fields):
            continue  # gone and back, or changed and back
        changes.append(Change(sid, before, after, fields, tuple(sorted(track.causes))))
    return sorted(changes, key=lambda c: c.statement_id)


def trusted_ids(root: Path, run_id: str, f: Filters) -> frozenset[str]:
    """The ids of the findings `query.query` keeps at `f` in the snapshot `run_id`."""
    return frozenset(x["id"] for x in findings(root, run_id, f))


def findings(root: Path, run_id: str, f: Filters) -> list[dict[str, Any]]:
    """The findings `query.query` keeps at `f` in the snapshot `run_id`, read as an
    in-memory store with the store's schema: each table a view of the schema's empty table
    and the snapshot's file, by column name (a column an older schema lacked is null)."""
    with duckdb.connect() as con:
        con.execute("ATTACH ':memory:' AS shape")
        con.execute("USE shape")
        con.execute(files("genokeeper.store").joinpath("schema.sql").read_text())
        con.execute("USE memory")
        for table in TABLES:
            path = root / run_id / f"{table}.parquet"
            rows = (
                f" UNION ALL BY NAME SELECT * FROM read_parquet({_literal(str(path))})"
                if path.is_file()
                else ""
            )
            con.execute(f"CREATE VIEW {table} AS SELECT * FROM shape.{table}{rows}")
        return query(con, run_id, f).findings


def kind_of(
    change: Change, f: Filters, before_trusted: bool, after_trusted: bool
) -> ChangeKind:
    """The kind of `change`, given whether it was and is trusted at the filters `f`. A
    research association is never trusted here: the feed lists one only when a filter
    names it, so it neither enters nor leaves. Without a field change (both states there)
    the change is its publications'."""
    state = change.after or change.before
    if state is not None and is_research_layer(state.source):
        before_trusted = after_trusted = False
    if change.before is None:
        return ChangeKind.ENTERED if after_trusted else ChangeKind.NEW
    if change.after is None or change.after.kind == Kind.NON_FINDING:
        return ChangeKind.LEFT if before_trusted else ChangeKind.GONE
    if before_trusted != after_trusted:
        return ChangeKind.ENTERED if after_trusted else ChangeKind.LEFT
    if not change.fields:
        return ChangeKind.PUBLICATIONS
    if _ASSESSMENT & set(change.fields):
        return ChangeKind.LEVEL
    return ChangeKind.DETAILS


@dataclass
class _Record:
    first_found: date | None = None
    last_changed: date | None = None
    assessment_changed: date | None = None
    in_oldest: bool = False


def dates(root: Path, until: str | None = None) -> dict[str, Dates]:
    """Each statement's three dates over the usable snapshots under `root` up to the one
    `until` names (every one when None; none when it names no usable snapshot), by its id
    in the latest snapshot that held it. A non-finding is read as absent."""
    snaps = snapshots(root)
    names = [s.run_id for s in snaps]
    if until is not None:
        snaps = snaps[: names.index(until) + 1] if until in names else []
    if not snaps:
        return {}
    oldest = local_day(snaps[0])
    records: dict[str, _Record] = {}
    with duckdb.connect() as con:
        for k, snap in enumerate(snaps):
            day = local_day(snap)
            _load(con, root, snap.run_id)
            if k == 0:
                for sid, kind in con.execute(
                    "SELECT statement_id, kind FROM after"
                ).fetchall():
                    records[sid] = _Record(
                        first_found=day if kind == Kind.FINDING else None,
                        last_changed=day,
                        assessment_changed=day,
                        in_oldest=True,
                    )
                continue
            for before, after in _pairs(con):
                if after is None:
                    assert before is not None
                    records[before.statement_id].last_changed = day
                    continue
                if (
                    before is not None
                    and before.statement_id != after.statement_id
                    and after.statement_id not in records
                ):
                    records[after.statement_id] = records.pop(before.statement_id)
                record = records.setdefault(after.statement_id, _Record())
                changed = set(COMPARED if before is None else _fields(before, after))
                if changed:
                    record.last_changed = day
                if changed & _ASSESSMENT:
                    record.assessment_changed = day
                if after.kind == Kind.FINDING and record.first_found is None:
                    record.first_found = day
    return {
        sid: Dates(
            first_found=r.first_found,
            last_changed=r.last_changed,
            assessment_changed=r.assessment_changed,
            on_or_before=r.in_oldest
            and oldest in (r.first_found, r.last_changed, r.assessment_changed),
        )
        for sid, r in records.items()
    }


def read_publications(path: Path | None) -> dict[str, dict[str, list[str]]]:
    """A refresh's `publications.json` as `{statement_id: {"litvar": [pmid], "clinvar":
    [pmid]}}`, PMIDs of digits only; empty without one."""
    if path is None or not path.is_file():
        return {}
    document = json.loads(path.read_text())
    return {
        str(sid): {
            name: [str(p) for p in lists.get(name) or [] if _PMID.fullmatch(str(p))]
            for name in ("litvar", "clinvar")
        }
        for sid, lists in document.items()
        if isinstance(lists, dict)
    }
