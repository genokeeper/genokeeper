"""The DuckDB findings store: schema, runs, row accounting and Parquet snapshots.

A store records the schema version that created it (`meta.schema_version`, one row). This
code writes `SCHEMA_VERSION` into a new store and refuses, before touching its schema, one
another version wrote: a newer one (`StoreTooNew`) and an older one (`StoreTooOld`).
There are no migrations (track design, "Schema changes"): an older store is moved aside and
the run repeated; every run's Parquet snapshot keeps the history.
"""

import collections.abc
import json
import re
import secrets
import shutil
import weakref
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from enum import StrEnum
from importlib.resources import files
from pathlib import Path

import duckdb


class Level(StrEnum):
    """Display level of a call's genotype reliability."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class Validity(StrEnum):
    """Clinical validity of a statement, mapped from the source's own level."""

    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    CONFLICTING = "conflicting"


class RowType(StrEnum):
    """Where a call came from: a row type of our own gVCF (a reference block among them; on
    chrM, mutserve's call, `genokeeper.mito`), or no row."""

    PASS = "PASS"
    REFCALL = "RefCall"
    REF_BLOCK = "ref_block"
    NOCALL = "NoCall"
    NONE = "none"


class Kind(StrEnum):
    """What a statement is for the user."""

    FINDING = "finding"
    WATCH = "watch"
    NON_FINDING = "non_finding"
    # A claim that would be a finding, at a site whose genotype is unknown (or has no call).
    UNASSESSED = "unassessed"


class Topic(StrEnum):
    """The dashboard topic a finding sits under."""

    RISKS = "risks"
    CARRIER = "carrier"
    DRUGS = "drugs"
    TRAITS = "traits"
    WELLNESS = "wellness"
    ANCESTRY = "ancestry"
    RESEARCH = "research"


class AccountingError(ValueError):
    """A tool's row accounting does not add up to its input rows."""


class StoreTooNew(RuntimeError):
    """The store was written by a newer genokeeper, whose schema this code does not know."""


class StoreTooOld(RuntimeError):
    """The store was written by an older genokeeper, whose schema this code does not migrate."""


# The schema version this code writes and reads: bump with every change to schema.sql.
# 2: `variant` per run (`run_id` in its key, no `vrs_id`); kind `unassessed`, its view and
# its `accounting` column.
# 3: GRCh38 only: `variant` without its hg19 record, `call` without `mapping`, row type
# `ref_block` (and no `secondary`), `run` without `oc_freeze`; no `oc_clinvar_match` or
# `clinvar_crosscheck`.
# 4: `call.p_wrong` measured on HG003 (stage 2b), with its interval `p_wrong_lo` and
# `p_wrong_hi`.
# 5: `gene_call` (slice 1b), one pharmacogenomic gene call per gene and run; an unassessed
# statement names its vid or its gene.
# 6: an unassessed statement with neither vid nor gene carries its reason (slice 1d).
# 7: `coverage` (slice 3c), one row per run and gene; `statement.exclusion`.
SCHEMA_VERSION = 7


TABLES = ("variant", "call", "statement", "run", "accounting", "gene_call", "coverage")
_RUN_ID = re.compile(r"\d{8}T\d{6}Z-[0-9a-f]{6}")


def _schema_version(con: duckdb.DuckDBPyConnection) -> int | None:
    """The store's recorded schema version; None for a store that records none (new)."""
    has_meta = con.execute(
        "SELECT count(*) FROM information_schema.tables WHERE table_name = 'meta'"
    ).fetchone()
    if has_meta is None or has_meta[0] == 0:
        return None
    row = con.execute("SELECT max(schema_version) FROM meta").fetchone()
    return None if row is None or row[0] is None else int(row[0])


def check_schema_version(con: duckdb.DuckDBPyConnection, path: Path) -> int | None:
    """The schema version of the store `con` opened at `path`, None for a new one; refuse a
    store another version wrote, newer (`StoreTooNew`) or older (`StoreTooOld`). Only
    reads."""
    found = _schema_version(con)
    if found is not None and found > SCHEMA_VERSION:
        raise StoreTooNew(
            f"{path} has schema version {found}; this genokeeper reads up to"
            f" {SCHEMA_VERSION}: update genokeeper"
        )
    if found is not None and found < SCHEMA_VERSION:
        raise StoreTooOld(
            f"{path} has schema version {found}, older than this genokeeper's"
            f" {SCHEMA_VERSION}, and is not migrated: move it aside and rerun"
            " `genokeeper run`"
        )
    return found


def connect(path: Path) -> duckdb.DuckDBPyConnection:
    """Open (or create) the store at `path`, create its schema if missing and record
    `SCHEMA_VERSION` in a new one; refuse, untouched, a store another schema version
    wrote (`StoreTooNew`, `StoreTooOld`)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(path))
    try:
        found = check_schema_version(con, path)
        con.execute(files("genokeeper.store").joinpath("schema.sql").read_text())
        if found is None:
            con.execute(
                "INSERT INTO meta (schema_version) VALUES (?)", [SCHEMA_VERSION]
            )
    except BaseException:
        con.close()
        raise
    return con


# Connections inside a `transaction` block. DuckDB has no way to ask, and a BEGIN inside a
# transaction aborts that transaction, so the blocks keep track themselves.
_IN_TRANSACTION: weakref.WeakSet[duckdb.DuckDBPyConnection] = weakref.WeakSet()


@contextmanager
def transaction(con: duckdb.DuckDBPyConnection) -> Iterator[None]:
    """Run the block in one transaction: committed when it ends, rolled back when it raises.

    A block inside another `transaction(con)` block joins it, and the outermost one alone
    commits or rolls back. A transaction begun with a bare `con.begin()` is not seen.
    """
    if con in _IN_TRANSACTION:
        yield
        return
    con.begin()
    _IN_TRANSACTION.add(con)
    try:
        yield
    except BaseException:
        con.rollback()
        raise
    else:
        con.commit()
    finally:
        _IN_TRANSACTION.discard(con)


def start_run(
    con: duckdb.DuckDBPyConnection,
    inputs: collections.abc.Mapping[str, object],
    tools: collections.abc.Mapping[str, str],
) -> str:
    """Record a new run and return its id, `<UTC %Y%m%dT%H%M%SZ>-<6 hex>`."""
    started = datetime.now(UTC)
    run_id = f"{started:%Y%m%dT%H%M%SZ}-{secrets.token_hex(3)}"
    con.execute(
        "INSERT INTO run (run_id, started_at, inputs, tools) VALUES (?, ?, ?, ?)",
        [
            run_id,
            started.replace(tzinfo=None),
            json.dumps(inputs, sort_keys=True),
            json.dumps(tools, sort_keys=True),
        ],
    )
    return run_id


def record_accounting(
    con: duckdb.DuckDBPyConnection,
    run_id: str,
    tool: str,
    *,
    input_rows: int,
    findings: int,
    watch: int,
    non_findings: int,
    unassessed: int,
    dropped: dict[str, int],
) -> None:
    """Record where a tool's input rows went; refuse a negative count or an unbalanced tally:
    `input_rows = findings + watch + non_findings + unassessed + sum(dropped)`."""
    counts = {
        "input_rows": input_rows,
        "findings": findings,
        "watch": watch,
        "non_findings": non_findings,
        "unassessed": unassessed,
        **{f"dropped[{reason!r}]": n for reason, n in dropped.items()},
    }
    negative = {name: n for name, n in counts.items() if n < 0}
    if negative:
        raise ValueError(f"{tool} in run {run_id}: negative counts {negative}")
    accounted = findings + watch + non_findings + unassessed + sum(dropped.values())
    if input_rows != accounted:
        raise AccountingError(
            f"{tool} in run {run_id}: {input_rows} input rows, but {findings} findings"
            f" + {watch} watch + {non_findings} non-findings + {unassessed} unassessed"
            f" + {dropped} dropped = {accounted}"
        )
    con.execute(
        "INSERT INTO accounting (run_id, tool, input_rows, findings, watch, non_findings,"
        " unassessed, dropped) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        [
            run_id,
            tool,
            input_rows,
            findings,
            watch,
            non_findings,
            unassessed,
            json.dumps(dropped, sort_keys=True),
        ],
    )


def export_snapshot(con: duckdb.DuckDBPyConnection, run_id: str, out_dir: Path) -> Path:
    """Write the run's rows of every table to `<out_dir>/<run_id>/<table>.parquet`.

    `meta` has no run: the snapshot holds the store's schema version.
    """
    if not _RUN_ID.fullmatch(run_id):
        raise ValueError(f"not a run id: {run_id!r}")
    if con.execute("SELECT 1 FROM run WHERE run_id = ?", [run_id]).fetchone() is None:
        raise KeyError(f"no run {run_id} in the store")
    snapshot = out_dir / run_id
    # Written beside it and renamed into place: a snapshot folder is a whole one.
    partial = out_dir / f".{run_id}.partial"
    shutil.rmtree(partial, ignore_errors=True)
    partial.mkdir(parents=True)
    queries = {
        "meta": "SELECT * FROM meta",
        **{table: f"SELECT * FROM {table} WHERE run_id = $run_id" for table in TABLES},
    }
    for table, query in queries.items():
        params = {"run_id": run_id} if "$run_id" in query else None
        con.sql(query, params=params).write_parquet(str(partial / f"{table}.parquet"))
    shutil.rmtree(snapshot, ignore_errors=True)
    partial.rename(snapshot)
    return snapshot


def export_snapshots(con: duckdb.DuckDBPyConnection, out_dir: Path) -> list[Path]:
    """`export_snapshot` for each run of the store with no snapshot under `out_dir` (the
    new run's, and any an export that failed after its commit left out), oldest first."""
    runs = con.execute("SELECT run_id FROM run ORDER BY started_at, run_id").fetchall()
    return [
        export_snapshot(con, run_id, out_dir)
        for (run_id,) in runs
        if not (out_dir / run_id).is_dir()
    ]
