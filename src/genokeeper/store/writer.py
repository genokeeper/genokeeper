"""A run's row writes, taken in bulk.

DuckDB keeps about a third of a MiB per `INSERT ... ON CONFLICT` statement until the
transaction commits, so one statement per panel site ran a 270,000-site run out of memory.
`StoreWriter` buffers the rows instead and flushes each table in one statement: the rows go
to a newline-delimited JSON file in a private temporary directory under the run's work
directory, read back by `read_json` with every column typed from the schema, never inferred.

A table that keeps the first row written at its key (`variant`, `call`, each keyed within
one run) keeps the first buffered one, and a row the store already holds wins over both
(`ON CONFLICT DO NOTHING`). `statement`, `gene_call` and `coverage` are inserted plainly: a
second row at a key raises, as it did row by row. A column the schema gives a default
(`statement.exclusion`) takes it where the row leaves it out.

A statement flagged `exclusion` (an important negative) whose claim names no confirming test
(`confirm`) is refused when it is added. Its anchor comes from this writer, the one path a
run's statements take into the store, which a source flagging an exclusion cannot move
because the check reads the claim the source wrote, not a flag it sets.
"""

import json
import shutil
import tempfile
from collections.abc import Mapping
from pathlib import Path
from types import TracebackType
from typing import Self

import duckdb

# The tables a run writes row by row, in flush order: a variant before the rows naming it.
BUFFERED = ("variant", "call", "statement", "gene_call", "coverage")
# The tables whose rows were written `ON CONFLICT DO NOTHING`: the first row at a key wins.
_KEEP_FIRST = frozenset({"variant", "call"})
# Flush once this many rows are buffered, so a run's rows never all sit in memory.
MAX_ROWS = 100_000
# The temporary row directories under the work directory start with this. A killed run
# (SIGKILL, out of memory) leaves its own behind, holding the user's rows as plain text.
ROWS_PREFIX = ".store-rows-"


class StoreWriter:
    """Buffered inserts into `con`'s run tables, flushed in bulk inside the caller's
    transaction; the temporary files go under `work`.

    `flush()` writes every buffered row; reads that must see them come after it. As a
    `with` block, the writer starts by flushing what was buffered before, flushes the
    block's rows when it ends, and discards them when it raises: the caller's transaction
    then rolls back what was flushed.
    """

    def __init__(
        self, con: duckdb.DuckDBPyConnection, work: Path, *, max_rows: int = MAX_ROWS
    ) -> None:
        self.con = con
        self._work = work
        self._max_rows = max_rows
        self._rows: dict[str, list[dict[str, object]]] = {t: [] for t in BUFFERED}
        self._keys: dict[str, set[tuple[object, ...]]] = {t: set() for t in BUFFERED}
        self._columns: dict[str, dict[str, str]] = {}
        self._defaults: dict[str, dict[str, str]] = {}
        self._primary: dict[str, tuple[str, ...]] = {}

    def __enter__(self) -> Self:
        self.flush()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if exc_type is None:
            self.flush()
        else:
            self.discard()

    @property
    def pending(self) -> int:
        """How many rows are buffered and not yet flushed."""
        return sum(len(rows) for rows in self._rows.values())

    def _schema(self, table: str) -> dict[str, str]:
        """`table`'s columns and their DuckDB types, in schema order."""
        if table not in self._columns:
            found = self.con.execute(
                "SELECT column_name, data_type, column_default"
                " FROM information_schema.columns"
                " WHERE table_schema = 'main' AND table_name = ?"
                " ORDER BY ordinal_position",
                [table],
            ).fetchall()
            key = self.con.execute(
                "SELECT constraint_column_names FROM duckdb_constraints()"
                " WHERE schema_name = 'main' AND table_name = ?"
                " AND constraint_type = 'PRIMARY KEY'",
                [table],
            ).fetchone()
            if not found or key is None:
                raise ValueError(f"the store has no table {table} with a primary key")
            self._columns[table] = {str(c): str(t) for c, t, _ in found}
            self._defaults[table] = {
                str(c): str(d) for c, _, d in found if d is not None
            }
            self._primary[table] = tuple(key[0])
        return self._columns[table]

    def add(self, table: str, row: Mapping[str, object]) -> None:
        """Buffer one row of `table`; a column it leaves out is NULL."""
        if table not in self._rows:
            raise ValueError(f"{table} is not a buffered table: {BUFFERED}")
        columns = self._schema(table)
        unknown = set(row) - set(columns)
        if unknown:
            raise ValueError(f"{table} has no columns {sorted(unknown)}")
        if table == "statement" and row.get("exclusion"):
            _check_exclusion(row)
        if table in _KEEP_FIRST:
            key = tuple(row.get(c) for c in self._primary[table])
            if key in self._keys[table]:
                return
            self._keys[table].add(key)
        self._rows[table].append({c: row.get(c) for c in columns})
        if self.pending >= self._max_rows:
            self.flush()

    def discard(self) -> None:
        """Drop every buffered row."""
        for table in BUFFERED:
            self._rows[table] = []
            self._keys[table] = set()

    def flush(self) -> None:
        """Insert every buffered row, table by table in `BUFFERED` order, one statement
        each; the buffer is empty afterwards, also when an insert raises."""
        batches = {t: self._rows[t] for t in BUFFERED if self._rows[t]}
        self.discard()
        if not batches:
            return
        self._work.mkdir(parents=True, exist_ok=True)
        # mkdtemp creates the directory mode 700.
        folder = Path(tempfile.mkdtemp(prefix=ROWS_PREFIX, dir=self._work))
        try:
            for table, rows in batches.items():
                self._insert(table, rows, folder / f"{table}.jsonl")
        finally:
            shutil.rmtree(folder, ignore_errors=True)

    def _insert(self, table: str, rows: list[dict[str, object]], path: Path) -> None:
        with path.open("w", encoding="utf-8") as out:
            for row in rows:
                # repr floats: 1.8e-21 and 1e-300 come back as the same double.
                out.write(json.dumps(row) + "\n")
        columns = self._schema(table)
        names = ", ".join(f'"{c}"' for c in columns)
        # A column with a schema default takes it where the row left it out (null).
        defaults = self._defaults[table]
        values = ", ".join(
            f'COALESCE("{c}", {defaults[c]})' if c in defaults else f'"{c}"'
            for c in columns
        )
        # A JSON column is read as its text and cast on insert, as a bound string was: read
        # as JSON, the text would be stored as one JSON string.
        spec = ", ".join(
            f"'{c}': '{'VARCHAR' if t == 'JSON' else t}'" for c, t in columns.items()
        )
        conflict = " ON CONFLICT DO NOTHING" if table in _KEEP_FIRST else ""
        self.con.execute(
            f'INSERT INTO "{table}" ({names}) SELECT {values} FROM read_json(?,'
            f" format = 'newline_delimited', columns = {{{spec}}}){conflict}",
            [str(path)],
        )


def _check_exclusion(row: Mapping[str, object]) -> None:
    """Refuse an exclusion whose claim names no confirming test (a non-empty `confirm`)."""
    claim = row.get("claim")
    if isinstance(claim, str):
        claim = json.loads(claim)
    confirm = claim.get("confirm") if isinstance(claim, Mapping) else None
    if not isinstance(confirm, str) or not confirm.strip():
        raise ValueError(
            f"statement {row.get('statement_id')}: an exclusion needs its claim's"
            " confirming test (`confirm`)"
        )


def sweep_rows(work: Path) -> list[str]:
    """Remove the row directories a killed run left under `work`; the ones that could not
    be removed, each with why. Only while no other run can be writing (under its lock)."""
    failed: list[str] = []
    for folder in sorted(work.glob(f"{ROWS_PREFIX}*")):
        try:
            shutil.rmtree(folder)
        except OSError as e:
            failed.append(f"{folder}: {e}")
    return failed
