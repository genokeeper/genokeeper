import re
from enum import StrEnum
from pathlib import Path

import duckdb
import pytest

from genokeeper import history
from genokeeper.query import Filters
from genokeeper.store import db
from genokeeper.store.db import (
    SCHEMA_VERSION,
    AccountingError,
    Kind,
    Level,
    RowType,
    StoreTooNew,
    StoreTooOld,
    Topic,
    Validity,
    connect,
    export_snapshot,
    export_snapshots,
    record_accounting,
    start_run,
    transaction,
)

TABLES = (
    "variant",
    "call",
    "statement",
    "run",
    "accounting",
    "gene_call",
    "coverage",
)
VID = "GRCh38:chr1:169549811:C:T"
COUNTS = {
    "input_rows": 12,
    "findings": 3,
    "watch": 2,
    "non_findings": 4,
    "unassessed": 2,
}


@pytest.mark.parametrize(
    ("table", "column", "enum"),
    [
        ("statement", "kind", Kind),
        ("statement", "validity", Validity),
        ("statement", "topic", Topic),
        ("call", "level", Level),
        ("call", "row_type", RowType),
        ("gene_call", "level", Level),
    ],
)
def test_enumerated_columns_check_the_enum_values(
    tmp_path: Path, table: str, column: str, enum: type[StrEnum]
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    rows = con.execute(
        "SELECT expression FROM duckdb_constraints()"
        " WHERE table_name = ? AND constraint_type = 'CHECK'"
        " AND constraint_column_names = [?]",
        [table, column],
    ).fetchall()
    assert len(rows) == 1, rows
    assert set(re.findall(r"'([^']*)'", rows[0][0])) == {e.value for e in enum}


def test_schema_3_holds_grch38_runs_only(tmp_path: Path) -> None:
    """No hg19 record, call mapping, OpenCRAVAT freeze or OpenCRAVAT table; a reference
    block is a row type of its own."""
    con = connect(tmp_path / "findings.duckdb")

    def columns(table: str) -> list[str]:
        return [
            row[0]
            for row in con.execute(
                "SELECT column_name FROM information_schema.columns"
                " WHERE table_name = ? ORDER BY ordinal_position",
                [table],
            ).fetchall()
        ]

    assert columns("variant") == ["run_id", "vid", "rsid"]
    assert "mapping" not in columns("call")
    assert columns("run") == ["run_id", "started_at", "inputs", "tools"]
    names = {
        row[0]
        for row in con.execute(
            "SELECT table_name FROM information_schema.tables"
        ).fetchall()
    }
    assert names.isdisjoint({"oc_clinvar_match", "clinvar_crosscheck"})
    insert = (
        "INSERT INTO call (run_id, vid, dosage, level, row_type)"
        " VALUES ('r', ?, '0', 'medium', ?)"
    )
    con.execute(insert, ["v1", "ref_block"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["v2", "secondary"])


@pytest.mark.parametrize(
    ("column", "value"), [("validity", "Strong"), ("topic", "risk")]
)
def test_a_bad_validity_or_topic_is_refused(
    tmp_path: Path, column: str, value: str
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    insert = (
        f"INSERT INTO statement (run_id, statement_id, kind, {column})"
        " VALUES ('r', ?, 'watch', ?)"
    )
    con.execute(insert, ["s0", None])  # a statement without one is allowed
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["s1", value])


def test_a_store_records_its_schema_version_once(tmp_path: Path) -> None:
    db = tmp_path / "findings.duckdb"
    connect(db).close()
    con = connect(db)
    assert con.execute("SELECT schema_version FROM meta").fetchall() == [
        (SCHEMA_VERSION,)
    ]
    assert SCHEMA_VERSION == 7


def test_a_store_older_than_the_code_is_refused_untouched(tmp_path: Path) -> None:
    db = tmp_path / "findings.duckdb"
    con = connect(db)
    con.execute("UPDATE meta SET schema_version = ?", [SCHEMA_VERSION - 1])
    con.execute("DROP VIEW finding")  # an older schema's own shape
    con.close()
    with pytest.raises(StoreTooOld) as refused:
        connect(db)
    message = str(refused.value)
    assert len(message.splitlines()) == 1
    assert str(db) in message
    assert f"schema version {SCHEMA_VERSION - 1}" in message
    assert "move it aside and rerun `genokeeper run`" in message
    con = duckdb.connect(str(db))
    names = {
        row[0]
        for row in con.execute(
            "SELECT table_name FROM information_schema.tables"
        ).fetchall()
    }
    assert "finding" not in names  # this code's schema was never applied over it
    assert con.execute("SELECT schema_version FROM meta").fetchall() == [
        (SCHEMA_VERSION - 1,)
    ]


def test_variant_is_keyed_per_run_without_a_vrs_id(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    key = con.execute(
        "SELECT constraint_column_names FROM duckdb_constraints()"
        " WHERE table_name = 'variant' AND constraint_type = 'PRIMARY KEY'"
    ).fetchone()
    assert key == (["run_id", "vid"],)
    columns = [
        row[0]
        for row in con.execute(
            "SELECT column_name FROM information_schema.columns"
            " WHERE table_name = 'variant'"
        ).fetchall()
    ]
    assert "vrs_id" not in columns
    # Two runs may each read one allele their own way (another release's rsID).
    for run_id, rsid in (("r1", "rs1"), ("r2", "rs2")):
        con.execute(
            "INSERT INTO variant (run_id, vid, rsid) VALUES (?, ?, ?)",
            [run_id, VID, rsid],
        )
    assert con.execute("SELECT run_id, rsid FROM variant ORDER BY ALL").fetchall() == [
        ("r1", "rs1"),
        ("r2", "rs2"),
    ]


def test_a_store_newer_than_the_code_is_refused_untouched(tmp_path: Path) -> None:
    db = tmp_path / "findings.duckdb"
    con = connect(db)
    con.execute("UPDATE meta SET schema_version = ?", [SCHEMA_VERSION + 1])
    con.execute("DROP VIEW finding")  # a newer schema's own shape
    con.close()
    with pytest.raises(StoreTooNew, match=f"schema version {SCHEMA_VERSION + 1}"):
        connect(db)
    con = duckdb.connect(str(db))
    names = {
        row[0]
        for row in con.execute(
            "SELECT table_name FROM information_schema.tables"
        ).fetchall()
    }
    assert "finding" not in names  # this code's schema was never applied over it


def test_bad_kind_is_refused(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    insert = "INSERT INTO statement (run_id, statement_id, kind) VALUES ('r', ?, ?)"
    con.execute(insert, ["s1", "watch"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["s2", "findings"])


def test_accounting_is_keyed_by_run_and_tool(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    record_accounting(con, run_id, "panel", **COUNTS, dropped={"x": 1})
    record_accounting(con, run_id, "mito", **COUNTS, dropped={"x": 1})
    with pytest.raises(duckdb.ConstraintException):
        record_accounting(con, run_id, "panel", **COUNTS, dropped={"x": 1})


# Every case balances, so only a sign can be what is refused.
@pytest.mark.parametrize(
    ("counts", "dropped"),
    [
        (dict.fromkeys(COUNTS, 0) | {"input_rows": -1}, {"x": -1}),
        ({**COUNTS, "input_rows": 8, "findings": -1}, {"x": 1}),
        ({**COUNTS, "input_rows": 8, "watch": -2}, {"x": 1}),
        ({**COUNTS, "input_rows": 4, "non_findings": -4}, {"x": 1}),
        ({**COUNTS, "input_rows": 8, "unassessed": -2}, {"x": 1}),
        (COUNTS, {"x": 3, "y": -2}),
    ],
)
def test_accounting_refuses_negative_counts(
    tmp_path: Path, counts: dict[str, int], dropped: dict[str, int]
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    with pytest.raises(ValueError, match="negative"):
        record_accounting(con, run_id, "panel", **counts, dropped=dropped)
    assert con.execute("SELECT count(*) FROM accounting").fetchone() == (0,)


def test_finding_view_joins_calls_and_keeps_gene_findings(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    con.execute("INSERT INTO variant (run_id, vid) VALUES (?, ?)", [run_id, VID])
    con.execute(
        "INSERT INTO call (run_id, vid, gt, gq, level, row_type)"
        " VALUES (?, ?, '0/1', 50, 'high', 'PASS')",
        [run_id, VID],
    )
    statement = (
        "INSERT INTO statement (run_id, statement_id, kind, vid, gene, validity,"
        " source_level, source, source_version, topic)"
        " VALUES (?, ?, ?, ?, ?, 'strong', 'level', 'src', 'v1', ?)"
    )
    con.execute(statement, [run_id, "clinvar:642", "finding", VID, "F5", "risks"])
    con.execute(
        statement, [run_id, "pharmcat:CYP2C19", "finding", None, "CYP2C19", "drugs"]
    )
    con.execute(statement, [run_id, "clinvar:999", "watch", VID, "F5", "risks"])
    rows = con.execute(
        "SELECT statement_id, gt, gq, level FROM finding ORDER BY statement_id"
    ).fetchall()
    assert rows == [
        ("clinvar:642", "0/1", 50, "high"),
        ("pharmcat:CYP2C19", None, None, None),
    ]


def test_schema_is_idempotent(tmp_path: Path) -> None:
    db = tmp_path / "findings.duckdb"
    connect(db).close()
    con = connect(db)
    names = {
        row[0]
        for row in con.execute(
            "SELECT table_name FROM information_schema.tables"
        ).fetchall()
    }
    assert set(TABLES) | {"finding", "unassessed", "meta"} <= names


def test_unassessed_view_holds_only_unassessed_claims_with_their_call(
    tmp_path: Path,
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    con.execute(
        "INSERT INTO call (run_id, vid, gt, dosage, level, row_type)"
        " VALUES (?, ?, './.', 'unknown', 'low', 'NoCall')",
        [run_id, VID],
    )
    statement = (
        "INSERT INTO statement (run_id, statement_id, kind, vid, gene, validity,"
        " source_level, source, source_version, topic)"
        " VALUES (?, ?, ?, ?, 'MUTYH', 'strong', 'level', 'src', 'v1', 'risks')"
    )
    con.execute(statement, [run_id, "clinvar:1", "unassessed", VID])
    con.execute(statement, [run_id, "clinvar:2", "finding", VID])
    con.execute(statement, [run_id, "clinvar:3", "watch", VID])
    assert con.execute(
        "SELECT statement_id, kind, dosage FROM unassessed"
    ).fetchall() == [("clinvar:1", "unassessed", "unknown")]
    assert con.execute("SELECT statement_id FROM finding").fetchall() == [
        ("clinvar:2",)
    ]
    finding_columns = con.execute("SELECT * FROM finding LIMIT 0").description
    unassessed_columns = con.execute("SELECT * FROM unassessed LIMIT 0").description
    assert [c[0] for c in finding_columns] == [c[0] for c in unassessed_columns]


@pytest.mark.parametrize(
    "missing", ["validity", "source_level", "source", "source_version", "topic", "vid"]
)
def test_an_unassessed_claim_needs_a_findings_attributes_and_a_variant(
    tmp_path: Path, missing: str
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    values = {
        "vid": VID,
        "validity": "strong",
        "source_level": "2 stars",
        "source": "clinvar",
        "source_version": "2026-10",
        "topic": "risks",
    }
    insert = (
        "INSERT INTO statement (run_id, statement_id, kind, vid, validity, source_level,"
        " source, source_version, topic) VALUES ('r', ?, 'unassessed', ?, ?, ?, ?, ?, ?)"
    )
    con.execute(insert, ["clinvar:1", *values.values()])
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["clinvar:2", *(values | {missing: None}).values()])


def test_a_gene_level_unassessed_claim_names_its_gene_or_its_reason(
    tmp_path: Path,
) -> None:
    """Schema 5: an unresolved gene call is a gene-level unassessed statement (no vid),
    allowed when it names its gene; its gene call has a row of its own. Schema 6: one with
    neither vid nor gene (a polygenic score below its overlap) is allowed when it carries
    its reason."""
    con = connect(tmp_path / "findings.duckdb")
    insert = (
        "INSERT INTO statement (run_id, statement_id, kind, vid, gene, reason, validity,"
        " source_level, source, source_version, topic)"
        " VALUES ('r', ?, 'unassessed', NULL, ?, ?, 'moderate', 'CPIC Moderate',"
        " 'pharmcat', '3.4.0', 'drugs')"
    )
    con.execute(insert, ["pharmcat:CYP2D6", "CYP2D6", None])
    con.execute(insert, ["pgsc_calc:PGS000018", None, "score overlap 0.60, below 0.75"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["pharmcat:none", None, None])
    con.execute(
        "INSERT INTO gene_call (run_id, gene, diplotype, phenotype, activity_score,"
        " callers, agreed, level, warnings) VALUES ('r', 'CYP2D6', NULL, NULL, NULL,"
        ' \'{"BCyrius": "*1/*4", "Aldy": "*2/*4"}\', false, \'low\', [\'unresolved\'])'
    )
    with pytest.raises(duckdb.ConstraintException):  # one gene call per gene and run
        con.execute(
            "INSERT INTO gene_call (run_id, gene, level) VALUES ('r', 'CYP2D6', 'low')"
        )
    with pytest.raises(duckdb.ConstraintException):
        con.execute(
            "INSERT INTO gene_call (run_id, gene, level) VALUES ('r', 'CYP2B6', 'Low')"
        )


def test_finding_requires_source_version(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    insert = (
        "INSERT INTO statement (run_id, statement_id, kind, gene, validity, source_level,"
        " source, source_version, topic) VALUES ('r', ?, 'finding', 'F5', 'strong',"
        " '2 stars', 'clinvar', ?, 'risks')"
    )
    con.execute(insert, ["clinvar:1", "2026-10"])
    with pytest.raises(duckdb.ConstraintException):
        con.execute(insert, ["clinvar:2", None])


def test_accounting_must_balance(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {"vcf": "sha256:0"}, {"bcftools": "1.24"})
    record_accounting(
        con,
        run_id,
        "panel",
        input_rows=12,
        findings=3,
        watch=2,
        non_findings=4,
        unassessed=2,
        dropped={"x": 1},
    )
    for input_rows, unassessed in ((13, 2), (12, 0)):
        # One row too many, or the unassessed rows left uncounted.
        with pytest.raises(AccountingError, match="unassessed"):
            record_accounting(
                con,
                run_id,
                "mito",
                input_rows=input_rows,
                findings=3,
                watch=2,
                non_findings=4,
                unassessed=unassessed,
                dropped={"x": 1},
            )
    assert con.execute(
        "SELECT input_rows, findings, watch, non_findings, unassessed FROM accounting"
    ).fetchall() == [(12, 3, 2, 4, 2)]


def test_run_id_shape(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    assert re.fullmatch(r"\d{8}T\d{6}Z-[0-9a-f]{6}", run_id)
    assert con.execute("SELECT run_id FROM run").fetchall() == [(run_id,)]


def test_export_snapshot_writes_one_parquet_per_table(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    other = start_run(con, {}, {})
    for rid, rsid in ((run_id, "rs6025"), (other, "rs-other-release")):
        con.execute(
            "INSERT INTO variant (run_id, vid, rsid)"
            " VALUES (?, 'GRCh38:chr1:169549811:C:T', ?)",
            [rid, rsid],
        )
    for rid in (run_id, other):
        con.execute(
            "INSERT INTO call (run_id, vid, level, row_type)"
            " VALUES (?, 'GRCh38:chr1:169549811:C:T', 'high', 'PASS')",
            [rid],
        )
    snapshot = export_snapshot(con, run_id, tmp_path / "snapshots")
    assert sorted(p.name for p in snapshot.iterdir()) == sorted(
        f"{t}.parquet" for t in (*TABLES, "meta")
    )
    calls = con.execute(
        "SELECT run_id FROM read_parquet(?)", [str(snapshot / "call.parquet")]
    ).fetchall()
    assert calls == [(run_id,)]
    variants = con.execute(
        "SELECT run_id, rsid FROM read_parquet(?)", [str(snapshot / "variant.parquet")]
    ).fetchall()
    assert variants == [(run_id, "rs6025")]
    meta = con.execute(
        "SELECT schema_version FROM read_parquet(?)", [str(snapshot / "meta.parquet")]
    ).fetchall()
    assert meta == [(SCHEMA_VERSION,)]


def test_a_failed_export_is_repaired_by_the_next(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    con = connect(tmp_path / "findings.duckdb")
    root = tmp_path / "snapshots"
    first = start_run(con, {}, {})
    real = duckdb.DuckDBPyRelation.write_parquet

    def full_disk(self: duckdb.DuckDBPyRelation, path: str, **kw: object) -> None:
        if path.endswith("call.parquet"):
            raise OSError("No space left on device")
        real(self, path, **kw)  # pyright: ignore[reportArgumentType]

    with monkeypatch.context() as m:
        m.setattr(duckdb.DuckDBPyRelation, "write_parquet", full_disk)
        with pytest.raises(OSError, match="No space"):
            export_snapshots(con, root)
    assert not (root / first).exists()  # no partial snapshot reads as whole
    second = start_run(con, {}, {})
    assert set(export_snapshots(con, root)) == {root / first, root / second}
    assert sorted(p.name for p in (root / first).iterdir()) == sorted(
        f"{t}.parquet" for t in (*TABLES, "meta")
    )
    assert export_snapshots(con, root) == []


def _runs(con: duckdb.DuckDBPyConnection) -> int:
    row = con.execute("SELECT count(*) FROM run").fetchone()
    assert row is not None
    return row[0]


def test_transaction_commits_or_rolls_back(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    with transaction(con):
        start_run(con, {}, {})
    assert _runs(con) == 1
    with pytest.raises(RuntimeError), transaction(con):
        start_run(con, {}, {})
        raise RuntimeError("abort")
    assert _runs(con) == 1


def test_a_nested_transaction_joins_the_outer_one(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    with pytest.raises(RuntimeError), transaction(con):
        with transaction(con):
            start_run(con, {}, {})
        # The inner block's end committed nothing: the outer one rolls it back.
        raise RuntimeError("abort")
    assert _runs(con) == 0
    with transaction(con):
        with transaction(con):
            start_run(con, {}, {})
        start_run(con, {}, {})
    assert _runs(con) == 2


def test_store_has_coverage_v7(tmp_path: Path) -> None:
    assert SCHEMA_VERSION == 7
    assert "coverage" in db.TABLES
    con = connect(tmp_path / "findings.duckdb")
    key = con.execute(
        "SELECT constraint_column_names FROM duckdb_constraints()"
        " WHERE table_name = 'coverage' AND constraint_type = 'PRIMARY KEY'"
    ).fetchone()
    assert key == (["run_id", "gene"],)
    con.execute("UPDATE meta SET schema_version = 6")
    con.close()
    with pytest.raises(StoreTooOld, match="schema version 6"):
        connect(tmp_path / "findings.duckdb")


def test_statement_exclusion_defaults_false(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    con.execute("INSERT INTO variant (run_id, vid) VALUES (?, ?)", [run_id, VID])
    con.execute(
        "INSERT INTO statement (run_id, statement_id, kind, vid, gene, validity,"
        " source_level, source, source_version, topic)"
        " VALUES (?, 'clinvar:1', 'finding', ?, 'BRCA1', 'strong', 'Pathogenic (2★)',"
        " 'clinvar', '20261011', 'risks')",
        [run_id, VID],
    )
    assert con.execute("SELECT exclusion FROM statement").fetchall() == [(False,)]
    snapshot = export_snapshot(con, run_id, tmp_path / "snapshots")
    # A schema-6 snapshot: its statement file has no `exclusion`, and no coverage file.
    statements = snapshot / "statement.parquet"
    con.execute(
        f"COPY (SELECT * EXCLUDE (exclusion) FROM read_parquet('{statements}'))"
        f" TO '{statements}' (FORMAT parquet)"
    )
    (snapshot / "coverage.parquet").unlink()
    found = history.findings(tmp_path / "snapshots", run_id, Filters())
    assert [f["id"] for f in found] == ["clinvar:1"]
