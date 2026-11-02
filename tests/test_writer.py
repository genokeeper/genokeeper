import json
import math
import random
import time
from pathlib import Path

import duckdb
import pytest

from genokeeper.sources.panel import clinvar_claim
from genokeeper.store import writer as writer_module
from genokeeper.store.db import connect, start_run, transaction
from genokeeper.store.writer import BUFFERED, StoreWriter, sweep_rows

VID = "GRCh38:chr1:169549811:C:T"


def _call(run_id: str, vid: str = VID, **values: object) -> dict[str, object]:
    return {
        "run_id": run_id,
        "vid": vid,
        "gt": "0/1",
        "dosage": "1",
        "p_wrong": 1e-5,
        "level": "high",
        "gq": 50,
        "dp": 30,
        "ad": [15, 15],
        "vaf": 0.5,
        "row_type": "PASS",
        "flags": [],
        "warnings": [],
        **values,
    }


def _statement(run_id: str, statement_id: str, **values: object) -> dict[str, object]:
    return {
        "run_id": run_id,
        "statement_id": statement_id,
        "kind": "non_finding",
        "vid": VID,
        "gene": "MUTYH",
        "claim": clinvar_claim(
            variation_id=5294,
            clnsig="Pathogenic",
            clnrevstat="reviewed_by_expert_panel",
            condition="Familial adenomatous polyposis 2",
        ),
        "source": "clinvar",
        **values,
    }


def _store(tmp_path: Path) -> tuple[duckdb.DuckDBPyConnection, str]:
    con = connect(tmp_path / "findings.duckdb")
    return con, start_run(con, {}, {})


def _all(con: duckdb.DuckDBPyConnection, table: str, run_id: str) -> list[tuple]:
    return con.execute(
        f"SELECT * EXCLUDE (run_id) FROM {table} WHERE run_id = ? ORDER BY ALL",
        [run_id],
    ).fetchall()


def test_the_first_buffered_row_at_a_key_wins(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    writer.add("call", _call(run_id, gt="0/1"))
    writer.add("call", _call(run_id, gt="1/1"))
    writer.add("variant", {"run_id": run_id, "vid": VID, "rsid": "rs1"})
    writer.add("variant", {"run_id": run_id, "vid": VID, "rsid": "rs2"})
    assert writer.pending == 2
    writer.flush()
    assert con.execute("SELECT gt FROM call").fetchall() == [("0/1",)]
    assert con.execute("SELECT rsid FROM variant").fetchall() == [("rs1",)]


def test_each_run_keeps_its_own_variant_row(tmp_path: Path) -> None:
    # A later ClinVar release may name the allele otherwise: its run keeps that.
    con, first = _store(tmp_path)
    second = start_run(con, {}, {})
    for run_id, rsid in ((first, "rs1"), (second, "rs2")):
        writer = StoreWriter(con, tmp_path / "work")
        writer.add("variant", {"run_id": run_id, "vid": VID, "rsid": rsid})
        writer.add("variant", {"run_id": run_id, "vid": VID, "rsid": f"{rsid}-later"})
        writer.flush()
    assert con.execute("SELECT run_id, rsid FROM variant ORDER BY rsid").fetchall() == [
        (first, "rs1"),
        (second, "rs2"),
    ]


def test_a_row_the_store_holds_wins_over_a_buffered_one(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    con.execute(
        "INSERT INTO variant (run_id, vid, rsid) VALUES (?, ?, 'stored')", [run_id, VID]
    )
    writer = StoreWriter(con, tmp_path / "work")
    writer.add("variant", {"run_id": run_id, "vid": VID, "rsid": "buffered"})
    writer.add("variant", {"run_id": run_id, "vid": "GRCh38:chr1:2:A:G", "rsid": "new"})
    writer.flush()
    assert con.execute("SELECT vid, rsid FROM variant ORDER BY vid").fetchall() == [
        ("GRCh38:chr1:169549811:C:T", "stored"),
        ("GRCh38:chr1:2:A:G", "new"),
    ]


def test_the_first_row_wins_across_a_flush_the_buffer_size_forces(
    tmp_path: Path,
) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work", max_rows=2)
    writer.add("call", _call(run_id, gt="0/1"))
    writer.add("call", _call(run_id, vid="GRCh38:chr1:2:A:G"))
    assert writer.pending == 0  # the second row filled the buffer
    writer.add("call", _call(run_id, gt="1/1"))
    writer.flush()
    assert con.execute("SELECT gt FROM call WHERE vid = ?", [VID]).fetchall() == [
        ("0/1",)
    ]


def test_a_second_statement_at_a_key_still_raises(tmp_path: Path) -> None:
    # Statements were inserted without ON CONFLICT: a duplicate is a bug, not a no-op.
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    writer.add("statement", _statement(run_id, "clinvar:1"))
    writer.add("statement", _statement(run_id, "clinvar:1", gene="APC"))
    with pytest.raises(duckdb.ConstraintException):
        writer.flush()


def test_tables_flush_in_order_with_the_variant_first(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    order: list[str] = []
    insert = StoreWriter._insert  # pyright: ignore[reportPrivateUsage]

    def recording(self: StoreWriter, table: str, rows: list, path: Path) -> None:
        order.append(table)
        insert(self, table, rows, path)

    monkeypatch.setattr(StoreWriter, "_insert", recording)
    writer.add("coverage", {"run_id": run_id, "gene": "CBS", "bases": 1004})
    writer.add("gene_call", {"run_id": run_id, "gene": "CYP2C19", "level": "medium"})
    writer.add("statement", _statement(run_id, "clinvar:1"))
    writer.add("call", _call(run_id))
    writer.add("variant", {"run_id": run_id, "vid": VID})
    writer.flush()
    assert order == list(BUFFERED)
    assert BUFFERED[0] == "variant"
    # The schema declares no foreign keys; the order is the variant before what names it.
    assert con.execute(
        "SELECT count(*) FROM duckdb_constraints()"
        " WHERE constraint_type = 'FOREIGN KEY'"
    ).fetchone() == (0,)


def _bound(con: duckdb.DuckDBPyConnection, table: str, row: dict[str, object]) -> None:
    """Insert `row` with bound parameters, as each row was written before the writer."""
    columns = ", ".join(row)
    marks = ", ".join("?" for _ in row)
    con.execute(f"INSERT INTO {table} ({columns}) VALUES ({marks})", list(row.values()))


def test_rows_come_back_as_a_bound_insert_stores_them(tmp_path: Path) -> None:
    con = connect(tmp_path / "findings.duckdb")
    bulk, bound = start_run(con, {}, {}), start_run(con, {}, {})
    rng = random.Random(7)
    doubles = [1.8e-21, 1e-300, 5e-324, 0.1, 1 / 3, math.nan, math.inf]
    doubles += [rng.random() * 10.0 ** rng.randint(-300, 300) for _ in range(50)]
    calls = [
        _call(
            run_id="",
            vid=f"GRCh38:chr1:{n}:A:G",
            p_wrong=p,
            vaf=None if n % 3 == 0 else p,
            gt=None if n % 4 == 0 else "0/1",
            ad=None if n % 5 == 0 else [n, 2**31 - 1],
            flags=["low_gq", "Ünïcødé"] if n % 2 else [],
            warnings=["近い", "🧬 emoji"] if n % 2 else [],
        )
        for n, p in enumerate(doubles)
    ]
    statements = [
        _statement(
            "",
            f"clinvar:{n}",
            gene=gene,
            reason=reason,
            actionable=actionable,
            points=points,
            claim=claim,
            topic=None if n else "risks",
        )
        for n, (gene, reason, actionable, points, claim) in enumerate(
            [
                ("MUTYH", None, True, 1.8e-21, None),
                ("Ünïcødé", "Maladie de Ménière; 遺伝", False, None, None),
                (
                    "APC",
                    "",
                    None,
                    -0.0,
                    clinvar_claim(
                        variation_id=7,
                        clnsig="Pathogenic|risk_factor",
                        clnrevstat="no_assertion",
                        condition='Syndrome de Lynch 😀 "quoted" \\ back',
                    ),
                ),
            ]
        )
    ]
    variants = [{"vid": f"V{n}", "rsid": None if n else "rs7"} for n in range(3)]
    writer = StoreWriter(con, tmp_path / "work")
    with transaction(con), writer:
        for row in calls:
            writer.add("call", {**row, "run_id": bulk})
            _bound(con, "call", {**row, "run_id": bound})
        for row in statements:
            writer.add("statement", {**row, "run_id": bulk})
            _bound(con, "statement", {**row, "run_id": bound})
        for row in variants:
            writer.add("variant", {**row, "run_id": bulk})
    for table in ("call", "statement"):
        got, want = _all(con, table, bulk), _all(con, table, bound)
        # NaN != NaN: compare the representations.
        assert repr(got) == repr(want), table
        assert len(got) == len(calls if table == "call" else statements)
    assert con.execute(
        "SELECT p_wrong FROM call WHERE run_id = ? AND vid = 'GRCh38:chr1:0:A:G'",
        [bulk],
    ).fetchone() == (1.8e-21,)
    assert con.execute("SELECT vid, rsid FROM variant ORDER BY vid").fetchall() == [
        ("V0", "rs7"),
        ("V1", None),
        ("V2", None),
    ]
    # A claim keeps its text, not a JSON string holding it.
    assert con.execute(
        "SELECT claim->>'condition' FROM statement"
        " WHERE run_id = ? AND statement_id = 'clinvar:2'",
        [bulk],
    ).fetchone() == ('Syndrome de Lynch 😀 "quoted" \\ back',)


def test_a_check_violation_still_raises(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    writer.add("call", _call(run_id, level="certain"))
    with pytest.raises(duckdb.ConstraintException):
        writer.flush()
    assert writer.pending == 0


def test_rows_pass_through_a_private_directory_removed_after(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    con, run_id = _store(tmp_path)
    work = tmp_path / "work"
    writer = StoreWriter(con, work)
    seen: list[tuple[Path, int]] = []
    insert = StoreWriter._insert  # pyright: ignore[reportPrivateUsage]

    def recording(self: StoreWriter, table: str, rows: list, path: Path) -> None:
        seen.append((path.parent, path.parent.stat().st_mode & 0o777))
        insert(self, table, rows, path)

    monkeypatch.setattr(StoreWriter, "_insert", recording)
    writer.add("variant", {"run_id": run_id, "vid": VID})
    writer.add("call", _call(run_id))
    writer.flush()
    assert {(folder.parent, mode) for folder, mode in seen} == {(work, 0o700)}
    assert list(work.iterdir()) == []


def test_a_failed_flush_leaves_no_rows_and_no_directory(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    work = tmp_path / "work"
    writer = StoreWriter(con, work)
    with pytest.raises(duckdb.ConstraintException), transaction(con), writer:
        start_run(con, {}, {})
        writer.add("variant", {"run_id": run_id, "vid": VID})
        writer.add("call", _call(run_id))
        writer.flush()  # flushed rows, rolled back with the transaction
        writer.add("statement", _statement(run_id, "clinvar:1"))
        writer.add("call", _call(run_id, vid="GRCh38:chr1:2:A:G", row_type="Pass"))
    for table in ("variant", "call", "statement"):
        assert con.execute(f"SELECT count(*) FROM {table}").fetchone() == (0,)
    assert con.execute("SELECT count(*) FROM run").fetchone() == (1,)
    assert list(work.iterdir()) == []


def test_a_block_that_raises_discards_its_rows(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    writer.add("variant", {"run_id": run_id, "vid": "GRCh38:chr1:1:A:G"})
    with pytest.raises(RuntimeError), writer:
        # Entering flushed the row buffered before the block.
        assert writer.pending == 0
        writer.add("variant", {"run_id": run_id, "vid": VID})
        raise RuntimeError("abort")
    assert writer.pending == 0
    writer.flush()
    assert con.execute("SELECT vid FROM variant").fetchall() == [("GRCh38:chr1:1:A:G",)]


def test_an_unknown_table_or_column_is_refused(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    with pytest.raises(ValueError, match="not a buffered table"):
        writer.add("run", {"run_id": run_id})
    with pytest.raises(ValueError, match="no columns"):
        writer.add("call", _call(run_id, genotype="0/1"))


def test_a_panel_sized_run_writes_in_seconds(tmp_path: Path) -> None:
    # Row by row with ON CONFLICT, 20,000 calls took 56 s and 7 GiB in one transaction.
    n = 50_000
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    started = time.monotonic()
    with transaction(con), writer:
        for i in range(n):
            vid = f"GRCh38:chr1:{i + 1}:A:G"
            writer.add(
                "variant",
                {"run_id": run_id, "vid": vid, "rsid": f"rs{i}"},
            )
            writer.add("call", _call(run_id, vid=vid, p_wrong=i * 1e-9))
            writer.add("statement", _statement(run_id, f"clinvar:{i}", vid=vid))
    elapsed = time.monotonic() - started
    for table in ("variant", "call", "statement"):
        assert con.execute(f"SELECT count(*) FROM {table}").fetchone() == (n,)
    assert elapsed < 30, f"{elapsed:.1f} s"
    assert writer_module.MAX_ROWS < 3 * n  # the buffer flushed on its own on the way


def test_sweep_rows_removes_left_row_files_and_names_what_it_could_not(
    tmp_path: Path,
) -> None:
    gone = tmp_path / ".store-rows-a"
    gone.mkdir()
    (gone / "call.jsonl").write_text("{}\n")
    stuck = tmp_path / ".store-rows-b" / "inner"
    stuck.mkdir(parents=True)
    (stuck / "call.jsonl").write_text("{}\n")
    kept = tmp_path / "other"
    kept.mkdir()
    stuck.chmod(0o500)  # its file cannot be unlinked
    try:
        failed = sweep_rows(tmp_path)
    finally:
        stuck.chmod(0o700)
    assert not gone.exists() and kept.exists()
    assert len(failed) == 1 and failed[0].startswith(str(tmp_path / ".store-rows-b"))


def _exclusion(confirm: object) -> str:
    claim: dict[str, object] = {"type": "repeat_expansion", "gene": "HTT"}
    if confirm is not None:
        claim["confirm"] = confirm
    return json.dumps(claim)


def test_statement_exclusion_defaults_false(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path)
    con.execute("INSERT INTO variant (run_id, vid) VALUES (?, ?)", [run_id, VID])
    with StoreWriter(con, tmp_path / "work") as writer:
        writer.add("statement", _statement(run_id, "clinvar:1"))
        writer.add(
            "statement",
            _statement(
                run_id,
                "clinvar:2",
                claim=_exclusion("repeat-primed PCR"),
                exclusion=True,
            ),
        )
    assert con.execute(
        "SELECT statement_id, exclusion FROM statement ORDER BY ALL"
    ).fetchall() == [("clinvar:1", False), ("clinvar:2", True)]


@pytest.mark.parametrize("confirm", [None, "", "  "])
def test_exclusion_needs_a_confirm_test(tmp_path: Path, confirm: object) -> None:
    con, run_id = _store(tmp_path)
    writer = StoreWriter(con, tmp_path / "work")
    row = _statement(run_id, "expansionhunter:HTT", claim=_exclusion(confirm))
    with pytest.raises(ValueError, match="expansionhunter:HTT"):
        writer.add("statement", {**row, "exclusion": True})
    assert writer.pending == 0
    writer.flush()
    assert _all(con, "statement", run_id) == []
