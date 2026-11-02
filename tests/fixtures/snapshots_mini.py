"""Run snapshots written as `genokeeper run` writes them: a store with the real schema,
filled with synthetic rows and exported by `export_snapshot` to
`<root>/<run_id>/<table>.parquet`, in the store's own column order."""

import json
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime, timedelta
from importlib.resources import files
from pathlib import Path
from typing import Any

import duckdb

from genokeeper.refresh_log import RefreshRecord, append, publications_path
from genokeeper.store.db import SCHEMA_VERSION, export_snapshot

VID = "GRCh38:chr17:1000:C:T"


def run_id(day: int, month: int = 10, suffix: str = "aaaaaa") -> str:
    """A run id on 2026-<month>-<day> at noon UTC."""
    return f"2026{month:02d}{day:02d}T120000Z-{suffix}"


def started(rid: str) -> datetime:
    """The start time a run id names."""
    return datetime.strptime(rid[:16], "%Y%m%dT%H%M%SZ")


def file(sha: str) -> dict[str, str]:
    """A file input as `run.inputs` records one."""
    return {"path": f"/lab/{sha}", "sha256": sha}


def statement(
    sid: str,
    *,
    kind: str = "finding",
    vid: str | None = VID,
    gene: str | None = "BRCA1",
    rsid: str | None = None,
    validity: str | None = "strong",
    source_level: str | None = "Pathogenic (2★)",
    reason: str | None = None,
    effect: dict[str, Any] | None = None,
    actionable: bool = False,
    topic: str | None = "risks",
    source: str = "clinvar",
    source_version: str = "20261011",
    claim: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """A statement row (with `rsid` for its variant row)."""
    return {
        "statement_id": sid,
        "kind": kind,
        "vid": vid,
        "gene": gene,
        "rsid": rsid,
        "claim": claim
        or {
            "schema": "genomelab.claim/1",
            "type": "clinvar_germline",
            "classification": "Pathogenic",
            "condition": "Hereditary breast ovarian cancer",
        },
        "validity": validity,
        "source_level": source_level,
        "reason": reason,
        "effect": effect,
        "actionable": actionable,
        "topic": topic,
        "source": source,
        "source_id": sid.split(":", 1)[1],
        "source_version": source_version,
    }


def call(
    vid: str = VID,
    *,
    dosage: str = "1",
    level: str = "high",
    warnings: Sequence[str] = (),
) -> dict[str, Any]:
    """A call row."""
    unknown = dosage == "unknown"
    return {
        "vid": vid,
        "gt": {"unknown": "./.", "0": "0/0", "2": "1/1"}.get(dosage, "0/1"),
        "dosage": dosage,
        "p_wrong": None if unknown else 1e-4,
        "level": level,
        "gq": 40,
        "row_type": "NoCall" if unknown else "PASS",
        "flags": [],
        "warnings": list(warnings),
    }


def gene_call(gene: str, *, level: str | None = "medium") -> dict[str, Any]:
    """A gene call row."""
    return {
        "gene": gene,
        "diplotype": "*1/*2",
        "phenotype": "Intermediate Metabolizer",
        "activity_score": None,
        "callers": {"PharmCAT": "*1/*2"},
        "agreed": None,
        "level": level,
        "warnings": [],
    }


def _insert(con: duckdb.DuckDBPyConnection, table: str, row: Mapping[str, Any]) -> None:
    con.execute(
        f"INSERT INTO {table} ({', '.join(row)}) VALUES ({', '.join('?' for _ in row)})",
        list(row.values()),
    )


def fill(
    con: duckdb.DuckDBPyConnection,
    rid: str,
    started_at: datetime,
    tools: Mapping[str, Any],
    inputs: Mapping[str, Any],
    statements: Sequence[Mapping[str, Any]],
    calls: Sequence[Mapping[str, Any]] = (),
    gene_calls: Sequence[Mapping[str, Any]] = (),
) -> None:
    """Write one run's rows into a store `con` opened with the store's schema."""
    _insert(
        con,
        "run",
        {
            "run_id": rid,
            "started_at": started_at,
            "inputs": json.dumps(inputs),
            "tools": json.dumps(tools),
        },
    )
    variants: dict[str, str | None] = {}
    for s in statements:
        if s["vid"] is not None:
            variants.setdefault(s["vid"], s["rsid"])
    for c in calls:
        variants.setdefault(c["vid"], None)
    for vid, rsid in variants.items():
        _insert(con, "variant", {"run_id": rid, "vid": vid, "rsid": rsid})
    for c in calls:
        _insert(con, "call", {"run_id": rid, **c})
    for s in statements:
        row = {k: v for k, v in s.items() if k != "rsid"}
        row["claim"] = json.dumps(row["claim"])
        row["effect"] = None if row["effect"] is None else json.dumps(row["effect"])
        _insert(con, "statement", {"run_id": rid, **row})
    for g in gene_calls:
        _insert(
            con, "gene_call", {"run_id": rid, **g, "callers": json.dumps(g["callers"])}
        )


def store(path: Path | None = None) -> duckdb.DuckDBPyConnection:
    """A store with the real schema, in memory unless `path` is given."""
    con = duckdb.connect(":memory:" if path is None else str(path))
    con.execute(files("genokeeper.store").joinpath("schema.sql").read_text())
    con.execute("INSERT INTO meta (schema_version) VALUES (?)", [SCHEMA_VERSION])
    return con


def write(
    root: Path,
    rid: str,
    started_at: datetime,
    tools: Mapping[str, Any],
    inputs: Mapping[str, Any],
    statements: Sequence[Mapping[str, Any]],
    calls: Sequence[Mapping[str, Any]] = (),
    gene_calls: Sequence[Mapping[str, Any]] = (),
) -> Path:
    """Export one run's snapshot to `<root>/<rid>/` (`started(rid)` is the start time its
    id names)."""
    con = store()
    try:
        fill(con, rid, started_at, tools, inputs, statements, calls, gene_calls)
        return export_snapshot(con, rid, root)
    finally:
        con.close()


def refreshed(reports: Path, lab: Path, rid: str, pubs: Mapping[str, Any]) -> None:
    """A complete refresh of the run `rid`, an hour after it started, that wrote `pubs` as
    its publications."""
    at = started(rid).replace(tzinfo=UTC) + timedelta(hours=1)
    append(lab, RefreshRecord(at, at, "complete", (), None, None, rid))
    path = publications_path(reports, at)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(pubs))
