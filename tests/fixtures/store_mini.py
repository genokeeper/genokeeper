"""A findings store seeded with SQL inserts: one run, its variants, calls and statements."""

import json
from dataclasses import dataclass, field
from typing import Any

import duckdb

RUN = "20261018T120000Z-aaaaaa"
TOOLS = {"clinvar": "20261011", "deepvariant": "1.10.0"}


@dataclass
class Row:
    """One seeded statement, and (unless gene-level) its variant and call.

    `level=None` makes a gene-level statement (no vid); `call=False` keeps the variant and
    the statement's vid but writes no call row. `kind` is the statement's kind (`finding`,
    `watch`, `non_finding` or `unassessed`), `warnings` its call's; `gt=None` derives the
    call's `gt` from `dosage`. `gene_call` is a `gene_call` row's columns (callers as a
    dict), written for the gene; an id that is not `clinvar:<number>` gets `claim`.
    `reason` is the statement's.
    """

    id: str
    kind: str = "finding"
    gene: str | None = "BRCA1"
    level: str | None = "high"
    call: bool = True
    validity: str = "strong"
    dosage: str = "1"
    gt: str | None = None
    topic: str = "risks"
    actionable: bool = False
    effect: dict[str, Any] | None = None
    pos: int = 1000
    chrom: str = "chr17"
    rsid: str | None = None
    condition: str = "Hereditary breast ovarian cancer"
    classification: str = "Pathogenic"
    flags: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    gene_call: dict[str, Any] | None = None
    claim: dict[str, Any] | None = None
    reason: str | None = None


_GT = {"unknown": "./.", "0": "0/0", "2": "1/1"}


def seed(con: duckdb.DuckDBPyConnection, rows: list[Row], run_id: str = RUN) -> None:
    con.execute(
        "INSERT INTO run (run_id, started_at, inputs, tools) VALUES (?, ?, '{}', ?)"
        " ON CONFLICT DO NOTHING",
        [
            run_id,
            f"{run_id[:4]}-{run_id[4:6]}-{run_id[6:8]} 12:00:00",
            json.dumps(TOOLS),
        ],
    )
    for row in rows:
        vid = None
        if row.level is not None:
            vid = f"GRCh38:{row.chrom}:{row.pos}:C:T"
            con.execute(
                "INSERT INTO variant (run_id, vid, rsid) VALUES (?, ?, ?)"
                " ON CONFLICT DO NOTHING",
                [run_id, vid, row.rsid],
            )
        if vid is not None and row.call:
            con.execute(
                "INSERT INTO call (run_id, vid, gt, dosage, p_wrong, level, gq, dp, row_type,"
                " flags, warnings) VALUES (?, ?, ?, ?, ?, ?, 40, 30, ?, ?, ?)",
                [
                    run_id,
                    vid,
                    row.gt if row.gt is not None else _GT.get(row.dosage, "0/1"),
                    row.dosage,
                    None if row.dosage == "unknown" else 1e-4,
                    row.level,
                    "NoCall" if row.dosage == "unknown" else "PASS",
                    row.flags,
                    row.warnings,
                ],
            )
        if row.gene_call is not None:
            columns = {"run_id": run_id, "gene": row.gene, **row.gene_call}
            columns["callers"] = json.dumps(columns.get("callers"))
            con.execute(
                f"INSERT INTO gene_call ({', '.join(columns)})"
                f" VALUES ({', '.join('?' for _ in columns)})",
                list(columns.values()),
            )
        source_id = row.id.split(":")[1]
        claim = row.claim or {
            "schema": "genomelab.claim/1",
            "type": "clinvar_germline",
            "classification": row.classification,
            "review_status": "criteria_provided,_multiple_submitters,_no_conflicts",
            "condition": row.condition,
            "variation_id": int(source_id),
        }
        con.execute(
            "INSERT INTO statement (run_id, statement_id, kind, vid, gene, claim, validity,"
            " source_level, reason, effect, actionable, topic, source, source_id,"
            " source_version) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
            " 'clinvar', ?, '20261011')",
            [
                run_id,
                row.id,
                row.kind,
                vid,
                row.gene,
                json.dumps(claim),
                row.validity,
                f"{row.classification} (2★)",
                row.reason,
                None if row.effect is None else json.dumps(row.effect),
                row.actionable,
                row.topic,
                source_id,
            ],
        )


def seed_coverage(
    con: duckdb.DuckDBPyConnection,
    gene: str,
    callable_share: float,
    *,
    owner: str | None = None,
    run_id: str = RUN,
) -> None:
    """One `coverage` row of the run: reliable at the shipped 0.95 or more callable."""
    con.execute(
        "INSERT INTO coverage (run_id, gene, bases, callable, mean_depth, difficult,"
        " outside_benchmark, owner, reliable) VALUES (?, ?, 1200, ?, 34.5, 0.1, 0.2, ?, ?)",
        [run_id, gene, callable_share, owner, callable_share >= 0.95],
    )
