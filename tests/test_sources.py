"""The source seam: each source collects its rows, the promotion is a pure step across them,
and the run writes them with one accounting row per tool."""

import copy
import json
from dataclasses import dataclass, replace
from pathlib import Path

import duckdb
import pytest
from fixtures.clingen_mini import write_gene_validity
from fixtures.faf_mini import write_faf_table
from fixtures.promotion_mini import INHERITANCE, haploid, promotion_rows

from genokeeper import pipeline
from genokeeper.config import load_grading
from genokeeper.evidence.clingen import gene_level, load_curations
from genokeeper.evidence.clinvar import BA1_REASON, NOT_MATCHED_REASON
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim, ClaimSet
from genokeeper.sources import Collected, Row, RunContext, write_collected
from genokeeper.sources.panel import ClaimsSource, clinvar_claim
from genokeeper.sources.promote import promote_compound_hets
from genokeeper.store.db import Level, RowType, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

REASON = "two heterozygous variants; phase unknown"
# What `pipeline._promote_compound_hets` (the store-side promotion this seam replaced) made
# of `promotion_rows`, run once on it before any code moved: statement -> (kind, topic,
# reason).
FROZEN = {
    "clinvar:1": ("finding", "risks", REASON),
    "clinvar:2": ("finding", "risks", f"a reason of its own; {REASON}"),
    "clinvar:3": ("unassessed", "risks", None),
    "clinvar:4": ("non_finding", None, None),
    "clinvar:5": ("non_finding", None, None),
    "clinvar:6": ("non_finding", None, None),
    "clinvar:7": ("non_finding", None, None),
    "clinvar:8": ("non_finding", None, None),
    "clinvar:9": ("non_finding", None, None),
    "clinvar:10": ("non_finding", None, None),
    "clinvar:11": ("non_finding", None, None),
    "clinvar:12": ("finding", "risks", REASON),
    "clinvar:13": ("finding", "risks", REASON),
    "clinvar:14": ("non_finding", None, None),
    "clinvar:15": ("non_finding", None, None),
    "clinvar:16": ("non_finding", None, None),
    "clinvar:17": ("non_finding", None, None),
    "clinvar:18": ("non_finding", None, None),
    "clinvar:19": ("non_finding", None, None),
    "clinvar:20": ("finding", "risks", REASON),
    "clinvar:21": ("finding", "risks", REASON),
    "clinvar:22": ("finding", "risks", REASON),
    "clinvar:23": ("non_finding", None, None),
    "clinvar:24": ("non_finding", None, None),
    "clinvar:25": ("non_finding", None, None),
    "clinvar:26": ("finding", "risks", REASON),
    "clinvar:27": ("finding", "risks", REASON),
}


def _store(tmp_path: Path) -> tuple[duckdb.DuckDBPyConnection, str, StoreWriter]:
    con = connect(tmp_path / "findings.duckdb")
    return con, start_run(con, {}, {}), StoreWriter(con, tmp_path / "work")


def _ctx(tmp_path: Path, run_id: str) -> RunContext:
    return RunContext(
        run_id=run_id,
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance=INHERITANCE,
        male=True,
        haploid=haploid,
        source_version="20261012",
    )


def _collected(run_id: str, dropped: dict[str, int] | None = None) -> list[Collected]:
    """`promotion_rows` as one `Collected` per tool, in the order a run writes them;
    `dropped` is every tool's, its count added to the input rows."""
    by_tool: dict[str, list[Row]] = {}
    for tool, variant, call, statement in promotion_rows(run_id):
        by_tool.setdefault(tool, []).append(Row(variant, call, statement))
    gone = dropped or {}
    return [
        Collected(
            tool=tool,
            rows=rows,
            extra={},
            input_rows=len(rows) + sum(gone.values()),
            dropped=dict(gone),
        )
        for tool, rows in by_tool.items()
    ]


def _accounting(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, tuple]:
    rows = con.execute(
        "SELECT tool, input_rows, findings, watch, non_findings, unassessed, dropped"
        " FROM accounting WHERE run_id = ?",
        [run_id],
    ).fetchall()
    return {tool: (*counts, json.loads(dropped)) for tool, *counts, dropped in rows}


def test_promotion_is_pure() -> None:
    collected = _collected("run")
    before = copy.deepcopy(collected)
    out = promote_compound_hets(collected, load_gene_lists(), INHERITANCE, haploid)
    assert collected == before
    assert [c.tool for c in out] == [c.tool for c in collected]
    got = {
        row.statement["statement_id"]: (
            row.statement["kind"],
            row.statement["topic"],
            row.statement["reason"],
        )
        for c in out
        for row in c.rows
    }
    assert got == FROZEN
    # Only statements change: every variant and call is the one collected.
    assert [(r.variant, r.call) for c in out for r in c.rows] == [
        (r.variant, r.call) for c in before for r in c.rows
    ]


def test_accounting_from_rows(tmp_path: Path) -> None:
    con, run_id, writer = _store(tmp_path)
    collected = _collected(run_id, {"not_in_scope": 2})
    promoted = promote_compound_hets(collected, load_gene_lists(), INHERITANCE, haploid)
    with transaction(con), writer:
        write_collected(writer, con, run_id, promoted)
    # The kinds the rows carry after the promotion (FROZEN), plus what each tool dropped.
    assert _accounting(con, run_id) == {
        "first": (26, 7, 0, 16, 1, {"not_in_scope": 2}),
        "second": (5, 2, 0, 1, 0, {"not_in_scope": 2}),
    }
    kinds = dict(
        con.execute(
            "SELECT statement_id, kind FROM statement WHERE run_id = ?", [run_id]
        ).fetchall()
    )
    assert kinds == {sid: kind for sid, (kind, _, _) in FROZEN.items()}


def test_extra_tables_written(tmp_path: Path) -> None:
    # Schema 3 has no table only a source's extras fill (slice 1b's gene_call will be one):
    # a source's extra rows go through the writer as its own rows do.
    con, run_id, writer = _store(tmp_path)
    collected = Collected(
        tool="extras",
        rows=[],
        extra={
            "variant": [
                {"run_id": run_id, "vid": "GRCh38:chr1:1:A:G", "rsid": "rs1"},
                {"run_id": run_id, "vid": "GRCh38:chr1:2:A:G", "rsid": None},
                # The first row at a key wins, as the store has always kept it.
                {"run_id": run_id, "vid": "GRCh38:chr1:1:A:G", "rsid": "rs-later"},
            ],
        },
        input_rows=0,
        dropped={},
    )
    with transaction(con), writer:
        write_collected(writer, con, run_id, [collected])
    assert con.execute(
        "SELECT vid, rsid FROM variant WHERE run_id = ? ORDER BY 1", [run_id]
    ).fetchall() == [("GRCh38:chr1:1:A:G", "rs1"), ("GRCh38:chr1:2:A:G", None)]
    assert _accounting(con, run_id) == {"extras": (0, 0, 0, 0, 0, {})}


@dataclass(frozen=True)
class _OneHet:
    """A source the run has never heard of: one MUTYH het (a carrier finding: MUTYH is a
    Tier 4 gene), and its variant as an extra row."""

    tool: str
    n: int

    def collect(self, ctx: RunContext) -> Collected:
        vid = f"GRCh38:chr1:{1000 + self.n}:A:G"
        call: dict[str, object] = {
            "run_id": ctx.run_id,
            "vid": vid,
            "gt": "0/1",
            "dosage": "1",
            "level": "high",
            "row_type": "PASS",
        }
        statement: dict[str, object] = {
            "run_id": ctx.run_id,
            "statement_id": f"clinvar:{self.n}",
            "kind": "finding",
            "topic": "carrier",
            "vid": vid,
            "gene": "MUTYH",
            "claim": clinvar_claim(
                variation_id=self.n,
                clnsig="Pathogenic",
                clnrevstat="criteria_provided,_single_submitter",
                condition="not provided",
            ),
            "validity": "strong",
            "source_level": "criteria_provided,_single_submitter",
            "source": "clinvar",
            "source_id": str(self.n),
            "source_version": ctx.source_version,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, call, statement)],
            extra={"variant": [{"run_id": ctx.run_id, "vid": vid, "rsid": None}]},
            input_rows=2,
            dropped={"skipped": 1},
        )


def test_third_source_needs_no_run_change(tmp_path: Path) -> None:
    con, run_id, writer = _store(tmp_path)
    sources = [_OneHet("first", 1), _OneHet("second", 2), _OneHet("third", 3)]
    with transaction(con), writer:
        pipeline.store_sources(writer, con, _ctx(tmp_path, run_id), sources)
    # Each source's one het pairs with the others' in MUTYH: the promotion spans sources,
    # and each tool's accounting counts its promoted finding.
    assert _accounting(con, run_id) == {
        tool: (2, 1, 0, 0, 0, {"skipped": 1}) for tool in ("first", "second", "third")
    }
    assert con.execute(
        "SELECT statement_id, kind, reason FROM statement WHERE run_id = ? ORDER BY 1",
        [run_id],
    ).fetchall() == [(f"clinvar:{n}", "finding", REASON) for n in (1, 2, 3)]
    assert con.execute(
        "SELECT count(*) FROM variant WHERE run_id = ?", [run_id]
    ).fetchone() == (3,)


# OTC: Definitive for one disease (so Definitive gene-level), Disputed for the claim's own.
OWN_MONDO = "MONDO:0005001"


@pytest.mark.parametrize(
    ("disease_ids", "faf", "validity", "reason"),
    [
        # The curation of the claim's own disease grades it.
        (f"MONDO:{OWN_MONDO},MedGen:C1", None, "conflicting", "gene validity Disputed"),
        # No curation of its disease: the gene-level one, unmatched.
        ("MedGen:C2", None, "strong", NOT_MATCHED_REASON),
        # Common in gnomAD: BA1 makes it conflicting.
        ("MedGen:C3", 0.2, "conflicting", BA1_REASON),
    ],
)
def test_claims_are_graded_on_their_curation_and_frequency(
    tmp_path: Path,
    disease_ids: str,
    faf: float | None,
    validity: str,
    reason: str,
) -> None:
    validity_csv = write_gene_validity(
        tmp_path / "gene-validity.csv",
        [("OTC", "Definitive", "MONDO:0000077"), ("OTC", "Disputed", OWN_MONDO)],
    )
    curations = load_curations(validity_csv)
    vid = "GRCh38:chrX:38400000:A:G"
    claim = Claim(
        variation_id=1,
        gene="OTC",
        clnsig="Pathogenic",
        clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
        disease_ids=disease_ids,
        condition="Ornithine transcarbamylase deficiency",
        chrom="chrX",
        pos=38_400_000,
        ref="A",
        alt="G",
        vid=vid,
        rsid=None,
    )
    het = Called(
        gt="0/1",
        dosage="1",
        row_type=RowType.PASS,
        gq=50,
        dp=30,
        ad=(15, 15, 0),
        vaf=0.5,
        warnings=(),
    )
    graded = Reliability(p_wrong=1e-5, level=Level.HIGH, flags=(), warnings=())
    ctx = replace(
        _ctx(tmp_path, "20261020T000000Z-abcdef"),
        curations=curations,
        gene_level={gene: gene_level(found) for gene, found in curations.items()},
        faf=FafTable.open(
            write_faf_table(tmp_path / "faf", {} if faf is None else {vid: faf})
        ),
        male=False,
        haploid=lambda chrom, pos: False,
    )
    source = ClaimsSource(
        ClaimSet(claims=[claim], input_rows=1, dropped={}), [(claim, het, graded)]
    )
    (row,) = source.collect(ctx).rows
    assert row.statement["validity"] == validity
    assert str(row.statement["reason"]).startswith(reason)
