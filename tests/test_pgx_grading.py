"""PharmCAT's gene results as graded gene-level statements, the outside-call file PharmCAT
reads, and CPIC's gene-drug levels (track design, slice 1b, "The store" and "Gene-call
reliability")."""

import json
from dataclasses import replace
from pathlib import Path

import duckdb
import pytest
from fixtures.faf_mini import write_faf_table

from genokeeper.config import load_grading
from genokeeper.evidence.clinvar import grade_clinvar
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.pgx import (
    GENE_UNRESOLVED,
    NO_STRENGTH,
    cpic_actionable,
    grade_gene_result,
    load_cpic_alleles,
    load_cpic_levels,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.pgx.outside import OutsideCall, outside_lines
from genokeeper.pgx.pharmcat import DOSING, DrugRec, GeneResult, parse_report
from genokeeper.query import Filters, query
from genokeeper.sources import RunContext, write_collected
from genokeeper.sources.panel import clinvar_statement
from genokeeper.sources.pgx import (
    G6PD_GAP,
    IFNL3_RETIRED,
    PGX_TOOL,
    UGT1A1_28_UNREAD,
    PgxSource,
)
from genokeeper.store.db import Level, Validity, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

REPORT = Path(__file__).parent / "fixtures" / "pharmcat_report_mini.json"
LEVELS = {("CYP2C19", "clopidogrel"): "A", ("CYP2D6", "codeine"): "A"}


def _result(**values: object) -> GeneResult:
    """CYP2C19 *1/*1, a normal metabolizer, from PharmCAT's matcher, no position missing,
    with CPIC's Strong clopidogrel recommendation; then `values`."""
    base = GeneResult(
        gene="CYP2C19",
        diplotype="*1/*1",
        phenotype="Normal Metabolizer",
        activity_score=None,
        call_source="MATCHER",
        missing_positions=(),
        drugs=(DrugRec("clopidogrel", "CPIC", "Strong", ()),),
    )
    return replace(base, **values)  # type: ignore[arg-type]


def _cyp2d6(**values: object) -> OutsideCall:
    base = OutsideCall(
        gene="CYP2D6",
        diplotype="*1/*4",
        phenotype=None,
        activity_score="1.0",
        callers={"BCyrius": "*1/*4", "Aldy": "*2/*4"},
        agreed=True,
        unresolved_reason=None,
    )
    return replace(base, **values)  # type: ignore[arg-type]


# --- grading -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("classifications", "validity"),
    [
        (["Strong", "Optional"], Validity.STRONG),
        (["Optional", "Moderate", "No recommendation"], Validity.MODERATE),
        (["Optional"], Validity.LIMITED),
    ],
)
def test_validity_from_classification(
    classifications: list[str], validity: Validity
) -> None:
    """The strongest recommendation PharmCAT reports for the gene: Strong, Moderate and
    Optional are strong, moderate and limited."""
    drugs = tuple(
        DrugRec(f"drug{i}", "CPIC", c, ()) for i, c in enumerate(classifications)
    )
    got, _, _, warnings = grade_gene_result(
        _result(drugs=drugs), cpic_levels={}, outside=None
    )
    assert got == validity
    assert NO_STRENGTH not in warnings


def test_no_graded_recommendation_is_limited() -> None:
    drugs = (DrugRec("atomoxetine", "DPWG", "Unspecified", ()),)
    validity, _, _, warnings = grade_gene_result(
        _result(drugs=drugs), cpic_levels={}, outside=None
    )
    assert validity == Validity.LIMITED
    assert NO_STRENGTH in warnings


@pytest.mark.parametrize(
    ("levels", "drugs", "actionable"),
    [
        ({("CYP2C19", "clopidogrel"): "A"}, [("clopidogrel", "CPIC")], True),
        ({("CYP2C19", "clopidogrel"): "B"}, [("clopidogrel", "CPIC")], True),
        ({("CYP2C19", "clopidogrel"): "A/B"}, [("clopidogrel", "CPIC")], True),
        ({("CYP2C19", "clopidogrel"): "C"}, [("clopidogrel", "CPIC")], False),
        # Another gene's level does not count.
        ({("CYP2D6", "clopidogrel"): "A"}, [("clopidogrel", "CPIC")], False),
        ({}, [("clopidogrel", "CPIC")], False),
        # A DPWG guideline does.
        ({}, [("clopidogrel", "DPWG")], True),
        ({}, [("clopidogrel", "FDA label")], False),
    ],
)
def test_actionable_from_cpic_level(
    levels: dict[tuple[str, str], str],
    drugs: list[tuple[str, str]],
    actionable: bool,
) -> None:
    recs = tuple(DrugRec(d, s, "Strong", ()) for d, s in drugs)
    _, got, _, _ = grade_gene_result(
        _result(drugs=recs), cpic_levels=levels, outside=None
    )
    assert got is actionable


@pytest.mark.parametrize(
    ("result", "outside", "level"),
    [
        # PharmCAT alone is one caller.
        (_result(), None, Level.MEDIUM),
        (_result(missing_positions=("chr10:94762706",)), None, Level.LOW),
        # Two outside callers agreeing.
        (_result(gene="CYP2D6", call_source="OUTSIDE"), _cyp2d6(), Level.HIGH),
        # One outside caller.
        (
            _result(gene="CYP2B6", call_source="OUTSIDE"),
            _cyp2d6(gene="CYP2B6", callers={"Aldy": "*1/*6"}, agreed=False),
            Level.MEDIUM,
        ),
        # PharmCAT did not take the outside call.
        (_result(gene="CYP2D6", call_source="MATCHER"), _cyp2d6(), Level.LOW),
        (
            _result(gene="CYP2D6", diplotype=None, call_source="NONE"),
            _cyp2d6(agreed=False, unresolved_reason="activity scores differ"),
            Level.LOW,
        ),
    ],
)
def test_gene_call_reliability(
    result: GeneResult, outside: OutsideCall | None, level: Level
) -> None:
    _, _, got, _ = grade_gene_result(result, cpic_levels={}, outside=outside)
    assert got == level


def test_one_callers_score_is_medium_and_says_so() -> None:
    """Two callers, one score: not agreed, so medium, never high, and its warning says the
    score rests on one caller."""
    note = "one caller's score: Aldy stated no activity score for *4/*68"
    outside = _cyp2d6(diplotype=None, agreed=False, note=note)
    _, _, level, warnings = grade_gene_result(
        _result(gene="CYP2D6", call_source="OUTSIDE"), cpic_levels={}, outside=outside
    )
    assert level is Level.MEDIUM
    assert note in warnings


# --- the source and the store --------------------------------------------------------------


def _context(tmp_path: Path, run_id: str) -> RunContext:
    return RunContext(
        run_id=run_id,
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=True,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
    )


def _store(tmp_path: Path, source: PgxSource) -> tuple[duckdb.DuckDBPyConnection, str]:
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {"pharmcat": "3.4.0"})
    collected = source.collect(_context(tmp_path, run_id))
    with transaction(con):
        write_collected(StoreWriter(con, tmp_path / "work"), con, run_id, [collected])
    return con, run_id


def _statements(con: duckdb.DuckDBPyConnection) -> dict[str, dict[str, object]]:
    cursor = con.execute("SELECT * FROM statement")
    names = [c[0] for c in cursor.description]
    return {row[1]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()}


def test_normal_result_is_still_a_finding(tmp_path: Path) -> None:
    """A normal metabolizer is a result, not an absence: a finding under drugs."""
    source = PgxSource(
        results=[_result()], outside={}, cpic_levels=LEVELS, version="3.4.0"
    )
    con, _ = _store(tmp_path, source)
    row = _statements(con)["pharmcat:CYP2C19"]
    assert (row["kind"], row["topic"], row["vid"], row["gene"]) == (
        "finding",
        "drugs",
        None,
        "CYP2C19",
    )
    assert (row["validity"], row["actionable"], row["source"]) == (
        "strong",
        True,
        "pharmcat",
    )
    assert row["source_level"] == "CPIC Strong"
    gene_call = con.execute(
        "SELECT diplotype, phenotype, level, callers FROM gene_call"
    ).fetchall()
    assert gene_call == [
        ("*1/*1", "Normal Metabolizer", "medium", json.dumps({"PharmCAT": "*1/*1"}))
    ]
    accounting = con.execute(
        "SELECT tool, input_rows, findings, unassessed FROM accounting"
    ).fetchall()
    assert accounting == [(PGX_TOOL, 2, 1, 1)]  # the G6PD gap


def test_unresolved_gene_is_gene_level_unassessed(tmp_path: Path) -> None:
    """Outside callers that disagree leave the gene unresolved: no line goes to PharmCAT,
    and the store holds a gene-level unassessed statement naming the gene, with both calls
    in its gene call."""
    unresolved = _cyp2d6(
        diplotype=None,
        activity_score=None,
        agreed=False,
        unresolved_reason="BCyrius and Aldy differ on the activity score",
    )
    assert outside_lines([unresolved]) == ""
    result = _result(
        gene="CYP2D6",
        diplotype=None,
        phenotype=None,
        call_source="NONE",
        drugs=(DrugRec("codeine", "CPIC", "Moderate", ()),),
    )
    source = PgxSource(
        results=[result],
        outside={"CYP2D6": unresolved},
        cpic_levels=LEVELS,
        version="3.4.0",
    )
    con, run_id = _store(tmp_path, source)
    row = _statements(con)["pharmcat:CYP2D6"]
    assert (row["kind"], row["vid"], row["gene"], row["topic"]) == (
        "unassessed",
        None,
        "CYP2D6",
        "drugs",
    )
    level, callers, warnings = con.execute(
        "SELECT level, callers, warnings FROM gene_call WHERE gene = 'CYP2D6'"
    ).fetchone() or (None, None, None)
    assert level == "low"
    assert json.loads(str(callers)) == {"BCyrius": "*1/*4", "Aldy": "*2/*4"}
    reason = f"{GENE_UNRESOLVED}: BCyrius and Aldy differ on the activity score"
    assert reason in list(warnings or [])
    result_json = query(con, run_id, Filters(), list_unassessed=True).to_json()
    assert result_json["unassessed"]["by_gene"]["CYP2D6"]["reasons"] == {reason: 1}
    (site,) = [s for s in result_json["unassessed_sites"] if s["gene"] == "CYP2D6"]
    assert site["gene_call"]["callers"] == {"BCyrius": "*1/*4", "Aldy": "*2/*4"}
    assert site["marker"] == {"rsid": None, "grch38": None}


def test_pharmcat_no_call_is_unassessed(tmp_path: Path) -> None:
    result = _result(gene="HLA-A", diplotype=None, phenotype=None, call_source="NONE")
    con, _ = _store(
        tmp_path,
        PgxSource(results=[result], outside={}, cpic_levels={}, version="3.4.0"),
    )
    assert _statements(con)["pharmcat:HLA-A"]["kind"] == "unassessed"


def test_stored_claim_has_no_dosing_text(tmp_path: Path) -> None:
    """The claims written from PharmCAT's report hold no recommendation text and no word of
    a dose: only the gene's call and each drug's name, source and classification."""
    report = json.loads(REPORT.read_text())
    recommendations: list[str] = []
    implications: list[str] = []
    for drugs in report["drugs"].values():
        for drug in drugs.values():
            for guideline in drug["guidelines"]:
                for annotation in guideline["annotations"]:
                    recommendations.append(annotation["drugRecommendation"])
                    implications.extend(annotation["implications"])
    assert recommendations and any(DOSING.search(i) for i in implications)
    source = PgxSource(
        results=parse_report(REPORT), outside={}, cpic_levels=LEVELS, version="3.4.0"
    )
    con, _ = _store(tmp_path, source)
    claims = [str(c) for (c,) in con.execute("SELECT claim FROM statement").fetchall()]
    assert (
        len(claims) == 7
    )  # five genes, F2 (an unannotated gene call) and the G6PD gap
    for claim in claims:
        assert not DOSING.search(claim), claim
        assert not [t for t in recommendations + implications if t in claim]
        parsed = json.loads(claim)
        if parsed["type"] == "pharmcat_gap":
            assert parsed == {
                "schema": parsed["schema"],
                "type": "pharmcat_gap",
                "gene": "G6PD",
            }
            continue
        assert parsed["type"] == "pharmcat_gene"
        assert all(
            set(d) == {"drug", "source", "classification"} for d in parsed["drugs"]
        )
    reasons = [r for (r,) in con.execute("SELECT reason FROM statement").fetchall()]
    assert not [r for r in reasons if r is not None and DOSING.search(r)]


# --- the outside-call file ------------------------------------------------------------------


def test_outside_tsv_format() -> None:
    """PharmCAT's four columns (gene, diplotype, phenotype, activity score), tab-separated;
    an empty column keeps its tab, trailing empty columns are left off, one line per gene,
    and an unresolved gene has none."""
    calls = [
        _cyp2d6(diplotype="*1/*4", activity_score="1.0"),
        OutsideCall(
            gene="HLA-B",
            diplotype=None,
            phenotype="*57:01 positive",
            activity_score=None,
            callers={"OptiType": "*57:01/*08:01", "T1K": "*57:01/*08:01"},
            agreed=True,
            unresolved_reason=None,
        ),
        OutsideCall(
            gene="CYP2C9",
            diplotype=None,
            phenotype=None,
            activity_score="1.5",
            callers={"x": "y"},
            agreed=True,
            unresolved_reason=None,
        ),
        _cyp2d6(gene="CYP2B6", diplotype="*1/*6", activity_score=None),
        _cyp2d6(gene="HLA-A", unresolved_reason="OptiType and T1K differ on *31:01"),
    ]
    assert outside_lines(calls) == (
        "CYP2D6\t*1/*4\t\t1.0\n"
        "HLA-B\t\t*57:01 positive\n"
        "CYP2C9\t\t\t1.5\n"
        "CYP2B6\t*1/*6\n"
    )
    with pytest.raises(ValueError, match="CYP2D6"):
        outside_lines([_cyp2d6(), _cyp2d6()])


def test_phenotype_left_off_beside_a_diplotype() -> None:
    """Given a diplotype and a phenotype, PharmCAT looks recommendations up by the
    phenotype as written, and found none for HLA-B's `*15:02 negative; *57:01 positive;
    *58:01 negative` beside `*44:02/*57:01` (its run gave an empty lookup key), while the
    diplotype alone gave it each allele's status: the line carries the diplotype only."""
    hla = OutsideCall(
        gene="HLA-B",
        diplotype="*44:02/*57:01",
        phenotype="*15:02 negative; *57:01 positive; *58:01 negative",
        activity_score=None,
        callers={"OptiType": "*44:02/*57:01", "T1K": "*44:02/*57:01"},
        agreed=True,
        unresolved_reason=None,
    )
    assert outside_lines([hla]) == "HLA-B\t*44:02/*57:01\n"
    assert outside_lines([_cyp2d6(phenotype="Intermediate Metabolizer")]) == (
        "CYP2D6\t*1/*4\t\t1.0\n"
    )


@pytest.mark.parametrize(
    ("diplotype", "written"),
    [
        ("*1x3/*4", "*1x≥3/*4"),
        ("*2x4/*1x3", "*2x≥3/*1x≥3"),
        ("*1x≥3/*4", "*1x≥3/*4"),
        ("*1x2/*4", "*1x2/*4"),  # two copies keep their count
        ("*41x5/*4", "*41x≥3/*4"),
    ],
)
def test_three_copies_written(diplotype: str, written: str) -> None:
    """Three or more copies are CPIC's and PharmCAT's `x≥3` bin (PharmCAT 3.4.0 reads
    `*1≥3/*4` and `*1x3/*4` as undocumented alleles, Indeterminate), and the line carries
    its score: PharmCAT phenotypes `*1x≥3/*4` with 3.0 Ultrarapid."""
    line = outside_lines([_cyp2d6(diplotype=diplotype, activity_score="3.0")])
    assert line == f"CYP2D6\t{written}\t\t3.0\n"


# --- CPIC's levels ---------------------------------------------------------------------------


def test_cpic_levels_read_from_the_pair_list(tmp_path: Path) -> None:
    """CPIC's API rows, by their keys; a page, an error object or no rows is refused."""
    path = tmp_path / "cpic_pairs.json"
    rows = [
        {"genesymbol": "CYP2C19", "drugname": "clopidogrel", "cpiclevel": "A"},
        {"genesymbol": "CYP2D6", "drugname": "Codeine", "cpiclevel": "A"},
        {"genesymbol": "CYP2C9", "drugname": "celecoxib", "cpiclevel": "A/B"},
    ]
    path.write_text(json.dumps(rows))
    assert load_cpic_levels(path) == {
        ("CYP2C19", "clopidogrel"): "A",
        ("CYP2D6", "codeine"): "A",
        ("CYP2C9", "celecoxib"): "A/B",
    }
    alleles = tmp_path / "MT-RNR1_alleles.json"
    alleles.write_text(json.dumps([{"name": "m.1555A>G"}, {"name": "Reference"}]))
    assert load_cpic_alleles(alleles) == ["m.1555A>G", "Reference"]
    for bad in ('{"code": "PGRST", "message": "x"}', "[]", '[{"name": "x"}]'):
        path.write_text(bad)
        with pytest.raises(ValueError, match="cpiclevel"):
            load_cpic_levels(path)


def test_drug_response_actionable_from_the_drug_named() -> None:
    """A ClinVar drug_response claim is actionable when its condition names a drug whose
    pair with the gene is CPIC level A or B."""
    assert cpic_actionable("CYP2C19", "Clopidogrel response - Efficacy", LEVELS)
    assert not cpic_actionable("CYP2C19", "Warfarin response", LEVELS)
    assert not cpic_actionable("CYP2D6", "Clopidogrel response", LEVELS)
    assert not cpic_actionable(
        "CYP2C9", "celecoxib response", {("CYP2C9", "celecoxib"): "C"}
    )
    lists = load_gene_lists()

    def statement(clnsig: str, levels: dict[tuple[str, str], str]) -> object:
        grade = grade_clinvar(
            clnsig=clnsig,
            clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
            gene="CYP2C19",
            gene_validity=None,
            lists=lists,
            faf_grpmax=None,
            dosage="1",
            grading=load_grading(),
        )
        row = clinvar_statement(
            "r",
            variation_id=1,
            vid="GRCh38:chr10:94761900:C:T",
            gene="CYP2C19",
            clnsig=clnsig,
            clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
            condition="Clopidogrel response - Efficacy",
            consequence="",
            grade=grade,
            source_version="20261011",
            lists=lists,
            cpic_levels=levels,
        )
        return row["actionable"]

    assert statement("drug_response", LEVELS) is True
    assert statement("drug_response", {}) is False
    assert (
        statement("Pathogenic", LEVELS) is False
    )  # CPIC levels are for drug responses


# --- slice 3f: UGT1A1*28, IFNL3 and the G6PD gap --------------------------------------------


def test_untrusted_ugt1a1_warned(tmp_path: Path) -> None:
    """rs3064744 not trusted (Review Focus 4, the source half): UGT1A1's gene call carries
    the warning rather than a silent missing position; no other gene does, and a trusted
    run carries none."""
    ugt1a1 = _result(gene="UGT1A1", diplotype="*1/*1", drugs=())

    def warnings(unread: bool) -> dict[str, list[str]]:
        source = PgxSource(
            results=[ugt1a1, _result()],
            outside={},
            cpic_levels=LEVELS,
            version="3.4.0",
            ugt1a1_unread=unread,
        )
        collected = source.collect(_context(tmp_path, "r"))
        return {
            str(c["gene"]): list(c["warnings"])  # type: ignore[arg-type]
            for c in collected.extra["gene_call"]
        }

    unread = warnings(True)
    assert UGT1A1_28_UNREAD in unread["UGT1A1"]
    assert UGT1A1_28_UNREAD not in unread["CYP2C19"]
    assert UGT1A1_28_UNREAD not in warnings(False)["UGT1A1"]
    con, _ = _store(
        tmp_path,
        PgxSource(
            results=[ugt1a1],
            outside={},
            cpic_levels=LEVELS,
            version="3.4.0",
            ugt1a1_unread=True,
        ),
    )
    (stored,) = con.execute(
        "SELECT warnings FROM gene_call WHERE gene = 'UGT1A1'"
    ).fetchone() or (None,)
    assert UGT1A1_28_UNREAD in list(stored or [])


def test_ifnl3_retired_and_g6pd_gap(tmp_path: Path) -> None:
    """IFNL3's result is a non-finding retired by CPIC; one unassessed G6PD-gap statement
    every run names what PharmCAT does not call."""
    ifnl3 = _result(
        gene="IFNL3",
        diplotype="rs12979860 reference (C)/rs12979860 reference (C)",
        phenotype=None,
        drugs=(DrugRec("peginterferon alfa-2a", "CPIC", "Strong", ()),),
    )
    source = PgxSource(
        results=[ifnl3, _result()], outside={}, cpic_levels=LEVELS, version="3.4.0"
    )
    con, _ = _store(tmp_path, source)
    rows = _statements(con)
    row = rows["pharmcat:IFNL3"]
    assert (row["kind"], row["reason"], row["actionable"]) == (
        "non_finding",
        IFNL3_RETIRED,
        False,
    )
    gap = rows["pharmcat:G6PD-gap"]
    assert (gap["kind"], gap["gene"], gap["topic"], gap["reason"], gap["source"]) == (
        "unassessed",
        "G6PD",
        "drugs",
        G6PD_GAP,
        PGX_TOOL,
    )
    assert (gap["actionable"], gap["exclusion"]) == (False, False)
    accounting = con.execute(
        "SELECT input_rows, findings, non_findings, unassessed FROM accounting"
    ).fetchall()
    assert accounting == [(3, 1, 1, 1)]
    # With no gene result at all, the gap is still there.
    empty = PgxSource(results=[], outside={}, cpic_levels={}, version="3.4.0")
    (only,) = empty.collect(_context(tmp_path, "r")).rows
    assert only.statement["statement_id"] == "pharmcat:G6PD-gap"
