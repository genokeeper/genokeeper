import json
import re
from dataclasses import replace
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures import snapshots_mini as sm
from fixtures.store_mini import RUN, TOOLS, Row, seed, seed_coverage

from genokeeper.carrier import genotype
from genokeeper.cli import main
from genokeeper.insertions import INSERTIONS_TOOL
from genokeeper.query import HAPLOID_DOSAGE_NOTE, Filters, RunNotFound, query, show
from genokeeper.sources.cnv import CNV_TOOL
from genokeeper.sources.deletions import DELETIONS_TOOL
from genokeeper.sources.paralogs import LPA_TOOL, PARALOGS_TOOL
from genokeeper.sources.predicted import PREDICTED_TOOL
from genokeeper.sources.repeats import REPEATS_TOOL
from genokeeper.sources.smn import SMN_TOOL
from genokeeper.sources.yleaf import YLEAF_TOOL
from genokeeper.store.db import SCHEMA_VERSION, Level, Topic, Validity, connect

OLD_RUN = "20261001T120000Z-bbbbbb"


@pytest.fixture
def con(tmp_path: Path) -> duckdb.DuckDBPyConnection:
    return connect(tmp_path / "genomelab.duckdb")


def ids(con: duckdb.DuckDBPyConnection, f: Filters | None = None) -> list[str]:
    return [row["id"] for row in query(con, None, f or Filters()).findings]


def test_default_filters_hide_low_and_limited(con: duckdb.DuckDBPyConnection) -> None:
    seed(
        con,
        [
            Row("clinvar:1", level="high", validity="strong", pos=1001),
            Row("clinvar:2", level="medium", validity="moderate", pos=1002),
            Row("clinvar:3", level="low", validity="strong", pos=1003),
            Row("clinvar:4", level="high", validity="limited", pos=1004),
            Row("clinvar:5", level="low", validity="limited", pos=1005),
            Row("clinvar:6", level="high", validity="conflicting", pos=1006),
        ],
    )
    result = query(con, None, Filters())
    assert sorted(f["id"] for f in result.findings) == ["clinvar:1", "clinvar:2"]
    # clinvar:5 fails both filters and counts once, under the first (reliability).
    assert result.hidden == {
        "reliability_below:medium": 2,
        "validity_below:moderate": 2,
    }


def test_lowest_thresholds_hide_nothing(con: duckdb.DuckDBPyConnection) -> None:
    seed(
        con,
        [
            Row("clinvar:3", level="low", validity="strong", pos=1003),
            Row("clinvar:6", level="high", validity="conflicting", pos=1006),
        ],
    )
    result = query(
        con, None, Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    )
    assert sorted(f["id"] for f in result.findings) == ["clinvar:3", "clinvar:6"]
    assert result.hidden == {}


def test_high_threshold_names_the_threshold(con: duckdb.DuckDBPyConnection) -> None:
    seed(con, [Row("clinvar:2", level="medium", validity="moderate")])
    result = query(
        con, None, Filters(min_reliability=Level.HIGH, min_validity=Validity.STRONG)
    )
    assert result.findings == []
    assert result.hidden == {"reliability_below:high": 1}


NOCALL = "genotype unknown: nocall"
NO_READS = "genotype unknown: no_reads"


def _unassessed(id: str, pos: int, **values: Any) -> Row:
    """An unassessed claim at a NoCall site."""
    row = Row(
        id, kind="unassessed", level="low", dosage="unknown", pos=pos, warnings=[NOCALL]
    )
    return replace(row, **values)


def test_an_unassessed_claim_is_counted_never_a_finding(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            _unassessed("clinvar:7", 1007),
            _unassessed("clinvar:8", 1008, validity="limited"),
            Row("clinvar:9", pos=1009),
        ],
    )
    result = query(con, None, Filters(min_reliability=Level.HIGH))
    # Never a finding, whatever the reliability threshold: it is counted per gene.
    assert [f["id"] for f in result.findings] == ["clinvar:9"]
    assert result.findings[0]["kind"] == "finding"
    assert result.unassessed == {
        "total": 1,
        "by_gene": {
            "BRCA1": {
                "sites": 1,
                "withheld": 0,
                "topics": {"risks": 1},
                "reasons": {NOCALL: 1},
            }
        },
    }
    # The unassessed row still faces the validity filter, counted under its own prefix.
    assert result.hidden == {"unassessed:validity_below:moderate": 1}
    assert result.unassessed_sites is None
    assert "unassessed_sites" not in result.to_json()


def test_the_unassessed_summary_counts_genes_topics_and_reasons(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            _unassessed("clinvar:11", 1011),
            _unassessed(
                "clinvar:12", 1012, warnings=["het on a haploid contig", NOCALL]
            ),
            _unassessed("clinvar:13", 1013, warnings=[NO_READS, NOCALL]),
            _unassessed("clinvar:14", 1014, topic="carrier"),
            _unassessed("clinvar:15", 1015, call=False),
            _unassessed("clinvar:16", 1016, gene="CFTR", warnings=[NO_READS]),
        ],
    )
    result = query(con, None, Filters())
    assert result.findings == []
    assert result.unassessed == {
        "total": 6,
        "by_gene": {
            "BRCA1": {
                "sites": 5,
                "withheld": 0,
                "topics": {"risks": 4, "carrier": 1},
                # Each row once, under its first "genotype unknown" warning.
                "reasons": {NOCALL: 3, NO_READS: 1, "no call row": 1},
            },
            "CFTR": {
                "sites": 1,
                "withheld": 0,
                "topics": {"risks": 1},
                "reasons": {NO_READS: 1},
            },
        },
    }
    assert result.hidden == {}


def test_filters_drop_unassessed_claims_into_prefixed_hidden_counts(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            _unassessed("clinvar:21", 1021, validity="limited"),
            _unassessed("clinvar:22", 1022, topic="carrier"),
            _unassessed("clinvar:23", 1023, gene="CFTR"),
            _unassessed("clinvar:24", 1024, rsid="rs24"),
            _unassessed("clinvar:25", 1025, rsid="rs25", condition="x"),
            Row("clinvar:26", pos=1026, level="low"),
            Row("clinvar:27", pos=1027, rsid="rs27", condition="x"),
        ],
    )
    f = Filters(topic=Topic.RISKS, gene="BRCA1", marker="rs25", text="ovarian")
    result = query(con, None, f, list_unassessed=True)
    assert result.findings == []
    assert result.unassessed == {"total": 0, "by_gene": {}}
    assert result.unassessed_sites == []
    # No row disappears uncounted; the findings' reasons come first, each in filter order.
    assert list(result.hidden.items()) == [
        ("reliability_below:medium", 1),
        ("marker", 1),
        ("unassessed:validity_below:moderate", 1),
        ("unassessed:topic", 1),
        ("unassessed:gene", 1),
        ("unassessed:marker", 1),
        ("unassessed:text", 1),
    ]
    assert sum(result.hidden.values()) == 7


def test_list_unassessed_adds_the_sites_in_the_findings_shape_and_order(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            _unassessed("clinvar:32", 1032, validity="moderate"),
            _unassessed("clinvar:31", 1031, actionable=True),
            _unassessed("clinvar:33", 1033),
            Row("clinvar:34", pos=1034),
        ],
    )
    result = query(con, None, Filters(), list_unassessed=True)
    assert result.unassessed_sites is not None
    assert [s["id"] for s in result.unassessed_sites] == [
        "clinvar:31",
        "clinvar:33",
        "clinvar:32",
    ]
    assert {s["kind"] for s in result.unassessed_sites} == {"unassessed"}
    assert [f["kind"] for f in result.findings] == ["finding"]
    assert set(result.unassessed_sites[0]) == set(result.findings[0])
    assert result.unassessed_sites[0]["genotype"] == {
        "gt": "./.",
        "dosage": "unknown",
        "allele": "T",
    }
    assert result.to_json()["unassessed_sites"] == result.unassessed_sites


def test_show_finds_an_unassessed_claim(con: duckdb.DuckDBPyConnection) -> None:
    seed(con, [_unassessed("clinvar:41", 1041), Row("clinvar:42", pos=1042)])
    _run, site = show(con, None, "clinvar:41")
    assert site is not None and site["kind"] == "unassessed"
    _run, finding = show(con, None, "clinvar:42")
    assert finding is not None and finding["kind"] == "finding"


def test_gene_level_finding_skips_reliability_and_has_null_marker(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(con, [Row("clinvar:9", gene="LDLR", level=None)])
    result = query(con, None, Filters(min_reliability=Level.HIGH))
    [finding] = result.findings
    assert finding["marker"] == {"rsid": None, "grch38": None}
    assert finding["genotype"] == {"gt": None, "dosage": None, "allele": None}
    assert set(finding["reliability"].values()) == {None}
    assert finding["gene_call"] is None
    assert result.hidden == {}


PGX_CLAIM = {"schema": "genomelab.claim/1", "type": "pharmcat_gene", "drugs": []}


def _gene_call(level: str, **values: Any) -> dict[str, Any]:
    return {
        "diplotype": "*1/*2",
        "phenotype": "Intermediate Metabolizer",
        "activity_score": 1.5,
        "callers": {"BCyrius": "*1/*2", "Aldy": "*1/*2"},
        "agreed": True,
        "level": level,
        "warnings": [],
    } | values


def test_a_gene_level_row_shows_its_gene_call(con: duckdb.DuckDBPyConnection) -> None:
    """A pharmacogenomic gene result carries its gene call, and the call's reliability
    level faces the reliability filter as a variant's does."""
    seed(
        con,
        [
            Row(
                "pharmcat:CYP2D6",
                gene="CYP2D6",
                level=None,
                topic="drugs",
                gene_call=_gene_call("high"),
                claim=PGX_CLAIM,
            ),
            Row(
                "pharmcat:CYP2C19",
                gene="CYP2C19",
                level=None,
                topic="drugs",
                gene_call=_gene_call("low", warnings=["positions missing: chr10:1"]),
                claim=PGX_CLAIM,
            ),
        ],
    )
    result = query(con, None, Filters())
    [finding] = result.findings
    assert finding["id"] == "pharmcat:CYP2D6"
    assert finding["gene_call"] == {
        "diplotype": "*1/*2",
        "phenotype": "Intermediate Metabolizer",
        "activity_score": 1.5,
        "callers": {"BCyrius": "*1/*2", "Aldy": "*1/*2"},
        "agreed": True,
        "level": "high",
        "warnings": [],
    }
    assert finding["genotype"] == {"gt": None, "dosage": None, "allele": None}
    assert result.hidden == {"reliability_below:medium": 1}
    _run, shown = show(con, None, "pharmcat:CYP2C19")
    assert shown is not None
    assert shown["gene_call"]["warnings"] == ["positions missing: chr10:1"]


def test_a_gene_level_unassessed_claim_counts_why(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """An unresolved gene call is counted under its own reason; a gene-level unassessed
    claim without a gene call under `no gene call`; one with neither marker nor gene call
    (a polygenic score below its overlap, schema 6) under the reason it carries."""
    unresolved = "gene call unresolved: BCyrius and Aldy differ"
    seed(
        con,
        [
            Row(
                "pharmcat:CYP2D6",
                kind="unassessed",
                gene="CYP2D6",
                level=None,
                topic="drugs",
                gene_call=_gene_call("low", diplotype=None, warnings=[unresolved]),
                claim=PGX_CLAIM,
            ),
            Row(
                "pharmcat:HLA-A",
                kind="unassessed",
                gene="HLA-A",
                level=None,
                topic="drugs",
                claim=PGX_CLAIM,
            ),
            Row(
                "pgsc_calc:PGS000018",
                kind="unassessed",
                gene=None,
                level=None,
                claim=PGX_CLAIM,
                reason="score overlap 0.60, below 0.75",
            ),
        ],
    )
    result = query(con, None, Filters())
    assert result.unassessed["by_gene"]["CYP2D6"]["reasons"] == {unresolved: 1}
    assert result.unassessed["by_gene"]["HLA-A"]["reasons"] == {"no gene call": 1}
    assert result.unassessed["by_gene"][""]["reasons"] == {
        "score overlap 0.60, below 0.75": 1
    }


def test_stage3_after_classified_before_predicted(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """The stage-3 sources come after every other source and before the predicted ones,
    whatever their validity, actionability or gene."""
    claim = {"schema": "genomelab.claim/1", "type": "stage3"}
    seed(
        con,
        [
            Row(
                "predicted:1",
                validity="strong",
                actionable=True,
                gene="AAA",
                pos=1001,
                claim=claim,
            ),
            Row(
                "yleaf:haplogroup",
                validity="strong",
                gene=None,
                level=None,
                claim=claim,
            ),
            Row("smn:SMN1", validity="strong", gene="AAA", pos=1003, claim=claim),
            Row(
                "expansionhunter:HTT",
                validity="strong",
                gene="AAA",
                pos=1004,
                claim=claim,
            ),
            Row("clinvar:2", validity="limited", gene="ZZZ", pos=1002),
        ],
    )
    for statement_id, source in (
        ("predicted:1", PREDICTED_TOOL),
        ("yleaf:haplogroup", YLEAF_TOOL),
        ("smn:SMN1", SMN_TOOL),
        ("expansionhunter:HTT", REPEATS_TOOL),
    ):
        con.execute(
            "UPDATE statement SET source = ? WHERE statement_id = ?",
            [source, statement_id],
        )
    found = query(con, None, Filters(min_validity=Validity.LIMITED)).findings
    assert [f["id"] for f in found] == [
        "clinvar:2",
        "expansionhunter:HTT",
        "smn:SMN1",
        "yleaf:haplogroup",
        "predicted:1",
    ]


def test_deletion_after_classified_before_predicted(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """A deletion (slice 3b) is a stage-3 screening result: after a ClinVar-classified
    finding and before a predicted one, whatever their validity or gene."""
    claim = {"schema": "genomelab.claim/1", "type": "deletion"}
    seed(
        con,
        [
            Row("predicted:1", validity="strong", gene="AAA", pos=1001, claim=claim),
            Row("deletions:AAA:chr17:1003-3000", gene="AAA", pos=1003, claim=claim),
            Row("clinvar:2", validity="limited", gene="ZZZ", pos=1002),
        ],
    )
    for statement_id, source in (
        ("predicted:1", PREDICTED_TOOL),
        ("deletions:AAA:chr17:1003-3000", DELETIONS_TOOL),
    ):
        con.execute(
            "UPDATE statement SET source = ? WHERE statement_id = ?",
            [source, statement_id],
        )
    found = query(con, None, Filters(min_validity=Validity.LIMITED)).findings
    assert [f["id"] for f in found] == [
        "clinvar:2",
        "deletions:AAA:chr17:1003-3000",
        "predicted:1",
    ]


def test_3d_after_classified_before_predicted(con: duckdb.DuckDBPyConnection) -> None:
    """Slice 3d's sources (depth copy number, insertions, the paralog loci, Lp(a)) are
    stage-3 results: after a ClinVar-classified finding and before a predicted one,
    whatever their validity or gene."""
    claim = {"schema": "genomelab.claim/1", "type": "stage3"}
    sources = {
        "cnv:AAA:chr17:1003-3000": CNV_TOOL,
        "insertions:AAA:chr17:1004-1005": INSERTIONS_TOOL,
        "paralogs:STRC": PARALOGS_TOOL,
        "lpa:LPA": LPA_TOOL,
        "predicted:1": PREDICTED_TOOL,
    }
    seed(
        con,
        [
            Row("predicted:1", validity="strong", gene="AAA", pos=1001, claim=claim),
            Row("cnv:AAA:chr17:1003-3000", gene="AAA", pos=1003, claim=claim),
            Row("insertions:AAA:chr17:1004-1005", gene="AAA", pos=1004, claim=claim),
            Row("paralogs:STRC", gene="AAA", level=None, claim=claim),
            Row("lpa:LPA", gene="AAA", level=None, claim=claim),
            Row("clinvar:2", validity="limited", gene="ZZZ", pos=1002),
        ],
    )
    for statement_id, source in sources.items():
        con.execute(
            "UPDATE statement SET source = ? WHERE statement_id = ?",
            [source, statement_id],
        )
    found = query(con, None, Filters(min_validity=Validity.LIMITED)).findings
    assert [f["id"] for f in found] == [
        "clinvar:2",
        "cnv:AAA:chr17:1003-3000",
        "insertions:AAA:chr17:1004-1005",
        "lpa:LPA",
        "paralogs:STRC",
        "predicted:1",
    ]


def test_ordering_actionable_first(con: duckdb.DuckDBPyConnection) -> None:
    def rel(value: float) -> dict[str, Any]:
        return {
            "relative": {"measure": "OR", "value": value, "ci": [1.0, 9.0]},
            "absolute": None,
        }

    def absolute(risk: float) -> dict[str, Any]:
        return {
            "relative": None,
            "absolute": {"risk": risk, "baseline": 0.01, "population": "EUR"},
        }

    seed(
        con,
        [
            Row("clinvar:11", gene="ZZZ", validity="moderate", pos=1011),
            Row("clinvar:12", gene="AAA", validity="moderate", pos=1012),
            Row("clinvar:13", gene="MMM", validity="strong", pos=1013),
            Row(
                "clinvar:14", gene="ZZZ", validity="moderate", pos=1014, effect=rel(3.0)
            ),
            Row(
                "clinvar:15", gene="ZZZ", validity="moderate", pos=1015, effect=rel(1.5)
            ),
            Row(
                "clinvar:16",
                gene="ZZZ",
                validity="moderate",
                pos=1016,
                effect=absolute(0.2),
            ),
            Row(
                "clinvar:17", gene="ZZZ", validity="moderate", pos=1017, actionable=True
            ),
            Row(
                "clinvar:18",
                gene="ZZZ",
                validity="moderate",
                pos=1018,
                effect=absolute(0.5),
            ),
        ],
    )
    assert ids(con) == [
        "clinvar:17",  # actionable first, whatever its validity and effect
        "clinvar:13",  # then validity: strong before moderate
        "clinvar:18",  # then absolute risk, larger first
        "clinvar:16",
        "clinvar:14",  # then relative effect, larger first
        "clinvar:15",
        "clinvar:12",  # then gene name
        "clinvar:11",
    ]


def test_json_has_exact_keys(con: duckdb.DuckDBPyConnection) -> None:
    effect = {
        "relative": {"measure": "OR", "value": 2.0, "ci": [1.5, 2.5]},
        "absolute": None,
    }
    seed(
        con, [Row("clinvar:21", rsid="rs123", effect=effect, flags=["lowmappability"])]
    )
    result = query(con, None, Filters())
    document = json.loads(json.dumps(result.to_json()))
    assert set(document) == {"run_id", "versions", "findings", "hidden", "unassessed"}
    assert document["unassessed"] == {"total": 0, "by_gene": {}}
    assert document["run_id"] == RUN
    assert document["versions"] == TOOLS
    [finding] = document["findings"]
    assert set(finding) == {
        "id",
        "kind",
        "topic",
        "gene",
        "marker",
        "genotype",
        "reliability",
        "validity",
        "effect",
        "actionable",
        "claim",
        "sources",
        "gene_call",
        "dates",
    }
    assert finding["gene_call"] is None  # a variant finding has none
    assert set(finding["marker"]) == {"rsid", "grch38"}
    assert set(finding["genotype"]) == {"gt", "dosage", "allele"}
    assert set(finding["reliability"]) == {
        "level",
        "p_wrong",
        "p_wrong_ci",
        "gq",
        "flags",
        "warnings",
        "row_type",
    }
    assert set(finding["validity"]) == {"level", "source_level", "reason"}
    assert [set(source) for source in finding["sources"]] == [
        {"source", "id", "version"}
    ]
    assert finding["marker"] == {"rsid": "rs123", "grch38": "chr17:1000:C:T"}
    assert finding["reliability"]["flags"] == ["lowmappability"]
    assert finding["reliability"]["p_wrong"] == pytest.approx(1e-4)
    assert finding["effect"] == effect
    assert finding["claim"]["variation_id"] == 21
    assert finding["sources"] == [
        {"source": "clinvar", "id": "21", "version": "20261011"}
    ]


def test_filters_count_under_the_first_failing_filter(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            Row("clinvar:31", gene="BRCA1", topic="risks", pos=1031, rsid="rs31"),
            Row("clinvar:32", gene="BRCA1", topic="carrier", pos=1032),
            Row("clinvar:33", gene="CFTR", topic="risks", pos=1033),
            Row("clinvar:34", gene="BRCA1", topic="risks", pos=1034, level="low"),
            Row("clinvar:35", gene="BRCA1", topic="risks", pos=1035, rsid="rs35"),
        ],
    )
    result = query(
        con,
        None,
        Filters(topic=Topic.RISKS, gene="brca1", marker="RS31"),
    )
    assert [f["id"] for f in result.findings] == ["clinvar:31"]
    assert result.hidden == {
        "reliability_below:medium": 1,
        "topic": 1,
        "gene": 1,
        "marker": 1,
    }


@pytest.mark.parametrize(
    ("marker", "found"),
    [
        ("rs31", True),
        ("RS31", True),
        ("GRCh38:chr17:1031:C:T", True),
        ("grch38:chr17:1031:c:t", True),
        ("chr17:1031:C:T", True),
        ("17:1031:C:T", True),
        ("GRCh38:17:1031:C:T", True),
        ("chr17:1031:C:G", False),
        ("chr17:1031", False),
        ("rs3", False),
    ],
)
def test_query_marker_by_grch38_vid_and_rsid(
    con: duckdb.DuckDBPyConnection, marker: str, found: bool
) -> None:
    """`--marker` finds a site by its rsID or its GRCh38 vid, with or without the
    `GRCh38:` prefix and the `chr`; a position alone names no allele."""
    seed(con, [Row("clinvar:31", pos=1031, rsid="rs31"), Row("clinvar:32", pos=1032)])
    assert ids(con, Filters(marker=marker)) == (["clinvar:31"] if found else [])
    sites = genotype(con, None, marker=marker).sites
    assert [s["id"] for s in sites] == (["clinvar:31"] if found else [])


@pytest.mark.parametrize(("chrom", "noted"), [("chrX", True), ("chr17", False)])
def test_grch38_row_on_a_sex_chromosome_is_noted(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], chrom: str, noted: bool
) -> None:
    """The dosage note follows the vid's GRCh38 chromosome."""
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(con, [Row("clinvar:94", chrom=chrom, pos=1094)])
    con.close()
    assert main(["show", "clinvar:94", "--db", str(path)]) == 0
    assert (HAPLOID_DOSAGE_NOTE in capsys.readouterr().out) is noted


@pytest.mark.parametrize("text", ["ovarian", "LIKELY_PATH", "clinvar:41", "brca"])
def test_text_searches_gene_condition_classification_and_id(
    con: duckdb.DuckDBPyConnection, text: str
) -> None:
    seed(
        con,
        [
            Row(
                "clinvar:41",
                gene="BRCA2",
                pos=1041,
                condition="Breast-ovarian cancer",
                classification="Likely_pathogenic",
            ),
            Row(
                "clinvar:42",
                gene="LDLR",
                pos=1042,
                condition="Familial hypercholesterolemia",
            ),
        ],
    )
    result = query(con, None, Filters(text=text))
    assert [f["id"] for f in result.findings] == ["clinvar:41"]
    assert result.hidden == {"text": 1}


@pytest.mark.parametrize("text", ["hereditary hemochromatosis", "Hereditary_hemo"])
def test_text_reads_an_underscore_as_a_space(
    con: duckdb.DuckDBPyConnection, text: str
) -> None:
    """ClinVar stores a condition with `_` for a space; the user types spaces."""
    condition = "Hereditary_hemochromatosis|not_provided"
    seed(con, [Row("clinvar:45", gene="HFE", pos=1045, condition=condition)])
    result = query(con, None, Filters(text=text))
    assert [f["id"] for f in result.findings] == ["clinvar:45"]


def test_text_searches_a_drug_claims_drug_names(con: duckdb.DuckDBPyConnection) -> None:
    """A drug question names a drug: `--text` finds the pharmacogenomic claims about it."""
    claim = {
        "schema": "genomelab.claim/1",
        "type": "pharmgkb_clinical_annotation",
        "annotation_id": "655384799",
        "gene": "TNF",
        "category": "Efficacy",
        "condition": "Psoriasis",
        "drugs": [{"drug": "etanercept"}],
    }
    seed(
        con,
        [
            Row("pharmgkb:655384799", gene="TNF", topic="drugs", pos=1043, claim=claim),
            Row("clinvar:44", gene="LDLR", pos=1044),
        ],
    )
    result = query(con, None, Filters(text="Etanercept"))
    assert [f["id"] for f in result.findings] == ["pharmgkb:655384799"]


def test_latest_run_by_default_and_named_run_on_request(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(con, [Row("clinvar:51")], run_id=OLD_RUN)
    seed(con, [Row("clinvar:52", pos=1052)], run_id=RUN)
    latest = query(con, None, Filters())
    assert (latest.run_id, ids(con)) == (RUN, ["clinvar:52"])
    old = query(con, OLD_RUN, Filters())
    assert (old.run_id, [f["id"] for f in old.findings]) == (OLD_RUN, ["clinvar:51"])
    with pytest.raises(RunNotFound):
        query(con, "20260101T000000Z-cccccc", Filters())


def test_each_run_reads_its_own_variant_row(con: duckdb.DuckDBPyConnection) -> None:
    # One allele, named by each run's own ClinVar release.
    seed(con, [Row("clinvar:53", rsid="rs-old")], run_id=OLD_RUN)
    seed(con, [Row("clinvar:53", rsid="rs-new")], run_id=RUN)
    for run_id, rsid in ((OLD_RUN, "rs-old"), (RUN, "rs-new")):
        [finding] = query(con, run_id, Filters()).findings
        assert finding["marker"]["rsid"] == rsid
        _run, shown = show(con, run_id, "clinvar:53")
        assert shown is not None and shown["marker"]["rsid"] == rsid


def test_empty_store_has_no_latest_run(con: duckdb.DuckDBPyConnection) -> None:
    with pytest.raises(RunNotFound):
        query(con, None, Filters())


def test_show_returns_one_finding_whatever_the_filters(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(con, [Row("clinvar:61", level="low", validity="limited")])
    run, finding = show(con, None, "clinvar:61")
    assert run == RUN
    assert finding is not None and finding["id"] == "clinvar:61"
    assert show(con, None, "clinvar:62") == (RUN, None)


def test_hidden_counts_come_in_filter_order(con: duckdb.DuckDBPyConnection) -> None:
    # Seeded in the reverse of the filter order, so insertion order cannot pass for it.
    seed(
        con,
        [
            # Shares the kept row's variant (no second call row), so only text fails.
            Row("clinvar:81", pos=1087, call=False, condition="x"),
            Row("clinvar:82", pos=1082, rsid="rs1"),
            Row("clinvar:83", pos=1083, gene="CFTR"),
            Row("clinvar:84", pos=1084, topic="carrier"),
            Row("clinvar:85", pos=1085, validity="limited"),
            Row("clinvar:86", pos=1086, level="low"),
            Row("clinvar:87", pos=1087),
        ],
    )
    f = Filters(
        topic=Topic.RISKS, gene="BRCA1", marker="chr17:1087:C:T", text="ovarian"
    )
    result = query(con, None, f)
    assert [x["id"] for x in result.findings] == ["clinvar:87"]
    assert list(result.hidden) == [
        "reliability_below:medium",
        "validity_below:moderate",
        "topic",
        "gene",
        "marker",
        "text",
    ]


def test_an_unassessed_claim_without_a_call_row_is_unknown(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(con, [_unassessed("clinvar:91", 1091, rsid="rs91", call=False)])
    result = query(con, None, Filters(), list_unassessed=True)
    con.close()
    assert result.findings == []
    [site] = result.unassessed_sites or []
    assert site["marker"]["grch38"] == "chr17:1091:C:T"
    # The allele a call would count is known without the call.
    assert site["genotype"] == {"gt": None, "dosage": None, "allele": "T"}
    assert set(site["reliability"].values()) == {None}
    assert result.unassessed["by_gene"]["BRCA1"]["reasons"] == {"no call row": 1}
    assert main(["show", "clinvar:91", "--db", str(path)]) == 0
    out = capsys.readouterr().out
    assert "genotype: no call row - unknown, not a result" in out
    assert "kind: unassessed" in out
    assert "gene-level" not in out


def test_cli_query_prints_unassessed_counts_and_lists_them_on_request(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(con, [_unassessed("clinvar:95", 1095), Row("clinvar:96", pos=1096)])
    con.close()
    assert main(["query", "--db", str(path)]) == 0
    out = capsys.readouterr().out
    assert "run 20261018T120000Z-aaaaaa: 1 findings" in out
    assert "unassessed (genotype unknown or not covered, not a result): 1 sites" in out
    assert f"  BRCA1: 1 (risks 1) - {NOCALL} 1" in out
    assert "clinvar:95" not in out
    assert main(["query", "--db", str(path), "--list-unassessed"]) == 0
    out = capsys.readouterr().out
    assert "id: clinvar:95  kind: unassessed" in out
    assert main(["query", "--db", str(path), "--json"]) == 0
    document = json.loads(capsys.readouterr().out)
    assert "unassessed_sites" not in document
    assert document["unassessed"]["total"] == 1
    assert main(["query", "--db", str(path), "--json", "--list-unassessed"]) == 0
    document = json.loads(capsys.readouterr().out)
    assert [s["id"] for s in document["unassessed_sites"]] == ["clinvar:95"]
    assert main(["show", "clinvar:95", "--db", str(path), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["finding"]["kind"] == "unassessed"


@pytest.mark.parametrize(
    ("chrom", "noted"), [("chrX", True), ("chrY", True), ("chr17", False)]
)
def test_text_output_explains_dosage_on_a_sex_chromosome(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], chrom: str, noted: bool
) -> None:
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(con, [Row("clinvar:93", chrom=chrom, pos=1093)])
    con.close()
    assert main(["show", "clinvar:93", "--db", str(path)]) == 0
    out = capsys.readouterr().out
    assert (HAPLOID_DOSAGE_NOTE in out) is noted
    assert "1 = hemizygous carrier of the allele" in HAPLOID_DOSAGE_NOTE


@pytest.fixture
def db(tmp_path: Path) -> str:
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(
        con,
        [
            Row("clinvar:71", rsid="rs71", pos=1071),
            Row("clinvar:72", level="low", pos=1072),
        ],
    )
    con.close()
    return str(path)


def test_cli_query_json_roundtrip(db: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["query", "--json", "--db", db]) == 0
    document = json.loads(capsys.readouterr().out)
    assert document["run_id"] == RUN
    assert [f["id"] for f in document["findings"]] == ["clinvar:71"]
    assert document["hidden"] == {"reliability_below:medium": 1}


def test_cli_query_passes_filters(db: str, capsys: pytest.CaptureFixture[str]) -> None:
    argv = [
        "query",
        "--json",
        "--db",
        db,
        "--min-reliability",
        "low",
        "--marker",
        "rs71",
    ]
    assert (
        main([*argv, "--run", RUN, "--topic", "risks", "--min-validity", "strong"]) == 0
    )
    document = json.loads(capsys.readouterr().out)
    assert [f["id"] for f in document["findings"]] == ["clinvar:71"]
    assert document["hidden"] == {"marker": 1}


def test_cli_query_text_shows_findings_and_hidden_counts(
    db: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["query", "--db", db]) == 0
    out = capsys.readouterr().out
    assert "BRCA1" in out and "rs71" in out and "clinvar:71" in out
    assert "p_wrong 1.0e-04" in out
    assert "(measured on GIAB reference samples, DNBSEQ-T7)" in out
    assert "reliability_below:medium: 1" in out


def test_cli_query_prints_hidden_counts_in_filter_order(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "genomelab.duckdb"
    con = connect(path)
    seed(
        con,
        [
            Row("clinvar:84", pos=1084, topic="carrier"),
            Row("clinvar:86", pos=1086, level="low"),
        ],
    )
    con.close()
    assert main(["query", "--db", str(path), "--topic", "risks", "--json"]) == 0
    assert list(json.loads(capsys.readouterr().out)["hidden"]) == [
        "reliability_below:medium",
        "topic",
    ]
    assert main(["query", "--db", str(path), "--topic", "risks"]) == 0
    out = capsys.readouterr().out
    assert out.index("reliability_below:medium") < out.index("topic: 1")


def test_cli_query_nothing_matches_exits_zero(
    db: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["query", "--json", "--db", db, "--gene", "NOPE"]) == 0
    assert json.loads(capsys.readouterr().out)["findings"] == []


def test_cli_show(db: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["show", "clinvar:72", "--json", "--db", db]) == 0
    document = json.loads(capsys.readouterr().out)
    assert set(document) == {"run_id", "finding"}
    assert (document["run_id"], document["finding"]["id"]) == (RUN, "clinvar:72")
    assert main(["show", "clinvar:72", "--db", db]) == 0
    out = capsys.readouterr().out
    assert out.startswith(f"run {RUN}") and "clinvar:72" in out
    assert main(["show", "clinvar:99", "--db", db]) == 1
    assert "clinvar:99" in capsys.readouterr().err


def test_cli_missing_store_exits_2(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = tmp_path / "nowhere.duckdb"
    assert main(["query", "--db", str(missing)]) == 2
    assert str(missing) in capsys.readouterr().err
    assert not missing.exists()


@pytest.mark.parametrize("command", ["query", "show"])
@pytest.mark.parametrize("offset", [-1, 1])
def test_cli_a_store_of_another_schema_version_exits_2_in_one_line(
    db: str, capsys: pytest.CaptureFixture[str], command: str, offset: int
) -> None:
    with duckdb.connect(db) as con:
        con.execute("UPDATE meta SET schema_version = ?", [SCHEMA_VERSION + offset])
    argv = [command, *(["clinvar:71"] if command == "show" else []), "--db", db]
    assert main(argv) == 2
    err = capsys.readouterr().err
    assert err.startswith(f"genokeeper {command}: ")
    assert db in err and f"schema version {SCHEMA_VERSION + offset}" in err
    assert len(err.strip().splitlines()) == 1


def test_cli_missing_run_exits_2(db: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["query", "--db", db, "--run", "20260101T000000Z-cccccc"]) == 2
    assert "20260101T000000Z-cccccc" in capsys.readouterr().err
    assert main(["show", "clinvar:71", "--db", db, "--run", "nope"]) == 2


@pytest.fixture
def seeded_store_with_snapshots(tmp_path: Path) -> Path:
    """A store holding the latest of three snapshots beside it: `clinvar:1` absent on
    2026-10-20, limited on 2026-10-24 and strong on 2026-11-30."""
    root = tmp_path / "snapshots"
    tools, inputs = {"clinvar": "20261011"}, {"clinvar38": sm.file("cv1")}
    first, second, latest = (
        sm.run_id(20),
        sm.run_id(24),
        sm.run_id(30, month=11),
    )
    sm.write(root, first, sm.started(first), tools, inputs, [], [sm.call()])
    for rid, validity in ((second, "limited"), (latest, "strong")):
        sm.write(
            root,
            rid,
            sm.started(rid),
            tools,
            inputs,
            [sm.statement("clinvar:1", validity=validity)],
            [sm.call()],
        )
    db = tmp_path / "genomelab.duckdb"
    con = sm.store(db)
    sm.fill(
        con,
        latest,
        sm.started(latest),
        tools,
        inputs,
        [sm.statement("clinvar:1")],
        [sm.call()],
    )
    con.close()
    return db


def test_query_and_show_json_carry_dates(
    seeded_store_with_snapshots: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    db = str(seeded_store_with_snapshots)
    expected = {
        "first_found": "2026-10-24",
        "last_changed": "2026-11-30",
        "assessment_changed": "2026-11-30",
        "on_or_before": False,
    }
    assert main(["query", "--db", db, "--json"]) == 0
    (finding,) = json.loads(capsys.readouterr().out)["findings"]
    assert finding["dates"] == expected
    assert main(["show", "clinvar:1", "--db", db, "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["finding"]["dates"] == expected


def test_an_older_runs_dates_come_from_the_snapshots_up_to_it(
    seeded_store_with_snapshots: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The 2026-10-24 run's `limited` was last assessed that day, not on 2026-11-30."""
    db, second = seeded_store_with_snapshots, sm.run_id(24)
    con = duckdb.connect(str(db))
    sm.fill(
        con,
        second,
        sm.started(second),
        {"clinvar": "20261011"},
        {"clinvar38": sm.file("cv1")},
        [sm.statement("clinvar:1", validity="limited")],
        [sm.call()],
    )
    con.close()
    on_the_day = ["2026-10-24", "2026-10-24", "2026-10-24", False]
    for command in (["query", "--min-validity", "limited"], ["show", "clinvar:1"]):
        assert main([*command, "--db", str(db), "--run", second, "--json"]) == 0
        document = json.loads(capsys.readouterr().out)
        (finding,) = document.get("findings") or [document["finding"]]
        assert list(finding["dates"].values()) == on_the_day, command


def test_dates_are_null_without_snapshots(con: duckdb.DuckDBPyConnection) -> None:
    seed(con, [Row("clinvar:1")])
    (finding,) = query(con, None, Filters()).findings
    assert finding["dates"] is None


# --- slice 3c: a gene's coverage --------------------------------------------------------------

COVERAGE_KEYS = {
    *("bases", "callable", "mean_depth", "difficult", "outside_benchmark", "owner"),
    *("reliable", "note"),
}


def test_query_gene_carries_coverage(
    con: duckdb.DuckDBPyConnection, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    seed(con, [Row("clinvar:1", gene="CBS"), Row("clinvar:2", gene="PMS2", pos=1002)])
    seed_coverage(con, "CBS", 0.99)
    seed_coverage(con, "PMS2", 0.80)
    cbs = query(con, None, Filters(gene="cbs")).to_json()["coverage"]
    assert set(cbs) == COVERAGE_KEYS
    assert (cbs["reliable"], cbs["note"], cbs["callable"]) == (True, None, 0.99)
    pms2 = query(con, None, Filters(gene="PMS2")).to_json()["coverage"]
    assert (pms2["reliable"], pms2["note"]) == (
        False,
        "not reliably assessed: 80% callable",
    )
    assert "coverage" not in query(con, None, Filters()).to_json()
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["query", "--gene", "PMS2", "--db", db]) == 0
    assert "coverage of PMS2: not reliably assessed: 80% callable" in (
        capsys.readouterr().out
    )


def test_query_gene_without_coverage_row(con: duckdb.DuckDBPyConnection) -> None:
    seed(con, [Row("clinvar:1")])
    found = query(con, None, Filters(gene="SLC37A4")).to_json()["coverage"]
    assert set(found) == COVERAGE_KEYS
    assert found["note"] == "not measured" and found["bases"] is None


def test_owned_gene_is_not_called_unreliable(
    con: duckdb.DuckDBPyConnection, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    """SMN1's CDS is a paralog the masked reference cannot tell from SMN2's: its callable
    share is near zero, but its answer is the stage-3 source's, measured by its own
    calibration."""
    seed(con, [Row("clinvar:21", gene="SMN1", kind="non_finding", dosage="0")])
    seed_coverage(con, "SMN1", 0.03, owner=SMN_TOOL)
    found = query(con, None, Filters(gene="SMN1")).to_json()["coverage"]
    assert (found["reliable"], found["owner"], found["note"]) == (
        False,
        "smn",
        "assessed by smn",
    )
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["genotype", "--gene", "SMN1", "--db", db]) == 0
    out = capsys.readouterr().out
    assert "SMN1:" in out and "not reliably assessed" not in out


# --- slice 3f: research layers ---------------------------------------------------------------


def _research_row(con: duckdb.DuckDBPyConnection) -> None:
    """A carried level-3 PharmGKB annotation, stored as its research source writes it."""
    claim = {
        "schema": "genomelab.claim/1",
        "type": "pharmgkb_clinical_annotation",
        "annotation_id": "655384602",
        "gene": "RGS4",
        "category": "Efficacy",
        "condition": "Schizophrenia",
        "drugs": [{"drug": "risperidone"}],
    }
    seed(
        con,
        [
            Row(
                "pharmgkb_research:655384602",
                gene="RGS4",
                validity="limited",
                topic="research",
                rsid="rs951439",
                claim=claim,
            )
        ],
    )
    con.execute(
        "UPDATE statement SET source = 'pharmgkb_research',"
        " source_level = 'ClinPGx level 3'"
    )


def test_research_scores_hidden_until_named(con: duckdb.DuckDBPyConnection) -> None:
    """A research score is a research layer's: counted, not returned, until a gene, marker
    or text names it; `--text <condition>` finds it at `--min-validity limited`."""
    seed(
        con,
        [
            Row(
                "pgsc_calc:research:PGS002746",
                validity="limited",
                topic="research",
                claim={"type": "pgs_research", "condition": "ADHD"},
            )
        ],
    )
    con.execute("UPDATE statement SET source = 'pgsc_calc_research'")
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    for f in (Filters(), loosest, replace(loosest, topic=Topic.RESEARCH)):
        result = query(con, None, f)
        assert (result.findings, result.hidden) == ([], {"research_association": 1})
    limited = Filters(min_validity=Validity.LIMITED)
    for text in ("adhd", "pgsc_calc:research:", "PGS002746"):
        assert ids(con, replace(limited, text=text)) == ["pgsc_calc:research:PGS002746"]


def test_panel_variants_hidden_until_named(con: duckdb.DuckDBPyConnection) -> None:
    """A PanelApp panel statement is a research layer's: counted, not returned, until a
    gene, marker or text names it; `--text panelapp:` lists them all at
    `--min-validity limited`."""
    sid = "panelapp:GRCh38:chr17:1000:C:T:ENST01"
    seed(
        con,
        [
            Row(
                sid,
                gene="GENEZ",
                validity="limited",
                topic="research",
                claim={"type": "panelapp_variant", "classification": "predicted"},
            )
        ],
    )
    con.execute("UPDATE statement SET source = 'panelapp'")
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    for f in (Filters(), loosest, replace(loosest, topic=Topic.RESEARCH)):
        result = query(con, None, f)
        assert (result.findings, result.hidden) == ([], {"research_association": 1})
    limited = Filters(min_validity=Validity.LIMITED)
    for named in (
        replace(limited, gene="genez"),
        replace(limited, text="panelapp:"),
        replace(limited, marker="GRCh38:chr17:1000:C:T"),
    ):
        assert ids(con, named) == [sid], named


def test_level_3_hidden_by_default(con: duckdb.DuckDBPyConnection) -> None:
    """A level-3 annotation is a research association: counted, not returned, until a gene,
    marker or text (its drug or its phenotype) names it."""
    _research_row(con)
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    for f in (Filters(), loosest, replace(loosest, topic=Topic.RESEARCH)):
        result = query(con, None, f)
        assert (result.findings, result.hidden) == ([], {"research_association": 1})
    for narrowed in (
        replace(loosest, gene="rgs4"),
        replace(loosest, marker="rs951439"),
        replace(loosest, text="risperidone"),
        replace(loosest, text="schizophrenia"),
    ):
        assert ids(con, narrowed) == ["pharmgkb_research:655384602"], narrowed


def test_one_research_predicate() -> None:
    """Every research test goes through `is_research_layer`: no module compares a source
    with one research tool's name, so a fifth copy of the test fails here."""
    package = Path(__file__).parents[1] / "src" / "genokeeper"
    assert package.is_dir()
    pattern = re.compile(r"[!=]=\s*(GWAS_TOOL|PHARMGKB_RESEARCH_TOOL)\b")
    # Slice 3d's GBA1 rule rewrites a GWAS claim's own `risk_allele_copies`: a test of the
    # GWAS claim's shape, not of a research layer.
    gwas_claim = 'gwas = row.statement.get("source") == GWAS_TOOL'
    hits = [
        f"{path.relative_to(package)}:{n}"
        for path in sorted(package.rglob("*.py"))
        for n, line in enumerate(path.read_text().splitlines(), 1)
        if pattern.search(line) and line.strip() != gwas_claim
    ]
    assert hits == []


def test_a_gene_call_joins_only_the_statement_standing_for_its_gene(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """A gene call is the result of `<tool>:<gene>` alone: an HLA allele's unassessed
    statement and PharmCAT's G6PD gap on genes PharmCAT called (slice 3f) carry their own
    reasons, not the gene call's, and the query reads the store."""
    claim = {"schema": "genomelab.claim/1", "type": "hla_allele"}
    seed(
        con,
        [
            Row(
                "pharmcat:HLA-B",
                gene="HLA-B",
                level=None,
                topic="drugs",
                gene_call=_gene_call("high"),
                claim=PGX_CLAIM,
            ),
            Row(
                "pharmcat:G6PD",
                gene="G6PD",
                level=None,
                topic="drugs",
                gene_call=_gene_call("high"),
                claim=PGX_CLAIM,
            ),
            Row(
                "hla:HLA-B*27",
                kind="unassessed",
                gene="HLA-B",
                level=None,
                claim=claim,
                reason="not calibrated: run `genokeeper calibrate`",
            ),
            Row(
                "pharmcat:G6PD-gap",
                kind="unassessed",
                gene="G6PD",
                level=None,
                topic="drugs",
                claim=PGX_CLAIM,
                reason="not covered: PharmCAT does not call Villeurbanne",
            ),
        ],
    )
    result = query(con, None, Filters(), list_unassessed=True)
    assert result.unassessed_sites is not None
    sites = {s["id"]: s for s in result.unassessed_sites}
    assert {i: s["gene_call"] for i, s in sites.items()} == {
        "hla:HLA-B*27": None,
        "pharmcat:G6PD-gap": None,
    }
    assert result.unassessed["by_gene"]["HLA-B"]["reasons"] == {
        "not calibrated: run `genokeeper calibrate`": 1
    }
    found = {f["id"]: f["gene_call"] for f in result.findings}
    assert found["pharmcat:HLA-B"] is not None and found["pharmcat:G6PD"] is not None


def test_two_tools_on_one_gene_read_their_own_gene_calls(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """PharmCAT's CFTR result and CFTR's poly-T/TG tracts are two statements on one gene:
    each reads its own call, by the name its statement id ends with."""
    seed(
        con,
        [
            Row(
                "pharmcat:CFTR",
                gene="CFTR",
                level=None,
                topic="drugs",
                gene_call=_gene_call("high", diplotype="drug response"),
                claim=PGX_CLAIM,
            ),
            Row(
                "expansionhunter:CFTR-polyT-TG",
                gene="CFTR",
                level=None,
                gene_call=_gene_call(
                    "high", diplotype="tract lengths", gene="CFTR-polyT-TG"
                ),
                claim={"schema": "genomelab.claim/1", "type": "repeat_expansion"},
            ),
        ],
    )
    found = {f["id"]: f for f in query(con, None, Filters()).findings}
    assert found["pharmcat:CFTR"]["gene_call"]["diplotype"] == "drug response"
    assert (
        found["expansionhunter:CFTR-polyT-TG"]["gene_call"]["diplotype"]
        == "tract lengths"
    )


@pytest.mark.parametrize("kind", ["unassessed", "finding"])
def test_a_statement_without_a_gene_reads_no_gene_call(
    con: duckdb.DuckDBPyConnection, kind: str
) -> None:
    """RBCeq2's ABCG2 system (`bloodgroup:ABCG2`, no gene) shares its last word with
    PharmCAT's ABCG2 and is not its result: read as PharmCAT's, an unassessed one raises
    and a finding shows another tool's call."""
    seed(
        con,
        [
            Row(
                "pharmcat:ABCG2",
                gene="ABCG2",
                level=None,
                topic="drugs",
                gene_call=_gene_call("high"),
                claim=PGX_CLAIM,
            ),
            Row(
                "bloodgroup:ABCG2",
                kind=kind,
                gene=None,
                level=None,
                validity="limited",
                claim={"schema": "genomelab.claim/1", "type": "blood_group"},
                reason="not calibrated",
            ),
        ],
    )
    result = query(
        con,
        None,
        Filters(min_reliability=Level.LOW, min_validity=Validity.LIMITED),
        list_unassessed=True,
    )
    assert result.unassessed_sites is not None
    shown = {s["id"]: s for s in [*result.findings, *result.unassessed_sites]}
    assert shown["bloodgroup:ABCG2"]["gene_call"] is None
    assert shown["pharmcat:ABCG2"]["gene_call"] is not None
