"""`genokeeper genotype`: an answer per site from its dosage and reliability, a gene's answer
as counts with what was not looked at, and the genes whose main mechanism is not assessed."""

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures.run_context import run_context
from fixtures.store_mini import RUN, Row, seed, seed_coverage

from genokeeper.carrier import (
    NOT_COVERED,
    Answer,
    _gene_level_site,
    _gene_lines,
    _site_block,
    genotype,
    load_mechanisms,
    possible_carrier,
)
from genokeeper.cli import main
from genokeeper.evidence.predicted import NOT_HAPLOINSUFFICIENT
from genokeeper.insertions import INSERTIONS_TOOL
from genokeeper.paralogs import LOCI, GauchianCall, PsCall
from genokeeper.query import show
from genokeeper.sources import write_collected
from genokeeper.sources.cnv import CNV_TOOL, WHOLE_GAIN
from genokeeper.sources.deletions import CONFIRM_BY
from genokeeper.sources.paralogs import LpaSource, ParalogSource
from genokeeper.sources.repeats import REPEATS_TOOL
from genokeeper.sources.smn import SMN_TOOL
from genokeeper.store.db import SCHEMA_VERSION, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

NOCALL = "genotype unknown: nocall"


@pytest.fixture
def con(tmp_path: Path) -> duckdb.DuckDBPyConnection:
    return connect(tmp_path / "genomelab.duckdb")


def test_not_covered_names_the_repeat_loci_of_3a_and_3e() -> None:
    assert (
        "repeat expansions outside the 24 loci, RFC1 and CFTR's poly-T/TG of stages 3a"
        " and 3e (FMR1 and AFF2 follow the X dosage; the loci short reads cannot assess"
        " are stated as not covered, with the reason)"
    ) in NOT_COVERED
    assert "repeat expansions" not in NOT_COVERED


@pytest.mark.parametrize("kind", ["finding", "non_finding", "unassessed"])
def test_cftr_5t_is_not_cf_carrier(
    con: duckdb.DuckDBPyConnection,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    kind: str,
) -> None:
    """CFTR's poly-T/TG statement (`expansionhunter:CFTR`) is never cystic-fibrosis carrier
    status: with no CF-causing call in the store, the CF sites and the gene's counts are
    what they are without it, and the statement answers neither carrier nor not a carrier;
    the text prints it as its own line."""
    cf = Row("clinvar:41", gene="CFTR", kind="non_finding", dosage="0", level="high")
    repeat = Row(
        "expansionhunter:CFTR",
        kind=kind,
        gene="CFTR",
        level=None,
        validity="limited",
        topic="carrier",
        claim={
            "schema": "genomelab.claim/1",
            "type": "repeat_expansion",
            "gene": "CFTR",
            "locus_id": "CFTR",
            "alleles": {
                "tg": [{"repeats": 12}, {"repeats": 11}],
                "t": [{"repeats": 5}, {"repeats": 7}],
            },
        },
    )
    seed(con, [cf, repeat])
    con.execute(
        "UPDATE statement SET source = 'expansionhunter' WHERE source_id = 'CFTR'"
    )
    result = genotype(con, None, gene="CFTR").to_json()
    [lead, site] = result["sites"]
    assert lead["id"] == "expansionhunter:CFTR" and site["id"] == "clinvar:41"
    assert site["answer"] == Answer.NOT_A_CARRIER
    assert lead["answer"] == Answer.UNKNOWN
    assert "not cystic-fibrosis carrier status" in str(lead["not_a_finding_because"])
    summary = result["gene"]
    assert summary["carried"] == 0
    assert summary["by_answer"] == {
        "carrier": 0,
        "not_a_carrier": 1,
        "probably_not_a_carrier": 0,
        "unknown": 0,
        "needs_confirmation": 0,
        "possible_carrier": 0,
        "information": 0,
    }
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["genotype", "--gene", "CFTR", "--db", db]) == 0
    out = capsys.readouterr().out
    assert "CFTR poly-T/TG: not cystic-fibrosis carrier status" in out
    assert "answer: carrier of the allele" not in out
    assert "answer: not a carrier of the allele" in out  # the CF site's own answer


@pytest.mark.parametrize(
    ("kind", "reason", "confirm", "expected"),
    [
        ("finding", "possible carrier; needs confirmation", None, "unknown"),
        ("non_finding", "an expansion is not seen", "PCR", "unknown"),
        ("non_finding", None, None, "not_a_carrier"),
        ("finding", "screening result, needs confirmation", None, "carrier"),
    ],
)
def test_screening_repeat_statements_answer_unknown(
    con: duckdb.DuckDBPyConnection,
    kind: str,
    reason: str | None,
    confirm: str | None,
    expected: str,
) -> None:
    """A slice-3e "possible" finding never answers `carrier`, nor an exclusion-flagged normal
    (its claim names a confirming test) `not_a_carrier`; 3a's statements answer as ever."""
    repeat = Row(
        "expansionhunter:RFC1",
        kind=kind,
        gene="RFC1",
        level=None,
        validity="limited",
        reason=reason,
        claim={
            "schema": "genomelab.claim/1",
            "type": "repeat_expansion",
            "gene": "RFC1",
            "confirm": confirm,
        },
    )
    seed(con, [repeat])
    con.execute(
        "UPDATE statement SET source = 'expansionhunter' WHERE source_id = 'RFC1'"
    )
    [lead] = genotype(con, None, gene="RFC1").to_json()["sites"]
    assert lead["answer"] == expected


def test_not_covered_names_depth_limits() -> None:
    """Slice 3d screens copy number and insertions: what stays out is named."""
    assert NOT_COVERED[0] == (
        "copy-number variants under about 1 kb (most single-exon events), target exons in"
        " segmental duplications without a Parascopy model, and where an extra whole-gene"
        " copy lies (deletions of 50 bp or more, and depth losses and gains of 1 kb or"
        " more, in the screen genes are screened and can miss one; SLC37A4 is not"
        " covered)"
    )
    assert NOT_COVERED[1] == (
        "structural variants: inversions and translocations (short-read inversion recall"
        " 0.12-0.17 at 30x)"
    )


def _deletion(
    con: duckdb.DuckDBPyConnection,
    gene: str,
    kind: str,
    topic: str,
    reason: str | None = None,
) -> dict[str, Any]:
    """`genotype --gene` of a store holding one het deletion statement in `gene` (source
    `deletions`, at its symbolic vid, its validity reason `reason` when given)."""
    claim = {"schema": "genomelab.claim/1", "type": "deletion", "gene": gene}
    seed(
        con,
        [
            Row(f"deletions:{gene}:chr7:5000-7000", kind=kind, gene=gene, topic=topic,
                level="low", validity="limited", pos=5000, chrom="chr7", claim=claim),
        ],
    )  # fmt: skip
    vid = "GRCh38:chr7:5000:A:<DEL>"
    for table in ("variant", "call", "statement"):
        con.execute(f"UPDATE {table} SET vid = ?", [vid])
    con.execute("UPDATE statement SET source = 'deletions'")
    if reason is not None:
        con.execute("UPDATE statement SET reason = ?", [reason])
    return genotype(con, None, gene=gene).to_json()


def test_possible_carrier_3b_deletion(con: duckdb.DuckDBPyConnection) -> None:
    """A deletion under carrier is a possible carrier, never counted as carried; the
    same deletion under risks needs confirmation."""
    result = _deletion(con, "CFTR", "finding", "carrier", "likely present, …")
    (deletion,) = result["sites"]
    assert deletion["marker"]["grch38"] == "chr7:5000:A:<DEL>"
    assert deletion["genotype"]["dosage"] == "1"
    assert deletion["answer"] == Answer.POSSIBLE_CARRIER
    assert deletion["not_a_finding_because"] is None
    gene = result["gene"]
    assert (gene["deletions"], gene["carried"], gene["possible_carriers"]) == (1, 0, 1)
    assert gene["pathogenic_sites"] == 0
    con.execute("UPDATE statement SET topic = 'risks', gene = 'BRCA1'")
    (risk,) = genotype(con, None, gene="BRCA1").to_json()["sites"]
    assert risk["answer"] == Answer.NEEDS_CONFIRMATION


def test_possible_carrier_text(con: duckdb.DuckDBPyConnection) -> None:
    """The text output says a possible carrier needs confirmation, and its deletion line
    names each way a deletion is confirmed."""
    result = _deletion(con, "CFTR", "finding", "carrier")
    confirm = "confirmation by MLPA, a microarray, qPCR or long reads"
    (site,) = result["sites"]
    assert (
        "possibly a carrier: a screening result that needs confirmation; see its reason"
    ) in _site_block(site)
    lines = _gene_lines("CFTR", result["gene"])
    assert f"  deletions of 50 bp or more: 1, each needing {confirm}" in lines


def _store(
    tmp_path: Path, *sources: ParalogSource | LpaSource
) -> duckdb.DuckDBPyConnection:
    """A store holding one run of `sources`' statements, written as the run writes them."""
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {})
    ctx = replace(run_context(tmp_path), run_id=run_id)
    with transaction(con):
        collected = [source.collect(ctx) for source in sources]
        write_collected(StoreWriter(con, tmp_path / "work"), con, run_id, collected)
    return con


def test_possible_carrier_3d_strc(tmp_path: Path) -> None:
    """ParalogSource's STRC carrier finding (validity limited) answers possible carrier;
    an information copy number (RHD at 0 copies) and Lp(a) answer information, never
    carrier."""
    calls = {locus: PsCall(cn=2, qual=40.0, passed=True, rows=1) for locus in LOCI}
    calls["STRC"] = PsCall(cn=1, qual=40.0, passed=True, rows=1)
    calls["RHD"] = PsCall(cn=0, qual=40.0, passed=True, rows=1)
    calls["GBA1"] = PsCall(cn=4, qual=40.0, passed=True, rows=1)
    source = ParalogSource(
        calls=calls,
        unique_ratio=0.5,
        gauchian=GauchianCall(4, False, False, (), (), None),
        aldy_cyp2a6=2,
        withheld={},
        not_run=None,
    )
    lpa = LpaSource(kiv2=24.3, withheld=None, not_run=None)
    con = _store(tmp_path, source, lpa)
    result = genotype(con, None, gene="STRC").to_json()
    (strc,) = result["sites"]
    assert (strc["id"], strc["kind"], strc["validity"]["level"]) == (
        "paralogs:STRC",
        "finding",
        "limited",
    )
    assert strc["answer"] == Answer.POSSIBLE_CARRIER
    assert (result["gene"]["carried"], result["gene"]["possible_carriers"]) == (0, 1)
    for gene, statement_id in (("RHD", "paralogs:RHD"), ("LPA", "lpa:LPA")):
        (site,) = genotype(con, None, gene=gene).to_json()["sites"]
        assert (site["id"], site["kind"]) == (statement_id, "finding")
        assert site["answer"] == Answer.INFORMATION


def test_possible_carrier_3e_rfc1_shape(con: duckdb.DuckDBPyConnection) -> None:
    """3e's RFC1 carrier call (a repeat locus's gene-level finding at validity limited)
    answers possible carrier; from a source outside POSSIBLE_SOURCES it is not one, and its
    "possible" reason still keeps it from `carrier` (slice 3e's screening rule)."""
    seed(con, [Row("expansionhunter:RFC1", gene="RFC1", level=None, topic="carrier",
                   validity="limited", reason="possible carrier of …",
                   claim={"schema": "genomelab.claim/1", "type": "repeat"})])  # fmt: skip
    con.execute(f"UPDATE statement SET source = '{REPEATS_TOOL}'")
    (rfc1,) = genotype(con, None, gene="RFC1").to_json()["sites"]
    assert rfc1["answer"] == Answer.POSSIBLE_CARRIER
    con.execute("UPDATE statement SET source = 'repeats'")
    _run, finding = show(con, None, "expansionhunter:RFC1")
    assert finding is not None
    assert not possible_carrier(finding)
    assert _gene_level_site(finding)["answer"] == Answer.UNKNOWN


def test_repeat_carrier_at_strong_still_carrier(con: duckdb.DuckDBPyConnection) -> None:
    """3a's FXN carrier finding is at validity strong: it still answers carrier."""
    seed(con, [Row("expansionhunter:FXN", gene="FXN", level=None, topic="carrier",
                   claim={"schema": "genomelab.claim/1", "type": "repeat"})])  # fmt: skip
    con.execute(f"UPDATE statement SET source = '{REPEATS_TOOL}'")
    (fxn,) = genotype(con, None, gene="FXN").to_json()["sites"]
    assert (fxn["validity"]["level"], fxn["answer"]) == ("strong", Answer.CARRIER)


def test_smn1_finding_still_carrier(con: duckdb.DuckDBPyConnection) -> None:
    """3a's SMN1 finding is at validity strong: it still answers carrier."""
    seed(con, [Row("smn:SMN1", gene="SMN1", level=None, topic="carrier",
                   claim={"schema": "genomelab.claim/1", "type": "smn"})])  # fmt: skip
    con.execute(f"UPDATE statement SET source = '{SMN_TOOL}'")
    (smn1,) = genotype(con, None, gene="SMN1").to_json()["sites"]
    assert (smn1["validity"]["level"], smn1["answer"]) == ("strong", Answer.CARRIER)


def test_clinvar_limited_carrier_unchanged(con: duckdb.DuckDBPyConnection) -> None:
    """A ClinVar finding under carrier at validity limited answers as before: carrier."""
    seed(con, [Row("clinvar:41", gene="CFTR", topic="carrier", validity="limited")])
    (site,) = genotype(con, None, gene="CFTR").to_json()["sites"]
    assert site["answer"] == Answer.CARRIER
    assert site["genotype"]["dosage"] == "1"


def _sv_sites(
    con: duckdb.DuckDBPyConnection,
    gene: str,
    events: list[tuple[str, str, str, str]],
    start: int = 5000,
) -> dict[str, Any]:
    """`genotype --gene` of a store holding one het statement per event in `gene`, each
    (source, ALT, kind, reason) at its own position from `start`, under risks at validity
    limited."""
    rows = [
        Row(f"{source}:{gene}:chr7:{start + i}-7000", kind=kind, gene=gene,
            level="low", validity="limited", pos=start + i, chrom="chr7",
            reason=reason, claim={"schema": "genomelab.claim/1", "type": source})
        for i, (source, _alt, kind, reason) in enumerate(events)
    ]  # fmt: skip
    seed(con, rows)
    for i, (source, alt, _kind, _reason) in enumerate(events):
        old, new = f"GRCh38:chr7:{start + i}:C:T", f"GRCh38:chr7:{start + i}:A:{alt}"
        for table in ("variant", "call", "statement"):
            con.execute(f"UPDATE {table} SET vid = ? WHERE vid = ?", [new, old])
        con.execute("UPDATE statement SET source = ? WHERE vid = ?", [source, new])
    return genotype(con, None, gene=gene).to_json()


def test_cnv_and_insertion_answers(con: duckdb.DuckDBPyConnection) -> None:
    """A depth loss and an insertion under risks need confirmation; a whole-gene gain's
    non-finding gives WHOLE_GAIN; the gene counts its depth events and insertions."""
    gain = f"possible copy-number gain, {CONFIRM_BY}; {WHOLE_GAIN}; predicted, not classified"
    result = _sv_sites(
        con,
        "LDLR",
        [
            (CNV_TOOL, "<DEL>", "finding", "possible copy-number loss, …"),
            (CNV_TOOL, "<DUP>", "non_finding", gain),
            (INSERTIONS_TOOL, "<INS>", "finding", "possible insertion, …"),
        ],
    )
    loss, whole, insertion = result["sites"]
    assert (loss["answer"], insertion["answer"]) == (
        Answer.NEEDS_CONFIRMATION,
        Answer.NEEDS_CONFIRMATION,
    )
    assert (whole["answer"], whole["not_a_finding_because"]) == (
        Answer.NEEDS_CONFIRMATION,
        WHOLE_GAIN,
    )
    gene = result["gene"]
    assert (gene["copy_number"], gene["insertions"], gene["deletions"]) == (2, 1, 0)
    assert (gene["carried"], gene["possible_carriers"]) == (0, 0)
    assert "a structural variant that needs confirmation; see its reason" in (
        _site_block(loss)
    )
    gated = f"possible; {NOT_HAPLOINSUFFICIENT}; predicted, not classified"
    other = _sv_sites(
        con,
        "MUTYH",
        [
            (CNV_TOOL, "<DEL>", "non_finding", gated),
            (CNV_TOOL, "<DUP>", "non_finding", "possible copy-number gain, …"),
            (INSERTIONS_TOOL, "<INS>", "non_finding", "possible insertion, …"),
        ],
        start=6000,
    )
    assert [s["not_a_finding_because"] for s in other["sites"]] == [
        f"a copy-number loss is a {NOT_HAPLOINSUFFICIENT}",
        "one copy-number gain in a gene whose condition needs two copies",
        "one insertion in a gene whose condition needs two copies",
    ]


def test_deletion_non_finding_reason(con: duckdb.DuckDBPyConnection) -> None:
    result = _deletion(con, "MUTYH", "non_finding", "risks")
    (deletion,) = (s for s in result["sites"] if s["id"].startswith("deletions:"))
    assert deletion["answer"] == Answer.NEEDS_CONFIRMATION
    assert deletion["not_a_finding_because"] == (
        "one deletion in a gene whose condition needs two copies"
    )


def test_deletion_not_haploinsufficient_reason(con: duckdb.DuckDBPyConnection) -> None:
    reason = f"possible; {NOT_HAPLOINSUFFICIENT}; predicted, not classified"
    result = _deletion(con, "PCSK9", "non_finding", "risks", reason)
    (deletion,) = (s for s in result["sites"] if s["id"].startswith("deletions:"))
    assert deletion["not_a_finding_because"] == (
        "a deletion is a loss of function in a gene neither haploinsufficient nor Tier 3"
    )


def _site(con: duckdb.DuckDBPyConnection, row: Row) -> dict[str, Any]:
    seed(con, [row])
    result = genotype(con, None, marker=f"GRCh38:{row.chrom}:{row.pos}:C:T")
    assert result.assessed
    [site] = result.to_json()["sites"]
    return site


def test_dosage_zero_medium_is_not_a_carrier(con: duckdb.DuckDBPyConnection) -> None:
    site = _site(con, Row("clinvar:1", kind="non_finding", dosage="0", level="medium"))
    assert site["answer"] == Answer.NOT_A_CARRIER
    assert site["genotype"] == {"gt": "0/0", "dosage": "0", "allele": "T"}
    assert site["kind"] == "non_finding"
    assert site["classification"] == "Pathogenic"
    assert site["condition"] == "Hereditary breast ovarian cancer"
    assert "not carried" in str(site["not_a_finding_because"])


def test_dosage_zero_high_is_not_a_carrier_and_one_copy_is_a_carrier(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(
        con,
        [
            Row("clinvar:2", kind="non_finding", dosage="0", level="high", pos=1002),
            Row("clinvar:3", dosage="1", level="low", pos=1003),
        ],
    )
    sites = genotype(con, None, gene="brca1").to_json()["sites"]
    answers = {s["id"]: (s["answer"], s["not_a_finding_because"]) for s in sites}
    assert answers == {
        "clinvar:2": (Answer.NOT_A_CARRIER, "the allele is not carried (dosage 0)"),
        "clinvar:3": (Answer.CARRIER, None),
    }


def test_predicted_one_copy_says_why(con: duckdb.DuckDBPyConnection) -> None:
    """A stored one-copy predicted het (slice 1e) is not a finding because it is one
    predicted copy in a recessive gene, not for a classification outside disease risk."""
    claim = {"type": "predicted_variant", "classification": "predicted, not classified"}
    row = Row("predicted:1", kind="non_finding", gene="MUTYH", claim=claim)
    seed(con, [row])
    con.execute("UPDATE statement SET source = 'vep', topic = NULL")
    [site] = genotype(con, None, gene="MUTYH").to_json()["sites"]
    why = "one copy in an autosomal recessive gene outside the carrier list"
    assert (site["answer"], site["not_a_finding_because"]) == (
        Answer.CARRIER,
        f"predicted, not classified: {why}",
    )


def test_dosage_zero_low_is_probably_not(con: duckdb.DuckDBPyConnection) -> None:
    site = _site(con, Row("clinvar:4", kind="non_finding", dosage="0", level="low"))
    assert site["answer"] == Answer.PROBABLY_NOT_A_CARRIER


def test_unknown_dosage_is_unknown(con: duckdb.DuckDBPyConnection) -> None:
    row = Row(
        "clinvar:5", kind="unassessed", dosage="unknown", level="low", warnings=[NOCALL]
    )
    site = _site(con, row)
    assert site["answer"] == Answer.UNKNOWN
    assert site["not_a_finding_because"] == NOCALL


def test_no_call_row_is_unknown(con: duckdb.DuckDBPyConnection) -> None:
    site = _site(con, Row("clinvar:6", kind="unassessed", call=False))
    assert site["answer"] == Answer.UNKNOWN
    assert site["genotype"]["dosage"] is None
    assert site["not_a_finding_because"] == "no call row"


def test_marker_without_site_is_not_assessed(
    con: duckdb.DuckDBPyConnection, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    seed(con, [Row("clinvar:7", rsid="rs7")])
    result = genotype(con, None, marker="rs999")
    assert result.to_json() == {
        "run_id": RUN,
        "assessed": False,
        "sites": [],
        "gene": None,
    }
    assert genotype(con, None, marker="RS7").assessed
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["genotype", "--marker", "rs999", "--db", db]) == 0
    assert "rs999: not assessed" in capsys.readouterr().out


def test_gene_summary_has_no_gene_level_answer(
    con: duckdb.DuckDBPyConnection, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    seed(
        con,
        [
            Row("clinvar:11", kind="non_finding", dosage="0", level="high", pos=1011),
            Row("clinvar:12", kind="non_finding", dosage="0", level="high", pos=1012),
            Row("clinvar:13", kind="non_finding", dosage="0", level="low", pos=1013),
            Row("clinvar:14", kind="unassessed", call=False, pos=1014),
            # A benign site is no known pathogenic variant: not counted.
            Row(
                "clinvar:15",
                kind="non_finding",
                dosage="0",
                pos=1015,
                classification="Benign",
            ),
            # Another gene's site stays out of BRCA1's answer.
            Row("clinvar:16", gene="BRCA2", dosage="0", kind="non_finding", pos=1016),
        ],
    )
    result = genotype(con, None, gene="BRCA1").to_json()
    assert [s["id"] for s in result["sites"]] == [f"clinvar:{n}" for n in range(11, 16)]
    summary = result["gene"]
    assert "answer" not in summary
    assert summary == {
        "pathogenic_sites": 4,
        "carried": 0,
        "by_answer": {
            "carrier": 0,
            "not_a_carrier": 2,
            "probably_not_a_carrier": 1,
            "unknown": 1,
            "needs_confirmation": 0,
            "possible_carrier": 0,
            "information": 0,
        },
        "by_level": {"high": 2, "medium": 0, "low": 1},
        "possible_carriers": 0,
        "deletions": 0,
        "copy_number": 0,
        "insertions": 0,
        "not_covered": list(NOT_COVERED),
        "dominant_mechanism_not_assessed": None,
        "coverage": summary["coverage"],
    }
    assert summary["coverage"]["note"] == "not measured"
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["genotype", "--gene", "BRCA1", "--db", db]) == 0
    out = capsys.readouterr().out
    gene_lines = [line for line in out.splitlines() if line.startswith("BRCA1:")]
    assert gene_lines
    assert not [line for line in gene_lines if "not a carrier" in line.lower()]
    assert "4 known pathogenic sites" in gene_lines[0]


def test_smn1_flagged(con: duckdb.DuckDBPyConnection) -> None:
    mechanisms = load_mechanisms()
    assert {"SMN1", "HBA1", "HBA2", "CYP21A2", "GBA1", "F8"} <= set(mechanisms)
    seed(con, [Row("clinvar:21", gene="SMN1", kind="non_finding", dosage="0")])
    summary = genotype(con, None, gene="smn1").to_json()["gene"]
    assert summary["dominant_mechanism_not_assessed"] == mechanisms["SMN1"]
    assert "exon 7" in mechanisms["SMN1"]
    # A gene the store holds no site for is still flagged, never answered.
    empty = genotype(con, None, gene="F8").to_json()
    assert empty["assessed"] is False
    assert empty["gene"]["pathogenic_sites"] == 0
    assert "inversion" in empty["gene"]["dominant_mechanism_not_assessed"]


def test_every_mechanism_row_names_its_source() -> None:
    text = (
        Path(__file__).resolve().parents[1]
        / "src/genokeeper/resources/mechanism_not_assessed.tsv"
    ).read_text()
    rows = [line.split("\t") for line in text.splitlines() if not line.startswith("#")][
        1:
    ]
    assert rows and all(len(r) == 3 and r[2].startswith("GeneReviews") for r in rows)


def test_cli_exit_codes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = tmp_path / "genomelab.duckdb"
    with connect(path) as con:
        seed(con, [Row("clinvar:31", rsid="rs31", kind="non_finding", dosage="0")])
    db = str(path)
    assert main(["genotype", "--marker", "rs31", "--json", "--db", db]) == 0
    document = json.loads(capsys.readouterr().out)
    assert document["sites"][0]["answer"] == "not_a_carrier"
    assert main(["genotype", "--db", db]) == 2
    assert "--gene or --marker" in capsys.readouterr().err
    missing = tmp_path / "nowhere.duckdb"
    assert main(["genotype", "--gene", "BRCA1", "--db", str(missing)]) == 2
    assert str(missing) in capsys.readouterr().err
    assert not missing.exists()
    assert main(["genotype", "--gene", "BRCA1", "--db", db, "--run", "nope"]) == 2
    assert "nope" in capsys.readouterr().err
    with duckdb.connect(db) as con:
        con.execute("UPDATE meta SET schema_version = ?", [SCHEMA_VERSION + 1])
    assert main(["genotype", "--gene", "BRCA1", "--db", db]) == 2
    err = capsys.readouterr().err
    assert (
        err.startswith("genokeeper genotype: ") and len(err.strip().splitlines()) == 1
    )


def test_gene_summary_says_not_reliably_assessed(
    con: duckdb.DuckDBPyConnection, capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    """A gene below `[coverage] min_callable` with no stage-3 owner: its head starts with
    the coverage note, so "none carried" never reads as reassurance."""
    seed(con, [Row("clinvar:11", gene="PMS2", kind="non_finding", dosage="0")])
    seed_coverage(con, "PMS2", 0.82)
    summary = genotype(con, None, gene="pms2").to_json()["gene"]
    assert summary["coverage"]["note"] == "not reliably assessed: 82% callable"
    assert _gene_lines("PMS2", summary)[0].startswith(
        "not reliably assessed: 82% callable; PMS2: 1 known pathogenic sites"
    )
    con.close()
    db = str(tmp_path / "genomelab.duckdb")
    assert main(["genotype", "--gene", "PMS2", "--db", db]) == 0
    assert "not reliably assessed: 82% callable; PMS2:" in capsys.readouterr().out
