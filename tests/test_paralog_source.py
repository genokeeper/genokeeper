"""The paralog loci and Lp(a) as sources (`sources.paralogs`, slice 3d): STRC as the one
paralog carrier locus, RHD, AMY1, C4 and CYP2A6 as information, GBA1 from Gauchian in its
four handled shapes, the named genes not covered, and the LPA KIV-2 count from KILDA."""

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from fixtures.run_context import run_context

from genokeeper.deletions.exons import GBA1_HOMOLOGY, load_spans, locus_span
from genokeeper.paralogs import GauchianCall, PsCall
from genokeeper.sources import Collected, Fields
from genokeeper.sources.paralogs import (
    BIALLELIC_WITHHELD,
    GBA_REASON,
    LPA_REASON,
    LPA_TOOL,
    NOT_COVERED,
    PARALOGS_TOOL,
    RHD_NOTE,
    STRC_REASON,
    LpaSource,
    ParalogSource,
    gba1_contradicted,
)
from genokeeper.store.db import Kind, Topic

WILD = GauchianCall(4, False, False, (), (), None)
L483R = "c.1448T>G(p.Leu483Arg)"
CARRIER = GauchianCall(4, True, False, (L483R,), (), 5)
GOOD = PsCall(cn=2, qual=40.0, passed=True, rows=1)
EXONS_9_11 = locus_span(load_spans(), (GBA1_HOMOLOGY[0],))


def _ps(
    cn: int | None, *, qual: float | None = 40.0, rows: int = 1, passed: bool = True
) -> PsCall:
    """Parascopy's call; `cn` None with `passed`: rows that disagree or one unknown."""
    return PsCall(cn=cn, qual=qual, passed=passed, rows=rows)


def _collect(
    tmp_path: Path,
    calls: Mapping[str, PsCall] | None = None,
    *,
    unique_ratio: float | None = 1.0,
    gauchian: GauchianCall | None = WILD,
    aldy_cyp2a6: int | None = 2,
    withheld: Mapping[str, str] | None = None,
    not_run: str | None = None,
) -> Collected:
    every = {locus: GOOD for locus in ("STRC", "RHD", "AMY1", "C4", "CYP2A6")}
    every["GBA1"] = _ps(4)
    source = ParalogSource(
        calls={**every, **(calls or {})},
        unique_ratio=unique_ratio,
        gauchian=gauchian,
        aldy_cyp2a6=aldy_cyp2a6,
        withheld=withheld or {},
        not_run=not_run,
    )
    return source.collect(run_context(tmp_path))


def _statement(collected: Collected, locus: str) -> dict[str, Any]:
    (statement,) = (
        r.statement
        for r in collected.rows
        if r.statement["statement_id"] == f"{PARALOGS_TOOL}:{locus}"
    )
    return statement


def _gene_call(collected: Collected, locus: str) -> dict[str, Any]:
    (call,) = (c for c in collected.extra["gene_call"] if c["gene"] == locus)
    return call


def test_strc_carrier_needs_unique_depth(tmp_path: Path) -> None:
    one = _collect(tmp_path, {"STRC": _ps(1, qual=40.0)}, unique_ratio=0.5)
    statement = _statement(one, "STRC")
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == str(Topic.CARRIER)
    assert (statement["reason"], statement["validity"]) == (STRC_REASON, "limited")
    assert statement["actionable"] is False
    assert _gene_call(one, "STRC")["level"] == "medium"
    assert _gene_call(one, "STRC")["diplotype"] == "STRC 1"
    high = _statement(_collect(tmp_path, {"STRC": _ps(1)}, unique_ratio=0.95), "STRC")
    assert high["kind"] == str(Kind.UNASSESSED)
    assert "CATSPER2" in high["reason"]
    two = _statement(_collect(tmp_path, {"STRC": _ps(2)}), "STRC")
    assert two["kind"] == str(Kind.NON_FINDING)
    low = _statement(_collect(tmp_path, {"STRC": _ps(1, qual=10.0)}), "STRC")
    assert low["kind"] == str(Kind.UNASSESSED)
    assert "quality" in low["reason"]


def test_strc_rows_disagree_unassessed(tmp_path: Path) -> None:
    collected = _collect(tmp_path, {"STRC": _ps(None, rows=2)}, unique_ratio=0.5)
    statement = _statement(collected, "STRC")
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "2 rows" in statement["reason"]
    assert _gene_call(collected, "STRC")["level"] == "low"
    failed = _collect(tmp_path, {"STRC": _ps(1, passed=False)}, unique_ratio=0.5)
    assert _statement(failed, "STRC")["kind"] == str(Kind.UNASSESSED)
    no_row = _collect(
        tmp_path, {"STRC": _ps(None, qual=None, rows=0)}, unique_ratio=0.5
    )
    assert "no row" in _statement(no_row, "STRC")["reason"]


def test_strc_zero_copies_unassessed(tmp_path: Path) -> None:
    statement = _statement(_collect(tmp_path, {"STRC": _ps(0)}), "STRC")
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert statement["reason"] == "no copy in a healthy adult: a calling error"


def test_ncf1_sord_not_covered(tmp_path: Path) -> None:
    collected = _collect(tmp_path)
    reason = (
        "their common carrier alleles, the GT deletion and c.757delG, change no copy"
        " number, so a copy-number rule cannot see them"
    )
    for gene in ("NCF1", "SORD"):
        (statement,) = (
            r.statement for r in collected.rows if r.statement["gene"] == gene
        )
        assert statement["kind"] == str(Kind.UNASSESSED)
        assert statement["reason"] == reason
        assert statement["topic"] == str(Topic.CARRIER)
        assert "copy_number" not in json.loads(str(statement["claim"]))
    genes = {c["gene"] for c in collected.extra["gene_call"]}
    assert not genes & {"NCF1", "SORD"}


def test_c4_total_only(tmp_path: Path) -> None:
    collected = _collect(tmp_path, {"C4": _ps(4)})
    statement = _statement(collected, "C4")
    assert statement["kind"] == str(Kind.FINDING)
    assert json.loads(str(statement["claim"]))["copy_number"] == 4
    ids = {r.statement["statement_id"] for r in collected.rows}
    genes = {r.statement["gene"] for r in collected.rows}
    assert not {"paralogs:C4A", "paralogs:C4B"} & ids
    assert not {"C4A", "C4B"} & genes


def test_rhd_copy_number_only(tmp_path: Path) -> None:
    for copies in (0, 1, 2):
        collected = _collect(tmp_path, {"RHD": _ps(copies)})
        statement = _statement(collected, "RHD")
        assert statement["kind"] == str(Kind.FINDING)
        assert statement["topic"] == str(Topic.TRAITS)
        assert statement["validity"] == "limited"
        assert statement["actionable"] is False
        assert statement["reason"] == RHD_NOTE
        assert json.loads(str(statement["claim"]))["copy_number"] == copies
    unassessed = _statement(_collect(tmp_path, {"RHD": _ps(None)}), "RHD")
    assert unassessed["kind"] == str(Kind.UNASSESSED)
    assert RHD_NOTE in unassessed["reason"]
    for withheld in (
        _collect(tmp_path, withheld={"RHD": "not covered: no call in calibration"}),
        _collect(tmp_path, not_run="the paralog callers did not run"),
    ):
        assert RHD_NOTE in _statement(withheld, "RHD")["reason"]
    for collected in (_collect(tmp_path, {"RHD": _ps(0)}), _collect(tmp_path)):
        for row in collected.rows:
            text = f"{row.statement['reason']} {row.statement['claim']}"
            for word in ("RhD negative", "RhD positive", "blood type"):
                assert word not in text


def _every_statement(tmp_path: Path) -> list[Fields]:
    ctx = run_context(tmp_path)
    collected = [
        _collect(tmp_path, {"STRC": _ps(1)}, unique_ratio=0.5, gauchian=CARRIER),
        _collect(tmp_path, gauchian=None, aldy_cyp2a6=1),
        _collect(tmp_path, not_run="the paralog callers did not run"),
        LpaSource(kiv2=24.3, withheld=None, not_run=None).collect(ctx),
        LpaSource(kiv2=None, withheld=None, not_run=None).collect(ctx),
    ]
    return [r.statement for c in collected for r in c.rows]


def test_no_watch_no_exclusion(tmp_path: Path) -> None:
    statements = _every_statement(tmp_path)
    assert statements
    for statement in statements:
        assert statement["kind"] != str(Kind.WATCH)
        assert statement["actionable"] is False
        # 3c's schema 7 adds `exclusion`; 3d leaves it at its default, false.
        assert not statement.get("exclusion", False)


def test_information_loci(tmp_path: Path) -> None:
    collected = _collect(tmp_path, {"AMY1": _ps(6), "CYP2A6": _ps(2)})
    amy1 = _statement(collected, "AMY1")
    assert (amy1["kind"], amy1["topic"]) == (str(Kind.FINDING), str(Topic.TRAITS))
    cyp2a6 = _statement(collected, "CYP2A6")
    assert (cyp2a6["kind"], cyp2a6["topic"]) == (str(Kind.FINDING), str(Topic.DRUGS))
    assert cyp2a6["reason"] == "no CPIC guideline; information"
    assert _gene_call(collected, "CYP2A6")["agreed"] is True
    # UGT2B17 is not covered (pre-flight): no copy-number locus of its own.
    assert "UGT2B17" not in {c["gene"] for c in collected.extra["gene_call"]}


def test_cyp2a6_disagrees_with_aldy(tmp_path: Path) -> None:
    collected = _collect(tmp_path, {"CYP2A6": _ps(2)}, aldy_cyp2a6=1)
    statement = _statement(collected, "CYP2A6")
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "Parascopy and Aldy disagree" in statement["reason"]
    assert _gene_call(collected, "CYP2A6")["agreed"] is False
    none = _statement(_collect(tmp_path, aldy_cyp2a6=None), "CYP2A6")
    assert none["kind"] == str(Kind.UNASSESSED)


def test_gba_wild_type_non_finding(tmp_path: Path) -> None:
    statement = _statement(_collect(tmp_path), "GBA1")
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert statement["topic"] == str(Topic.CARRIER)


def test_gba_biallelic_withheld(tmp_path: Path) -> None:
    biallelic = GauchianCall(4, True, True, (L483R,), (), 9)
    statement = _statement(_collect(tmp_path, gauchian=biallelic), "GBA1")
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert statement["reason"] == BIALLELIC_WITHHELD


def test_gba_no_parascopy_unassessed(tmp_path: Path) -> None:
    for call in (_ps(None), _ps(3)):
        collected = _collect(tmp_path, {"GBA1": call}, gauchian=CARRIER)
        statement = _statement(collected, "GBA1")
        assert statement["kind"] == str(Kind.UNASSESSED)
        assert "Parascopy" in statement["reason"]
    assert _gene_call(collected, "GBA1")["agreed"] is False
    gone = ParalogSource(
        calls={},
        unique_ratio=None,
        gauchian=CARRIER,
        aldy_cyp2a6=None,
        withheld={},
        not_run=None,
    ).collect(run_context(tmp_path))
    assert _statement(gone, "GBA1")["kind"] == str(Kind.UNASSESSED)


def test_gba_carrier_finding(tmp_path: Path) -> None:
    floor = run_context(tmp_path).grading.paralogs.gba_min_reads
    at_floor = GauchianCall(4, True, False, (L483R,), (), floor)
    collected = _collect(tmp_path, {"GBA1": _ps(4)}, gauchian=at_floor)
    statement = _statement(collected, "GBA1")
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == str(Topic.CARRIER)
    assert (statement["reason"], statement["validity"]) == (GBA_REASON, "limited")
    assert statement["actionable"] is False
    assert _gene_call(collected, "GBA1")["level"] == "medium"
    assert _gene_call(collected, "GBA1")["agreed"] is True
    below = GauchianCall(4, True, False, (L483R,), (), floor - 1)
    low = _statement(_collect(tmp_path, gauchian=below), "GBA1")
    assert low["kind"] == str(Kind.UNASSESSED)
    assert "reads" in low["reason"]
    unknown = GauchianCall(4, True, False, (L483R,), (), None)
    none = _statement(_collect(tmp_path, gauchian=unknown), "GBA1")
    assert none["kind"] == str(Kind.UNASSESSED)
    assert "reads" in none["reason"]


def test_gauchian_unhandled_unassessed(tmp_path: Path) -> None:
    for call in (
        GauchianCall(3, False, False, (), (), None),
        GauchianCall(4, False, False, (), ("c.1A>T",), None),
        GauchianCall(4, None, False, (), (), None),
        GauchianCall(4, True, False, (L483R,), ("c.1A>T",), 7),
        GauchianCall(None, True, False, (L483R,), (), 7),
        GauchianCall(4, False, False, (), (), None, "intron9"),
        GauchianCall(4, True, False, (L483R,), (), 7, "intron9"),
    ):
        statement = _statement(_collect(tmp_path, gauchian=call), "GBA1")
        assert statement["kind"] == str(Kind.UNASSESSED), call
        assert statement["reason"].startswith("Gauchian output not handled: "), call


def _clinvar_at(pos: int, *, source: str = "clinvar") -> Fields:
    return {
        "statement_id": "clinvar:1",
        "kind": str(Kind.FINDING),
        "vid": f"GRCh38:chr1:{pos}:C:A",
        "gene": "GBA1",
        "topic": str(Topic.CARRIER),
        "source": source,
    }


def test_gba1_exon_call_contradicted() -> None:
    exon10 = _clinvar_at(155235252)  # c.1448T>G, GBA1 exon 10
    assert gba1_contradicted(exon10, WILD, EXONS_9_11)
    assert gba1_contradicted(_clinvar_at(155235252, source="vep"), WILD, EXONS_9_11)
    assert not gba1_contradicted(exon10, CARRIER, EXONS_9_11)
    named_other = GauchianCall(4, False, False, (), (L483R,), None)
    assert not gba1_contradicted(exon10, named_other, EXONS_9_11)
    # A variant Gauchian names elsewhere in exons 9-11 does not cover this position.
    elsewhere = GauchianCall(4, True, False, ("c.1604G>A(p.Arg535His)",), (), 6)
    assert gba1_contradicted(exon10, elsewhere, EXONS_9_11)
    exon8 = _clinvar_at(155236300)
    assert not gba1_contradicted(exon8, WILD, EXONS_9_11)
    assert not gba1_contradicted(exon10, None, EXONS_9_11)
    assert not gba1_contradicted(_clinvar_at(155235252, source="smn"), WILD, EXONS_9_11)


def test_withheld_locus_unassessed(tmp_path: Path) -> None:
    reason = "not covered: HG002's libraries disagree"
    collected = _collect(tmp_path, withheld={"AMY1": reason, "GBA1": reason})
    assert _statement(collected, "AMY1")["kind"] == str(Kind.UNASSESSED)
    assert _statement(collected, "AMY1")["reason"] == reason
    assert _statement(collected, "GBA1")["reason"] == reason
    assert _statement(collected, "C4")["kind"] == str(Kind.FINDING)
    not_run = _collect(tmp_path, not_run="the paralog callers did not run")
    for locus in ("STRC", "RHD", "AMY1", "C4", "CYP2A6", "GBA1"):
        statement = _statement(not_run, locus)
        assert statement["kind"] == str(Kind.UNASSESSED)
        assert "the paralog callers did not run" in statement["reason"]


def test_not_covered_statements(tmp_path: Path) -> None:
    collected = _collect(tmp_path)
    assert set(NOT_COVERED) == {
        "HBA1", "HBA2", "CYP21A2", "TNXB", "PMS2", "F8", "OTOA", "GSTT1", "GSTM1",
        "NCF1", "SORD", "UGT2B17",
    }  # fmt: skip
    assert NOT_COVERED["UGT2B17"] == (
        Topic.TRAITS,
        "Parascopy's UGT2B15 model covers 21% of UGT2B17's MANE span",
    )
    for gene, (topic, reason) in NOT_COVERED.items():
        (statement,) = (
            r.statement for r in collected.rows if r.statement["gene"] == gene
        )
        assert statement["kind"] == str(Kind.UNASSESSED)
        assert (statement["topic"], statement["reason"]) == (str(topic), reason)
        assert statement["source"] == PARALOGS_TOOL
    assert collected.input_rows == len(collected.rows) == 6 + len(NOT_COVERED)
    assert len(collected.extra["gene_call"]) == 6


def test_lpa_information(tmp_path: Path) -> None:
    ctx = run_context(tmp_path)
    found = LpaSource(kiv2=24.34, withheld=None, not_run=None).collect(ctx)
    assert found.tool == LPA_TOOL
    (row,) = found.rows
    statement = row.statement
    assert statement["statement_id"] == "lpa:LPA"
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == str(Topic.WELLNESS)
    assert statement["validity"] == "limited"
    assert statement["reason"] == LPA_REASON
    assert json.loads(str(statement["claim"]))["kiv2"] == 24.3
    reason = "not covered: KILDA off HG002's count by 9.0 (30%)"
    (withheld,) = LpaSource(kiv2=24.3, withheld=reason, not_run=None).collect(ctx).rows
    assert withheld.statement["kind"] == str(Kind.UNASSESSED)
    assert withheld.statement["reason"] == reason
    assert json.loads(str(withheld.statement["claim"]))["kiv2"] is None
    (no_call,) = LpaSource(kiv2=None, withheld=None, not_run=None).collect(ctx).rows
    assert no_call.statement["kind"] == str(Kind.UNASSESSED)
