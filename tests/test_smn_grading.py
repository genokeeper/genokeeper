"""SMN1 grading (`sources.smn`): the statement and gene call from the two callers' results
built by hand."""

import json
from pathlib import Path
from typing import Any

from fixtures.run_context import run_context

from genokeeper.evidence.pgx import GENE_UNRESOLVED
from genokeeper.smn import ParascopyCall, SmnCaller
from genokeeper.sources.smn import SmnSource
from genokeeper.store.db import Kind


def _caller(
    smn1: int | None, *, smn2: int | None = 2, marker: int | None = 0
) -> SmnCaller:
    return SmnCaller(smn1, smn2, marker, "high")


def _parascopy(
    smn1: int | None, *, qual: float = 40.0, passed: bool = True
) -> ParascopyCall:
    return ParascopyCall(smn1, 2, qual, passed)


def _grade(
    tmp_path: Path,
    caller: SmnCaller | None,
    parascopy: ParascopyCall | None,
    *,
    not_run: str | None = None,
    withheld: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    source = SmnSource(caller, parascopy, not_run, withheld, 20.0)
    collected = source.collect(run_context(tmp_path))
    return collected.rows[0].statement, collected.extra["gene_call"][0]


def test_one_copy_both_agree_finding(tmp_path: Path) -> None:
    statement, gene_call = _grade(tmp_path, _caller(1), _parascopy(1))
    assert statement["statement_id"] == "smn:SMN1"
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == "carrier"
    assert statement["validity"] == "strong"
    assert statement["actionable"] is False
    assert "screening result" in statement["reason"]
    assert "MLPA or qPCR" in statement["reason"]
    assert gene_call["level"] == "medium"
    assert gene_call["agreed"] is True
    assert gene_call["diplotype"] == "SMN1 1, SMN2 2"
    assert set(json.loads(gene_call["callers"])) == {
        "SMNCopyNumberCaller 1.1.2",
        "Parascopy 1.19.0",
    }


def test_two_copies_reduced_risk(tmp_path: Path) -> None:
    statement, _ = _grade(tmp_path, _caller(2), _parascopy(2))
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert "2+0 carrier" in statement["reason"]
    three, _ = _grade(tmp_path, _caller(3), _parascopy(3))
    assert three["kind"] == str(Kind.NON_FINDING)


def test_silent_carrier_marker_watch(tmp_path: Path) -> None:
    statement, _ = _grade(tmp_path, _caller(2, marker=1), _parascopy(2))
    assert statement["kind"] == str(Kind.WATCH)
    assert "possible silent carrier" in statement["reason"]
    assert "screening result, needs confirmation" in statement["reason"]


def test_smn_disagreement_unassessed(tmp_path: Path) -> None:
    statement, gene_call = _grade(tmp_path, _caller(1), _parascopy(2))
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "disagree" in statement["reason"]
    assert gene_call["agreed"] is False


def test_parascopy_low_quality_unassessed(tmp_path: Path) -> None:
    low, _ = _grade(tmp_path, _caller(1), _parascopy(1, qual=19.0))
    assert low["kind"] == str(Kind.UNASSESSED)
    failed, _ = _grade(tmp_path, _caller(1), _parascopy(1, passed=False))
    assert failed["kind"] == str(Kind.UNASSESSED)
    edge, _ = _grade(tmp_path, _caller(1), _parascopy(1, qual=20.0))
    assert edge["kind"] == str(Kind.FINDING)


def test_zero_copies_unassessed(tmp_path: Path) -> None:
    statement, _ = _grade(tmp_path, _caller(0), _parascopy(0))
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "calling error" in statement["reason"]


def test_not_run_and_withheld_unassessed(tmp_path: Path) -> None:
    not_run, call = _grade(tmp_path, None, None, not_run="stage 3 not calibrated")
    assert not_run["kind"] == str(Kind.UNASSESSED)
    assert not_run["reason"] == "stage 3 not calibrated"
    # `query` refuses an unassessed gene call whose warning lacks the prefix.
    assert call["warnings"] == [f"{GENE_UNRESOLVED}: stage 3 not calibrated"]
    withheld, _ = _grade(
        tmp_path, _caller(1), _parascopy(1), withheld="false call in calibration"
    )
    assert withheld["kind"] == str(Kind.UNASSESSED)
    assert withheld["reason"] == "false call in calibration"
    none, _ = _grade(tmp_path, _caller(None), _parascopy(1))
    assert none["kind"] == str(Kind.UNASSESSED)
