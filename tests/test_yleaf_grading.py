"""Y haplogroup grading (`sources.yleaf`): the statement from Yleaf's result built by hand."""

import json
from pathlib import Path
from typing import Any

from fixtures.run_context import run_context

from genokeeper.sources import Collected
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.sources.yleaf import YleafSource
from genokeeper.store.db import Kind
from genokeeper.yhaplogroup import YHaplogroup


def _collect(
    tmp_path: Path,
    found: YHaplogroup | None,
    *,
    not_run: str | None = None,
    withheld: str | None = None,
    male: bool = True,
) -> Collected:
    source = YleafSource(found, not_run, withheld, 0.95, male)
    return source.collect(run_context(tmp_path))


def _statement(collected: Collected) -> dict[str, Any]:
    (row,) = collected.rows
    return row.statement


def test_haplogroup_finding(tmp_path: Path) -> None:
    strong = _statement(_collect(tmp_path, YHaplogroup("R-L21", "L21", 0.995, 120)))
    assert strong["kind"] == str(Kind.FINDING)
    assert strong["statement_id"] == "yleaf:haplogroup"
    assert (strong["vid"], strong["gene"], strong["topic"]) == (None, None, "ancestry")
    assert strong["validity"] == "strong"
    assert strong["actionable"] is False
    assert json.loads(strong["claim"]) == {
        "schema": CLAIM_SCHEMA,
        "type": "y_haplogroup",
        "haplogroup": "R-L21",
        "marker": "L21",
        "qc_score": 0.995,
        "valid_markers": 120,
        "tree": "Yleaf 4.1.6 bundled YFull",
    }
    moderate = _statement(_collect(tmp_path, YHaplogroup("R-L21", "L21", 0.95, 120)))
    assert (moderate["kind"], moderate["validity"]) == ("finding", "moderate")


def test_low_qc_unassessed(tmp_path: Path) -> None:
    statement = _statement(_collect(tmp_path, YHaplogroup("R-L21", "L21", 0.94, 40)))
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "0.94" in statement["reason"]
    assert json.loads(statement["claim"])["haplogroup"] is None


def test_female_no_statement(tmp_path: Path) -> None:
    collected = _collect(tmp_path, YHaplogroup("R-L21", "L21", 0.99, 120), male=False)
    assert collected.rows == []
    assert collected.input_rows == 0


def test_not_run_unassessed(tmp_path: Path) -> None:
    not_run = _statement(_collect(tmp_path, None, not_run="stage 3 not calibrated"))
    assert not_run["kind"] == str(Kind.UNASSESSED)
    assert not_run["reason"] == "stage 3 not calibrated"
    found = YHaplogroup("R-L21", "L21", 0.99, 120)
    withheld = _statement(_collect(tmp_path, found, withheld="L816 not derived"))
    assert withheld["kind"] == str(Kind.UNASSESSED)
    assert _statement(_collect(tmp_path, None))["kind"] == str(Kind.UNASSESSED)
