"""The HLA statements (`sources.hla.HlaSource`; slice 3f, "HLA"): class I alleles where
OptiType and T1K agree, class II research rows from T1K alone, the coeliac rule-out where T1K
and the six DQ tag SNPs agree, the class I/II table, and KIR named as not covered.

The callers' calls start from HG002's (`fixtures/hg002/`), with single genes changed."""

import json
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from fixtures.run_context import run_context

from genokeeper.pgx.hla import HlaCall, read_optitype, read_t1k
from genokeeper.pgx.hla_typing import (
    ASSOCIATION_NOTE,
    B27_NOTE,
    CONFIRM_DQ,
    EXCLUDED,
    TAG_RSIDS,
    TAGS_EUROPEAN,
    DqStatus,
    carried_reason,
)
from genokeeper.pgx.pharmcat import DOSING
from genokeeper.query import Filters, query
from genokeeper.sources import INFORMATION_NOTE, write_collected
from genokeeper.sources.hla import (
    CLASS_II_ONE_CALLER,
    HLA_ASSOCIATION_TOOL,
    HLA_TOOL,
    KIR_REASON,
    HlaSource,
)
from genokeeper.store.db import Kind, Topic, Validity, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

HG002 = Path(__file__).parent / "fixtures" / "hg002"
OPT = read_optitype(HG002 / "hg002.optitype_result.tsv.txt")
T1K = read_t1k(HG002 / "hg002.t1k_genotype.tsv.txt", genes=None)
NO_TAGS = {rsid: "0" for rsid in TAG_RSIDS}
# HG002's DQ8 tag (rs7454108 1|0) beside its homozygous rs2395182.
HG002_TAGS = NO_TAGS | {"rs2395182": "2", "rs7454108": "1"}
# T1K's class II with neither DQ2, DQ8 nor a half-heterodimer.
CLEAN_DQ = {"HLA-DQA1": ("*01:01", "*01:02"), "HLA-DQB1": ("*05:01", "*06:03")}


def collect(
    tmp_path: Path,
    *,
    opt: Mapping[str, HlaCall] | None = None,
    t1k: Mapping[str, HlaCall] | None = None,
    tags: dict[str, str] | None = None,
    not_run: str | None = None,
    withheld: dict[str, str] | None = None,
    raw: tuple[Mapping[str, HlaCall] | None, Mapping[str, HlaCall] | None]
    | None = None,
) -> dict[str, dict[str, Any]]:
    """Each statement by id, its claim parsed; the callers are HG002's with `opt` and `t1k`
    laid over them (`raw` instead: the two mappings as given)."""
    optitype, typed = (
        raw if raw is not None else ({**OPT, **(opt or {})}, {**T1K, **(t1k or {})})
    )
    source = HlaSource(
        optitype,
        typed,
        HG002_TAGS if tags is None else tags,
        not_run,
        withheld or {},
    )
    found = source.collect(run_context(tmp_path))
    assert found.tool == HLA_TOOL
    assert found.input_rows == len(found.rows)
    statements: dict[str, dict[str, Any]] = {}
    for row in found.rows:
        assert row.variant is None and row.call is None
        s = dict(row.statement)
        s["claim"] = json.loads(str(s["claim"]))
        statements[str(s["statement_id"])] = s
    return statements


def test_class_i_allele_needs_both_callers(tmp_path: Path) -> None:
    both = {"HLA-B": ("*27:05", "*35:08")}
    s = collect(tmp_path, opt=both, t1k=both)["hla:HLA-B*27"]
    assert (s["kind"], s["topic"], s["validity"], s["reason"]) == (
        "finding",
        "risks",
        "moderate",
        B27_NOTE,
    )
    assert s["actionable"] is False
    # T1K's equivalents all in the group: still positive.
    s = collect(tmp_path, opt=both, t1k={"HLA-B": ("*27:05,*27:02", "*35:08")})
    assert s["hla:HLA-B*27"]["kind"] == "finding"
    s = collect(tmp_path, opt=both)["hla:HLA-B*27"]
    assert s["kind"] == "unassessed"
    assert "OptiType positive" in s["reason"] and "T1K negative" in s["reason"]
    # HG002 carries no B*27: both callers negative.
    assert collect(tmp_path)["hla:HLA-B*27"]["kind"] == "non_finding"


def test_drug_allele_wording(tmp_path: Path) -> None:
    both = {"HLA-B": ("*13:01", "*35:08")}
    s = collect(tmp_path, opt=both, t1k=both)["hla:HLA-B*13:01"]
    assert (s["kind"], s["topic"], s["validity"]) == ("finding", "drugs", "limited")
    assert s["claim"]["drugs"] == ["dapsone", "co-trimoxazole"]
    assert s["actionable"] is False
    assert not DOSING.search(json.dumps(s))
    # HG002's C*04:01, the 2B nevirapine allele, agreed by both callers.
    c04 = collect(tmp_path)["hla:HLA-C*04:01"]
    assert (c04["kind"], c04["claim"]["drugs"]) == ("finding", ["nevirapine"])


def test_pharmcat_alleles_left_out(tmp_path: Path) -> None:
    ids = set(collect(tmp_path))
    for allele in ("HLA-B*15:02", "HLA-B*57:01", "HLA-B*58:01", "HLA-A*31:01"):
        assert f"hla:{allele}" not in ids


def coeliac(tmp_path: Path, t1k: Mapping[str, HlaCall], tags: dict[str, str]):
    return collect(tmp_path, t1k=t1k, tags=tags)["hla:coeliac"]


def test_coeliac_rule_out_needs_both(tmp_path: Path) -> None:
    s = coeliac(tmp_path, CLEAN_DQ, NO_TAGS)
    assert (s["kind"], s["reason"], s["exclusion"]) == (
        "non_finding",
        f"{EXCLUDED}; {TAGS_EUROPEAN}",
        True,
    )
    assert s["claim"]["confirm"] == CONFIRM_DQ == "HLA-DQ typing in a laboratory"
    dq25 = {"HLA-DQA1": ("*05:01", "*01:01"), "HLA-DQB1": ("*02:01", "*05:01")}
    s = coeliac(tmp_path, dq25, NO_TAGS | {"rs2187668": "1"})
    assert (s["kind"], s["topic"], s["validity"], s["exclusion"]) == (
        "finding",
        "risks",
        "limited",
        False,
    )
    assert s["reason"] == f"{carried_reason(DqStatus(1, 0, 0, True))}; {TAGS_EUROPEAN}"
    assert s["reason"].startswith("possibly carries HLA-DQ2:")
    half = {"HLA-DQA1": ("*05:05", "*01:01"), "HLA-DQB1": ("*03:01", "*05:01")}
    s = coeliac(tmp_path, half, NO_TAGS)
    assert (s["kind"], s["exclusion"]) == ("unassessed", False)


def test_half_heterodimer_finding_wording(tmp_path: Path) -> None:
    half = {"HLA-DQA1": ("*05:05", "*01:01"), "HLA-DQB1": ("*03:01", "*05:01")}
    s = coeliac(tmp_path, half, NO_TAGS | {"rs4639334": "1"})
    assert s["kind"] == "finding"
    assert s["reason"] == f"{carried_reason(DqStatus(0, 0, 0, True))}; {TAGS_EUROPEAN}"
    assert "DQ2 or" not in s["reason"] and "DQ8" not in s["reason"]
    dq8 = {"HLA-DQA1": ("*03:01", "*01:01"), "HLA-DQB1": ("*03:02", "*05:01")}
    s = coeliac(tmp_path, dq8, NO_TAGS | {"rs7454108": "1"})
    assert (s["kind"], s["reason"]) == (
        "finding",
        f"{carried_reason(DqStatus(0, 0, 1, False))}; {TAGS_EUROPEAN}",
    )
    # HG002 itself: DQ8 heterozygous by both.
    assert collect(tmp_path)["hla:coeliac"]["reason"] == s["reason"]


def test_tags_see_half_t1k_misses_it(tmp_path: Path) -> None:
    none = {"HLA-DQA1": ("*01:01", "*01:02"), "HLA-DQB1": ("*05:01", "*06:02")}
    s = coeliac(tmp_path, none, NO_TAGS | {"rs4639334": "1"})
    assert s["kind"] == "unassessed"
    assert EXCLUDED not in s["reason"] and s["exclusion"] is False


def test_dropped_class_ii_allele_never_excludes(tmp_path: Path) -> None:
    """Review Focus 2: OptiType and T1K agree on class I, but T1K dropped a DQB1 allele
    (quality 0): the coeliac statement is unassessed, never "practically excluded"."""
    text = (HG002 / "hg002.t1k_genotype.tsv.txt").read_text()
    clean = (
        "HLA-DQA1\t2\tHLA-DQA1*01:01:01\t65.1\t40\tHLA-DQA1*01:02:01\t63.7\t39\t\n"
        "HLA-DQB1\t2\tHLA-DQB1*05:01:01\t43.1\t30\tHLA-DQB1*06:03:01\t34.5\t0\t\n"
    )
    lines = [
        line
        for line in text.splitlines(keepends=True)
        if not line.startswith(("HLA-DQA1\t", "HLA-DQB1\t"))
    ]
    path = tmp_path / "t1k_genotype.tsv"
    path.write_text("".join(lines) + clean)
    typed = read_t1k(path, genes=None)
    assert typed["HLA-DQB1"] == ("*05:01", None)
    s = collect(tmp_path, raw=(OPT, typed), tags=NO_TAGS)
    assert s["hla:HLA-B*27"]["kind"] == "non_finding"  # class I agreed
    rule = s["hla:coeliac"]
    assert (rule["kind"], rule["exclusion"]) == ("unassessed", False)
    assert EXCLUDED not in rule["reason"]
    # T1K's DQB1 missing altogether, or T1K not run at all: the same.
    del typed["HLA-DQB1"]
    for raw in ((OPT, typed), (OPT, None)):
        rule = collect(tmp_path, raw=raw, tags=NO_TAGS)["hla:coeliac"]
        assert (rule["kind"], rule["exclusion"]) == ("unassessed", False)


def test_dq_type_disagreement_unassessed(tmp_path: Path) -> None:
    dq25 = {"HLA-DQA1": ("*05:01", "*01:01"), "HLA-DQB1": ("*02:01", "*05:01")}
    s = coeliac(tmp_path, dq25, NO_TAGS | {"rs7454108": "1"})
    assert s["kind"] == "unassessed"
    assert "DQ2" not in s["reason"] and "DQ8" not in s["reason"]


def test_equivalents_spanning_risk_unknown(tmp_path: Path) -> None:
    for t1k in (
        {"HLA-DQA1": ("*03:01", "*01:01"), "HLA-DQB1": ("*03:02,*03:01", "*05:01")},
        {"HLA-DQA1": ("*05:05,*06:01", "*01:01"), "HLA-DQB1": ("*05:01", "*06:02")},
    ):
        for tags in (NO_TAGS, HG002_TAGS, NO_TAGS | {"rs4639334": "1"}):
            assert coeliac(tmp_path, t1k, tags)["kind"] == "unassessed"


def test_hla_associations_research(tmp_path: Path) -> None:
    c06 = {"HLA-C": ("*06:02", "*04:01")}
    found = collect(tmp_path, opt=c06, t1k=c06)
    s = found["hla:HLA-C*06:02:association"]
    assert (s["kind"], s["topic"], s["validity"], s["actionable"]) == (
        "finding",
        "research",
        "limited",
        False,
    )
    assert s["reason"] == ASSOCIATION_NOTE
    assert s["claim"]["type"] == "hla_association"
    assert "10.4238/2011.December.15.2" in s["source_id"]
    # The same allele's co-trimoxazole statement is a separate row (ruling R3).
    assert found["hla:HLA-C*06:02"]["topic"] == "drugs"
    b51 = {"HLA-B": ("*51:01", "*35:08")}
    s = collect(tmp_path, opt=b51, t1k=b51)["hla:HLA-B*51:association"]
    assert (s["kind"], s["topic"], s["reason"]) == (
        "finding",
        "research",
        ASSOCIATION_NOTE,
    )
    s = collect(tmp_path, opt=b51)["hla:HLA-B*51:association"]
    assert s["kind"] == "unassessed"
    assert "OptiType positive" in s["reason"] and "T1K negative" in s["reason"]


def test_hla_associations_are_a_research_layer(tmp_path: Path) -> None:
    """Ruling R19: the associations' source is a research layer's, hidden from a query
    that names no gene, kept for `gene=HLA-C`."""
    con = connect(tmp_path / "s.duckdb")
    run_id = start_run(con, {}, {})
    c06 = {**T1K, "HLA-C": ("*06:02", "*04:01")}
    found = HlaSource({**OPT, **c06}, c06, HG002_TAGS, None, {}).collect(
        replace(run_context(tmp_path), run_id=run_id)
    )
    with transaction(con):
        write_collected(StoreWriter(con, tmp_path / "w"), con, run_id, [found])
    loose = Filters(min_validity=Validity.LIMITED)
    sid = "hla:HLA-C*06:02:association"
    assert sid not in {x["id"] for x in query(con, run_id, loose).findings}
    narrowed = query(con, run_id, replace(loose, gene="HLA-C")).findings
    kept = {x["id"]: x for x in narrowed}[sid]
    assert kept["sources"][0]["source"] == HLA_ASSOCIATION_TOOL


@pytest.mark.parametrize(
    ("t1k", "tags"),
    [
        (CLEAN_DQ, NO_TAGS),
        ({}, HG002_TAGS),
        ({"HLA-DQB1": ("*02:01", "*05:01")}, NO_TAGS),
    ],
)
def test_only_the_rule_out_is_flagged(
    tmp_path: Path, t1k: Mapping[str, HlaCall], tags: dict[str, str]
) -> None:
    for sid, s in collect(tmp_path, t1k=t1k, tags=tags).items():
        assert isinstance(s["exclusion"], bool)
        assert s["exclusion"] is (sid == "hla:coeliac" and s["kind"] == "non_finding")
        assert s["kind"] != "watch" and s["actionable"] is False
        assert s["source"] == (
            HLA_ASSOCIATION_TOOL if sid.endswith(":association") else HLA_TOOL
        )


def test_withheld_and_not_run(tmp_path: Path) -> None:
    reason = "not covered: OptiType and T1K differ from HG002's truth on HLA-B"
    found = collect(tmp_path, withheld={"HLA-B": reason})
    for sid, s in found.items():
        if sid.startswith("hla:HLA-B"):
            assert (s["kind"], s["reason"]) == ("unassessed", reason)
    assert found["hla:HLA-C*04:01"]["kind"] == "finding"
    dq = collect(tmp_path, t1k=CLEAN_DQ, tags=NO_TAGS, withheld={"DQ": "tags"})
    assert (dq["hla:coeliac"]["kind"], dq["hla:coeliac"]["reason"]) == (
        "unassessed",
        "tags",
    )
    not_run = "not calibrated: run `genokeeper calibrate`"
    found = collect(tmp_path, not_run=not_run)
    for sid, s in found.items():
        assert s["kind"] == "unassessed"
        assert s["reason"] == (KIR_REASON if sid == "hla:KIR" else not_run)
    for run in (collect(tmp_path), found, collect(tmp_path, raw=(None, None))):
        kir = run["hla:KIR"]
        assert (kir["kind"], kir["reason"]) == ("unassessed", KIR_REASON)


def test_callers_missing_unassessed(tmp_path: Path) -> None:
    for sid, s in collect(tmp_path, raw=(None, None)).items():
        assert s["kind"] == "unassessed", sid
    s = collect(tmp_path, raw=(None, T1K))["hla:HLA-C*04:01"]
    assert (s["kind"], s["reason"]) == ("unassessed", "OptiType made no call")


def test_table_is_research(tmp_path: Path) -> None:
    found = collect(tmp_path)
    s = found["hla:table"]
    assert (s["kind"], s["topic"], s["validity"], s["exclusion"]) == (
        "finding",
        "research",
        "limited",
        False,
    )
    genes = {g["gene"]: g for g in s["claim"]["genes"]}
    assert genes["HLA-A"]["alleles"] == ["*01:01", "*26:01"]
    for gene in ("HLA-A", "HLA-B", "HLA-C"):
        assert genes[gene]["callers"] == "OptiType and T1K"
    for gene in ("HLA-DQA1", "HLA-DQB1", "HLA-DRB1"):
        assert (genes[gene]["callers"], genes[gene]["note"]) == (
            "T1K",
            CLASS_II_ONE_CALLER,
        )
    assert genes["HLA-DQB1"]["alleles"] == ["*03:02", "*05:01"]
    # The class II genes no calibration measure covers are named as information.
    uncalibrated = "HLA-DRB4, HLA-DPA1, HLA-DPB1"
    assert s["reason"] == f"{CLASS_II_ONE_CALLER}; {uncalibrated}: {INFORMATION_NOTE}"
    # A class I gene the callers disagree on is listed without alleles.
    s = collect(tmp_path, opt={"HLA-B": ("*35:08", "*38:02")})["hla:table"]
    b = {g["gene"]: g for g in s["claim"]["genes"]}["HLA-B"]
    assert b["alleles"] is None and "*38:02" in b["callers"]
    # DQB1*06:02 is a research row from T1K alone.
    s = collect(tmp_path, t1k={"HLA-DQB1": ("*06:02", "*05:01")})
    row = s["hla:HLA-DQB1*06:02:association"]
    assert (row["kind"], row["topic"], row["validity"]) == (
        "finding",
        str(Topic.RESEARCH),
        str(Validity.LIMITED),
    )
    assert CLASS_II_ONE_CALLER in row["reason"]
    assert found["hla:HLA-DQB1*06:02:association"]["kind"] == str(Kind.NON_FINDING)
