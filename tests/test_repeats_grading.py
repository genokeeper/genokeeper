"""Repeat grading (`genokeeper.repeats.loci`, `sources.repeats`): the thresholds, and the
statement and gene call of a locus from ExpansionHunter calls built by hand."""

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from fixtures.run_context import run_context

from genokeeper.config import load_grading
from genokeeper.repeats.expansionhunter import EH_VERSION, read_calls
from genokeeper.repeats.loci import (
    Allele,
    RepeatCall,
    RepeatClass,
    RepeatLocus,
    Rfc1Call,
    Rule,
    classify,
)
from genokeeper.repeats.not_covered import NOT_COVERED
from genokeeper.repeats.rfc1 import RFC1_PAIRS
from genokeeper.smn import ParascopyCall, SmnCaller
from genokeeper.sources import Collected
from genokeeper.sources.repeats import (
    BETWEEN,
    BEYOND_FRAGMENT,
    CFTR_UNREAD,
    FINDING_REASON,
    FINDING_REASONS,
    NO_SPANNING,
    NOT_SET,
    RFC1_NORMAL,
    RFC1_NOT_RUN,
    RFC1_UNCLEAR,
    TWO_ALLELES_ON_X,
    X_DOSAGE,
    X_NOT_CALIBRATED,
    RepeatSource,
    XDosage,
    cftr_grade,
    grade,
    rfc1_grade,
)
from genokeeper.sources.smn import SmnSource
from genokeeper.store.db import Kind, Level, Validity

LOCI = load_grading().repeats.loci
HTT, AR, ATXN7 = LOCI["HTT"], LOCI["AR"], LOCI["ATXN7"]
# The run's X dosage under a calibrated [dosage]: one X with a Y, two X without one.
ONE_X = XDosage(1, 0.5, True)
TWO_X = XDosage(2, 1.0, False)


def _allele(
    repeats: int,
    ci: tuple[int, int] | None = None,
    *,
    spanning: int = 20,
    flanking: int = 0,
    inrepeat: int = 0,
) -> Allele:
    return Allele(repeats, ci or (repeats, repeats), spanning, flanking, inrepeat)


def _call(
    locus: RepeatLocus,
    *alleles: Allele,
    low_depth: bool = False,
    longest_flanking: int = 0,
) -> RepeatCall:
    return RepeatCall(locus.variant_id, alleles, 30.0, low_depth, 150, longest_flanking)


def _source(*calls: RepeatCall, **kw: Any) -> RepeatSource:
    fields: dict[str, Any] = {
        "calls": {c.variant_id: c for c in calls},
        "not_run": None,
        "loci": LOCI,
        "withheld": {},
        "caps": {},
        "images": {},
        "pcr_like": False,
        "insert": 400.0,
        "male": True,
        "min_allele_reads": 5,
        "min_spanning": 3,
        "x_dosage": ONE_X,
        "rfc1_not_run": None,
    }
    return RepeatSource(**(fields | kw))


def _collect(tmp_path: Path, source: RepeatSource) -> Collected:
    return source.collect(run_context(tmp_path))


def _statement(collected: Collected, gene: str) -> dict[str, Any]:
    (row,) = [
        r.statement
        for r in collected.rows
        if r.statement["statement_id"] == f"expansionhunter:{gene}"
    ]
    return row


def _gene_call(collected: Collected, gene: str) -> dict[str, Any]:
    (row,) = [r for r in collected.extra["gene_call"] if r["gene"] == gene]
    return row


def test_classify_ranges() -> None:
    assert classify(HTT, 26) is RepeatClass.NORMAL
    assert classify(HTT, 27) is RepeatClass.INTERMEDIATE
    assert classify(HTT, 35) is RepeatClass.INTERMEDIATE
    assert classify(HTT, 36) is RepeatClass.REDUCED_PENETRANCE
    assert classify(HTT, 39) is RepeatClass.REDUCED_PENETRANCE
    assert classify(HTT, 40) is RepeatClass.PATHOGENIC
    assert classify(AR, 35) is RepeatClass.GAP
    assert classify(ATXN7, 36) is RepeatClass.GAP


def test_expanded_lower_bound_is_finding(tmp_path: Path) -> None:
    call = _call(HTT, _allele(18), _allele(42, (41, 44), spanning=8))
    collected = _collect(tmp_path, _source(call))
    statement = _statement(collected, "HTT")
    assert statement["kind"] == str(Kind.FINDING)
    assert "needs confirmation" in statement["reason"]
    assert statement["actionable"] is False
    assert statement["topic"] == "risks"
    gene_call = _gene_call(collected, "HTT")
    assert gene_call["diplotype"] == "CAG 18/42"
    assert gene_call["phenotype"] == "expanded (screening)"
    callers = json.loads(gene_call["callers"])
    assert callers == {
        f"ExpansionHunter {EH_VERSION}": "18/42, CI 18-18/41-44, reads 20/0/0; 8/0/0"
    }
    claim = json.loads(statement["claim"])
    assert claim["type"] == "repeat_expansion"
    assert claim["class"] == "pathogenic"
    assert claim["thresholds"]["pathogenic_min"] == 40


def test_interval_crossing_threshold_is_watch(tmp_path: Path) -> None:
    call = _call(HTT, _allele(18), _allele(40, (38, 41)))
    statement = _statement(_collect(tmp_path, _source(call)), "HTT")
    assert statement["kind"] == str(Kind.WATCH)
    assert "borderline" in statement["reason"]
    assert "needs confirmation" in statement["reason"]


def test_normal_is_non_finding(tmp_path: Path) -> None:
    call = _call(HTT, _allele(17, (16, 18)), _allele(20, (19, 21)))
    collected = _collect(tmp_path, _source(call))
    assert _statement(collected, "HTT")["kind"] == str(Kind.NON_FINDING)
    assert _gene_call(collected, "HTT")["phenotype"] == "normal"
    # A normal estimate whose interval reaches the intermediate range is not a non-finding.
    upper = _call(HTT, _allele(18), _allele(26, (25, 29)))
    assert _statement(_collect(tmp_path, _source(upper)), "HTT")["kind"] == "watch"


def test_a_watch_never_reads_class_normal(tmp_path: Path) -> None:
    """`class` is the highest class the longest allele's interval reaches, the bound the
    grading reads: a borderline watch whose estimate is normal reads the class its
    interval reaches, a normal the class of its estimate."""
    watch = _call(HTT, _allele(17), _allele(25, (23, 29)))
    statement = _statement(_collect(tmp_path, _source(watch)), "HTT")
    assert statement["kind"] == str(Kind.WATCH)
    assert json.loads(statement["claim"])["class"] == "intermediate"
    normal = _call(HTT, _allele(17, (16, 18)), _allele(20, (19, 21)))
    claim = json.loads(_statement(_collect(tmp_path, _source(normal)), "HTT")["claim"])
    assert claim["class"] == "normal"


def test_a_watch_beyond_the_fragment_claims_no_reliable_expansion(
    tmp_path: Path,
) -> None:
    """An allele past the fragment on too few reads is a watch: its warning says the size
    is not measured, never that the expansion is detected reliably."""
    few = _call(HTT, _allele(18), _allele(150, (120, 160), spanning=0, flanking=2))
    collected = _collect(tmp_path, _source(few))
    assert _statement(collected, "HTT")["kind"] == str(Kind.WATCH)
    warnings = _gene_call(collected, "HTT")["warnings"]
    assert BEYOND_FRAGMENT in warnings
    assert not any("reliabl" in w for w in warnings)


def test_intermediate_is_watch_with_the_note(tmp_path: Path) -> None:
    call = _call(HTT, _allele(18), _allele(30))
    statement = _statement(_collect(tmp_path, _source(call)), "HTT")
    assert statement["kind"] == str(Kind.WATCH)
    assert "intermediate" in statement["reason"]
    assert "father" in statement["reason"]


def test_too_few_reads_not_finding(tmp_path: Path) -> None:
    call = _call(HTT, _allele(18), _allele(45, (42, 48), spanning=3))
    statement = _statement(_collect(tmp_path, _source(call)), "HTT")
    assert statement["kind"] == str(Kind.WATCH)


def _level(tmp_path: Path, locus: RepeatLocus, *alleles: Allele, **kw: Any) -> Any:
    return _gene_call(
        _collect(tmp_path, _source(_call(locus, *alleles), **kw)), locus.gene
    )


def test_male_x_two_alleles_low(tmp_path: Path) -> None:
    gene_call = _level(tmp_path, AR, _allele(20), _allele(22))
    assert gene_call["level"] == str(Level.LOW)
    assert TWO_ALLELES_ON_X in gene_call["warnings"]
    assert _level(tmp_path, AR, _allele(20))["level"] == str(Level.HIGH)
    assert _level(tmp_path, AR, _allele(20), _allele(22), male=False)["level"] == "high"


def test_beyond_insert_medium(tmp_path: Path) -> None:
    # 150 CAG repeats are 450 bp, past a 400 bp fragment: the expansion is reliable, its
    # size is not (design decision, 2026-10-28).
    long = _allele(150, (120, 160))
    gene_call = _level(tmp_path, HTT, _allele(18), long)
    assert (gene_call["level"], gene_call["warnings"]) == ("medium", [BEYOND_FRAGMENT])
    # Every other cause of low still holds.
    low_depth = _call(HTT, _allele(18), long, low_depth=True)
    assert _gene_call(_collect(tmp_path, _source(low_depth)), "HTT")["level"] == "low"
    two_on_x = _level(tmp_path, AR, _allele(20), _allele(150, (120, 160)))
    assert two_on_x["level"] == "low"


def test_beyond_insert_findings_at_the_default_filters(tmp_path: Path) -> None:
    from genokeeper.query import Filters, query
    from genokeeper.sources import write_collected
    from genokeeper.store.db import connect, start_run, transaction
    from genokeeper.store.writer import StoreWriter

    c9, fxn = LOCI["C9orf72"], LOCI["FXN"]
    source = _source(
        _call(c9, _allele(5), _allele(200, (150, 300), spanning=0, inrepeat=30)),
        _call(fxn, _allele(10), _allele(300, (250, 400), spanning=0, flanking=26)),
    )
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {})
    ctx = replace(run_context(tmp_path), run_id=run_id)
    with transaction(con):
        write_collected(
            StoreWriter(con, tmp_path / "work"), con, run_id, [source.collect(ctx)]
        )
    findings = {f["gene"]: f for f in query(con, run_id, Filters()).findings}
    for gene in ("C9orf72", "FXN"):
        found = findings[gene]
        assert found["kind"] == str(Kind.FINDING)
        assert found["gene_call"]["level"] == "medium"
        assert BEYOND_FRAGMENT in found["gene_call"]["warnings"]


def test_beyond_read_medium(tmp_path: Path) -> None:
    # 60 CAG repeats are 180 bp: past a 150 bp read, inside the fragment.
    assert _level(tmp_path, HTT, _allele(18), _allele(60))["level"] == "medium"
    no_span = _allele(18, spanning=0, flanking=10)
    assert _level(tmp_path, HTT, no_span)["level"] == "medium"
    # Spanning reads must support every allele, not one of the two.
    one = _level(
        tmp_path, HTT, _allele(18), _allele(20, (19, 21), spanning=0, flanking=9)
    )
    assert (one["level"], one["warnings"]) == ("medium", [NO_SPANNING])


def test_cap_lowers_level(tmp_path: Path) -> None:
    gene_call = _level(tmp_path, HTT, _allele(18), _allele(19), caps={"HTT": Level.LOW})
    assert gene_call["level"] == str(Level.LOW)


def test_withheld_and_not_run_unassessed(tmp_path: Path) -> None:
    call = _call(HTT, _allele(18), _allele(19))
    withheld = _collect(
        tmp_path, _source(call, withheld={"HTT": "false call in calibration"})
    )
    assert _statement(withheld, "HTT")["kind"] == str(Kind.UNASSESSED)
    assert _statement(withheld, "HTT")["reason"] == "false call in calibration"
    assert _statement(withheld, "AR")["kind"] == str(Kind.UNASSESSED)  # no call
    not_run = _collect(tmp_path, _source(call, not_run="stage 3 not calibrated"))
    assert {r.statement["kind"] for r in not_run.rows} == {str(Kind.UNASSESSED)}
    assert _statement(not_run, "HTT")["reason"] == "stage 3 not calibrated"
    low_depth = _collect(tmp_path, _source(replace(call, low_depth=True)))
    assert _statement(low_depth, "HTT")["kind"] == str(Kind.UNASSESSED)
    assert _gene_call(low_depth, "HTT")["level"] == "low"


def test_collect_writes_every_locus(tmp_path: Path) -> None:
    collected = _collect(tmp_path, _source(_call(HTT, _allele(18))))
    # The 24 loci, RFC1, CFTR and the 22 not covered; a gene call for each graded one.
    assert collected.input_rows == len(collected.rows) == 48
    assert len(collected.extra["gene_call"]) == 26


def test_fxn_topic_carrier(tmp_path: Path) -> None:
    fxn = LOCI["FXN"]
    call = _call(fxn, _allele(10), _allele(70, (60, 80)))
    statement = _statement(_collect(tmp_path, _source(call)), "FXN")
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == "carrier"


def test_pcr_like_warns_c9orf72(tmp_path: Path) -> None:
    locus = LOCI["C9orf72"]
    gene_call = _level(tmp_path, locus, _allele(5), _allele(8), pcr_like=True)
    assert "library looks PCR-amplified; a miss is possible" in gene_call["warnings"]
    other = _level(tmp_path, HTT, _allele(18), _allele(19), pcr_like=True)
    assert other["warnings"] == []


def test_rows_satisfy_the_store(tmp_path: Path) -> None:
    from genokeeper.sources import write_collected
    from genokeeper.store.db import connect, start_run, transaction
    from genokeeper.store.writer import StoreWriter

    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {})
    ctx = replace(run_context(tmp_path), run_id=run_id)
    source = _source(_call(HTT, _allele(18), _allele(42, (41, 44))))
    with transaction(con):
        write_collected(
            StoreWriter(con, tmp_path / "work"), con, run_id, [source.collect(ctx)]
        )
    assert con.execute("SELECT count(*) FROM statement").fetchone() == (48,)
    assert con.execute("SELECT count(*) FROM gene_call").fetchone() == (26,)
    assert con.execute(
        "SELECT count(*) FROM statement WHERE actionable"
    ).fetchone() == (0,)


def _gene_of(call: RepeatCall) -> str:
    return next(x.gene for x in LOCI.values() if x.variant_id == call.variant_id)


def test_every_non_normal_stage3_reason_needs_confirmation(tmp_path: Path) -> None:
    """Every finding and watch item of the repeat loci, RFC1, CFTR's 5T with TG12 and SMN1
    says "screening result, needs confirmation" (the Y haplogroup's finding is ancestry, not
    a non-normal result)."""
    calls = [
        _call(HTT, _allele(18), _allele(42, (41, 44))),  # finding
        _call(HTT, _allele(18), _allele(40, (38, 41))),  # borderline
        _call(HTT, _allele(18), _allele(30)),  # intermediate
        _call(HTT, _allele(18), _allele(37)),  # reduced penetrance
        _call(HTT, _allele(18), _allele(45, (42, 48), spanning=3)),  # too few reads
        _call(LOCI["DMPK"], _allele(12), _allele(40)),  # premutation
    ]
    statements = [
        _statement(_collect(tmp_path, _source(c)), _gene_of(c)) for c in calls
    ]
    carrier = _rfc1("carrier", "1 expanded pathogenic allele, 1 short reference allele")
    statements.append(_statement(_collect(tmp_path, _source(rfc1=carrier)), "RFC1"))
    statements.append(_cftr_statement(tmp_path, (11, 12), (5, 7)))
    for caller, marker in ((1, 0), (2, 1)):
        smn = SmnSource(
            SmnCaller(caller, 2, marker, "high"),
            ParascopyCall(caller, 2, 40.0, True),
            None,
            None,
            20.0,
        )
        statements.append(smn.collect(run_context(tmp_path)).rows[0].statement)
    assert {s["kind"] for s in statements} == {str(Kind.FINDING), str(Kind.WATCH)}
    for statement in statements:
        assert "screening result, needs confirmation" in statement["reason"], statement


# Slice 3e: the thirteen loci's rules (`grade`, `excluded`).

CNBP, TBP, THAP11 = LOCI["CNBP"], LOCI["TBP"], LOCI["THAP11"]
FMR1, AFF2, FGF14 = LOCI["FMR1"], LOCI["AFF2"], LOCI["FGF14"]
TCF4, PABPN1 = LOCI["TCF4"], LOCI["PABPN1"]
NON_STANDARD = [g for g, locus in LOCI.items() if locus.rule is not Rule.STANDARD]
ELEVEN = {
    "HTT", "AR", "DMPK", "C9orf72", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A",
    "ATN1", "FXN",
}  # fmt: skip


def _graded(locus: RepeatLocus, call: RepeatCall, **kw: Any) -> Any:
    fields: dict[str, Any] = {
        "male": True,
        "x_dosage": ONE_X,
        "min_allele_reads": 5,
        "min_spanning": 3,
    }
    return grade(locus, call, **(fields | kw))


def _one(tmp_path: Path, call: RepeatCall, **kw: Any) -> dict[str, Any]:
    """The statement `collect` writes for the locus of `call`."""
    return _statement(_collect(tmp_path, _source(call, **kw)), _gene_of(call))


def test_exclusion_two_distinct_normal_is_non_finding(tmp_path: Path) -> None:
    call = _call(
        CNBP,
        _allele(18, (18, 19), spanning=6),
        _allele(22, (21, 22), spanning=5),
        longest_flanking=20,
    )
    statement = _one(tmp_path, call)
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert (
        statement["reason"] == "normal: 18 and 22 CAGG, each allele read across in full"
    )
    assert statement["exclusion"] is True
    assert _gene_call(_collect(tmp_path, _source(call)), "CNBP")["level"] == "high"


def test_exclusion_homozygous_looking_call_unassessed(tmp_path: Path) -> None:
    """Two equal alleles are never normal at an exclusion-only locus: the homozygous case
    is not measured (the closing design, amended 2026-10-30 after the final review)."""
    call = _call(CNBP, _allele(20, spanning=8), _allele(20, spanning=8))
    statement = _one(tmp_path, call)
    assert statement["kind"] == "unassessed"
    assert statement["reason"] == (
        "not assessed (the two alleles are the same length: a longer one may be hidden);"
        " needs PCR sizing"
    )
    assert statement["exclusion"] is False
    fgf14 = _call(FGF14, _allele(15, spanning=8), _allele(15, spanning=8))
    assert _one(tmp_path, fgf14)["kind"] == "unassessed"


def test_longer_flanking_guards_the_one_x_call(tmp_path: Path) -> None:
    # Two distinct spanned alleles leave no room for a hidden third: the flanking read
    # does not matter (the spec's exclusion rule, as amended 2026-10-30).
    distinct = _call(
        CNBP, _allele(18, spanning=8), _allele(22, spanning=8), longest_flanking=30
    )
    statement = _one(tmp_path, distinct)
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert statement["exclusion"] is True
    # A flanking read past the one allele of a one-X call can be size mosaicism.
    fmr1 = _call(FMR1, _allele(30, spanning=8), longest_flanking=40)
    statement = _one(tmp_path, fmr1, male=True, x_dosage=ONE_X)
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "a flanking read is longer" in statement["reason"]
    assert statement["reason"].endswith("needs a targeted PCR test")


HG002 = Path(__file__).parent / "fixtures" / "hg002"


def test_inrepeat_reads_from_the_json(tmp_path: Path) -> None:
    """HG002's own ExpansionHunter output with in-repeat reads at CNBP, FMR1 and CFTR's T
    tract: `ADIR` stays 0 for alleles shorter than a read, so the rule reads the JSON's
    `CountsOfInrepeatReads`, and all three are unassessed; without them, normal."""
    data = json.loads((HG002 / "hg002.eh3e.json").read_text())
    vcf = (HG002 / "hg002.eh3e.vcf.txt").read_text()
    # FMR1 with enough spanning reads that only the in-repeat reads can refuse it.
    vcf = vcf.replace("1:SPANNING:31:31-31:1:10:0:", "1:SPANNING:31:31-31:8:10:0:")
    (tmp_path / "eh.vcf").write_text(vcf)
    (tmp_path / "eh.json").write_text(json.dumps(data))
    before = read_calls(tmp_path / "eh.vcf", tmp_path / "eh.json")
    for locus, variant in (("CNBP", "CNBP"), ("FMR1", "FMR1"), ("CFTR", "CFTR_T")):
        data["LocusResults"][locus]["Variants"][variant]["CountsOfInrepeatReads"] = (
            "(37, 2)"
        )
    (tmp_path / "eh.json").write_text(json.dumps(data))
    after = read_calls(tmp_path / "eh.vcf", tmp_path / "eh.json")
    clean = _collect(tmp_path, _source(*before.values()))
    hidden = _collect(tmp_path, _source(*after.values()))
    for gene in ("CNBP", "FMR1"):
        assert _statement(clean, gene)["kind"] == str(Kind.NON_FINDING), gene
        statement = _statement(hidden, gene)
        assert statement["kind"] == str(Kind.UNASSESSED), gene
        assert "reads inside the repeat" in statement["reason"], gene
        assert statement["exclusion"] is False
    assert _statement(hidden, "FMR1")["reason"].endswith("needs a targeted PCR test")
    cnbp = _statement(hidden, "CNBP")
    assert cnbp["reason"].endswith("needs repeat-primed PCR")
    # The rule's unassessed keeps the alleles in the claim.
    assert [a["repeats"] for a in json.loads(cnbp["claim"])["alleles"]] == [15, 16]
    gene_call = _gene_call(hidden, "CNBP")
    assert (gene_call["level"], gene_call["diplotype"], gene_call["phenotype"]) == (
        "low",
        None,
        "not assessed",
    )
    assert _statement(clean, "CFTR")["kind"] == str(Kind.NON_FINDING)
    assert _statement(hidden, "CFTR")["reason"] == CFTR_UNREAD


def test_exclusion_too_few_spanning_unassessed(tmp_path: Path) -> None:
    call = _call(CNBP, _allele(18, spanning=2), _allele(22))
    statement = _one(tmp_path, call)
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert "an allele has fewer than 3 spanning reads" in statement["reason"]


def test_exclusion_between_ranges_unassessed(tmp_path: Path) -> None:
    tbp = _call(TBP, _allele(36), _allele(42, (41, 44), spanning=1))
    assert _one(tmp_path, tbp)["kind"] == str(Kind.UNASSESSED)
    assert _one(tmp_path, tbp)["reason"] == BETWEEN
    thap11 = _call(THAP11, _allele(20), _allele(41))
    assert _one(tmp_path, thap11)["reason"] == BETWEEN


def test_exclusion_estimate_in_disease_range_is_possible_finding(
    tmp_path: Path,
) -> None:
    call = _call(CNBP, _allele(18), _allele(80, (78, 90), spanning=0, flanking=12))
    collected = _collect(tmp_path, _source(call))
    statement = _statement(collected, "CNBP")
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["validity"] == "limited"
    assert statement["reason"].startswith(
        "possible expansion; screening result, needs confirmation"
    )
    assert statement["actionable"] is False
    assert statement["exclusion"] is False
    assert _gene_call(collected, "CNBP")["level"] == "low"
    few = _call(CNBP, _allele(18), _allele(80, (78, 90), spanning=0, flanking=1))
    assert _one(tmp_path, few)["reason"] == (
        "not assessed: an allele estimated in the disease range on too few reads; needs"
        " repeat-primed PCR"
    )


def _calls_at_twice_the_minimum(locus: RepeatLocus) -> list[RepeatCall]:
    big = locus.pathogenic_min * 2
    expanded = _allele(big, (big, big + 10), spanning=10, flanking=10, inrepeat=10)
    normal = _allele(locus.normal_max - 1)
    return [_call(locus, expanded), _call(locus, normal, expanded)]


def test_exclusion_finding_is_only_possible(tmp_path: Path) -> None:
    """The exclusion-only boundary: a non-standard locus never yields a watch, a finding
    only at validity limited worded "possible", and FMR1, AFF2 and FGF14 never one."""
    seen: set[str] = set()
    for gene in NON_STANDARD:
        locus = LOCI[gene]
        for call in _calls_at_twice_the_minimum(locus):
            for male in (True, False):
                for x_dosage in (
                    None,
                    *(
                        XDosage(n, 0.74 if n is None else n / 2, y)
                        for n in (None, 0, 1, 2, 3)
                        for y in (False, True)
                    ),
                ):
                    statement = _one(tmp_path, call, male=male, x_dosage=x_dosage)
                    seen.add(statement["kind"])
                    assert statement["kind"] != str(Kind.WATCH), gene
                    if statement["kind"] == str(Kind.FINDING):
                        assert gene not in {"FMR1", "AFF2", "FGF14"}
                        assert statement["validity"] == "limited", gene
                        assert statement["reason"].startswith("possible"), gene
    assert seen == {str(Kind.FINDING), str(Kind.UNASSESSED)}


def test_exclusion_flag_only_on_rule_normals(tmp_path: Path) -> None:
    """`exclusion` is true on the normals of the three non-standard rules alone, each with
    its confirming test, and false (with `confirm` null) on every other statement."""
    normal = {
        gene: _call(
            locus,
            _allele(max(locus.normal_max - 4, 1), spanning=6),
            *(() if locus.rule is Rule.X_LINKED else (_allele(locus.normal_max - 2),)),
        )
        for gene, locus in LOCI.items()
    }
    collected = _collect(tmp_path, _source(*normal.values()))
    for gene, locus in LOCI.items():
        statement = _statement(collected, gene)
        claim = json.loads(statement["claim"])
        flagged = locus.rule is not Rule.STANDARD
        assert statement["kind"] == str(Kind.NON_FINDING), gene
        assert statement["exclusion"] is flagged, gene
        if not flagged:
            assert claim["confirm"] is None, gene
    confirm = {
        g: json.loads(_statement(collected, g)["claim"])["confirm"]
        for g in ("FMR1", "AFF2", "CNBP", "FGF14")
    }
    assert confirm == {
        "FMR1": "a targeted PCR test",
        "AFF2": "a targeted PCR test",
        "CNBP": "PCR sizing",
        "FGF14": "PCR sizing",
    }
    # A non-normal statement is never flagged.
    other = _collect(
        tmp_path, _source(_call(CNBP, _allele(18), _allele(22, spanning=1)))
    )
    assert _statement(other, "CNBP")["exclusion"] is False
    assert json.loads(_statement(other, "CNBP")["claim"])["confirm"] is None


def test_x_linked_male_normal(tmp_path: Path) -> None:
    call = _call(FMR1, _allele(30, (30, 30), spanning=8))
    collected = _collect(tmp_path, _source(call, male=True, x_dosage=ONE_X))
    statement = _statement(collected, "FMR1")
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert statement["reason"] == (
        "normal, 30 CGG on the one X, read across in full; says nothing about anyone"
        " else's FMR1"
    )
    assert _gene_call(collected, "FMR1")["level"] == "high"


def test_x_linked_male_upper_bound_45_unassessed(tmp_path: Path) -> None:
    high = _one(tmp_path, _call(FMR1, _allele(44, (43, 45), spanning=8)))
    assert high["kind"] == str(Kind.UNASSESSED)
    assert high["reason"].endswith("needs a targeted PCR test")
    at = _one(tmp_path, _call(FMR1, _allele(44, (43, 44), spanning=8)))
    assert at["kind"] == str(Kind.NON_FINDING)
    assert at["reason"].startswith("normal, 44 CGG")


def test_x_linked_female_two_distinct_normal(tmp_path: Path) -> None:
    two = _call(FMR1, _allele(29), _allele(31))
    statement = _one(tmp_path, two, male=False, x_dosage=TWO_X)
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert statement["reason"] == (
        "normal: 29 and 31 CGG, each X's allele read across in full; says nothing about"
        " anyone else's FMR1"
    )
    same = _one(
        tmp_path, _call(FMR1, _allele(30), _allele(30)), male=False, x_dosage=TWO_X
    )
    assert same["kind"] == str(Kind.UNASSESSED)
    assert same["reason"] == (
        "not assessed (the two alleles are the same length: a longer one may be hidden);"
        " needs a targeted PCR test"
    )


def test_x_linked_dosage_mismatch_unassessed(tmp_path: Path) -> None:
    """FMR1 is graded only when the run's X dosage matches the call: a mismatch, a ratio
    between copies or two X with a Y (XXY) says the measured dosage and the PCR test; an
    uncalibrated [dosage] names `genokeeper calibrate`."""
    call = _call(FMR1, _allele(30, spanning=8))
    for male, x_dosage, copies in (
        (True, XDosage(2, 1.0, True), "2 X and a Y"),
        (False, XDosage(1, 0.5, False), "1 X"),
        (True, XDosage(None, 0.74, True), "X ratio 0.74, between copies"),
        (False, XDosage(2, 1.0, True), "2 X and a Y"),
    ):
        statement = _one(tmp_path, call, male=male, x_dosage=x_dosage)
        assert statement["kind"] == str(Kind.UNASSESSED)
        assert statement["reason"] == X_DOSAGE.format(copies=copies)
        assert statement["exclusion"] is False
    for male in (True, False):
        statement = _one(tmp_path, call, male=male, x_dosage=None)
        assert (statement["kind"], statement["reason"]) == (
            "unassessed",
            "not assessed: the X dosage is not calibrated; run `genokeeper calibrate`",
        )
        assert statement["reason"] == X_NOT_CALIBRATED
    # The run's own checks come first, as at every locus.
    for kw, reason in (
        ({"calls": {}}, "no call at this locus"),
        ({"not_run": "stage 3 not calibrated"}, "stage 3 not calibrated"),
        ({"withheld": {"FMR1": "a withheld reason"}}, "a withheld reason"),
    ):
        source = _source(call, x_dosage=None, **kw)
        assert _statement(_collect(tmp_path, source), "FMR1")["reason"] == reason


def test_aff2_normal_max_39(tmp_path: Path) -> None:
    assert AFF2.normal_max == 39
    statement = _one(tmp_path, _call(AFF2, _allele(38, (37, 40), spanning=8)))
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert _one(tmp_path, _call(AFF2, _allele(38, (37, 39), spanning=8)))[
        "kind"
    ] == str(Kind.NON_FINDING)


def test_distinct_large_estimate_unassessed_not_watch(tmp_path: Path) -> None:
    call = _call(FGF14, _allele(300, (250, 400), spanning=0, flanking=30))
    statement = _one(tmp_path, call)
    assert statement["kind"] == str(Kind.UNASSESSED)
    assert statement["reason"].endswith("needs long-range PCR or long reads")


def test_min_spanning_unset_unassessed(tmp_path: Path) -> None:
    calls = [
        _call(
            locus,
            _allele(max(locus.normal_max - 4, 1), spanning=6),
            *(() if locus.rule is Rule.X_LINKED else (_allele(locus.normal_max - 2),)),
        )
        for locus in LOCI.values()
    ]
    collected = _collect(tmp_path, _source(*calls, min_spanning=None))
    for gene in NON_STANDARD:
        statement = _statement(collected, gene)
        assert (statement["kind"], statement["reason"]) == ("unassessed", NOT_SET), gene
    for gene in ELEVEN:
        assert _statement(collected, gene)["kind"] == str(Kind.NON_FINDING), gene


def test_tcf4_finding_limited(tmp_path: Path) -> None:
    call = _call(TCF4, _allele(18), _allele(55, (51, 60), spanning=8))
    statement = _one(tmp_path, call)
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["reason"] == FINDING_REASONS["TCF4"]
    assert statement["validity"] == "limited"
    assert statement["actionable"] is False
    assert statement["exclusion"] is False


def test_tcf4_intermediate_watch(tmp_path: Path) -> None:
    statement = _one(tmp_path, _call(TCF4, _allele(18), _allele(45)))
    assert statement["kind"] == str(Kind.WATCH)
    assert "intermediate" in statement["reason"]
    assert statement["exclusion"] is False


def test_pabpn1_single_eleven_is_watch(tmp_path: Path) -> None:
    statement = _one(tmp_path, _call(PABPN1, _allele(10), _allele(11)))
    assert statement["kind"] == str(Kind.WATCH)
    assert "11 is a modifier allele" in statement["reason"]
    normal = _one(tmp_path, _call(PABPN1, _allele(10), _allele(10)))
    assert normal["kind"] == str(Kind.NON_FINDING)
    assert normal["exclusion"] is False


def test_pabpn1_eleven_eleven_possible_recessive(tmp_path: Path) -> None:
    call = _call(PABPN1, _allele(11), _allele(11))
    collected = _collect(tmp_path, _source(call))
    statement = _statement(collected, "PABPN1")
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["validity"] == "limited"
    assert statement["reason"].startswith(
        "possible recessive oculopharyngeal muscular dystrophy"
    )
    assert "screening result, needs confirmation" in statement["reason"]
    gene_call = _gene_call(collected, "PABPN1")
    assert (gene_call["level"], gene_call["phenotype"]) == (
        "low",
        "possible recessive (screening)",
    )
    dominant = _one(tmp_path, _call(PABPN1, _allele(11), _allele(12, spanning=8)))
    assert dominant["kind"] == str(Kind.FINDING)
    assert dominant["reason"] == FINDING_REASON
    assert dominant["validity"] == "strong"


# Slice 3e: RFC1's motif call, CFTR's poly-T/TG and the not-covered loci.

RFC1_CONFIRM = "repeat-primed PCR or a Southern blot"
LABELS = ("normal", "potential carrier", "carrier", "affected")


def _rfc1(call: str, description: str, *, low_depth: bool = False) -> Rfc1Call:
    alleles = (_allele(11, spanning=6), _allele(12, spanning=6))
    return Rfc1Call(call, description, alleles, 30.0, low_depth)


def _rfc1_statement(tmp_path: Path, rfc1: Rfc1Call | None, **kw: Any) -> Any:
    collected = _collect(tmp_path, _source(rfc1=rfc1, **kw))
    statement = _statement(collected, "RFC1")
    return statement, _gene_call(collected, "RFC1")


def test_rfc1_normal_seen_is_non_finding(tmp_path: Path) -> None:
    for description, level in (
        ("no expanded allele", Level.HIGH),
        ("1 expanded benign allele, 1 short reference allele", Level.MEDIUM),
    ):
        graded = rfc1_grade(_rfc1("normal", description))
        assert (graded.kind, graded.reason, graded.level) == (
            Kind.NON_FINDING,
            RFC1_NORMAL,
            level,
        )
        statement, gene_call = _rfc1_statement(tmp_path, _rfc1("normal", description))
        assert (statement["kind"], statement["reason"]) == ("non_finding", RFC1_NORMAL)
        assert gene_call["level"] == str(level)
        assert gene_call["diplotype"] == "AARRG 11/12"
    _, benign = _rfc1_statement(
        tmp_path,
        _rfc1("normal", "1 expanded benign allele, 1 short reference allele"),
    )
    assert "an expanded allele of a benign motif" in benign["warnings"]


def test_rfc1_unclear_calls_unassessed(tmp_path: Path) -> None:
    for description in (
        "expanded allele may exist but not observed",
        "2 expanded alleles may exist but not observed",
        "2 expanded alleles (possibly reference)",
    ):
        graded = rfc1_grade(_rfc1("normal", description))
        assert (graded.kind, graded.reason) == (Kind.UNASSESSED, RFC1_UNCLEAR)
        statement, _ = _rfc1_statement(tmp_path, _rfc1("normal", description))
        assert (statement["kind"], statement["reason"]) == ("unassessed", RFC1_UNCLEAR)


def test_rfc1_carrier_and_affected_are_possible_findings(tmp_path: Path) -> None:
    for call, description in RFC1_PAIRS:
        if call == "normal":
            continue
        statement, gene_call = _rfc1_statement(tmp_path, _rfc1(call, description))
        assert statement["kind"] == str(Kind.FINDING), call
        assert statement["validity"] == "limited"
        assert gene_call["level"] == "low"
        assert gene_call["phenotype"] == f"possible {call} (screening)"
        assert statement["reason"].startswith(f"possible {call}")
        assert "needs confirmation" in statement["reason"]
        assert RFC1_CONFIRM in statement["reason"]
        assert statement["actionable"] is False


def test_rfc1_finding_is_only_possible(tmp_path: Path) -> None:
    descriptions = {d for _, d in RFC1_PAIRS} | {"something new"}
    for call in (*LABELS, "unknown"):
        for description in descriptions:
            statement, _ = _rfc1_statement(tmp_path, _rfc1(call, description))
            assert statement["kind"] != str(Kind.WATCH)
            if statement["kind"] == str(Kind.FINDING):
                assert statement["validity"] == "limited"
                assert statement["reason"].startswith("possible")


def test_rfc1_normal_flagged(tmp_path: Path) -> None:
    normals = {
        "no expanded allele",
        "1 expanded benign allele, 1 short reference allele",
    }
    for call, description in [*RFC1_PAIRS, ("normal", "something new")]:
        statement, _ = _rfc1_statement(tmp_path, _rfc1(call, description))
        claim = json.loads(statement["claim"])
        flagged = description in normals
        assert statement["exclusion"] is flagged, description
        assert claim["confirm"] == (RFC1_CONFIRM if flagged else None), description
        assert claim["rfc1"] == {"call": call, "description": description}
    # Every other claim's `rfc1` is null.
    collected = _collect(tmp_path, _source(rfc1=_rfc1("normal", "no expanded allele")))
    for row in collected.rows:
        claim = json.loads(str(row.statement["claim"]))
        if row.statement["gene"] != "RFC1" and "locus_id" in claim:
            assert claim["rfc1"] is None


def test_rfc1_vocabulary() -> None:
    """The nine (label, description) pairs ehunter/locus/RFC1MotifAnalysis.cpp's
    `getRFC1Status` writes at commit 9e93322f4a, in source order, and their grades."""
    assert RFC1_PAIRS == (
        ("normal", "no expanded allele"),
        ("normal", "expanded allele may exist but not observed"),
        ("carrier", "1 expanded pathogenic allele, 1 short reference allele"),
        ("normal", "1 expanded benign allele, 1 short reference allele"),
        ("normal", "2 expanded alleles may exist but not observed"),
        ("affected", "2 expanded pathogenic alleles"),
        ("carrier", "1 expanded pathogenic allele, 1 expanded benign allele"),
        ("potential carrier", "2 expanded alleles with >30% pathogenic kmers"),
        ("normal", "2 expanded alleles (possibly reference)"),
    )
    graded = [rfc1_grade(_rfc1(call, description)) for call, description in RFC1_PAIRS]
    kinds = [g.kind for g in graded]
    assert kinds.count(Kind.NON_FINDING) == 2
    assert kinds.count(Kind.UNASSESSED) == 3
    assert kinds.count(Kind.FINDING) == 4
    assert {g.validity for g in graded if g.kind is Kind.FINDING} == {Validity.LIMITED}
    other = rfc1_grade(_rfc1("normal", "1 expanded pathogenic allele"))
    assert (other.kind, other.reason) == (
        Kind.UNASSESSED,
        "not assessed: unknown RFC1 call 'normal' (1 expanded pathogenic allele)",
    )
    # A positive label beside a description it never takes is unknown too.
    positive = rfc1_grade(_rfc1("carrier", "no expanded allele"))
    assert (positive.kind, positive.reason) == (
        Kind.UNASSESSED,
        "not assessed: unknown RFC1 call 'carrier' (no expanded allele)",
    )


def test_rfc1_not_run_withheld_unassessed(tmp_path: Path) -> None:
    normal = _rfc1("normal", "no expanded allele")
    # Unless told otherwise, the source has no motif call because the analysis did not run.
    built = {k: v for k, v in vars(_source()).items() if k != "rfc1_not_run"}
    assert RepeatSource(**built).rfc1_not_run == RFC1_NOT_RUN
    cases = [
        ({"rfc1": None, "rfc1_not_run": RFC1_NOT_RUN}, RFC1_NOT_RUN),
        ({"rfc1": normal, "rfc1_not_run": RFC1_NOT_RUN}, RFC1_NOT_RUN),
        (
            {"rfc1": normal, "not_run": "stage 3 not calibrated"},
            "stage 3 not calibrated",
        ),
        (
            {"rfc1": None, "rfc1_not_run": "rfc1 failed: no analysis"},
            "rfc1 failed: no analysis",
        ),
        (
            {"rfc1": normal, "withheld": {"RFC1": "a finding on HG003"}},
            "a finding on HG003",
        ),
        ({"rfc1": None, "rfc1_not_run": None}, "no call at this locus"),
        ({"rfc1": replace(normal, low_depth=True)}, "low depth at this locus"),
    ]
    for kw, reason in cases:
        statement = _statement(_collect(tmp_path, _source(**kw)), "RFC1")
        assert (statement["kind"], statement["reason"]) == ("unassessed", reason), kw
        assert statement["exclusion"] is False
    capped = _rfc1_statement(tmp_path, normal, caps={"RFC1": Level.MEDIUM})[1]
    assert capped["level"] == "medium"


def _cftr_calls(
    tg: tuple[int, int], t: tuple[int, int], **kw: Any
) -> tuple[RepeatCall, RepeatCall]:
    return (
        RepeatCall("CFTR_TG", tuple(_allele(n) for n in tg), 30.0, False, 150),
        RepeatCall("CFTR_T", tuple(_allele(n, **kw) for n in t), 30.0, False, 150),
    )


def _cftr_statement(
    tmp_path: Path, tg: tuple[int, int], t: tuple[int, int], **kw: Any
) -> dict[str, Any]:
    return _statement(_collect(tmp_path, _source(*_cftr_calls(tg, t, **kw))), "CFTR")


def test_cftr_grades(tmp_path: Path) -> None:
    clear = _cftr_statement(tmp_path, (11, 12), (7, 9))
    assert (clear["kind"], clear["reason"]) == (
        "non_finding",
        "no 5T allele: T 7/9, TG 11/12",
    )
    assert clear["exclusion"] is False
    five = _cftr_statement(tmp_path, (11, 12), (5, 7))
    assert (five["kind"], five["topic"], five["validity"]) == (
        "finding",
        "carrier",
        "limited",
    )
    assert five["reason"].startswith(
        "possible variable-penetrance CFTR allele: 5T with TG12"
    )
    assert "not determined" in five["reason"]
    assert "on both copies" not in five["reason"]
    collected = _collect(tmp_path, _source(*_cftr_calls((11, 12), (5, 7))))
    gene_call = _gene_call(collected, "CFTR")
    assert (gene_call["level"], gene_call["diplotype"]) == ("low", "TG 11/12 T 5/7")
    claim = json.loads(five["claim"])
    assert [a["repeats"] for a in claim["alleles"]["tg"]] == [11, 12]
    assert [a["repeats"] for a in claim["alleles"]["t"]] == [5, 7]
    assert (claim["class"], claim["thresholds"], claim["image"]) == (None, None, None)
    low = _cftr_statement(tmp_path, (11, 11), (5, 7))
    assert (low["kind"], low["validity"]) == ("finding", "limited")
    assert low["reason"].startswith("information: 5T with TG11")
    assert "alone it causes no disease" in five["reason"]
    both = _cftr_statement(tmp_path, (11, 12), (5, 5))
    assert "on both copies" in both["reason"]
    assert "alone it causes no disease" not in both["reason"]
    assert (
        "5T on both copies beside a long TG can itself cause a CFTR-related disorder"
        " (absent vas deferens in a man)"
    ) in both["reason"]
    assert both["reason"].endswith("needs confirmation by CFTR poly-T/TG testing")
    unread = _cftr_statement(tmp_path, (11, 12), (5, 7), spanning=2)
    assert (unread["kind"], unread["reason"]) == ("unassessed", CFTR_UNREAD)
    graded = [
        cftr_grade(*_cftr_calls(tg, t), min_spanning=3)
        for tg in ((11, 11), (11, 12), (12, 13))
        for t in ((5, 5), (5, 7), (7, 9))
    ]
    assert {g.kind for g in graded} == {Kind.FINDING, Kind.NON_FINDING}
    statements = [
        _cftr_statement(tmp_path, tg, t)
        for tg in ((11, 11), (11, 12))
        for t in ((5, 7), (7, 9))
    ]
    assert not any(s["actionable"] for s in statements)
    assert {s["kind"] for s in statements} == {"finding", "non_finding"}


def test_cftr_unassessed(tmp_path: Path) -> None:
    calls = _cftr_calls((11, 12), (7, 9))
    for kw, reason in (
        ({}, "no call at this locus"),
        ({"not_run": "stage 3 not calibrated"}, "stage 3 not calibrated"),
        ({"withheld": {"CFTR": "libraries disagree"}}, "libraries disagree"),
        ({"low_depth": True}, "low depth at this locus"),
    ):
        found = (
            ()
            if not kw
            else (calls[0], replace(calls[1], low_depth=True))
            if "low_depth" in kw
            else calls
        )
        kw = {k: v for k, v in kw.items() if k != "low_depth"}
        statement = _statement(_collect(tmp_path, _source(*found, **kw)), "CFTR")
        assert (statement["kind"], statement["reason"]) == ("unassessed", reason)
    capped = _gene_call(
        _collect(tmp_path, _source(*calls, caps={"CFTR": Level.MEDIUM})), "CFTR"
    )
    assert capped["level"] == "medium"


def test_not_covered_statements(tmp_path: Path) -> None:
    assert len(NOT_COVERED) == 22
    for kw in ({}, {"not_run": "stage 3 not calibrated"}):
        collected = _collect(tmp_path, _source(**kw))
        called = {g["gene"] for g in collected.extra["gene_call"]}
        for entry in NOT_COVERED:
            statement = _statement(collected, entry.gene)
            assert statement["kind"] == str(Kind.UNASSESSED)
            assert statement["reason"].startswith("not covered:"), entry.gene
            assert statement["source_level"] == "not covered"
            assert statement["exclusion"] is False
            assert set(json.loads(statement["claim"])) == {
                "schema",
                "type",
                "gene",
                "condition",
            }
            assert entry.gene not in called
    cstb = _statement(_collect(tmp_path, _source()), "CSTB")
    assert cstb["topic"] == "carrier"
    samd12 = _statement(_collect(tmp_path, _source()), "SAMD12")
    assert "(TTTCA)" in samd12["reason"]
