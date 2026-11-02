"""HLA alleles beyond PharmCAT's, and the coeliac DQ types from T1K's class II alleles and
from Monsuur et al. 2008's six tag SNPs (slice 3f, "HLA")."""

import re

import pytest

from genokeeper.pgx.hla import PGX_HLA, HlaCall
from genokeeper.pgx.hla_typing import (
    B27_NOTE,
    EXCLUDED,
    TAG_ALLELES,
    TAG_RSIDS,
    TAGS_EUROPEAN,
    DqStatus,
    carried_reason,
    coeliac,
    dq_from_t1k,
    dq_from_tags,
    load_hg002_truth,
    load_hla_alleles,
)
from genokeeper.pgx.pharmcat import DOSING
from genokeeper.store.db import Kind, Topic, Validity

NONE = {rsid: "0" for rsid in TAG_RSIDS}


def tags(**dosages: str) -> dict[str, str]:
    """The six tag SNPs' dosages, reference except those named."""
    return NONE | dosages


def t1k(dqa1: HlaCall, dqb1: HlaCall) -> dict[str, HlaCall]:
    return {"HLA-DQA1": dqa1, "HLA-DQB1": dqb1}


def test_tag_alleles_are_monsuurs_six_on_grch38() -> None:
    """The six tags at their GRCh38 positions, the tag allele the ALT in each (pre-flight
    record, dbSNP 157)."""
    assert dict(zip(TAG_RSIDS, TAG_ALLELES, strict=True)) == {
        "rs2187668": ("chr6", 32638107, "C", "T"),
        "rs2395182": ("chr6", 32445540, "G", "T"),
        "rs7775228": ("chr6", 32690302, "T", "C"),
        "rs4713586": ("chr6", 32691805, "A", "G"),
        "rs7454108": ("chr6", 32713706, "T", "C"),
        "rs4639334": ("chr6", 32634437, "G", "A"),
    }


def test_dq_from_tags_monsuurs_rules() -> None:
    """DQ2.5, DQ8 and DQ7 copies are their tag's dosage (DQ7 the half-heterodimer); DQ2.2 is
    the T,G,A haplotype, T,G,G being DQ4; an unknown dosage is no status."""
    assert dq_from_tags(NONE) == DqStatus(dq25=0, dq22=0, dq8=0, half=False)
    assert dq_from_tags(tags(rs2187668="2")) == DqStatus(2, 0, 0, False)
    assert dq_from_tags(tags(rs7454108="1")) == DqStatus(0, 0, 1, False)
    assert dq_from_tags(tags(rs4639334="1")) == DqStatus(0, 0, 0, True)
    # HG002 (v5.0q): rs2395182 1|1, rs7454108 1|0, the rest reference.
    assert dq_from_tags(tags(rs2395182="2", rs7454108="1")) == DqStatus(0, 0, 1, False)
    # DQ2.2: T at rs2395182 on both chromosomes, one C at rs7775228, A at rs4713586.
    assert dq_from_tags(tags(rs2395182="2", rs7775228="1")) == DqStatus(0, 1, 0, False)
    assert dq_from_tags(tags(rs2395182="2", rs7775228="2")) == DqStatus(0, 2, 0, False)
    # T,G,G on both chromosomes is DQ4, not DQ2.2.
    assert dq_from_tags(tags(rs2395182="2", rs7775228="2", rs4713586="2")) == DqStatus(
        0, 0, 0, False
    )
    assert dq_from_tags(tags(rs2187668="unknown")) is None
    assert dq_from_tags({k: v for k, v in NONE.items() if k != "rs4639334"}) is None


@pytest.mark.parametrize(
    ("a", "b", "c"),
    [
        # One T and one C: on the same chromosome or not.
        ("1", "1", "0"),
        # Two Cs and one G: the G on the T chromosome or the other.
        ("1", "2", "1"),
        # One C and one G: together (DQ4) or apart (DQ2.2 beside a G).
        ("2", "1", "1"),
    ],
)
def test_dq22_unphased_ambiguity_is_no_status(a: str, b: str, c: str) -> None:
    """Unphased, DQ2.2's copies are known only when every phasing gives the same count."""
    assert dq_from_tags(tags(rs2395182=a, rs7775228=b, rs4713586=c)) is None


def test_dq_from_t1k() -> None:
    """T1K's two DQA1 and two DQB1 alleles give the full types and the half-heterodimer
    (DQA1*05 or DQB1*02); an allele dropped is no status."""
    assert dq_from_t1k(t1k(("*05:01", "*01:01"), ("*02:01", "*05:01"))) == DqStatus(
        1, 0, 0, True
    )
    assert dq_from_t1k(t1k(("*02:01", "*01:01"), ("*02:02", "*05:01"))) == DqStatus(
        0, 1, 0, True
    )
    # HG002: DQA1*03:01/*01:05, DQB1*03:02/*05:01, DQ8 heterozygous.
    assert dq_from_t1k(t1k(("*03:01", "*01:05"), ("*03:02", "*05:01"))) == DqStatus(
        0, 0, 1, False
    )
    assert dq_from_t1k(t1k(("*05:05", "*01:01"), ("*03:01", "*05:01"))) == DqStatus(
        0, 0, 0, True
    )
    assert dq_from_t1k(t1k(("*01:01", "*01:02"), ("*05:01", "*06:02"))) == DqStatus(
        0, 0, 0, False
    )
    assert dq_from_t1k(t1k(("*01:01", "*01:02"), (None, "*05:01"))) is None
    assert dq_from_t1k({"HLA-DQA1": ("*01:01", "*01:02")}) is None


def test_equivalents_spanning_risk_unknown() -> None:
    """T1K's equivalent names of one allele that span a risk and a non-risk name: no
    status, so the coeliac statement is unassessed whatever the tags say."""
    spans = [
        t1k(("*03:01", "*01:01"), ("*03:02,*03:01", "*05:01")),
        t1k(("*05:05,*06:01", "*01:01"), ("*05:01", "*06:02")),
    ]
    for calls in spans:
        assert dq_from_t1k(calls) is None
        for tag in (NONE, tags(rs7454108="1"), tags(rs4639334="1")):
            kind, _ = coeliac(dq_from_t1k(calls), dq_from_tags(tag))
            assert kind is Kind.UNASSESSED
    # Equivalents that agree on the risk (both DQB1*02:01 at two fields) are one allele.
    assert dq_from_t1k(t1k(("*05:01", "*01:01"), ("*02:01", "*05:01,*05:02"))) == (
        DqStatus(1, 0, 0, True)
    )


def test_coeliac_rule_out_needs_both_callers_clean() -> None:
    nothing = DqStatus(0, 0, 0, False)
    assert coeliac(nothing, nothing) == (
        Kind.NON_FINDING,
        f"{EXCLUDED}; {TAGS_EUROPEAN}",
    )
    half = DqStatus(0, 0, 0, True)
    for t, g in ((half, nothing), (nothing, half), (None, nothing), (nothing, None)):
        kind, reason = coeliac(t, g)
        assert kind is Kind.UNASSESSED
        assert EXCLUDED not in reason


def test_dq_type_disagreement_unassessed() -> None:
    """Both see a risk type but not the same full types: unassessed, naming neither."""
    kind, reason = coeliac(DqStatus(1, 0, 0, True), DqStatus(0, 0, 1, False))
    assert kind is Kind.UNASSESSED
    assert "DQ2" not in reason and "DQ8" not in reason


def test_carried_reason_wording() -> None:
    common = ": common (about 30–40% of Europeans), not a diagnosis"
    assert (
        carried_reason(DqStatus(1, 0, 0, True)) == f"possibly carries HLA-DQ2{common}"
    )
    assert (
        carried_reason(DqStatus(0, 1, 0, True)) == f"possibly carries HLA-DQ2{common}"
    )
    assert (
        carried_reason(DqStatus(0, 0, 2, False)) == f"possibly carries HLA-DQ8{common}"
    )
    assert carried_reason(DqStatus(1, 0, 1, True)) == (
        f"possibly carries HLA-DQ2 and -DQ8{common}"
    )
    half = carried_reason(DqStatus(0, 0, 0, True))
    assert half == (
        "possibly carries half of HLA-DQ2 (DQA1*05 or DQB1*02 alone): coeliac disease is"
        " rare with it alone; not a diagnosis"
    )
    with pytest.raises(ValueError, match="carries"):
        carried_reason(DqStatus(0, 0, 0, False))


def test_alleles_table() -> None:
    """The class I alleles PharmCAT does not read, one allele per row (ruling R4), and the
    research associations; none of PharmCAT's own."""
    rows = {(a.gene, a.allele, a.topic): a for a in load_hla_alleles()}
    pharmcat = {(g, a) for g, alleles in PGX_HLA.items() for a in alleles}
    assert not pharmcat & {(g, a) for g, a, _ in rows}
    b13 = rows["HLA-B", "*13:01", Topic.DRUGS]
    assert (b13.drugs, b13.validity, b13.level) == (
        ("dapsone", "co-trimoxazole"),
        Validity.LIMITED,
        "PharmGKB 2A",
    )
    assert "1183681792" in b13.source and "1451704000" in b13.source
    shared = [a for a in rows.values() if "1451703980" in a.source]
    assert sorted(a.allele for a in shared) == ["*06:02", "*07:27", "*08:01"]
    b27 = rows["HLA-B", "*27", Topic.RISKS]
    assert (b27.validity, b27.note, b27.drugs) == (Validity.MODERATE, B27_NOTE, ())
    research = {(g, a) for g, a, t in rows if t is Topic.RESEARCH}
    assert research == {("HLA-C", "*06:02"), ("HLA-B", "*51"), ("HLA-DQB1", "*06:02")}
    for a in rows.values():
        assert not DOSING.search(" ".join((a.condition, a.note, a.level, *a.drugs)))
        assert re.fullmatch(r"\*\d+(:\d+)?", a.allele)


def test_hg002_truth() -> None:
    """HG002's published types (Chin et al. 2020), maternal first; the maternal DQA1 written
    so that the G group's 01:01 and the sequence's 01:05 both match."""
    truth = load_hg002_truth()
    assert truth["HLA-B"] == ("*35:08", "*38:01")
    assert truth["HLA-DQA1"] == ("*01:05,*01:01", "*03:01")
    assert set(truth) == {
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-DQA1",
        "HLA-DQB1",
        "HLA-DRB1",
    }
