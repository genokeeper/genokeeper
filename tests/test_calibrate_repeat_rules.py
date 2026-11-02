"""The repeat calibration's trio checks (`genokeeper.calibrate.repeat_rules`), on hand-built
calls: HG002's two libraries, HG003 (the father) and HG004 (the mother), each normal at every
locus unless a test says otherwise."""

from collections.abc import Mapping
from dataclasses import replace
from typing import Any

import pytest

from genokeeper.calibrate.repeat_rules import (
    CHILD,
    FATHER,
    LIBRARIES_DISAGREE,
    MOTHER,
    NOT_INHERITED,
    NOT_VALIDATED,
    PARENT_FINDING,
    RFC1_WITHHELD,
    RULE_DID_NOT_HOLD,
    X_RULE_DID_NOT_HOLD,
    SampleRepeats,
    measure_repeat_rules,
    mendelian,
    merge_withheld,
)
from genokeeper.calibrate.stage3 import (
    FALSE_CALL,
    NO_CALL,
    RULES,
    SMN_KEY,
    SMN_WITHHELD,
    Stage3Table,
    report,
)
from genokeeper.config import load_grading
from genokeeper.repeats.loci import Allele, RepeatCall, Rfc1Call
from genokeeper.repeats.rfc1 import RFC1_LOCUS
from genokeeper.sources.repeats import CFTR_GENE, CFTR_T, CFTR_TG
from genokeeper.store.db import Level

LOCI = load_grading().repeats.loci
CNBP, FMR1, FGF14 = LOCI["CNBP"], LOCI["FMR1"], LOCI["FGF14"]
HG002 = (CHILD, "pcr")
READ = 150


def _call(variant: str, *counts: int, spanning: int = 20) -> RepeatCall:
    alleles = tuple(Allele(n, (n, n), spanning, 0, 0) for n in counts)
    return RepeatCall(variant, alleles, 30.0, False, READ)


def _rfc1(
    call: str = "normal",
    description: str = "1 expanded benign allele, 1 short reference allele",
) -> Rfc1Call:
    return Rfc1Call(call, description, (), 30.0, False)


def _sample(male: bool) -> SampleRepeats:
    """Normal everywhere: 5/7 on an autosome, 7 on a male X, 6/8 on a female's; CFTR TG
    11/12 and T 7/9; RFC1 with HG002's benign expanded allele seen."""
    calls = {
        locus.variant_id: (
            _call(locus.variant_id, 7)
            if male and locus.chrom == "chrX"
            else _call(locus.variant_id, 6, 8)
            if locus.chrom == "chrX"
            else _call(locus.variant_id, 5, 7)
        )
        for locus in LOCI.values()
    }
    calls[CFTR_TG] = _call(CFTR_TG, 11, 12)
    calls[CFTR_T] = _call(CFTR_T, 7, 9)
    return SampleRepeats(sex="male" if male else "female", calls=calls, rfc1=_rfc1())


def _trio() -> dict[str, SampleRepeats]:
    return {
        CHILD: _sample(True),
        "pcr": _sample(True),
        FATHER: _sample(True),
        MOTHER: _sample(False),
    }


def _with(
    samples: Mapping[str, SampleRepeats], keys: tuple[str, ...], *calls: RepeatCall
) -> dict[str, SampleRepeats]:
    """`samples` with `calls` replacing the calls of each of `keys`."""
    changed = {c.variant_id: c for c in calls}
    return {
        key: replace(s, calls={**s.calls, **changed}) if key in keys else s
        for key, s in samples.items()
    }


def _without(
    samples: Mapping[str, SampleRepeats], key: str, variant: str
) -> dict[str, SampleRepeats]:
    s = samples[key]
    calls = {v: c for v, c in s.calls.items() if v != variant}
    return {**samples, key: replace(s, calls=calls)}


def _measure(
    samples: Mapping[str, SampleRepeats],
) -> tuple[dict[str, str], dict[str, Level], dict[str, Any]]:
    return measure_repeat_rules(samples, LOCI, min_spanning=3, min_allele_reads=4)


def test_a_normal_trio_withholds_nothing() -> None:
    withheld, caps, measured = _measure(_trio())
    assert (withheld, caps) == ({}, {})
    assert set(measured) == {CHILD, "pcr", FATHER, MOTHER}


def test_measure_needs_the_four_samples() -> None:
    samples = _trio()
    del samples[MOTHER]
    with pytest.raises(ValueError, match="hg004"):
        _measure(samples)


def test_merge_withheld_precedence() -> None:
    assert merge_withheld(
        {"CNBP": RULE_DID_NOT_HOLD, "HTT": PARENT_FINDING},
        {"CNBP": NO_CALL, "HTT": FALSE_CALL, "TBP": LIBRARIES_DISAGREE},
    ) == {"CNBP": NO_CALL, "HTT": FALSE_CALL, "TBP": LIBRARIES_DISAGREE}
    # A key's only reason passes through: RFC1's, and 3a's SMN and Y haplogroup.
    assert merge_withheld({RFC1_LOCUS: RFC1_WITHHELD}, {SMN_KEY: SMN_WITHHELD}) == {
        RFC1_LOCUS: RFC1_WITHHELD,
        SMN_KEY: SMN_WITHHELD,
    }
    # A size off the pedigree before a parent's finding; no validating normal last.
    assert merge_withheld(
        {"CNBP": PARENT_FINDING, "TBP": NOT_VALIDATED},
        {"CNBP": NOT_INHERITED, "TBP": PARENT_FINDING},
    ) == {"CNBP": NOT_INHERITED, "TBP": PARENT_FINDING}
    with pytest.raises(ValueError, match="not covered: something else"):
        merge_withheld({"CNBP": "not covered: something else"})
    with pytest.raises(ValueError, match=RFC1_LOCUS):
        merge_withheld({RFC1_LOCUS: RFC1_WITHHELD}, {RFC1_LOCUS: NO_CALL})


def test_rfc1_needs_hg002s_benign_expansion_seen() -> None:
    """RFC1 is kept only when both HG002 libraries see its benign ~115-unit AAAAG allele
    (the closing design, slice 3e) and HG003 and HG004 are non-findings; a sample with no
    RFC1 call withholds RFC1 alone."""
    samples = _trio()
    carrier = _rfc1("carrier", "1 expanded pathogenic allele, 1 short reference allele")
    missed = _rfc1(description="no expanded allele")
    for key, rfc1 in (
        (MOTHER, carrier),
        (FATHER, _rfc1(description="2 expanded alleles (possibly reference)")),
        ("pcr", None),
        (FATHER, None),
        (CHILD, missed),
        ("pcr", missed),
    ):
        changed = {**samples, key: replace(samples[key], rfc1=rfc1)}
        withheld, _, measured = _measure(changed)
        assert withheld == {RFC1_LOCUS: RFC1_WITHHELD}, (key, rfc1)
    assert measured["pcr"]["rfc1"]["kind"] == "non_finding"
    # The parents need a normal, not HG002's benign expansion.
    plain = {k: replace(s, rfc1=missed) for k, s in samples.items() if k not in HG002}
    assert _measure({**samples, **plain})[0] == {}


def test_hg002s_libraries_must_agree() -> None:
    """HG002's two libraries, one DNA sequenced twice, must agree: a normal beside an
    unassessed call, or sizes more than `MENDEL_TOLERANCE` apart, withholds the locus."""
    unread = _call(CNBP.variant_id, 5, 7, spanning=1)
    withheld, _, _ = _measure(_with(_trio(), ("pcr",), unread))
    assert withheld == {"CNBP": LIBRARIES_DISAGREE}
    withheld, _, _ = _measure(_with(_trio(), ("pcr",), _call(CNBP.variant_id, 5, 10)))
    assert withheld == {"CNBP": LIBRARIES_DISAGREE}
    # Within the tolerance the libraries agree.
    assert _measure(_with(_trio(), ("pcr",), _call(CNBP.variant_id, 5, 9)))[0] == {}
    # A parent the rule does not hold on is no finding: nothing withheld.
    assert _measure(_with(_trio(), (FATHER, MOTHER), unread))[0] == {}
    withheld, _, _ = _measure(_without(_trio(), "pcr", CNBP.variant_id))
    assert withheld == {"CNBP": NO_CALL}
    withheld, _, _ = _measure(_with(_trio(), (CHILD,), _call(FGF14.variant_id, 7, 7)))
    assert withheld == {"FGF14": LIBRARIES_DISAGREE}


def test_a_finding_on_a_healthy_sample_withholds() -> None:
    """A finding on HG002 at a non-standard locus withholds it `RULE_DID_NOT_HOLD`, on a
    parent `PARENT_FINDING`."""
    expanded = _call(CNBP.variant_id, 5, 100)
    withheld, _, _ = _measure(_with(_trio(), HG002, expanded))
    assert withheld == {"CNBP": RULE_DID_NOT_HOLD}
    withheld, _, _ = _measure(_with(_trio(), (FATHER,), expanded))
    assert withheld == {"CNBP": PARENT_FINDING}


def test_a_locus_no_sample_gave_a_normal_is_withheld() -> None:
    """All four unassessed: the rule was never seen to give a normal."""
    unread = _call(CNBP.variant_id, 5, 7, spanning=1)
    withheld, _, _ = _measure(_with(_trio(), (*HG002, FATHER, MOTHER), unread))
    assert withheld == {"CNBP": NOT_VALIDATED}
    # One parent's normal validates it.
    assert _measure(_with(_trio(), (*HG002, FATHER), unread))[0] == {}


def test_x_linked_needs_hg003_and_the_two_x_branch_needs_hg004() -> None:
    fmr1 = FMR1.variant_id
    samples = _with(_trio(), (*HG002, FATHER), _call(fmr1, 30))
    withheld, _, _ = _measure(_with(samples, (MOTHER,), _call(fmr1, 30, 30)))
    assert withheld == {"FMR1:two_x": X_RULE_DID_NOT_HOLD}
    samples = _with(_trio(), (*HG002, FATHER), _call(fmr1, 31))
    assert _measure(_with(samples, (MOTHER,), _call(fmr1, 29, 31)))[0] == {}
    samples = _with(_trio(), (*HG002, FATHER), _call(fmr1, 40))
    withheld, _, _ = _measure(_with(samples, (MOTHER,), _call(fmr1, 29, 31)))
    assert withheld == {"FMR1:two_x": X_RULE_DID_NOT_HOLD}
    # HG003 unassessed beside HG002's one-X normals: kept.
    assert _measure(_with(_trio(), (FATHER,), _call(fmr1, 7, spanning=1)))[0] == {}
    hidden = replace(_call(fmr1, 7), inrepeat=2)
    withheld, _, _ = _measure(_with(_trio(), ("pcr",), hidden))
    assert withheld == {"FMR1": LIBRARIES_DISAGREE}


def test_x_linked_needs_one_one_x_normal() -> None:
    """HG002's libraries unassessed (one spanning read, as at FMR1) beside HG003's one-X
    normal keep FMR1; all three unassessed withhold it, whatever HG004 shows."""
    fmr1 = FMR1.variant_id
    samples = _with(_trio(), HG002, _call(fmr1, 31, spanning=1))
    samples = _with(samples, (FATHER,), _call(fmr1, 31))
    samples = _with(samples, (MOTHER,), _call(fmr1, 29, 31))
    withheld, _, measured = _measure(samples)
    assert withheld == {}
    assert measured[CHILD]["repeats"]["FMR1"]["kind"] == "unassessed"
    withheld, _, _ = _measure(_with(samples, (FATHER,), _call(fmr1, 31, spanning=1)))
    assert withheld == {"FMR1": NOT_VALIDATED}


def test_parent_finding_withholds() -> None:
    withheld, _, _ = _measure(
        _with(_trio(), (FATHER,), _call(LOCI["TCF4"].variant_id, 5, 60))
    )
    assert withheld == {"TCF4": PARENT_FINDING}
    five_t = _with(_trio(), (MOTHER,), _call(CFTR_TG, 11, 12), _call(CFTR_T, 5, 7))
    assert _measure(five_t)[0] == {}


def test_trio_checks_apply_to_the_eleven() -> None:
    htt = LOCI["HTT"].variant_id
    withheld, _, _ = _measure(_with(_trio(), (MOTHER,), _call(htt, 17, 45)))
    assert withheld == {"HTT": PARENT_FINDING}
    withheld, caps, _ = _measure(_with(_trio(), (CHILD,), _call(htt, 5, 10)))
    assert withheld == {} and caps == {"HTT": Level.LOW}


def test_cftr_libraries_disagree_withheld() -> None:
    withheld, _, _ = _measure(_with(_trio(), ("pcr",), _call(CFTR_T, 7, 7)))
    assert withheld == {CFTR_GENE: LIBRARIES_DISAGREE}
    # Other TG counts disagree too.
    withheld, _, _ = _measure(_with(_trio(), ("pcr",), _call(CFTR_TG, 11, 13)))
    assert withheld == {CFTR_GENE: LIBRARIES_DISAGREE}
    withheld, _, _ = _measure(_without(_trio(), "pcr", CFTR_TG))
    assert withheld == {CFTR_GENE: NO_CALL}
    # Both libraries unread compare equal, but the rule did not hold on HG002.
    unread = _call(CFTR_T, 7, 9, spanning=1)
    withheld, _, _ = _measure(_with(_trio(), HG002, unread))
    assert withheld == {CFTR_GENE: RULE_DID_NOT_HOLD}


def test_mendelian() -> None:
    def m(chrom: str, child: tuple[int, ...], father: tuple[int, ...] | None,
          mother: tuple[int, ...] | None) -> bool | None:  # fmt: skip
        return mendelian(
            "CAG",
            chrom,
            _call("X", *child),
            None if father is None else _call("X", *father),
            None if mother is None else _call("X", *mother),
        )

    assert m("chr4", (5, 9), (5, 6), (7, 12)) is True  # 9 is 2 from the mother's 7
    assert m("chr4", (9, 5), (5, 6), (7, 12)) is True
    assert m("chr4", (5, 10), (5, 6), (5, 7)) is False  # 10 is 3 from every allele
    assert m("chr4", (5, 7), (5, 7), None) is None
    assert m("chrX", (7,), (40,), (6, 20)) is True
    assert m("chrX", (10,), (10,), (6, 7)) is False
    assert m("chr4", (5, 60), (5, 60), (5, 60)) is None  # 180 bases, beyond the read
    assert m("chr4", (5, 50), (5, 50), (5, 50)) is None  # 150 bases, not under it


def test_mendelian_failure_caps_low() -> None:
    withheld, caps, measured = _measure(
        _with(_trio(), HG002, _call(CNBP.variant_id, 5, 10))
    )
    assert (withheld, caps) == ({}, {"CNBP": Level.LOW})
    assert measured[CHILD]["repeats"]["CNBP"]["mendelian"] is False
    _, caps, _ = _measure(_with(_trio(), (CHILD,), _call(CFTR_TG, 11, 15)))
    assert caps == {CFTR_GENE: Level.LOW}


def test_hg002_homozygous_rests_on_the_parents() -> None:
    """HG002 homozygous-looking at an EXCLUSION locus is graded as the run grades it,
    unassessed; the locus is kept on the parents' normals when HG002's size is within
    `MENDEL_TOLERANCE` of an allele of each parent, and withheld otherwise or when a parent
    has no call."""
    cnbp = CNBP.variant_id
    samples = _with(_trio(), HG002, _call(cnbp, 10, 10))
    samples = _with(samples, (FATHER,), _call(cnbp, 10, 12))
    samples = _with(samples, (MOTHER,), _call(cnbp, 10, 11))
    withheld, caps, measured = _measure(samples)
    assert (withheld, caps) == ({}, {})
    assert measured[CHILD]["repeats"]["CNBP"]["kind"] == "unassessed"
    # HG002's 10 is more than 2 from each of HG004's alleles.
    withheld, _, _ = _measure(_with(samples, (MOTHER,), _call(cnbp, 13, 15)))
    assert withheld == {"CNBP": NOT_INHERITED}
    withheld, _, _ = _measure(_without(samples, MOTHER, cnbp))
    assert withheld == {"CNBP": NOT_INHERITED}
    # Against the default parents (5/7 each) HG002's 20/20 follows neither.
    withheld, _, _ = _measure(_with(_trio(), HG002, _call(cnbp, 20, 20)))
    assert withheld == {"CNBP": NOT_INHERITED}


def test_report_lists_the_rules() -> None:
    withheld, caps, measured = _measure(
        _with(_trio(), HG002, _call(CNBP.variant_id, 5, 10))
    )
    table = Stage3Table(
        versions={},
        withheld=withheld | {"FMR1:two_x": X_RULE_DID_NOT_HOLD},
        caps=caps,
        flank={},
        measured={RULES: measured},
    )
    text = report(table)
    assert "trial" not in text and "homozygous_ok" not in text
    (cnbp,) = [line for line in text.splitlines() if line.startswith("  CNBP\t")]
    assert "pcrfree 5/10 non_finding, 20 spanning, Mendelian no" in cnbp
    assert cnbp.endswith("capped at low")
    assert f"FMR1:two_x\twithheld: {X_RULE_DID_NOT_HOLD}" in text
    assert (
        "RFC1\tpcrfree normal (1 expanded benign allele, 1 short reference allele)"
        in text
    )
    assert "CFTR\tpcrfree TG 11/12 T 7/9 non_finding" in text
