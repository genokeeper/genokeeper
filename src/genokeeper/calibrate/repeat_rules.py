"""Slice 3e's repeat calibration on the trio: HG002's two libraries (`pcrfree`, `pcr`), its
father HG003 and its mother HG004, all healthy, each graded as a run would grade it (spec,
slice 3e "Calibration"; the plan's "Values and builder decisions").

`measure_repeat_rules` gives what the trio withholds (key -> one reason, `merge_withheld`)
and caps:

- a non-standard locus, each sample graded as the run grades it, is withheld only on
  evidence that the rule misleads or for lack of any validation (amended 2026-10-31): a
  finding on HG002; HG002's two libraries disagreeing (a normal beside an unassessed call,
  or sizes more than `MENDEL_TOLERANCE` repeats apart); a homozygous-looking HG002, which
  is unassessed, whose size is not within `MENDEL_TOLERANCE` of an allele of each parent
  (a parent with no call withholds too); or no sample graded normal there (for an X-linked
  locus, none of HG002's libraries and HG003 a one-X normal). An unassessed sample alone
  withholds nothing: an unassessed result never misleads;
- an X-linked locus's two-X branch (`two_x_key`) is withheld unless HG004 is a two-X normal
  whose alleles include HG002's within `MENDEL_TOLERANCE` repeats;
- RFC1 is withheld unless both HG002 libraries see its benign expanded allele
  (`HG002_RFC1`) and HG003 and HG004 are normals; CFTR unless its rule holds on both HG002
  libraries and they agree;
- a finding on HG003 or HG004 withholds its locus (CFTR excepted: its alleles are common in
  healthy people);
- an HG002 allele shorter than the read with no parental allele within `MENDEL_TOLERANCE`
  repeats (`mendelian`) caps the locus at low.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from genokeeper.calibrate.records import RECORDED
from genokeeper.calibrate.stage3 import (
    FALSE_CALL,
    NO_CALL,
    SMN_WITHHELD,
    YLEAF_WITHHELD,
    two_x_key,
)
from genokeeper.repeats.loci import RepeatCall, RepeatLocus, Rfc1Call, Rule
from genokeeper.repeats.rfc1 import RFC1_LOCUS
from genokeeper.sex import Sex
from genokeeper.sources.repeats import (
    CFTR_GENE,
    CFTR_T,
    CFTR_TG,
    Graded,
    XDosage,
    cftr_grade,
    grade,
    rfc1_grade,
)
from genokeeper.store.db import Kind, Level

CHILD = "pcrfree"
FATHER = "hg003"
MOTHER = "hg004"
_HG002 = (CHILD, "pcr")
# 3a's size bound (`stage3.MAX_SIZE_ERROR`): an HG002 allele further than this from every
# parental allele is not inherited as called.
MENDEL_TOLERANCE = 2
# CFTR's two tracts (`sources.repeats`), each with its motif, on chr7.
_CFTR = ((CFTR_TG, "TG"), (CFTR_T, "T"))
_CFTR_CHROM = "chr7"

RFC1_WITHHELD = "not covered: the RFC1 motif call was not a clear normal on a healthy calibration sample"
# HG002's RFC1 call on both libraries: its benign ~115-unit AAAAG allele seen (the closing
# design, slice 3e; GIAB's v5.0q truth, R5). Its anchor comes from this constant, which an
# ExpansionHunter build that misses the allele cannot move because the pair is code, not
# read from the build's output.
HG002_RFC1 = ("normal", "1 expanded benign allele, 1 short reference allele")
LIBRARIES_DISAGREE = "not covered: HG002's two libraries disagree in calibration"
RULE_DID_NOT_HOLD = (
    "not covered: the exclusion rule did not hold on HG002 in calibration"
)
X_RULE_DID_NOT_HOLD = "not covered: the X rule did not hold on the trio in calibration"
PARENT_FINDING = "not covered: a finding on a healthy calibration parent"
NOT_INHERITED = "not covered: HG002's sizes do not follow its parents' in calibration"
NOT_VALIDATED = (
    "not covered: the rule never gave a normal on a healthy calibration sample"
)
# One reason per withheld key, the most basic failure first.
WITHHOLD_ORDER = (
    NO_CALL,
    FALSE_CALL,
    LIBRARIES_DISAGREE,
    RULE_DID_NOT_HOLD,
    X_RULE_DID_NOT_HOLD,
    NOT_INHERITED,
    PARENT_FINDING,
    NOT_VALIDATED,
)
# Reasons each its key's only one: RFC1's, and 3a's SMN and Y haplogroup.
_ONLY = frozenset({RFC1_WITHHELD, SMN_WITHHELD, YLEAF_WITHHELD})


def merge_withheld(*maps: Mapping[str, str]) -> dict[str, str]:
    """`maps` merged: per key the reason earliest in `WITHHOLD_ORDER`, or its key's only
    reason (`RFC1_WITHHELD`, 3a's `SMN_WITHHELD` and `YLEAF_WITHHELD`). `ValueError` for any
    other reason, or an only reason beside another."""
    # A reason outside the order is refused, never ranked by a guess. Its anchor comes from
    # `WITHHOLD_ORDER` and `_ONLY`, this module's constants, which a caller's map cannot
    # move because the ranks are code, not read from the maps.
    rank = {reason: i for i, reason in enumerate(WITHHOLD_ORDER)}
    given: dict[str, set[str]] = {}
    for found in maps:
        for key, reason in found.items():
            given.setdefault(key, set()).add(reason)
    merged: dict[str, str] = {}
    for key, reasons in given.items():
        if reasons & _ONLY:
            if len(reasons) > 1:
                raise ValueError(
                    f"{key}: {sorted(reasons)}: an only reason with another"
                )
            (merged[key],) = reasons
            continue
        unknown = reasons - set(rank)
        if unknown:
            raise ValueError(f"{key}: no withholding order for {sorted(unknown)}")
        merged[key] = min(reasons, key=rank.__getitem__)
    return merged


@dataclass(frozen=True)
class SampleRepeats:
    """One trio sample's repeat calls (by `VARID`) and its RFC1 motif call."""

    sex: Sex
    calls: Mapping[str, RepeatCall]
    rfc1: Rfc1Call | None


def _near(repeats: int, parent: RepeatCall) -> bool:
    return any(abs(repeats - a.repeats) <= MENDEL_TOLERANCE for a in parent.alleles)


def mendelian(
    motif: str,
    chrom: str,
    child: RepeatCall | None,
    father: RepeatCall | None,
    mother: RepeatCall | None,
) -> bool | None:
    """Whether `child`'s alleles are inherited within `MENDEL_TOLERANCE` repeats: on chrX
    the son's one allele from a mother's allele, elsewhere one allele from a father's and the
    other from a mother's. None when a call it needs is missing or has another allele count
    (the father's is not needed on chrX), or a child allele is not shorter than the read."""
    if child is None or any(
        a.repeats * len(motif) >= child.read_length for a in child.alleles
    ):
        return None
    if chrom == "chrX":
        if mother is None or not mother.alleles or len(child.alleles) != 1:
            return None
        return _near(child.alleles[0].repeats, mother)
    if (
        father is None
        or mother is None
        or not father.alleles
        or not mother.alleles
        or len(child.alleles) != 2
    ):
        return None
    a, b = (allele.repeats for allele in child.alleles)
    return (_near(a, father) and _near(b, mother)) or (
        _near(b, father) and _near(a, mother)
    )


def _alleles(call: RepeatCall | None) -> list[int] | None:
    return (
        None if call is None or not call.alleles else [a.repeats for a in call.alleles]
    )


def _kind(graded: Graded | None) -> str | None:
    return None if graded is None else str(graded.kind)


def _cftr_grade(sample: SampleRepeats, min_spanning: int) -> Graded | None:
    tg, t = sample.calls.get(CFTR_TG), sample.calls.get(CFTR_T)
    if tg is None or t is None or not tg.alleles or not t.alleles:
        return None
    return cftr_grade(tg, t, min_spanning=min_spanning)


def _homozygous(call: RepeatCall | None) -> bool:
    return (
        call is not None
        and len(call.alleles) == 2
        and call.alleles[0].repeats == call.alleles[1].repeats
    )


def _misleads(
    locus: RepeatLocus,
    calls: Mapping[str, RepeatCall | None],
    grades: Mapping[str, Graded | None],
) -> set[str]:
    """A non-standard locus's reasons from the trio (`PARENT_FINDING` aside): a finding on
    HG002, HG002's libraries disagreeing, a homozygous-looking HG002 size no parent pair
    carries (`NOT_INHERITED`), or no healthy sample graded normal (`NOT_VALIDATED`). An
    unassessed sample is no reason: an unassessed result never misleads."""
    kinds = {key: None if g is None else g.kind for key, g in grades.items()}
    reasons: set[str] = set()
    if any(kinds[key] is Kind.FINDING for key in _HG002):
        reasons.add(
            X_RULE_DID_NOT_HOLD if locus.rule is Rule.X_LINKED else RULE_DID_NOT_HOLD
        )
    # HG002's two libraries are one DNA sequenced twice. Its anchor comes from the other
    # library's call, which an artefact of one library's preparation cannot move because
    # the two were prepared and sequenced apart.
    sizes = [sorted(_alleles(calls[key]) or []) for key in _HG002]
    if {kinds[key] for key in _HG002} == {Kind.NON_FINDING, Kind.UNASSESSED} or (
        all(sizes)
        and (
            len(sizes[0]) != len(sizes[1])
            or any(abs(a - b) > MENDEL_TOLERANCE for a, b in zip(*sizes, strict=True))
        )
    ):
        reasons.add(LIBRARIES_DISAGREE)
    # A homozygous-looking HG002 is unassessed, so a normal at the locus rests on the
    # parents' normals; its size must be one each parent carries, or the rule's reading of
    # the trio is not the pedigree's. Its anchor comes from HG003's and HG004's own calls,
    # which an allele dropping out of HG002's reads cannot move because each parent was
    # sequenced apart; a parent with no call is no evidence, so it withholds too.
    if locus.rule is not Rule.X_LINKED and any(
        _homozygous(calls[key])
        and mendelian(
            locus.motif, locus.chrom, calls[key], calls[FATHER], calls[MOTHER]
        )
        is not True
        for key in _HG002
    ):
        reasons.add(NOT_INHERITED)
    # A rule never seen to give a normal on a healthy sample is unvalidated. Its anchor
    # comes from the calibration samples' grades, which a rule that withholds every normal
    # in calibration cannot move because a normal must appear among them; an X-linked
    # locus's one-X branch counts only the one-X samples, HG004 checking the two-X branch.
    validators = (*_HG002, FATHER) if locus.rule is Rule.X_LINKED else RECORDED
    if not any(kinds[key] is Kind.NON_FINDING for key in validators):
        reasons.add(NOT_VALIDATED)
    return reasons


def measure_repeat_rules(
    samples: Mapping[str, SampleRepeats],
    loci: Mapping[str, RepeatLocus],
    *,
    min_spanning: int,
    min_allele_reads: int,
) -> tuple[dict[str, str], dict[str, Level], dict[str, Any]]:
    """What the trio withholds (key -> one reason) and caps (key -> level), and what was
    measured per sample and locus (alleles, kind, fewest spanning reads, the Mendelian
    verdict on `CHILD`; RFC1's call; CFTR's counts). `samples` holds every `RECORDED` key;
    `ValueError` otherwise."""
    # A missing parent would pass the parent and X checks unseen. Its anchor comes from
    # `records.RECORDED`, 3c's constant, which a caller's mapping cannot move because the
    # keys are checked against the code's tuple.
    if set(samples) != set(RECORDED):
        raise ValueError(
            f"the repeat rules need {', '.join(RECORDED)}; given {', '.join(samples)}"
        )
    pairs: list[tuple[str, str]] = []
    caps: dict[str, Level] = {}
    measured: dict[str, Any] = {key: {"repeats": {}} for key in RECORDED}

    def graded(key: str, call: RepeatCall | None, locus: RepeatLocus) -> Graded | None:
        if call is None or not call.alleles:
            return None
        male = samples[key].sex == "male"
        return grade(
            locus,
            call,
            male=male,
            x_dosage=XDosage(1, 0.5, True) if male else XDosage(2, 1.0, False),
            min_allele_reads=min_allele_reads,
            min_spanning=min_spanning,
        )

    for gene, locus in loci.items():
        calls = {key: samples[key].calls.get(locus.variant_id) for key in RECORDED}
        # Each sample is graded as the run grades it: a rule that does not hold on HG002
        # is the rule the run applies.
        grades = {key: graded(key, calls[key], locus) for key in RECORDED}
        normal = {
            key: g is not None and g.kind is Kind.NON_FINDING
            for key, g in grades.items()
        }
        if any(grades[key] is None for key in _HG002):
            pairs.append((gene, NO_CALL))
        if locus.rule is not Rule.STANDARD:
            pairs += [(gene, reason) for reason in _misleads(locus, calls, grades)]
        if any(
            (g := grades[key]) is not None and g.kind is Kind.FINDING
            for key in (FATHER, MOTHER)
        ):
            pairs.append((gene, PARENT_FINDING))
        if locus.rule is Rule.X_LINKED:
            mother = calls[MOTHER]
            hg002 = [_alleles(calls[key]) or [] for key in _HG002]
            if not (
                normal[MOTHER]
                and mother is not None
                and all(_near(n, mother) for found in hg002 for n in found)
            ):
                pairs.append((two_x_key(gene), X_RULE_DID_NOT_HOLD))
        verdict = mendelian(
            locus.motif, locus.chrom, calls[CHILD], calls[FATHER], calls[MOTHER]
        )
        if verdict is False:
            caps[gene] = Level.LOW
        for key in RECORDED:
            call = calls[key]
            entry: dict[str, Any] = {
                "alleles": _alleles(call),
                "kind": _kind(grades[key]),
                "spanning": None
                if call is None or not call.alleles
                else min(a.spanning for a in call.alleles),
            }
            if key == CHILD:
                entry["mendelian"] = verdict
            measured[key]["repeats"][gene] = entry

    rfc1 = {key: samples[key].rfc1 for key in RECORDED}
    rfc1_grades = {key: None if c is None else rfc1_grade(c) for key, c in rfc1.items()}
    seen = all(
        (c := rfc1[key]) is not None and (c.call, c.description) == HG002_RFC1
        for key in _HG002
    )
    parents = all(
        (g := rfc1_grades[key]) is not None and g.kind is Kind.NON_FINDING
        for key in (FATHER, MOTHER)
    )
    if not (seen and parents):
        pairs.append((RFC1_LOCUS, RFC1_WITHHELD))
    for key, call in rfc1.items():
        measured[key]["rfc1"] = (
            None
            if call is None
            else {
                "call": call.call,
                "description": call.description,
                "kind": _kind(rfc1_grades[key]),
            }
        )

    cftr = {key: _cftr_grade(samples[key], min_spanning) for key in RECORDED}
    counts = {
        key: [_alleles(samples[key].calls.get(variant)) for variant, _ in _CFTR]
        for key in RECORDED
    }
    if any(cftr[key] is None for key in _HG002):
        pairs.append((CFTR_GENE, NO_CALL))
    elif any(
        (g := cftr[key]) is not None and g.kind is Kind.UNASSESSED for key in _HG002
    ):
        pairs.append((CFTR_GENE, RULE_DID_NOT_HOLD))
    elif cftr[CHILD] != cftr["pcr"] or [sorted(c or []) for c in counts[CHILD]] != [
        sorted(c or []) for c in counts["pcr"]
    ]:
        pairs.append((CFTR_GENE, LIBRARIES_DISAGREE))
    cftr_mendel = [
        mendelian(
            motif,
            _CFTR_CHROM,
            *(samples[key].calls.get(variant) for key in (CHILD, FATHER, MOTHER)),
        )
        for variant, motif in _CFTR
    ]
    if False in cftr_mendel:
        caps[CFTR_GENE] = Level.LOW
    for key in RECORDED:
        tg, t = counts[key]
        measured[key]["cftr"] = {"tg": tg, "t": t, "kind": _kind(cftr[key])}
    measured[CHILD]["cftr"]["mendelian"] = cftr_mendel

    withheld = merge_withheld(*({key: reason} for key, reason in pairs))
    return withheld, caps, measured
