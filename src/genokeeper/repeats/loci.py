"""The repeat loci's thresholds and ExpansionHunter's calls as the grading reads them.

The thresholds are STRchive v2.28.0's, written into `[repeat_loci]` of config/grading.toml;
a repeat count between two of a locus's classes is a `GAP`, which the grading treats as
borderline. Each locus is graded by its `Rule`: slice 3a's eleven and TCF4 and PABPN1 by
the standard one, the rest (slice 3e) exclusion-only.
"""

from dataclasses import dataclass
from enum import StrEnum

from genokeeper.store.db import Topic, Validity


class RepeatClass(StrEnum):
    """Where a repeat count falls among a locus's ranges."""

    NORMAL = "normal"
    INTERMEDIATE = "intermediate"
    REDUCED_PENETRANCE = "reduced penetrance"
    PATHOGENIC = "pathogenic"
    GAP = "gap"


class Rule(StrEnum):
    """How a locus is graded: `STANDARD` sizes it (slice 3a); `EXCLUSION` says only
    "normal" or unassessed, or "possible" above the disease minimum; `X_LINKED` is that by
    the run's X dosage and never a finding; `DISTINCT` (FGF14) says only "normal", from two
    distinct alleles."""

    STANDARD = "standard"
    EXCLUSION = "exclusion"
    X_LINKED = "x_linked"
    DISTINCT = "distinct"


@dataclass(frozen=True)
class RepeatLocus:
    """One locus: its catalog ids, motif and chromosome, and its ranges in repeat units
    (inclusive); `evidence` is STRchive's `evidence` field. `rule` grades it, `validity`
    overrides the evidence's mapping, and `recessive_at` is the count from which both
    alleles below the disease minimum are a possible recessive finding."""

    gene: str
    locus_id: str
    variant_id: str
    chrom: str
    motif: str
    normal_max: int
    intermediate: tuple[int, int] | None
    reduced_penetrance: tuple[int, int] | None
    pathogenic_min: int
    evidence: str
    topic: Topic
    note: str | None
    rule: Rule = Rule.STANDARD
    validity: Validity | None = None
    recessive_at: int | None = None


def _within(span: tuple[int, int] | None, repeats: int) -> bool:
    return span is not None and span[0] <= repeats <= span[1]


def classify(locus: RepeatLocus, repeats: int) -> RepeatClass:
    """The class of a repeat count; a count in none of the ranges is a `GAP`."""
    if repeats <= locus.normal_max:
        return RepeatClass.NORMAL
    if _within(locus.intermediate, repeats):
        return RepeatClass.INTERMEDIATE
    if _within(locus.reduced_penetrance, repeats):
        return RepeatClass.REDUCED_PENETRANCE
    if repeats >= locus.pathogenic_min:
        return RepeatClass.PATHOGENIC
    return RepeatClass.GAP


@dataclass(frozen=True)
class Allele:
    """One allele: ExpansionHunter's repeat count and its interval, and the reads that
    support it (spanning, flanking, in-repeat)."""

    repeats: int
    ci: tuple[int, int]
    spanning: int
    flanking: int
    inrepeat: int


@dataclass(frozen=True)
class RepeatCall:
    """ExpansionHunter's call at one variant: one or two alleles (one for a hemizygous
    locus), the locus's coverage, ExpansionHunter's `LowDepth` filter, the read length, the
    most units any flanking read holds and the variant's in-repeat reads (its JSON's
    `CountsOfInrepeatReads`, which, unlike an allele's `ADIR`, counts every one)."""

    variant_id: str
    alleles: tuple[Allele, ...]
    coverage: float
    low_depth: bool
    read_length: int
    longest_flanking: int = 0
    inrepeat: int = 0


@dataclass(frozen=True)
class Rfc1Call:
    """ExpansionHunter's RFC1 motif analysis: its label (`call`) and `description`, the
    locus's alleles, its coverage and ExpansionHunter's `LowDepth` filter."""

    call: str
    description: str
    alleles: tuple[Allele, ...]
    coverage: float
    low_depth: bool
