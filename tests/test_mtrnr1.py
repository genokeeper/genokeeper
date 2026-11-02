"""MT-RNR1 from mutserve's heteroplasmy levels at CPIC's MT-RNR1 variants."""

from genokeeper.mito import MitoCall
from genokeeper.pgx.mtrnr1 import mtrnr1_outside
from genokeeper.pgx.outside import outside_lines

# CPIC's MT-RNR1 allele names: `Reference` and one per variant.
ALLELES = ["Reference", "m.827A>G", "m.961T>del", "m.1494C>T", "m.1555A>G"]
MIN_HOMOPLASMY = 0.9
# A haplogroup H variant outside MT-RNR1 and one inside it CPIC does not list.
OTHERS = [MitoCall(73, "A", "G", 1.0, 3100), MitoCall(750, "A", "G", 0.99, 2900)]


def test_mtrnr1_homoplasmic_1555() -> None:
    calls = [*OTHERS, MitoCall(1555, "A", "G", 0.998, 3050)]
    outside = mtrnr1_outside(calls, ALLELES, min_homoplasmy=MIN_HOMOPLASMY)
    assert (outside.gene, outside.diplotype, outside.unresolved_reason) == (
        "MT-RNR1",
        "m.1555A>G",
        None,
    )
    assert (outside.agreed, outside.callers) == (
        False,
        {"mutserve": "m.1555A>G at 0.998"},
    )
    assert outside_lines([outside]) == "MT-RNR1\tm.1555A>G\n"
    # No CPIC variant called: the reference, PharmCAT's `Reference`.
    reference = mtrnr1_outside(OTHERS, ALLELES, min_homoplasmy=MIN_HOMOPLASMY)
    assert (reference.diplotype, reference.unresolved_reason) == ("Reference", None)
    assert reference.callers == {"mutserve": "no CPIC MT-RNR1 variant"}


def test_mtrnr1_heteroplasmic_unresolved() -> None:
    """A defining variant below min_homoplasmy is neither carried nor absent: the gene is
    unresolved with its level stated, and PharmCAT gets no line for it."""
    calls = [*OTHERS, MitoCall(1555, "A", "G", 0.42, 2980)]
    outside = mtrnr1_outside(calls, ALLELES, min_homoplasmy=MIN_HOMOPLASMY)
    assert outside.diplotype is None
    assert outside.unresolved_reason == ("m.1555A>G at heteroplasmy 0.42, below 0.9")
    assert outside.callers == {"mutserve": "m.1555A>G at 0.42"}
    assert outside_lines([outside]) == ""


def test_mtrnr1_two_variants_or_no_coverage_unresolved() -> None:
    """MT-RNR1 is one allele: two CPIC variants carried cannot be written as one. Without
    reads at a CPIC position (mutserve's raw coverage), the reference is not asserted."""
    calls = [MitoCall(827, "A", "G", 1.0, 3000), MitoCall(1555, "A", "G", 0.99, 3000)]
    two = mtrnr1_outside(calls, ALLELES, min_homoplasmy=MIN_HOMOPLASMY)
    assert two.unresolved_reason == "more than one CPIC variant: m.827A>G, m.1555A>G"
    coverage = {827: 3000, 1494: 3000, 1555: 4}
    thin = mtrnr1_outside(
        OTHERS, ALLELES, min_homoplasmy=MIN_HOMOPLASMY, coverage=coverage, min_dp=10
    )
    assert thin.unresolved_reason == "fewer than 10 reads at m.1555"
    covered = {827: 3000, 1494: 3000, 1555: 3000}
    assert mtrnr1_outside(
        OTHERS, ALLELES, min_homoplasmy=MIN_HOMOPLASMY, coverage=covered, min_dp=10
    ).diplotype == ("Reference")


def test_mtrnr1_variant_on_too_few_reads_unresolved() -> None:
    calls = [*OTHERS, MitoCall(1555, "A", "G", 1.0, 3)]
    outside = mtrnr1_outside(calls, ALLELES, min_homoplasmy=MIN_HOMOPLASMY, min_dp=10)
    assert outside.diplotype is None
    assert outside.unresolved_reason == "m.1555A>G on 3 reads, fewer than 10"
