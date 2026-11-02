"""CYP2B6 from Aldy alone: its call, read from its output, as PharmCAT's outside call.

`fixtures/aldy_cyp2b6.aldy` is written in the shape Aldy 4.8.3 writes (`aldy/genotype.py`
and `diplotype.write_decomposition`; CYP2B6 has a CPIC function table and no scores, so its
solution line states `cpic=` alone): a *1/*6 sample with *6's two defining variants. Aldy's
own test data holds only chr22, so no CYP2B6 output of a real run is at hand.
"""

from dataclasses import replace
from pathlib import Path

from genokeeper.pgx.cyp2b6 import cyp2b6_outside
from genokeeper.pgx.cyp2d6 import read_aldy
from genokeeper.pgx.outside import outside_lines

FIXTURE = Path(__file__).parent / "fixtures" / "aldy_cyp2b6.aldy"


def test_aldy_cyp2b6_parsed() -> None:
    """Aldy's one solution is the gene's call, from one caller: never agreed (medium in
    the gene call's grading), written to PharmCAT as its diplotype."""
    call = read_aldy(FIXTURE)
    assert (call.caller, call.diplotype, call.phenotype, call.activity_score) == (
        "Aldy",
        "*1/*6",
        "intermediate",
        None,
    )
    outside = cyp2b6_outside(call)
    assert (outside.gene, outside.diplotype, outside.activity_score) == (
        "CYP2B6",
        "*1/*6",
        None,
    )
    assert (outside.agreed, outside.unresolved_reason) == (False, None)
    assert outside.callers == {"Aldy": "*1/*6"}
    assert outside_lines([outside]) == "CYP2B6\t*1/*6\n"


def test_aldy_cyp2b6_sub_alleles_and_alternatives() -> None:
    """Sub-allele numbers go (PharmCAT names star alleles); solutions naming one diplotype
    are one call; solutions naming several, or none, leave the gene unresolved."""
    call = replace(read_aldy(FIXTURE), diplotype="*6.001/*1.002 or *1.001/*6.002")
    assert cyp2b6_outside(call).diplotype == "*1/*6"
    split = cyp2b6_outside(replace(call, diplotype="*1/*6 or *4/*9"))
    assert (split.diplotype, split.unresolved_reason) == (
        None,
        "Aldy's solutions name different diplotypes: *1/*6 or *4/*9",
    )
    assert split.callers == {"Aldy": "*1/*6 or *4/*9"}
    nocall = cyp2b6_outside(replace(call, diplotype=None))
    assert (nocall.unresolved_reason, nocall.callers) == (
        "Aldy made no call",
        {"Aldy": "no call"},
    )
    assert outside_lines([split]) == outside_lines([nocall]) == ""
