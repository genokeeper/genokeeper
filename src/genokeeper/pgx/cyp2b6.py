"""CYP2B6 from Aldy alone (track design, slice 1b, "Outside calls"; design decision,
2026-10-19: Aldy called 87% of CYP2B6 correctly against about 48% for the others).

Aldy runs as for CYP2D6 (`cyp2d6.run_aldy(..., "CYP2B6")`). Its licence (IURTC) permits
non-commercial research only and excludes diagnostic use (README). CPIC gives CYP2B6 no
activity scores, so PharmCAT gets the diplotype, named as PharmCAT names star alleles
(`cyp2d6.canonical_diplotype`: no sub-allele numbers), and looks its phenotype up. One
caller is never agreed: the gene call's reliability is medium at best.
"""

from genokeeper.pgx.cyp2d6 import alternatives, canonical_diplotype
from genokeeper.pgx.outside import CallerCall, OutsideCall

GENE = "CYP2B6"


def cyp2b6_outside(call: CallerCall) -> OutsideCall:
    """Aldy's CYP2B6 call as PharmCAT's outside call: its diplotype where every solution
    names the same one; otherwise, or with no call, unresolved with Aldy's call shown."""
    callers = {call.caller: call.diplotype or "no call"}
    reason = None
    diplotype = None
    if call.diplotype is None:
        reason = f"{call.caller} made no call"
    else:
        names = {canonical_diplotype(d) for d in alternatives(call.diplotype)}
        if None in names:
            reason = f"{call.caller}'s call names no diplotype PharmCAT reads: {call.diplotype}"
        elif len(names) > 1:
            reason = (
                f"{call.caller}'s solutions name different diplotypes: {call.diplotype}"
            )
        else:
            diplotype = names.pop()
    return OutsideCall(
        gene=GENE,
        diplotype=diplotype,
        phenotype=None,
        activity_score=None,
        callers=callers,
        agreed=False,
        unresolved_reason=reason,
    )
