"""MT-RNR1 from mutserve's heteroplasmy levels (track design, slice 1b, "Outside calls").

PharmCAT reads no MT-RNR1 position from a VCF; its outside call names one allele, as CPIC
names it (`m.1555A>G`; `evidence.pgx.load_cpic_alleles`), or `Reference`. mutserve (2a, `mito.call_chrm`)
gives each chrM variant's level:

- a CPIC variant at a level of at least [mito] min_homoplasmy, on at least `min_dp` reads,
  is the allele carried;
- one below either is neither carried nor absent: the gene is unresolved, its level or reads
  stated;
- two CPIC variants carried cannot be one allele: unresolved, both named;
- none is `Reference`, and, given mutserve's raw coverage, only where every CPIC
  single-base position has at least `min_dp` reads; otherwise unresolved.

mutserve calls no indels, so CPIC's deletions at m.960 and m.961 (all of uncertain risk) are
never seen: `Reference` means none of CPIC's single-base variants. One caller is never agreed.
"""

import re
from collections.abc import Collection, Mapping, Sequence

from genokeeper.mito import MitoCall
from genokeeper.pgx.outside import OutsideCall

GENE = "MT-RNR1"
CALLER = "mutserve"
REFERENCE = "Reference"
# A single-base variant as CPIC names it.
_SNV = re.compile(r"^m\.(\d+)([ACGT])>([ACGT])$")


def mtrnr1_outside(
    calls: Sequence[MitoCall],
    alleles: Collection[str],
    *,
    min_homoplasmy: float,
    coverage: Mapping[int, int] | None = None,
    min_dp: int = 1,
) -> OutsideCall:
    """MT-RNR1's outside call from mutserve's `calls`, by CPIC's MT-RNR1 allele names
    `alleles` (see the module docstring); `coverage` is mutserve's reads per position, if
    known."""
    defined = {
        (int(m[1]), m[2], m[3]): name
        for name in alleles
        if (m := _SNV.match(name)) is not None
    }
    found = [
        (defined[(c.pos, c.ref, c.alt)], c)
        for c in calls
        if (c.pos, c.ref, c.alt) in defined
    ]
    callers = {
        CALLER: ", ".join(f"{name} at {c.level:.3g}" for name, c in found)
        or "no CPIC MT-RNR1 variant"
    }
    low = [(n, c) for n, c in found if c.level < min_homoplasmy or c.dp < min_dp]
    thin = sorted(
        {
            pos
            for pos, _, _ in defined
            if coverage is not None and coverage.get(pos, 0) < min_dp
        }
    )
    reason = None
    if low:
        reason = "; ".join(
            f"{name} on {c.dp} reads, fewer than {min_dp}"
            if c.dp < min_dp
            else f"{name} at heteroplasmy {c.level:.3g}, below {min_homoplasmy:g}"
            for name, c in low
        )
    elif len(found) > 1:
        reason = f"more than one CPIC variant: {', '.join(name for name, _ in found)}"
    elif not found and thin:
        reason = f"fewer than {min_dp} reads at {', '.join(f'm.{p}' for p in thin)}"
    allele = None if reason is not None else (found[0][0] if found else REFERENCE)
    return OutsideCall(
        gene=GENE,
        diplotype=allele,
        phenotype=None,
        activity_score=None,
        callers=callers,
        agreed=False,
        unresolved_reason=reason,
    )
