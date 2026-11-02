"""Outside calls: a gene's call from callers other than PharmCAT (CYP2D6 from BCyrius and
Aldy, CYP2B6 from Aldy, HLA-A and HLA-B from OptiType and T1K, MT-RNR1 from mutserve), and
the file PharmCAT reads them from (`-po`), overriding its own (track design, slice 1b).

`OutsideCall` is the one type every outside caller returns. An unresolved call (the callers
disagree, or a heteroplasmy too low to call) keeps every caller's call in `callers` and its
reason in `unresolved_reason`, and writes no line: PharmCAT then has no call for the gene,
and the store holds it as a gene-level unassessed statement.
"""

import re
from collections.abc import Sequence
from dataclasses import dataclass

# A gene copy number written `x<N>`: CPIC and PharmCAT bin three or more copies as `x≥3`.
_COPIES = re.compile(r"x(\d+)\b")
_BIN = "x≥3"
# CPIC's activity score of three or more copies of a function allele (`*1x≥3`): one bin.
BIN_SCORE = "≥3.0"
# The core gate's label for every outside caller's run, as for PharmCAT's.
LABEL = "dna:pgx"


class CallerFailed(RuntimeError):
    """An outside caller ran but wrote no result."""


@dataclass(frozen=True)
class OutsideCall:
    """One gene's call from outside PharmCAT: its diplotype (or single allele), phenotype
    and CPIC activity score as CPIC writes it (`cpic_score`), each optional; every caller's
    own call (caller -> call); whether two callers agreed; why it is unresolved, None when it
    is not; and a warning a resolved call carries for its reader (`note`), such as resting on
    one caller's score."""

    gene: str
    diplotype: str | None
    phenotype: str | None
    activity_score: str | None
    callers: dict[str, str]
    agreed: bool
    unresolved_reason: str | None
    note: str | None = None


@dataclass(frozen=True)
class CallerCall:
    """One outside caller's own call of a gene: the caller, the diplotype as it writes it
    (None for its no-call), and its phenotype and CPIC activity score as it states them
    (the score as `cpic_score` writes it). `abstained` is true when it names a diplotype
    and states no score for it (BCyrius's `-` or `n/a`, Aldy's solutions all without
    `cpic_score`), not when it states several. `filtered` is the caller's own non-PASS filter on its call."""

    caller: str
    diplotype: str | None
    phenotype: str | None
    activity_score: str | None
    abstained: bool = False
    filtered: str | None = None


def _binned(diplotype: str) -> str:
    """`diplotype` with each copy number of three or more written CPIC's way (`*1x3` ->
    `*1x≥3`); PharmCAT 3.4.0 knows the bin, not the count (`*1x3` and `*1≥3` are
    undocumented alleles to it, phenotyped Indeterminate)."""
    return _COPIES.sub(
        lambda m: _BIN if int(m.group(1)) >= 3 else m.group(0), diplotype
    )


def cpic_score(stated: str) -> str | None:
    """A caller's stated CPIC activity score as CPIC writes it (`1.0`, `0.25`), a score of 3
    or more as CPIC's bin `BIN_SCORE` (`3.0`, `4.0`, `≥3.0`): three or more copies of a
    function allele are one bin, so 3 and 4 copies agree. None for no score (`-`, `n/a`,
    several)."""
    text = stated.strip()
    try:
        value = float(text.removeprefix("≥"))
    except ValueError:
        return None
    if value >= 3:
        return BIN_SCORE
    if text.startswith("≥"):
        return None
    text = f"{value:.2f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def outside_lines(calls: Sequence[OutsideCall]) -> str:
    """PharmCAT's outside-call file for `calls`: one line per resolved gene, tab-separated
    gene, diplotype, phenotype and activity score, an empty column keeping its tab and
    trailing empty columns left off. A diplotype of three or more copies is written in its
    bin (`_binned`), with its score: given a score, PharmCAT phenotypes even a diplotype it
    does not know. A phenotype is written only for a call without a diplotype: given both,
    PharmCAT looks recommendations up by the phenotype as written, and its run found none
    for an HLA gene's (an empty lookup key) where the diplotype alone gave each allele's
    status. Two calls of one gene are refused (`ValueError`). Write it as UTF-8."""
    genes = [c.gene for c in calls]
    repeated = sorted({g for g in genes if genes.count(g) > 1})
    if repeated:
        raise ValueError(f"more than one outside call for {', '.join(repeated)}")
    lines: list[str] = []
    for call in calls:
        if call.unresolved_reason is not None:
            continue
        diplotype = None if call.diplotype is None else _binned(call.diplotype)
        phenotype = None if diplotype else call.phenotype
        columns = [
            call.gene,
            diplotype or "",
            phenotype or "",
            call.activity_score or "",
        ]
        while len(columns) > 2 and not columns[-1]:
            columns.pop()
        lines.append("\t".join(columns) + "\n")
    return "".join(lines)
