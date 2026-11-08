"""What a calibration table withheld, kept for the agent (the owner's ruling, 2026-11-06:
nothing is hidden from the agent; each result says how reliable and how significant it is).

A table that withholds an item leaves its statement unassessed with the table's reason, and
`reveal` leaves that as it is: no kind, reason, validity, topic or flag moves, no gene call,
no call and no accounting, so what it keeps is never a finding, a watch item, a normal or an
exclusion, and the claim's own keys keep their nulls for every reader of them. It only adds
the claim's `WITHHELD` list, which no grading reader opens.

Each source hands over what its tools measured as typed data beside the statement
(`Row.measured`, per caller where there are several, a caller with nothing absent), so whether
a value exists is never read off display text, and that data reaches the store only inside
the entry. No measured value (no call, a tool not run): the item stays as it was. `reveal`
grades the source a second time with nothing withheld, so the entry says what the statement
would be were the value true, in the source's own words, whatever its rules. Five
values have no unassessed statement to reveal, and their sources write the entry themselves:
an HLA gene in the `hla:table` finding, a large depth event in the `cnv:screen` statement,
a call of a band no table measured (deletions, depth events, insertions, named loci),
UGT1A1*28's call beside PharmCAT's gene result, which is never given it (`sources.pgx`), and
a polygenic score below the overlap (`sources.prs`), which no table withheld.
"""

import json
from collections.abc import Mapping
from dataclasses import replace
from typing import Any, cast

from genokeeper.pgx.hla_typing import CONFIRM_DQ
from genokeeper.sources import Collected, Row, RunContext, Source
from genokeeper.sources.hla import CONFIRM_TYPING
from genokeeper.sources.panel import WITHHELD, withheld_entry
from genokeeper.store.db import Kind

# The fields a source carries a table's reasons in: one reason, or item -> reason.
_REASONS = ("withheld", "heteroplasmy_withheld")
# The test that confirms what each claim type measures, where a claim names none (`confirm`;
# none for ancestry: no test exists).
_CONFIRM = {
    "repeat_expansion": "a targeted PCR test",
    "smn_copy_number": "MLPA or qPCR",
    "y_haplogroup": "Y-chromosome SNP testing in a laboratory",
    "paralog_copy_number": "a targeted copy-number test such as MLPA",
    "lpa_kiv2": "a one-off Lp(a) blood test",
    "hla_allele": CONFIRM_TYPING,
    "hla_association": CONFIRM_TYPING,
    "hla_coeliac": CONFIRM_DQ,
    "hhv6_ici": "ddPCR",
    "roh_summary": "a SNP microarray",
    "blood_group": "serology",
    "mitomap": "mtDNA testing in a clinical laboratory",
}


def _revealed(row: Row, ungated: Row, reason: str | None) -> Row:
    """`row`, an unassessed statement `ungated` is otherwise and with a value, with what
    was withheld in its claim; `reason`: the table's, where the statement names another."""
    statement, graded = row.statement, ungated.statement
    if (
        not row.measured
        or statement["kind"] != Kind.UNASSESSED
        or (graded["kind"], graded["reason"])
        == (statement["kind"], statement["reason"])
    ):
        return row
    claim = json.loads(str(statement["claim"]))
    entry = withheld_entry(
        reason or str(statement["reason"]),
        row.measured,
        would_be=Kind(str(graded["kind"])),
        significance=cast(str | None, graded["reason"]),
        confirm=json.loads(str(graded["claim"])).get("confirm")
        or _CONFIRM.get(claim["type"]),
    )
    held = json.dumps({**claim, WITHHELD: [entry]}, ensure_ascii=False)
    return replace(row, statement={**statement, "claim": held})


def reveal(source: Source, ctx: RunContext) -> Collected:
    """What `source` collects, with what a table withheld in the claims it left unassessed
    (the module docstring)."""
    shown = source.collect(ctx)
    name = next((n for n in _REASONS if getattr(source, n, None)), None)
    if name is None:
        return shown
    reasons = getattr(source, name)
    nothing = {} if isinstance(reasons, Mapping) else None
    ungated = cast(Any, replace)(source, **{name: nothing}).collect(ctx).rows
    reason = reasons if isinstance(reasons, str) else None
    rows = [
        _revealed(row, graded, reason)
        for row, graded in zip(shown.rows, ungated, strict=True)
    ]
    return replace(shown, rows=rows)
