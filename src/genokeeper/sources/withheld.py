"""What a calibration table withheld, kept for the agent (the owner's ruling, 2026-11-06:
nothing is hidden from the agent; each result says how reliable and how significant it is).

A table that withholds an item leaves its statement unassessed with the table's reason, and
`reveal` leaves that as it is: no kind, reason, validity, topic or flag moves, no gene call,
no call and no accounting, so what it keeps is never a finding, a watch item, a normal or an
exclusion, and the claim's own keys keep their nulls for every reader of them. It only adds
the claim's `WITHHELD` list, which no grading reader opens.

It grades the source a second time with nothing withheld, so the entry says what the
statement would be were the value true, in the source's own words, whatever its rules. An
item with no measured value (`_MEASURED`: no call, a tool not run) stays as it was. Four
values have no unassessed statement to reveal, and their sources write the entry themselves:
an HLA gene in the `hla:table` finding, a large depth event in the `cnv:screen` statement,
a call of a band no table measured (deletions, depth events, insertions, named loci) and
UGT1A1*28's call beside PharmCAT's gene result, which is never given it (`sources.pgx`).
"""

import json
from collections.abc import Mapping
from dataclasses import replace
from typing import Any, cast

from genokeeper.pgx.hla_typing import CONFIRM_DQ
from genokeeper.sources import (
    WITHHELD,
    Collected,
    Fields,
    Row,
    RunContext,
    Source,
    withheld_entry,
)
from genokeeper.sources.hla import CONFIRM_TYPING
from genokeeper.store.db import Kind

# The fields a source carries a table's reasons in: one reason, or item -> reason.
_REASONS = ("withheld", "heteroplasmy_withheld")
_COPIES = "a targeted copy-number test such as MLPA"
# What each claim type measures, and the test that confirms it where a claim names none
# (`confirm`; None where no test exists).
_MEASURED: dict[str, tuple[tuple[str, ...], str | None]] = {
    "repeat_expansion": (("alleles", "range", "rfc1"), "a targeted PCR test"),
    "smn_copy_number": (("smn1", "smn2", "g27134tg", "callers"), "MLPA or qPCR"),
    "y_haplogroup": (
        ("haplogroup", "marker", "qc_score", "valid_markers"),
        "Y-chromosome SNP testing in a laboratory",
    ),
    "paralog_copy_number": (("called", "callers", "variants"), _COPIES),
    "lpa_kiv2": (("kiv2",), "a one-off Lp(a) blood test"),
    "hla_allele": (("callers", "status"), CONFIRM_TYPING),
    "hla_association": (("callers", "status"), CONFIRM_TYPING),
    "hla_coeliac": (("t1k", "tags"), CONFIRM_DQ),
    "hhv6_ici": (("species", "ratio"), "ddPCR"),
    "roh_summary": (
        ("total_long", "n_long", "longest", "isodisomy"),
        "a SNP microarray",
    ),
    "blood_group": (("genotype", "phenotype"), "serology"),
    "ancestry_similarity": (("population", "low_confidence", "probabilities"), None),
    "ancestry_admixture": (("proportions",), None),
    "ancestry_neanderthal": (("neanderthal_bp", "denisovan_bp", "share"), None),
    "mitomap": (("level",), "mtDNA testing in a clinical laboratory"),
}
# What a claim holds for a caller that made no call (`hla.status`).
_NO_VALUE = (None, "unknown")
# Each caller's line, placeholders included: shown beside a value, never one alone.
_LINES = "callers"


def _measured(value: object) -> bool:
    """Whether a claim's key holds a measurement, at any depth."""
    if isinstance(value, Mapping):
        return any(_measured(v) for v in cast(Mapping[str, object], value).values())
    if isinstance(value, list):
        return any(_measured(v) for v in cast(list[object], value))
    return value not in _NO_VALUE


def _revealed(row: Row, ungated: Fields | None, reason: str | None) -> Row:
    """`row`, an unassessed statement `ungated` is otherwise and with a value, with what
    was withheld in its claim; `reason`: the table's, where the statement names another."""
    statement = row.statement
    if (
        ungated is None
        or statement["kind"] != Kind.UNASSESSED
        or (ungated["kind"], ungated["reason"])
        == (statement["kind"], statement["reason"])
    ):
        return row
    claim = json.loads(str(statement["claim"]))
    graded = json.loads(str(ungated["claim"]))
    keys, confirm = _MEASURED.get(claim["type"], ((), None))
    value = {key: graded[key] for key in keys if _measured(graded.get(key))}
    if not value.keys() - {_LINES}:
        return row
    entry = withheld_entry(
        reason or str(statement["reason"]),
        value,
        would_be=Kind(str(ungated["kind"])),
        significance=cast(str | None, ungated["reason"]),
        confirm=graded.get("confirm") or confirm,
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
    ungated = {
        str(row.statement["statement_id"]): row.statement
        for row in cast(Any, replace)(source, **{name: nothing}).collect(ctx).rows
    }
    reason = reasons if isinstance(reasons, str) else None
    rows = [
        _revealed(row, ungated.get(str(row.statement["statement_id"])), reason)
        for row in shown.rows
    ]
    return replace(shown, rows=rows)
