"""SMN1 copy number as a source (`SmnSource`; slice 3a, "Grading", SMN1/SMN2).

One statement (`smn:SMN1`) under carrier and one `gene_call`. A finding needs one SMN1 copy
from both SMNCopyNumberCaller and Parascopy, worded "screening result, needs confirmation";
the tools disagreeing, a call missing, a tool not run or withheld, or no SMN1 copy at all
(a calling error in a healthy adult) is unassessed, never normal. SMN2's copy number is
shown as information. Nothing here is actionable.
"""

import json
from dataclasses import dataclass

from genokeeper.evidence.pgx import GENE_UNRESOLVED
from genokeeper.smn import (
    PARASCOPY_VERSION,
    SMNCALLER_VERSION,
    ParascopyCall,
    SmnCaller,
)
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.sources.repeats import CONFIRM
from genokeeper.store.db import Kind, Level, Topic, Validity

SMN_TOOL = "smn"
SMN_GENE = "SMN1"
SMN_CONDITION = "spinal muscular atrophy"
CALLER_NAME = f"SMNCopyNumberCaller {SMNCALLER_VERSION}"
PARASCOPY_NAME = f"Parascopy {PARASCOPY_VERSION}"
FINDING_REASON = f"{CONFIRM} by MLPA or qPCR: one SMN1 copy, SMA carrier"
TWO_COPIES_REASON = (
    "reduced risk; a 2+0 carrier (both copies on one chromosome) is not excluded"
)
SILENT_CARRIER_REASON = f"{CONFIRM}: possible silent carrier (g.27134T>G)"
ZERO_REASON = "SMN1 = 0 in a healthy adult: a calling error"


def _shown(value: float | None) -> str:
    return "?" if value is None else str(value)


@dataclass(frozen=True)
class SmnSource:
    """Both callers' results (None: the tool gave none), why the tools did not run
    (`not_run`) or the calibration withholds SMN (`withheld`), and the Parascopy quality
    below which its call does not count."""

    caller: SmnCaller | None
    parascopy: ParascopyCall | None
    not_run: str | None
    withheld: str | None
    min_qual: float
    tool: str = SMN_TOOL

    def _parascopy_smn1(self) -> int | None:
        p = self.parascopy
        if p is None or not p.passed or p.qual is None or p.qual < self.min_qual:
            return None
        return p.smn1

    def _grade(self) -> tuple[Kind, str | None, str | None]:
        """The statement's kind, its reason and the gene call's phenotype."""
        if self.not_run is not None:
            return Kind.UNASSESSED, self.not_run, None
        if self.withheld is not None:
            return Kind.UNASSESSED, self.withheld, None
        a = None if self.caller is None else self.caller.smn1
        b = self._parascopy_smn1()
        if a is None:
            return Kind.UNASSESSED, f"{CALLER_NAME} gave no SMN1 copy number", None
        if b is None:
            low = f"failed, or quality below {self.min_qual:g}"
            reason = f"{PARASCOPY_NAME} gave no usable SMN1 copy number ({low})"
            return Kind.UNASSESSED, reason, None
        if a != b:
            return (
                Kind.UNASSESSED,
                f"the two callers disagree on SMN1 ({a} and {b})",
                None,
            )
        if a == 0:
            return Kind.UNASSESSED, ZERO_REASON, None
        if a == 1:
            return Kind.FINDING, FINDING_REASON, "SMA carrier (screening)"
        marker = self.caller.g27134tg if self.caller is not None else None
        if a == 2:
            if marker is not None and marker >= 1:
                return Kind.WATCH, SILENT_CARRIER_REASON, "possible silent carrier"
            return Kind.NON_FINDING, TWO_COPIES_REASON, "reduced risk"
        return Kind.NON_FINDING, None, "not a carrier"

    def collect(self, ctx: RunContext) -> Collected:
        """The SMN1 statement and its gene call."""
        kind, reason, phenotype = self._grade()
        caller, parascopy = self.caller, self.parascopy
        smn1 = None if caller is None else caller.smn1
        smn2 = None if caller is None else caller.smn2
        assessed = kind is not Kind.UNASSESSED
        callers: dict[str, str] = {}
        if caller is not None:
            callers[CALLER_NAME] = (
                f"SMN1 {_shown(caller.smn1)}, SMN2 {_shown(caller.smn2)},"
                f" g.27134T>G {_shown(caller.g27134tg)}"
            )
        if parascopy is not None:
            callers[PARASCOPY_NAME] = (
                f"SMN1 {_shown(parascopy.smn1)}, SMN2 {_shown(parascopy.smn2)},"
                f" qual {_shown(parascopy.qual)}, {'pass' if parascopy.passed else 'fail'}"
            )
        other = self._parascopy_smn1()
        agreed = None if smn1 is None or other is None else smn1 == other
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "smn_copy_number",
            "gene": SMN_GENE,
            "condition": SMN_CONDITION,
            "smn1": smn1 if assessed else None,
            "smn2": smn2,
            "g27134tg": None if caller is None else caller.g27134tg,
            "callers": callers,
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{SMN_TOOL}:{SMN_GENE}",
            "kind": str(kind),
            "vid": None,
            "gene": SMN_GENE,
            "claim": json.dumps(claim),
            "validity": str(Validity.STRONG),
            "source_level": f"{CALLER_NAME} and {PARASCOPY_NAME}",
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.CARRIER),
            "source": SMN_TOOL,
            "source_id": SMN_GENE,
            "source_version": f"{SMNCALLER_VERSION}; {PARASCOPY_VERSION}",
            "points": None,
        }
        gene_call: Fields = {
            "run_id": ctx.run_id,
            "gene": SMN_GENE,
            "diplotype": f"SMN1 {_shown(smn1)}, SMN2 {_shown(smn2)}"
            if assessed
            else None,
            "phenotype": phenotype,
            "activity_score": None,
            "callers": json.dumps(callers),
            "agreed": agreed,
            "level": str(Level.MEDIUM if assessed else Level.LOW),
            "warnings": []
            if assessed or reason is None
            else [f"{GENE_UNRESOLVED}: {reason}"],
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement)],
            extra={"gene_call": [gene_call]},
            input_rows=1,
            dropped={},
        )
