"""Inherited chromosomally integrated HHV-6 as a source (`Hhv6Source`; slice 3f, "HHV-6
band"). Tool `hhv6`, one statement `hhv6:ici` under wellness, never actionable, never
flagged `exclusion`.

- A ratio in the band is a finding of information at validity moderate: `ICI` with the
  species, `ICI_CONFIRM` (ddPCR, also the claim's `confirm`) and `INFORMATION_NOTE`, since
  the band's width is arithmetic (the trio offers no positive).
- A ratio below the calibration's background ceiling is a non-finding.
- A ratio between the ceiling and the band is unassessed (`BETWEEN`), so is one above the
  band (`ABOVE`) and a call without a human depth (`hhv6.NO_DEPTH`).
- `not_run` (no table, a failed step) and `withheld` (the table's reason) make it unassessed.
"""

import json
import math
from dataclasses import dataclass

from genokeeper.hhv6 import NO_DEPTH, Hhv6Call
from genokeeper.sources import INFORMATION_NOTE, Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

HHV6_TOOL = "hhv6"
ICI = (
    "inherited chromosomally integrated {species}: blood tests may read it as an active"
    " infection; each child has a 50% chance of inheriting it"
)
CONFIRM = "ddPCR"
ICI_CONFIRM = "a ddPCR test confirms it"
BETWEEN = "viral depth between background and one copy: needs a ddPCR test"
ABOVE = "viral depth above one integrated copy: needs a ddPCR test"
BELOW = "viral depth below the calibration's background"
SOURCE_ID = "NC_001664.4+NC_000898.1"


@dataclass(frozen=True)
class Hhv6Source:
    """The run's call (`hhv6.classify`; None when it has none), why HHV-6 was not run, and
    the calibration table's withholding reason."""

    call: Hhv6Call | None
    not_run: str | None
    withheld: str | None
    tool: str = HHV6_TOOL

    def _kind(self, ctx: RunContext) -> tuple[Kind, str | None]:
        call = self.call
        if self.not_run is not None:
            return Kind.UNASSESSED, self.not_run
        if self.withheld is not None:
            return Kind.UNASSESSED, self.withheld
        if call is None:
            return Kind.UNASSESSED, "not run"
        if math.isnan(call.ratio):
            return Kind.UNASSESSED, NO_DEPTH
        if call.kind is Kind.FINDING:
            ici = ICI.format(species=call.species)
            return Kind.FINDING, f"{ici}; {ICI_CONFIRM}; {INFORMATION_NOTE}"
        if call.kind is Kind.NON_FINDING:
            return Kind.NON_FINDING, BELOW
        above = call.ratio > ctx.grading.hhv6.max_ratio
        return Kind.UNASSESSED, ABOVE if above else BETWEEN

    def collect(self, ctx: RunContext) -> Collected:
        """The one `hhv6:ici` statement."""
        kind, reason = self._kind(ctx)
        call = self.call
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "hhv6_ici",
            "species": None if call is None else call.species,
            "ratio": None if call is None or math.isnan(call.ratio) else call.ratio,
            "confirm": CONFIRM,
        }
        ratio = claim["ratio"]
        measured = (
            None if ratio is None else {"species": claim["species"], "ratio": ratio}
        )
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{HHV6_TOOL}:ici",
            "kind": str(kind),
            "vid": None,
            "gene": None,
            "claim": json.dumps(claim, ensure_ascii=False),
            "validity": str(Validity.MODERATE),
            "source_level": "viral depth over the haploid human depth",
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.WELLNESS),
            "source": HHV6_TOOL,
            "source_id": SOURCE_ID,
            "source_version": SOURCE_ID,
            "points": None,
            "exclusion": False,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement, measured)],
            extra={},
            input_rows=1,
            dropped={},
        )
