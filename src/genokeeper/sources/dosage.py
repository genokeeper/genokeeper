"""Sex-chromosome dosage as a source (`DosageSource`; slice 3c, the closing design's
"Sex-chromosome dosage").

One statement with neither vid nor gene (`dosage:sex-chromosomes`): the pattern
`sex.karyotype` names from the X and Y depth ratios and the chrX allele balance, a finding
or a non-finding, never actionable (every finding but the mosaic loss of Y is a "screening
result, needs a karyotype"). The positives are shown, not calibrated: no truth sample
carries one. The statement is unassessed with `not_run`: the run passes the dosage step's
failure, or the `foundation` table's own reason when that table is not current
(`tables.open_table`), since the `[dosage]` bounds are set from it; `[dosage] calibrated_on`
is only the statement's `source_version`.
"""

import json
from dataclasses import dataclass

from genokeeper.config import DosageThresholds
from genokeeper.sex import Dosage, karyotype
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

DOSAGE_TOOL = "dosage"
STATEMENT_ID = f"{DOSAGE_TOOL}:sex-chromosomes"
SOURCE_LEVEL = "mosdepth depth ratios and chrX allele balance"


@dataclass(frozen=True)
class DosageSource:
    """The measured dosage (None: not measured), why it was not (`not_run`), and the
    `[dosage]` thresholds that grade it."""

    dosage: Dosage | None
    not_run: str | None
    thresholds: DosageThresholds
    tool: str = DOSAGE_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """The dosage's one statement."""
        d = self.dosage
        reason = self.not_run
        if reason is None and d is None:
            reason = "sex-chromosome dosage not measured"
        found = (
            None if reason is not None or d is None else karyotype(d, self.thresholds)
        )
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "sex_chromosome_dosage",
            "label": None if found is None else found.label,
            "x_ratio": None if d is None else d.x_ratio,
            "y_ratio": None if d is None else d.y_ratio,
            "x_het_share": None if d is None else d.x_het_share,
            "x_sites": None if d is None else d.x_sites,
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": STATEMENT_ID,
            "kind": str(Kind.UNASSESSED if found is None else found.kind),
            "vid": None,
            "gene": None,
            "claim": json.dumps(claim),
            "validity": str(Validity.MODERATE if found is None else found.validity),
            "source_level": SOURCE_LEVEL,
            "reason": reason if found is None else found.reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.RISKS if found is None else found.topic),
            "source": DOSAGE_TOOL,
            "source_id": "sex-chromosomes",
            "source_version": self.thresholds.calibrated_on,
            "points": None,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement)],
            extra={},
            input_rows=1,
            dropped={},
        )
