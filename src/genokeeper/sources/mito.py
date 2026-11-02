"""The mitochondrial haplogroup as a source (`MitoSource`; track design, slice 1d,
"Haplogroup").

One statement with neither vid nor gene (`haplogrep3:haplogroup`) under ancestry: a finding
with Haplogrep3's haplogroup, its quality and tree, and the chrM positions it read as rCRS
for want of reads. Its validity follows the quality (builder's call, the plan's cut-offs):
strong at 0.9 and above, moderate at 0.8, limited below. When too much of chrM is uncovered
(`haplogroup.too_uncovered`) no haplogroup is stated: the claim is unassessed, its reason
the coverage.
"""

import json
from collections.abc import Sequence
from dataclasses import dataclass

from genokeeper.haplogroup import HAPLOGREP_VERSION, RCRS_LENGTH, TREE, Haplogroup
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

MITO_TOOL = "haplogrep3"


def haplogroup_validity(quality: float) -> Validity:
    """Strong at a Haplogrep3 quality of 0.9 and above, moderate at 0.8, else limited."""
    if quality >= 0.9:
        return Validity.STRONG
    return Validity.MODERATE if quality >= 0.8 else Validity.LIMITED


@dataclass(frozen=True)
class MitoSource:
    """Haplogrep3's haplogroup (None: not stated, chrM too uncovered), the chrM positions
    below `min_dp` reads, and that depth."""

    haplogroup: Haplogroup | None
    uncovered: Sequence[int]
    min_dp: int
    tool: str = MITO_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """The haplogroup's one statement."""
        found = self.haplogroup
        holes = len(self.uncovered)
        if found is None:
            validity, level = Validity.LIMITED, "Haplogrep3 not run"
            reason = (
                f"{holes} of {RCRS_LENGTH} chrM positions below {self.min_dp} reads;"
                " no haplogroup stated"
            )
        else:
            validity = haplogroup_validity(found.quality)
            level = f"Haplogrep3 quality {found.quality:.2f}"
            reason = (
                f"{holes} chrM position{'s' * (holes != 1)} below {self.min_dp} reads"
                " read as rCRS"
                if holes
                else None
            )
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "mito_haplogroup",
            "haplogroup": None if found is None else found.name,
            "quality": None if found is None else found.quality,
            "tree": TREE,
            "uncovered_positions": None if found is None else list(self.uncovered),
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{MITO_TOOL}:haplogroup",
            "kind": str(Kind.UNASSESSED if found is None else Kind.FINDING),
            "vid": None,
            "gene": None,
            "claim": json.dumps(claim),
            "validity": str(validity),
            "source_level": level,
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.ANCESTRY),
            "source": MITO_TOOL,
            "source_id": TREE,
            "source_version": HAPLOGREP_VERSION,
            "points": None,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement)],
            extra={},
            input_rows=1,
            dropped={},
        )
