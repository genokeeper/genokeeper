"""The Y haplogroup as a source (`YleafSource`; slice 3a, "Grading", Y haplogroup).

One statement with neither vid nor gene (`yleaf:haplogroup`) under ancestry, as
`MitoSource` writes the mitochondrial one: a finding with Yleaf's haplogroup, its marker
and QC score. Validity follows the QC score: strong at 0.99 and above, moderate at the
configured minimum. A sample that is not male has no Y haplogroup and no statement. When
Yleaf did not run, the calibration withholds it, or its QC score is below the minimum, no
haplogroup is stated: the claim is unassessed, its reason why.
"""

import json
from dataclasses import asdict, dataclass

from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity
from genokeeper.yhaplogroup import YLEAF_VERSION, YHaplogroup

YLEAF_TOOL = "yleaf"
TREE = f"Yleaf {YLEAF_VERSION} bundled YFull"
STRONG_QC = 0.99


@dataclass(frozen=True)
class YleafSource:
    """Yleaf's haplogroup (None: not stated), why Yleaf did not run (`not_run`) or the
    calibration withholds it (`withheld`), the QC score below which no haplogroup is
    stated, and whether the sample is male."""

    found: YHaplogroup | None
    not_run: str | None
    withheld: str | None
    min_qc: float
    male: bool
    tool: str = YLEAF_TOOL

    def _reason(self) -> str | None:
        """Why no haplogroup is stated, if none is."""
        if self.not_run is not None:
            return self.not_run
        if self.withheld is not None:
            return self.withheld
        if self.found is None:
            return "Yleaf gave no haplogroup"
        if self.found.qc_score < self.min_qc:
            return (
                f"Yleaf QC score {self.found.qc_score:.2f} below {self.min_qc:g};"
                " no haplogroup stated"
            )
        return None

    def collect(self, ctx: RunContext) -> Collected:
        """The haplogroup's one statement; none for a sample that is not male."""
        if not self.male:
            return Collected(
                tool=self.tool, rows=[], extra={}, input_rows=0, dropped={}
            )
        reason = self._reason()
        found = None if reason is not None else self.found
        if found is None:
            # The validity a stated haplogroup has at its lowest, so the default filters
            # name the gap (as FMR1's not-covered statement is named).
            validity, level = Validity.MODERATE, "Yleaf not stated"
        else:
            validity = (
                Validity.STRONG if found.qc_score >= STRONG_QC else Validity.MODERATE
            )
            level = f"Yleaf QC score {found.qc_score:.2f}"
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "y_haplogroup",
            "haplogroup": None if found is None else found.name,
            "marker": None if found is None else found.marker,
            "qc_score": None if self.found is None else self.found.qc_score,
            "valid_markers": None if self.found is None else self.found.valid_markers,
            "tree": TREE,
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{YLEAF_TOOL}:haplogroup",
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
            "source": YLEAF_TOOL,
            "source_id": TREE,
            "source_version": YLEAF_VERSION,
            "points": None,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement, self.found and asdict(self.found))],
            extra={},
            input_rows=1,
            dropped={},
        )
