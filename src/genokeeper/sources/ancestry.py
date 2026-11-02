"""Ancestry as a source (`AncestrySource`; slice 3f, "Admixture"). Tool `ancestry`, three
statements under ancestry, each a finding of information at validity limited, never
actionable, never flagged `exclusion`, its source level the tool and its version:

- `ancestry:similarity`: pgsc_calc's random-forest similarity to the HGDP+1kGP
  superpopulations (`pgsc.parse_similarity`), `SIMILARITY`;
- `ancestry:admixture`: fastmixture's projection on the panel's superpopulations
  (`ancestry.admixture.run_projection`), `ADMIXTURE`;
- `ancestry:neanderthal`: hmmix's archaic share (`ancestry.archaic.archaic_share`),
  `ARCHAIC` (R6).

Each carries `INFORMATION_NOTE`: none of them is checked against a reference genome here.
A result missing, a step in `not_run` or in `withheld` (keyed `similarity`, `admixture`,
`neanderthal`) leaves its statement unassessed with the reason.
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass

from genokeeper.ancestry.archaic import ARCHAIC_GENOME_BP, ArchaicShare
from genokeeper.prs.pgsc import PGSC_CALC_VERSION, PopSimilarity
from genokeeper.sources import INFORMATION_NOTE, Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity
from genokeeper.typing_files import HMMIX_VERSION

ANCESTRY_TOOL = "ancestry"
FASTMIXTURE_VERSION = "1.3.0"
SIMILARITY = "similarity to the HGDP+1kGP populations, not ancestry proportions"
ADMIXTURE = (
    "proportions against the panel's superpopulations; country-level ancestry needs"
    " other panels"
)
ARCHAIC = "every non-African genome carries some; about 1–2% by this method"
LOW_CONFIDENCE = "closest population assigned with low confidence"
NOT_RUN = "not run"
_TOOLS = {
    "similarity": f"pgsc_calc {PGSC_CALC_VERSION}",
    "admixture": f"fastmixture {FASTMIXTURE_VERSION}",
    "neanderthal": f"hmmix {HMMIX_VERSION}",
}


@dataclass(frozen=True)
class AncestrySource:
    """The run's similarity, admixture proportions and archaic share (each None when it
    has none), why a step was not run, and the calibration's withholding reasons."""

    similarity: PopSimilarity | None
    admixture: Mapping[str, float] | None
    archaic: ArchaicShare | None
    not_run: Mapping[str, str]
    withheld: Mapping[str, str]
    tool: str = ANCESTRY_TOOL

    def _results(self) -> dict[str, tuple[dict[str, object] | None, list[str]]]:
        """Per statement: its claim's own fields (None without a result) and its
        reasons."""
        sim, admix, arch = self.similarity, self.admixture, self.archaic
        return {
            "similarity": (
                None
                if sim is None
                else {
                    "population": sim.population,
                    "low_confidence": sim.low_confidence,
                    "probabilities": sim.probabilities,
                },
                [SIMILARITY]
                + ([LOW_CONFIDENCE] if sim is not None and sim.low_confidence else []),
            ),
            "admixture": (
                None if admix is None else {"proportions": dict(admix)},
                [ADMIXTURE],
            ),
            "neanderthal": (
                None
                if arch is None
                else {
                    "neanderthal_bp": arch.neanderthal_bp,
                    "denisovan_bp": arch.denisovan_bp,
                    "share": arch.share,
                    "genome_bp": ARCHAIC_GENOME_BP,
                },
                [ARCHAIC],
            ),
        }

    def collect(self, ctx: RunContext) -> Collected:
        """The three statements."""
        rows: list[Row] = []
        for name, (fields, reasons) in self._results().items():
            kind, reason = Kind.FINDING, "; ".join([*reasons, INFORMATION_NOTE])
            if name in self.not_run:
                kind, reason = Kind.UNASSESSED, self.not_run[name]
            elif name in self.withheld:
                kind, reason = Kind.UNASSESSED, self.withheld[name]
            elif fields is None:
                kind, reason = Kind.UNASSESSED, NOT_RUN
            claim = {
                "schema": CLAIM_SCHEMA,
                "type": f"ancestry_{name}",
                **(fields if kind is Kind.FINDING and fields is not None else {}),
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"{ANCESTRY_TOOL}:{name}",
                "kind": str(kind),
                "vid": None,
                "gene": None,
                "claim": json.dumps(claim, ensure_ascii=False),
                "validity": str(Validity.LIMITED),
                "source_level": _TOOLS[name],
                "reason": reason,
                "effect": None,
                "actionable": False,
                "topic": str(Topic.ANCESTRY),
                "source": ANCESTRY_TOOL,
                "source_id": name,
                "source_version": _TOOLS[name],
                "points": None,
                "exclusion": False,
            }
            rows.append(Row(None, None, statement))
        return Collected(
            tool=self.tool, rows=rows, extra={}, input_rows=len(rows), dropped={}
        )
