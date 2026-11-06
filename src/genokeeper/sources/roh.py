"""Runs of homozygosity as a source (`RohSource`; slice 3f, "ROH"). Tool `roh`, one
statement `roh:summary` under ancestry, a finding of information, never actionable, never
flagged `exclusion`.

The reason states the total length in runs of `[roh] min_run` or more, their count and the
longest run, and names each run of `isodisomy_min` or more as one that may suggest isodisomy;
heterodisomy leaves no run. Nothing speaks of the parents (R3 section 4). `not_run` (no
table, no population from pgsc_calc, a failed step) and `withheld` (the calibration's reason)
make it unassessed.
"""

import json
from dataclasses import asdict, dataclass

from genokeeper.roh import RohSummary, Run
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

ROH_TOOL = "roh"
ISODISOMY = (
    "a single run of {length} Mb on {chrom} may suggest isodisomy;"
    " heterodisomy cannot be seen"
)
SOURCE_VERSION = "bcftools roh 1.25"


def _mb(bases: int) -> str:
    return f"{bases / 1_000_000:.1f}"


def _run(run: Run | None) -> dict[str, object] | None:
    if run is None:
        return None
    return {
        "chrom": run.chrom,
        "start": run.start,
        "end": run.end,
        "length": run.length,
    }


@dataclass(frozen=True)
class RohSource:
    """The run's summary (`roh.summarise`; None when it has none), why ROH was not run,
    and the calibration table's withholding reason."""

    summary: RohSummary | None
    not_run: str | None
    withheld: str | None
    tool: str = ROH_TOOL

    def _reason(self, ctx: RunContext) -> tuple[Kind, str]:
        if self.not_run is not None:
            return Kind.UNASSESSED, self.not_run
        if self.withheld is not None:
            return Kind.UNASSESSED, self.withheld
        s = self.summary
        if s is None:
            return Kind.UNASSESSED, "not run"
        floor = _mb(ctx.grading.roh.min_run).removesuffix(".0")
        parts = [f"{_mb(s.total_long)} Mb in {s.n_long} runs of {floor} Mb or more"]
        if s.longest is not None:
            parts.append(f"the longest {_mb(s.longest.length)} Mb on {s.longest.chrom}")
        parts += [
            ISODISOMY.format(length=_mb(r.length), chrom=r.chrom) for r in s.isodisomy
        ]
        return Kind.FINDING, "; ".join(parts)

    def collect(self, ctx: RunContext) -> Collected:
        """The one `roh:summary` statement."""
        kind, reason = self._reason(ctx)
        s = self.summary
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "roh_summary",
            "min_run": ctx.grading.roh.min_run,
            "total_long": None if s is None else s.total_long,
            "n_long": None if s is None else s.n_long,
            "longest": None if s is None else _run(s.longest),
            "isodisomy": [] if s is None else [_run(r) for r in s.isodisomy],
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{ROH_TOOL}:summary",
            "kind": str(kind),
            "vid": None,
            "gene": None,
            "claim": json.dumps(claim, ensure_ascii=False),
            "validity": str(Validity.MODERATE),
            "source_level": "bcftools roh on our genotypes",
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.ANCESTRY),
            "source": ROH_TOOL,
            "source_id": "summary",
            "source_version": SOURCE_VERSION,
            "points": None,
            "exclusion": False,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement, s and asdict(s))],
            extra={},
            input_rows=1,
            dropped={},
        )
