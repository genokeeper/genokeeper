"""PharmCAT's gene results as a source (`PgxSource`): per gene a gene-level statement under
drugs and a `gene_call` row (track design, slice 1b, "The store").

A resolved gene is a finding, a normal result included: a normal metabolizer is a result,
not an absence. An unresolved one (`evidence.pgx.unresolved`) is unassessed: the schema
allows a gene-level unassessed statement because it names its gene.

The claim holds the gene's call and, per drug, its name, source and classification; no
implication and no recommendation text. The control, against dosing text in the store: its
anchor comes from `pgx_claim`'s fixed keys, filled from `GeneResult` fields that
`pharmcat.parse_report` takes from named report keys, which a report's recommendation text
cannot reach because neither function reads a free-text field into the claim.

Slice 3f: UGT1A1's gene call says so when UGT1A1*28 was not read (`UGT1A1_28_UNREAD`: the
calibration does not trust rs3064744's low-graded call); IFNL3's result is a non-finding
retired by CPIC (`IFNL3_RETIRED`); and every run holds one unassessed `pharmcat:G6PD-gap`
naming the G6PD haplotypes PharmCAT does not call (`G6PD_GAP`).
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.evidence.pgx import grade_gene_result, source_level, unresolved
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.pgx.outside import OutsideCall
from genokeeper.pgx.pharmcat import UGT1A1_28, GeneResult
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA, WITHHELD, withheld_entry
from genokeeper.store.db import Kind, Topic, Validity

PGX_TOOL = "pharmcat"
UGT1A1_28_UNREAD = (
    "UGT1A1*28 (rs3064744, TA repeat) not read: graded low, not calibrated"
)
UGT1A1_28_MEANS = (
    "UGT1A1*28 is the TA7 allele (ALT CATAT; reference TA6): PharmCAT's diplotype here did"
    " not read this call, so it may differ with it"
)
UGT1A1_28_CONFIRM = "UGT1A1 genotyping in a laboratory"
IFNL3_RETIRED = "retired by CPIC (peginterferon pairs)"
G6PD_GAP = "not covered: PharmCAT does not call the G6PD Mediterranean haplotype or Villeurbanne"


def ugt1a1_read(
    graded: Mapping[tuple[str, int, str, str], tuple[Called, Reliability]],
) -> list[Fields]:
    """UGT1A1*28's calls, the ALT each dosage counts, where it is known."""
    found = ((a, *graded[a]) for a in UGT1A1_28 if a in graded)
    return [
        {"alt": a[3], "gt": c.gt, "dosage": c.dosage, "level": r.level}
        for a, c, r in found
        if c.dosage != "unknown"
    ]


def pgx_claim(result: GeneResult) -> str:
    """The claim JSON of a gene result (schema `genomelab.claim/1`, type
    `pharmcat_gene`)."""
    return json.dumps(
        {
            "schema": CLAIM_SCHEMA,
            "type": "pharmcat_gene",
            "gene": result.gene,
            "diplotype": result.diplotype,
            "phenotype": result.phenotype,
            "activity_score": result.activity_score,
            "call_source": result.call_source,
            "drugs": [
                {
                    "drug": d.drug,
                    "source": d.source,
                    "classification": d.classification,
                }
                for d in result.drugs
            ],
        },
        ensure_ascii=False,
    )


def _score(text: str | None) -> float | None:
    """A plain activity score as a number; a bin (`≥3.0`) or alternatives stay in the
    claim's text only."""
    try:
        return None if text is None else float(text)
    except ValueError:
        return None


@dataclass(frozen=True)
class PgxSource:
    """PharmCAT's gene results (`pharmcat.parse_report`), the outside calls it was given
    (gene -> call), CPIC's levels (`evidence.pgx.load_cpic_levels`) and PharmCAT's version,
    the statements' source version. `ugt1a1_unread`: UGT1A1*28 was written as a no-call
    because the calibration does not trust it; `ugt1a1_calls` (`ugt1a1_read`): those calls,
    kept under the gene result's claim's `withheld`."""

    results: Sequence[GeneResult]
    outside: Mapping[str, OutsideCall]
    cpic_levels: Mapping[tuple[str, str], str]
    version: str
    tool: str = PGX_TOOL
    ugt1a1_unread: bool = False
    ugt1a1_calls: Sequence[Mapping[str, object]] = ()

    def _g6pd_gap(self, ctx: RunContext) -> Row:
        claim = {"schema": CLAIM_SCHEMA, "type": "pharmcat_gap", "gene": "G6PD"}
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{PGX_TOOL}:G6PD-gap",
            "kind": str(Kind.UNASSESSED),
            "vid": None,
            "gene": "G6PD",
            "claim": json.dumps(claim, ensure_ascii=False),
            "validity": str(Validity.STRONG),
            "source_level": "not covered",
            "reason": G6PD_GAP,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.DRUGS),
            "source": PGX_TOOL,
            "source_id": "G6PD-gap",
            "source_version": self.version,
            "points": None,
            "exclusion": False,
        }
        return Row(None, None, statement)

    def collect(self, ctx: RunContext) -> Collected:
        """One gene-level statement and one gene call per gene result, and the G6PD
        gap."""
        rows: list[Row] = []
        gene_calls: list[Fields] = []
        for result in self.results:
            outside = self.outside.get(result.gene)
            validity, actionable, level, warnings = grade_gene_result(
                result, cpic_levels=self.cpic_levels, outside=outside
            )
            kind = (
                Kind.FINDING if unresolved(result, outside) is None else Kind.UNASSESSED
            )
            claim = pgx_claim(result)
            if result.gene == "UGT1A1" and self.ugt1a1_unread:
                warnings = (*warnings, UGT1A1_28_UNREAD)
                if self.ugt1a1_calls:
                    entry = withheld_entry(
                        UGT1A1_28_UNREAD,
                        {"rs3064744": list(self.ugt1a1_calls)},
                        would_be=kind,
                        significance=UGT1A1_28_MEANS,
                        confirm=UGT1A1_28_CONFIRM,
                    )
                    held = {**json.loads(claim), WITHHELD: [entry]}
                    claim = json.dumps(held, ensure_ascii=False)
            reason = "; ".join(warnings) or None
            if result.gene == "IFNL3":
                kind, actionable, reason = Kind.NON_FINDING, False, IFNL3_RETIRED
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"pharmcat:{result.gene}",
                "kind": str(kind),
                "vid": None,
                "gene": result.gene,
                "claim": claim,
                "validity": str(validity),
                "source_level": source_level(result),
                "reason": reason,
                "effect": None,
                "actionable": actionable,
                "topic": str(Topic.DRUGS),
                "source": PGX_TOOL,
                "source_id": result.gene,
                "source_version": self.version,
                "points": None,
            }
            rows.append(Row(None, None, statement))
            callers = (
                dict(outside.callers)
                if outside is not None
                else {"PharmCAT": result.diplotype or "no call"}
            )
            gene_calls.append(
                {
                    "run_id": ctx.run_id,
                    "gene": result.gene,
                    "diplotype": result.diplotype,
                    "phenotype": result.phenotype,
                    "activity_score": _score(result.activity_score),
                    "callers": json.dumps(callers, ensure_ascii=False),
                    "agreed": None if outside is None else outside.agreed,
                    "level": str(level),
                    "warnings": list(warnings),
                }
            )
        rows.append(self._g6pd_gap(ctx))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={"gene_call": gene_calls},
            input_rows=len(self.results) + 1,
            dropped={},
        )
