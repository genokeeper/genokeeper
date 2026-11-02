"""Variants ClinVar has not classified, annotated by VEP and graded as predictions, as a
source (`PredictedSource`; track design, slice 1e).

Each VEP record that has a rule (`evidence.predicted.predicted_rule`) comes with the user's
call at its allele and the call's grade; `collect` grades it (`grade_predicted`) with
gnomAD's frequency at its vid, ClinGen's haploinsufficient genes, where its new stop
escapes NMD on its MANE transcript (UCSC's track) and ClinVar's P/LP truncations of the
release, and writes one statement for each that is a finding, unassessed, or one copy in a
recessive gene (a non-finding the promotion reads), its id the vid and the unversioned
transcript, as a variant can have records on two; the rest are counted under the reason
they are none. Every graded record, and every other VEP input variant off chrM (mutserve's)
counted once as `NO_PREDICTED_EFFECT`, is one of these, so the accounting balances.

Two hets in one recessive gene, Tier 3 carrier findings or one-copy non-findings, are a
candidate pair for read-backed phasing (`promote`), with each other or with a
heterozygous ClinVar claim.
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.annotate.vep import VEP_VERSION, VepRecord, locus
from genokeeper.evidence.clinvar import BA1_REASON, BS1_REASON, clnsig_class
from genokeeper.evidence.predicted import (
    BELOW_THRESHOLDS,
    NOT_CARRIED,
    NOT_HAPLOINSUFFICIENT,
    STOP_CONSEQUENCES,
    TRUNCATING,
    EscapeRegion,
    Site,
    escape,
    grade_predicted,
    truncations_3prime,
)
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim
from genokeeper.sources import Collected, Fields, Row, RunContext, panel

PREDICTED_TOOL = "vep"
NO_PREDICTED_EFFECT = "no_predicted_effect"
# A non-finding's reason -> the drop reason it is counted under.
DROP_REASONS = {
    BELOW_THRESHOLDS: NO_PREDICTED_EFFECT,
    BA1_REASON: "population_frequency",
    BS1_REASON: "population_frequency",
    NOT_HAPLOINSUFFICIENT: "not_haploinsufficient",
    NOT_CARRIED: "not_carried",
}


def clinvar_truncations(claims: Sequence[Claim]) -> list[tuple[str, str, int]]:
    """The release's P/LP claims whose molecular consequence truncates (nonsense or
    frameshift), as (gene, chrom, pos)."""
    return [
        (c.gene, c.chrom, c.pos)
        for c in claims
        if clnsig_class(c.clnsig) == "pathogenic"
        and any(part.split("|")[0] in TRUNCATING for part in c.consequence.split(","))
    ]


@dataclass(frozen=True)
class PredictedSource:
    """The graded candidates (record, call, its grade, homopolymer run) of one run's VEP
    annotation, how many variants VEP was given (each counted once, its records apart), UCSC's escape regions per transcript,
    ClinVar's truncations (`clinvar_truncations`) and ClinGen's haploinsufficient genes."""

    sites: Sequence[tuple[VepRecord, Called, Reliability, int]]
    input_rows: int
    escapes: Mapping[str, Sequence[EscapeRegion]]
    truncations: Sequence[tuple[str, str, int]]
    haploinsufficient: frozenset[str]
    version: str = VEP_VERSION
    tool: str = PREDICTED_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One variant, call and statement per site graded a statement (one recessive
        copy included); the rest, and the input variants without a candidate, counted by
        reason."""
        faf = ctx.faf.lookup([r.vid for r, _, _, _ in self.sites])
        dropped = dict.fromkeys(sorted(set(DROP_REASONS.values())), 0)
        # A variant can have records on two transcripts or genes (`parse_vep`): each graded
        # record is counted once, each variant with none once.
        dropped[NO_PREDICTED_EFFECT] = self.input_rows - len(
            {r.vid for r, _, _, _ in self.sites}
        )
        rows: list[Row] = []
        for record, called, reliability, run in self.sites:
            chrom, pos, _, _ = locus(record.vid)
            transcript = record.transcript.split(".")[0]
            regions = self.escapes.get(transcript, [])
            stops = record.terms & STOP_CONSEQUENCES
            found = escape(record, regions) if stops else None
            site = Site(
                record=record,
                dosage=called.dosage,
                faf=faf.get(record.vid),
                homopolymer=run,
                flags=reliability.flags,
                escape=found,
                truncations_3prime=0
                if found is None
                else truncations_3prime(found, regions, record.gene, self.truncations),
            )
            grade = grade_predicted(
                site,
                haploinsufficient=self.haploinsufficient,
                lists=ctx.lists,
                grading=ctx.grading,
                inheritance=ctx.inheritance.get(record.gene),
                male=ctx.male,
                haploid=ctx.haploid(chrom, pos),
            )
            if grade.reason in DROP_REASONS:
                dropped[DROP_REASONS[grade.reason]] += 1
                continue
            claim = {
                "schema": panel.CLAIM_SCHEMA,
                "type": "predicted_variant",
                "classification": "predicted, not classified",
                "consequence": record.consequence,
                "transcript": record.transcript,
                "hgvsc": record.hgvsc,
                "hgvsp": record.hgvsp,
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"predicted:{record.vid}:{transcript}",
                "kind": str(grade.kind),
                "vid": record.vid,
                "gene": record.gene,
                "claim": json.dumps(claim),
                "validity": None if grade.level is None else str(grade.level),
                "source_level": grade.source_level,
                "reason": grade.reason,
                "effect": None,
                # A prediction is no classification: never actionable (`query`'s order).
                "actionable": False,
                "topic": None if grade.topic is None else str(grade.topic),
                "source": PREDICTED_TOOL,
                "source_id": record.transcript,
                "source_version": self.version,
                "points": None,
            }
            variant: Fields = {"run_id": ctx.run_id, "vid": record.vid, "rsid": None}
            call = panel.call_row(ctx.run_id, record.vid, called, reliability)
            rows.append(Row(variant, call, statement))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=dropped[NO_PREDICTED_EFFECT] + len(self.sites),
            dropped=dropped,
        )
