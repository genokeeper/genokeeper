"""ClinVar's genome-wide GRCh38 claims, each on our own call and graded, as a source
(`ClaimsSource`), and the rows a ClinVar statement and a call are stored as."""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.evidence.clinvar import (
    ValidityGrade,
    clnsig_class,
    grade_clinvar,
    sf_actionable,
)
from genokeeper.evidence.disease import match_curation, parse_clndisdb
from genokeeper.evidence.pgx import cpic_actionable
from genokeeper.genelists import GeneList
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim, ClaimSet
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.store.db import Kind, Level

CLAIMS_TOOL = "clinvar"
# Written into every stored claim since before the rename to genokeeper; a new name would be
# a new schema, read beside this one.
CLAIM_SCHEMA = "genomelab.claim/1"
# The claim key of what a calibration table withheld (`withheld_entry`, `sources.withheld`).
WITHHELD = "withheld"


def withheld_entry(
    reason: str,
    value: Fields,
    *,
    would_be: Kind,
    significance: str | None,
    confirm: str | None,
) -> Fields:
    """One element of a claim's `WITHHELD` list: what the tool measured (`value`), at the
    lowest reliability, why it is withheld, and were it true the kind and reason the
    statement would have (`would_be`, `significance`) and the test that confirms it."""
    return {
        "reliability": str(Level.LOW),
        "reason": reason,
        "value": value,
        "would_be": str(would_be),
        "significance": significance,
        "confirm": confirm,
    }


def clinvar_claim(
    *, variation_id: int, clnsig: str, clnrevstat: str, condition: str
) -> str:
    """The claim JSON of a ClinVar germline statement (schema `genomelab.claim/1`)."""
    return json.dumps(
        {
            "schema": CLAIM_SCHEMA,
            "type": "clinvar_germline",
            "classification": clnsig,
            "review_status": clnrevstat,
            "condition": condition,
            "variation_id": variation_id,
        }
    )


def clinvar_statement(
    run_id: str,
    *,
    variation_id: int,
    vid: str,
    gene: str,
    clnsig: str,
    clnrevstat: str,
    condition: str,
    consequence: str,
    grade: ValidityGrade,
    source_version: str,
    lists: Mapping[str, GeneList],
    cpic_levels: Mapping[tuple[str, str], str] | None = None,
) -> Fields:
    """One ClinVar statement row, `clinvar:<variation_id>`; an unmapped CLNSIG's alert
    joins the reason. `actionable` is a finding ACMG SF v3.3 would report
    (`evidence.clinvar.sf_actionable`, from the grade's kind and topic and the claim's
    molecular `consequence`), or a drug_response claim whose gene a drug its condition names
    pairs with at CPIC level A or B (`cpic_levels`, `evidence.pgx.cpic_actionable`)."""
    reason = "; ".join(r for r in (grade.reason, grade.alert) if r) or None
    actionable = sf_actionable(
        clnsig=clnsig,
        gene=gene,
        lists=dict(lists),
        kind=grade.kind,
        topic=grade.topic,
        vid=vid,
        consequence=consequence,
    ) or (
        clnsig_class(clnsig) == "drug_response"
        and cpic_actionable(gene, condition, cpic_levels or {})
    )
    return {
        "run_id": run_id,
        "statement_id": f"clinvar:{variation_id}",
        "kind": str(grade.kind),
        "vid": vid,
        "gene": gene,
        "claim": clinvar_claim(
            variation_id=variation_id,
            clnsig=clnsig,
            clnrevstat=clnrevstat,
            condition=condition,
        ),
        "validity": None if grade.level is None else str(grade.level),
        "source_level": grade.source_level,
        "reason": reason,
        "effect": None,
        "actionable": actionable,
        "topic": None if grade.topic is None else str(grade.topic),
        "source": "clinvar",
        "source_id": str(variation_id),
        "source_version": source_version,
        "points": None,
    }


def call_row(run_id: str, vid: str, called: Called, reliability: Reliability) -> Fields:
    """The user's call at `vid` with its grade, as a row; the store keeps the first call
    written at a vid (`write_collected`)."""
    return {
        "run_id": run_id,
        "vid": vid,
        "gt": called.gt,
        "dosage": called.dosage,
        "p_wrong": reliability.p_wrong,
        "p_wrong_lo": None
        if reliability.p_wrong_ci is None
        else reliability.p_wrong_ci[0],
        "p_wrong_hi": None
        if reliability.p_wrong_ci is None
        else reliability.p_wrong_ci[1],
        "level": str(reliability.level),
        "gq": called.gq,
        "dp": called.dp,
        "ad": None if called.ad is None else list(called.ad),
        "vaf": called.vaf,
        "row_type": str(called.row_type),
        "flags": list(reliability.flags),
        "warnings": list(reliability.warnings),
    }


def grade_claim(
    claim: Claim, called: Called, ctx: RunContext, faf_grpmax: float | None
) -> ValidityGrade:
    """One claim graded on the user's call: on the ClinGen curation of its own disease,
    else the gene-level `ctx.gene_level` (`match_curation`), with gnomAD's filtering allele
    frequency at its vid; `ctx.haploid` takes the claim's GRCh38 locus. A validity GenCC
    gave is named so in the reasons (`ctx.validity_source`); a drug_response claim names
    the ClinPGx annotation at its own allele (`ctx.clinpgx`) or its gene's call
    (`ctx.gene_calls`)."""
    curation = match_curation(
        claim.gene,
        parse_clndisdb(claim.disease_ids),
        ctx.curations,
        ctx.gene_level.get(claim.gene),
    )
    return grade_clinvar(
        clnsig=claim.clnsig,
        clnrevstat=claim.clnrevstat,
        gene=claim.gene,
        gene_validity=curation.classification,
        validity_matched=curation.matched,
        validity_source=ctx.validity_source.get(claim.gene),
        lists=dict(ctx.lists),
        faf_grpmax=faf_grpmax,
        vid=claim.vid,
        ba1_exceptions=ctx.ba1_exceptions,
        dosage=called.dosage,
        grading=ctx.grading,
        inheritance=ctx.inheritance.get(claim.gene),
        male=ctx.male,
        haploid=ctx.haploid(claim.chrom, claim.pos),
        clinpgx=ctx.clinpgx,
        gene_calls=ctx.gene_calls,
    )


@dataclass(frozen=True)
class ClaimsSource:
    """The genome-wide claims of one ClinVar GRCh38 release (`claims`, for its accounting),
    each with the user's call from our own gVCF (`genotype.gvcf.calls_at`) or, on chrM,
    from mutserve (`genokeeper.mito`), and its grade (`graded`, one per claim)."""

    claims: ClaimSet
    graded: Sequence[tuple[Claim, Called, Reliability]]
    tool: str = CLAIMS_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """Each claim's variant, call and ClinVar statement, graded by `grade_claim` (one
        `ctx.faf` lookup over the claims' vids; no frequency when the table lacks the
        vid)."""
        faf = ctx.faf.lookup([claim.vid for claim, _, _ in self.graded])
        rows: list[Row] = []
        for claim, called, reliability in self.graded:
            grade = grade_claim(claim, called, ctx, faf.get(claim.vid))
            variant: Fields = {
                "run_id": ctx.run_id,
                "vid": claim.vid,
                "rsid": claim.rsid,
            }
            statement = clinvar_statement(
                ctx.run_id,
                variation_id=claim.variation_id,
                vid=claim.vid,
                gene=claim.gene,
                clnsig=claim.clnsig,
                clnrevstat=claim.clnrevstat,
                condition=claim.condition,
                consequence=claim.consequence,
                grade=grade,
                source_version=ctx.source_version,
                lists=ctx.lists,
                cpic_levels=ctx.cpic_levels,
            )
            rows.append(
                Row(
                    variant,
                    call_row(ctx.run_id, claim.vid, called, reliability),
                    statement,
                )
            )
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=self.claims.input_rows,
            dropped=self.claims.dropped,
        )
