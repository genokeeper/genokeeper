"""Polygenic scores and APOE as a source (`PrsSource`; track design, slice 1d, "Findings").

Each chosen score is one statement with neither vid nor gene (`pgsc_calc:<PGS ID>`): a
finding with its percentile, decile, Z and `Z_norm2` against the population pgsc_calc finds
the user closest to, in `effect.score` beside a null relative and absolute effect (a Z is
no odds ratio, and the order would rank it as one). The decile is what the user is shown
(slice 3f, "PGS"; `DECILE_REASON` in a scored statement's reasons); the percentile stays
for the change feed. Its validity is moderate with a PGS Catalog evaluation in
that population and limited otherwise, or when pgsc_calc assigned the population with low
confidence. A score our target matches below `MIN_OVERLAP` of its variants is unassessed,
its reason the rate, and what pgsc_calc computed is kept in the claim's `WITHHELD` list
for the agent (`sources.withheld`); one pgsc_calc gave no percentile at a passing rate is
unassessed for that reason. A weak score (slice 3g) adds `WEAK_REASON` last, scored or not,
and each condition no score qualifies for is one unassessed `pgsc_calc:not-covered:<slug>`.

A research score (`ResearchScoreSource`, a research layer: `sources.research`) is one
`pgsc_calc:research:<PGS ID>` statement, limited and under the research topic, whatever its
overlap (given first in its reason): nothing is hidden from the agent, and its reliability
is stated.

APOE is one gene-level statement (`apoe:APOE`, gene APOE, under risks) from our calls at
the ClinVar claims of rs429358 C and rs7412 T, with a gene call (the diplotype, at the
weaker call's level). Its validity is the weaker of the two SNPs' ClinVar grades, each
SNP's strongest claim first (builder's call). An unknown call, or a double heterozygote no
read phased, leaves it unassessed with its reason; a release without a claim at either SNP
drops it, counted as `apoe_no_clinvar_claim`.
"""

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim
from genokeeper.phasing import Phase
from genokeeper.prs.pgsc import (
    APOE_E2,
    APOE_E4,
    MIN_OVERLAP,
    PGSC_CALC_VERSION,
    ScoreResult,
    apoe_genotype,
)
from genokeeper.prs.scores import NotCovered, ResearchScore, Score
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import (
    CLAIM_SCHEMA,
    WITHHELD,
    grade_claim,
    withheld_entry,
)
from genokeeper.store.db import Kind, Level, Topic, Validity

PRS_TOOL = "pgsc_calc"
PRS_RESEARCH_TOOL = "pgsc_calc_research"
RESEARCH_SOURCE_LEVEL = "research score: did not pass the selection rule"
RESEARCH_REASON = (
    "this score says little about one person: it did not pass the selection rule and is"
    " shown for research only"
)
WITHHELD_SIGNIFICANCE = (
    "Too few of the score's variants matched for its decile to be reliable."
)
DECILE_REASON = (
    "shown in deciles: one genome's score is less precise than its percentile"
    " (Ding et al., Nat Genet 2022)"
)
APOE_DROPPED = "apoe_no_clinvar_claim"
_SNPS = {APOE_E4: "rs429358", APOE_E2: "rs7412"}
_STRENGTH = {
    Validity.CONFLICTING: 0,
    Validity.LIMITED: 1,
    Validity.MODERATE: 2,
    Validity.STRONG: 3,
}
_LEVEL = {Level.LOW: 0, Level.MEDIUM: 1, Level.HIGH: 2}
WEAK_REASON = (
    "weak predictive power: the score explains little of the differences between people,"
    " so the percentile says little about this one person"
)


def decile(percentile: float) -> int:
    """The decile (1..10) a percentile (0..100) falls in; 100 is the tenth."""
    return min(10, math.floor(percentile / 10) + 1)


def prs_validity(score: Score, result: ScoreResult) -> tuple[Validity, str | None]:
    """Moderate when the score has a PGS Catalog evaluation in the population pgsc_calc
    finds the user closest to, assigned with confidence; limited, with why, otherwise."""
    if result.population not in score.evaluation_ancestry:
        return Validity.LIMITED, f"no PGS Catalog evaluation in {result.population}"
    if result.low_confidence:
        return Validity.LIMITED, (
            f"closest population {result.population} assigned with low confidence"
        )
    return Validity.MODERATE, None


def _score_effect(result: ScoreResult, scored: bool) -> dict[str, object]:
    """The statement's effect: a null relative and absolute effect, and the score, whose
    values are null unless `scored`."""
    return {
        "relative": None,
        "absolute": None,
        "score": {
            "percentile": result.percentile if scored else None,
            "z": result.z if scored else None,
            "population": result.population,
            "match_rate": result.match_rate,
            "z_norm2": result.z_norm2 if scored else None,
            "decile": decile(result.percentile)
            if scored and result.percentile is not None
            else None,
        },
    }


def _score_row(ctx: RunContext, score: Score, result: ScoreResult) -> Row:
    validity, why = prs_validity(score, result)
    unscored = (
        f"score overlap {result.match_rate:.2f}, below {MIN_OVERLAP}"
        if result.match_rate < MIN_OVERLAP
        else "pgsc_calc gave no percentile"
        if result.percentile is None
        else None
    )
    scored = unscored is None
    reasons = [
        r
        for r in (
            unscored,
            why,
            None if unscored else DECILE_REASON,
            WEAK_REASON if score.weak else None,
        )
        if r
    ]
    claim: Fields = {
        "schema": CLAIM_SCHEMA,
        "type": "pgs",
        "pgs_id": score.pgs_id,
        "condition": score.condition,
        "trait_efo": score.trait_efo,
        "publication": score.publication,
        "development_ancestry": score.development_ancestry,
        "evaluations": [list(e) for e in score.evaluations],
        "licence": score.licence,
        "weak": score.weak,
    }
    if unscored is not None and result.percentile is not None:
        claim[WITHHELD] = [
            withheld_entry(
                unscored,
                {
                    "decile": decile(result.percentile),
                    "percentile": result.percentile,
                    "match_rate": result.match_rate,
                    "population": result.population,
                },
                would_be=Kind.FINDING,
                significance=WITHHELD_SIGNIFICANCE,
                confirm=None,
            )
        ]
    effect = _score_effect(result, scored)
    statement: Fields = {
        "run_id": ctx.run_id,
        "statement_id": f"{PRS_TOOL}:{score.pgs_id}",
        "kind": str(Kind.FINDING if scored else Kind.UNASSESSED),
        "vid": None,
        "gene": None,
        "claim": json.dumps(claim, ensure_ascii=False),
        "validity": str(validity),
        "source_level": "PGS Catalog evaluations in "
        + ", ".join(sorted(score.evaluation_ancestry)),
        "reason": "; ".join(reasons) or None,
        "effect": json.dumps(effect),
        "actionable": False,
        "topic": str(score.topic),
        "source": PRS_TOOL,
        "source_id": score.pgs_id,
        "source_version": PGSC_CALC_VERSION,
        "points": None,
    }
    return Row(None, None, statement)


def _research_row(ctx: RunContext, score: ResearchScore, result: ScoreResult) -> Row:
    scored = result.percentile is not None
    overlap = (
        f"score overlap {result.match_rate:.2f}, below {MIN_OVERLAP}"
        if result.match_rate < MIN_OVERLAP
        else None
    )
    reasons = [
        r
        for r in (
            overlap,
            None if scored else "pgsc_calc gave no percentile",
            score.why,
            DECILE_REASON if scored else None,
            RESEARCH_REASON,
        )
        if r
    ]
    claim = {
        "schema": CLAIM_SCHEMA,
        "type": "pgs_research",
        "pgs_id": score.pgs_id,
        "condition": score.condition,
        "trait_efo": score.trait_efo,
        "publication": score.publication,
        "development_ancestry": score.development_ancestry,
        "evaluations": [list(e) for e in score.evaluations],
        "licence": score.licence,
        "why": score.why,
    }
    statement: Fields = {
        "run_id": ctx.run_id,
        "statement_id": f"{PRS_TOOL}:research:{score.pgs_id}",
        "kind": str(Kind.FINDING if scored else Kind.UNASSESSED),
        "vid": None,
        "gene": None,
        "claim": json.dumps(claim, ensure_ascii=False),
        "validity": str(Validity.LIMITED),
        "source_level": RESEARCH_SOURCE_LEVEL,
        "reason": "; ".join(reasons),
        "effect": json.dumps(_score_effect(result, scored)),
        "actionable": False,
        "topic": str(Topic.RESEARCH),
        "source": PRS_RESEARCH_TOOL,
        "source_id": score.pgs_id,
        "source_version": PGSC_CALC_VERSION,
        "points": None,
    }
    return Row(None, None, statement)


def _not_covered_row(ctx: RunContext, entry: NotCovered) -> Row:
    claim = {
        "schema": CLAIM_SCHEMA,
        "type": "pgs_not_covered",
        "condition": entry.condition,
    }
    statement: Fields = {
        "run_id": ctx.run_id,
        "statement_id": f"{PRS_TOOL}:not-covered:{entry.slug}",
        "kind": str(Kind.UNASSESSED),
        "vid": None,
        "gene": None,
        "claim": json.dumps(claim, ensure_ascii=False),
        "validity": str(Validity.LIMITED),
        "source_level": "no qualifying PGS Catalog score",
        "reason": f"not covered: no score meets the criteria ({entry.why})",
        "effect": None,
        "actionable": False,
        "topic": str(entry.topic),
        "source": PRS_TOOL,
        "source_id": entry.slug,
        "source_version": PGSC_CALC_VERSION,
        "points": None,
    }
    return Row(None, None, statement)


@dataclass(frozen=True)
class ResearchScoreSource:
    """pgsc_calc's result for each research score (`pgsc.parse_results`, pgs_id -> result,
    one for every score)."""

    scores: Sequence[ResearchScore]
    results: Mapping[str, ScoreResult]
    tool: str = PRS_RESEARCH_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One statement per research score."""
        rows = [_research_row(ctx, s, self.results[s.pgs_id]) for s in self.scores]
        return Collected(
            tool=self.tool, rows=rows, extra={}, input_rows=len(rows), dropped={}
        )


@dataclass(frozen=True)
class PrsSource:
    """pgsc_calc's result for each chosen score (`pgsc.parse_results`, pgs_id -> result,
    one for every score), and APOE: the graded ClinVar claims at its two SNPs (`apoe`, as
    the run grades every claim) and the read-backed phase of a double heterozygote; and
    the conditions no score qualifies for (`not_covered`)."""

    scores: Sequence[Score]
    results: Mapping[str, ScoreResult]
    apoe: Sequence[tuple[Claim, Called, Reliability]]
    apoe_phase: Phase | None = None
    not_covered: Sequence[NotCovered] = ()
    tool: str = PRS_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One statement per score and per not-covered condition, and APOE's statement
        and gene call."""
        rows = [_score_row(ctx, s, self.results[s.pgs_id]) for s in self.scores]
        rows += [_not_covered_row(ctx, n) for n in self.not_covered]
        apoe, gene_calls = self._apoe(ctx)
        return Collected(
            tool=self.tool,
            rows=rows + apoe,
            extra={"gene_call": gene_calls},
            input_rows=len(self.scores) + 1 + len(self.not_covered),
            dropped={APOE_DROPPED: 1 - len(apoe)},
        )

    def _apoe(self, ctx: RunContext) -> tuple[list[Row], list[Fields]]:
        at: dict[Allele, list[tuple[Claim, Called, Reliability]]] = {
            snp: [
                g
                for g in self.apoe
                if (g[0].chrom, g[0].pos, g[0].ref, g[0].alt) == snp
            ]
            for snp in _SNPS
        }
        if not all(at.values()):
            return [], []
        faf = ctx.faf.lookup([claim.vid for claim, _, _ in self.apoe])
        strongest = {
            snp: max(
                (
                    (grade_claim(claim, called, ctx, faf.get(claim.vid)), claim)
                    for claim, called, _ in graded
                ),
                key=lambda g: -1 if g[0].level is None else _STRENGTH[g[0].level],
            )
            for snp, graded in at.items()
        }
        validity = min(
            (g.level or Validity.LIMITED for g, _ in strongest.values()),
            key=_STRENGTH.__getitem__,
        )
        calls = {snp: graded[0][1] for snp, graded in at.items()}
        level = min(
            (graded[0][2].level for graded in at.values()), key=_LEVEL.__getitem__
        )
        diplotype = apoe_genotype(
            {snp: c.dosage for snp, c in calls.items()}, self.apoe_phase
        )
        unknown = [_SNPS[snp] for snp, c in calls.items() if c.dosage == "unknown"]
        reason = None
        if diplotype is None:
            reason = (
                f"genotype unknown at {' and '.join(unknown)}"
                if unknown
                else "ε2/ε4 or ε1/ε3: phase unknown"
            )
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": "apoe:APOE",
            "kind": str(Kind.FINDING if diplotype else Kind.UNASSESSED),
            "vid": None,
            "gene": "APOE",
            "claim": json.dumps(
                {
                    "schema": CLAIM_SCHEMA,
                    "type": "apoe_diplotype",
                    "gene": "APOE",
                    "diplotype": diplotype,
                },
                ensure_ascii=False,
            ),
            "validity": str(validity),
            "source_level": "; ".join(
                f"{_SNPS[snp]} {g.source_level}" for snp, (g, _) in strongest.items()
            ),
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.RISKS),
            "source": "clinvar",
            "source_id": ";".join(str(c.variation_id) for _, c in strongest.values()),
            "source_version": ctx.source_version,
            "points": None,
        }
        gene_calls: list[Fields] = []
        if diplotype is not None:
            gene_calls.append(
                {
                    "run_id": ctx.run_id,
                    "gene": "APOE",
                    "diplotype": diplotype,
                    "phenotype": None,
                    "activity_score": None,
                    # The caller's name from before the rename, as the stored runs have it.
                    "callers": json.dumps({"genomelab": diplotype}, ensure_ascii=False),
                    "agreed": None,
                    "level": str(level),
                    "warnings": [],
                }
            )
        return [Row(None, None, statement)], gene_calls
