"""PharmGKB's clinical annotations as a source (`PharmGkbSource`): each annotation located
at its allele (`evidence.pharmgkb.locate`), graded on the user's call there, under drugs
(track design, slice 1b, "PharmGKB clinical annotations").

Carried (one or two copies of the allele) is a finding, not carried a non-finding, an
unknown genotype an unassessed claim; validity is strong for level 1B and moderate for 2A
and 2B (ClinPGx's own words; A3s); `actionable` when one of the annotation's drugs is a CPIC level A or B pair with its
gene (`evidence.pgx.cpic_actionable`).

Built with `tool=PHARMGKB_RESEARCH_TOOL` over the level-3 annotations
(`ClinicalAnnotations.research`; slice 3f, "Research layers, one predicate"), it writes
`pharmgkb_research:<id>` statements under research: limited, never actionable, with the
reason `RESEARCH_REASON`; the query, the change feed and the dashboard hide them as they hide
GWAS associations (`sources.research.is_research_layer`).

The claim holds the annotation's id, gene, categories, phenotypes (as its `condition`) and
drugs; no annotation text. The control, against dosing text in the store: its anchor comes
from `pharmgkb_claim`'s fixed keys, filled from `Annotation` fields that
`evidence.pharmgkb.read_clinical_annotations` takes from named columns, which an
annotation's text cannot reach because no free-text column is read and every kept item
that speaks of a dose is left out.
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.evidence.pgx import cpic_actionable
from genokeeper.evidence.pharmgkb import Annotation, validity
from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.sources import Collected, Fields, Row, RunContext, panel
from genokeeper.store.db import Kind, Topic
from genokeeper.variants.key import make_vid

PHARMGKB_TOOL = "pharmgkb"
PHARMGKB_RESEARCH_TOOL = "pharmgkb_research"
RESEARCH_REASON = "one study or an unconfirmed association"


def pharmgkb_claim(annotation: Annotation) -> str:
    """The claim JSON of an annotation (schema `genomelab.claim/1`, type
    `pharmgkb_clinical_annotation`)."""
    return json.dumps(
        {
            "schema": panel.CLAIM_SCHEMA,
            "type": "pharmgkb_clinical_annotation",
            "annotation_id": annotation.annotation_id,
            "gene": annotation.gene,
            "category": "; ".join(annotation.categories) or None,
            "condition": "; ".join(annotation.phenotypes) or None,
            "drugs": [{"drug": drug} for drug in annotation.drugs],
        },
        ensure_ascii=False,
    )


def _kind(dosage: str) -> Kind:
    if dosage == "unknown":
        return Kind.UNASSESSED
    return Kind.NON_FINDING if dosage == "0" else Kind.FINDING


@dataclass(frozen=True)
class PharmGkbSource:
    """The located annotations, each with its allele (normalised GRCh38), the user's call
    there and its grade; the download's rows and where the rest went; and the download's
    CREATED date, the statements' source version."""

    located: Sequence[tuple[Annotation, Allele, Called, Reliability]]
    input_rows: int
    dropped: Mapping[str, int]
    version: str
    tool: str = PHARMGKB_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One variant, call and statement per located annotation."""
        from genokeeper.sources.research import is_research_layer  # it imports this

        rows: list[Row] = []
        research = is_research_layer(self.tool)
        for annotation, allele, called, reliability in self.located:
            vid = make_vid(*allele)
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"{self.tool}:{annotation.annotation_id}",
                "kind": str(_kind(called.dosage)),
                "vid": vid,
                "gene": annotation.gene,
                "claim": pharmgkb_claim(annotation),
                "validity": str(validity(annotation.level)),
                "source_level": f"ClinPGx level {annotation.level}",
                "reason": RESEARCH_REASON if research else None,
                "effect": None,
                "actionable": not research
                and cpic_actionable(
                    annotation.gene, "; ".join(annotation.drugs), ctx.cpic_levels
                ),
                "topic": str(Topic.RESEARCH if research else Topic.DRUGS),
                "source": self.tool,
                "source_id": annotation.annotation_id,
                "source_version": self.version,
                "points": None,
                "exclusion": False,
            }
            variant: Fields = {
                "run_id": ctx.run_id,
                "vid": vid,
                "rsid": annotation.variant,
            }
            call = panel.call_row(ctx.run_id, vid, called, reliability)
            rows.append(Row(variant, call, statement))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=self.input_rows,
            dropped=dict(self.dropped),
        )
