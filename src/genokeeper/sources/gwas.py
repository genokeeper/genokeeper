"""Single-SNP research associations as a source (`GwasSource`; track design, slice 1d,
"Research associations"; stage 1 design, "Clinical validity").

Each located association (`evidence.gwas`) is one statement at its panel allele's vid,
graded on the user's call there by its risk allele's copies (`graded`, from
`genotype.sites.risk_copies`): one or two a finding, none a non-finding, unknown (only on
an unknown genotype) unassessed; not countable on a known call (another allele beside a
reference risk allele; ruling R26) no statement, counted as `other_allele_at_site`. The
claim names the risk allele and its copies, so an answer can say "two copies of risk
allele C" whichever allele the vid counts; `check_calls` recounts both from the stored
call with the same function.

Validity is limited for every one, whatever its p-value: the GWAS Catalog has no
replication field, and a single SNP's association is research-grade (stage 1 design). The
effect is the catalog's odds ratio or beta, never an absolute risk. The topic is that of the
first of its traits the EFO table maps (`prs.scores.topic_for`) to a topic other than
research, else research.
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.evidence.gwas import Association
from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.genotype.sites import NOT_COUNTABLE, risk_copies
from genokeeper.prs.scores import topic_for
from genokeeper.sources import Collected, Fields, Row, RunContext, panel
from genokeeper.store.db import Kind, Topic, Validity
from genokeeper.variants.key import make_vid

GWAS_TOOL = "gwas_catalog"
OTHER_ALLELE = "other_allele_at_site"  # the drop reason of an uncountable association


def graded(
    vid: str, risk_allele: str, gt: str | None, dosage: str | None
) -> tuple[int | None, Kind] | None:
    """An association's risk-allele copies on the call at `vid` (`risk_copies`) and the
    kind they give; None when they are not countable."""
    copies = risk_copies(vid, risk_allele, gt, dosage)
    if copies == NOT_COUNTABLE:
        return None
    if copies is None:
        return None, Kind.UNASSESSED
    return copies, Kind.FINDING if copies else Kind.NON_FINDING


@dataclass(frozen=True)
class GwasSource:
    """The located associations, each with its allele (the panel's, normalised GRCh38),
    the user's call there and its grade; the download's rows and where the rest went;
    and its release, the statements' source version."""

    located: Sequence[tuple[Association, Allele, Called, Reliability]]
    input_rows: int
    dropped: Mapping[str, int]
    version: str
    tool: str = GWAS_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One variant, call and statement per located association whose copies are
        countable; the rest are dropped as `OTHER_ALLELE`."""
        topics: dict[str, Topic] = {}
        rows: list[Row] = []
        dropped = dict(self.dropped)
        for a, allele, called, reliability in self.located:
            vid = make_vid(*allele)
            grade = graded(vid, a.risk_allele, called.gt, called.dosage)
            if grade is None:
                dropped[OTHER_ALLELE] = dropped.get(OTHER_ALLELE, 0) + 1
                continue
            copies, kind = grade
            for efo in a.efo_ids:
                topics.setdefault(efo, topic_for(efo))
            topic = next(
                (topics[e] for e in a.efo_ids if topics[e] is not Topic.RESEARCH),
                Topic.RESEARCH,
            )
            relative = (
                None
                if a.measure is None
                else {
                    "measure": a.measure,
                    "value": a.value,
                    "ci": None if a.ci is None else list(a.ci),
                }
            )
            claim = {
                "schema": panel.CLAIM_SCHEMA,
                "type": "gwas_association",
                "condition": a.trait,
                "trait_efo": ", ".join(a.efo_ids),
                "risk_allele": a.risk_allele,
                "risk_allele_copies": copies,
                "p_value": a.p_value,
                "pubmed_id": a.pubmed_id,
                "study_accession": a.study,
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"gwas:{a.rsid}-{a.risk_allele}:"
                + ("+".join(a.efo_ids) or a.trait),
                "kind": str(kind),
                "vid": vid,
                "gene": a.gene,
                "claim": json.dumps(claim, ensure_ascii=False),
                "validity": str(Validity.LIMITED),
                "source_level": f"GWAS Catalog p = {a.p_value}",
                "reason": None,
                "effect": json.dumps({"relative": relative, "absolute": None}),
                "actionable": False,
                "topic": str(topic),
                "source": GWAS_TOOL,
                "source_id": a.study,
                "source_version": self.version,
                "points": None,
            }
            variant: Fields = {"run_id": ctx.run_id, "vid": vid, "rsid": a.rsid}
            call = panel.call_row(ctx.run_id, vid, called, reliability)
            rows.append(Row(variant, call, statement))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=self.input_rows,
            dropped=dropped,
        )
