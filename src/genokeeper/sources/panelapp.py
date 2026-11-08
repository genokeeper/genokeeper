"""Rare, probably damaging variants the user carries in the genes of Genomics England
PanelApp's Intellectual disability panel, as research statements (`PanelAppSource`).

The panel (285, version `PANEL_VERSION`, one JSON file the bootstrap pins by sha256) rates
how well each gene's link to the condition is established (green, amber or red: evidence for
the gene, never for a variant) and names its mode of inheritance. The layer reads no tool of its own: it reads the records `PredictedSource`
graded, so VEP's and SpliceAI's steps and their keys do not depend on it, and SpliceAI scored
only the listed genes (`genelists.LISTED`), which a splice change's statement says.

A statement (`panelapp:<vid>:<transcript>`, as a variant can have records on two, `predicted:`'s
convention) stands for a record in a panel gene that `predicted_rule` calls a loss of function
that passes every LoF gate, or a missense at its calibrated REVEL or AlphaMissense score;
that `grade_predicted` found carried, and rare by this layer's own limit (`MAX_FAF`); every
other record is counted under why. A variant a `predicted:` statement already states keeps it alone
(`without_repeats`, which runs after the compound-heterozygote promotion, as that can make a
statement of one). Limited, under research and never actionable: a computer's prediction on
a broad panel, which most people carry something on.

The layer has one boundary (`pipeline._panelapp_layer`), inside which its statements are
computed once (`collect_panel`): anything wrong there leaves one unassessed statement
(`FAILED_STATEMENT`) and every other statement as it is. The store phase only puts the run's
id on rows already made (`PanelAppSource`).
"""

import json
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path

from genokeeper.annotate.vep import VepRecord, locus
from genokeeper.evidence.clinvar import GENOTYPE_REASON
from genokeeper.evidence.predicted import (
    GENETICIST_NOTE,
    LOF,
    MISSENSE,
    NOT_CARRIED,
    PREDICTED_REASON,
    SPLICE_CONSEQUENCES,
    Site,
    failed_gate,
    predicted_rule,
)
from genokeeper.genelists import LISTED
from genokeeper.sources import Collected, Fields, Row, RunContext, panel
from genokeeper.sources.predicted import PREDICTED_TOOL, PredictedSource
from genokeeper.store.db import Kind, Topic, Validity

PANELAPP_TOOL = "panelapp"
PANEL_ID = 285
PANEL_VERSION = "11.93"
CLASSIFICATION = "predicted, not classified"
PANEL_SOURCE_LEVEL = "Genomics England PanelApp"
# The one statement of a layer that could not run.
FAILED_STATEMENT = "panelapp:panel"
# This layer's one rarity limit, whatever the gene's inheritance: gnomAD's filtering allele
# frequency below 1%. `grade_predicted`'s rule (BS1, and BA1's 0.05 alone for a recessive
# gene, where a pathogenic allele may be common) answers another question, and would call
# an allele at 4% rare.
MAX_FAF = 0.01
NOT_A_DIAGNOSIS = "this is not a diagnosis"
ONE_COPY_NOTE = "one copy of a variant in a recessive gene is not a diagnosis"
BIALLELIC = "BIALLELIC"
SPLICEAI_NOTE = (
    "SpliceAI scores only the ACMG SF and carrier genes, so it did not score this splice"
    " change"
)
NO_INHERITANCE = "not stated by the panel"
# Why a record is no statement.
NOT_IN_PANEL = "not_in_panel"
NOT_A_PANEL_RULE = "not_a_panel_rule"
ABSENT = "not_carried"
POPULATION_FREQUENCY = "population_frequency"
FAILED_LOF_GATE = "failed_lof_gate"
ALREADY_PREDICTED = "already_predicted"
_DROPS = (
    NOT_IN_PANEL,
    NOT_A_PANEL_RULE,
    ABSENT,
    POPULATION_FREQUENCY,
    FAILED_LOF_GATE,
    ALREADY_PREDICTED,
)
# PanelApp's confidence level -> its colour.
_LEVELS = {"3": "green", "2": "amber", "1": "red", "0": "red"}
# What each colour says, of the gene and never of the user's variant.
_GENE_EVIDENCE = {
    "green": "its link to the condition is well established",
    "amber": "the evidence for its link to the condition is moderate",
    "red": "the evidence for its link to the condition is weak",
}


def panelapp_path(lab: Path) -> Path:
    """Where the bootstrap files the panel (`scripts/bootstrap.sh`, `step_panelapp`)."""
    return lab / "panelapp" / str(PANEL_ID) / PANEL_VERSION / "panel.json"


@dataclass(frozen=True)
class PanelGene:
    """One gene's confidence colour and its inheritance as the panel writes it."""

    confidence: str
    inheritance: str


@dataclass(frozen=True)
class GenePanel:
    """The panel's name and version, its genes by symbol, and its genes' GRCh38 Ensembl ids
    (a panel's HGNC symbol can be older than VEP's)."""

    name: str
    version: str
    genes: dict[str, PanelGene]
    by_ensembl: dict[str, str] = field(default_factory=dict)

    def gene(self, record: VepRecord) -> PanelGene | None:
        """The panel's entry for `record`'s gene: by VEP's symbol, else by its Ensembl id."""
        found = self.genes.get(record.gene)
        if found is not None:
            return found
        symbol = self.by_ensembl.get(record.gene_id)
        return None if symbol is None else self.genes[symbol]


def load_panel(path: Path) -> GenePanel:
    """The panel file; refuse a missing file, another panel or version, no genes, a
    duplicate gene or an unknown confidence level."""
    if not path.is_file():
        raise FileNotFoundError(f"missing input: no {path}; run scripts/bootstrap.sh")
    data = json.loads(path.read_text())
    if data["id"] != PANEL_ID:
        raise ValueError(f"{path}: panel {data['id']}, expected {PANEL_ID}")
    if data["version"] != PANEL_VERSION:
        raise ValueError(f"{path}: version {data['version']}, expected {PANEL_VERSION}")
    genes: dict[str, PanelGene] = {}
    by_ensembl: dict[str, str] = {}
    for entry in data["genes"]:
        gene_data = entry["gene_data"]
        symbol = gene_data["gene_symbol"]
        if symbol in genes:
            raise ValueError(f"{path}: duplicate gene {symbol}")
        level = entry["confidence_level"]
        if level not in _LEVELS:
            raise ValueError(f"{path}: unknown confidence level {level!r} of {symbol}")
        inheritance = entry["mode_of_inheritance"] or NO_INHERITANCE
        genes[symbol] = PanelGene(_LEVELS[level], inheritance)
        for release in gene_data["ensembl_genes"]["GRch38"].values():
            by_ensembl[release["ensembl_id"]] = symbol
    if not genes:
        raise ValueError(f"{path}: no genes")
    return GenePanel(data["name"], data["version"], genes, by_ensembl)


def _dropped_reason(
    rule: str, site: Site, grade_reason: str, ctx: RunContext
) -> str | None:
    """Why a record under a panel rule is no statement, or None when it qualifies:
    `grade_predicted`'s carried test, this layer's rarity limit, and the LoF filter."""
    if grade_reason == NOT_CARRIED:
        return ABSENT
    if site.faf is not None and site.faf >= MAX_FAF:
        return POPULATION_FREQUENCY
    if rule == LOF and failed_gate(site, ctx.grading.predicted) is not None:
        return FAILED_LOF_GATE
    return None


def _statement(
    record: VepRecord,
    site: Site,
    found: tuple[str, str],
    gene_panel: GenePanel,
    panel_gene: PanelGene,
    ctx: RunContext,
) -> Fields:
    rule, source_level = found
    unknown = site.dosage == "unknown"
    listed = any(record.gene in ctx.lists[name].genes for name in LISTED)
    chrom, pos, _, _ = locus(record.vid)
    # A het, on two chromosomes, in a gene the panel calls biallelic: not a homozygote, not
    # a hemizygous male, not a monoallelic or X-linked gene, where one copy can matter.
    one_copy = (
        site.dosage == "1"
        and not ctx.haploid(chrom, pos)
        and panel_gene.inheritance.startswith(BIALLELIC)
    )
    rating = (
        f"gene rated {panel_gene.confidence} on the panel:"
        f" {_GENE_EVIDENCE[panel_gene.confidence]}; this says nothing about this variant"
    )
    reasons = [
        PREDICTED_REASON,
        rule,
        f"in Genomics England PanelApp's {gene_panel.name} panel ({PANEL_ID}, v{gene_panel.version})",
        rating,
        f"inheritance: {panel_gene.inheritance}",
        *([SPLICEAI_NOTE] if record.terms & SPLICE_CONSEQUENCES and not listed else []),
        *([GENOTYPE_REASON] if unknown else []),
        NOT_A_DIAGNOSIS,
        *([ONE_COPY_NOTE] if one_copy else []),
        GENETICIST_NOTE,
    ]
    claim = {
        "schema": panel.CLAIM_SCHEMA,
        "type": "panelapp_variant",
        "classification": CLASSIFICATION,
        "consequence": record.consequence,
        "transcript": record.transcript,
        "hgvsc": record.hgvsc,
        "hgvsp": record.hgvsp,
        "panel": gene_panel.name,
        "panel_version": gene_panel.version,
        "panel_confidence": panel_gene.confidence,
        "panel_inheritance": panel_gene.inheritance,
    }
    transcript = record.transcript.split(".")[0]
    return {
        "run_id": ctx.run_id,
        "statement_id": f"panelapp:{record.vid}:{transcript}",
        "kind": str(Kind.UNASSESSED if unknown else Kind.FINDING),
        "vid": record.vid,
        "gene": record.gene,
        "claim": json.dumps(claim),
        "validity": str(Validity.LIMITED),
        "source_level": source_level,
        "reason": "; ".join(reasons),
        "effect": None,
        "actionable": False,  # a prediction on a broad panel: never actionable
        "topic": str(Topic.RESEARCH),
        "source": PANELAPP_TOOL,
        "source_id": record.transcript,
        "source_version": gene_panel.version,
        "points": None,
    }


def failed_panel(ctx: RunContext, failed: str) -> Collected:
    """The one unassessed statement of a layer that could not run: what failed."""
    claim = {
        "schema": panel.CLAIM_SCHEMA,
        "type": "panelapp_variant",
        "classification": CLASSIFICATION,
    }
    statement: Fields = {
        "run_id": ctx.run_id,
        "statement_id": FAILED_STATEMENT,
        "kind": str(Kind.UNASSESSED),
        "vid": None,
        "gene": None,
        "claim": json.dumps(claim),
        "validity": str(Validity.LIMITED),
        "source_level": PANEL_SOURCE_LEVEL,
        "reason": f"the panel's research layer could not run: {failed}",
        "effect": None,
        "actionable": False,
        "topic": str(Topic.RESEARCH),
        "source": PANELAPP_TOOL,
        "source_id": str(PANEL_ID),
        "source_version": PANEL_VERSION,
        "points": None,
    }
    return Collected(PANELAPP_TOOL, [Row(None, None, statement)], {}, 1, {})


def collect_panel(
    gene_panel: GenePanel, predicted: PredictedSource, ctx: RunContext
) -> Collected:
    """One statement per qualifying record; the rest counted by reason (see the module
    docstring)."""
    dropped: dict[str, int] = dict.fromkeys(_DROPS, 0)
    rows: list[Row] = []
    graded = predicted.graded(ctx)
    for record, called, reliability, site, grade in graded:
        panel_gene = gene_panel.gene(record)
        if panel_gene is None:
            dropped[NOT_IN_PANEL] += 1
            continue
        found = predicted_rule(record, ctx.grading.predicted)
        if found is None or found[0] not in (LOF, MISSENSE):
            dropped[NOT_A_PANEL_RULE] += 1
            continue
        why = _dropped_reason(found[0], site, grade.reason or "", ctx)
        if why is not None:
            dropped[why] += 1
            continue
        variant: Fields = {"run_id": ctx.run_id, "vid": record.vid, "rsid": None}
        call = panel.call_row(ctx.run_id, record.vid, called, reliability)
        statement = _statement(record, site, found, gene_panel, panel_gene, ctx)
        rows.append(Row(variant, call, statement))
    return Collected(PANELAPP_TOOL, rows, {}, len(graded), dropped)


def _for_run(fields: Fields | None, run_id: str) -> Fields | None:
    return None if fields is None else {**fields, "run_id": run_id}


@dataclass(frozen=True)
class PanelAppSource:
    """What the layer's boundary computed (`collect_panel`, or `failed_panel`), under its
    placeholder run id."""

    collected: Collected
    tool: str = PANELAPP_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """The rows the boundary made, with the run's id: nothing here can fail."""
        rows = [
            replace(
                row,
                variant=_for_run(row.variant, ctx.run_id),
                call=_for_run(row.call, ctx.run_id),
                statement={**row.statement, "run_id": ctx.run_id},
            )
            for row in self.collected.rows
        ]
        return replace(self.collected, rows=rows)


def without_repeats(collected: Sequence[Collected]) -> list[Collected]:
    """`collected` with every panel statement dropped whose variant a `predicted:`
    statement states: after the compound-heterozygote promotion, which can make one. A
    `predicted:` non-finding (one copy in a recessive gene) states nothing to the user, so
    it does not count."""
    stated = {
        row.statement["vid"]
        for c in collected
        if c.tool == PREDICTED_TOOL
        for row in c.rows
        if row.statement["kind"] != Kind.NON_FINDING
    }
    result: list[Collected] = []
    for c in collected:
        if c.tool != PANELAPP_TOOL:
            result.append(c)
            continue
        kept = [row for row in c.rows if row.statement["vid"] not in stated]
        gone = c.dropped.get(ALREADY_PREDICTED, 0) + len(c.rows) - len(kept)
        result.append(
            replace(c, rows=kept, dropped={**c.dropped, ALREADY_PREDICTED: gone})
        )
    return result
