"""The key paths of the query's JSON contract, which the genome-ask skill documents."""

from typing import Any

# Every key path of `query --json` (with `--list-unassessed`) and `show --json`, nested ones
# included. A finding's keys are rooted at the top, whether it sits in `findings`, in
# `unassessed_sites` or in `show`'s `finding`. `versions` (tool -> version), `hidden`
# (reason -> count), `topics` (topic -> count) and `reasons` (reason -> count) are maps:
# their keys are data. `by_gene` maps a gene (data) to an object whose keys are the
# contract's.
CONTRACT = {
    *("run_id", "versions", "findings", "hidden", "finding"),
    *("unassessed", "unassessed_sites"),
    *("unassessed.total", "unassessed.by_gene"),
    *("unassessed.by_gene.sites", "unassessed.by_gene.topics"),
    "unassessed.by_gene.reasons",
    *("id", "kind", "topic", "gene", "marker", "genotype", "reliability", "validity"),
    *("effect", "actionable", "claim", "sources"),
    *("marker.rsid", "marker.grch38"),
    *("genotype.gt", "genotype.dosage", "genotype.allele"),
    *("reliability.level", "reliability.p_wrong", "reliability.p_wrong_ci"),
    *("reliability.gq", "reliability.flags"),
    *("reliability.warnings", "reliability.row_type"),
    *("validity.level", "validity.source_level", "validity.reason"),
    *("effect.relative", "effect.absolute"),
    *("effect.relative.measure", "effect.relative.value", "effect.relative.ci"),
    *("effect.absolute.risk", "effect.absolute.baseline", "effect.absolute.population"),
    *("claim.schema", "claim.type", "claim.classification", "claim.review_status"),
    *("claim.condition", "claim.variation_id"),
    *("sources.source", "sources.id", "sources.version"),
    "gene_call",
    # Stage 4: the dates from the run snapshots, null without them.
    *("dates", "dates.first_found", "dates.last_changed", "dates.assessment_changed"),
    "dates.on_or_before",
}
# Slice 1b: the keys a pharmacogenomic gene result adds, its claim's and its gene call's.
PGX_CONTRACT = {
    *("claim.gene", "claim.diplotype", "claim.phenotype", "claim.activity_score"),
    *("claim.call_source", "claim.drugs"),
    *("claim.drugs.drug", "claim.drugs.source", "claim.drugs.classification"),
    *("gene_call.diplotype", "gene_call.phenotype", "gene_call.activity_score"),
    *("gene_call.callers", "gene_call.agreed", "gene_call.level", "gene_call.warnings"),
}
# Slice 1b: the keys a PharmGKB clinical annotation's claim adds (its `gene`, `condition`
# and `drugs` are the keys above).
PHARMGKB_CONTRACT = {"claim.annotation_id", "claim.category"}
# Slice 1d: the keys a polygenic score adds, its claim's and its effect's (APOE's claim has
# the pharmacogenomic `gene` and `diplotype`).
PRS_CONTRACT = {
    *("claim.pgs_id", "claim.trait_efo", "claim.publication"),
    *("claim.development_ancestry", "claim.evaluations", "claim.licence"),
    "claim.weak",  # slice 3g
    *("effect.score", "effect.score.percentile", "effect.score.z"),
    *("effect.score.population", "effect.score.match_rate"),
    # Slice 3f: pgsc_calc's Z_norm2 and the decile shown.
    *("effect.score.z_norm2", "effect.score.decile"),
}
# Slice 1d: the keys the mitochondrial haplogroup's claim adds.
HAPLOGROUP_CONTRACT = {
    *("claim.haplogroup", "claim.quality", "claim.tree", "claim.uncovered_positions"),
}
# Slice 1d: the keys a research association's claim adds (its `condition` and `trait_efo`
# are the keys above).
GWAS_CONTRACT = {
    *("claim.risk_allele", "claim.risk_allele_copies", "claim.p_value"),
    *("claim.pubmed_id", "claim.study_accession"),
}
# Slice 1e: the keys a predicted finding's claim adds.
PREDICTED_CONTRACT = {
    *("claim.consequence", "claim.transcript", "claim.hgvsc", "claim.hgvsp"),
}
# Slice 3a: the keys a repeat locus's, SMN1's and the Y haplogroup's claims add (their
# `gene`, `haplogroup` and `tree` are the keys above).
STAGE3_CONTRACT = {
    *("claim.locus_id", "claim.motif", "claim.class", "claim.image", "claim.alleles"),
    *("claim.alleles.repeats", "claim.alleles.ci", "claim.alleles.spanning"),
    *("claim.alleles.flanking", "claim.alleles.inrepeat", "claim.thresholds"),
    *("claim.thresholds.normal_max", "claim.thresholds.intermediate"),
    *("claim.thresholds.reduced_penetrance", "claim.thresholds.pathogenic_min"),
    *("claim.smn1", "claim.smn2", "claim.g27134tg", "claim.callers"),
    *("claim.marker", "claim.qc_score", "claim.valid_markers"),
    # Slice 3e: the test that would confirm a repeat normal flagged `exclusion`, and
    # RFC1's motif call.
    "claim.confirm",
    "claim.rfc1",
}
# Slice 3c: the keys the sex-chromosome dosage statement's claim adds, and a gene's
# coverage, which `query --gene` prints at the top.
DOSAGE_CONTRACT = {
    *("claim.label", "claim.x_ratio", "claim.y_ratio", "claim.x_het_share"),
    "claim.x_sites",
}
_COVERAGE = (
    *("bases", "callable", "mean_depth", "difficult", "outside_benchmark", "owner"),
    *("reliable", "note"),
)
COVERAGE_CONTRACT = {"coverage", *(f"coverage.{key}" for key in _COVERAGE)}
# Slice 3f: the keys its claims add to a run's query (HLA, MITOMAP, HHV-6, ROH, blood
# groups, ancestry, traits; `gene`, `condition`, `drugs`, `callers`, `classification` and
# `haplogroup` are keys above), whatever the statement's kind.
_DQ = ("dq25", "dq22", "dq8", "half")
_RUN = ("chrom", "start", "end", "length")
SLICE3F_CONTRACT = {
    *("claim.allele", "claim.status", "claim.confirm", "claim.note"),
    *("claim.t1k", "claim.tags", "claim.tag_dosages", "claim.genes"),
    *("claim.species", "claim.ratio"),
    *("claim.min_run", "claim.total_long", "claim.n_long", "claim.longest"),
    "claim.isodisomy",
    *("claim.system", "claim.genotype"),
    *("claim.trait", "claim.result", "claim.sites", "claim.sites.rsid"),
    *("claim.sites.vid", "claim.sites.effect_allele", "claim.sites.copies"),
    *("claim.reliability", "claim.sites.genotype", "claim.sites.level"),
}
# Slice 3f: the keys only a result writes (a carried MITOMAP variant, the coeliac readings,
# the HLA table's genes, the runs of homozygosity, the ancestry statements), which a run
# whose 3f statements are all unassessed does not show.
SLICE3F_FINDING_CONTRACT = {
    *("claim.locus", "claim.level"),
    *(f"claim.{side}.{key}" for side in ("t1k", "tags") for key in _DQ),
    *("claim.genes.gene", "claim.genes.alleles", "claim.genes.callers"),
    "claim.genes.note",
    *(f"claim.{run}.{key}" for run in ("longest", "isodisomy") for key in _RUN),
    *("claim.population", "claim.low_confidence", "claim.probabilities"),
    *("claim.proportions", "claim.neanderthal_bp", "claim.denisovan_bp"),
    *("claim.share", "claim.genome_bp"),
}
CONTRACT |= PGX_CONTRACT | PHARMGKB_CONTRACT | PRS_CONTRACT
CONTRACT |= HAPLOGROUP_CONTRACT | GWAS_CONTRACT | PREDICTED_CONTRACT | STAGE3_CONTRACT
CONTRACT |= DOSAGE_CONTRACT | COVERAGE_CONTRACT | SLICE3F_CONTRACT
CONTRACT |= SLICE3F_FINDING_CONTRACT
# The key paths `genotype --json` adds to those: `run_id` and the finding keys it shares.
GENOTYPE_CONTRACT = {
    *("assessed", "sites"),
    *("answer", "classification", "condition", "not_a_finding_because"),
    *("gene.pathogenic_sites", "gene.carried", "gene.by_answer", "gene.by_level"),
    *("gene.deletions", "gene.not_covered", "gene.dominant_mechanism_not_assessed"),
    *(f"gene.{path}" for path in COVERAGE_CONTRACT),
    *("gene.possible_carriers", "gene.copy_number", "gene.insertions"),
}
# Stage 4: every key path of `changes --json`. A change's `before` and `after` hold the
# compared fields and the source's release.
_SIDE = (
    *("kind", "validity", "source_level", "reason", "effect", "actionable", "topic"),
    *("dosage", "level", "result", "source_version"),
)
CHANGES_CONTRACT = {
    *("base", "run", "refresh", "failed", "check", "changes", "research", "note"),
    *("refresh.started", "refresh.ended", "refresh.outcome", "refresh.sources"),
    *("refresh.base_run", "refresh.pre_run", "refresh.run", "refresh.error"),
    *("refresh.sources.name", "refresh.sources.before", "refresh.sources.after"),
    *("refresh.sources.status", "refresh.sources.error"),
    *("failed.name", "failed.before"),
    *("check.ran_at", "check.newer", "check.failed", "check.days", "check.stale"),
    *("check.newer.name", "check.newer.kind", "check.newer.local"),
    *("check.newer.upstream",),
    *("changes.id", "changes.kind", "changes.gene", "changes.marker", "changes.topic"),
    *("changes.before", "changes.after", "changes.fields", "changes.causes"),
    *("changes.dates", "changes.publications"),
    *("changes.marker.rsid", "changes.marker.grch38"),
    *(f"changes.{side}.{key}" for side in ("before", "after") for key in _SIDE),
    *("changes.dates.first_found", "changes.dates.last_changed"),
    *("changes.dates.assessment_changed", "changes.dates.on_or_before"),
    *("changes.publications.new", "changes.publications.attached_by_lab"),
}
# `genotype --json`'s `by_answer` (answer -> count) and `by_level` (level -> count) are maps
# too, and so is a gene call's `callers` (caller -> call), and `changes --json`'s `research`
# (change kind -> count); its `sites` are findings with four keys more, rooted at the top
# like the rest. Slice 3f's HLA claims hold two more, `status` (caller -> status) and
# `tag_dosages` (rsid -> dosage), and its ancestry claims two: `probabilities` and
# `proportions` (population -> share).
_MAPS = frozenset(
    {
        *("versions", "hidden", "topics", "reasons", "by_answer", "by_level"),
        *("callers", "research", "status", "tag_dosages"),
        *("probabilities", "proportions"),
    }
)
# CFTR's repeat claim keys its `alleles` by tract (`tg`, `t`), each an allele list.
_KEYED = frozenset({"by_gene", "alleles"})
_FINDING_ROOTS = frozenset({"findings", "finding", "unassessed_sites", "sites"})


def key_paths(value: Any, prefix: str = "") -> set[str]:
    """Every key path in a JSON value: a map's keys are data, a keyed map's values carry
    the path of the map itself, and a finding's keys are rooted at the top."""
    if isinstance(value, list):
        return set().union(*(key_paths(item, prefix) for item in value))
    if not isinstance(value, dict):
        return set()
    paths: set[str] = set()
    for key, inner in value.items():
        path = f"{prefix}{key}"
        paths.add(path)
        if key in _MAPS:
            continue
        if key in _KEYED and isinstance(inner, dict):
            for item in inner.values():
                paths |= key_paths(item, f"{path}.")
            continue
        paths |= key_paths(
            inner, "" if not prefix and key in _FINDING_ROOTS else f"{path}."
        )
    return paths
