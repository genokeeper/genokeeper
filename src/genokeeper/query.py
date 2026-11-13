"""Findings of one run, filtered by reliability and validity, in the spec's order, and the
run's unassessed claims, counted per gene.

`query` returns the findings a filter keeps and counts the ones it hides, each once under the
first filter it fails: `research_association`, `reliability_below:<threshold>`,
`validity_below:<threshold>`, `topic`, `gene`, `marker`, `text`. A single-SNP research
association (source `gwas_catalog`, slice 1d) is kept only by a filter that names a gene, a
marker or a text: a user carries the common risk allele of most of the catalog's hundreds
of thousands, so unnarrowed they would bury every other finding (design decision,
2026-10-23). A gene-level finding faces the reliability filter with its gene
call's level (`gene_call`, a pharmacogenomic gene result, slice 1b), and without one never.

A claim at a site whose genotype is unknown (dosage `unknown`, or no call row) is not a
finding but `unassessed` (the store's kind). The unassessed claims the same filters keep,
reliability aside, are counted per gene with their topics and why the genotype is unknown
(`unassessed`) and, on request, listed in the findings' shape (`unassessed_sites`); those
the filters drop are counted in `hidden` under the same reason prefixed `unassessed:`, so
no claim disappears uncounted. Every listed element names its `kind`.

With a gene, `query` also returns the gene's coverage in the run (`gene_coverage`, slice
3c): its `coverage` row with a `note` — "not reliably assessed: N% callable" below
`[coverage] min_callable`, "assessed by <owner>" for a gene a stage-3 source states,
whatever its callable share, "not measured" for a gene the run has no row for — so that "no
finding" in a gene the reads barely cover never reads as reassurance.

A finding's `genotype.allele` is the allele its `dosage` counts: the GRCh38 ALT. The JSON
shape is the contract the genome-ask skill reads; the text output adds a note on chrX and chrY
findings that there, for a male outside the PARs, dosage counts copies on a single-copy
chromosome.
"""

import argparse
import json
import sys
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper.config import load_paths
from genokeeper.evidence.pgx import GENE_UNRESOLVED
from genokeeper.insertions import INSERTIONS_TOOL
from genokeeper.sources.cnv import CNV_TOOL
from genokeeper.sources.deletions import DELETIONS_TOOL
from genokeeper.sources.panel import WITHHELD
from genokeeper.sources.paralogs import LPA_TOOL, PARALOGS_TOOL
from genokeeper.sources.predicted import PREDICTED_TOOL
from genokeeper.sources.repeats import REPEATS_TOOL
from genokeeper.sources.research import is_research_layer
from genokeeper.sources.smn import SMN_TOOL
from genokeeper.sources.yleaf import YLEAF_TOOL
from genokeeper.store.db import (
    Kind,
    Level,
    StoreTooNew,
    StoreTooOld,
    Topic,
    Validity,
    check_schema_version,
)

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers
    from genokeeper.history import Dates

# A male's chrX and chrY outside the PARs hold one copy: the text output says what the
# dosage means there (the query knows neither the sample's sex nor the PARs).
HAPLOID_DOSAGE_NOTE = (
    "on a male's chrX or chrY outside the PARs, dosage counts copies on a single-copy"
    " chromosome (1 = hemizygous carrier of the allele)"
)
_SEX_CONTIGS = frozenset({"chrX", "chrY"})

_LEVEL_RANK = {Level.LOW: 0, Level.MEDIUM: 1, Level.HIGH: 2}
_VALIDITY_RANK = {
    Validity.CONFLICTING: 0,
    Validity.LIMITED: 1,
    Validity.MODERATE: 2,
    Validity.STRONG: 3,
}

# Hidden-count reasons in filter order; a threshold reason carries `:<threshold>`. An
# unassessed claim's reason carries the prefix, and those come after the findings'.
_REASONS = (
    "research_association",
    "reliability_below",
    "validity_below",
    "topic",
    "gene",
    "marker",
    "text",
)
RESEARCH_ASSOCIATION = "research_association"
UNASSESSED_PREFIX = "unassessed:"
# Why an unassessed claim's genotype is unknown: its call's first warning that starts with
# this, else `NO_CALL_REASON` when it has no call row.
GENOTYPE_UNKNOWN = "genotype unknown"
NO_CALL_REASON = "no call row"
# A gene-level unassessed claim's reason: its gene call's first warning that starts with
# `GENE_UNRESOLVED`; without a gene call, the reason it carries, else `NO_GENE_CALL`.
NO_GENE_CALL = "no gene call"
# A gene call is the result of the one gene-level statement whose id is `<tool>:<its name>`
# (PharmCAT's gene result, APOE, SMN1, a repeat locus): another gene-level statement on that
# gene (an HLA allele's, PharmCAT's G6PD gap, slice 3f) is not joined to it. The name is the
# gene's, or starts with it where two tools answer one gene (`sources.repeats.CFTR_CALL`), so
# a statement without a gene owns none (RBCeq2's `bloodgroup:ABCG2` is not PharmCAT's ABCG2).
# `{s}` is the statement's alias, `g` the gene call's.
OWNS_GENE_CALL = (
    "ends_with({s}.statement_id, ':' || g.gene) AND starts_with(g.gene, {s}.gene)"
)

_COLUMNS = (
    "f.statement_id, f.kind, f.topic, f.gene, f.vid, f.claim, f.validity, f.source_level,"
    " f.reason, f.effect, f.actionable, f.source, f.source_id, f.source_version, f.gt,"
    " f.dosage, f.p_wrong, f.p_wrong_lo, f.p_wrong_hi, f.gq, f.level, f.flags, f.warnings,"
    " f.row_type, v.rsid, g.gene, g.diplotype, g.phenotype, g.activity_score, g.callers,"
    " g.agreed, g.level, g.warnings"
)


def _select(kinds: Collection[Kind]) -> str:
    """The columns of the statements of `kinds` with their call (the `finding` and
    `unassessed` views' join, over any kinds), the run's variant row and, for a gene-level
    statement, its gene call."""
    listed = ", ".join(f"'{Kind(kind)}'" for kind in kinds)
    return (
        f"SELECT {_COLUMNS} FROM (SELECT s.*, c.* EXCLUDE (run_id, vid)"
        " FROM statement AS s"
        " LEFT JOIN call AS c ON c.run_id = s.run_id AND c.vid = s.vid"
        f" WHERE s.kind IN ({listed})) AS f"
        " LEFT JOIN variant AS v ON v.run_id = f.run_id AND v.vid = f.vid"
        " LEFT JOIN gene_call AS g"
        " ON f.vid IS NULL AND g.run_id = f.run_id"
        f" AND {OWNS_GENE_CALL.format(s='f')}"
    )


class RunNotFound(LookupError):
    """The store is missing or unreadable (another schema version wrote it), or has no such
    run, or no run at all."""


@dataclass(frozen=True)
class Filters:
    """What a query keeps; the defaults are the dashboard's default view."""

    topic: Topic | None = None
    min_reliability: Level = Level.MEDIUM
    min_validity: Validity = Validity.MODERATE
    gene: str | None = None
    marker: str | None = None
    text: str | None = None


@dataclass(frozen=True)
class QueryResult:
    """The selected run, its tool versions, the kept findings in order, the hidden counts,
    the kept unassessed claims per gene, when listed those claims in order and, for a
    query by gene, the gene's coverage (`gene_coverage`)."""

    run_id: str
    versions: dict[str, str]
    findings: list[dict[str, Any]]
    hidden: dict[str, int]
    unassessed: dict[str, Any]
    unassessed_sites: list[dict[str, Any]] | None = None
    coverage: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        document = {
            "run_id": self.run_id,
            "versions": self.versions,
            "findings": self.findings,
            "hidden": self.hidden,
            "unassessed": self.unassessed,
        }
        if self.unassessed_sites is not None:
            document["unassessed_sites"] = self.unassessed_sites
        if self.coverage is not None:
            document["coverage"] = self.coverage
        return document


COVERAGE_COLUMNS = (
    *("bases", "callable", "mean_depth", "difficult", "outside_benchmark", "owner"),
    "reliable",
)
NOT_MEASURED = "not measured"


def coverage_note(row: Mapping[str, Any]) -> str | None:
    """What a gene's coverage row means for its answer: "assessed by <owner>" for a gene a
    stage-3 source states (its own calibration measures it, whatever the callable share),
    "not reliably assessed: N% callable" below `[coverage] min_callable`, else None."""
    if row["owner"] is not None:
        return f"assessed by {row['owner']}"
    if not row["reliable"]:
        return f"not reliably assessed: {row['callable']:.0%} callable"
    return None


def gene_coverage(
    con: duckdb.DuckDBPyConnection, run_id: str, gene: str
) -> dict[str, Any]:
    """The run's `coverage` row of `gene` (case-insensitive) with its `note`
    (`coverage_note`); every column null and the note "not measured" without a row."""
    row = con.execute(
        f"SELECT {', '.join(COVERAGE_COLUMNS)} FROM coverage"
        " WHERE run_id = ? AND lower(gene) = lower(?)",
        [run_id, gene],
    ).fetchone()
    if row is None:
        return dict.fromkeys(COVERAGE_COLUMNS) | {"note": NOT_MEASURED}
    found = dict(zip(COVERAGE_COLUMNS, row, strict=True))
    return found | {"note": coverage_note(found)}


def _run(
    con: duckdb.DuckDBPyConnection, run_id: str | None
) -> tuple[str, dict[str, str]]:
    if run_id is None:
        row = con.execute(
            "SELECT run_id, tools FROM run ORDER BY started_at DESC, run_id DESC LIMIT 1"
        ).fetchone()
        if row is None:
            raise RunNotFound("the store has no runs")
    else:
        row = con.execute(
            "SELECT run_id, tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
        if row is None:
            raise RunNotFound(f"no run {run_id} in the store")
    return row[0], json.loads(row[1]) if row[1] is not None else {}


def _load(value: str | None) -> Any:
    return json.loads(value) if value is not None else None


def _finding(
    row: tuple[Any, ...], dates: Mapping[str, Dates] | None = None
) -> dict[str, Any]:
    """A row of `_select` in the JSON contract's shape; `dates` are the snapshots' dates
    by statement id (`dates` null without them, or for a statement they lack)."""
    (
        statement_id,
        kind,
        topic,
        gene,
        vid,
        claim,
        validity,
        source_level,
        reason,
        effect,
        actionable,
        source,
        source_id,
        source_version,
        gt,
        dosage,
        p_wrong,
        p_wrong_lo,
        p_wrong_hi,
        gq,
        level,
        flags,
        warnings,
        row_type,
        rsid,
        call_gene,
        diplotype,
        phenotype,
        activity_score,
        callers,
        agreed,
        gene_level,
        gene_warnings,
    ) = row
    has_call = level is not None
    grch38 = vid.removeprefix("GRCh38:") if vid is not None else None
    gene_call = (
        None
        if call_gene is None
        else {
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": activity_score,
            "callers": _load(callers),
            "agreed": agreed,
            "level": gene_level,
            "warnings": list(gene_warnings or []),
        }
    )
    return {
        "id": statement_id,
        "kind": kind,
        "topic": topic,
        "gene": gene,
        "marker": {"rsid": rsid, "grch38": grch38},
        # `dosage` counts copies of the GRCh38 ALT.
        "genotype": {
            "gt": gt,
            "dosage": dosage,
            "allele": grch38.rsplit(":", 1)[1] if grch38 is not None else None,
        },
        "reliability": {
            "level": level,
            "p_wrong": p_wrong,
            "p_wrong_ci": None if p_wrong_lo is None else [p_wrong_lo, p_wrong_hi],
            "gq": gq,
            "flags": list(flags or []) if has_call else None,
            "warnings": list(warnings or []) if has_call else None,
            "row_type": row_type,
        },
        "validity": {"level": validity, "source_level": source_level, "reason": reason},
        "effect": _load(effect),
        "actionable": actionable,
        "claim": _load(claim),
        "sources": [{"source": source, "id": source_id, "version": source_version}],
        "gene_call": gene_call,
        "dates": found.to_json()
        if dates is not None and (found := dates.get(statement_id)) is not None
        else None,
    }


def _strip_chr(position: str) -> str:
    return position[3:] if position.lower().startswith("chr") else position


def _marker_matches(finding: dict[str, Any], marker: str) -> bool:
    m = finding["marker"]
    wanted = marker.lower()
    if m["rsid"] is not None and m["rsid"].lower() == wanted:
        return True
    # A GRCh38 vid, with or without its `GRCh38:` prefix and the `chr` of its contig.
    vid = _strip_chr(wanted.removeprefix("grch38:"))
    return m["grch38"] is not None and vid == _strip_chr(m["grch38"].lower())


def text_fields(finding: dict[str, Any]) -> list[str]:
    """What `--text` searches in a finding (the dashboard's search reads the same): its
    gene, the claim's condition and classification, its id, and each drug's name."""
    claim = finding["claim"] or {}
    haystack = [
        finding["gene"],
        claim.get("condition"),
        claim.get("classification"),
        finding["id"],
        # A pharmacogenomic claim's drugs (slice 1b).
        *(d.get("drug") for d in claim.get("drugs") or [] if isinstance(d, dict)),
    ]
    return [str(field) for field in haystack if field is not None]


def searchable(text: str) -> str:
    """`text` as both searches match it: lower-cased, and ClinVar's `_` and `|` as the
    dashboard shows them, a space and ", "."""
    return text.lower().replace("_", " ").replace("|", ", ")


def _text_matches(finding: dict[str, Any], text: str) -> bool:
    wanted = searchable(text)
    return any(wanted in searchable(field) for field in text_fields(finding))


def _hidden_by(
    finding: dict[str, Any], f: Filters, *, reliability: bool = True
) -> str | None:
    """The first filter `finding` fails, as its hidden-count reason, or None if it is kept;
    `reliability=False` skips the reliability filter (an unassessed claim has no genotype
    to grade). A gene-level finding's level is its gene call's, or its claim's
    `reliability` (a trait, slice 3f)."""
    if is_research_layer(finding["sources"][0]["source"]) and not (
        f.gene or f.marker or f.text
    ):
        return RESEARCH_ASSOCIATION
    claim = finding["claim"] if isinstance(finding["claim"], dict) else {}
    level = (
        finding["reliability"]["level"]
        or (finding["gene_call"] or {}).get("level")
        or claim.get("reliability")
    )
    if (
        reliability
        and level is not None
        and _LEVEL_RANK[Level(level)] < _LEVEL_RANK[f.min_reliability]
    ):
        return f"reliability_below:{f.min_reliability}"
    validity = finding["validity"]["level"]
    if _VALIDITY_RANK[Validity(validity)] < _VALIDITY_RANK[f.min_validity]:
        return f"validity_below:{f.min_validity}"
    if f.topic is not None and finding["topic"] != f.topic:
        return "topic"
    if f.gene is not None and (finding["gene"] or "").lower() != f.gene.lower():
        return "gene"
    if f.marker is not None and not _marker_matches(finding, f.marker):
        return "marker"
    if f.text is not None and not _text_matches(finding, f.text):
        return "text"
    return None


def _descending(value: float | None) -> tuple[bool, float]:
    """A sort key putting larger values first and None last."""
    return (value is None, -value if value is not None else 0.0)


def _effect_value(effect: dict[str, Any] | None, part: str, key: str) -> float | None:
    inner = (effect or {}).get(part)
    return inner.get(key) if inner else None


# Slice 3d's sources: depth copy number, insertions, the paralog loci and Lp(a).
_SLICE3D_TOOLS = frozenset({CNV_TOOL, INSERTIONS_TOOL, PARALOGS_TOOL, LPA_TOOL})
_STAGE3_TOOLS = (
    frozenset({REPEATS_TOOL, SMN_TOOL, YLEAF_TOOL, DELETIONS_TOOL}) | _SLICE3D_TOOLS
)


def _tier(source: str) -> int:
    """0 for every source but the stage-3 tools and the predicted one, 1 for the stage-3
    screening tools (3a's, 3b's deletions and 3d's), 2 for predicted."""
    if source == PREDICTED_TOOL:
        return 2
    return 1 if source in _STAGE3_TOOLS else 0


def _order(finding: dict[str, Any]) -> tuple[Any, ...]:
    """The feed's order: a tier first (a stage-3 screening result after every other source
    and a predicted finding after those, so neither outranks a ClinVar-classified one); then
    actionable first, validity, absolute and relative effect, gene and id."""
    effect = finding["effect"]
    return (
        _tier(finding["sources"][0]["source"]),
        finding["actionable"] is not True,
        -_VALIDITY_RANK[Validity(finding["validity"]["level"])],
        _descending(_effect_value(effect, "absolute", "risk")),
        _descending(_effect_value(effect, "relative", "value")),
        (finding["gene"] is None, finding["gene"] or ""),
        finding["id"],
    )


def unassessed_reason(site: dict[str, Any]) -> str:
    """Why an unassessed claim's genotype is unknown: its call's first warning starting with
    `GENOTYPE_UNKNOWN`, or `NO_CALL_REASON` without a call row.

    A run stores an unassessed claim only beside an unknown call (`pipeline.check_calls`),
    and every unknown call carries such a warning (`genotype.gvcf`, `genokeeper.mito`); a call
    without one is a broken store, refused (`ValueError`) rather than counted under a
    reason the contract does not name.

    A gene-level claim (no marker) waits on its gene call instead: its first warning
    starting with `GENE_UNRESOLVED` (`ValueError` for a gene call without one, as for a
    call). Without a gene call, it says itself why it is unassessed (a polygenic score below
    its overlap, an APOE double heterozygote unphased): its stored reason, else
    `NO_GENE_CALL`."""
    if site["marker"]["grch38"] is None:
        gene_call = site["gene_call"]
        if gene_call is None:
            return site["validity"]["reason"] or NO_GENE_CALL
        unresolved = next(
            (w for w in gene_call["warnings"] if w.startswith(GENE_UNRESOLVED)), None
        )
        if unresolved is None:
            raise ValueError(
                f"unassessed claim {site['id']} has a gene call without a"
                f" {GENE_UNRESOLVED!r} warning: the store is inconsistent"
            )
        return unresolved
    warnings = site["reliability"]["warnings"]
    if warnings is None:
        return NO_CALL_REASON
    found = next((w for w in warnings if w.startswith(GENOTYPE_UNKNOWN)), None)
    if found is None:
        raise ValueError(
            f"unassessed claim {site['id']} has a call without a"
            f" {GENOTYPE_UNKNOWN!r} warning: the store is inconsistent"
        )
    return found


def _counted(counts: dict[str, int]) -> dict[str, int]:
    """Counts largest first, ties by key."""
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _unassessed_summary(sites: list[dict[str, Any]]) -> dict[str, Any]:
    """`{"total", "by_gene": {gene: {"sites", "withheld", "topics", "reasons"}}}`, genes
    in order, `withheld` counting the claims holding what a table withheld; a claim without
    a gene counts under `""`."""
    by_gene: dict[str, dict[str, Any]] = {}
    for site in sites:
        entry = by_gene.setdefault(
            site["gene"] or "", {"sites": 0, "withheld": 0, "topics": {}, "reasons": {}}
        )
        entry["sites"] += 1
        entry["withheld"] += WITHHELD in (site["claim"] or {})
        for part, key in (
            ("topics", site["topic"]),
            ("reasons", unassessed_reason(site)),
        ):
            entry[part][key] = entry[part].get(key, 0) + 1
    return {
        "total": len(sites),
        "by_gene": {
            gene: {
                "sites": entry["sites"],
                "withheld": entry["withheld"],
                "topics": _counted(entry["topics"]),
                "reasons": _counted(entry["reasons"]),
            }
            for gene, entry in sorted(by_gene.items())
        },
    }


def _hidden_order(reason: str) -> tuple[bool, int]:
    unassessed = reason.startswith(UNASSESSED_PREFIX)
    base = reason.removeprefix(UNASSESSED_PREFIX).split(":")[0]
    return unassessed, _REASONS.index(base)


def query(
    con: duckdb.DuckDBPyConnection,
    run_id: str | None,
    f: Filters,
    *,
    list_unassessed: bool = False,
    dates: Mapping[str, Dates] | None = None,
) -> QueryResult:
    """The findings of `run_id` (None: the latest run) that `f` keeps, what it hid, and the
    unassessed claims it keeps (reliability aside) counted per gene; `list_unassessed`
    also lists those claims. Each carries its `dates` entry (`history.dates`). With
    `f.gene`, the gene's coverage too."""
    run, versions = _run(con, run_id)
    hidden: dict[str, int] = {}
    kept: dict[Kind, list[dict[str, Any]]] = {Kind.FINDING: [], Kind.UNASSESSED: []}
    for kind, prefix in ((Kind.FINDING, ""), (Kind.UNASSESSED, UNASSESSED_PREFIX)):
        rows = con.execute(f"{_select([kind])} WHERE f.run_id = ?", [run]).fetchall()
        for row in (_finding(r, dates) for r in rows):
            reason = _hidden_by(row, f, reliability=kind is Kind.FINDING)
            if reason is None:
                kept[kind].append(row)
            else:
                hidden[prefix + reason] = hidden.get(prefix + reason, 0) + 1
    sites = sorted(kept[Kind.UNASSESSED], key=_order)
    return QueryResult(
        run_id=run,
        versions=versions,
        findings=sorted(kept[Kind.FINDING], key=_order),
        hidden=dict(sorted(hidden.items(), key=lambda kv: _hidden_order(kv[0]))),
        unassessed=_unassessed_summary(sites),
        unassessed_sites=sites if list_unassessed else None,
        coverage=None if f.gene is None else gene_coverage(con, run, f.gene),
    )


def show(
    con: duckdb.DuckDBPyConnection,
    run_id: str | None,
    statement_id: str,
    *,
    dates: Mapping[str, Dates] | None = None,
) -> tuple[str, dict[str, Any] | None]:
    """The selected run and its statement `statement_id` (None if it has none), whatever
    the filters and whatever its kind; its `kind` says which.

    `run_id=None` selects the latest run.
    """
    run, _versions = _run(con, run_id)
    row = con.execute(
        f"{_select(list(Kind))} WHERE f.run_id = ? AND f.statement_id = ?",
        [run, statement_id],
    ).fetchone()
    return run, _finding(row, dates) if row is not None else None


# --- command line ---------------------------------------------------------------------------


def _text_block(finding: dict[str, Any]) -> str:
    m, g, r, v = (
        finding["marker"],
        finding["genotype"],
        finding["reliability"],
        finding["validity"],
    )
    marker = " ".join(
        f"{label}{value}"
        for label, value in (
            ("", m["rsid"]),
            ("GRCh38 ", m["grch38"]),
        )
        if value is not None
    )
    lines = [f"{finding['gene'] or '(no gene)'}  {marker or '(gene-level)'}"]
    lines.append(
        f"  id: {finding['id']}  kind: {finding['kind']}  topic: {finding['topic']}"
    )
    gene_call = finding["gene_call"]
    if m["grch38"] is None and gene_call is not None:
        callers = "; ".join(f"{k} {v}" for k, v in (gene_call["callers"] or {}).items())
        lines.append(
            f"  gene call: {gene_call['diplotype'] or 'unresolved'}"
            f" ({gene_call['phenotype'] or 'no phenotype'}),"
            f" activity score {gene_call['activity_score']}; reliability"
            f" {gene_call['level']}; callers: {callers or 'none'}"
        )
        if gene_call["warnings"]:
            lines.append(f"  warnings: {', '.join(gene_call['warnings'])}")
    elif m["grch38"] is None:
        lines.append("  genotype: none (gene-level finding)")
    elif r["level"] is None:
        lines.append("  genotype: no call row - unknown, not a result")
    else:
        if g["dosage"] == "unknown":
            lines.append(
                f"  genotype: {g['gt']} (dosage unknown, allele {g['allele']})"
                " - unknown, not a result"
            )
        else:
            lines.append(
                f"  genotype: {g['gt']} (dosage {g['dosage']} = copies of {g['allele']})"
            )
        if m["grch38"].split(":")[0] in _SEX_CONTIGS:
            lines.append(f"  note: {HAPLOID_DOSAGE_NOTE}")
        if r["p_wrong"] is None:
            p_wrong = "n/a (not measured)"
        else:
            ci = r["p_wrong_ci"]
            interval = f" [{ci[0]:.1e}, {ci[1]:.1e}]" if ci is not None else ""
            p_wrong = (
                f"{r['p_wrong']:.1e}{interval}"
                " (measured on GIAB reference samples, DNBSEQ-T7)"
            )
        gq = f", GQ {r['gq']}" if r["gq"] is not None else ""
        flags = ", ".join(r["flags"]) or "none"
        lines.append(
            f"  reliability: {r['level']}, p_wrong {p_wrong}{gq},"
            f" flags: {flags}; {r['row_type']}"
        )
        if r["warnings"]:
            lines.append(f"  warnings: {', '.join(r['warnings'])}")
    reason = f"; {v['reason']}" if v["reason"] else ""
    lines.append(f"  validity: {v['level']} (source: {v['source_level']}){reason}")
    if finding["effect"] is not None:
        lines.append(f"  effect: {json.dumps(finding['effect'])}")
    lines.append(f"  actionable: {'yes' if finding['actionable'] else 'no'}")
    sources = "; ".join(
        f"{s['source']} {s['id']} ({s['version']})" for s in finding["sources"]
    )
    lines.append(f"  sources: {sources}")
    return "\n".join(lines)


def open_store(path: Path) -> duckdb.DuckDBPyConnection:
    """Open an existing store read-only; never create one. Refuse (`RunNotFound`) a store
    another schema version wrote."""
    if not path.is_file():
        raise RunNotFound(f"no findings store at {path}")
    try:
        con = duckdb.connect(str(path), read_only=True)
    except duckdb.Error as e:
        raise RunNotFound(f"cannot open the findings store at {path}: {e}") from e
    try:
        check_schema_version(con, path)
    except (StoreTooNew, StoreTooOld) as e:
        con.close()
        raise RunNotFound(str(e)) from e
    return con


def _snapshot_dates(db: Path, run: str | None) -> dict[str, Dates]:
    """The dates from the snapshots beside the store (`<db's folder>/snapshots`), up to
    the run `run` (None: every one)."""
    from genokeeper.history import dates  # history imports this module

    return dates(Path(db).parent / "snapshots", run)


def _coverage_text(coverage: Mapping[str, Any]) -> str:
    """A gene's coverage in one line: its note first, then its callable share."""
    if coverage["bases"] is None:
        return str(coverage["note"])
    shares = (
        f"{coverage['callable']:.0%} of {coverage['bases']:,} CDS bases callable, mean"
        f" depth {coverage['mean_depth']:.1f}, {coverage['difficult']:.0%} in difficult"
        f" regions, {coverage['outside_benchmark']:.0%} outside the calibrated benchmark"
    )
    note = coverage["note"]
    return shares if note is None else f"{note}; {shares}"


def _query_command(args: argparse.Namespace) -> int:
    filters = Filters(
        topic=args.topic,
        min_reliability=args.min_reliability,
        min_validity=args.min_validity,
        gene=args.gene,
        marker=args.marker,
        text=args.text,
    )
    try:
        with open_store(args.db) as con:
            result = query(
                con,
                args.run,
                filters,
                list_unassessed=args.list_unassessed,
                dates=_snapshot_dates(args.db, args.run),
            )
    except RunNotFound as e:
        print(f"genokeeper query: {e}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(result.to_json(), indent=2, ensure_ascii=False))
        return 0
    print(f"run {result.run_id}: {len(result.findings)} findings")
    for finding in result.findings:
        print()
        print(_text_block(finding))
    print()
    if result.coverage is not None:
        print(f"coverage of {args.gene}: {_coverage_text(result.coverage)}")
    summary = result.unassessed
    print(
        "unassessed (genotype unknown or not covered, not a result):"
        f" {summary['total']} sites under these filters"
    )
    for gene, entry in summary["by_gene"].items():
        topics = ", ".join(f"{t} {n}" for t, n in entry["topics"].items())
        reasons = "; ".join(f"{r} {n}" for r, n in entry["reasons"].items())
        print(f"  {gene or '(no gene)'}: {entry['sites']} ({topics}) - {reasons}")
    for site in result.unassessed_sites or []:
        print()
        print(_text_block(site))
    print()
    if result.hidden:
        print("hidden by the filters:")
        for reason, count in result.hidden.items():
            print(f"  {reason}: {count}")
    else:
        print("hidden by the filters: none")
    return 0


def _show_command(args: argparse.Namespace) -> int:
    try:
        with open_store(args.db) as con:
            run, finding = show(
                con,
                args.run,
                args.statement_id,
                dates=_snapshot_dates(args.db, args.run),
            )
    except RunNotFound as e:
        print(f"genokeeper show: {e}", file=sys.stderr)
        return 2
    if finding is None:
        print(
            f"genokeeper show: no statement {args.statement_id} in this run",
            file=sys.stderr,
        )
        return 1
    if args.json:
        document = {"run_id": run, "finding": finding}
        print(json.dumps(document, indent=2, ensure_ascii=False))
    else:
        print(f"run {run}")
        print(_text_block(finding))
    return 0


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--db",
        type=Path,
        default=load_paths().reports / "genomelab.duckdb",
        help="the findings store (default: reports/genomelab.duckdb)",
    )
    parser.add_argument(
        "--run", metavar="RUN_ID", help="a run id (default: the latest run)"
    )
    parser.add_argument("--json", action="store_true", help="print the JSON contract")


def register(subparsers: Subparsers) -> None:
    """Add the `query` and `show` sub-commands."""
    q = subparsers.add_parser("query", help="list findings, filtered and ordered")
    q.add_argument("--topic", type=Topic, choices=list(Topic))
    q.add_argument(
        "--min-reliability", type=Level, choices=list(Level), default=Level.MEDIUM
    )
    q.add_argument(
        "--min-validity",
        type=Validity,
        choices=list(Validity),
        default=Validity.MODERATE,
    )
    q.add_argument("--gene", help="a gene symbol (case-insensitive)")
    q.add_argument(
        "--marker",
        help="an rsID or a GRCh38 vid (chr1:123:A:G)",
    )
    q.add_argument(
        "--text",
        help="a substring of the gene, condition, classification, id or a drug's name",
    )
    q.add_argument(
        "--list-unassessed",
        action="store_true",
        help="also list the unassessed claims (genotype unknown or not covered) the filters keep",
    )
    _common(q)
    q.set_defaults(handler=_query_command)
    s = subparsers.add_parser("show", help="one statement in full, whatever its kind")
    s.add_argument("statement_id")
    _common(s)
    s.set_defaults(handler=_show_command)
