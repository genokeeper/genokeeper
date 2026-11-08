"""`genokeeper genotype`: the user's genotype at the sites of a gene or a marker, over the
statements of every kind, each with an answer (track design, A4).

A site's answer comes from its dosage and its reliability alone (the level judged on the
error rate measured on a GIAB sample, `genotype.reliability`): `carrier` (one or two copies of the
claim's allele), `not_a_carrier` (dosage 0 at medium or high reliability),
`probably_not_a_carrier` (dosage 0 at low) or `unknown` (dosage `unknown`, or no call row);
a structural variant (`SV_NOUNS`: slice 3b's deletions, 3d's depth losses and gains and its
insertions) and a PanelApp panel statement `needs_confirmation`, whatever its dosage. A
screening finding under carrier at
validity limited (`possible_carrier`: 3b's deletions, 3d's, 3e's RFC1) answers
`possible_carrier`, never `carrier` (the closing design's "Where each kind of result goes").
A research association counts the copies of its own risk allele, which may be the
reference rather than the ALT the dosage counts: they are recounted from its call
(`genotype.sites.risk_copies`, as `sources.gwas` grades it) and stand for the dosage;
copies that are not countable (another allele beside a reference one) answer unknown.
A marker the store holds no site for is "not assessed" (`assessed` false).

A gene gets no answer of its own, but for one a stage-3 tool (`STAGE3_SOURCES`: the repeat
loci, SMN1, slice 3d's paralog loci and Lp(a)) states at the gene level: that statement leads its sites, answering from its kind
(`GENE_LEVEL_ANSWER`: a finding is only a possible carrier) with its claim (SMN1's copy
numbers, a locus's alleles). Its summary
counts the gene's possible carriers, its deletion, depth and insertion sites and its known
pathogenic sites (the
ClinVar pathogenic class: Pathogenic, Likely pathogenic) by answer and, where the genotype
is known, by reliability, and lists what the run never looked at (`NOT_COVERED`). A gene
whose pathogenic alleles are mostly events SNV calling cannot see names that mechanism
(`dominant_mechanism_not_assessed`, from `resources/mechanism_not_assessed.tsv`), so "none
of the N variants carried" is never read as reassurance, unless the run's own stage-3
statement on the gene is assessed (SMN1's copy number). The control: its anchor comes from
the summary's keys, written in `_gene_summary` from counts and that one assessed-or-not bit,
which no statement's text can move, because the text reaches the summary only as a count and
a stage-3 statement only by its kind.

A gene's summary also carries its coverage in the run (`query.gene_coverage`, slice 3c); a
gene below `[coverage] min_callable` that no stage-3 source states leads its text with
"not reliably assessed: N% callable", so its counts are no reassurance.

The JSON shape (`GenotypeResult.to_json`) is the contract the genome-ask skill reads: each
site is `query`'s finding shape plus `answer`, `classification`, `condition` and
`not_a_finding_because`.
"""

import argparse
import json
import sys
from dataclasses import dataclass
from enum import StrEnum
from importlib.resources import files
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper.deletions.screen import DELETION_SCREEN_LIMITS
from genokeeper.evidence.clinvar import parse_clnsig
from genokeeper.evidence.predicted import NOT_HAPLOINSUFFICIENT
from genokeeper.genotype.sites import risk_copies
from genokeeper.insertions import INSERTIONS_TOOL
from genokeeper.query import (
    RunNotFound,
    _common,
    _finding,
    _marker_matches,
    _run,
    _select,
    _text_block,
    gene_coverage,
    open_store,
    unassessed_reason,
)
from genokeeper.sources.cnv import CNV_TOOL, GAIN_WORDING, LOSS_WORDING, WHOLE_GAIN
from genokeeper.sources.deletions import DELETIONS_TOOL
from genokeeper.sources.panelapp import PANELAPP_TOOL
from genokeeper.sources.paralogs import LPA_TOOL, PARALOGS_TOOL
from genokeeper.sources.predicted import PREDICTED_TOOL
from genokeeper.sources.repeats import CFTR_CALL, REPEATS_TOOL
from genokeeper.sources.smn import SMN_TOOL
from genokeeper.store.db import Kind, Level, Topic, Validity

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

MECHANISMS_FILE = "mechanism_not_assessed.tsv"


class Answer(StrEnum):
    """Whether the user carries a site's allele, from its dosage and reliability."""

    CARRIER = "carrier"
    NOT_A_CARRIER = "not_a_carrier"
    PROBABLY_NOT_A_CARRIER = "probably_not_a_carrier"
    UNKNOWN = "unknown"
    NEEDS_CONFIRMATION = "needs_confirmation"
    POSSIBLE_CARRIER = "possible_carrier"
    INFORMATION = "information"


# What no run of this track looks at, in any gene: a gene's answer always lists it.
NOT_COVERED = (
    (
        "copy-number variants under about 1 kb (most single-exon events), target exons in"
        " segmental duplications without a Parascopy model, and where an extra whole-gene"
        " copy lies (deletions of 50 bp or more, and depth losses and gains of 1 kb or"
        " more, in the screen genes are screened and can miss one; SLC37A4 is not"
        " covered)"
    ),
    (
        "structural variants: inversions and translocations (short-read inversion recall"
        " 0.12-0.17 at 30x)"
    ),
    (
        "repeat expansions outside the 24 loci, RFC1 and CFTR's poly-T/TG of stages 3a"
        " and 3e (FMR1 and AFF2 follow the X dosage; the loci short reads cannot assess"
        " are stated as not covered, with the reason)"
    ),
    "variants in paralogs and pseudogene regions short reads cannot place",
    "variants ClinVar does not classify as pathogenic or likely pathogenic",
)

# CFTR's poly-T/TG statement (slice 3e) is a variable-penetrance allele, never cystic-fibrosis
# carrier status: whatever its kind, it answers neither carrier nor not a carrier, so a
# "no 5T allele" never reads as a negative CF test. The anchor is the statement id, which
# `sources.repeats` writes from `REPEATS_TOOL` and `CFTR_CALL` and no claim text can move.
CFTR_REPEAT_ID = f"{REPEATS_TOOL}:{CFTR_CALL}"
CFTR_REPEAT_LINE = (
    "CFTR poly-T/TG: not cystic-fibrosis carrier status (a variable-penetrance allele"
    " or its absence); see its reason"
)

# The sources whose gene-level statement answers for its gene, and the answer each kind
# gives: a finding is a screening result that needs confirmation, never `carrier` (BR-009:
# validity strong is the gene's link to the disease, not the call, which no calibration
# sample has measured on a positive); a watch (borderline, intermediate, a possible silent
# carrier) is not known either way.
STAGE3_SOURCES = (REPEATS_TOOL, SMN_TOOL, PARALOGS_TOOL, LPA_TOOL)
GENE_LEVEL_ANSWER = {
    Kind.FINDING: Answer.POSSIBLE_CARRIER,
    Kind.NON_FINDING: Answer.NOT_A_CARRIER,
    Kind.WATCH: Answer.UNKNOWN,
    Kind.UNASSESSED: Answer.UNKNOWN,
}

# Slice 3d: the structural-variant sources, each site's noun (a depth gain's is
# `GAIN_WORDING`'s, by its ALT), and the screening sources whose carrier finding at validity
# limited is only a possible carrier (3e's RFC1 calls are `RepeatSource` rows).
SV_NOUNS = {
    DELETIONS_TOOL: "deletion",
    INSERTIONS_TOOL: "insertion",
    CNV_TOOL: LOSS_WORDING.noun,
}
POSSIBLE_SOURCES = (
    DELETIONS_TOOL,
    CNV_TOOL,
    INSERTIONS_TOOL,
    PARALOGS_TOOL,
    REPEATS_TOOL,
)
# A gene-level finding under these topics is information (slice 3d's copy numbers, Lp(a)),
# never carrier status.
INFORMATION_TOPICS = (Topic.TRAITS, Topic.DRUGS, Topic.WELLNESS)

_KNOWN_LEVELS = (Level.HIGH, Level.MEDIUM, Level.LOW)
_ANSWER_TEXT = {
    Answer.CARRIER: "carrier of the allele",
    Answer.NOT_A_CARRIER: "not a carrier of the allele",
    Answer.PROBABLY_NOT_A_CARRIER: "probably not a carrier of the allele",
    Answer.UNKNOWN: "unknown - the genotype could not be determined",
    Answer.NEEDS_CONFIRMATION: (
        "a structural variant that needs confirmation; see its reason"
    ),
    Answer.POSSIBLE_CARRIER: (
        "possibly a carrier: a screening result that needs confirmation; see its reason"
    ),
    Answer.INFORMATION: "information, not a carrier status; see its reason",
}
PANEL_ANSWER_TEXT = "a computer's prediction that needs confirmation; see its reason"


def load_mechanisms() -> dict[str, str]:
    """Gene to the mechanism this track does not assess, from the packaged file."""
    text = files("genokeeper").joinpath("resources", MECHANISMS_FILE).read_text()
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    mechanisms: dict[str, str] = {}
    for number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != 3 or not all(fields):
            raise ValueError(f"{MECHANISMS_FILE}: row {number}: expected 3 columns")
        gene, mechanism, _source = fields
        mechanisms[gene] = mechanism
    return mechanisms


def _dosage(site: dict[str, Any]) -> str | None:
    """The copies of the claim's allele: the call's dosage, or a research association's
    copies of its risk allele on that call (`unknown` when they cannot be known or
    counted)."""
    claim, call = site["claim"] or {}, site["genotype"]
    if "risk_allele" not in claim:
        return call["dosage"]
    vid = site["marker"]["grch38"]
    copies = risk_copies(vid, claim["risk_allele"], call["gt"], call["dosage"])
    return str(copies) if isinstance(copies, int) else "unknown"


def _source(statement: dict[str, Any]) -> str:
    return statement["sources"][0]["source"]


def sv_noun(site: dict[str, Any]) -> str | None:
    """The site's structural-variant noun (`SV_NOUNS`), None for any other source; a depth
    event whose ALT is `<DUP>` is a copy-number gain."""
    if _source(site) == CNV_TOOL and site["marker"]["grch38"].endswith(":<DUP>"):
        return GAIN_WORDING.noun
    return SV_NOUNS.get(_source(site))


def possible_carrier(statement: dict[str, Any]) -> bool:
    """Whether a statement is a doubtful carrier result: a finding under carrier at
    validity limited from one of `POSSIBLE_SOURCES` (a ClinVar or predicted one is not)."""
    return (
        statement["kind"] == Kind.FINDING
        and statement["topic"] == Topic.CARRIER
        and statement["validity"]["level"] == Validity.LIMITED
        and _source(statement) in POSSIBLE_SOURCES
    )


def answer(site: dict[str, Any]) -> Answer:
    """A site's answer from its dosage (`_dosage`) and its reliability level; a possible
    carrier (`possible_carrier`) is one whatever its dosage, and any other structural
    variant (`sv_noun`) or PanelApp panel statement needs confirmation; none is ever
    `carrier`."""
    if possible_carrier(site):
        return Answer.POSSIBLE_CARRIER
    if sv_noun(site) is not None or _source(site) == PANELAPP_TOOL:
        return Answer.NEEDS_CONFIRMATION
    dosage = _dosage(site)
    if dosage is None or dosage == "unknown":
        return Answer.UNKNOWN
    if dosage in ("1", "2"):
        return Answer.CARRIER
    if dosage != "0":
        raise ValueError(f"site {site['id']}: dosage {dosage!r}")
    if site["reliability"]["level"] == Level.LOW:
        return Answer.PROBABLY_NOT_A_CARRIER
    return Answer.NOT_A_CARRIER


def _class(site: dict[str, Any]) -> str | None:
    classification = (site["claim"] or {}).get("classification")
    return parse_clnsig(classification).cls if classification else None


def not_a_finding_because(site: dict[str, Any]) -> str | None:
    """Why a site's statement is not a finding (None for a finding), from its kind, its
    dosage and its ClinVar class; a predicted one (slice 1e) is stored as a non-finding
    only as one copy in a recessive gene (`evidence.predicted`), and so is a structural
    variant (`sv_noun`: slice 3b's `sources.deletions`, 3d's), or as a loss of function in
    a gene neither haploinsufficient nor Tier 3 (its reason names the gate), or as a
    whole-gene gain outside the triplosensitive genes (`WHOLE_GAIN`)."""
    kind = site["kind"]
    if kind == Kind.FINDING:
        return None
    if kind == Kind.UNASSESSED:
        return unassessed_reason(site)
    if kind == Kind.WATCH:
        return "a variant of uncertain significance in a listed gene: watched, not a finding"
    cls = _class(site)
    if _dosage(site) == "0":
        return "the allele is not carried (dosage 0)"
    noun = sv_noun(site)
    if noun is not None:
        reason = site["validity"]["reason"] or ""
        if NOT_HAPLOINSUFFICIENT in reason:
            article = "an" if noun[0] in "aeiou" else "a"
            return f"{article} {noun} is a {NOT_HAPLOINSUFFICIENT}"
        if WHOLE_GAIN in reason:
            return WHOLE_GAIN
        return f"one {noun} in a gene whose condition needs two copies"
    if site["sources"][0]["source"] == PREDICTED_TOOL:
        return (
            "predicted, not classified: one copy in an autosomal recessive gene outside"
            " the carrier list"
        )
    if cls == "benign":
        return "classified benign or likely benign"
    if cls == "vus":
        return "a variant of uncertain significance outside the listed genes"
    if cls in ("pathogenic", "conflicting"):
        return "one copy in an autosomal recessive gene outside the carrier list"
    return "a classification outside disease risk"


def _site(finding: dict[str, Any]) -> dict[str, Any]:
    claim = finding["claim"] or {}
    return {
        **finding,
        "answer": answer(finding),
        "classification": claim.get("classification"),
        "condition": claim.get("condition"),
        "not_a_finding_because": not_a_finding_because(finding),
    }


def _gene_level_answer(finding: dict[str, Any], kind: Kind) -> Answer:
    """A gene-level statement's answer from its kind and topic; a possible carrier
    (`possible_carrier`) is only that, a finding under `INFORMATION_TOPICS` is information,
    and a non-finding at low reliability, or SMN1's copies (a 2+0 carrier is never
    excluded), is only probably not. Any other slice-3e "possible" finding and an
    exclusion-flagged normal (its claim names the test that would confirm it) are screening
    results: unknown, as CFTR's poly-T/TG is."""
    if possible_carrier(finding):
        return Answer.POSSIBLE_CARRIER
    if kind is Kind.FINDING and finding["topic"] in INFORMATION_TOPICS:
        return Answer.INFORMATION
    screening = (
        kind is Kind.FINDING
        and (finding["validity"]["reason"] or "").startswith("possible")
    ) or (kind is Kind.NON_FINDING and (finding["claim"] or {}).get("confirm"))
    if screening:
        return Answer.UNKNOWN
    if kind is Kind.NON_FINDING:
        level = finding["reliability"]["level"] or (finding["gene_call"] or {}).get(
            "level"
        )
        if level == Level.LOW or finding["sources"][0]["source"] == SMN_TOOL:
            return Answer.PROBABLY_NOT_A_CARRIER
    return GENE_LEVEL_ANSWER[kind]


def _gene_level_site(finding: dict[str, Any]) -> dict[str, Any]:
    """A stage-3 gene-level statement as a site: its answer (`_gene_level_answer`), and
    why it is not a finding from its own reason (or its gene call's phenotype)."""
    kind = Kind(finding["kind"])
    if finding["id"] == CFTR_REPEAT_ID:
        return {
            **finding,
            "answer": Answer.UNKNOWN,
            "classification": None,
            "condition": None,
            "not_a_finding_because": CFTR_REPEAT_LINE,
        }
    if kind is Kind.FINDING:
        because = None
    elif kind is Kind.UNASSESSED:
        because = unassessed_reason(finding)
    else:
        phenotype = (finding["gene_call"] or {}).get("phenotype")
        because = finding["validity"]["reason"] or phenotype or str(kind)
    return {
        **finding,
        "answer": _gene_level_answer(finding, kind),
        "classification": None,
        "condition": None,
        "not_a_finding_because": because,
    }


def _position(site: dict[str, Any]) -> tuple[str, str, int, str]:
    chrom, pos, _ref, _alt = site["marker"]["grch38"].split(":")
    return site["gene"] or "", chrom, int(pos), site["id"]


def _gene_summary(
    sites: list[dict[str, Any]],
    gene: str,
    *,
    mechanism_assessed: bool = False,
    coverage: dict[str, Any] | None = None,
) -> dict[str, Any]:
    pathogenic = [s for s in sites if _class(s) == "pathogenic"]
    by_answer = {a.value: 0 for a in Answer}
    by_level = {level.value: 0 for level in _KNOWN_LEVELS}
    for site in pathogenic:
        by_answer[site["answer"]] += 1
        if site["answer"] != Answer.UNKNOWN:
            by_level[site["reliability"]["level"]] += 1

    def count(source: str) -> int:
        return sum(1 for s in sites if _source(s) == source)

    return {
        "pathogenic_sites": len(pathogenic),
        "carried": by_answer[Answer.CARRIER],
        "by_answer": by_answer,
        "by_level": by_level,
        "possible_carriers": sum(
            1 for s in sites if s["answer"] == Answer.POSSIBLE_CARRIER
        ),
        "deletions": count(DELETIONS_TOOL),
        "copy_number": count(CNV_TOOL),
        "insertions": count(INSERTIONS_TOOL),
        "not_covered": list(NOT_COVERED),
        "dominant_mechanism_not_assessed": None
        if mechanism_assessed
        else load_mechanisms().get(gene.upper()),
        "deletion_screen_limit": DELETION_SCREEN_LIMITS.get(gene.upper()),
        "coverage": coverage,
    }


@dataclass(frozen=True)
class GenotypeResult:
    """The selected run, the matching sites in order and, for a gene, its summary."""

    run_id: str
    sites: list[dict[str, Any]]
    gene: dict[str, Any] | None

    @property
    def assessed(self) -> bool:
        return bool(self.sites)

    def to_json(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "assessed": self.assessed,
            "sites": self.sites,
            "gene": self.gene,
        }


def genotype(
    con: duckdb.DuckDBPyConnection,
    run_id: str | None,
    *,
    gene: str | None = None,
    marker: str | None = None,
) -> GenotypeResult:
    """The sites of `gene` (case-insensitive) and/or `marker` (an rsID or a GRCh38 vid, with
    or without its `GRCh38:` prefix and `chr`) in `run_id` (None: the latest run), over the
    statements of every kind, with their answers; with `gene`, its stage-3 gene-level
    statement first, and the gene's summary over all its sites."""
    if gene is None and marker is None:
        raise ValueError("give a gene or a marker")
    run, _versions = _run(con, run_id)
    rows = con.execute(
        f"{_select(list(Kind))} WHERE f.run_id = ? AND f.vid IS NOT NULL", [run]
    ).fetchall()
    sites = sorted(map(_site, map(_finding, rows)), key=_position)
    summary = None
    if gene is not None:
        sites = [s for s in sites if (s["gene"] or "").lower() == gene.lower()]
        listed = ", ".join("?" * len(STAGE3_SOURCES))
        stated = con.execute(
            f"{_select(list(Kind))} WHERE f.run_id = ? AND f.vid IS NULL"
            f" AND f.source IN ({listed}) AND lower(f.gene) = lower(?)"
            " ORDER BY f.statement_id",
            [run, *STAGE3_SOURCES, gene],
        ).fetchall()
        level = [_gene_level_site(_finding(r)) for r in stated]
        assessed = any(s["kind"] != Kind.UNASSESSED for s in level)
        summary = _gene_summary(
            level + sites,
            gene,
            mechanism_assessed=assessed,
            coverage=gene_coverage(con, run, gene),
        )
        sites = level + sites
    if marker is not None:
        sites = [s for s in sites if _marker_matches(s, marker)]
    return GenotypeResult(run_id=run, sites=sites, gene=summary)


# --- command line ---------------------------------------------------------------------------


def _gene_lines(gene: str, summary: dict[str, Any]) -> list[str]:
    """The gene's counts, never a gene-level answer; led by the coverage note when the
    gene is not reliably assessed, so a count of none carried never reads as
    reassurance."""
    answers, levels = summary["by_answer"], summary["by_level"]
    known = sum(levels.values())
    level_counts = ", ".join(f"{level} {n}" for level, n in levels.items())
    coverage = summary.get("coverage") or {}
    unreliable = coverage.get("reliable") is False and coverage.get("owner") is None
    head = (
        (f"{coverage['note']}; " if unreliable else "")
        + f"{gene.upper()}: {summary['pathogenic_sites']} known pathogenic sites in the"
        f" store; the allele carried at {summary['carried']}; genotype known at {known}"
        f" ({level_counts}), unknown at {answers[Answer.UNKNOWN]}"
    )
    lines = [head]
    if summary["dominant_mechanism_not_assessed"] is not None:
        lines.append(
            "  main mechanism not assessed: "
            + summary["dominant_mechanism_not_assessed"]
        )
    if summary["deletion_screen_limit"] is not None:
        lines.append(f"  deletion screen limit: {summary['deletion_screen_limit']}")
    if summary["deletions"]:
        lines.append(
            f"  deletions of 50 bp or more: {summary['deletions']}, each needing"
            " confirmation by MLPA, a microarray, qPCR or long reads"
        )
    for key, what in (
        ("copy_number", "copy-number changes from read depth"),
        ("insertions", "insertions of 50 bp or more"),
    ):
        if summary[key]:
            lines.append(
                f"  {what}: {summary[key]}, each needing confirmation; see each reason"
            )
    if summary["possible_carriers"]:
        lines.append(
            f"  possible carrier at {summary['possible_carriers']}: screening results"
            " that need confirmation"
        )
    lines.append("  not looked at: " + "; ".join(summary["not_covered"]))
    return lines


def _site_block(site: dict[str, Any]) -> str:
    lines = [_text_block(site)]
    if site["id"] == CFTR_REPEAT_ID:
        return "\n".join([*lines, f"  {CFTR_REPEAT_LINE}"])
    lines.append(f"  claim: {site['classification']}; condition: {site['condition']}")
    panel = _source(site) == PANELAPP_TOOL
    text = PANEL_ANSWER_TEXT if panel else _ANSWER_TEXT[Answer(site["answer"])]
    lines.append(f"  answer: {text}")
    if site["not_a_finding_because"] is not None:
        lines.append(f"  not a finding: {site['not_a_finding_because']}")
    return "\n".join(lines)


def _genotype_command(args: argparse.Namespace) -> int:
    if args.gene is None and args.marker is None:
        print("genokeeper genotype: give --gene or --marker", file=sys.stderr)
        return 2
    try:
        with open_store(args.db) as con:
            result = genotype(con, args.run, gene=args.gene, marker=args.marker)
    except RunNotFound as e:
        print(f"genokeeper genotype: {e}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(result.to_json(), indent=2, ensure_ascii=False))
        return 0
    print(f"run {result.run_id}")
    if args.gene is not None and result.gene is not None:
        print("\n".join(_gene_lines(args.gene, result.gene)))
    if not result.assessed:
        asked = args.marker if args.marker is not None else args.gene
        print(f"{asked}: not assessed - the store holds no site for it in this run")
    for site in result.sites:
        print()
        print(_site_block(site))
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `genotype` sub-command."""
    g = subparsers.add_parser(
        "genotype", help="the genotype and an answer at a gene's or a marker's sites"
    )
    g.add_argument("--gene", help="a gene symbol (case-insensitive)")
    g.add_argument(
        "--marker",
        help="an rsID or a GRCh38 vid (chr1:123:A:G)",
    )
    _common(g)
    g.set_defaults(handler=_genotype_command)
