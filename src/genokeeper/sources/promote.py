"""Compound-heterozygote promotion: a pure step over every source's collected rows, with
the phase read-backed phasing found for each candidate pair (`genokeeper.phasing`)."""

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from itertools import combinations

from genokeeper.annotate.vep import locus
from genokeeper.evidence.clinvar import parse_clnsig
from genokeeper.evidence.predicted import FAILED_GATE
from genokeeper.genelists import REPORT_BIALLELIC, GeneList, carrier_genes, gene_mode
from genokeeper.phasing import IN_CIS_REASON, IN_TRANS_REASON, Pair, Phase
from genokeeper.sources import Collected, Fields, Row
from genokeeper.sources.predicted import PREDICTED_TOOL
from genokeeper.store.db import Kind, Topic

COMPOUND_HET_REASON = "two heterozygous variants; phase unknown"
_SOURCES = frozenset({"clinvar", PREDICTED_TOOL})

type _Het = tuple[int, int, object]  # (source, row, vid)


def _hets(
    collected: Sequence[Collected],
    lists: Mapping[str, GeneList],
    inheritance: Mapping[str, str],
    haploid: Callable[[str, int], bool],
) -> dict[str, list[_Het]]:
    """Per recessive gene, its heterozygous pathogenic-class statements, ClinVar's and the
    predicted ones (slice 1e, `sources.predicted`), that their single-het grading left as
    one copy: a carrier finding in a carrier gene (Tier 3 or 4, `carrier_genes`), a
    non-finding in an AR gene outside them,
    or a predicted research statement in an AR gene that failed no LoF gate (design
    decision, 2026-10-28); each judged on the call the store keeps at its vid."""
    carrier_list = carrier_genes(lists)
    calls: dict[object, Fields] = {}
    for c in collected:
        for row in c.rows:
            if row.call is not None:
                calls.setdefault(row.call["vid"], row.call)
    by_gene: dict[str, list[_Het]] = {}
    for i, c in enumerate(collected):
        for j, row in enumerate(c.rows):
            s = row.statement
            gene = str(s.get("gene") or "")
            mode = gene_mode(gene, lists, inheritance.get(gene))
            carrier = s.get("kind") == Kind.FINDING and s.get("topic") == Topic.CARRIER
            single = s.get("kind") == Kind.NON_FINDING and gene not in carrier_list
            source = s.get("source")
            research = (
                source == PREDICTED_TOOL
                and s.get("kind") == Kind.FINDING
                and s.get("topic") == Topic.RESEARCH
                and FAILED_GATE not in str(s.get("reason"))
            )
            if not gene or not (carrier or ((single or research) and mode == "AR")):
                continue
            if source not in _SOURCES:
                continue
            if source == "clinvar":
                claim = json.loads(str(s.get("claim")))
                if parse_clnsig(claim["classification"]).cls != "pathogenic":
                    continue
            vid = s.get("vid")
            call = calls.get(vid) if vid is not None else None
            if call is None or call.get("dosage") != "1":
                continue
            chrom, pos, _, _ = locus(str(vid))
            if not haploid(chrom, pos):
                by_gene.setdefault(gene, []).append((i, j, vid))
    return {
        gene: hets
        for gene, hets in by_gene.items()
        if len({vid for _, _, vid in hets}) >= 2
    }


def _pair(a: object, b: object) -> Pair:
    x, y = sorted((locus(str(a)), locus(str(b))))
    return x, y


def candidate_pairs(
    collected: Sequence[Collected],
    lists: Mapping[str, GeneList],
    inheritance: Mapping[str, str],
    haploid: Callable[[str, int], bool],
) -> list[Pair]:
    """The pairs read-backed phasing is asked about: every two distinct alleles of a gene
    `promote_compound_hets` would promote, as GRCh38 `(chrom, pos, ref, alt)` from their
    vids, each pair in order, sorted."""
    pairs: set[Pair] = set()
    for hets in _hets(collected, lists, inheritance, haploid).values():
        vids = sorted({str(vid) for _, _, vid in hets})
        pairs |= {_pair(a, b) for a, b in combinations(vids, 2)}
    return sorted(pairs)


def promote_compound_hets(
    collected: Sequence[Collected],
    lists: Mapping[str, GeneList],
    inheritance: Mapping[str, str],
    haploid: Callable[[str, int], bool],
    phases: Mapping[Pair, Phase] | None = None,
) -> list[Collected]:
    """`collected` with possible compound heterozygotes promoted to findings; the input is
    left as it is.

    Two or more heterozygous pathogenic-class statements at different variants in one
    recessive gene, ClinVar's or predicted (slice 1e), each graded as one copy (`_hets`: a
    carrier finding in a carrier gene, a non-finding, or predicted research), are judged by their
    pairs' `phases` (`candidate_pairs`; a pair without one is unphased): a statement in
    trans with any partner becomes a finding under risks, "in trans (read-backed)"; one in
    cis with every partner keeps its single-het kind and topic, "in cis (read-backed)", as
    two variants on one haplotype still make one copy; any other becomes a finding under
    risks, "two heterozygous variants; phase unknown". Without `phases` every pair is
    unphased. A promoted ClinVar statement in an ACMG SF gene whose report rule is
    `biallelic` becomes actionable: SF v3.3 reports two P/LP variants "even when phase is
    undetermined" (Table 1 footnote a); a predicted one never does. Each statement is judged on the call the store keeps at its vid, the first
    one written (`write_collected`), whichever source collected the statement. A copy on a haploid locus (`haploid(chrom,
    pos)`, at the vid's GRCh38 locus) is no het and never enters the promotion; nor does
    an unknown genotype (only a dosage of 1 counts), so an unassessed statement never
    counts toward a pair and is never promoted.
    """
    known = phases or {}
    outcome: dict[tuple[int, int], tuple[bool, str, bool]] = {}
    for gene, hets in _hets(collected, lists, inheritance, haploid).items():
        sf = lists["acmg_sf"].report.get(gene) == REPORT_BIALLELIC
        for i, j, vid in hets:
            partners = {str(other) for _, _, other in hets if other != vid}
            found = {known.get(_pair(vid, p), Phase.UNPHASED) for p in partners}
            reported = sf and collected[i].rows[j].statement.get("source") == "clinvar"
            if Phase.IN_TRANS in found:
                outcome[(i, j)] = (True, IN_TRANS_REASON, reported)
            elif found == {Phase.IN_CIS}:
                outcome[(i, j)] = (False, IN_CIS_REASON, False)
            else:
                outcome[(i, j)] = (True, COMPOUND_HET_REASON, reported)
    return [
        replace(
            c,
            rows=[
                _judged(row, *outcome[(i, j)]) if (i, j) in outcome else row
                for j, row in enumerate(c.rows)
            ],
        )
        for i, c in enumerate(collected)
    ]


def _judged(row: Row, finding: bool, why: str, reported: bool) -> Row:
    reason = "; ".join(str(r) for r in (row.statement.get("reason"), why) if r)
    if not finding:
        return replace(row, statement={**row.statement, "reason": reason})
    return replace(
        row,
        statement={
            **row.statement,
            "kind": str(Kind.FINDING),
            "topic": str(Topic.RISKS),
            "reason": reason,
            "actionable": row.statement.get("actionable") is True or reported,
        },
    )
