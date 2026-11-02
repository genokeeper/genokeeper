"""Eleven traits as a source (`TraitSource`): slice 3f's five ("Traits") and slice 3g's six
single-gene traits. Tool `traits`, one statement per trait, not per rsid (`trait:<name>` for
each of `load_traits()`), under traits at validity moderate, never actionable, never flagged
`exclusion`, on our genotypes at `resources/traits/traits.tsv`'s alleles (`TRAIT_ALLELES`,
graded as every genotype is). Each site's effect-allele copies come from
`genotype.sites.risk_copies`, so an effect allele that is the REF (ADH1B's T) counts the
reference copies.

- `any` (lactase persistence): "likely" when any enhancer allele is carried, "unlikely" only
  when every one is called not carried, unassessed otherwise.
- `tas2r38`: PAV copies at rs713598, rs1726866 and rs10246939, assuming only PAV and AVI:
  all three 2, 1 or 0 give PAV/PAV, PAV/AVI or AVI/AVI; any other combination is a rare
  haplotype (`RARE_TAS2R38`), unassessed.
- `single`: the effect allele's copies at one site.
- `mc1r`: R, the effect copies over the `strong` sites, and r over the `weak` ones; R of two
  or more, else (every strong site known) R of one, r alone (`MC1R_WEAK_ONLY`) or none.

A known result is a finding of information: its wording, the trait's note, and
`INFORMATION_NOTE` (the phenotype is the literature's, not checked here); an unknown or
uncountable genotype leaves it unassessed. A statement is gene-level (ruling R24): each
site's genotype and level are in its claim, and its `reliability` is the lowest level over
the sites its result reads.
"""

import csv
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Literal, cast

from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.genotype.sites import NOT_COUNTABLE, risk_copies
from genokeeper.sources import INFORMATION_NOTE, Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Level, Topic, Validity
from genokeeper.variants.key import make_vid

TRAITS_TOOL = "traits"
TRAITS_FILE = (
    Path(__file__).resolve().parents[1] / "resources" / "traits" / "traits.tsv"
)
TRAITS_VERSION = "2026-10-30"
RARE_TAS2R38 = "rare TAS2R38 haplotype"
MC1R_WEAK_ONLY = "weak variants only: little effect expected"
STRONG, WEAK = "strong", "weak"
LIKELY = "likely"
UNLIKELY = "unlikely"
_DIPLOTYPES = {0: "AVI/AVI", 1: "PAV/AVI", 2: "PAV/PAV"}
_RANK = (Level.LOW, Level.MEDIUM, Level.HIGH)

type Rule = Literal["any", "single", "tas2r38", "mc1r"]


@dataclass(frozen=True)
class Trait:
    """One trait: its name, gene, sites (alleles, their rsids and effect alleles), rule,
    wordings by effect-allele copies (or PAV haplotypes) 0, 1 and 2, note and sources;
    `classes` is each site's MC1R class (`strong` or `weak`), empty for every other rule."""

    name: str
    gene: str
    alleles: tuple[Allele, ...]
    rsids: tuple[str, ...]
    effects: tuple[str, ...]
    classes: tuple[str, ...]
    rule: Rule
    wordings: tuple[str, str, str]
    note: str | None
    source: str


@cache
def load_traits(path: Path = TRAITS_FILE) -> tuple[Trait, ...]:
    """The packaged traits, in the file's order; a trait whose rows disagree on a
    trait-level column, or an `mc1r` site without a class (any other site with one), is
    refused."""
    with path.open(newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    rows: dict[str, list[dict[str, str]]] = {}
    for row in csv.DictReader(lines, delimiter="\t"):
        rows.setdefault(row["trait"], []).append(row)
    traits: list[Trait] = []
    shared = ("gene", "rule", "none", "one", "two", "note", "source")
    for name, sites in rows.items():
        first = sites[0]
        if any(site[c] != first[c] for site in sites for c in shared):
            raise ValueError(f"{path}: {name}'s rows disagree")
        if first["rule"] not in ("any", "single", "tas2r38", "mc1r"):
            raise ValueError(f"{path}: {name}'s rule {first['rule']!r}")
        classes = tuple(s["class"] for s in sites)
        allowed = (STRONG, WEAK) if first["rule"] == "mc1r" else ("",)
        if any(c not in allowed for c in classes):
            raise ValueError(f"{path}: {name}'s class {classes!r}")
        traits.append(
            Trait(
                name=name,
                gene=first["gene"],
                alleles=tuple(
                    (s["chrom"], int(s["pos"]), s["ref"], s["alt"]) for s in sites
                ),
                rsids=tuple(s["rsid"] for s in sites),
                effects=tuple(s["effect"] for s in sites),
                classes=classes,
                rule=cast(Rule, first["rule"]),
                wordings=(first["none"], first["one"], first["two"]),
                note=first["note"] or None,
                source=first["source"],
            )
        )
    return tuple(traits)


# Every allele a trait reads, for the run's genotypes step.
TRAIT_ALLELES: tuple[Allele, ...] = tuple(a for t in load_traits() for a in t.alleles)


def _copies(
    allele: Allele, effect: str, graded: Mapping[Allele, tuple[Called, Reliability]]
) -> int | Literal["not_countable"] | None:
    got = graded.get(allele)
    if got is None:
        return None
    called = got[0]
    return risk_copies(make_vid(*allele), effect, called.gt, called.dosage)


def _why(rsids: Sequence[str], copies: Sequence[int | str | None]) -> str:
    """Why these sites leave a result unknown ("" when every one is counted)."""
    unknown = [r for r, c in zip(rsids, copies, strict=True) if c is None]
    other = [r for r, c in zip(rsids, copies, strict=True) if c == NOT_COUNTABLE]
    return "; ".join(
        part
        for part in (
            f"genotype unknown at {', '.join(unknown)}" if unknown else "",
            f"another allele at {', '.join(other)}" if other else "",
        )
        if part
    )


def _mc1r(
    trait: Trait, copies: Sequence[int | str | None]
) -> tuple[str | None, str | None, str | None]:
    """MC1R's result from R, the effect copies at its strong sites, and r at its weak ones:
    R of two or more whatever else is unknown; then, every strong site known, R of one, r
    alone, or (every weak site known too) none. Two strong copies at different sites count
    as two: strong alleles are not seen on one haplotype."""
    sites = {
        k: (
            [r for r, c in zip(trait.rsids, trait.classes, strict=True) if c == k],
            [n for n, c in zip(copies, trait.classes, strict=True) if c == k],
        )
        for k in (STRONG, WEAK)
    }
    strong, weak = (
        sum(n for n in sites[k][1] if isinstance(n, int)) for k in (STRONG, WEAK)
    )
    none, one, two = trait.wordings
    if strong >= 2:
        return "two strong", two, None
    if why := _why(*sites[STRONG]):
        return None, None, why
    if strong == 1:
        return "one strong", one, None
    if weak > 0:
        return "weak only", MC1R_WEAK_ONLY, None
    if why := _why(*sites[WEAK]):
        return None, None, why
    return "none", none, None


def _result(
    trait: Trait, copies: Sequence[int | str | None]
) -> tuple[int | str | None, str | None, str | None]:
    """The trait's result for the claim, its wording, and why it is unknown."""
    if trait.rule == "mc1r":
        return _mc1r(trait, copies)
    counted = [c for c in copies if isinstance(c, int)]
    why = _why(trait.rsids, copies)
    if trait.rule == "any":
        if any(c > 0 for c in counted):
            return LIKELY, trait.wordings[1], None
        return (None, None, why) if why else (UNLIKELY, trait.wordings[0], None)
    if why:
        return None, None, why
    if trait.rule == "single":
        return counted[0], trait.wordings[counted[0]], None
    if len(set(counted)) != 1:
        return None, None, RARE_TAS2R38
    return _DIPLOTYPES[counted[0]], trait.wordings[counted[0]], None


def _reads(
    trait: Trait, result: object, copies: Sequence[int | str | None]
) -> list[bool]:
    """The sites a result reads: "likely" the carried ones, "two strong" the carried strong
    ones, "one strong" every strong one, every other result them all."""
    carried = [isinstance(c, int) and c > 0 for c in copies]
    if result == LIKELY:
        return carried
    if result == "two strong":
        return [k and c == STRONG for k, c in zip(carried, trait.classes, strict=True)]
    if result == "one strong":
        return [c == STRONG for c in trait.classes]
    return [True] * len(copies)


@dataclass(frozen=True)
class TraitSource:
    """Our graded genotypes at `TRAIT_ALLELES` (an allele missing is unknown)."""

    graded: Mapping[Allele, tuple[Called, Reliability]]
    tool: str = TRAITS_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One gene-level statement per trait (ruling R24)."""
        rows: list[Row] = []
        for trait in load_traits():
            got = [self.graded.get(a) for a in trait.alleles]
            copies = [
                _copies(a, e, self.graded)
                for a, e in zip(trait.alleles, trait.effects, strict=True)
            ]
            result, wording, why = _result(trait, copies)
            kind = Kind.UNASSESSED if wording is None else Kind.FINDING
            reason, reliability = why, None
            if wording is not None:
                note = [trait.note] if trait.note else []
                reason = "; ".join([wording, *note, INFORMATION_NOTE])
                read = [
                    g[1].level
                    for g, r in zip(got, _reads(trait, result, copies), strict=True)
                    if g is not None and r
                ]
                reliability = str(min(read, key=_RANK.index))
            claim = {
                "schema": CLAIM_SCHEMA,
                "type": "trait",
                "trait": trait.name,
                "gene": trait.gene,
                "result": result,
                "reliability": reliability,
                "sites": [
                    {
                        "rsid": r,
                        "vid": make_vid(*a),
                        "effect_allele": e,
                        "copies": c,
                        "genotype": None if g is None else g[0].gt,
                        "level": None if g is None else str(g[1].level),
                    }
                    for r, a, e, c, g in zip(
                        trait.rsids,
                        trait.alleles,
                        trait.effects,
                        copies,
                        got,
                        strict=True,
                    )
                ],
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"trait:{trait.name}",
                "kind": str(kind),
                "vid": None,
                "gene": trait.gene,
                "claim": json.dumps(claim, ensure_ascii=False),
                "validity": str(Validity.MODERATE),
                "source_level": trait.source,
                "reason": reason,
                "effect": None,
                "actionable": False,
                "topic": str(Topic.TRAITS),
                "source": TRAITS_TOOL,
                "source_id": trait.name,
                "source_version": TRAITS_VERSION,
                "points": None,
                "exclusion": False,
            }
            rows.append(Row(None, None, statement))
        return Collected(
            tool=self.tool, rows=rows, extra={}, input_rows=len(rows), dropped={}
        )
