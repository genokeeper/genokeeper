"""The polygenic scores a run computes (`load_scores`, the packaged `scores.toml`), and the
topic of a trait (`topic_for`, the packaged `efo_topics.tsv`).

Each score's evaluations carry their sample's ancestry in pgsc_calc's population codes
(`POPULATIONS`, the HGDP+1kGP panel's superpopulations), so a run compares them with the
population pgsc_calc finds the user closest to; a mixed sample set is `MIXED`, which
matches none. The control, against a score graded moderate on an evaluation in another
ancestry: its anchor comes from `parse_scores`, which refuses any ancestry outside those
codes, so a catalog category written as it stands ("European") cannot reach the
comparison, where it would never match and leave the grade silently limited.

A score may carry `weak` (slice 3g), and the same file names the conditions no score
qualifies for (`load_not_covered`), each refused if a chosen score also names it.

The same file's `research_score` tables (`load_research_scores`) are scores outside the
chosen set (not reviewed under the selection rule, or not passing it) and are computed anyway:
each has the fields of a score and a required `why`, none repeats a chosen score, and a
condition can hold several.
"""

import functools
import re
import tomllib
from dataclasses import dataclass
from importlib.resources import files
from typing import Any

from genokeeper.store.db import Topic

# The HGDP+1kGP panel's superpopulations, as pgsc_calc names the closest one.
POPULATIONS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")
MIXED = "mixed"
_FIELDS = (
    "pgs_id",
    "condition",
    "trait_efo",
    "publication",
    "development_ancestry",
    "evaluations",
    "licence",
)


@dataclass(frozen=True)
class Score:
    """One chosen PGS Catalog score: the condition it stands for, its trait's EFO and topic,
    its publication, development ancestry and licence, the evaluations (PPM ID,
    population code) the choice rests on, and whether its predictive power is weak."""

    pgs_id: str
    evaluations: tuple[tuple[str, str], ...]
    condition: str
    trait_efo: str
    topic: Topic
    publication: str
    development_ancestry: str
    licence: str
    weak: bool = False

    @property
    def evaluation_ancestry(self) -> set[str]:
        """The population codes it has an evaluation in."""
        return {ancestry for _, ancestry in self.evaluations}


def _resource(name: str) -> str:
    return files("genokeeper").joinpath("resources", "pgs", name).read_text()


@functools.cache
def _efo_topics() -> dict[str, Topic]:
    lines = _resource("efo_topics.tsv").splitlines()
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    return {efo: Topic(topic) for efo, topic, _label in rows[1:]}  # after the header


def topic_for(efo_id: str) -> Topic:
    """The topic of a trait by its EFO (or MONDO) ID; `Topic.RESEARCH` for a trait the
    table does not map (stage 1 design, "Topics and ordering")."""
    return _efo_topics().get(efo_id, Topic.RESEARCH)


def _score(entry: dict[str, Any]) -> Score:
    if set(entry) - {"weak"} != set(_FIELDS):
        raise ValueError(f"score fields {sorted(entry)}: expected {sorted(_FIELDS)}")
    weak = entry.get("weak", False)
    if weak is not True and weak is not False:
        raise ValueError(f"{entry['pgs_id']}: weak {weak!r} is not a boolean")
    evaluations = tuple((str(p), str(a)) for p, a in entry["evaluations"])
    wrong = [a for _, a in evaluations if a not in (*POPULATIONS, MIXED)]
    if not evaluations or wrong:
        raise ValueError(f"{entry['pgs_id']}: evaluation ancestry {wrong}")
    return Score(
        pgs_id=entry["pgs_id"],
        evaluations=evaluations,
        condition=entry["condition"],
        trait_efo=entry["trait_efo"],
        topic=topic_for(entry["trait_efo"]),
        publication=entry["publication"],
        development_ancestry=entry["development_ancestry"],
        licence=entry["licence"],
        weak=weak,
    )


def parse_scores(text: str) -> tuple[Score, ...]:
    """The scores of a `scores.toml` text; `ValueError` for a missing field, a duplicate
    PGS ID or condition, an evaluation outside `POPULATIONS` and `MIXED`, or a `weak` that
    is not a boolean."""
    scores = [_score(entry) for entry in tomllib.loads(text)["score"]]
    for key in ("pgs_id", "condition"):
        values = [getattr(s, key) for s in scores]
        if len(set(values)) != len(values):
            raise ValueError(f"a duplicate {key} among {values}")
    return tuple(scores)


def load_scores() -> tuple[Score, ...]:
    """The packaged score list, in its order."""
    return parse_scores(_resource("scores.toml"))


@dataclass(frozen=True, kw_only=True)
class ResearchScore(Score):
    """A score outside the chosen set, and `why`."""

    why: str


def parse_research_scores(text: str) -> tuple[ResearchScore, ...]:
    """The research scores of a `scores.toml` text; `ValueError` for what a score refuses,
    a `why` missing or blank, a PGS ID repeated or a chosen score's."""
    data = tomllib.loads(text)
    chosen = {e["pgs_id"] for e in data.get("score", [])}
    scores: list[ResearchScore] = []
    for entry in data.get("research_score", []):
        why = entry.get("why")
        if not isinstance(why, str) or not why.strip():
            raise ValueError(f"{entry['pgs_id']}: a research score needs a why")
        base = _score({k: v for k, v in entry.items() if k != "why"})
        if base.pgs_id in chosen:
            raise ValueError(f"{base.pgs_id}: a chosen score, not a research one")
        scores.append(ResearchScore(**vars(base), why=why))
    ids = [s.pgs_id for s in scores]
    if len(set(ids)) != len(ids):
        raise ValueError(f"a duplicate pgs_id among {ids}")
    return tuple(scores)


def load_research_scores() -> tuple[ResearchScore, ...]:
    """The packaged research scores, in their order."""
    return parse_research_scores(_resource("scores.toml"))


@dataclass(frozen=True)
class NotCovered:
    """A condition no PGS Catalog score qualifies for, its topic, and why."""

    condition: str
    topic: Topic
    why: str

    @property
    def slug(self) -> str:
        """The condition in lower case, each run of other characters one hyphen."""
        return re.sub(r"[^a-z0-9]+", "-", self.condition.lower()).strip("-")


def parse_not_covered(text: str) -> tuple[NotCovered, ...]:
    """The not-covered conditions of a `scores.toml` text; `ValueError` for fields other
    than condition, topic and why, an unknown topic, a duplicate condition or slug, or a
    condition whose slug a score in the same text names."""
    data = tomllib.loads(text)
    chosen = {
        NotCovered(e["condition"], Topic.RISKS, "").slug for e in data.get("score", [])
    }
    entries: list[NotCovered] = []
    for entry in data.get("not_covered", []):
        if set(entry) != {"condition", "topic", "why"}:
            raise ValueError(f"not-covered fields {sorted(entry)}")
        entries.append(
            NotCovered(entry["condition"], Topic(entry["topic"]), entry["why"])
        )
        if entries[-1].slug in chosen:
            raise ValueError(f"{entry['condition']}: a chosen score names it")
    slugs = [e.slug for e in entries]
    if len(set(slugs)) != len(slugs):
        raise ValueError(f"a duplicate condition or slug among {slugs}")
    return tuple(entries)


def load_not_covered() -> tuple[NotCovered, ...]:
    """The packaged not-covered conditions, in their order."""
    return parse_not_covered(_resource("scores.toml"))
