"""MedlinePlus Genetics summaries: what a gene does and what a condition is, in the National
Library of Medicine's own plain-language text (public domain), which grounds the generated
explanations (spec item 11, `genokeeper.explain`).

The bootstrap fetches the whole `ghr-summaries.xml` into `<lab>/medlineplus/<fetch date>/`.
A gene is found by its exact symbol; a condition by its name or a synonym after one
normalisation (case, `_` as a space, runs of whitespace), so ClinVar's
`Hereditary_hemochromatosis` finds `Hereditary hemochromatosis` and nothing fuzzier does. A
normalised name that two conditions share names neither.
"""

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

MEDLINEPLUS_URL = "https://medlineplus.gov/download/ghr-summaries.xml"
MEDLINEPLUS_FILE = "ghr-summaries.xml"


@dataclass(frozen=True)
class Summary:
    """One summary: `id` is `medlineplus:gene:<SYMBOL>` or `medlineplus:condition:<page
    name>`; `text` the gene's function or the condition's description, its paragraphs
    joined by a blank line; `published` the summary's own date, its source version."""

    id: str
    name: str
    text: str
    published: str


def _normalised(name: str) -> str:
    return " ".join(name.replace("_", " ").split()).casefold()


@dataclass(frozen=True)
class MedlinePlus:
    """The gene summaries by symbol, the condition summaries by every normalised name and
    synonym."""

    genes: dict[str, Summary]
    conditions: dict[str, Summary]

    def condition(self, name: str) -> Summary | None:
        return self.conditions.get(_normalised(name))


def _tag(element: ET.Element) -> str:
    """The element's name without its namespace."""
    return element.tag.rpartition("}")[2]


def _child(element: ET.Element, name: str) -> str:
    """The text of `element`'s first direct child called `name`, or ""."""
    return next((c.text or "" for c in element if _tag(c) == name), "").strip()


def _text(summary: ET.Element, role: str) -> str:
    """The summary's text of `role`: its paragraphs with the inline markup dropped."""
    for text in summary.iter():
        if _tag(text) == "text" and _child(text, "text-role") == role:
            paragraphs = (
                "".join(p.itertext()).strip() for p in text.iter() if _tag(p) == "p"
            )
            return "\n\n".join(paragraphs)
    return ""


def _summary(element: ET.Element) -> Summary:
    """A `gene-summary` or `health-condition-summary` element as a `Summary`."""
    if _tag(element) == "gene-summary":
        id, role = f"medlineplus:gene:{_child(element, 'gene-symbol')}", "function"
    else:
        page = _child(element, "ghr-page").rpartition("/")[2]
        id, role = f"medlineplus:condition:{page}", "description"
    return Summary(
        id=id,
        name=_child(element, "name"),
        text=_text(element, role),
        published=_child(element, "published"),
    )


def load_medlineplus(path: Path) -> MedlinePlus:
    """The summaries in the `ghr-summaries.xml` at `path`; raises `ValueError` for a file
    whose root is not `summaries` or that holds no gene summary."""
    genes: dict[str, Summary] = {}
    conditions: dict[str, Summary] = {}
    shared: set[str] = set()
    try:
        events = ET.iterparse(path, events=("start", "end"))
        _, root = next(events)
        if _tag(root) != "summaries":
            raise ValueError(f"{path}: not MedlinePlus Genetics summaries ({root.tag})")
        for event, element in events:
            kind = _tag(element)
            if event != "end" or kind not in (
                "gene-summary",
                "health-condition-summary",
            ):
                continue
            summary = _summary(element)
            if kind == "gene-summary":
                genes[_child(element, "gene-symbol")] = summary
            else:
                names = [summary.name]
                for synonyms in (c for c in element if _tag(c) == "synonym-list"):
                    names.extend(s.text or "" for s in synonyms)
                for key in {_normalised(n) for n in names} - {""}:
                    if key in conditions:
                        shared.add(key)
                    conditions[key] = summary
            # Each summary is read once: drop it, so the whole file is never held.
            root.clear()
    except (ET.ParseError, OSError) as e:
        raise ValueError(
            f"{path}: not a readable MedlinePlus Genetics file ({e})"
        ) from e
    if not genes:
        raise ValueError(f"{path}: holds no MedlinePlus Genetics gene summary")
    for key in shared:
        del conditions[key]
    return MedlinePlus(genes=genes, conditions=conditions)
