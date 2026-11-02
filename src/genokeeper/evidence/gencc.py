"""GenCC's gene-disease validity submissions, beside ClinGen's curations (slice 3f, "GenCC").

`scripts/bootstrap.sh` (step_gencc) fetches GenCC's export in its new format,
https://thegencc.org/download/action/submissions-export-tsv?format=new, whole into
`<lab>/gencc/<its newest submitted_run_date as YYYYMMDD>/gencc-submissions.tsv`: one tab-
separated row per submission (a submitter's classification of one gene and one disease),
its columns found by name.

A gene ClinGen curates keeps ClinGen's curations only (`merge`); for every other gene its
GenCC submissions are aggregated per disease (`aggregate`): the strongest classification
stands only when no submitter's is Disputed or Refuted, which otherwise caps it at Disputed,
so one lab's Definitive never outweighs another's dispute. GenCC's titles are read in
ClinGen's words (`TITLES`); Supportive and Animal are mapped to Limited and No Known Disease
Relationship, neither of them positive. A title not in `TITLES` refuses the file
(`ValueError`): it could be positive or negative, and no reading of it is safe.
"""

import csv
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

from genokeeper.evidence.clingen import NEGATIVE_VALIDITY, Curation, validity_rank

GENCC_FILE = "gencc-submissions.tsv"
GENCC_SOURCE = "GenCC"
_GENE, _DISEASE, _TITLE = "gene_symbol", "disease_curie", "classification_title"
_RUN_DATE = "submitted_run_date"
# GenCC's classification title (lower-cased) -> ClinGen's word for it.
TITLES = {
    "definitive": "Definitive",
    "strong": "Strong",
    "moderate": "Moderate",
    "limited": "Limited",
    "supportive": "Limited",
    "disputed evidence": "Disputed",
    "refuted evidence": "Refuted",
    "no known disease relationship": "No Known Disease Relationship",
    "animal": "No Known Disease Relationship",
    "animal model only": "No Known Disease Relationship",
}
_DATE = re.compile(r"([0-9]{4})-([0-9]{2})-([0-9]{2})")
_DATED = re.compile(r"[0-9]{8}")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames or []
        missing = [c for c in (_GENE, _DISEASE, _TITLE, _RUN_DATE) if c not in columns]
        if missing:
            raise ValueError(f"{path}: not GenCC's export (no {', '.join(missing)})")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: no submission rows")
    return rows


def load_gencc(path: Path) -> dict[str, list[Curation]]:
    """Gene symbol -> its submissions, in file order, each a `Curation` in ClinGen's
    words. Raises ValueError on a file without the columns read, without a
    submission, or with a title `TITLES` does not know."""
    curations: dict[str, list[Curation]] = {}
    for row in _rows(path):
        gene, title = (row[_GENE] or "").strip(), (row[_TITLE] or "").strip()
        if not gene:
            continue
        word = TITLES.get(title.lower())
        if word is None:
            raise ValueError(f"{path}: unknown GenCC classification {title!r}")
        mondo = (row[_DISEASE] or "").strip()
        curations.setdefault(gene, []).append(Curation(mondo, word))
    return curations


def gencc_date(path: Path) -> str:
    """The export's newest `submitted_run_date`, as YYYYMMDD: the folder it is filed
    under."""
    days = [
        "".join(m.groups())
        for row in _rows(path)
        if (m := _DATE.match((row[_RUN_DATE] or "").strip()))
    ]
    if not days:
        raise ValueError(f"{path}: no {_RUN_DATE}")
    return max(days)


def aggregate(curations: Sequence[Curation]) -> list[Curation]:
    """One curation per disease, in order of first appearance: Disputed when any
    submission is Disputed or Refuted, else the strongest (`clingen.validity_rank`)."""
    by_disease: dict[str, list[Curation]] = {}
    for c in curations:
        by_disease.setdefault(c.mondo, []).append(c)
    found: list[Curation] = []
    for mondo, group in by_disease.items():
        if any(c.classification.lower() in NEGATIVE_VALIDITY for c in group):
            found.append(Curation(mondo, "Disputed"))
        else:
            best = max(group, key=lambda c: validity_rank(c.classification))
            found.append(Curation(mondo, best.classification))
    return found


def merge(
    clingen: Mapping[str, list[Curation]], gencc: Mapping[str, list[Curation]]
) -> tuple[dict[str, list[Curation]], dict[str, str]]:
    """ClinGen's curations, and GenCC's aggregated for every gene ClinGen does not curate;
    with gene -> "GenCC" for each gene taken from GenCC."""
    curations = dict(clingen)
    sources: dict[str, str] = {}
    for gene, submissions in gencc.items():
        if gene not in clingen:
            curations[gene] = aggregate(submissions)
            sources[gene] = GENCC_SOURCE
    return curations, sources


def newest_gencc(lab: Path) -> Path | None:
    """The GenCC export of the newest `<lab>/gencc/<YYYYMMDD>/` holding one, or None."""
    folder = lab / "gencc"
    dated = sorted(
        (d for d in folder.glob("*") if _DATED.fullmatch(d.name)), reverse=True
    )
    return next((d / GENCC_FILE for d in dated if (d / GENCC_FILE).is_file()), None)
