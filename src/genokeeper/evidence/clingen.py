"""ClinGen gene-disease validity, from ClinGen's own public CSV export.

`scripts/bootstrap.sh` downloads
https://search.clinicalgenome.org/kb/gene-validity/download into
`<lab>/clingen/<YYYYMMDD>/gene-validity.csv`. The file opens with a preamble
(`"FILE CREATED: YYYY-MM-DD"` among it), then a `+++` rule, the header row (with
`"GENE SYMBOL"`, `"DISEASE ID (MONDO)"` and `"CLASSIFICATION"` among its columns), another
rule and one row per curation: one gene and one disease. Columns are found by name, never by
position.

Slice 1e reads ClinGen's dosage sensitivity curations too (`load_dosage`), which the
bootstrap downloads into `<lab>/clingen-dosage/<YYYYMMDD>/`; slice 3d reads their
triplosensitivity scores and genomic locations as well, and ClinGen's region curations
(`load_regions`) and recurrent CNV regions (`load_recurrent`), which the bootstrap downloads
into `<lab>/clingen-regions/<YYYYMMDD>/`. A genomic location `chrN:start-end` is 1-based
inclusive; it is read as 0-based half-open.
"""

import csv
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

GENE_COLUMN = "GENE SYMBOL"
CLASSIFICATION_COLUMN = "CLASSIFICATION"
DISEASE_COLUMN = "DISEASE ID (MONDO)"
FILE_NAME = "gene-validity.csv"
# The one definition of a positive gene-disease validity, lower-cased; Limited is not one.
POSITIVE_VALIDITY = frozenset({"definitive", "strong", "moderate"})
NEGATIVE_VALIDITY = frozenset({"disputed", "refuted"})

# Several curations of one gene: the strongest positive one wins, then Refuted, then Disputed;
# any other value only when nothing else is known.
_VALIDITY_RANK = {
    "definitive": ("Definitive", 6),
    "strong": ("Strong", 5),
    "moderate": ("Moderate", 4),
    "limited": ("Limited", 3),
    "refuted": ("Refuted", 2),
    "refuted evidence": ("Refuted", 2),
    "disputed": ("Disputed", 1),
    "disputed evidence": ("Disputed", 1),
}
_CREATED = re.compile(r"FILE CREATED:\s*(\d{4})-(\d{2})-(\d{2})")
# ClinGen's dosage sensitivity curations (slice 1e), which the bootstrap fetches whole.
DOSAGE_FILE = "ClinGen_gene_curation_list_GRCh38.tsv"
DOSAGE_GENE = "Gene Symbol"
DOSAGE_SCORE = "Haploinsufficiency Score"
TRIPLO_SCORE = "Triplosensitivity Score"
LOCATION = "Genomic Location"
SUFFICIENT = "3"  # "Sufficient evidence for dosage pathogenicity"
# ClinGen's region curations and recurrent CNV regions (slice 3d), fetched whole together.
REGIONS_FILE = "ClinGen_region_curation_list_GRCh38.tsv"
RECURRENT_FILE = "ClinGen_recurrent_CNV_V2.1-hg38.bed"
REGION_NAME = "ISCA Region Name"
# The BED's colour of a recurrent region; its black (0,0,0) rows are breakpoint clusters.
RECURRENT_RGB = "255,111,0"
_LOCATION = re.compile(r"(chr[^:\s]+):(\d+)-(\d+)")
_DOSAGE_DATE = re.compile(r"#\s*(\d{1,2}) ([A-Z][a-z]{2}),\s*(\d{4})")
# The line ClinGen rewrites at every build of either file: the validity export's request
# day, the dosage list's build day.
_BUILD_DATE = re.compile(rf"\"?FILE CREATED:|{_DOSAGE_DATE.pattern}\s*$")


@dataclass(frozen=True)
class Curation:
    """One gene-disease curation: the disease's MONDO ID (`MONDO:0011450`; empty when the
    export has no MONDO column, which matches no disease) and the classification in
    ClinGen's own words (`Definitive`, `Disputed`, ...). Who curated a gene is
    `evidence.gencc.merge`'s second result."""

    mondo: str
    classification: str


def validity_rank(classification: str) -> int:
    """How strongly a classification ranks among several: Definitive 6 ... Disputed 1, and
    any other value (No Known Disease Relationship) 0."""
    return _VALIDITY_RANK.get(classification.lower(), ("", 0))[1]


def _rows(path: Path) -> Iterator[list[str]]:
    with path.open(newline="", encoding="utf-8") as f:
        yield from csv.reader(f)


def same_curations(a: Path, b: Path) -> bool:
    """Whether two downloads of one ClinGen file hold the same lines but for the date line
    of their build."""

    def kept(path: Path) -> list[str]:
        lines = path.read_text(encoding="utf-8").splitlines()
        return [line for line in lines if not _BUILD_DATE.match(line)]

    return kept(a) == kept(b)


def gene_validity_date(path: Path) -> str:
    """The date the export was created, as YYYYMMDD, from its `FILE CREATED:` line."""
    for row in _rows(path):
        if row and (m := _CREATED.fullmatch(row[0].strip())):
            return "".join(m.groups())
        if GENE_COLUMN in row:
            break
    raise ValueError(f"{path}: no 'FILE CREATED: YYYY-MM-DD' line before the header")


def load_curations(path: Path) -> dict[str, list[Curation]]:
    """Gene symbol -> its curations, in file order, one per disease row.

    Raises ValueError on a file without the header, without the gene or classification
    column, or without a curation. A header without `DISEASE ID (MONDO)` gives curations
    with an empty `mondo`.
    """
    rows = _rows(path)
    header = next((row for row in rows if GENE_COLUMN in row), None)
    if header is None:
        raise ValueError(f"{path}: no header row with a {GENE_COLUMN!r} column")
    if CLASSIFICATION_COLUMN not in header:
        raise ValueError(f"{path}: the header has no {CLASSIFICATION_COLUMN!r} column")
    gene_at, value_at = header.index(GENE_COLUMN), header.index(CLASSIFICATION_COLUMN)
    mondo_at = header.index(DISEASE_COLUMN) if DISEASE_COLUMN in header else None
    width = max(gene_at, value_at, -1 if mondo_at is None else mondo_at)
    curations: dict[str, list[Curation]] = {}
    for row in rows:
        if not row or row[0].startswith("+") or len(row) <= width:
            continue
        gene, value = row[gene_at].strip(), row[value_at].strip()
        if not gene or not value:
            continue
        word = _VALIDITY_RANK.get(value.lower(), (value, 0))[0]
        mondo = "" if mondo_at is None else row[mondo_at].strip()
        curations.setdefault(gene, []).append(Curation(mondo, word))
    if not curations:
        raise ValueError(f"{path}: no curation rows under the header")
    return curations


def gene_level(curations: list[Curation]) -> str:
    """One classification for several curations: the strongest positive one
    (`POSITIVE_VALIDITY`); else Refuted, then Disputed; else the strongest other (Limited,
    then No Known Disease Relationship)."""
    for group in (POSITIVE_VALIDITY, NEGATIVE_VALIDITY):
        if found := [c for c in curations if c.classification.lower() in group]:
            curations = found
            break
    return max(curations, key=lambda c: validity_rank(c.classification)).classification


@dataclass(frozen=True)
class DosageList:
    """ClinGen's dosage sensitivity curations (slice 1e): the date the list was made, as
    YYYYMMDD, the genes its haploinsufficiency score rates 3, sufficient evidence, those its
    triplosensitivity score rates 3 (slice 3d), and each gene's genomic location as
    (chrom, start, end), 0-based half-open (a gene whose location is not `chrN:start-end`,
    as `tbd`, has none)."""

    date: str
    haploinsufficient: frozenset[str]
    triplosensitive: frozenset[str] = frozenset()
    spans: dict[str, tuple[str, int, int]] = field(
        default_factory=dict[str, tuple[str, int, int]]
    )


def _location(text: str) -> tuple[str, int, int] | None:
    """`chrN:start-end` (1-based inclusive) as 0-based half-open, or None."""
    m = _LOCATION.fullmatch(text.strip())
    return None if m is None else (m[1], int(m[2]) - 1, int(m[3]))


def _curation_table(
    path: Path, first: str, needed: Sequence[str]
) -> tuple[str, list[dict[str, str]]]:
    """A ClinGen curation list (genes or regions): `#` lines, one of the first five its
    date (`#23 Oct,2026`), then a `#<first>` header row and one row per curation. Returns
    the date as YYYYMMDD and the rows keyed by column name. Raises ValueError on a file
    without the date, the header, a `needed` column or a row: its anchor comes from these
    column names, which a maintenance page or another file cannot pass because they are
    code, read against the bytes on disk."""
    lines = path.read_text(encoding="utf-8").splitlines()
    date = next((m for line in lines[:5] if (m := _DOSAGE_DATE.fullmatch(line))), None)
    if date is None:
        raise ValueError(f"{path}: no '#DD Mon,YYYY' date line")
    made = datetime.strptime(f"{date[1]} {date[2]} {date[3]}", "%d %b %Y")
    header = next((line for line in lines if line.startswith(f"#{first}")), None)
    if header is None:
        raise ValueError(f"{path}: no '#{first}' header row")
    columns = header.removeprefix("#").split("\t")
    absent = [c for c in needed if c not in columns]
    if absent:
        raise ValueError(f"{path}: the header lacks {absent}")
    width = max(columns.index(c) for c in needed)
    rows = [
        dict(zip(columns, cells, strict=False))
        for line in lines
        if line and not line.startswith("#")
        if len(cells := line.split("\t")) > width
    ]
    if not rows:
        raise ValueError(f"{path}: no rows under the header")
    return made.strftime("%Y%m%d"), rows


def load_dosage(path: Path) -> DosageList:
    """ClinGen's gene curation list (`DOSAGE_FILE`, `_curation_table`'s shape under a
    `#Gene Symbol` header): the genes scored 3 for haploinsufficiency and for
    triplosensitivity, and each gene's location."""
    date, rows = _curation_table(
        path, DOSAGE_GENE, (DOSAGE_GENE, DOSAGE_SCORE, TRIPLO_SCORE, LOCATION)
    )
    spans = {r[DOSAGE_GENE]: _location(r[LOCATION]) for r in rows}
    return DosageList(
        date=date,
        haploinsufficient=frozenset(
            r[DOSAGE_GENE] for r in rows if r[DOSAGE_SCORE].strip() == SUFFICIENT
        ),
        triplosensitive=frozenset(
            r[DOSAGE_GENE] for r in rows if r[TRIPLO_SCORE].strip() == SUFFICIENT
        ),
        spans={gene: span for gene, span in spans.items() if span is not None},
    )


def regions_date(path: Path) -> str:
    """The date ClinGen's region curation list (`REGIONS_FILE`) was made, as YYYYMMDD."""
    return _curation_table(path, "ISCA ID", (REGION_NAME, LOCATION))[0]


def load_regions(path: Path) -> list[tuple[str, int, int, str]]:
    """ClinGen's region curations (`REGIONS_FILE`, `_curation_table`'s shape under an
    `#ISCA ID` header) whose haploinsufficiency or triplosensitivity score is 3, with a
    location: (chrom, start, end, region name), 0-based half-open; what a large event must
    overlap to be shown."""
    _, rows = _curation_table(
        path, "ISCA ID", (REGION_NAME, LOCATION, DOSAGE_SCORE, TRIPLO_SCORE)
    )
    kept: list[tuple[str, int, int, str]] = []
    for r in rows:
        span = _location(r[LOCATION])
        scored = SUFFICIENT in (r[DOSAGE_SCORE].strip(), r[TRIPLO_SCORE].strip())
        if span is not None and scored:
            kept.append((*span, r[REGION_NAME]))
    return kept


def load_recurrent(path: Path) -> list[tuple[str, int, int, str]]:
    """ClinGen's recurrent CNV regions (`RECURRENT_FILE`): the BED's rows coloured
    `RECURRENT_RGB`, as (chrom, start, end, name), 0-based half-open as BED. Raises
    ValueError on a row without nine columns or a file with no recurrent row: its anchor
    comes from that colour and shape, which a maintenance page cannot pass because they are
    checked on the bytes on disk."""
    kept: list[tuple[str, int, int, str]] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cells = line.split("\t")
        if len(cells) < 9 or not (cells[1].isdigit() and cells[2].isdigit()):
            raise ValueError(f"{path}:{number}: not a nine-column BED row")
        if cells[8].strip() == RECURRENT_RGB:
            kept.append((cells[0], int(cells[1]), int(cells[2]), cells[3]))
    if not kept:
        raise ValueError(f"{path}: no recurrent region ({RECURRENT_RGB}) rows")
    return kept


def load_gene_validity(path: Path) -> dict[str, str]:
    """Gene symbol -> ClinGen's gene-level classification (`gene_level`), in ClinGen's own
    words. Raises ValueError as `load_curations` does."""
    return {gene: gene_level(found) for gene, found in load_curations(path).items()}
