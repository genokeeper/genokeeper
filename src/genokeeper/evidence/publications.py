"""Publications about a variant: the PubMed IDs two public files list, read in one streaming
pass each (the refresh keeps only the asked statements' PMIDs, `genokeeper.refresh`).

- LitVar2's whole variants file (`litvar2_variants.json.gz`, NCBI, public domain, about 2 GB):
  one JSON object per line, text-mined mentions. An rsID can span several records (the
  dbSNP record `litvar@rs738409##` and the ClinGen-allele records `litvar@CA...#rs738409##`),
  so a variant's PMIDs are the union over its records. `rsid` is `"rs"` and digits, or null
  for a protein-level or gene-level record, which names no rsID and is skipped.
- ClinVar's `var_citations.txt` (tab-delimited, one row per citation): the papers submitters
  attached to a classification, keyed by VariationID. Only `PubMed` rows are PMIDs; the file
  also has `PubMedCentral` (a PMC number) and `NCBIBookShelf` rows, which are not.

The bootstrap fetches each file whole into `<lab>/litvar/<YYYYMMDD>/` and
`<lab>/clinvar-citations/<YYYYMMDD>/`; nothing here asks a network for a single variant.
Every PMID is a string of digits.
"""

import gzip
import json
import re
from collections.abc import Collection, Iterator
from pathlib import Path

LITVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz"
LITVAR_FILE = "litvar2_variants.json.gz"
CITATIONS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"
)
CITATIONS_FILE = "var_citations.txt"

CHECKED_RECORDS = 1000
_CITATION_COLUMNS = ("VariationID", "citation_source", "citation_id")
_PUBMED = "PubMed"
# The line's rsID, found without parsing the line: the whole file is read for a few thousand
# asked rsIDs, and most of its lines name none. `check_litvar` proves this finds the rsID of
# every checked record whose JSON has one, so a change of the file's spacing is refused there.
_RSID = re.compile(rb'"rsid":\s*"(rs[0-9]+)"')


def _pmid(value: object, where: str) -> str:
    text = str(value)
    if isinstance(value, bool) or not text.isascii() or not text.isdigit():
        raise ValueError(f"{where + ': ' if where else ''}{value!r} is not a PMID")
    return str(int(text))


def _litvar_problem(line: bytes) -> str | None:
    """What is wrong with this line as a LitVar2 record, or None; `rsid` is read as the
    reader finds it, so a record the reader would miss is a problem."""
    try:
        record = json.loads(line)
    except ValueError:
        return "not JSON"
    if not isinstance(record, dict) or not isinstance(record.get("_id"), str):
        return "no _id"
    rsid = record.get("rsid", False)
    if rsid is not None and not (
        isinstance(rsid, str) and re.fullmatch("rs[0-9]+", rsid)
    ):
        return f"rsid is {rsid!r}"
    pmids = record.get("pmids")
    if not isinstance(pmids, list):
        return "no pmids list"
    try:
        for value in pmids:
            _pmid(value, "")
    except ValueError as error:
        return str(error)
    found = _RSID.search(line)
    if (found.group(1).decode() if found else None) != rsid:
        return "the reader would not find its rsid"
    return None


def check_litvar(path: Path) -> None:
    """Raise `ValueError` unless the first 1,000 lines of the gzip file `path` are LitVar2
    records: a JSON object with `_id`, `rsid` (null or `rs` and digits, found by the
    reader's own pattern) and a `pmids` list, and at least one of them names an rsID."""
    named = 0
    seen = 0
    try:
        with gzip.open(path, "rb") as handle:
            for number, line in zip(
                range(1, CHECKED_RECORDS + 1), handle, strict=False
            ):
                problem = _litvar_problem(line)
                if problem:
                    raise ValueError(f"{path.name} line {number}: {problem}")
                seen += 1
                named += _RSID.search(line) is not None
    except (OSError, EOFError) as error:
        raise ValueError(f"{path.name} is not a readable gzip file: {error}") from error
    if not seen:
        raise ValueError(f"{path.name} holds no record")
    if not named:
        raise ValueError(f"{path.name}: none of its first {seen} records names an rsID")


def litvar_pmids(path: Path, rsids: Collection[str]) -> dict[str, frozenset[str]]:
    """The PMIDs LitVar2 lists for each of `rsids` (`rs` and digits), the union of every
    record that names it; one pass over the file. An rsID with no PMID is absent."""
    wanted = frozenset(rsids)
    found: dict[str, set[str]] = {}
    if not wanted:
        return {}
    with gzip.open(path, "rb") as handle:
        for line in handle:
            match = _RSID.search(line)
            if match is None or match.group(1).decode() not in wanted:
                continue
            record = json.loads(line)
            rsid = record["rsid"]
            found.setdefault(rsid, set()).update(
                _pmid(value, f"{path.name} {record['_id']}")
                for value in record["pmids"]
            )
    return {rsid: frozenset(pmids) for rsid, pmids in found.items() if pmids}


def _rows(path: Path, limit: int | None) -> Iterator[tuple[str, str, str]]:
    """Each row's (VariationID, citation_source, citation_id), the columns found by the
    header's names; at most `limit` rows. A header without them, or a short row, raises."""
    with path.open(encoding="utf-8", newline="") as handle:
        names = handle.readline().removeprefix("#").rstrip("\r\n").split("\t")
        missing = [name for name in _CITATION_COLUMNS if name not in names]
        if missing:
            raise ValueError(f"{path.name}: no column {', '.join(missing)}")
        vid, source, cid = (names.index(name) for name in _CITATION_COLUMNS)
        for number, line in enumerate(handle, start=2):
            if limit is not None and number > limit + 1:
                return
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= max(vid, source, cid):
                raise ValueError(f"{path.name} line {number}: too few columns")
            yield fields[vid], fields[source], fields[cid]


def check_citations(path: Path) -> None:
    """Raise `ValueError` unless `path` is ClinVar's `var_citations.txt`: a header naming
    VariationID, citation_source and citation_id, and first 1,000 rows each with a numeric
    VariationID, a source and an id, at least one row."""
    try:
        count = 0
        for number, (vid, source, cid) in enumerate(
            _rows(path, CHECKED_RECORDS), start=2
        ):
            if not vid.isdigit() or not source or not cid:
                raise ValueError(f"{path.name} line {number}: not a citation row")
            if source == _PUBMED:
                _pmid(cid, f"{path.name} line {number}")
            count += 1
    except (OSError, UnicodeDecodeError) as error:
        raise ValueError(f"{path.name} is not a readable text file: {error}") from error
    if not count:
        raise ValueError(f"{path.name} holds no citation row")


def cited_pmids(
    path: Path, variation_ids: Collection[str]
) -> dict[str, frozenset[str]]:
    """The PubMed IDs `var_citations.txt` attaches to each of `variation_ids` (digit
    strings); one pass. PubMedCentral and bookshelf rows are not PMIDs and are skipped. A
    VariationID with no PubMed row is absent."""
    wanted = frozenset(variation_ids)
    found: dict[str, set[str]] = {}
    if not wanted:
        return {}
    for vid, source, cid in _rows(path, None):
        if source == _PUBMED and vid in wanted:
            found.setdefault(vid, set()).add(_pmid(cid, f"{path.name} {vid}"))
    return {vid: frozenset(pmids) for vid, pmids in found.items()}
