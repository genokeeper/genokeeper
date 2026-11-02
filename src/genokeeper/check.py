"""`genokeeper check-updates`: ask each source's newest version at a fixed URL and say, per
source, whether it is newer than the one on this machine (spec, "The daily check").

A request is for a version probe only (a checksum file, a header, a listing, a small JSON
document) at a URL written in `PROBES`; no request depends on the user's data or on what the
lab holds. An upstream string is untrusted: each probe reads one shape out of the response and
the result must `fullmatch` the probe's pattern, so a hostile tag, a captive portal's page or a
garbled header is "could not check" (`failed`), never `same`. The result is
`<lab>/refresh/check.json`; a notification text built from it carries source names, versions
and days only and reaches the notifier as an argument (`genokeeper.schedule.notify`).
"""

import argparse
import csv
import fcntl
import hashlib
import io
import json
import os
import re
import sys
import urllib.error
import urllib.request
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

from genokeeper import refresh_log
from genokeeper.config import load_paths
from genokeeper.evidence.clingen import DOSAGE_FILE, FILE_NAME
from genokeeper.evidence.gencc import GENCC_FILE
from genokeeper.lang import Language, language
from genokeeper.refresh_log import SourceStatus

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

CHECK = "check.json"  # under <lab>/refresh/
CSPEC_FILE = "cspec_versions.tsv"  # gene, SVI id and version, beside the BS1 thresholds
TIMEOUT = 30  # seconds per request
MAX_BODY = (
    2 * 1024 * 1024
)  # a probe reads at most this much: the big files are never fetched
VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
# ClinGen's two files carry their version in their rows (about 1 MB and 0.25 MB), so they
# are read whole, up to `WHOLE_MAX`; the cap above is for the multi-GB files.
WHOLE = frozenset({VALIDITY_URL, DOSAGE_URL})
WHOLE_MAX = 64 * 1024 * 1024
# GenCC builds its export per request and serves 20 a day; a HEAD is exempt from that limit,
# so its probe sends a HEAD (`head`), never a GET, and compares the validators it answers.
GENCC_URL = "https://thegencc.org/download/action/submissions-export-tsv?format=new"
HEADS = frozenset({GENCC_URL})
NOT_DAILY = (
    "not checked daily: GenCC gives no change date; genokeeper refresh fetches it"
)
REFRESH_DAYS = (
    30  # a refreshed source's newer release is news once the refresh is this old
)
USER_AGENT = "genokeeper-check-updates (version probes only)"

type Http = Callable[[str], tuple[bytes, Mapping[str, str]]]
type Kind = Literal["refreshed", "pinned", "cspec"]
type Status = Literal["newer", "same", "failed"]

_DATE = r"[0-9]{8}"
# At most 9 digits a part: a longer one is no version, and could not be read as a number.
_SEMVER = r"v?[0-9]{1,9}(\.[0-9]{1,9}){0,3}"
_CSPEC_VERSION = r"[0-9]+(\.[0-9]+){0,2}"  # CSpec's own "1.0.0" included
_CSPEC_LIST = rf"GN[0-9]{{3}}={_CSPEC_VERSION}(,GN[0-9]{{3}}={_CSPEC_VERSION})*"
# A validator a HEAD answers, named: GenCC's ETag (its file's md5), else the Last-Modified
# day, else the length. The local version lists all three of the newest download.
_VALIDATOR = r"(etag:[0-9a-f]{32}|last-modified:[0-9]{8}|content-length:[0-9]{1,15})"
_VALIDATORS = rf"{_VALIDATOR}(,{_VALIDATOR}){{0,2}}"


@dataclass(frozen=True)
class Probe:
    """One source: where its newest version is asked, the shape it must have, how the
    response yields it, how the version on this machine is read (the same scheme), and
    whether the upstream one is newer (`compare`; by default later by part, `newer_than`)."""

    name: str
    kind: Kind
    url: str
    pattern: re.Pattern[str]
    upstream: Callable[[bytes, Mapping[str, str]], str]
    local: Callable[[Path], str | None]
    compare: Callable[[str, str], bool] | None = None


@dataclass(frozen=True)
class Found:
    """What one probe found: the versions and `newer`, `same` or `failed` (`error` is then
    a class or a code, never text)."""

    name: str
    kind: str
    upstream: str | None
    local: str | None
    status: Status
    error: str | None = None


class ShapeError(ValueError):
    """A response that does not carry the version shape its probe reads."""


class NotCheckedDaily(ShapeError):
    """A HEAD that answers none of the validators its probe compares (`NOT_DAILY`)."""


# --- reading a response ---------------------------------------------------------------------


def _day(moment: datetime) -> str:
    return moment.astimezone(UTC).strftime("%Y%m%d")


def _last_modified(body: bytes, headers: Mapping[str, str]) -> str:
    """The `Last-Modified` date (UTC, YYYYMMDD) of a response typed as a file, not as a
    web page, whose body does not open as one either."""
    kind = headers.get("content-type", "").lower()
    page = body[:512].lstrip().lower().startswith((b"<!doctype html", b"<html"))
    if not kind or kind.startswith("text/html") or page:
        raise ShapeError("not the file")
    return _day(parsedate_to_datetime(headers["last-modified"]))


def _md5_name(body: bytes, headers: Mapping[str, str]) -> str:
    """The date in the file name ClinVar's `.md5` names (`clinvar_YYYYMMDD.vcf.gz`)."""
    del headers
    found = re.fullmatch(
        r"[0-9a-f]{32}\s+\S*clinvar_([0-9]{8})\.vcf\.gz\s*", body.decode("ascii")
    )
    if found is None:
        raise ShapeError("not ClinVar's .md5")
    return found[1]


def _json(body: bytes) -> dict[str, Any]:
    document = json.loads(body)
    if not isinstance(document, dict):
        raise ShapeError("not a JSON object")
    return document  # pyright: ignore[reportUnknownVariableType]


def _release_date(body: bytes, headers: Mapping[str, str]) -> str:
    """A GitHub release's `published_at` day (CPIC's data releases)."""
    del headers
    published = str(_json(body).get("published_at"))
    found = re.fullmatch(
        r"([0-9]{4})-([0-9]{2})-([0-9]{2})T[0-9]{2}:[0-9]{2}:[0-9]{2}Z", published
    )
    if found is None:
        raise ShapeError("no published_at")
    return "".join(found.groups())


def _tag(tag: str) -> Callable[[bytes, Mapping[str, str]], str]:
    """A GitHub `releases/latest` tag: `tag` is its whole shape with the version as its one
    group (`v([0-9]+[.][0-9]+[.][0-9]+)`), so nothing else survives."""

    def upstream(body: bytes, headers: Mapping[str, str]) -> str:
        del headers
        release = _json(body)
        found = re.fullmatch(tag, str(release.get("tag_name")))
        if found is None or release.get("draft") or release.get("prerelease"):
            raise ShapeError("an unexpected release tag")
        return found[1]

    return upstream


def _pgs_release(body: bytes, headers: Mapping[str, str]) -> str:
    del headers
    latest = _json(body).get("latest_release")
    found = re.fullmatch(
        r"([0-9]{4})-([0-9]{2})-([0-9]{2})",
        str(latest.get("date") if isinstance(latest, dict) else None),  # pyright: ignore[reportUnknownMemberType]
    )
    if found is None:
        raise ShapeError("no latest_release.date")
    return "".join(found.groups())


def _gwas_release(body: bytes, headers: Mapping[str, str]) -> str:
    """The date of the associations zip in the `latest/` listing: its release's build date."""
    del headers
    found = re.findall(
        r'href="gwas-catalog-associations_ontology-annotated-full\.zip">[^<]*</a></td>'
        r'<td align="right">([0-9]{4})-([0-9]{2})-([0-9]{2}) ',
        body.decode("utf-8", "replace"),
    )
    if len(found) != 1:
        raise ShapeError("not the GWAS Catalog listing")
    return "".join(found[0])


def _gnomad_release(body: bytes, headers: Mapping[str, str]) -> str:
    """The newest `release/<major>.<minor>[.<patch>]/` prefix of the bucket, as
    `<major>.<minor>` (a patch release is not a new data release)."""
    del headers
    prefixes = _json(body).get("prefixes")
    if not isinstance(prefixes, list):
        raise ShapeError("no prefixes")
    found = [
        (int(m[1]), int(m[2]))
        for p in prefixes  # pyright: ignore[reportUnknownVariableType]
        if (m := re.fullmatch(r"release/([0-9]+)\.([0-9]+)(?:\.[0-9]+)?/", str(p)))
    ]
    if not found:
        raise ShapeError("no release prefix")
    return ".".join(map(str, max(found)))


def _newest_day(rows: list[list[str]], header: list[str] | None, column: str) -> str:
    """The newest `YYYY-MM-DD` date (as YYYYMMDD) in `column` of `rows`."""
    if header is None or column not in header:
        raise ShapeError(f"no {column} column")
    at = header.index(column)
    days = [
        found[1] + found[2] + found[3]
        for r in rows
        if len(r) > at
        and (found := re.match(r"([0-9]{4})-([0-9]{2})-([0-9]{2})", r[at]))
    ]
    if not days:
        raise ShapeError(f"no {column}")
    return max(days)


def validity_date(text: str) -> str:
    """The newest `CLASSIFICATION DATE` (YYYYMMDD) of a ClinGen gene-disease validity
    export. Its `FILE CREATED` line is the day of the request, so it says nothing of the
    curations; the newest classification does."""
    rows = list(csv.reader(io.StringIO(text)))
    header = next((r for r in rows if "GENE SYMBOL" in r), None)
    return _newest_day(rows, header, "CLASSIFICATION DATE")


def dosage_date(text: str) -> str:
    """The newest `Date Last Evaluated` (YYYYMMDD) of ClinGen's dosage sensitivity list,
    whose own date line is the day of the nightly build."""
    rows = [line.removeprefix("#").split("\t") for line in text.splitlines()]
    header = next((r for r in rows if r[0] == "Gene Symbol"), None)
    return _newest_day(rows, header, "Date Last Evaluated")


def _whole(read: Callable[[str], str]) -> Callable[[bytes, Mapping[str, str]], str]:
    """The version `read` finds in a ClinGen file read whole (one cut at `WHOLE_MAX` is no
    whole file)."""

    def upstream(body: bytes, headers: Mapping[str, str]) -> str:
        if len(body) >= WHOLE_MAX:
            raise ShapeError("longer than a probe reads")
        return read(body.decode("utf-8"))

    return upstream


def cspec_rows() -> list[tuple[str, str, str]]:
    """(gene, SVI id, version) of every specification recorded beside a BS1 threshold."""
    text = files("genokeeper").joinpath("resources", CSPEC_FILE).read_text("utf-8")
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    rows = [tuple(line.split("\t")) for line in lines[1:]]
    return [(gene, svi, version) for gene, svi, version in rows]


def _validator(body: bytes, headers: Mapping[str, str]) -> str:
    """The first validator a HEAD of a file (not a web page) answers: its ETag, else its
    Last-Modified day, else its length; none of them is `NotCheckedDaily`."""
    del body
    kind = headers.get("content-type", "").lower()
    if not kind or kind.startswith("text/html"):
        raise ShapeError("not the file")
    if etag := headers.get("etag"):
        return "etag:" + etag.removeprefix("W/").strip('"')
    if modified := headers.get("last-modified"):
        return f"last-modified:{_day(parsedate_to_datetime(modified))}"
    if length := headers.get("content-length"):
        return f"content-length:{length}"
    raise NotCheckedDaily(NOT_DAILY)


def _differs(upstream: str, local: str) -> bool:
    """Whether the upstream validator is not among the local file's."""
    return upstream not in local.split(",")


def _cspec(body: bytes, headers: Mapping[str, str]) -> str:
    """`GN002=2.0,...`: the version CSpec's one list gives for each recorded
    specification the list names (one request, not one per gene)."""
    del headers
    wanted = {svi for _, svi, _ in cspec_rows()}
    listed = _json(body).get("data")
    if not isinstance(listed, list):
        raise ShapeError("no data")
    versions: dict[str, str] = {}
    for item in listed:  # pyright: ignore[reportUnknownVariableType]
        if not isinstance(item, dict):
            continue
        svi = str(item.get("@id", "")).rsplit("/", 1)[-1]  # pyright: ignore[reportUnknownMemberType]
        if svi not in wanted:
            continue
        version = str(item.get("version"))  # pyright: ignore[reportUnknownMemberType]
        if not re.fullmatch(_CSPEC_VERSION, version):
            raise ShapeError("not a specification version")
        if versions.setdefault(svi, version) != version:
            raise ShapeError("one specification listed twice")
    return ",".join(f"{svi}={versions[svi]}" for svi in sorted(versions))


# --- the version on this machine ------------------------------------------------------------


def _newest_dated(folder: Path) -> Path | None:
    found = (
        sorted(
            d for d in folder.iterdir() if d.is_dir() and re.fullmatch(_DATE, d.name)
        )
        if folder.is_dir()
        else []
    )
    return found[-1] if found else None


def _folder(name: str) -> Callable[[Path], str | None]:
    """The name of the newest dated folder under `<lab>/<name>/`."""

    def local(lab: Path) -> str | None:
        newest = _newest_dated(lab / name)
        return None if newest is None else newest.name

    return local


def _from_file(
    folder: str, name: str, read: Callable[[str], str]
) -> Callable[[Path], str | None]:
    """The version `read` finds in `name` of the newest dated folder under `<lab>/<folder>/`."""

    def local(lab: Path) -> str | None:
        newest = _newest_dated(lab / folder)
        path = None if newest is None else newest / name
        return read(path.read_text("utf-8")) if path and path.is_file() else None

    return local


def _validators(folder: str, name: str) -> Callable[[Path], str | None]:
    """The validators of `name` in the newest dated folder under `<lab>/<folder>/`: its md5
    (GenCC's ETag), its modification day (the Last-Modified `curl -R` put on it) and its
    length."""

    def local(lab: Path) -> str | None:
        newest = _newest_dated(lab / folder)
        path = None if newest is None else newest / name
        if path is None or not path.is_file():
            return None
        md5 = hashlib.md5()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                md5.update(chunk)
        stat = path.stat()
        day = _day(datetime.fromtimestamp(stat.st_mtime, UTC))
        return (
            f"etag:{md5.hexdigest()},last-modified:{day},content-length:{stat.st_size}"
        )

    return local


def _pinned(read: Callable[[], str]) -> Callable[[Path], str | None]:
    return lambda lab: read()


def _pharmcat() -> str:
    from genokeeper.pgx.pharmcat import PHARMCAT_VERSION

    return PHARMCAT_VERSION


def _hla() -> str:
    from genokeeper.pgx.hla import IMGT_HLA_VERSION

    return IMGT_HLA_VERSION


def _vep() -> str:
    from genokeeper.annotate.vep import VEP_VERSION

    return VEP_VERSION


def _haplogrep() -> str:
    from genokeeper.containers import IMAGES

    found = re.fullmatch(
        r"genomelab/haplogrep3:([0-9]+\.[0-9]+\.[0-9]+)-[0-9a-f]+", IMAGES["haplogrep3"]
    )
    if found is None:
        raise ShapeError("no haplogrep3 version in the image name")
    return found[1]


def _strchive() -> str:
    """The STRchive release the repeat thresholds were read from, named in the SOURCE.md
    beside the catalog (the catalog itself is Illumina's, probed as `expansionhunter`)."""
    text = files("genokeeper").joinpath("resources", "repeats", "SOURCE.md")
    found = re.search(r"STRchive v([0-9]+\.[0-9]+\.[0-9]+)", text.read_text("utf-8"))
    if found is None:
        raise ShapeError("no STRchive version in resources/repeats/SOURCE.md")
    return found[1]


def _expansionhunter() -> str:
    from genokeeper.repeats.expansionhunter import EH_VERSION

    # The update probe compares Illumina's release tags; the commit is named by the image digest.
    return EH_VERSION.split(" ", 1)[0]


def _yleaf() -> str:
    from genokeeper.yhaplogroup import YLEAF_VERSION

    return YLEAF_VERSION


def _gnomad() -> str:
    from genokeeper.evidence.frequency import GNOMAD_RELEASE

    return GNOMAD_RELEASE


def _scores(lab: Path) -> str | None:
    """The day the PGS scoring files were fetched: the newest file under `<lab>/pgs/scores`."""
    folder = lab / "pgs" / "scores"
    stamps = (
        [p.stat().st_mtime for p in folder.iterdir() if p.is_file()]
        if folder.is_dir()
        else []
    )
    return _day(datetime.fromtimestamp(max(stamps), UTC)) if stamps else None


def _cspec_local(lab: Path) -> str | None:
    return ",".join(f"{svi}={version}" for _, svi, version in cspec_rows())


def _re(pattern: str) -> re.Pattern[str]:
    return re.compile(pattern)


def _github(repo: str) -> str:
    return f"https://api.github.com/repos/{repo}/releases/latest"


_DAY, _VERSION = _re(_DATE), _re(_SEMVER)
PROBES: tuple[Probe, ...] = (
    Probe("clinvar", "refreshed",
          "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz.md5",
          _DAY, _md5_name, _folder("clinvar")),
    Probe("clingen", "refreshed", VALIDITY_URL,
          _DAY, _whole(validity_date),
          _from_file("clingen", FILE_NAME, validity_date)),
    Probe("clingen_dosage", "refreshed", DOSAGE_URL, _DAY, _whole(dosage_date),
          _from_file("clingen-dosage", DOSAGE_FILE, dosage_date)),
    Probe("cpic", "refreshed", _github("cpicpgx/cpic-data"),
          _DAY, _release_date, _folder("cpic")),
    Probe("pharmgkb", "refreshed",
          "https://api.clinpgx.org/v1/download/file/data/summaryAnnotations.zip",
          _DAY, _last_modified, _folder("pharmgkb")),
    Probe("gwas", "refreshed",
          "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/",
          _DAY, _gwas_release, _folder("gwas")),
    Probe("medlineplus", "refreshed", "https://medlineplus.gov/download/ghr-summaries.xml",
          _DAY, _last_modified, _folder("medlineplus")),
    Probe("litvar", "refreshed",
          "https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz",
          _DAY, _last_modified, _folder("litvar")),
    Probe("clinvar_citations", "refreshed",
          "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt",
          _DAY, _last_modified, _folder("clinvar-citations")),
    Probe("gencc", "refreshed", GENCC_URL, _re(_VALIDATORS), _validator,
          _validators("gencc", GENCC_FILE), _differs),
    Probe("cgd", "pinned", "https://research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz",
          _DAY, _last_modified, _folder("cgd")),
    Probe("pharmcat", "pinned", _github("PharmGKB/PharmCAT"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)"), _pinned(_pharmcat)),
    Probe("ipd_imgt_hla", "pinned", _github("ANHIG/IMGTHLA"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)-alpha"), _pinned(_hla)),
    Probe("ensembl_vep", "pinned", _github("Ensembl/ensembl-vep"),
          _VERSION, _tag(r"release/([0-9]+\.[0-9]+)"), _pinned(_vep)),
    Probe("haplogrep3", "pinned", _github("genepi/haplogrep3"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)"), _pinned(_haplogrep)),
    Probe("strchive", "pinned", _github("dashnowlab/STRchive"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)"), _pinned(_strchive)),
    Probe("expansionhunter", "pinned", _github("Illumina/ExpansionHunter"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)"), _pinned(_expansionhunter)),
    Probe("yleaf", "pinned", _github("genid/Yleaf"),
          _VERSION, _tag(r"v([0-9]+\.[0-9]+\.[0-9]+)"), _pinned(_yleaf)),
    Probe("pgs_catalog", "pinned", "https://www.pgscatalog.org/rest/info",
          _DAY, _pgs_release, _scores),
    Probe("gnomad", "pinned",
          "https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o"
          "?prefix=release/&delimiter=/&fields=prefixes,nextPageToken",
          _VERSION, _gnomad_release, _pinned(_gnomad)),
    Probe("cspec", "cspec", "https://cspec.genome.network/cspec/api/svis",
          _re(_CSPEC_LIST), _cspec, _cspec_local),
)  # fmt: skip
LABELS: dict[str, str] = {
    "clinvar": "ClinVar",
    "clingen": "ClinGen (validity)",
    "clingen_dosage": "ClinGen (dosage sensitivity)",
    "cpic": "CPIC",
    "pharmgkb": "ClinPGx",
    "gwas": "GWAS Catalog",
    "medlineplus": "MedlinePlus",
    "litvar": "LitVar",
    "clinvar_citations": "ClinVar (citations)",
    "gencc": "GenCC",
    "cgd": "CGD",
    "pharmcat": "PharmCAT",
    "ipd_imgt_hla": "IPD-IMGT/HLA",
    "ensembl_vep": "Ensembl VEP",
    "haplogrep3": "Haplogrep3",
    "strchive": "STRchive",
    "expansionhunter": "ExpansionHunter",
    "yleaf": "Yleaf",
    "pgs_catalog": "PGS Catalog",
    "gnomad": "gnomAD",
}
# The Russian names that differ from `LABELS`.
_LABELS_RU: dict[str, str] = {
    "clingen": "ClinGen (валидность)",
    "clingen_dosage": "ClinGen (дозовая чувствительность)",
    "clinvar_citations": "ClinVar (статьи)",
}
# A notice's sentences: a heading per kind of source, then when the last refresh was.
_NOTICE: dict[Language, dict[str, str]] = {
    "en": {
        "refreshed": "New source releases",
        "pinned": "New releases of pinned components (not applied by themselves)",
        "cspec": "ClinGen CSpec specifications changed",
        "never": "No refresh has run yet: genokeeper refresh.",
        "last": "The last refresh was {days} ago: genokeeper refresh.",
    },
    "ru": {
        "refreshed": "Новые версии источников",
        "pinned": "Вышли новые версии закреплённых компонентов (сами не применяются)",
        "cspec": "Изменились спецификации ClinGen CSpec",
        "never": "Обновление ещё не запускалось: genokeeper refresh.",
        "last": "Последнее обновление было {days} назад: genokeeper refresh.",
    },
}
# Recorded, never news: MedlinePlus is rebuilt daily and is no grading input; the PGS
# Catalog's release is catalog-wide and weekly, while the scoring files are sha-pinned.
QUIET = frozenset({"medlineplus", "pgs_catalog"})


# --- the probes -----------------------------------------------------------------------------


def head(url: str) -> Mapping[str, str]:
    """HEAD `url`, 30 s at most: its headers, and no body."""
    request = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT}, method="HEAD"
    )
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        return {k.lower(): v for k, v in response.headers.items()}


def http(url: str) -> tuple[bytes, Mapping[str, str]]:
    """GET `url`, 30 s at most, and read at most `MAX_BODY` bytes (`WHOLE_MAX` of a ClinGen
    file): the connection is closed after that, so a probe of a 2 GB file reads its headers
    and nothing more. A URL in `HEADS` gets a HEAD only (`head`), with an empty body."""
    if url in HEADS:
        return b"", head(url)
    whole = url in WHOLE
    asked = {"User-Agent": USER_AGENT}
    if not whole:
        asked["Range"] = f"bytes=0-{MAX_BODY - 1}"
    request = urllib.request.Request(url, headers=asked)
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        body = response.read(WHOLE_MAX if whole else MAX_BODY)
        return body, {k.lower(): v for k, v in response.headers.items()}


def _get(fetch: Http, url: str) -> tuple[bytes, Mapping[str, str]]:
    try:
        return fetch(url)
    except Exception:  # noqa: BLE001 - one retry, then the error is the probe's
        return fetch(url)


def _error_class(error: Exception) -> str:
    if isinstance(error, urllib.error.HTTPError):
        return f"http:{error.code}"
    name = type(error).__name__
    return name if refresh_log.ERROR_CLASS.fullmatch(name) else "Error"


def _key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.removeprefix("v").split("."))


def newer_than(upstream: str, local: str) -> bool:
    """Whether the version `upstream` is later than `local` (numbers by part, zero-padded)."""
    a, b = _key(upstream), _key(local)
    width = max(len(a), len(b))
    return a + (0,) * (width - len(a)) > b + (0,) * (width - len(b))


def _compare(probe: Probe, up: str, local: str | None) -> Found:
    if local is None:
        return Found(probe.name, probe.kind, up, None, "failed", "no-local-version")
    if not probe.pattern.fullmatch(local):
        return Found(probe.name, probe.kind, up, local, "failed", "local-shape")
    newer = (probe.compare or newer_than)(up, local)
    return Found(probe.name, probe.kind, up, local, "newer" if newer else "same")


def _cspec_found(up: str, local: str) -> list[Found]:
    """One entry per recorded specification: `newer` when its version differs."""
    now = dict(pair.split("=") for pair in up.split(","))
    found: list[Found] = []
    for svi, version in (pair.split("=") for pair in local.split(",")):
        name = f"cspec:{svi}"
        if svi not in now:
            found.append(Found(name, "cspec", None, version, "failed", "missing"))
            continue
        status: Status = "same" if now[svi] == version else "newer"
        found.append(Found(name, "cspec", now[svi], version, status))
    return found


def probe_one(probe: Probe, fetch: Http, lab: Path) -> list[Found]:
    """What `probe` found, as one entry (`cspec`: one per specification). Any error, and any
    response that does not have the probe's shape, is `failed`."""
    try:
        body, headers = _get(fetch, probe.url)
        up = probe.upstream(body, headers)
        if not probe.pattern.fullmatch(up):
            raise ShapeError("not the probe's version shape")
    except Exception as e:  # noqa: BLE001 - every failure reads "could not check"
        error = NOT_DAILY if isinstance(e, NotCheckedDaily) else _error_class(e)
        return [Found(probe.name, probe.kind, None, None, "failed", error)]
    try:
        local = probe.local(lab)
    except Exception as e:  # noqa: BLE001
        return [
            Found(
                probe.name, probe.kind, up, None, "failed", f"local:{_error_class(e)}"
            )
        ]
    if probe.kind == "cspec":
        return _cspec_found(up, local or "")
    return [_compare(probe, up, local)]


# --- check.json -----------------------------------------------------------------------------


def read_state(lab: Path) -> dict[str, Any] | None:
    """`<lab>/refresh/check.json`, or None when there is none or it is not one."""
    try:
        state = json.loads((lab / "refresh" / CHECK).read_text("utf-8"))
        datetime.fromisoformat(state["ran_at"])
        assert isinstance(state["sources"], list) and isinstance(
            state["notified"], dict
        )
    except OSError, ValueError, KeyError, TypeError, AssertionError:
        return None
    return state  # pyright: ignore[reportUnknownVariableType]


def _write_state(lab: Path, state: Mapping[str, Any]) -> None:
    folder = lab / "refresh"
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = folder / f"{CHECK}.{os.getpid()}.part"
    temporary.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, folder / CHECK)


@contextmanager
def _one_check(lab: Path) -> Iterator[None]:
    """Hold `<lab>/refresh/check.lock`: a check started meanwhile waits, then reads what
    this one announced."""
    folder = lab / "refresh"
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (folder / "check.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def check(lab: Path, fetch: Http, now: datetime) -> dict[str, Any]:
    """Probe every source, write `<lab>/refresh/check.json` and return it:
    `{"ran_at", "notified", "sources": [Found as a dict]}`. `notified` is kept from the
    last check."""
    previous = read_state(lab)
    found = [f for probe in PROBES for f in probe_one(probe, fetch, lab)]
    state: dict[str, Any] = {
        "ran_at": now.isoformat(),
        "notified": previous["notified"] if previous else {},
        "sources": [asdict(f) for f in found],
    }
    _write_state(lab, state)
    return state


def summary(state: Mapping[str, Any] | None) -> dict[str, Any]:
    """The part of a check that the dashboard and `changes` show: when it ran, what is
    newer, and what could not be checked. Never an age: the viewer computes it."""
    if state is None:
        return {"ran_at": None, "newer": [], "failed": []}
    sources = state["sources"]
    return {
        "ran_at": state["ran_at"],
        "newer": [
            {k: s[k] for k in ("name", "kind", "local", "upstream")}
            for s in sources
            if s["status"] == "newer" and s["name"] not in QUIET
        ],
        "failed": [s["name"] for s in sources if s["status"] == "failed"],
    }


# --- the notification -----------------------------------------------------------------------

_SHAPES = {p.name: p.pattern for p in PROBES}


def _valid(source: Mapping[str, Any]) -> bool:
    """Whether a recorded newer version still has its probe's shape (a text built from
    check.json carries nothing else)."""
    name = str(source["name"])
    pattern = _SHAPES["cspec" if name.startswith("cspec:") else name]
    version = str(source["upstream"])
    if name.startswith("cspec:"):
        return bool(
            re.fullmatch(r"GN[0-9]{3}", name[6:])
            and re.fullmatch(_CSPEC_VERSION, version)
        )
    return bool(pattern.fullmatch(version))


def pending(
    state: Mapping[str, Any], last_refresh: datetime | None, now: datetime
) -> dict[str, str]:
    """The `{name: version}` set a notification would announce: the newer pinned and CSpec
    versions, and the newer refreshed sources once the last refresh is `REFRESH_DAYS` old."""
    stale = last_refresh is None or (now - last_refresh).days >= REFRESH_DAYS
    return {
        s["name"]: s["upstream"]
        for s in state["sources"]
        if s["status"] == "newer"
        and _valid(s)
        and s["name"] not in QUIET
        and (stale or s["kind"] != "refreshed")
    }


def notice(
    state: Mapping[str, Any],
    last_refresh: datetime | None,
    now: datetime,
    *,
    lang: Language = "en",
) -> str | None:
    """The text in `lang` announcing what is newer, or None when there is nothing to say or
    this very set was already announced (`state["notified"]`). Source names, versions and
    days only: never a gene or a variant."""
    found = pending(state, last_refresh, now)
    if not found or found.items() <= state["notified"].items():
        return None
    kind = {s["name"]: s["kind"] for s in state["sources"]}
    text = _NOTICE[lang]
    labels = {**LABELS, **_LABELS_RU} if lang == "ru" else LABELS
    parts: list[str] = []
    for wanted in ("refreshed", "pinned", "cspec"):
        names = sorted(n for n in found if kind[n] == wanted)
        if names:
            listed = ", ".join(f"{labels.get(n, n)} {found[n]}" for n in names)
            parts.append(f"{text[wanted]}: {listed}.")
    if last_refresh is None:
        parts.append(text["never"])
    else:
        days = (now - last_refresh).days
        unit = "дн." if lang == "ru" else "day" if days == 1 else "days"
        parts.append(text["last"].format(days=f"{days} {unit}"))
    return " ".join(parts)


def mark_notified(lab: Path, state: dict[str, Any], names: Mapping[str, str]) -> None:
    """Record that `names` was announced, so the same set is not announced again. Earlier
    announcements stay: a probe that failed today leaves its version out of `names`."""
    state["notified"] = {**state["notified"], **names}
    _write_state(lab, state)


def last_refresh(lab: Path) -> datetime | None:
    """When the last refresh that did not fail ended, one in which every source failed
    (offline) not counted; None when there is none."""
    ended = [
        r.ended
        for r in refresh_log.read(lab)
        if r.outcome != "failed"
        and not (r.sources and all(s.status is SourceStatus.FAILED for s in r.sources))
    ]
    return max(ended) if ended else None


# --- command line ---------------------------------------------------------------------------


def run(
    lab: Path,
    fetch: Http,
    now: datetime,
    *,
    notifier: Callable[[str], int] | None,
    lang: Language = "en",
) -> tuple[dict[str, Any], int]:
    """`check`, then (with a `notifier`) announce what is news once, in `lang`: exit code 0,
    or 1 when the notifier failed (the same text is then tried again at the next check)."""
    code = 0
    with _one_check(lab):
        state = check(lab, fetch, now)
        if notifier is not None:
            refreshed = last_refresh(lab)
            text = notice(state, refreshed, now, lang=lang)
            if text is not None:
                if notifier(text) == 0:
                    mark_notified(lab, state, pending(state, refreshed, now))
                else:
                    print(
                        "genokeeper check-updates: the notifier failed", file=sys.stderr
                    )
                    code = 1
    return state, code


def _command(args: argparse.Namespace) -> int:
    notifier: Callable[[str], int] | None = None
    if args.notify:
        from genokeeper.schedule import notify

        notifier = notify
    state, code = run(
        load_paths().lab, http, datetime.now(UTC), notifier=notifier, lang=language()
    )
    for s in state["sources"]:
        print(f"{s['name']}: {s['status']} ({s['local']} -> {s['upstream']})")
    return code


def register(subparsers: Subparsers) -> None:
    """Add the `check-updates` sub-command."""
    c = subparsers.add_parser(
        "check-updates",
        help="ask each source for its newest version; write <lab>/refresh/check.json",
    )
    c.add_argument(
        "--notify",
        action="store_true",
        help="send a macOS notification when there is news (once per set of versions)",
    )
    c.set_defaults(handler=_command)
