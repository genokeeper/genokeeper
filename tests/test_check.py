"""`genokeeper check-updates`: every probe against the response recorded from its live
endpoint (`fixtures/check_responses/`, recorded with curl on 2026-10-27), the strict version
shapes, and the notification's once-per-set rule."""

import hashlib
import json
import os
import re
import subprocess
import threading
import tomllib
import urllib.error
import urllib.request
from collections.abc import Callable, Mapping
from datetime import UTC, datetime, timedelta
from email.utils import parsedate_to_datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Self

import pytest

from genokeeper import check, cli, refresh, refresh_log, schedule
from genokeeper.check import PROBES, Probe, notice, probe_one
from genokeeper.check import check as run_check
from genokeeper.evidence.clingen import DOSAGE_FILE
from genokeeper.lang import Language
from genokeeper.refresh_log import RefreshRecord

RECORDED = Path(__file__).parent / "fixtures" / "check_responses"
URLS = {p.url: p for p in PROBES}
NOW = datetime(2026, 10, 27, 17, 0, tzinfo=UTC)
HOSTILE = '1.0" & do shell script "x'
GRADING = Path(__file__).parents[1] / "config" / "grading.toml"
BOOTSTRAP = Path(__file__).parents[1] / "scripts" / "bootstrap.sh"
REPEATS_SOURCE = (
    Path(__file__).parents[1] / "src/genokeeper/resources/repeats/SOURCE.md"
)
# A matching lab's GenCC export: the fixture's head, whose md5 a replayed HEAD answers as its
# ETag (GenCC's ETag is the file's md5; the recorded one is the real 25 MB export's).
GENCC_HEAD = Path(__file__).parent / "fixtures" / "gencc_head.tsv"
GENCC_MD5 = hashlib.md5(GENCC_HEAD.read_bytes()).hexdigest()
RECORDED_ETAG = "468125bc137a1c303d463b8b67a13632"

# What each probe reads out of its recorded response (the version it normalises to).
RECORDED_VERSION = {
    "clinvar": "20261025",
    "clingen": "20261025",
    "clingen_dosage": "20261023",
    "cpic": "20261026",
    "pharmgkb": "20261005",
    "gwas": "20261022",
    "medlineplus": "20261027",
    "litvar": "20261025",
    "clinvar_citations": "20261025",
    "cgd": "20250221",
    "pharmcat": "3.4.0",
    "ipd_imgt_hla": "3.66.0",
    "ensembl_vep": "116.3",
    "haplogrep3": "3.3.2",
    "strchive": "2.28.0",
    "expansionhunter": "5.0.0",
    "yleaf": "4.1.6",
    "pgs_catalog": "20261023",
    "gnomad": "4.1",
    "gencc": f"etag:{GENCC_MD5}",
}

type Reply = tuple[bytes, Mapping[str, str]]


def recorded(name: str) -> Reply:
    body = RECORDED / f"{name}.body"
    headers: dict[str, str] = {}
    head = RECORDED / f"{name}.headers"
    for line in head.read_text().splitlines() if head.is_file() else []:
        key, _, value = line.partition(":")
        if value:
            headers[key.strip().lower()] = value.strip()
    return (body.read_bytes() if body.is_file() else b""), headers


def replay(url: str) -> Reply:
    name = URLS[url].name
    body, headers = recorded(name)
    if name == "gencc":  # see GENCC_MD5
        headers = {**headers, "etag": f'"{GENCC_MD5}"'}
    return body, headers


def bootstrap_day(name: str, tmp_path: Path) -> str:
    """The folder the bootstrap files `name`'s download under, by its own `file_day` on a
    file dated as `curl -R` dates it (the recorded Last-Modified), in a zone west of UTC."""
    stamp = tmp_path / f"{name}.stamp"
    stamp.touch()
    at = parsedate_to_datetime(recorded(name)[1]["last-modified"]).timestamp()
    os.utime(stamp, (at, at))
    script = ["bash", "-c", 'source "$0"; file_day "$1"', str(BOOTSTRAP), str(stamp)]
    west = {**os.environ, "TZ": "Pacific/Pago_Pago"}
    return subprocess.check_output(script, env=west, text=True).strip()


def lab_matching(tmp_path: Path) -> Path:
    """A lab whose newest versions are the recorded ones: what the bootstrap and a refresh
    on 2026-10-27 would have filed. A folder dated by the server's Last-Modified is named
    by the bootstrap's own code; the others by the version their file carries."""
    lab = tmp_path / "lab"
    lab.mkdir(parents=True)
    for folder, day in (
        ("clinvar", "20261025"),
        ("cpic", "20261027"),
        ("pharmgkb", "20261005"),
        ("gwas", "20261022"),
        ("medlineplus", "20261027"),
        ("litvar", bootstrap_day("litvar", lab)),
        ("clinvar-citations", bootstrap_day("clinvar_citations", lab)),
        ("cgd", bootstrap_day("cgd", lab)),
    ):
        (lab / folder / day).mkdir(parents=True)
    for folder, name, body in (
        ("clingen", "gene-validity.csv", "clingen"),
        ("clingen-dosage", DOSAGE_FILE, "clingen_dosage"),
    ):
        (lab / folder / "20261027").mkdir(parents=True)
        (lab / folder / "20261027" / name).write_bytes(recorded(body)[0])
    (lab / "gencc" / "20261025").mkdir(parents=True)
    (lab / "gencc" / "20261025" / "gencc-submissions.tsv").write_bytes(
        GENCC_HEAD.read_bytes()
    )
    scores = lab / "pgs" / "scores"
    scores.mkdir(parents=True)
    (scores / "PGS000001.txt.gz").write_bytes(b"")
    fetched = datetime(2026, 10, 23, 12, tzinfo=UTC).timestamp()
    os.utime(scores / "PGS000001.txt.gz", (fetched, fetched))
    return lab


def by_name(state: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {s["name"]: s for s in state["sources"]}


def reply_with(
    name: str, body: bytes, headers: Mapping[str, str]
) -> Callable[[str], Reply]:
    """`replay`, except that the probe `name` gets `body` and `headers`."""

    def fetch(url: str) -> Reply:
        return (body, headers) if URLS[url].name == name else replay(url)

    return fetch


# --- each probe, against its recorded response -----------------------------------------------


@pytest.mark.parametrize("probe", PROBES, ids=lambda p: p.name)
def test_each_probe_reads_its_recorded_response(probe: Probe, tmp_path: Path) -> None:
    found = probe_one(probe, replay, lab_matching(tmp_path))
    if probe.kind == "cspec":
        assert len(found) == len(check.cspec_rows())
        assert all(f.kind == "cspec" and f.upstream for f in found)
        return
    [one] = found
    assert (one.name, one.upstream) == (probe.name, RECORDED_VERSION[probe.name])
    assert probe.pattern.fullmatch(one.upstream or "")


@pytest.mark.parametrize("probe", PROBES, ids=lambda p: p.name)
def test_matching_versions_read_same(probe: Probe, tmp_path: Path) -> None:
    """D4: the recorded response and the version the bootstrap would have filed under it
    are one scheme per probe, so they read `same`."""
    found = probe_one(probe, replay, lab_matching(tmp_path))
    assert found and all(f.status == "same" and f.error is None for f in found)


def test_matching_lab_reads_same_for_every_source(tmp_path: Path) -> None:
    state = run_check(lab_matching(tmp_path), replay, NOW)
    assert {s["status"] for s in state["sources"]} == {"same"}
    assert {p.name for p in PROBES if p.kind != "cspec"} == set(RECORDED_VERSION)


def test_an_older_release_on_disk_reads_newer(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    for folder in (
        "clinvar",
        "cpic",
        "pharmgkb",
        "gwas",
        "litvar",
        "clinvar-citations",
    ):
        (lab / folder / "20200101").mkdir(parents=True)
    state = by_name(run_check(lab, replay, NOW))
    assert state["clinvar"]["status"] == "newer"
    assert (state["clinvar"]["local"], state["clinvar"]["upstream"]) == (
        "20200101",
        "20261025",
    )
    assert all(
        state[n]["status"] == "newer" for n in ("cpic", "pharmgkb", "gwas", "litvar")
    )
    # A release newer than the upstream one (a fetch date after a release date) is not older.
    (lab / "cpic" / "20261230").mkdir()
    assert by_name(run_check(lab, replay, NOW))["cpic"]["status"] == "same"
    # No release on disk at all cannot be compared: could not check.
    assert state["medlineplus"]["status"] == "failed"
    assert state["medlineplus"]["error"] == "no-local-version"


@pytest.mark.parametrize(
    ("upstream", "local", "newer"),
    [
        ("3.5.0", "3.4.0", True),
        ("3.4.0", "3.4.0", False),
        ("3.10", "3.9", True),
        ("4.1", "4.1.0", False),
        ("116.4", "116.3", True),
        ("20261026", "20261027", False),
        ("v3.5.0", "3.4.0", True),
    ],
)
def test_versions_compare_by_part(upstream: str, local: str, newer: bool) -> None:
    assert check.newer_than(upstream, local) is newer


def test_pinned_newer_is_reported_never_applied(tmp_path: Path) -> None:
    lab = lab_matching(tmp_path)
    release = json.loads(recorded("pharmcat")[0])
    release["tag_name"] = "v3.5.0"
    state = by_name(
        run_check(lab, reply_with("pharmcat", json.dumps(release).encode(), {}), NOW)
    )
    assert state["pharmcat"]["kind"] == "pinned"
    assert (state["pharmcat"]["status"], state["pharmcat"]["upstream"]) == (
        "newer",
        "3.5.0",
    )
    assert state["pharmcat"]["local"] == "3.4.0"
    # The bucket's 4.1.2 is a patch release of 4.1: not a new data release; 4.2 is.
    assert state["gnomad"]["status"] == "same"
    bucket = {"prefixes": ["release/4.1.2/", "release/4.2/", "release/lr/"]}
    state = by_name(
        run_check(lab, reply_with("gnomad", json.dumps(bucket).encode(), {}), NOW)
    )
    assert (state["gnomad"]["status"], state["gnomad"]["upstream"]) == ("newer", "4.2")


def test_pinned_versions_are_the_analysers_pins() -> None:
    from genokeeper.annotate.vep import VEP_VERSION
    from genokeeper.pgx.hla import IMGT_HLA_VERSION
    from genokeeper.pgx.pharmcat import PHARMCAT_VERSION
    from genokeeper.repeats.expansionhunter import EH_VERSION
    from genokeeper.yhaplogroup import YLEAF_VERSION

    here = {p.name: p.local(Path("/nonexistent")) for p in PROBES if p.kind == "pinned"}
    assert here["pharmcat"] == PHARMCAT_VERSION
    assert here["ipd_imgt_hla"] == IMGT_HLA_VERSION
    assert here["ensembl_vep"] == VEP_VERSION
    assert here["gnomad"] == "4.1"
    assert here["yleaf"] == YLEAF_VERSION
    assert here["expansionhunter"] == EH_VERSION.split(" ", 1)[0]
    source = REPEATS_SOURCE.read_text("utf-8")
    assert f"STRchive v{here['strchive']}" in source
    assert re.fullmatch(r"[0-9]+(\.[0-9]+){2}", here["haplogrep3"] or "")


def test_refreshed_probes_are_the_refreshed_sources(tmp_path: Path) -> None:
    """A source the refresh fetches has a probe, in the folder the bootstrap files it under."""
    refreshed = {p.name: p for p in PROBES if p.kind == "refreshed"}
    assert set(refreshed) == set(refresh.REFRESHED)
    lab = tmp_path / "lab"
    for name, probe in refreshed.items():
        if name in ("clingen", "clingen_dosage", "gencc"):  # read from the newest file
            continue
        (lab / refresh.FOLDERS[name] / "20990102").mkdir(parents=True)
        (lab / refresh.FOLDERS[name] / "20990101").mkdir(parents=True)
        assert probe.local(lab) == "20990102", name


def test_gencc_probe_never_gets(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """GenCC's export is built per request and limited to 20 a day, so its probe sends one
    HEAD (exempt from that limit) and no GET: its ETag (the file's md5), else Last-Modified,
    else Content-Length, against the newest download's. Without any of the three it says it
    is not checked daily and asks nothing more."""
    gencc = URLS[check.GENCC_URL]
    assert (gencc.name, gencc.kind) == ("gencc", "refreshed")
    lab = lab_matching(tmp_path)
    local = lab / "gencc" / "20261025" / "gencc-submissions.tsv"
    sent: list[tuple[str, str]] = []
    answer: Mapping[str, str] = recorded("gencc")[1]

    class Response:
        headers = property(lambda self: answer)

        def __enter__(self) -> Self:
            return self

        def __exit__(self, *args: object) -> None:
            return None

        def read(self, *args: object) -> bytes:
            raise AssertionError("a HEAD has no body to read")

    def urlopen(request: urllib.request.Request, timeout: float) -> Response:
        sent.append((request.get_method(), request.full_url))
        return Response()

    monkeypatch.setattr(check.urllib.request, "urlopen", urlopen)
    [found] = probe_one(gencc, check.http, lab)
    assert sent == [("HEAD", check.GENCC_URL)]
    day = check._day(datetime.fromtimestamp(local.stat().st_mtime, UTC))
    mine = f"etag:{GENCC_MD5},last-modified:{day},content-length:{local.stat().st_size}"
    assert (found.status, found.upstream, found.local) == (
        "newer",
        f"etag:{RECORDED_ETAG}",
        mine,
    )
    # Without an ETag: the Last-Modified day, which `curl -R` put on the download ...
    modified = recorded("gencc")[1]["last-modified"]
    stamp = parsedate_to_datetime(modified).timestamp()
    os.utime(local, (stamp, stamp))
    answer = {"content-type": "text/tab-separated-values", "last-modified": modified}
    assert [f.status for f in probe_one(gencc, check.http, lab)] == ["same"]
    # ... else its length.
    answer = {"content-type": "text/plain", "content-length": "1"}
    [found] = probe_one(gencc, check.http, lab)
    assert (found.status, found.upstream) == ("newer", "content-length:1")
    # None of the three: not checked daily, and nothing but the one HEAD was sent.
    answer = {"content-type": "text/tab-separated-values", "server": "nginx"}
    [found] = probe_one(gencc, check.http, lab)
    assert (found.status, found.upstream, found.error) == (
        "failed",
        None,
        "not checked daily: GenCC gives no change date; genokeeper refresh fetches it",
    )
    assert [method for method, _ in sent] == ["HEAD"] * 4
    assert {url for _, url in sent} == {check.GENCC_URL}


# --- could not check, never same -------------------------------------------------------------


def hostile_reply(probe: Probe, text: str) -> Reply:
    """The probe's own response with `text` where its version is read."""
    by_kind: dict[str, Any] = {
        "cpic": {"published_at": text},
        "pgs_catalog": {"latest_release": {"date": text}},
        "gnomad": {"prefixes": [text]},
        "cspec": {"data": [{"@id": "x/GN002", "version": text}]},
    }
    document = by_kind.get(
        probe.name,
        {"tag_name": text, "draft": False, "prerelease": False},
    )
    if probe.url.startswith("https://api.github.com") or probe.name in by_kind:
        return json.dumps(document).encode(), {"content-type": "application/json"}
    return text.encode(), {"last-modified": text, "content-type": "text/plain"}


@pytest.mark.parametrize("probe", PROBES, ids=lambda p: p.name)
@pytest.mark.parametrize(
    "text",
    [HOSTILE, "3.4.0\n", "٣.٤.٠", "20261025 && reboot", "../../etc", ""],
    ids=["quote-and-script", "newline", "arabic-digits", "shell", "path", "empty"],
)
def test_hostile_version_is_could_not_check(
    probe: Probe, text: str, tmp_path: Path
) -> None:
    """Review focus 5: a `"` in a GitHub tag, a garbled header: failed, with no upstream
    string recorded and nothing in the state that came from the response."""
    lab = lab_matching(tmp_path)
    state = run_check(lab, lambda url: hostile_reply(probe, text), NOW)
    mine = [s for s in state["sources"] if s["name"].split(":")[0] == probe.name]
    assert mine and all(s["status"] == "failed" for s in mine)
    assert all(s["upstream"] is None for s in mine)
    raw = (lab / "refresh" / "check.json").read_text()
    assert "do shell script" not in raw and "reboot" not in raw


ARABIC = "٢٠٢٦١٠٢٥"  # eight digits that `\\d` would take for a date


@pytest.mark.parametrize(
    ("read", "body"),
    [
        (check._md5_name, f"{'a' * 32}  clinvar_{ARABIC}.vcf.gz\n"),
        (check._release_date, json.dumps({"published_at": "٢٠٢٦-١٠-٢٦T١٠:٠٠:٠٠Z"})),
        (
            next(p.upstream for p in PROBES if p.name == "pharmcat"),
            json.dumps({"tag_name": "v٣.٤.٠"}),
        ),
        (check._pgs_release, json.dumps({"latest_release": {"date": "٢٠٢٦-١٠-٢٣"}})),
        (check._gnomad_release, json.dumps({"prefixes": ["release/٤.١/"]})),
        (
            check._gwas_release,
            (
                '<a href="gwas-catalog-associations_ontology-annotated-full.zip">x</a>'
                '</td><td align="right">٢٠٢٦-١٠-٢٢ 11:40'
            ),
        ),
        (
            URLS[check.VALIDITY_URL].upstream,
            '"GENE SYMBOL","CLASSIFICATION DATE"\n"A","٢٠٢٦-١٠-٢٥T00:00:00"\n',
        ),
    ],
    ids=["clinvar", "cpic", "tag", "pgs", "gnomad", "gwas", "clingen"],
)
def test_unicode_digits_are_not_digits(read: Callable[..., str], body: str) -> None:
    with pytest.raises(ValueError):
        read(body.encode(), {})


def test_a_folder_named_in_unicode_digits_is_not_a_release(tmp_path: Path) -> None:
    (tmp_path / "clinvar" / ARABIC).mkdir(parents=True)
    assert check._folder("clinvar")(tmp_path) is None


def test_hostile_version_beside_good_ones_fails_only_its_source(tmp_path: Path) -> None:
    release = json.loads(recorded("pharmcat")[0])
    release["tag_name"] = f"v{HOSTILE}"
    state = by_name(
        run_check(
            lab_matching(tmp_path),
            reply_with("pharmcat", json.dumps(release).encode(), {}),
            NOW,
        )
    )
    assert state["pharmcat"]["status"] == "failed"
    assert {s["status"] for n, s in state.items() if n != "pharmcat"} == {"same"}


def test_captive_portal_html_is_could_not_check(tmp_path: Path) -> None:
    """A captive portal answers every URL with its own page, even with a Last-Modified."""
    page = b"<html><body>Sign in to the Wi-Fi network</body></html>"
    headers = {
        "content-type": "text/html; charset=utf-8",
        "last-modified": "Tue, 27 Oct 2026 07:31:57 GMT",
    }
    state = run_check(lab_matching(tmp_path), lambda url: (page, headers), NOW)
    assert state["sources"]
    assert {s["status"] for s in state["sources"]} == {"failed"}


def test_offline_is_could_not_check_for_every_probe(tmp_path: Path) -> None:
    lab = lab_matching(tmp_path)
    asked: list[str] = []

    def offline(url: str) -> Reply:
        asked.append(url)
        raise urllib.error.URLError("no route to host")

    state = run_check(lab, offline, NOW)
    assert {s["status"] for s in state["sources"]} == {"failed"}
    assert {s["error"] for s in state["sources"]} == {"URLError"}
    assert {s["name"] for s in state["sources"]} == {p.name for p in PROBES}
    assert state["ran_at"] == NOW.isoformat()
    assert len(asked) == 2 * len(PROBES)  # one retry each
    on_disk = json.loads((lab / "refresh" / "check.json").read_text())
    assert on_disk == state


def test_one_probe_failing_fails_only_that_source(tmp_path: Path) -> None:
    def fetch(url: str) -> Reply:
        if URLS[url].name == "gwas":
            raise urllib.error.HTTPError(url, 503, "busy", {}, None)  # pyright: ignore[reportArgumentType]
        return replay(url)

    state = by_name(run_check(lab_matching(tmp_path), fetch, NOW))
    assert (state["gwas"]["status"], state["gwas"]["error"]) == ("failed", "http:503")
    assert {s["status"] for n, s in state.items() if n != "gwas"} == {"same"}


def test_a_probe_that_works_on_the_retry_is_not_failed(tmp_path: Path) -> None:
    calls: list[str] = []

    def flaky(url: str) -> Reply:
        calls.append(url)
        if URLS[url].name == "cpic" and calls.count(url) == 1:
            raise TimeoutError
        return replay(url)

    state = by_name(run_check(lab_matching(tmp_path), flaky, NOW))
    assert state["cpic"]["status"] == "same"


def test_validity_probe_reads_the_newest_classification_not_the_request_day() -> None:
    """The export's FILE CREATED line is the day it was requested; comparing it would call
    every request newer. The newest classification date is the export's version."""
    body, _ = recorded("clingen")
    assert b"FILE CREATED: 2026-10-27" in body
    assert check.validity_date(body.decode()) == "20261025"


def test_dosage_probe_reads_the_newest_evaluation_not_the_build_day() -> None:
    """The dosage list names the day it was built, every night; its newest "Date Last
    Evaluated" is its version, as the bootstrap files it only when its rows change."""
    body, _ = recorded("clingen_dosage")
    assert body.splitlines()[1] == b"#27 Oct,2026"
    assert check.dosage_date(body.decode()) == "20261023"


def test_validity_probe_reads_an_export_past_the_probe_cap() -> None:
    body, _ = recorded("clingen")
    rows = body.splitlines(keepends=True)
    grown = b"".join([*rows, *rows[-1:] * (3 * check.MAX_BODY // len(rows[-1]))])
    assert len(grown) > 2 * check.MAX_BODY
    assert URLS[check.VALIDITY_URL].upstream(grown, {}) == "20261025"


def test_a_version_too_long_to_be_a_number_is_could_not_check(tmp_path: Path) -> None:
    release = json.loads(recorded("pharmcat")[0])
    release["tag_name"] = "v1.2." + "9" * 5000
    lab = lab_matching(tmp_path)
    state = by_name(
        run_check(lab, reply_with("pharmcat", json.dumps(release).encode(), {}), NOW)
    )
    assert state["pharmcat"]["status"] == "failed"
    assert (lab / "refresh" / "check.json").is_file()


def test_a_last_modified_that_is_not_the_file_is_could_not_check(
    tmp_path: Path,
) -> None:
    """Only a response typed as something other than a web page, and not one in its body,
    says when the file changed."""
    litvar, lab = next(p for p in PROBES if p.name == "litvar"), lab_matching(tmp_path)
    old = {"last-modified": "Wed, 01 Jan 2020 00:00:00 GMT"}
    typed = {**old, "content-type": "application/octet-stream"}
    for reply in (
        (b"<html>portal</html>", old),
        (b"", old),
        (b"<!DOCTYPE html>", typed),
        # Typed as a page, whose body does not open like one (a BOM, a comment, a prolog).
        (b"\xef\xbb\xbfSign in", {**old, "content-type": "text/html"}),
        (b"<!-- portal --><html>", {**old, "content-type": "text/html; charset=utf-8"}),
    ):
        [found] = probe_one(litvar, lambda url, reply=reply: reply, lab)
        assert found.status == "failed", reply


def test_quiet_sources_are_never_shown_as_newer() -> None:
    """MedlinePlus (rebuilt daily, no grading input) and the PGS Catalog (a weekly catalog
    release against sha-pinned files) are recorded, never announced nor listed as newer."""
    state = newer_state(
        medlineplus="20261027", pgs_catalog="20261030", clinvar="20261101"
    )
    assert set(check.pending(state, None, NOW)) == {"clinvar"}
    assert [s["name"] for s in check.summary(state)["newer"]] == ["clinvar"]


# --- ClinGen CSpec ---------------------------------------------------------------------------


def cspec_document(changes: Mapping[str, str], drop: tuple[str, ...] = ()) -> bytes:
    listed = json.loads(recorded("cspec")[0])
    data: list[dict[str, Any]] = []
    for item in listed["data"]:
        svi = item["@id"].rsplit("/", 1)[1]
        if svi in drop:
            continue
        data.append({**item, "version": changes.get(svi, item["version"])})
    return json.dumps({"data": data}).encode()


def test_cspec_version_change_is_reported(tmp_path: Path) -> None:
    lab = lab_matching(tmp_path)
    body = cspec_document({"GN002": "2.1"}, drop=("GN105",))
    state = by_name(run_check(lab, reply_with("cspec", body, {}), NOW))
    assert state["cspec:GN002"]["status"] == "newer"
    assert (state["cspec:GN002"]["local"], state["cspec:GN002"]["upstream"]) == (
        "2.0",
        "2.1",
    )
    assert state["cspec:GN105"]["status"] == "failed"
    assert state["cspec:GN105"]["error"] == "missing"
    others = {n: s for n, s in state.items() if n.startswith("cspec:")}
    assert {
        s["status"]
        for n, s in others.items()
        if n not in ("cspec:GN002", "cspec:GN105")
    } == {"same"}
    text = notice(run_check(lab, reply_with("cspec", body, {}), NOW), NOW, NOW)
    assert text is not None and "GN002 2.1" in text


def test_cspec_list_is_one_request_for_every_recorded_gene(tmp_path: Path) -> None:
    asked: list[str] = []

    def spy(url: str) -> Reply:
        asked.append(url)
        return replay(url)

    run_check(lab_matching(tmp_path), spy, NOW)
    assert asked.count(URLS["https://cspec.genome.network/cspec/api/svis"].url) == 1
    assert not [u for u in asked if "cspec" in u and not u.endswith("/svis")]


def test_cspec_resource_lists_the_vcep_genes_of_the_thresholds() -> None:
    rows = check.cspec_rows()
    genes = {gene for gene, _, _ in rows}
    thresholds = tomllib.loads(GRADING.read_text())["bs1_by_gene"]
    assert genes <= set(thresholds)
    assert {"BRCA1", "MYH7", "TP53", "LDLR"} <= genes
    assert all(re.fullmatch(r"GN[0-9]{3}", svi) for _, svi, _ in rows)
    assert all(re.fullmatch(r"[0-9]+(\.[0-9]+){1,2}", v) for _, _, v in rows)


def test_a_specification_listed_twice_fails_the_probe(tmp_path: Path) -> None:
    listed = json.loads(cspec_document({}))
    listed["data"].append({**listed["data"][0], "version": "9.9"})
    state = run_check(
        lab_matching(tmp_path),
        reply_with("cspec", json.dumps(listed).encode(), {}),
        NOW,
    )
    cspec = [s for s in state["sources"] if s["name"].startswith("cspec")]
    assert [(s["name"], s["status"]) for s in cspec] == [("cspec", "failed")]


def test_each_specification_version_has_its_shape(tmp_path: Path) -> None:
    """A version field that smuggles another specification in is no version."""
    listed = json.loads(cspec_document({"GN105": "1.0,GN101=99.0"}, drop=("GN101",)))
    state = run_check(
        lab_matching(tmp_path),
        reply_with("cspec", json.dumps(listed).encode(), {}),
        NOW,
    )
    cspec = [s for s in state["sources"] if s["name"].startswith("cspec")]
    assert [(s["name"], s["status"]) for s in cspec] == [("cspec", "failed")]


# --- the requests ----------------------------------------------------------------------------


def test_every_requested_url_is_a_probe_constant(tmp_path: Path) -> None:
    """O11: the same URLs, in the same order, whatever the lab holds; each is a constant of
    `PROBES`, over https, with nothing to fill in."""
    seen: list[list[str]] = []
    for lab in (lab_matching(tmp_path / "a"), tmp_path / "b" / "empty"):
        asked: list[str] = []

        def spy(url: str, asked: list[str] = asked) -> Reply:
            asked.append(url)
            return replay(url)

        run_check(lab, spy, NOW)
        seen.append(asked)
    assert seen[0] == seen[1] == [p.url for p in PROBES]
    assert len({p.url for p in PROBES}) == len(PROBES)
    for p in PROBES:
        assert p.url.startswith("https://") and not re.search(r"[{}$%]", p.url)


def test_the_real_http_reads_the_headers_and_not_the_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A probe of a huge file reads at most `MAX_BODY` bytes of it; a ClinGen file is read
    whole."""
    big = b"x" * (3 * check.MAX_BODY)

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            if self.path == "/missing":
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Last-Modified", "Sun, 25 Oct 2026 10:38:43 GMT")
            self.send_header("Content-Length", str(len(big)))
            self.end_headers()
            try:
                self.wfile.write(big)
            except OSError:
                pass  # the client hung up after its cap: what is being tested

        def log_message(self, format: str, *args: Any) -> None:
            pass

    for key in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("no_proxy", "127.0.0.1")
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        base = f"http://127.0.0.1:{server.server_address[1]}"
        body, headers = check.http(f"{base}/big")
        assert len(body) == check.MAX_BODY
        assert headers["last-modified"] == "Sun, 25 Oct 2026 10:38:43 GMT"
        # A ClinGen file is read whole: its version is in its rows.
        monkeypatch.setattr(check, "WHOLE", frozenset({f"{base}/big"}))
        assert check.http(f"{base}/big")[0] == big
        with pytest.raises(urllib.error.HTTPError) as caught:
            check.http(f"{base}/missing")
        assert caught.value.code == 404
    finally:
        server.shutdown()
        server.server_close()


# --- check.json and the notification ---------------------------------------------------------


def newer_state(**versions: str) -> dict[str, Any]:
    kinds = {p.name: p.kind for p in PROBES}
    return {
        "ran_at": NOW.isoformat(),
        "notified": {},
        "sources": [
            {
                "name": name,
                "kind": kinds.get(name, "cspec"),
                "upstream": version,
                "local": "1",
                "status": "newer",
                "error": None,
            }
            for name, version in versions.items()
        ],
    }


def test_check_json_keeps_what_was_notified(tmp_path: Path) -> None:
    lab = lab_matching(tmp_path)
    first = run_check(lab, replay, NOW)
    check.mark_notified(lab, first, {"clinvar": "20261025"})
    second = run_check(lab, replay, NOW + timedelta(days=1))
    assert second["notified"] == {"clinvar": "20261025"}
    assert check.read_state(lab) == second
    (lab / "refresh" / "check.json").write_text("{not json")
    assert check.read_state(lab) is None
    assert run_check(lab, replay, NOW)["notified"] == {}


def test_notice_needs_thirty_days_for_refreshed_sources() -> None:
    state = newer_state(clinvar="20261101", cpic="20261026")
    assert notice(state, NOW - timedelta(days=29), NOW) is None
    text = notice(state, NOW - timedelta(days=30), NOW, lang="ru")
    assert text is not None
    assert "ClinVar 20261101" in text and "CPIC 20261026" in text and "30 дн." in text
    # A refresh that never happened is as old as can be.
    assert notice(state, None, NOW) is not None
    # A pinned component or a specification needs no age.
    for name, version in (("pharmcat", "3.5.0"), ("cspec:GN002", "2.1")):
        text = notice(newer_state(**{name: version}), NOW - timedelta(days=1), NOW)
        assert text is not None and version in text and "1 day ago" in text


def test_notice_leaves_out_what_is_not_newer_or_is_quiet() -> None:
    state = newer_state(medlineplus="20261027")
    assert notice(state, None, NOW) is None  # rebuilt daily, no grading input
    state = newer_state(clinvar="20261101")
    state["sources"].append({**state["sources"][0], "name": "cpic", "status": "same"})
    state["sources"].append({**state["sources"][0], "name": "gwas", "status": "failed"})
    text = notice(state, None, NOW)
    assert text is not None and "CPIC" not in text and "GWAS" not in text


def test_notice_fires_once_per_newer_set(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    state = newer_state(clinvar="20261101", pharmcat="3.5.0")
    old = NOW - timedelta(days=40)
    first = notice(state, old, NOW)
    assert first is not None
    check.mark_notified(lab, state, check.pending(state, old, NOW))
    assert notice(state, old, NOW) is None
    assert notice(state, old, NOW + timedelta(days=3)) is None
    # What was announced and is no longer newer leaves nothing to say.
    fewer = newer_state(pharmcat="3.5.0")
    fewer["notified"] = state["notified"]
    assert notice(fewer, old, NOW) is None
    # A version not announced before is news, and the text names it.
    more = newer_state(clinvar="20261108", pharmcat="3.5.0")
    more["notified"] = state["notified"]
    again = notice(more, old, NOW)
    assert again is not None and "20261108" in again


def test_a_version_announced_earlier_stays_announced_after_its_probe_fails(
    tmp_path: Path,
) -> None:
    lab, old = tmp_path / "lab", NOW - timedelta(days=40)
    day1 = newer_state(pharmcat="3.5.0")
    check.mark_notified(lab, day1, check.pending(day1, old, NOW))
    # Day 2: the PharmCAT probe failed, so ClinVar alone is announced.
    day2 = newer_state(clinvar="20261101")
    day2["notified"] = day1["notified"]
    check.mark_notified(lab, day2, check.pending(day2, old, NOW))
    # Day 3: both are newer again; only a version not announced before is news.
    day3 = newer_state(clinvar="20261101", pharmcat="3.5.0")
    day3["notified"] = day2["notified"]
    assert notice(day3, old, NOW) is None


def test_a_hand_edited_version_never_reaches_the_text() -> None:
    """The text is built from check.json: a value that has no longer its probe's shape is
    left out."""
    state = newer_state(clinvar=HOSTILE, pharmcat="3.5.0", **{"cspec:GN002": '2.0"'})
    text = notice(state, None, NOW)
    assert text is not None and "3.5.0" in text
    assert '"' not in text and "shell" not in text


def test_a_notice_is_english_unless_russian_is_asked() -> None:
    state = newer_state(clinvar="20261101", clingen="20261101")
    english = notice(state, NOW - timedelta(days=30), NOW)
    assert english == (
        "New source releases: ClinGen (validity) 20261101, ClinVar 20261101."
        " The last refresh was 30 days ago: genokeeper refresh."
    )
    russian = notice(state, None, NOW, lang="ru")
    assert russian == (
        "Новые версии источников: ClinGen (валидность) 20261101, ClinVar 20261101."
        " Обновление ещё не запускалось: genokeeper refresh."
    )


@pytest.mark.parametrize("lang", ["en", "ru"])
def test_every_notice_template_is_in_its_language_and_names_no_gene_or_variant(
    lang: Language,
) -> None:
    genes = set(tomllib.loads(GRADING.read_text())["bs1_by_gene"])
    genes |= {gene for gene, _, _ in check.cspec_rows()} | {"CYP2D6", "APOE", "HLA-B"}
    everything = newer_state(
        **{
            p.name: "20261101" if p.pattern.pattern.startswith("[0-9]{8}") else "9.9.9"
            for p in PROBES
            if p.kind != "cspec"
        },
        **{"cspec:GN002": "2.1"},
    )
    cases = [
        (everything, NOW - timedelta(days=45)),
        (everything, None),
        (newer_state(clinvar="20261101"), NOW - timedelta(days=31)),
        (newer_state(pharmcat="3.5.0"), NOW - timedelta(days=2)),
        (newer_state(**{"cspec:GN002": "2.1"}), NOW - timedelta(days=2)),
    ]
    for state, refreshed in cases:
        text = notice(state, refreshed, NOW, lang=lang)
        assert text is not None
        assert bool(re.search("[А-Яа-я]", text)) is (lang == "ru") and len(text) < 900
        assert not re.search(r"""["'`\\$;&|<>\n]""", text)
        assert not re.search(r"\brs[0-9]+|\bchr[0-9XYM]|[0-9]+:[ACGT]+>?", text)
        assert not genes & set(re.findall(r"[A-Za-z0-9-]+", text)), text
    assert set(check.LABELS) == {p.name for p in PROBES if p.kind != "cspec"}


# --- the command -----------------------------------------------------------------------------


def command_lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    lab = tmp_path / "lab"
    (lab / "clinvar" / "20260901").mkdir(parents=True)
    old = NOW - timedelta(days=60)
    refresh_log.append(
        lab, RefreshRecord(old, old, "complete", (), None, None, "r", None)
    )
    monkeypatch.setattr(check, "http", replay)
    monkeypatch.setattr(check, "load_paths", lambda: type("P", (), {"lab": lab})())
    return lab


def test_command_notifies_once_and_retries_a_failed_notifier(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = command_lab(tmp_path, monkeypatch)
    sent: list[str] = []
    code = [1]

    def notifier(text: str) -> int:
        sent.append(text)
        return code[0]

    monkeypatch.setattr(schedule, "notify", notifier)
    assert cli.main(["check-updates", "--notify"]) == 1  # the notifier failed
    assert len(sent) == 1 and "ClinVar" in sent[0]
    assert check.read_state(lab)["notified"] == {}  # type: ignore[index]
    code[0] = 0
    assert cli.main(["check-updates", "--notify"]) == 0
    assert len(sent) == 2 and sent[0] == sent[1]
    assert check.read_state(lab)["notified"]["clinvar"] == "20261025"  # type: ignore[index]
    assert cli.main(["check-updates", "--notify"]) == 0
    assert len(sent) == 2  # the same set is not announced twice
    assert "clinvar: newer (20260901 -> 20261025)" in capsys.readouterr().out


@pytest.mark.parametrize(("env", "word"), [(None, "release"), ("ru", "версии")])
def test_command_notifies_in_the_language_of_genokeeper_lang(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env: str | None, word: str
) -> None:
    command_lab(tmp_path, monkeypatch)
    if env is None:
        monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    else:
        monkeypatch.setenv("GENOKEEPER_LANG", env)
    sent: list[str] = []
    monkeypatch.setattr(schedule, "notify", lambda text: sent.append(text) or 0)
    assert cli.main(["check-updates", "--notify"]) == 0
    assert len(sent) == 1 and word in sent[0]


def test_two_concurrent_checks_announce_a_set_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A check started while another is still announcing waits for it, then reads that
    it was announced."""
    lab = command_lab(tmp_path, monkeypatch)
    sent: list[str] = []
    entered, release = threading.Event(), threading.Event()

    def notifier(text: str) -> int:
        sent.append(text)
        entered.set()
        release.wait(5)
        return 0

    def one() -> None:
        check.run(lab, replay, NOW, notifier=notifier)

    first = threading.Thread(target=one)
    first.start()
    assert entered.wait(5)
    second = threading.Thread(target=one)
    second.start()
    second.join(0.5)  # without a lock it is announcing too, within this wait
    release.set()
    first.join(5)
    second.join(5)
    assert len(sent) == 1


def test_check_json_is_written_through_a_name_of_its_own_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    moved: list[str] = []
    real = os.replace
    monkeypatch.setattr(
        os, "replace", lambda a, b: (moved.append(Path(a).name), real(a, b))
    )
    run_check(tmp_path / "lab", replay, NOW)
    assert moved == [f"check.json.{os.getpid()}.part"]


def test_command_without_notify_only_writes_check_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = command_lab(tmp_path, monkeypatch)

    def refuse(text: str) -> int:
        raise AssertionError("no notification without --notify")

    monkeypatch.setattr(schedule, "notify", refuse)
    assert cli.main(["check-updates"]) == 0
    assert check.read_state(lab) is not None
