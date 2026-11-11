"""`python -m genokeeper.demo package OUT`: a fresh, path-scrubbed store and the site's pages,
refused (nothing under a final name) when a needle survives in any output."""

import base64
import hashlib
import html
import json
import os
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures.store_mini import Row, seed

import genokeeper
from genokeeper import changes, refresh_log
from genokeeper.calibrate.tables import REPORT, TABLES, folder
from genokeeper.demo import (
    MANIFEST,
    REFRESH_LOG_FILE,
    SITE_FILES,
    STORE_FILE,
    demo_paths,
)
from genokeeper.demo import package as pkg
from genokeeper.demo.frame import CALIBRATION_NOTE, FRAME, TABLE_SAMPLES
from genokeeper.demo.package import Leak, package, scan_bytes
from genokeeper.query import Filters, open_store, query
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus
from genokeeper.store.db import (
    SCHEMA_VERSION,
    Level,
    Validity,
    connect,
    export_snapshots,
)

RUN1, RUN2 = "20261012T120000Z-aaaaaa", "20261109T120000Z-bbbbbb"
# Four findings enter or leave the trusted tier between the runs; two stay as they were.
BEFORE = {"clinvar:1": "limited", "clinvar:2": "strong", "clinvar:3": "conflicting"}
BEFORE |= {"clinvar:4": "moderate", "clinvar:5": "strong", "clinvar:6": "limited"}
AFTER = {**BEFORE, "clinvar:1": "strong", "clinvar:2": "conflicting"}
AFTER |= {"clinvar:3": "moderate", "clinvar:4": "limited"}
LEVELS = {"clinvar:5": "low", "clinvar:6": "medium"}
# The page's slider order: `rel=<i>` and `ev=<i>` index these.
REL = (Level.LOW, Level.MEDIUM, Level.HIGH)
EV = (Validity.CONFLICTING, Validity.LIMITED, Validity.MODERATE, Validity.STRONG)
SITE = "https://genokeeper.github.io/genokeeper/social.png"


def _rows(validity: dict[str, str]) -> list[Row]:
    return [
        Row(
            sid,
            gene=f"GENE{sid[-1]}",
            pos=1000 + int(sid[-1]),
            validity=v,
            level=LEVELS.get(sid, "high"),
        )
        for sid, v in validity.items()
    ]


def _lab(tmp_path: Path) -> tuple[Path, Path]:
    """A home holding a lab whose demo root holds a two-run store with the maintainer's
    paths planted (JSON, `VARCHAR[]`), its refresh-log record, and one table's report."""
    home = tmp_path / "Users" / "someone"
    lab = home / "genome-lab"
    root = demo_paths(lab)
    con = connect(root.reports / STORE_FILE)
    try:
        seed(con, _rows(BEFORE), run_id=RUN1)
        seed(con, _rows(AFTER), run_id=RUN2)
        tools = {"clinvar": "20261108", "deepvariant": "1.10.0"}
        con.execute(
            "UPDATE run SET tools = ? WHERE run_id = ?", [json.dumps(tools), RUN2]
        )
        inputs = {
            "fastq": f"{root.repo}/raw/r1.fq.gz",
            "clinvar38": {"path": f"{lab}/clinvar/x.vcf.gz", "sha256": "ab"},
            "vep": f"{home}/.vep/cache",
        }
        con.execute("UPDATE run SET inputs = ?", [json.dumps(inputs)])
        con.execute(
            "UPDATE call SET warnings = ? WHERE run_id = ? AND vid LIKE '%:1001:%'",
            [["first", f"log at {lab}/work/x.log"], RUN2],
        )
    finally:
        con.close()
    at = datetime(2026, 11, 10, 9, tzinfo=UTC)
    clinvar = SourceOutcome("clinvar", "20261011", "20261108", SourceStatus.NEW)
    record = RefreshRecord(at, at, "complete", (clinvar,), RUN1, None, RUN2)
    refresh_log.append(root.repo / "lab", record)
    report = folder(lab, "stage3") / REPORT
    report.parent.mkdir(parents=True)
    report.write_text(f"measured <b>on</b> {lab}/calibration/hg002\n")
    return lab, home


def _package(tmp_path: Path, private_ids: Path | None = None) -> Path:
    lab, home = _lab(tmp_path)
    return package(lab, tmp_path / "out", home=home, private_ids=private_ids)


def _needles(tmp_path: Path) -> dict[str, str]:
    home = tmp_path / "Users" / "someone"
    lab = home / "genome-lab"
    return {
        str(i): re.escape(str(p))
        for i, p in enumerate((home, lab, demo_paths(lab).repo))
    }


def _read(out: Path, name: str) -> str:
    return (out / name).read_text(encoding="utf-8")


def _page_data(page: str) -> dict[str, Any]:
    [text] = re.findall(
        r'<script type="application/json" id="data">(.*?)</script>', page, re.DOTALL
    )
    return json.loads(text)


def _report(out: Path, tmp_path: Path) -> dict[str, Any]:
    """`changes.report` at the default filters over the shipped store, as the demo folder
    will have it: its snapshots exported, a lab holding the shipped refresh log."""
    reports, lab = tmp_path / "check" / "reports", tmp_path / "check" / "lab"
    with open_store(out / STORE_FILE) as con:
        export_snapshots(con, reports / "snapshots")
    (lab / "refresh").mkdir(parents=True)
    shutil.copyfile(out / REFRESH_LOG_FILE, lab / "refresh" / refresh_log.LOG)
    return changes.report(
        reports,
        lab,
        base=None,
        run=None,
        f=Filters(),
        every=False,
        gene=None,
        marker=None,
        text=None,
    )


def test_shipped_store_bytes_hold_no_planted_path(tmp_path: Path) -> None:
    out = _package(tmp_path)
    assert scan_bytes(out / STORE_FILE, _needles(tmp_path)) == []
    assert str(tmp_path).encode() not in (out / STORE_FILE).read_bytes()
    with open_store(out / STORE_FILE) as con:
        rows = con.execute("SELECT inputs FROM run").fetchall()
    assert {json.loads(inputs)["fastq"] for (inputs,) in rows} == {
        "<demo>/raw/r1.fq.gz"
    }
    assert json.loads(rows[0][0])["clinvar38"]["path"] == "<lab>/clinvar/x.vcf.gz"
    assert json.loads(rows[0][0])["vep"] == "~/.vep/cache"


def test_scrub_covers_list_columns(tmp_path: Path) -> None:
    out = _package(tmp_path)
    with open_store(out / STORE_FILE) as con:
        found = con.execute(
            "SELECT warnings FROM call WHERE run_id = ? AND vid LIKE '%:1001:%'", [RUN2]
        ).fetchone()
    assert found == (["first", "log at <lab>/work/x.log"],)


def test_scan_refuses_planted_path_in_any_table(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "Users" / "someone"
    scrub = pkg.scrub

    def scrub_then_plant(con: duckdb.DuckDBPyConnection, replacements: Any) -> None:
        scrub(con, replacements)
        con.execute("CREATE TABLE notes (body VARCHAR)")
        con.execute("INSERT INTO notes VALUES (?)", [f"kept {home}/x"])

    monkeypatch.setattr(pkg, "scrub", scrub_then_plant)
    with pytest.raises(Leak, match=r"notes\.body"):
        _package(tmp_path)
    assert os.listdir(tmp_path / "out") == []


def test_private_ids_skip_comments_and_missing_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    lab, home = _lab(tmp_path)
    out = tmp_path / "out"
    # Each of these would match the pages were it a pattern: "#" a link, "" and "  " any.
    ids = tmp_path / "genome-private-ids"
    ids.write_text("#\n\n  \nKIT-[0-9]{4}\n")
    package(lab, out, home=home, private_ids=ids)
    assert set(os.listdir(out)) == {STORE_FILE, *SITE_FILES, MANIFEST, REFRESH_LOG_FILE}
    # A pattern line is a needle, in the store's text as anywhere.
    with connect(demo_paths(lab).reports / STORE_FILE) as con:
        con.execute(
            "UPDATE statement SET reason = 'kit KIT-0042' WHERE run_id = ?", [RUN2]
        )
    with pytest.raises(Leak, match=r"statement\.reason: private ID 1"):
        package(lab, tmp_path / "out2", home=home, private_ids=ids)
    assert os.listdir(tmp_path / "out2") == []
    capsys.readouterr()
    package(lab, tmp_path / "out3", home=home, private_ids=None)
    err = capsys.readouterr().err.splitlines()
    assert [line for line in err if "private-ID" in line] == [
        "genokeeper demo package: no private-ID file: that check is skipped"
    ]


def test_package_writes_site_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GENOKEEPER_LANG", "ru")
    out = _package(tmp_path)
    assert set(os.listdir(out)) == {STORE_FILE, *SITE_FILES, MANIFEST, REFRESH_LOG_FILE}
    manifest = json.loads(_read(out, MANIFEST))
    assert set(manifest) == {"version", "schema_version", "runs", "date", "files"}
    assert (manifest["version"], manifest["schema_version"]) == (
        genokeeper.__version__,
        SCHEMA_VERSION,
    )
    assert manifest["runs"] == [
        {"run_id": RUN1, "clinvar": "20261011"},
        {"run_id": RUN2, "clinvar": "20261108"},
    ]
    assert manifest["files"] == {
        name: hashlib.sha256((out / name).read_bytes()).hexdigest()
        for name in (STORE_FILE, *SITE_FILES, REFRESH_LOG_FILE)
    }
    first = html.escape(FRAME.split(". ")[0])
    for name in SITE_FILES:
        assert first in _read(out, name), name
    assert _page_data(_read(out, "dashboard.html"))["lang"] == "en"


def test_dashboard_shows_a_refresh_and_no_check_warning(tmp_path: Path) -> None:
    out = _package(tmp_path)
    shown = _page_data(_read(out, "dashboard.html"))["changes"]
    assert "check" not in shown
    assert shown["default_base"] == RUN1  # the refresh's base run
    assert shown["bases"][RUN1]["rows"]


def test_landing_numbers(tmp_path: Path) -> None:
    out = _package(tmp_path)
    index = _read(out, "index.html")
    with open_store(out / STORE_FILE) as con:
        every = Filters(min_reliability=REL[0], min_validity=EV[0])
        loosest = len(query(con, None, every).findings)
        default = len(query(con, None, Filters()).findings)
    assert (loosest, default) == (6, 2)
    assert f"One genome, {loosest:,} candidate findings." in index
    assert f"{default:,} pass at the default thresholds." in index
    n = len(_report(out, tmp_path)["changes"])
    assert n == 4
    assert f", {n} findings changed at the default thresholds." in index


def test_landing_matrix(tmp_path: Path) -> None:
    out = _package(tmp_path)
    cells = re.findall(
        r'<td( class="default")?><a class="t[0-5]"'
        r' href="dashboard\.html#rel=(\d)&amp;ev=(\d)">([\d,]+)',
        _read(out, "index.html"),
    )
    assert len(cells) == 12
    with open_store(out / STORE_FILE) as con:
        for default, r, v, count in cells:
            f = Filters(min_reliability=REL[int(r)], min_validity=EV[int(v)])
            assert int(count.replace(",", "")) == len(query(con, None, f).findings)
    assert [(r, v) for default, r, v, _ in cells if default] == [("1", "2")]
    assert {(r, v) for _, r, v, _ in cells} == {
        (str(r), str(v)) for r in range(3) for v in range(4)
    }


def test_landing_diff_lists_changes(tmp_path: Path) -> None:
    out = _package(tmp_path)
    index = _read(out, "index.html")
    [diff] = re.findall(
        r'<div class="diff"[^>]*aria-label="[^"]+">(.*?)\n</div>', index, re.DOTALL
    )
    rows = _report(out, tmp_path)["changes"]
    for listed in rows[:3]:
        label = f"{listed['gene']} {listed['marker']['grch38']}"
        assert label in diff
        lines = diff.split(label, 1)[1].split("\n")
        assert lines[1].startswith('<div class="m">- ') and (
            f"{listed['before']['validity']} evidence" in lines[1]
        )
        assert lines[2].startswith('<div class="p">+ ') and (
            f"{listed['after']['validity']} evidence" in lines[2]
        )
    assert rows[3]["gene"] not in diff
    assert 'href="dashboard.html#tab=changes"' in index


def test_landing_has_share_tags(tmp_path: Path) -> None:
    index = _read(_package(tmp_path), "index.html")
    tags = dict(
        re.findall(r'<meta (?:property|name)="([^"]+)" content="([^"]*)">', index)
    )
    assert tags["og:title"] and tags["og:description"] and tags["twitter:card"]
    assert tags["og:image"] == SITE
    assert re.search(r"<title>[^<]+</title>", index)


def test_calibration_summary_covers_every_table(tmp_path: Path) -> None:
    assert set(TABLE_SAMPLES) == set(TABLES)
    out = _package(tmp_path)
    page, index = _read(out, "calibration.html"), _read(out, "index.html")
    for name, measured in TABLE_SAMPLES.items():
        assert measured.for_hg002 in ("independent", "same sample"), name
        row = (
            f"<td>{html.escape(measured.measures)}</td>"
            f"<td>{html.escape(measured.samples)}</td>"
        )
        assert row in page and row in index
    assert TABLE_SAMPLES["small_variants"].for_hg002 == "independent"
    assert html.escape(CALIBRATION_NOTE.split(". ")[0]) in page
    # A table's report, its paths rewritten, escaped in its own section.
    assert (
        "<pre>measured &lt;b&gt;on&lt;/b&gt; &lt;lab&gt;/calibration/hg002\n</pre>"
        in page
    )


def _sha256(text: str) -> str:
    return base64.b64encode(hashlib.sha256(text.encode()).digest()).decode()


def test_site_pages_load_nothing(tmp_path: Path) -> None:
    out = _package(tmp_path)
    for name in ("index.html", "calibration.html"):
        page = _read(out, name)
        assert "<script" not in page and "<link" not in page and "@import" not in page
        assert re.findall(r"\bsrc\s*=|url\(", page) == []
        assert re.findall(r'content="(https?://[^"]*)"', page) == (
            [SITE] if name == "index.html" else []
        )
        [style] = re.findall(r"<style>(.*?)</style>", page, re.DOTALL)
        csp = f"default-src 'none'; style-src 'sha256-{_sha256(style)}'"
        assert f'<meta http-equiv="Content-Security-Policy" content="{csp}">' in page
