"""`python -m genokeeper.demo package OUT`: a fresh, path-scrubbed store and the site's pages,
refused (nothing under a final name) when a needle survives in any output."""

import hashlib
import html
import json
import os
import re
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures.demo_site import RUN1, RUN2, seeded
from fixtures.medlineplus_mini import gene, write_medlineplus
from fixtures.store_mini import Row

import genokeeper
from genokeeper import changes, refresh_log
from genokeeper.calibrate.tables import REPORT, TABLES, folder
from genokeeper.config import GRADING_FILE, load_grading
from genokeeper.dashboard import _sha256
from genokeeper.demo import (
    MANIFEST,
    REFRESH_LOG_FILE,
    SHIPPED,
    SITE_FILES,
    STORE_FILE,
    demo_paths,
)
from genokeeper.demo import package as pkg
from genokeeper.demo.frame import (
    CALIBRATION_NOTE,
    CALIBRATION_PAGE,
    FRAME,
    TABLE_SAMPLES,
)
from genokeeper.demo.package import Leak, package, scan_bytes
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE
from genokeeper.explain import LOOSEST, explain, newest_medlineplus
from genokeeper.query import Filters, open_store, query
from genokeeper.store.db import (
    SCHEMA_VERSION,
    Level,
    Validity,
    connect,
    export_snapshots,
)

# Seven changes: 1 and 2 enter or leave the default list, 6 leaves and 7, 8 enter the list of
# everything called, 3 and 4 change their grading inside both; 5 stays as it was.
BEFORE = {"clinvar:1": "limited", "clinvar:2": "strong", "clinvar:3": "conflicting"}
BEFORE |= {"clinvar:4": "moderate", "clinvar:5": "strong", "clinvar:6": "limited"}
AFTER = {**BEFORE, "clinvar:1": "strong", "clinvar:2": "conflicting"}
AFTER |= {"clinvar:3": "limited", "clinvar:4": "strong", "clinvar:7": "limited"}
AFTER = {sid: v for sid, v in AFTER.items() if sid != "clinvar:6"} | {
    "clinvar:8": "limited"
}
LEVELS = {
    "clinvar:5": "low",
    "clinvar:6": "medium",
    "clinvar:7": "low",
    "clinvar:8": "low",
}
# The landing's three counts: every change, then the entered or left at each list.
TIERS = ((LOOSEST, True), (LOOSEST, False), (Filters(), False))
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
    with connect(seeded(lab, _rows(BEFORE), _rows(AFTER))) as con:
        tools = {"clinvar": "20261108", "deepvariant": "1.10.0"}
        con.execute(
            "UPDATE run SET tools = ? WHERE run_id = ?", [json.dumps(tools), RUN2]
        )
        inputs = {
            "fastq": f"{root.repo}/raw/r1.fq.gz",
            "clinvar38": {"path": f"{lab}/clinvar/x.vcf.gz", "sha256": "ab"},
            "vep": f"{lab}/vep/cache",
            "grading": {"path": str(GRADING_FILE)},
        }
        con.execute("UPDATE run SET inputs = ?", [json.dumps(inputs)])
        con.execute(
            "UPDATE call SET warnings = ? WHERE run_id = ? AND vid LIKE '%:1001:%'",
            [["first", f"log at {lab}/work/x.log"], RUN2],
        )
    report = folder(lab, "stage3") / REPORT
    report.parent.mkdir(parents=True)
    report.write_text(f"measured <b>on</b> {lab}/calibration/hg002\n")
    return lab, home


def _package(tmp_path: Path) -> Path:
    lab, home = _lab(tmp_path)
    return package(
        lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True
    )


def _needles(tmp_path: Path) -> dict[str, str]:
    home = tmp_path / "Users" / "someone"
    lab = home / "genome-lab"
    return {str(i): str(p) for i, p in enumerate((home, lab, demo_paths(lab).repo))}


def _read(out: Path, name: str) -> str:
    return (out / name).read_text(encoding="utf-8")


def _page_data(page: str) -> dict[str, Any]:
    [text] = re.findall(
        r'<script type="application/json" id="data">(.*?)</script>', page, re.DOTALL
    )
    return json.loads(text)


def _report(
    out: Path, tmp_path: Path, f: Filters = TIERS[2][0], every: bool = False
) -> dict[str, Any]:
    """`changes.report` at `f` (the default filters) over the shipped store, as the demo
    folder will have it: its snapshots exported, a lab holding the shipped refresh log."""
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
        f=f,
        every=every,
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
    assert json.loads(rows[0][0])["vep"] == "<lab>/vep/cache"
    grading = json.loads(rows[0][0])["grading"]["path"]
    assert grading == "<genokeeper>/config/grading.toml"


@pytest.mark.parametrize(
    "planted", ["{elsewhere}/config/grading.toml", "~/Dev/genokeeper/config/g.toml"]
)
def test_scan_refuses_any_absolute_path_left(tmp_path: Path, planted: str) -> None:
    """A path under neither root, lab nor checkout (the grading config of another checkout,
    or one the home rewrite already shortened) survives every rewrite: the scan names it."""
    lab, home = _lab(tmp_path)
    path = planted.format(elsewhere=tmp_path / "elsewhere")
    with connect(demo_paths(lab).reports / STORE_FILE) as con:
        con.execute("UPDATE run SET inputs = ?", [json.dumps({"grading": path})])
    with pytest.raises(Leak, match=r"^refusing to package: run\.inputs: a path$"):
        package(lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True)
    assert os.listdir(tmp_path / "out") == []


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


def test_scan_reads_the_refresh_log(tmp_path: Path) -> None:
    """The log is copied as it is: the scan of the shipped files is its only control."""
    lab, home = _lab(tmp_path)
    with (demo_paths(lab).repo / "lab" / "refresh" / refresh_log.LOG).open("a") as log:
        log.write(f"cut at {home}/x\n")
    with pytest.raises(
        Leak, match=f"^refusing to package: {REFRESH_LOG_FILE}: the home"
    ):
        package(lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True)
    assert os.listdir(tmp_path / "out") == []


def test_private_ids_are_read_as_the_guard_reads_them(tmp_path: Path) -> None:
    """`grep -E`'s POSIX classes, in a page alone; a pattern it cannot read refuses."""
    lab, home = _lab(tmp_path)
    (folder(lab, "stage3") / REPORT).write_text("kit KIT-0042\n")
    ids = tmp_path / "genome-private-ids"
    ids.write_text("KIT-[[:digit:]]{4}\n")
    with pytest.raises(Leak, match=r"calibration\.html: private ID on line 1$"):
        package(lab, tmp_path / "out", home=home, private_ids=ids)
    ids.write_text("# a comment\nKIT-[\n")
    with pytest.raises(Leak, match="cannot read the private ID on line 2"):
        package(lab, tmp_path / "out", home=home, private_ids=ids)
    assert os.listdir(tmp_path / "out") == []


def test_package_refuses_an_out_holding_other_files(tmp_path: Path) -> None:
    lab, home = _lab(tmp_path)
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "shot.png").write_bytes(b"png")
    with pytest.raises(FileExistsError, match="holds shot.png"):
        package(lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True)
    assert os.listdir(tmp_path / "out") == ["shot.png"]


def test_private_ids_skip_comments_and_missing_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    lab, home = _lab(tmp_path)
    out = tmp_path / "out"
    # Each of these would match the pages were it a pattern: "#" a link, "" and "  " any.
    ids = tmp_path / "genome-private-ids"
    ids.write_text("#\n\n  \nKIT-[0-9]{4}\n")
    package(lab, out, home=home, private_ids=ids)
    assert set(os.listdir(out)) == {*SHIPPED, MANIFEST}
    # A pattern line is a needle, in the store's text as anywhere.
    with connect(demo_paths(lab).reports / STORE_FILE) as con:
        con.execute(
            "UPDATE statement SET reason = 'kit KIT-0042' WHERE run_id = ?", [RUN2]
        )
    with pytest.raises(Leak, match=r"statement\.reason: private ID on line 4"):
        package(lab, tmp_path / "out2", home=home, private_ids=ids)
    assert os.listdir(tmp_path / "out2") == []
    # Without the file, nothing is written unless the waiver names it; then the run says so.
    with pytest.raises(Leak, match="no private-ID file.*--no-private-ids"):
        package(lab, tmp_path / "out3", home=home, private_ids=None)
    assert not (tmp_path / "out3").exists()
    capsys.readouterr()
    package(lab, tmp_path / "out3", home=home, private_ids=None, no_private_ids=True)
    err = capsys.readouterr().err.splitlines()
    assert [line for line in err if "private-ID" in line] == [
        "genokeeper demo package: --no-private-ids: the private-ID check is skipped"
    ]


def test_package_writes_site_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GENOKEEPER_LANG", "ru")
    out = _package(tmp_path)
    assert set(os.listdir(out)) == {*SHIPPED, MANIFEST}
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
        name: hashlib.sha256((out / name).read_bytes()).hexdigest() for name in SHIPPED
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
    # The shipped log is what makes `changes --json` in the demo folder show a refresh.
    assert _report(out, tmp_path)["refresh"]["run"] == RUN2


def test_landing_numbers(tmp_path: Path) -> None:
    out = _package(tmp_path)
    index = _read(out, "index.html")
    with open_store(out / STORE_FILE) as con:
        every = Filters(min_reliability=REL[0], min_validity=EV[0])
        loosest = len(query(con, None, every).findings)
        default = len(query(con, None, Filters()).findings)
    assert (loosest, default) == (7, 2)
    assert f"One genome, {loosest:,} candidate findings." in index
    assert f"{default:,} remain at the default thresholds." in index
    n = [
        len(_report(out, tmp_path / str(i), *t)["changes"]) for i, t in enumerate(TIERS)
    ]
    assert n == [7, 3, 2]
    assert (
        f"ClinVar brought {n[0]} changes; {n[1]} findings entered or left the list of"
        f" everything called, and {n[2]} the list at the default thresholds.</p>"
    ) in index


def test_landing_matrix(tmp_path: Path) -> None:
    out = _package(tmp_path)
    cells = re.findall(
        r'<td( class="default")?><a class="t[0-5]"'
        r' href="dashboard\.html#rel=(\d)&amp;ev=(\d)">([\d,]+)',
        _read(out, "index.html"),
    )
    assert len(cells) == 12
    # "high only" on top: the display order flips, each link keeps its own indices.
    assert [r for _, r, _, _ in cells[::4]] == ["2", "1", "0"]
    with open_store(out / STORE_FILE) as con:
        for default, r, v, count in cells:
            f = Filters(min_reliability=REL[int(r)], min_validity=EV[int(v)])
            assert int(count.replace(",", "")) == len(query(con, None, f).findings)
    assert [(r, v) for default, r, v, _ in cells if default] == [("1", "2")]
    assert {(r, v) for _, r, v, _ in cells} == {
        (str(r), str(v)) for r in range(3) for v in range(4)
    }


def test_landing_examples_are_the_everything_called_rows(tmp_path: Path) -> None:
    out = _package(tmp_path)
    reports = [_report(out, tmp_path / str(i), *t) for i, t in enumerate(TIERS)]
    reports[0][
        "changes"
    ].reverse()  # the first three of every change are none of the list's
    rows, tier = reports[0]["changes"], reports[1]["changes"]
    assert {r["id"] for r in rows[:3]}.isdisjoint(r["id"] for r in tier)
    with open_store(out / STORE_FILE) as con:
        index = pkg._landing(tmp_path, con, pkg._runs(con), reports)  # pyright: ignore[reportPrivateUsage]
    [diff] = re.findall(
        r'<div class="diff"[^>]*aria-label="[^"]+">(.*?)\n</div>', index, re.DOTALL
    )
    assert len(tier) == 3 and "Three of them:" not in index
    for listed in tier:
        label = f"{listed['gene']} {listed['marker']['grch38']}"
        assert label in diff
        lines = diff.split(label, 1)[1].split("\n")
        before, after = (
            "not reported" if side is None else f"{side['validity']} evidence"
            for side in (listed["before"], listed["after"])
        )
        assert lines[1].startswith('<div class="m">- ') and before in lines[1]
        assert lines[2].startswith('<div class="p">+ ') and after in lines[2]
    assert rows[0]["gene"] not in diff
    assert 'href="dashboard.html#tab=changes&amp;rel=0&amp;ev=0"' in index


def test_landing_diff_reads_clinvar_as_words() -> None:
    side = {"kind": "finding", "validity": "conflicting", "level": "high"}
    side["source_level"] = "Conflicting_classifications_of_pathogenicity|risk_factor"
    side["reason"] = "gnomAD FAF above BA1"  # the reason the store holds, shown
    assert pkg._side(side) == (  # pyright: ignore[reportPrivateUsage]
        "Conflicting classifications of pathogenicity, risk factor,"
        " conflicting evidence (gnomAD FAF above BA1), genotype high"
    )


def test_landing_says_only_what_is_measured(tmp_path: Path) -> None:
    lab, home = _lab(tmp_path)
    with connect(demo_paths(lab).reports / STORE_FILE) as con:
        con.execute(
            "UPDATE statement SET topic = 'drugs' WHERE statement_id = 'clinvar:4'"
        )
    index = _read(
        package(
            lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True
        ),
        "index.html",
    )
    assert "The 2 at the default thresholds by topic: drugs 1, risks 1." in index
    for said in ("Every class of finding", "GitHub and sends nothing", "Anthropic"):
        assert said not in index
    for said in (
        "Polygenic scores, GWAS associations and VEP predictions carry none of their own",
        "This run is HG002's PCR-free library",
        "HG002 may also be in DeepVariant's training data.",
        "DeepVariant's own genotype quality",
        "The download sends nothing.",
        "with a local model nothing leaves the machine",
        "No finding is not an all-clear",
    ):
        assert said in index, said
    assert "HG002's father: a different library of a related sample" in CALIBRATION_NOTE
    assert "one Ashkenazi Jewish family" in CALIBRATION_NOTE
    assert "outside the calibrated regions" in CALIBRATION_NOTE
    assert "a different genome" not in CALIBRATION_NOTE
    assert TABLE_SAMPLES["small_variants"].for_hg002 == "related sample"


# Report lines as `calibrate` prints them (`build.report`, `classes.report`).
SMALL = "SNP F1 0.9912 (the part stops below 0.995); below the limit\nINDEL F1 0.9934\n"
DELETIONS = """short band: h 0.5, prior 0.1
  pcrfree: 10 rare calls, 9 false; truth 30 of 70 found
  pcr: 12 rare calls, 11 false; truth 31 of 70 found
  delly|short: p_wrong 0.997 (0.98-1)
  both|short: p_wrong 0.004 (0.001-0.01)
long band: h none, prior 0.2
  pcrfree: 3 rare calls, 1 false; truth 56 of 127 found
  dysgu|long: p_wrong 0.83 (0.5-0.95)
"""


def test_landing_numbers_come_from_reports_store_and_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab, home = _lab(tmp_path)
    for name, text in (("small_variants", SMALL), ("deletions", DELETIONS)):
        folder(lab, name).mkdir(parents=True)
        (folder(lab, name) / REPORT).write_text(text)
    with connect(demo_paths(lab).reports / STORE_FILE) as con:
        withheld = json.dumps({"withheld": [{"reason": "a"}, {"reason": "b"}]})
        con.execute(
            "UPDATE statement SET claim = json_merge_patch(claim, ?)"
            " WHERE statement_id = 'clinvar:3'",
            [withheld],
        )
    grading = replace(load_grading(), p_high=0.002, p_low=0.02, ba1=0.07)
    monkeypatch.setattr(pkg, "load_grading", lambda: grading)
    index = _read(
        package(
            lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True
        ),
        "index.html",
    )
    for said in (
        "Small-variant genotypes, chr1-22, on HG003: SNP F1 0.9912, INDEL F1 0.9934.",
        "truth events found: pcrfree 86 of 197, pcr 31 of 70;",
        "2 of 3 classes wrong with a chance above 0.02, up to 0.997.",
        "This run holds 2 withheld values",
        "is at most 0.002 ",
        "low above 0.02 ",
        "above BA1 (0.07)",
        "A clean result is not an exclusion",
    ):
        assert said in index, said


def test_landing_built_on(tmp_path: Path) -> None:
    index = _read(_package(tmp_path), "index.html")
    built = index.split("<h2>Built on</h2>")[1].split("<h2>")[0]
    tools = pkg.citations()
    groups = list(dict.fromkeys(t["group"] for t in tools))  # in the file's order
    assert set(groups) == set(pkg.GROUP_MEASURED)
    assert re.findall(r"<h3>([^<]+)</h3>", built) == [html.escape(g) for g in groups]
    for t in tools:
        link = f'<a href="{html.escape(t["url"])}">{html.escape(t["name"])}</a>'
        assert link in built and html.escape(t["citation"]) in built, t["name"]
    # The run's recorded versions (the fixture's `tools`), where `version_key` names one.
    assert '<a href="https://doi.org/10.1038/nbt.4235">DeepVariant</a> 1.10.0.' in built
    assert "ClinVar</a> 20261108." in built and "BWA-MEM</a>. Li" in built
    assert "not genokeeper's own accuracy" in built and "endorses" in built
    for word in ("validated", "clinical-grade", "powered by", "<img"):
        assert word not in index.lower()


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
    for measured in TABLE_SAMPLES.values():
        row = (
            f"<td>{html.escape(measured.measures)}</td>"
            f"<td>{html.escape(measured.samples)}</td>"
        )
        assert row in page and row in index
    assert html.escape(CALIBRATION_NOTE.split(". ")[0]) in page
    # A table's report, its paths rewritten, escaped in its own section.
    assert (
        "<pre>measured &lt;b&gt;on&lt;/b&gt; &lt;lab&gt;/calibration/hg002\n</pre>"
        in page
    )
    # The note points at the calibration page from the others only: a link on the landing.
    assert (
        CALIBRATION_PAGE not in page and "outside the calibrated regions.</p>" in page
    )
    link = '<a href="calibration.html">the calibration page</a> says which table'
    assert link in index and CALIBRATION_PAGE in _read(out, "dashboard.html")
    assert '<ul class="plain methods">\n    <li>Tools and versions' in page
    assert ".methods { overflow-wrap: anywhere; }" in page  # md5 lists wrap at 375 px


def test_site_pages_load_nothing(tmp_path: Path) -> None:
    out = _package(tmp_path)
    for name in ("index.html", "calibration.html"):
        page = _read(out, name)
        assert re.search(
            r'</div>\n<main class="sheet">.*</main>\n</body>', page, re.DOTALL
        )
        assert "<script" not in page and "<link" not in page and "@import" not in page
        assert re.findall(r"\bsrc\s*=|url\(", page) == []
        assert re.findall(r'content="(https?://[^"]*)"', page) == (
            [SITE] if name == "index.html" else []
        )
        [style] = re.findall(r"<style>(.*?)</style>", page, re.DOTALL)
        csp = f"default-src 'none'; style-src 'sha256-{_sha256(style)}'"
        assert f'<meta http-equiv="Content-Security-Policy" content="{csp}">' in page


def test_readme_built_on_names_every_citation() -> None:
    readme = (Path(__file__).parents[1] / "README.md").read_text(encoding="utf-8")
    listed = readme.split("## Built on")[-1].split("\n## ")[0]
    missing = [
        t["name"] for t in pkg.citations() if f"[{t['name']}]({t['url']})" not in listed
    ]
    assert missing == []


def test_dashboard_shows_the_demo_roots_explanations(tmp_path: Path) -> None:
    """Written as the demo's are: another model and tone, grounded in the lab's MedlinePlus."""
    lab, home = _lab(tmp_path)
    write_medlineplus(lab / "medlineplus" / "20261020" / MEDLINEPLUS_FILE,
                      gene("GENE2", "The GENE2 gene makes a protein."))  # fmt: skip
    with open_store(demo_paths(lab).reports / STORE_FILE) as con:
        [finding] = [
            f for f in query(con, None, LOOSEST).findings if f["gene"] == "GENE2"
        ]

    def ask(text: str) -> str:
        return json.dumps(
            {"title": "T.", "text": "X.", "cites": ["medlineplus:gene:GENE2"]}
        )

    cache = demo_paths(lab).reports / "explanations"
    model = "claude-sonnet-5-5"
    explain(
        [finding], newest_medlineplus(lab), cache, ask, model=model, tone="engaging"
    )
    out = package(
        lab, tmp_path / "out", home=home, private_ids=None, no_private_ids=True
    )
    shown = _page_data(_read(out, "dashboard.html"))["items"]["clinvar:2"][
        "explanation"
    ]
    assert shown == {"title": "T.", "text": "X.", "cites": ["medlineplus:gene:GENE2"],
                     "model": model}  # fmt: skip


def test_built_on_shows_each_version_once() -> None:
    tools = {"haplogrep3": "3.2.2", "haplogrep3_tree": "phylotree-rcrs@17.3"}
    tools |= {
        "acmg_sf": "3.3",
        "carrier_tier3": "2021",
        "gnomad_faf": "4.1, floor 1e-06",
    }
    built = pkg._built_on(tools)
    assert "PhyloTree Build 17</a> phylotree-rcrs@17.3. " in built
    assert "Haplogrep 3</a> 3.2.2. " in built
    for name in ("ACMG SF v3.3", "ACMG carrier screening tiers (2021)", "gnomAD v4.1"):
        assert re.search(rf">{re.escape(name)}[^<]*</a>\. ", built), name
    assert "Filtering allele frequency (FAF)</a> 4.1, floor 1e-06. " in built
