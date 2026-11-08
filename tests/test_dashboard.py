"""The dashboard's data is the query's: one view per filter combination, the findings once.
Its page keeps every source text as data, under a CSP pinned to the inlined script and style."""

import base64
import hashlib
import json
import re
import shutil
import subprocess
from dataclasses import replace
from datetime import UTC, datetime
from importlib.resources import files
from itertools import product
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures import snapshots_mini as sm
from fixtures.medlineplus_mini import condition, gene, write_medlineplus
from fixtures.store_mini import RUN, TOOLS, Row, seed

from genokeeper import dashboard, refresh_log
from genokeeper.carrier import NOT_COVERED, load_mechanisms
from genokeeper.changes import NO_BASE, report
from genokeeper.cli import main
from genokeeper.config import Paths
from genokeeper.dashboard import payload, render
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE, load_medlineplus
from genokeeper.explain import MODEL, explain, lookup, watch_items
from genokeeper.history import Change, diff, kind_of, states, trusted_ids
from genokeeper.lang import LANGUAGES
from genokeeper.query import Filters, open_store, query
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus
from genokeeper.sources.gwas import GWAS_TOOL
from genokeeper.store.db import Level, Topic, Validity, connect, export_snapshot

PGX_CLAIM = {"schema": "genomelab.claim/1", "type": "pharmcat_gene", "drugs": []}
GENE_CALL = {
    "diplotype": "*1/*2",
    "phenotype": "Intermediate Metabolizer",
    "activity_score": 1.5,
    "callers": {"BCyrius": "*1/*2", "Aldy": "*1/*2"},
    "agreed": True,
    "level": "high",
    "warnings": [],
}
NOCALL = "genotype unknown: nocall"
NO_CACHE = Path(__file__).parent / "no-explanations"

# Findings across levels, validities and topics, a gene result, a chrX and a chrM row, and
# two unassessed claims (a NoCall site, a gene without a gene call).
ROWS = [
    Row("clinvar:1", rsid="rs1", pos=1001, actionable=True),
    Row("clinvar:2", level="medium", validity="moderate", pos=1002, topic="carrier"),
    Row("clinvar:3", level="low", validity="limited", pos=1003),
    Row("clinvar:4", validity="conflicting", pos=1004, topic="carrier"),
    Row("clinvar:5", gene="DMD", chrom="chrX", gt="1", validity="moderate"),
    Row(
        "clinvar:6",
        gene="MT-ND4",
        chrom="chrM",
        gt="1",
        level="medium",
        warnings=["heteroplasmy at level 0.23"],
    ),
    Row(
        "pharmcat:CYP2D6",
        gene="CYP2D6",
        level=None,
        topic="drugs",
        gene_call=GENE_CALL,
        claim=PGX_CLAIM,
    ),
    Row(
        "clinvar:7",
        kind="unassessed",
        level="low",
        dosage="unknown",
        pos=1007,
        warnings=[NOCALL],
    ),
    Row(
        "pharmcat:CYP2C19",
        kind="unassessed",
        gene="CYP2C19",
        level=None,
        topic="drugs",
        claim=PGX_CLAIM,
    ),
]


@pytest.fixture
def con(tmp_path: Path) -> duckdb.DuckDBPyConnection:
    con = connect(tmp_path / "genomelab.duckdb")
    seed(con, ROWS)
    return con


def test_every_view_is_the_query(con: duckdb.DuckDBPyConnection) -> None:
    data = payload(con, None, NO_CACHE)
    combinations = list(product(Level, Validity, (None, *Topic)))
    assert len(data["views"]) == len(combinations)
    for reliability, validity, topic in combinations:
        result = query(
            con,
            None,
            Filters(topic=topic, min_reliability=reliability, min_validity=validity),
        )
        assert data["views"][f"{reliability}|{validity}|{topic or ''}"] == {
            "findings": [f["id"] for f in result.findings],
            "hidden": result.hidden,
            "unassessed": result.unassessed,
        }


def test_items_hold_every_finding_once(con: duckdb.DuckDBPyConnection) -> None:
    data = payload(con, None, NO_CACHE)
    loosest = query(
        con, None, Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    )
    assert (data["run_id"], data["versions"]) == (RUN, TOOLS)
    assert [
        {k: v for k, v in item.items() if k not in ("explanation", "search")}
        for item in data["items"].values()
    ] == loosest.findings
    assert {item["explanation"] for item in data["items"].values()} == {None}
    assert {i for view in data["views"].values() for i in view["findings"]} == set(
        data["items"]
    )
    assert {item["kind"] for item in data["items"].values()} == {"finding"}


def test_default_view_is_the_spec_default(con: duckdb.DuckDBPyConnection) -> None:
    data = payload(con, None, NO_CACHE)
    assert data["default"] == "medium|moderate|"
    assert data["default"] in data["views"]


def _sha256(text: str) -> str:
    return base64.b64encode(hashlib.sha256(text.encode()).digest()).decode()


def _resource(name: str) -> str:
    folder = files("genokeeper").joinpath("resources", "dashboard")
    return folder.joinpath(name).read_text(encoding="utf-8")


def test_csp_first_in_head(con: duckdb.DuckDBPyConnection) -> None:
    page = render(payload(con, None, NO_CACHE))
    csp = re.match(
        r'<meta http-equiv="Content-Security-Policy" content="([^"]*)">',
        page.split("<head>", 1)[1].lstrip(),
    )
    assert csp is not None
    [script] = re.findall(r"<script>(.*?)</script>", page, re.DOTALL)
    [style] = re.findall(r"<style>(.*?)</style>", page, re.DOTALL)
    assert (script, style) == (_resource("dashboard.js"), _resource("dashboard.css"))
    assert csp[1] == (
        f"default-src 'none'; script-src 'sha256-{_sha256(script)}';"
        f" style-src 'sha256-{_sha256(style)}'"
    )


def test_hostile_text_stays_data(con: duckdb.DuckDBPyConnection) -> None:
    hostile = "</script><img src=x onerror=alert(1)>"
    seed(con, [Row("clinvar:9", gene="<b>X</b>", condition=hostile, pos=1009)])
    data = payload(con, None, NO_CACHE)
    page = render(data)
    assert "</script><img" not in page and "<b>X" not in page and "<b>x" not in page
    [text] = re.findall(
        r'<script type="application/json" id="data">(.*?)</script>', page, re.DOTALL
    )
    assert json.loads(text) == data
    assert data["items"]["clinvar:9"]["claim"]["condition"] == hostile
    # The search string carries the same markup, lower-cased, and stays escaped too.
    assert hostile in data["items"]["clinvar:9"]["search"]
    assert "<b>x</b>" in data["items"]["clinvar:9"]["search"]


def test_script_builds_no_markup_from_strings() -> None:
    script = _resource("dashboard.js")
    for sink in (
        "innerHTML",
        "outerHTML",
        "insertAdjacentHTML",
        "document.write",
        "eval(",
        "new Function",
    ):
        assert sink not in script


def test_dashboard_command_writes_the_page(
    con: duckdb.DuckDBPyConnection,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    con.close()
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    out = tmp_path / "reports" / "dashboard.html"
    db = str(tmp_path / "genomelab.duckdb")
    cache = str(NO_CACHE)
    assert main(["dashboard", "--db", db, "--cache", cache, "--out", str(out)]) == 0
    assert '<script type="application/json" id="data">' in out.read_text(
        encoding="utf-8"
    )
    assert list(out.parent.iterdir()) == [out]
    # No MedlinePlus file at all: the page is built, with one warning naming the lab.
    [warning] = capsys.readouterr().err.splitlines()
    assert str(tmp_path / "lab") in warning and "scripts/bootstrap.sh" in warning


def test_dashboard_command_refuses_missing_store(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    missing = tmp_path / "nowhere.duckdb"
    out = tmp_path / "dashboard.html"
    assert main(["dashboard", "--db", str(missing), "--out", str(out)]) == 2
    [err] = capsys.readouterr().err.splitlines()
    assert err.startswith("genokeeper dashboard: ") and str(missing) in err
    assert not out.exists() and not missing.exists()


# --- slice 1c: explanations on the page ----------------------------------------------------

HOSTILE_EXPLANATION = "</script><img src=x onerror=alert(1)> и <b>X</b>"


def _page_data(page: str) -> dict[str, Any]:
    [text] = re.findall(
        r'<script type="application/json" id="data">(.*?)</script>', page, re.DOTALL
    )
    return json.loads(text)


def _explained(
    con: duckdb.DuckDBPyConnection,
    tmp_path: Path,
    text: str = "Объяснение.",
    title: str | None = None,
) -> tuple[Path, Path]:
    """A lab's MedlinePlus file, and a cache `explain` filled through a stub `Ask`
    that answers `text` (headed `title`, or `text` again) citing each finding's own
    source."""
    ghr = write_medlineplus(
        tmp_path / "lab" / "medlineplus" / "20261020" / MEDLINEPLUS_FILE,
        gene("BRCA1", "The BRCA1 gene makes a protein."),
        condition("hboc", "Hereditary breast ovarian cancer", "A condition."),
    )
    cache = tmp_path / "explanations"
    findings = list(payload(con, None, NO_CACHE)["items"].values())

    def ask(prompt: str) -> str:
        own = re.search(r'"id": "([^"]+)"', prompt[prompt.index("\n{\n") :])
        assert own is not None
        answer = {"title": title or text, "text": text, "cites": [own[1]]}
        return json.dumps(answer, ensure_ascii=False)

    assert explain(findings, load_medlineplus(ghr), cache, ask)["written"] == len(
        findings
    )
    return ghr, cache


def test_cached_explanation_shown(
    con: duckdb.DuckDBPyConnection,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    ghr, cache = _explained(con, tmp_path)
    medlineplus = load_medlineplus(ghr)
    data = payload(con, None, cache, medlineplus)
    assert data["items"]["clinvar:1"]["explanation"] == {
        "title": "Объяснение.",
        "text": "Объяснение.",
        "cites": ["clinvar:1"],
        "model": MODEL,
    }
    # The command finds the lab's MedlinePlus and the cache, so the keys match.
    con.close()
    monkeypatch.setenv("GENOKEEPER_HOME", str(ghr.parents[2]))
    out = tmp_path / "dashboard.html"
    db = str(tmp_path / "genomelab.duckdb")
    assert (
        main(["dashboard", "--db", db, "--cache", str(cache), "--out", str(out)]) == 0
    )
    assert _page_data(out.read_text(encoding="utf-8")) == data
    # A MedlinePlus file that does not parse: a warning, and the page built without the
    # explanations it grounded, keeping one whose sources did not change.
    ghr.write_text("<html>maintenance</html>", encoding="utf-8")
    assert (
        main(["dashboard", "--db", db, "--cache", str(cache), "--out", str(out)]) == 0
    )
    assert str(ghr) in capsys.readouterr().err
    items = _page_data(out.read_text(encoding="utf-8"))["items"]
    assert items["clinvar:1"]["explanation"] is None
    assert items["pharmcat:CYP2D6"]["explanation"]["cites"] == ["pharmcat:CYP2D6"]


def test_stale_explanation_not_shown(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """An explanation is shown while its sources' text stands, across a new release of
    them, and not once the claim's text changed."""
    ghr, cache = _explained(con, tmp_path)
    medlineplus = load_medlineplus(ghr)
    finding = payload(con, None, NO_CACHE)["items"]["clinvar:1"]
    assert lookup(finding, medlineplus, cache) is not None
    newer = {**finding, "sources": [{**finding["sources"][0], "version": "20261101"}]}
    assert lookup(newer, medlineplus, cache) is not None
    reworded = {**finding, "claim": {**finding["claim"], "condition": "Another"}}
    assert lookup(reworded, medlineplus, cache) is None


def test_hostile_explanation_stays_data(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    ghr, cache = _explained(con, tmp_path, HOSTILE_EXPLANATION)
    medlineplus = load_medlineplus(ghr)
    data = payload(con, None, cache, medlineplus)
    page = render(data)
    assert "</script><img" not in page and "<b>X" not in page
    assert _page_data(page) == data
    explanation = data["items"]["clinvar:1"]["explanation"]
    assert explanation["title"] == explanation["text"] == HOSTILE_EXPLANATION


def test_item_search_text(con: duckdb.DuckDBPyConnection, tmp_path: Path) -> None:
    """The page's search string: the query's `--text` fields and the explanation's
    words, lower-cased, with ClinVar's `_` read as a space."""
    warfarin = {**PGX_CLAIM, "drugs": [{"drug": "Warfarin", "source": "CPIC"}]}
    seed(
        con,
        [
            Row("clinvar:9", condition="Hereditary_hemochromatosis", pos=1009),
            Row(
                "pharmcat:CYP2C9",
                gene="CYP2C9",
                level=None,
                topic="drugs",
                gene_call=GENE_CALL,
                claim=warfarin,
            ),
        ],
    )
    ghr, cache = _explained(
        con, tmp_path, "Наследственный Гемохроматоз.", "Избыток железа"
    )
    items = payload(con, None, cache, load_medlineplus(ghr))["items"]
    assert "hereditary hemochromatosis" in items["clinvar:9"]["search"]
    assert "warfarin" in items["pharmcat:CYP2C9"]["search"]
    assert "избыток железа" in items["clinvar:9"]["search"]
    assert "наследственный гемохроматоз." in items["clinvar:9"]["search"]


# --- stage 4: the Changes tab --------------------------------------------------------------

D1, D2, D3 = sm.run_id(1), sm.run_id(2), sm.run_id(3)
ANALYSER = "a" * 64
PUBMED = "https://pubmed.ncbi.nlm.nih.gov/"
T0 = datetime(2026, 10, 3, 13, 0, tzinfo=UTC)


def _vid(pos: int) -> str:
    return f"GRCh38:chr17:{pos}:C:T"


def _at(sid: str, pos: int, **kw: Any) -> dict[str, Any]:
    return sm.statement(sid, vid=_vid(pos), **kw)


def _snapshot(
    con: duckdb.DuckDBPyConnection,
    root: Path,
    rid: str,
    clinvar: str,
    statements: list[dict[str, Any]],
    calls: list[dict[str, Any]],
    gene_calls: list[dict[str, Any]],
) -> None:
    tools = {"clinvar": clinvar, "analyser": ANALYSER}
    inputs = {"clinvar38": sm.file(f"cv{clinvar}"), "gwas": sm.file(f"gw{clinvar}")}
    sm.fill(con, rid, sm.started(rid), tools, inputs, statements, calls, gene_calls)
    export_snapshot(con, rid, root)


def _history(reports: Path) -> None:
    """A store at `reports/genomelab.duckdb` holding three runs, each exported to
    `reports/snapshots/`: from D1 to D3 a statement enters, leaves (by validity, by
    reliability, by its gene call, by going, by becoming a non-finding), moves, appears
    below the line, changes a detail, merges its ClinVar ID, and a GWAS row appears."""
    unchanged = [_at("clinvar:1", 1001), _at("clinvar:10", 1010, kind="watch")]
    unchanged_calls = [sm.call(_vid(1001)), sm.call(_vid(1010))]
    gwas = _at(
        "gwas:rs1-T:EFO_1",
        1011,
        gene="PPARG",
        rsid="rs1",
        source="gwas_catalog",
        topic="research",
    )
    early = [
        _at("clinvar:2", 1002, validity="moderate"),
        _at("clinvar:3", 1003, validity="limited"),
        _at("clinvar:4", 1004),
        _at("clinvar:5", 1005),
        _at("clinvar:7", 1007, reason="old reason"),
        _at("clinvar:8", 1008),
        _at("clinvar:9", 1009),
        sm.statement("pharmcat:CYP2D6", vid=None, gene="CYP2D6", source="pharmcat"),
    ]
    late = [
        _at("clinvar:2", 1002, validity="limited"),
        _at("clinvar:3", 1003),
        _at("clinvar:4", 1004),
        _at("clinvar:6", 1006, validity="conflicting"),
        _at("clinvar:7", 1007, reason="new reason"),
        _at("clinvar:8", 1008, kind="non_finding"),
        _at("clinvar:99", 1009, validity="moderate"),
        sm.statement("pharmcat:CYP2D6", vid=None, gene="CYP2D6", source="pharmcat"),
        gwas,
    ]
    calls = [sm.call(_vid(p)) for p in (1002, 1003, 1005, 1006, 1007, 1008, 1009)]
    reports.mkdir(parents=True)
    con = sm.store(reports / "genomelab.duckdb")
    try:
        root = reports / "snapshots"
        _snapshot(
            con,
            root,
            D1,
            "20261001",
            unchanged + early,
            [*unchanged_calls, *calls, sm.call(_vid(1004), level="medium")],
            [sm.gene_call("CYP2D6", level="medium")],
        )
        # D2: the run in between, with the same ClinVar and one reliability moved.
        _snapshot(
            con,
            root,
            D2,
            "20261001",
            unchanged + early,
            [*unchanged_calls, *calls, sm.call(_vid(1004), level="low")],
            [sm.gene_call("CYP2D6", level="medium")],
        )
        _snapshot(
            con,
            root,
            D3,
            "20261101",
            unchanged + late,
            [
                *unchanged_calls,
                *calls,
                sm.call(_vid(1004), level="low"),
                sm.call(_vid(1011)),
            ],
            [sm.gene_call("CYP2D6", level="low")],
        )
    finally:
        con.close()


def _changes(reports: Path, lab: Path) -> dict[str, Any]:
    with open_store(reports / "genomelab.duckdb") as con:
        data = payload(con, None, NO_CACHE, reports=reports, lab=lab)
    return data["changes"]


def test_changes_views_match_history_for_every_slider_pair(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    sm.refreshed(reports, lab, D2, {"clinvar:1": {"litvar": ["1"], "clinvar": []}})
    sm.refreshed(reports, lab, D3, {"clinvar:1": {"litvar": ["1", "2"], "clinvar": []}})
    root = reports / "snapshots"
    changes = _changes(reports, lab)
    assert (changes["run"], changes["default_base"]) == (D3, D2)
    assert list(changes["bases"]) == [D2, D1]
    pairs = list(product(Level, Validity))
    for base, compared in changes["bases"].items():
        found = {c.statement_id: c for c in diff(root, base, D3)}
        before_states, after_states = states(root, base), states(root, D3)
        for sid in ("clinvar:1",):  # the statements D3's refresh lists
            if sid not in found and base == D2:  # D1 lists none: a baseline
                found[sid] = Change(sid, before_states[sid], after_states[sid], (), ())
        listed = {
            sid: c
            for sid, c in found.items()
            if getattr(c.after or c.before, "source", None) != GWAS_TOOL
        }
        assert set(compared["rows"]) == set(listed)
        assert len(compared["views"]) == len(pairs)
        for rel, val in pairs:
            f = Filters(min_reliability=rel, min_validity=val)
            trusted_before = trusted_ids(root, base, f)
            trusted_after = trusted_ids(root, D3, f)
            expected = {
                sid: str(
                    kind_of(
                        c,
                        f,
                        c.before is not None
                        and c.before.statement_id in trusted_before,
                        c.after is not None and sid in trusted_after,
                    )
                )
                for sid, c in listed.items()
            }
            assert compared["views"][f"{rel}|{val}"] == expected, (base, rel, val)
        # The rows are `changes --all`'s at the default filters, with their history.
        document = report(
            reports,
            lab,
            base=base,
            run=None,
            f=Filters(),
            every=True,
            gene=None,
            marker=None,
            text=None,
        )
        assert [
            {k: v for k, v in row.items() if k not in ("history", "source")}
            for row in compared["rows"].values()
        ] == document["changes"]
        assert {row["source"] for row in compared["rows"].values()} == {
            "clinvar",
            "pharmcat",
        }
        [gwas] = [c for c in found.values() if c not in listed.values()]
        assert compared["research"] == [
            {"kind": kind, "causes": list(gwas.causes), "count": n}
            for kind, n in document["research"].items()
        ]
        assert [c.split("@")[0] for c in gwas.causes] == ["source:gwas"]
    views = changes["bases"][D2]["views"]
    # The oracle discriminates: the sliders move a statement across the line.
    assert views["medium|moderate"]["clinvar:6"] == "new"
    assert views["low|conflicting"]["clinvar:6"] == "entered"
    assert views["medium|moderate"]["clinvar:2"] == "left"
    assert views["medium|limited"]["clinvar:2"] == "level"
    assert views["medium|moderate"]["clinvar:1"] == "publications"
    assert views["medium|moderate"]["clinvar:8"] == "left"  # became a non-finding
    assert changes["bases"][D1]["views"]["medium|moderate"]["clinvar:4"] == "left"
    assert changes["bases"][D1]["views"]["low|moderate"]["clinvar:4"] == "level"
    # A row's history spans every snapshot up to the run, under its earlier id too.
    rows = changes["bases"][D2]["rows"]
    assert rows["clinvar:99"]["before"]["validity"] == "strong"
    assert rows["clinvar:99"]["history"] == [
        {
            "run": rid,
            "date": f"2026-10-0{day}",
            "kind": "finding",
            "validity": validity,
            "level": "high",
        }
        for day, rid, validity in (
            (1, D1, "strong"),
            (2, D2, "strong"),
            (3, D3, "moderate"),
        )
    ]
    assert [h["kind"] for h in rows["clinvar:8"]["history"]] == [
        "finding",
        "finding",
        "non_finding",
    ]
    assert rows["clinvar:8"]["after"]["kind"] == "non_finding"
    assert [h["level"] for h in rows["pharmcat:CYP2D6"]["history"]] == [
        "medium",
        "medium",
        "low",
    ]


def test_level_3_never_enters_change_feed(tmp_path: Path) -> None:
    """The Changes tab counts a new level-3 PharmGKB annotation under research and lists
    no row for it."""
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    reports.mkdir(parents=True)
    research = _at(
        "pharmgkb_research:655384602",
        1002,
        gene="RGS4",
        rsid="rs951439",
        validity="limited",
        source="pharmgkb_research",
        topic="research",
    )
    unchanged = [_at("clinvar:1", 1001)]
    calls = [sm.call(_vid(1001)), sm.call(_vid(1002))]
    con = sm.store(reports / "genomelab.duckdb")
    try:
        root = reports / "snapshots"
        _snapshot(con, root, D1, "20261001", unchanged, calls[:1], [])
        _snapshot(con, root, D2, "20261001", [*unchanged, research], calls, [])
    finally:
        con.close()
    compared = _changes(reports, lab)["bases"][D1]
    assert compared["rows"] == {}
    assert [(r["kind"], r["count"]) for r in compared["research"]] == [("new", 1)]


def test_a_research_score_is_counted_never_listed_on_the_changes_tab(
    tmp_path: Path,
) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    reports.mkdir(parents=True)
    research = _at(
        "pgsc_calc:research:PGS002746",
        1002,
        gene=None,
        validity="limited",
        source="pgsc_calc_research",
        topic="research",
    )
    unchanged = [_at("clinvar:1", 1001)]
    calls = [sm.call(_vid(1001)), sm.call(_vid(1002))]
    con = sm.store(reports / "genomelab.duckdb")
    try:
        root = reports / "snapshots"
        _snapshot(con, root, D1, "20261001", unchanged, calls[:1], [])
        _snapshot(con, root, D2, "20261001", [*unchanged, research], calls, [])
    finally:
        con.close()
    compared = _changes(reports, lab)["bases"][D1]
    assert compared["rows"] == {}
    assert [(r["kind"], r["count"]) for r in compared["research"]] == [("new", 1)]


def test_banner_names_failed_sources_with_their_before(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    sources = (
        SourceOutcome("clinvar", "20261001", "20261101", SourceStatus.NEW),
        SourceOutcome("gwas", "2026-10-21", None, SourceStatus.FAILED, "http:503"),
        SourceOutcome("cpic", "20261001", "20261001", SourceStatus.UNCHANGED),
        SourceOutcome("litvar", None, None, SourceStatus.FAILED, "TimeoutError"),
    )
    refresh_log.append(
        lab, RefreshRecord(T0, T0, "partial", sources, D2, None, D3, None)
    )
    changes = _changes(reports, lab)
    assert changes["failed"] == [
        {"name": "gwas", "before": "2026-10-21"},
        {"name": "litvar", "before": None},
    ]
    # Another run's refresh is not this run's banner.
    refresh_log.append(lab, RefreshRecord(T0, T0, "complete", (), D3, None, D2, None))
    assert _changes(reports, lab)["failed"][0]["name"] == "gwas"


def test_a_refresh_that_made_no_run_is_compared_with_itself(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    refresh_log.append(lab, RefreshRecord(T0, T0, "complete", (), D2, None, D3, None))
    refresh_log.append(lab, RefreshRecord(T0, T0, "complete", (), D3, None, D3, None))
    changes = _changes(reports, lab)
    assert changes["default_base"] == D3
    assert changes["bases"][D3]["rows"] == {} and list(changes["bases"])[1:] == [D2, D1]


def test_tab_says_no_base_on_first_snapshot(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    reports.mkdir()
    con = sm.store(reports / "genomelab.duckdb")
    _snapshot(
        con,
        reports / "snapshots",
        D1,
        "20261001",
        [_at("clinvar:1", 1001)],
        [sm.call(_vid(1001))],
        [],
    )
    con.close()
    changes = _changes(reports, lab)
    assert changes["run"] == D1
    assert (changes["default_base"], changes["note"], changes["bases"]) == (
        None,
        NO_BASE,
        {},
    )
    # The feed's findings carry their dates from the snapshots.
    with open_store(reports / "genomelab.duckdb") as con:
        items = payload(con, None, NO_CACHE, reports=reports, lab=lab)["items"]
    assert items["clinvar:1"]["dates"] == {
        "first_found": "2026-10-01",
        "last_changed": "2026-10-01",
        "assessment_changed": "2026-10-01",
        "on_or_before": True,
    }


def test_non_digit_pmid_never_reaches_the_payload(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    hostile = [
        "123",
        "45x",
        "javascript:alert(1)",
        "١٢٣",
        " 7",
        "8\n",
        "//evil.example",
    ]
    sm.refreshed(reports, lab, D2, {"clinvar:1": {"litvar": [], "clinvar": []}})
    sm.refreshed(
        reports, lab, D3, {"clinvar:1": {"litvar": hostile, "clinvar": hostile[::-1]}}
    )
    with open_store(reports / "genomelab.duckdb") as con:
        data = payload(con, None, NO_CACHE, reports=reports, lab=lab)
    row = data["changes"]["bases"][D2]["rows"]["clinvar:1"]
    assert row["publications"] == {"new": ["123"], "attached_by_lab": ["123"]}
    pmids = [
        pmid
        for compared in data["changes"]["bases"].values()
        for listed in compared["rows"].values()
        for found in listed["publications"].values()
        for pmid in found
    ]
    assert pmids and all(re.fullmatch(r"[0-9]+", pmid) for pmid in pmids)
    page = render(data)
    for text in ("45x", "javascript:alert", "١٢٣", "//evil.example"):
        assert text not in page


def test_page_loads_nothing_from_the_network(tmp_path: Path) -> None:
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    sm.refreshed(reports, lab, D2, {"clinvar:1": {"litvar": ["1"], "clinvar": []}})
    sm.refreshed(
        reports, lab, D3, {"clinvar:1": {"litvar": ["1", "2"], "clinvar": ["3"]}}
    )
    with open_store(reports / "genomelab.duckdb") as con:
        data = payload(con, None, NO_CACHE, reports=reports, lab=lab)
    assert data["changes"]["bases"][D2]["rows"]["clinvar:1"]["publications"]["new"]
    page = render(data)
    js, css = _resource("dashboard.js"), _resource("dashboard.css")
    [csp] = re.findall(r'http-equiv="Content-Security-Policy" content="([^"]*)"', page)
    assert csp == (
        f"default-src 'none'; script-src 'sha256-{_sha256(js)}';"
        f" style-src 'sha256-{_sha256(css)}'"
    )
    # Outside the data (which may hold any text, inside application/json), the page's
    # only URL is the script's PubMed prefix.
    outside = re.sub(
        r'<script type="application/json" id="data">.*?</script>',
        "",
        page,
        flags=re.DOTALL,
    )
    assert re.findall(r"https?://[^\"'`\s]*", outside) == [PUBMED]
    assert re.findall(r"https?://[^\"'`\s]*", js) == [PUBMED]
    [name] = re.findall(rf'const (\w+) = "{re.escape(PUBMED)}";', js)
    uses = re.findall(rf"(.{{0,16}})\b{name}\b(.{{0,3}})", js)
    assert len(uses) >= 2
    for before, after in uses:
        assert before.endswith("const ") or (
            before.endswith(".href = ") and after == " + "
        )


def test_write_dashboard_writes_the_page_with_the_tab(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = Paths(
        repo=tmp_path,
        raw=tmp_path / "raw",
        reports=tmp_path / "reports",
        lab=tmp_path / "lab",
    )
    _history(paths.reports)
    out = dashboard.write_dashboard(paths)
    assert out == paths.reports / "dashboard.html"
    data = _page_data(out.read_text(encoding="utf-8"))
    assert data["run_id"] == D3
    assert data["changes"]["default_base"] == D2
    assert data["changes"]["bases"][D2]["views"]["medium|moderate"]["clinvar:3"] == (
        "entered"
    )
    assert not out.with_name(out.name + ".tmp").exists()
    # An older run: its changes compare with the snapshots before it.
    out = dashboard.write_dashboard(paths, D2)
    changes = _page_data(out.read_text(encoding="utf-8"))["changes"]
    assert (changes["run"], list(changes["bases"])) == (D2, [D1])


def _write_check(lab: Path, ran_at: datetime, *sources: dict[str, Any]) -> None:
    (lab / "refresh").mkdir(parents=True, exist_ok=True)
    (lab / "refresh" / "check.json").write_text(
        json.dumps({"ran_at": ran_at.isoformat(), "notified": {}, "sources": sources})
    )


def test_tab_carries_the_check_without_its_age(tmp_path: Path) -> None:
    """The payload holds `ran_at`, what is newer and what failed, never an age: the page
    outlives the check, and the script computes the age when it is viewed."""
    reports, lab = tmp_path / "reports", tmp_path / "lab"
    _history(reports)
    ran = datetime(2026, 8, 1, 10, tzinfo=UTC)  # months before the page is generated
    newer = {
        "name": "pharmcat",
        "kind": "pinned",
        "upstream": "3.5.0",
        "local": "3.4.0",
        "status": "newer",
        "error": None,
    }
    broken = {**newer, "name": "gwas", "status": "failed", "upstream": None}
    _write_check(lab, ran, newer, broken)
    block = _changes(reports, lab)["check"]
    assert block == {
        "ran_at": ran.isoformat(),
        "newer": [
            {
                "name": "pharmcat",
                "kind": "pinned",
                "local": "3.4.0",
                "upstream": "3.5.0",
            }
        ],
        "failed": ["gwas"],
        "stale_days": 35,
    }
    # No check at all: never ran, not "nothing newer".
    assert _changes(reports, tmp_path / "other")["check"] == {
        "ran_at": None,
        "newer": [],
        "failed": [],
        "stale_days": 35,
    }
    # The page is the data and a script that reads the age off the clock.
    with open_store(reports / "genomelab.duckdb") as con:
        data = payload(con, None, NO_CACHE, reports=reports, lab=lab)
    page = render(data)
    assert ran.isoformat() in page
    js = _resource("dashboard.js")
    assert "Date.now()" in js and "Date.parse(info.ran_at)" in js
    assert "info.stale_days" in js and "info.failed" in js and "info.newer" in js


# The page's script under a minimal DOM: each element's text and links, by id, after it
# rendered `data` (argv: the script, the data).
_DOM = r"""
const fs = require("fs");
const vm = require("vm");
class Text {
  constructor(t) { this.text = String(t); }
  get textContent() { return this.text; }
  links() { return []; }
}
class Node {
  constructor(tag) {
    Object.assign(this, { tag, kids: [], style: {}, className: "", hidden: false, value: "", href: "" });
  }
  append(...xs) { for (const x of xs) this.kids.push(typeof x === "string" ? new Text(x) : x); }
  replaceChildren(...xs) { this.kids = []; this.append(...xs); }
  set textContent(t) { this.kids = [new Text(t)]; }
  get textContent() { return this.kids.map((k) => k.textContent).join("\n"); }
  links() { return [...(this.tag === "a" ? [this.href] : []), ...this.kids.flatMap((k) => k.links())]; }
  setAttribute() {}
  addEventListener() {}
}
const byId = {};
const document = {
  createElement: (tag) => new Node(tag),
  getElementById: (id) => (byId[id] ??= new Node("div")),
};
document.getElementById("data").textContent = fs.readFileSync(process.argv[3], "utf8");
vm.runInNewContext(fs.readFileSync(process.argv[2], "utf8"), { document });
const out = {};
for (const [id, node] of Object.entries(byId)) out[id] = { text: node.textContent, links: node.links() };
console.log(JSON.stringify(out));
"""


def _rendered(data: dict[str, Any], tmp_path: Path) -> dict[str, dict[str, Any]]:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not installed")
    (tmp_path / "dom.js").write_text(_DOM)
    (tmp_path / "page.js").write_text(_resource("dashboard.js"))
    (tmp_path / "data.json").write_text(json.dumps(data, ensure_ascii=False))
    done = subprocess.run(
        [node, "dom.js", "page.js", "data.json"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(done.stdout)


IMAGE = "/lab/reports/work/repeats/k/rev.HTT.svg"
RESULT_ROWS = [
    Row(
        "yleaf:haplogroup",
        gene=None,
        level=None,
        topic="ancestry",
        claim={"type": "y_haplogroup", "haplogroup": "R1b-M269"},
    ),
    Row(
        "haplogrep3:haplogroup",
        gene=None,
        level=None,
        topic="ancestry",
        claim={"type": "mito_haplogroup", "haplogroup": "H1a"},
    ),
    Row(
        "pgsc_calc:PGS000018",
        gene=None,
        level=None,
        claim={"type": "pgs", "condition": "Coronary artery disease"},
        effect={
            "relative": None,
            "absolute": None,
            "score": {
                "percentile": 91.5,
                "z": 1.37,
                "population": "EUR",
                "match_rate": 0.98,
            },
        },
    ),
    Row(
        "predicted:GRCh38:chr17:1050:C:T:ENST1",
        pos=1050,
        claim={
            "type": "predicted_variant",
            "classification": "predicted, not classified",
            "consequence": "stop_gained",
            "hgvsc": "ENST1:c.121C>T",
            "hgvsp": "ENSP1:p.Arg41Ter",
        },
    ),
    Row(
        "repeats:HTT",
        gene="HTT",
        level=None,
        topic="carrier",
        gene_call={**GENE_CALL, "diplotype": "CAG 17/20", "phenotype": None},
        claim={"type": "repeat_expansion", "class": "normal", "image": IMAGE},
    ),
]


def test_rows_show_each_claim_types_result(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A haplogroup, a polygenic score, a predicted finding and a repeat locus show their
    result in the row, the search string holds it, and the REViewer image is linked."""
    seed(con, RESULT_ROWS)
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    feed = _rendered(data, tmp_path)["feed"]
    for shown in (
        "Y-гаплогруппа R1b-M269",
        "митохондриальная гаплогруппа H1a",
        "Coronary artery disease: 91.5-й перцентиль",
        "z 1.37 (популяция EUR",
        "stop_gained, ENSP1:p.Arg41Ter",
        "HTT: CAG 17/20, норма",
    ):
        assert shown in feed["text"], shown
    assert f"file://{IMAGE}" in feed["links"]
    search = {sid: item["search"] for sid, item in data["items"].items()}
    assert "r1b-m269" in search["yleaf:haplogroup"]
    assert "p.arg41ter" in search["predicted:GRCh38:chr17:1050:C:T:ENST1"]
    assert "cag 17/20" in search["repeats:HTT"] and "normal" in search["repeats:HTT"]


WEAK = "слабая предсказательная сила"


def test_pgs_headline_names_weak_power(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A weak score names its weakness beside the percentile, in the row's result."""
    (pgs,) = (r for r in RESULT_ROWS if r.id == "pgsc_calc:PGS000018")
    weak = replace(
        pgs,
        id="pgsc_calc:PGS004885",
        claim={"type": "pgs", "condition": "depression", "weak": True},
    )
    seed(con, [weak])
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    lines = _rendered(data, tmp_path)["feed"]["text"].split("\n")
    assert f"depression: 91.5-й перцентиль, {WEAK}" in lines


def test_pgs_headline_without_weak_is_unchanged(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A score from a store before slice 3g (no `weak` key) shows today's headline."""
    seed(con, [r for r in RESULT_ROWS if r.id == "pgsc_calc:PGS000018"])
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    text = _rendered(data, tmp_path)["feed"]["text"]
    assert "Coronary artery disease: 91.5-й перцентиль" in text.split("\n")
    assert WEAK not in text


def test_not_covered_is_not_shown_as_genotype_unknown(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A not-covered condition is unassessed: the heading allows for it, and its line
    carries its stored reason, which the skill keys on (`not covered:`)."""
    reason = (
        "not covered: no score meets the criteria (evaluated only in its own cohort)"
    )
    row = Row(
        "pgsc_calc:not-covered:migraine", kind="unassessed", gene=None, level=None
    )
    claim = {"type": "pgs_not_covered", "condition": "migraine"}
    seed(con, [replace(row, validity="limited", reason=reason, claim=claim)])
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    shown = _rendered(data, tmp_path)["unassessed"]["text"].split("\n")
    assert "генотип не определён или условие не покрыто, см. причину" in shown[0]
    assert f"без гена: 1; {reason} (1)" in shown


def test_pgs_line_shows_decile(con: duckdb.DuckDBPyConnection, tmp_path: Path) -> None:
    """A score stored with its decile shows the decile, not the percentile, and its
    score line adds `Z_norm2` when present (slice 3f, "PGS")."""
    score = {
        "percentile": 91.5,
        "z": 1.37,
        "population": "EUR",
        "match_rate": 0.98,
        "z_norm2": 1.2,
        "decile": 10,
    }
    row = Row(
        "pgsc_calc:PGS000018",
        gene=None,
        level=None,
        claim={"type": "pgs", "condition": "Coronary artery disease"},
        effect={"relative": None, "absolute": None, "score": score},
    )
    seed(con, [row])
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    text = _rendered(data, tmp_path)["feed"]["text"]
    assert "Coronary artery disease: 10-й дециль" in text
    assert "перцентиль" not in text
    assert "z 1.37, Z_norm2 1.2 (популяция EUR" in text


def test_carrier_topic_lists_what_is_not_covered(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """The carrier topic names what no run looks at and each gene whose mechanism this
    run did not assess, in Russian; SMN1 leaves the list once its copy number is."""
    smn = {"type": "smn_copy_number", "smn1": 2}
    seed(con, [Row("smn:SMN1", gene="SMN1", level=None, topic="carrier", claim=smn)])
    con.execute("UPDATE statement SET source = 'smn' WHERE statement_id = 'smn:SMN1'")
    data = payload(con, None, NO_CACHE, lang="ru")
    carrier = _rendered({**data, "default": "low|conflicting|carrier"}, tmp_path)
    shown = carrier["not-covered"]["text"]
    assert "структурные варианты" in shown and "не значит «не носитель»" in shown
    for symbol in ("HBA1", "HBA2", "CYP21A2", "F8", "GBA1", "DMD", "FMR1"):
        assert f"{symbol}: " in shown, symbol
    assert "SMN1: " not in shown
    every = _rendered({**data, "default": "low|conflicting|"}, tmp_path)
    assert every["not-covered"]["text"] == ""


def test_russian_not_covered_lists_match_the_carrier_ones() -> None:
    """Each Russian line keeps its English line's gene symbols and numbers, so an English
    line that changes without its translation fails here (3d's copy-number line once did)."""
    ru = dashboard.messages("ru")
    assert len(ru["carrier_not_covered"]) == len(NOT_COVERED)
    assert set(ru["carrier_mechanisms"]) == set(load_mechanisms())
    mechanisms = load_mechanisms()
    pairs = [*zip(NOT_COVERED, ru["carrier_not_covered"])]
    pairs += [
        (mechanisms[gene], text) for gene, text in ru["carrier_mechanisms"].items()
    ]
    for en, translated in pairs:
        tokens = set(re.findall(r"[A-Z][A-Z0-9]{2,}|\b\d+\b", en))
        assert tokens <= set(re.findall(r"[A-Z][A-Z0-9]{2,}|\b\d+\b", translated)), en


def test_english_not_covered_lists_are_the_carrier_ones(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    data = payload(con, None, NO_CACHE)
    assert data["not_covered"] == {"all": list(NOT_COVERED), "genes": load_mechanisms()}
    carrier = _rendered({**data, "default": "low|conflicting|carrier"}, tmp_path)
    shown = carrier["not-covered"]["text"]
    assert "does not mean “not a carrier”" in shown
    assert "Genes whose main mechanism is not assessed" in shown


def test_dashboard_names_the_reference_cause() -> None:
    """A change the masked reference moved is grouped under its own label."""
    js = _resource("dashboard.js")
    assert 'reference: "cause_reference"' in js
    assert dashboard.messages("ru")["cause_reference"] == "референсный геном"
    assert dashboard.messages("en")["cause_reference"] == "reference genome"


def _texts(message: Any) -> list[str]:
    """A message's texts: itself, each of its plural forms, or none for a table."""
    if isinstance(message, str):
        return [message]
    if isinstance(message, dict) and "other" in message:
        return [str(form) for form in message.values()]
    return []


def test_catalog_has_every_message_in_every_language() -> None:
    """Each message has exactly the languages `genokeeper.lang` names, with the same
    placeholders in each, and the script and page use every message and no other."""
    text = _resource("messages.json")
    catalog: dict[str, dict[str, Any]] = json.loads(text)
    en, ru = dashboard.messages("en"), dashboard.messages("ru")
    for key, entry in catalog.items():
        assert tuple(entry) == LANGUAGES, key
        holes = {
            lang: {h for form in _texts(own[key]) for h in re.findall(r"{(\w+)}", form)}
            for lang, own in (("en", en), ("ru", ru))
        }
        assert holes["en"] == holes["ru"], key
    js = _resource("dashboard.js")
    used = set(re.findall(r'"(\w+)"', js)) | set(re.findall(r"\bM\.(\w+)", js))
    page = {m[1] for m in re.finditer(r"\{\{(\w+)\}\}", _resource("page.html"))}
    carrier = {"carrier_not_covered", "carrier_mechanisms"}
    assert set(catalog) - carrier <= used | page
    called = set(re.findall(r'\bt\("(\w+)"', js)) | set(re.findall(r"\bM\.(\w+)", js))
    assert called <= set(catalog)


# --- slice 3f: the "Excluded by screening" block --------------------------------------------

RULED_OUT = (
    "coeliac disease practically excluded: no HLA-DQ2 or -DQ8;"
    " the tag SNPs are validated in Europeans"
)
CONFIRM = "HLA-DQ typing in a laboratory"
EXCLUSION_ROWS = [
    Row(
        "hla:coeliac",
        kind="non_finding",
        gene="HLA-DQ",
        level=None,
        claim={"type": "hla_coeliac", "confirm": CONFIRM},
        reason=RULED_OUT,
    ),
    Row(
        "hla:HLA-B*27",
        kind="non_finding",
        gene="HLA-B",
        level=None,
        claim={"type": "hla_allele", "confirm": "HLA-B typing"},
        reason="not carried",
    ),
    Row(
        "hla:HLA-B*51:association",
        kind="unassessed",
        gene="HLA-B",
        level=None,
        claim={"type": "hla_association", "confirm": "HLA-B typing"},
        reason="the callers do not agree",
    ),
    Row(
        "hla:HLA-C*06:02:association",
        gene="HLA-C",
        level=None,
        claim={"type": "hla_association", "confirm": "HLA-C typing"},
        reason="carried",
    ),
]
# Flagged: the rule-out, and (to prove the kind filter) the unassessed and finding rows.
FLAGGED = ("hla:coeliac", "hla:HLA-B*51:association", "hla:HLA-C*06:02:association")


def _flag(con: duckdb.DuckDBPyConnection, *ids: str) -> None:
    for sid in ids:
        con.execute(
            "UPDATE statement SET exclusion = true WHERE statement_id = ?", [sid]
        )


def test_excluded_lists_only_flagged_non_findings(
    con: duckdb.DuckDBPyConnection,
) -> None:
    seed(con, EXCLUSION_ROWS)
    _flag(con, *FLAGGED)
    expected = [
        {
            "id": "hla:coeliac",
            "gene": "HLA-DQ",
            "wording": RULED_OUT,
            "confirm": CONFIRM,
        }
    ]
    assert dashboard._excluded(con, RUN) == expected
    assert payload(con, None, NO_CACHE)["excluded"] == expected


def test_excluded_block_rendered(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """The page holds the block's container under the feed, the script names each item's
    gene, reason and confirming test, and the CSP still pins the script's bytes."""
    seed(con, EXCLUSION_ROWS)
    _flag(con, *FLAGGED)
    data = payload(con, None, NO_CACHE, lang="ru")
    page = render(data)
    assert page.index('id="feed"') < page.index('id="excluded"')
    [script] = re.findall(r"<script>(.*?)</script>", page, re.DOTALL)
    [csp] = re.findall(r'http-equiv="Content-Security-Policy" content="([^"]*)"', page)
    assert f"script-src 'sha256-{_sha256(script)}';" in csp
    assert _page_data(page)["excluded"][0]["gene"] == "HLA-DQ"
    shown = _rendered(data, tmp_path)["excluded"]["text"]
    assert "Исключено скринингом" in shown
    assert f"HLA-DQ: {RULED_OUT}; подтвердить: {CONFIRM}" in shown
    assert "HLA-B" not in shown and "HLA-C" not in shown


REPEAT = {"type": "repeat_expansion"}
SCREENED_ROWS = [
    Row(
        "expansionhunter:HTT",
        kind="watch",
        gene="HTT",
        level=None,
        claim=REPEAT,
        reason="screening result, needs confirmation: borderline, indeterminate",
        gene_call={"diplotype": "CAG 17/25", "phenotype": "borderline"},
    ),
    Row(
        "expansionhunter:DMPK",
        kind="non_finding",
        gene="DMPK",
        level=None,
        claim=REPEAT,
        gene_call={"diplotype": "CAG 5/12", "phenotype": "normal"},
    ),
    Row(
        "smn:SMN1",
        kind="non_finding",
        gene="SMN1",
        level=None,
        claim={"type": "smn_copy_number"},
        gene_call={"diplotype": "SMN1 3, SMN2 2", "phenotype": "not a carrier"},
    ),
    # Not listed: a flagged normal (the "Excluded by screening" block's), a stage-3
    # finding and unassessed claim (the feed's and its counts'), CFTR's poly-T/TG (no
    # cystic-fibrosis test), and a gene-level non-finding of another source.
    Row(
        "expansionhunter:CNBP",
        kind="non_finding",
        gene="CNBP",
        level=None,
        claim={"type": "repeat_expansion", "confirm": "PCR sizing"},
        reason="normal: 15 and 18 CAGG, each allele read across in full",
    ),
    Row("expansionhunter:ATXN1", gene="ATXN1", level=None, claim=REPEAT),
    Row(
        "expansionhunter:ATXN2",
        kind="unassessed",
        gene="ATXN2",
        level=None,
        claim=REPEAT,
    ),
    Row(
        "expansionhunter:CFTR-polyT-TG",
        kind="non_finding",
        gene="CFTR",
        level=None,
        claim=REPEAT,
        reason="no 5T allele: T 7/9, TG 11/12",
    ),
    Row("hhv6:ici", kind="non_finding", gene="U94", level=None, claim={"type": "x"}),
]


def test_screened_lists_stage3_watch_items_and_non_findings_apart(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A stage-3 watch item and a stage-3 non-finding reach neither the feed nor its
    counts: the page lists them in a block of their own, a watch apart from a normal."""
    seed(con, SCREENED_ROWS)
    _flag(con, "expansionhunter:CNBP")
    for prefix in ("expansionhunter", "smn", "hhv6"):
        con.execute(
            "UPDATE statement SET source = ? WHERE statement_id LIKE ?",
            [prefix, f"{prefix}:%"],
        )
    data = payload(con, None, NO_CACHE, lang="ru")
    page = render(data)
    assert page.index('id="excluded"') < page.index('id="screened"')
    assert data["screened"] == {
        "watch": [
            {
                "id": "expansionhunter:HTT",
                "label": "HTT CAG 17/25",
                "wording": SCREENED_ROWS[0].reason,
                "explanation": None,
            }
        ],
        "normal": [
            {"id": "expansionhunter:DMPK", "label": "DMPK CAG 5/12", "wording": None},
            {"id": "smn:SMN1", "label": "SMN1 3, SMN2 2", "wording": None},
        ],
    }
    shown = _rendered(data, tmp_path)["screened"]["text"]
    assert "Под наблюдением: ни норма, ни находка" in shown
    assert f"HTT CAG 17/25: {SCREENED_ROWS[0].reason}" in shown
    assert "Скрининг: не находка" in shown and "DMPK CAG 5/12" in shown
    assert "not a carrier" not in shown
    assert shown.index("HTT") < shown.index("Скрининг") < shown.index("DMPK")
    for name in ("CNBP", "ATXN1", "ATXN2", "CFTR", "U94"):
        assert name not in shown
    other = connect(tmp_path / "other.duckdb")
    seed(other, ROWS)
    empty = payload(other, None, NO_CACHE)
    assert empty["screened"] == {"watch": [], "normal": []}
    assert _rendered(empty, tmp_path)["screened"]["text"].count("none in this run") == 2


def test_a_watch_items_cached_explanation_shown(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A stage-3 watch item `explain` wrote an explanation for (BR-025) shows it in the
    stage-3 block; a non-finding has none."""
    seed(con, SCREENED_ROWS)
    con.execute(
        "UPDATE statement SET source = 'expansionhunter'"
        " WHERE statement_id LIKE 'expansionhunter:%'"
    )
    (watch,) = watch_items(con, RUN)
    cache = tmp_path / "explanations"

    def ask(prompt: str) -> str:
        answer = {"title": "Пограничный размер.", "text": "Не норма и не находка.",
                  "cites": [watch["id"]]}  # fmt: skip
        return json.dumps(answer, ensure_ascii=False)

    assert explain([watch], None, cache, ask, lang="ru")["written"] == 1
    data = payload(con, None, cache, None, lang="ru")
    (listed,) = data["screened"]["watch"]
    assert listed["explanation"]["text"] == "Не норма и не находка."
    assert all("explanation" not in x for x in data["screened"]["normal"])
    shown = _rendered(data, tmp_path)["screened"]["text"]
    assert "Не норма и не находка." in shown
    assert payload(con, None, NO_CACHE)["screened"]["watch"][0]["explanation"] is None


def test_excluded_empty_run(con: duckdb.DuckDBPyConnection, tmp_path: Path) -> None:
    data = payload(con, None, NO_CACHE, lang="ru")
    assert data["excluded"] == []
    shown = _rendered(data, tmp_path)["excluded"]["text"]
    assert "Исключено скринингом" in shown and "нет исключений в этом запуске" in shown


# --- slice 3f: a headline for every claim type ----------------------------------------------

SLICE3F_ROWS = [
    Row(
        "hla:HLA-B*27",
        gene="HLA-B",
        level=None,
        claim={
            "type": "hla_allele",
            "gene": "HLA-B",
            "allele": "*27",
            "condition": "ankylosing spondylitis",
        },
    ),
    Row(
        "hla:HLA-C*06:02:association",
        gene="HLA-C",
        level=None,
        claim={
            "type": "hla_association",
            "gene": "HLA-C",
            "allele": "*06:02",
            "condition": "psoriasis",
        },
    ),
    Row(
        "hla:coeliac",
        gene="HLA-DQ",
        level=None,
        claim={
            "type": "hla_coeliac",
            "t1k": {"dq25": 0, "dq22": 0, "dq8": 1, "half": False},
            "confirm": CONFIRM,
        },
    ),
    Row(
        "hla:table",
        gene=None,
        level=None,
        claim={"type": "hla_table", "genes": [{"gene": "HLA-A"}, {"gene": "HLA-B"}]},
    ),
    Row(
        "mitomap:3243A>G",
        gene="MT-TL1",
        chrom="chrM",
        pos=3243,
        claim={
            "type": "mitomap",
            "allele": "m.3243A>G",
            "condition": "MELAS",
            "level": 0.3,
        },
    ),
    Row(
        "hhv6:ici",
        gene=None,
        level=None,
        claim={"type": "hhv6_ici", "species": "HHV-6B", "ratio": 1.02},
    ),
    Row(
        "roh:summary",
        gene=None,
        level=None,
        claim={"type": "roh_summary", "total_long": 12_020_000, "n_long": 2},
    ),
    Row(
        "bloodgroup:FY",
        gene=None,
        level=None,
        claim={"type": "blood_group", "system": "FY", "phenotype": "Fy(a-b+)"},
    ),
    Row(
        "ancestry:similarity",
        gene=None,
        level=None,
        claim={"type": "ancestry_similarity", "population": "EUR"},
    ),
    Row(
        "ancestry:admixture",
        gene=None,
        level=None,
        claim={
            "type": "ancestry_admixture",
            "proportions": {"AFR": 0.03, "EUR": 0.85, "MID": 0.12},
        },
    ),
    Row(
        "ancestry:neanderthal",
        gene=None,
        level=None,
        claim={"type": "ancestry_neanderthal", "share": 0.011324},
    ),
    Row(
        "trait:lactase",
        gene="LCT",
        level=None,
        claim={"type": "trait", "trait": "lactase", "result": "likely"},
    ),
    Row(
        "trait:ACTN3",
        gene="ACTN3",
        level=None,
        claim={"type": "trait", "trait": "ACTN3", "result": 2},
    ),
    Row(
        "trait:TAS2R38",
        gene="TAS2R38",
        level=None,
        claim={"type": "trait", "trait": "TAS2R38", "result": "PAV/AVI"},
    ),
    Row(
        "trait:MC1R",
        gene="MC1R",
        level=None,
        claim={"type": "trait", "trait": "MC1R", "result": "one strong"},
    ),
    Row(
        "trait:HERC2",
        gene="HERC2",
        level=None,
        claim={"type": "trait", "trait": "HERC2", "result": 2},
    ),
]


def test_headline_for_every_3f_claim_type(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """Each slice 3f claim type has its own result in the row; none falls back to the
    claim's type or to a condition alone."""
    seed(con, SLICE3F_ROWS)
    data = {**payload(con, None, NO_CACHE, lang="ru"), "default": "low|conflicting|"}
    feed = _rendered(data, tmp_path)["feed"]["text"]
    for shown in (
        "HLA-B: аллель HLA-B*27, ankylosing spondylitis",
        "HLA-C: исследовательская ассоциация HLA-C*06:02, psoriasis",
        "HLA-DQ: целиакия: найден DQ8 — результат скрининга",
        "таблица HLA: 2 гена",
        "MT-TL1: MITOMAP m.3243A>G, MELAS, уровень 30%",
        "встроенный HHV-6B, отношение глубин 1.02",
        "участки гомозиготности: 12.02 Мб в 2 участках",
        "группа крови FY: Fy(a-b+)",
        "ближе всего к популяции EUR",
        "доли предков: EUR 85%, MID 12%, AFR 3%",
        "доля неандертальской ДНК 1.13%",
        "LCT: сохранение лактазы у взрослых: вероятно",
        "ACTN3: R577X, белок быстрых мышечных волокон: копий аллеля эффекта 2",
        "TAS2R38: восприятие горького вкуса: PAV/AVI",
        "MC1R: варианты рыжих волос: один сильный вариант",
        "HERC2: цвет глаз: копий аллеля эффекта 2",
    ):
        assert shown in feed, shown
    for row in SLICE3F_ROWS:
        assert row.claim is not None
        assert f": {row.claim['type']}\n" not in feed, row.id


def test_mc1r_two_strong_names_copies() -> None:
    """MC1R's "two strong" counts copies (one strong variant on both is two), so its label
    names copies, not variants, in each language."""
    en, ru = dashboard.messages("en"), dashboard.messages("ru")
    assert en["trait_results"]["two strong"] == "two copies with a strong variant"
    assert ru["trait_results"]["two strong"] == "две копии с сильным вариантом"


# --- the page's language --------------------------------------------------------------------


def test_english_page_names_each_view(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """The English page names the levels, the topics, a percentile, a decile, a gene call,
    a genotype's reliability, a weak score and the empty states, in English."""
    pgs = next(r for r in RESULT_ROWS if r.id == "pgsc_calc:PGS000018")
    weak = replace(
        pgs,
        id="pgsc_calc:PGS004885",
        claim={"type": "pgs", "condition": "depression", "weak": True},
    )
    assert pgs.effect is not None
    decile = replace(
        pgs,
        id="pgsc_calc:PGS000001",
        claim={"type": "pgs", "condition": "asthma"},
        effect={**pgs.effect, "score": {**pgs.effect["score"], "decile": 1}},
    )
    seed(con, [*RESULT_ROWS, *SLICE3F_ROWS, weak, decile])
    data = payload(con, None, NO_CACHE)
    assert data["lang"] == "en"
    shown = _rendered({**data, "default": "low|conflicting|"}, tmp_path)
    feed = shown["feed"]["text"]
    for text in (
        "Coronary artery disease: 91.5th percentile",
        "depression: 91.5th percentile, weak predictive power",
        "asthma: 1st decile",
        "CYP2D6: *1/*2, Intermediate Metabolizer",
        "reliability high, calls: BCyrius *1/*2; Aldy *1/*2",
        "0/1, reliability high, chance of error 1.0e-4 (measured on GIAB reference samples, DNBSEQ-T7)",
        "Y haplogroup R1b-M269",
        "HTT: CAG 17/20, normal",
        "HLA table: 2 genes",
        "runs of homozygosity: 12.02 Mb in 2 runs",
        "Strength of evidence: \nstrong",
        "\nmed.\n",
    ):
        assert text in feed, text
    assert not re.search("[а-яё]", feed, re.IGNORECASE)
    assert shown["reliability-label"]["text"] == "Genotype reliability: at least “low”"
    assert shown["validity-label"]["text"] == (
        "Strength of evidence: at least “conflicting”"
    )
    chips = shown["chips"]["text"].split("\n")
    assert chips[:4] == ["All topics", "Risks", "Carrier status", "Drugs"]
    assert (
        shown["excluded"]["text"] == "Excluded by screening\nno exclusions in this run"
    )
    assert shown["count"]["text"].endswith(" findings · actionable first")
    hidden = _rendered(data, tmp_path)["hidden"]["text"]
    assert "with a genotype below “medium”" in hidden
    assert hidden.endswith("Move the sliders to see them.")


def test_english_count_is_singular_for_one(
    tmp_path: Path, con: duckdb.DuckDBPyConnection
) -> None:
    con.execute("DELETE FROM statement WHERE statement_id <> 'clinvar:1'")
    shown = _rendered(payload(con, None, NO_CACHE), tmp_path)
    assert shown["count"]["text"] == "Showing 1 of 1 finding · actionable first"


@pytest.mark.parametrize(
    ("env", "flag", "lang", "tab"),
    [
        (None, None, "en", "Feed"),
        ("ru", None, "ru", "Лента"),
        ("ru", "en", "en", "Feed"),
        (None, "ru", "ru", "Лента"),
    ],
)
def test_page_language_is_flag_then_env_then_english(
    con: duckdb.DuckDBPyConnection,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    env: str | None,
    flag: str | None,
    lang: str,
    tab: str,
) -> None:
    con.close()
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    if env is None:
        monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    else:
        monkeypatch.setenv("GENOKEEPER_LANG", env)
    out = tmp_path / "dashboard.html"
    db = str(tmp_path / "genomelab.duckdb")
    args = ["dashboard", "--db", db, "--cache", str(NO_CACHE), "--out", str(out)]
    assert main([*args, *(["--lang", flag] if flag else [])]) == 0
    page = out.read_text(encoding="utf-8")
    assert page.startswith(f'<!doctype html>\n<html lang="{lang}">')
    assert f'id="tab-feed">{tab}</button>' in page
    assert _page_data(page)["lang"] == lang


def test_unknown_language_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GENOKEEPER_LANG", "fr")
    out = tmp_path / "dashboard.html"
    assert (
        main(["dashboard", "--db", str(tmp_path / "x.duckdb"), "--out", str(out)]) == 2
    )
    [err] = capsys.readouterr().err.splitlines()
    assert err.startswith("genokeeper dashboard: unknown language 'fr'")
    assert not out.exists()


def test_the_page_shows_the_explanations_of_its_own_language(
    con: duckdb.DuckDBPyConnection, tmp_path: Path
) -> None:
    """A Russian page shows the Russian cache's explanations, an English page none of them:
    `payload` looks each finding up in the page's language."""
    ghr = write_medlineplus(
        tmp_path / "lab" / "medlineplus" / "20261020" / MEDLINEPLUS_FILE,
        gene("BRCA1", "The BRCA1 gene makes a protein."),
        condition("hboc", "Hereditary breast ovarian cancer", "A condition."),
    )
    medlineplus = load_medlineplus(ghr)
    cache = tmp_path / "explanations"
    findings = list(payload(con, None, NO_CACHE)["items"].values())

    def ask(prompt: str) -> str:
        own = re.search(r'"id": "([^"]+)"', prompt[prompt.index("\n{\n") :])
        assert own is not None
        return json.dumps({"title": "Заголовок", "text": "Текст.", "cites": [own[1]]})

    explain(findings, medlineplus, cache, ask, lang="ru")
    ru = payload(con, None, cache, medlineplus, lang="ru")["items"]["clinvar:1"]
    en = payload(con, None, cache, medlineplus, lang="en")["items"]["clinvar:1"]
    assert (ru["explanation"]["text"], en["explanation"]) == ("Текст.", None)
