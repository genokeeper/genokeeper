"""`python -m genokeeper.demo package OUT`: what the HG002 demo publishes, written from the
demo root (`demo_paths(lab)`) into `OUT`.

The store is copied into a scratch root, every text value of every table (its VARCHAR,
JSON and VARCHAR[] columns) has its absolute paths rewritten (`scrub`: the demo root to
`<demo>`, the lab to `<lab>`, the home directory to `~`), and the shipped store is a fresh
database written from that copy (`COPY FROM DATABASE`): DuckDB keeps a replaced string in the
file's bytes after an UPDATE, and a fresh copy does not. From the scratch copy, with its
snapshots and a lab holding `DEMO_MARKER` and the demo's refresh-log record, come the
dashboard (English, the frame in its banner), the calibration page (each table's report
with the frame's calibration paragraph on top) and the landing page, its numbers from the
store. `MANIFEST` names each file's sha256.

The scan is the control: the adversary is the maintainer's own machine leaking into a public
file by mistake; the anchor is the outputs as written, read back. It covers every file
`package` ships, `MANIFEST` included, and every column of every table of the shipped store
(`_texts`). `scan_bytes` looks in each for the home directory, the lab and the demo root, as
exact strings, and `scan_private` runs each pattern of the private-ID file the genome-data
guard reads (`private_id_file`) through `grep -E`, the guard's engine. Any hit is a `Leak`, and
nothing is renamed from `OUT/.tmp-*` to its final name. The bare user name is not a needle: it
can match a condition or a drug name, and the home path covers it where it matters.
"""

import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from importlib.resources import files
from math import log1p
from pathlib import Path
from typing import Any

import duckdb

import genokeeper
from genokeeper import changes, refresh_log
from genokeeper.calibrate.tables import REPORT, TABLES, folder
from genokeeper.checksum import sha256_file
from genokeeper.config import _REPO, DEMO_MARKER, Paths
from genokeeper.dashboard import _banner, _sha256, write_dashboard
from genokeeper.demo import (
    MANIFEST,
    REFRESH_LOG_FILE,
    SHIPPED,
    STORE_FILE,
    demo_paths,
)
from genokeeper.demo.frame import (
    CALIBRATION_NOTE,
    FRAME,
    FRAME_SUMMARY,
    POINTER,
    TABLE_SAMPLES,
)
from genokeeper.explain import LOOSEST
from genokeeper.query import Filters, open_store, query
from genokeeper.store.db import (
    SCHEMA_VERSION,
    Level,
    Validity,
    connect,
    export_snapshots,
)

DOI = "10.5281/zenodo.24134245"
REPOSITORY = "https://github.com/genokeeper/genokeeper"
SOCIAL_IMAGE = "https://genokeeper.github.io/genokeeper/social.png"
PRIVATE_IDS = "genome-private-ids"
# The dashboard's slider order: its `#rel=<i>&ev=<i>` index these.
REL = (Level.LOW, Level.MEDIUM, Level.HIGH)
EV = (Validity.CONFLICTING, Validity.LIMITED, Validity.MODERATE, Validity.STRONG)
_REL_NAMES = ("any genotype", "medium or better", "high only")
_SAY = "genokeeper demo package"


class Leak(Exception):
    """A shipped output holds a needle (each hit names the table and column, or the file,
    and which needle, never its text), or `grep -E` cannot read a private-ID pattern."""


def private_id_file() -> Path | None:
    """The private-ID file the genome-data guard reads, `<git common dir>/info/
    genome-private-ids` of this checkout (`_REPO`); None outside a checkout or without one."""
    try:
        done = subprocess.run(
            [
                "git",
                "-C",
                _REPO,
                "rev-parse",
                "--path-format=absolute",
                "--git-common-dir",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except OSError, subprocess.CalledProcessError:
        return None
    path = Path(done.stdout.strip()) / "info" / PRIVATE_IDS
    return path if path.is_file() else None


def _private_patterns(path: Path) -> dict[int, str]:
    """The file's patterns by line number, as the guard reads them: blank and `#` lines
    skipped."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return {
        n: line
        for n, line in enumerate(lines, 1)
        if line.strip() and not line.lstrip().startswith("#")
    }


def _quote(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _string(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def _columns(con: duckdb.DuckDBPyConnection) -> list[tuple[str, str, str, str]]:
    """Each column of each table of the database, in any schema: (table as named in a
    query, its name, the column, its type)."""
    rows = con.execute(
        "SELECT c.table_schema, c.table_name, c.column_name, c.data_type"
        " FROM information_schema.columns AS c JOIN information_schema.tables AS t"
        " ON t.table_catalog = c.table_catalog AND t.table_schema = c.table_schema"
        " AND t.table_name = c.table_name"
        " WHERE t.table_type = 'BASE TABLE' AND c.table_catalog = current_database()"
        " ORDER BY c.table_schema, c.table_name, c.ordinal_position"
    ).fetchall()
    return [
        (f"{_quote(schema)}.{_quote(table)}", table, column, kind)
        for schema, table, column, kind in rows
    ]


def _longest_first(replacements: Mapping[str, str]) -> list[tuple[str, str]]:
    return sorted(replacements.items(), key=lambda kv: len(kv[0]), reverse=True)


def scrub(con: duckdb.DuckDBPyConnection, replacements: Mapping[str, str]) -> None:
    """Rewrite each key of `replacements` (longest first) in every VARCHAR, JSON and
    VARCHAR[] column of every table of `con`."""
    order = _longest_first(replacements)

    def rewritten(expr: str) -> str:
        for old, new in order:
            expr = f"replace({expr}, {_string(old)}, {_string(new)})"
        return expr

    for table, _, column, kind in _columns(con):
        c = _quote(column)
        value = {
            "VARCHAR": rewritten(c),
            "JSON": f"CAST({rewritten(f'CAST({c} AS VARCHAR)')} AS JSON)",
            "VARCHAR[]": f"list_transform({c}, lambda x: {rewritten('x')})",
        }.get(kind)
        if value is None:
            continue
        held = " OR ".join(
            f"contains(CAST({c} AS VARCHAR), {_string(old)})" for old, _ in order
        )
        con.execute(f"UPDATE {table} SET {c} = {value} WHERE {held}")


def _texts(con: duckdb.DuckDBPyConnection, folder: Path) -> list[Path]:
    """Each column of each table of `con`, its values as text one per line, in
    `folder/<table>.<column>`."""
    folder.mkdir()
    texts: list[Path] = []
    for table, name, column, _ in _columns(con):
        c = _quote(column)
        values = con.execute(
            f"SELECT CAST({c} AS VARCHAR) FROM {table} WHERE {c} IS NOT NULL"
        )
        texts.append(folder / f"{name}.{column}")
        with texts[-1].open("a", encoding="utf-8") as handle:  # a name in two schemas
            while batch := values.fetchmany(10_000):
                handle.writelines(f"{value}\n" for (value,) in batch)
    return texts


def scan_bytes(path: Path, needles: Mapping[str, str]) -> list[str]:
    """`<file name>: <needle's name>` for each string of `needles` in the file's bytes."""
    data = path.read_bytes()
    return [
        f"{path.name}: {name}"
        for name, text in needles.items()
        if text.encode() in data
    ]


def scan_private(paths: Sequence[Path], patterns: Mapping[int, str]) -> list[str]:
    """`<file name>: private ID on line <n>` for each pattern (line number -> ERE) `grep -E`
    finds in a file of `paths`."""
    hits: list[str] = []
    for n, pattern in patterns.items():
        found = subprocess.run(
            ["grep", "-E", "-a", "-l", "-e", pattern, "--", *paths],
            capture_output=True,
            text=True,
            check=False,
        )
        if found.returncode > 1:
            raise Leak(f"grep -E cannot read the private ID on line {n}")
        hits += [
            f"{Path(f).name}: private ID on line {n}" for f in found.stdout.splitlines()
        ]
    return hits


def _note(template: str) -> str:
    """The calibration note as HTML on the site page `template`: on the calibration page
    without its pointer to itself, on the others with a link to that page."""
    head, _, tail = CALIBRATION_NOTE.partition(POINTER)
    if template == "calibration.html":
        return html.escape(head) + "."
    link = '; <a href="calibration.html">the calibration page</a> says'
    return html.escape(head) + link + html.escape(tail)


def _fill(template: str, parts: Mapping[str, str]) -> str:
    """`resources/demo/<template>` with the site's style, its CSP and `parts` filled in, in
    one pass: a filled-in part is never searched for placeholders."""
    folder_ = files("genokeeper").joinpath("resources", "demo")
    css = folder_.joinpath("site.css").read_text(encoding="utf-8")
    frame = _banner((FRAME_SUMMARY, FRAME))
    frame = frame.replace("</details>", f"<p>{_note(template)}</p></details>")
    every = {
        "css": css,
        "csp": f"default-src 'none'; style-src 'sha256-{_sha256(css)}'",
        "frame": f'<div class="frame">{frame}</div>',
        "repository": REPOSITORY,
        **parts,
    }
    page = folder_.joinpath(template).read_text(encoding="utf-8")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: every[m[1]], page)


def _calibration_table() -> str:
    same = {"same sample": ' class="same"', "independent": ""}
    rows = "\n".join(
        f"<tr><td>{html.escape(m.measures)}</td><td>{html.escape(m.samples)}</td>"
        f"<td{same[m.for_hg002]}>{m.for_hg002}</td></tr>"
        for m in TABLE_SAMPLES.values()
    )
    return (
        '<div class="wrap"><table class="cal">\n<thead><tr><th>Class</th>'
        "<th>Measured on</th><th>For HG002</th></tr></thead>\n"
        f"<tbody>\n{rows}\n</tbody>\n</table></div>"
    )


def _day(release: str | None) -> str:
    """A release `YYYYMMDD` as `YYYY-MM-DD`; anything else as it is."""
    text = release or "unknown"
    return (
        f"{text[:4]}-{text[4:6]}-{text[6:]}" if re.fullmatch(r"\d{8}", text) else text
    )


def _count(n: int, noun: str) -> str:
    return f"{n:,} {noun}" + ("" if n == 1 else "s")


def _side(side: Mapping[str, Any] | None) -> str:
    if side is None:
        return "not reported"
    parts = [
        None if side["kind"] == "finding" else str(side["kind"]).replace("_", " "),
        None
        if side["source_level"] is None
        else str(side["source_level"]).replace("_", " ").replace("|", ", "),
        None if side["validity"] is None else f"{side['validity']} evidence",
        None if side["level"] is None else f"genotype {side['level']}",
    ]
    return ", ".join(str(p) for p in parts if p is not None)


def _diff(rows: Sequence[Mapping[str, Any]], before: str, after: str) -> str:
    """The first three changes as a unified diff, one before and one after line each."""
    if not rows:
        return ""
    lines = [f'<div class="h">ClinVar {_day(before)} → {_day(after)}</div>']
    for row in rows[:3]:
        marker = row["marker"]["rsid"] or row["marker"]["grch38"] or row["id"]
        label = f"{row['gene'] or '(no gene)'} {marker}"
        lines += [
            f'<div class="g h">{html.escape(label)} ({html.escape(row["kind"])})</div>',
            f'<div class="m">- {html.escape(_side(row["before"]))}</div>',
            f'<div class="p">+ {html.escape(_side(row["after"]))}</div>',
        ]
    shown = min(len(rows), 3)
    aria = (
        f"{_count(shown, 'finding')} whose grading changed between two ClinVar releases"
    )
    return (
        f'<div class="diff" role="img" aria-label="{aria}">\n'
        + "\n".join(lines)
        + "\n</div>"
    )


def _grid() -> list[list[Filters]]:
    """The matrix's thresholds, rows by `REL` and columns by `EV`."""
    return [[Filters(min_reliability=r, min_validity=v) for v in EV] for r in REL]


def _matrix(counts: Mapping[Filters, int]) -> str:
    """Figure 1's rows: each cell the findings `query` keeps at its pair of thresholds,
    shaded by count, linking to the dashboard at that view; the default pair outlined. "high
    only" is the top row: each link keeps its `REL` index."""
    top = max(counts.values())
    rows: list[str] = []
    for r, row in reversed(list(enumerate(_grid()))):
        cells: list[str] = []
        for v, f in enumerate(row):
            n = counts[f]
            shade = int(5 * log1p(n) / log1p(top)) if top else 0
            label = (
                "<small>default</small>"
                if f == Filters()
                else "<small>everything called</small>"
                if f == LOOSEST
                else ""
            )
            mark = ' class="default"' if f == Filters() else ""
            cells.append(
                f'<td{mark}><a class="t{shade}" href="dashboard.html#rel={r}&amp;ev={v}">'
                f"{n:,}{label}</a></td>"
            )
        rows.append(f"        <tr><th>{_REL_NAMES[r]}</th>{''.join(cells)}</tr>")
    return "\n".join(rows)


def _runs(con: duckdb.DuckDBPyConnection) -> list[dict[str, Any]]:
    """The store's runs, oldest first, each with its ClinVar release."""
    rows = con.execute(
        "SELECT run_id, tools->>'clinvar' FROM run ORDER BY started_at, run_id"
    ).fetchall()
    return [{"run_id": run_id, "clinvar": clinvar} for run_id, clinvar in rows]


def _landing(
    con: duckdb.DuckDBPyConnection,
    runs: Sequence[Mapping[str, Any]],
    reports: Sequence[Mapping[str, Any]],
) -> str:
    """The landing page; `reports` are every change, then those entering or leaving the
    list of everything called and the list at the default thresholds."""
    run = runs[-1]
    grid = [f for row in _grid() for f in row]
    counts = {f: len(query(con, run["run_id"], f).findings) for f in grid}
    loosest, default = counts[LOOSEST], counts[Filters()]
    clinvar = {r["run_id"]: r["clinvar"] for r in runs}
    base = clinvar.get(reports[0]["base"]) if reports[0]["base"] else None
    rows, everything, defaults = (r["changes"] for r in reports)
    changed = (
        f"The reads are the same. Between its releases of {_day(base)} and"
        f" {_day(run['clinvar'])}, ClinVar brought {_count(len(rows), 'change')};"
        f" {_count(len(everything), 'finding')} entered or left the list of everything"
        f" called, and {len(defaults):,} the list at the default thresholds."
    )
    if len(everything) > 3:
        changed += " Three of them:"
    description = (
        f"{loosest:,} candidate findings for GIAB's HG002, {default:,} at the default"
        " thresholds, each graded by a measured error rate and the strength of its evidence."
    )
    version = genokeeper.__version__
    return _fill(
        "index.html",
        {
            "description": html.escape(description),
            "og_title": html.escape(f"One genome, {loosest:,} candidate findings"),
            "og_image": SOCIAL_IMAGE,
            "doi": DOI,
            "loosest": f"{loosest:,}",
            "default": f"{default:,}",
            "matrix": _matrix(counts),
            "caption": html.escape(
                f"Run {run['run_id']}, ClinVar {_day(run['clinvar'])}, genokeeper {version}."
            ),
            "changed": html.escape(changed),
            "diff": _diff(everything, base or "", run["clinvar"] or ""),
            "calibration": _calibration_table(),
            "cite": (
                f"Genokeeper {html.escape(version)}. "
                f'<a href="https://doi.org/{DOI}">doi:{DOI}</a>. HG002 reads: ENA'
                " PRJNA1465120. Truth sets: Genome in a Bottle v4.2.1 and v5.0q. Reference:"
                " GIAB's GRCh38 with the GRC false duplications masked."
            ),
        },
    )


def _rewrite(text: str, replacements: Mapping[str, str]) -> str:
    for old, new in _longest_first(replacements):
        text = text.replace(old, new)
    return text


def _calibration(
    lab: Path,
    con: duckdb.DuckDBPyConnection,
    runs: Sequence[Mapping[str, Any]],
    replacements: Mapping[str, str],
) -> str:
    sections: list[str] = []
    for name in TABLES:
        report = folder(lab, name) / REPORT
        if not report.is_file():
            continue
        text = _rewrite(report.read_text(encoding="utf-8"), replacements)
        title = f"{name}: {TABLE_SAMPLES[name].measures}"
        sections.append(
            f'<details class="report"><summary>{html.escape(title)}</summary>'
            f"<pre>{html.escape(text)}</pre></details>"
        )
    found = con.execute(
        "SELECT tools FROM run WHERE run_id = ?", [runs[-1]["run_id"]]
    ).fetchone()
    tools: dict[str, Any] = json.loads(found[0]) if found is not None else {}
    versions = ", ".join(f"{k} {v}" for k, v in sorted(tools.items()))
    releases = " and ".join(_day(r["clinvar"]) for r in runs)
    methods = [
        f"Tools and versions of run {runs[-1]['run_id']}: {versions}.",
        "Reference: GIAB's GRCh38 with the GRC false duplications masked.",
        f"ClinVar releases: {releases}.",
    ]
    return _fill(
        "calibration.html",
        {
            "note": _note("calibration.html"),
            "calibration": _calibration_table(),
            "reports": "\n".join(sections) or "<p>No table's report was found.</p>",
            "methods": "\n".join(f"    <li>{html.escape(m)}</li>" for m in methods),
        },
    )


def _write_fresh(source: Path, target: Path) -> None:
    """A new database at `target` holding `source`'s schema and rows, none of its stale
    blocks."""
    con = duckdb.connect(str(target))
    try:
        name = con.execute("SELECT current_database()").fetchone()
        assert name is not None
        con.execute(f"ATTACH {_string(str(source))} AS scratch (READ_ONLY)")
        con.execute(f"COPY FROM DATABASE scratch TO {_quote(str(name[0]))}")
        con.execute("DETACH scratch")
        con.execute("CHECKPOINT")
    finally:
        con.close()


def package(lab: Path, out: Path, *, home: Path, private_ids: Path | None) -> Path:
    """Write the demo's files from `demo_paths(lab)` into `out` (see the module); return
    `out`. Raises `Leak`, with nothing under a final name, when an output holds the home
    directory, the lab, the demo root or a private-ID pattern, and `FileExistsError` when
    `out` holds anything but an earlier run's `.tmp-*` folder."""
    root = demo_paths(lab)
    if out.is_dir() and (
        held := sorted(p.name for p in out.iterdir() if not p.name.startswith(".tmp-"))
    ):
        raise FileExistsError(
            f"{out} holds {', '.join(held)}: package writes only into an empty folder"
        )
    needles = {"the home directory": str(home), "the lab": str(lab)}
    needles["the demo root"] = str(root.repo)
    if private_ids is None:
        print(f"{_SAY}: no private-ID file: that check is skipped", file=sys.stderr)
    patterns = {} if private_ids is None else _private_patterns(private_ids)
    replacements = {str(root.repo): "<demo>", str(lab): "<lab>", str(home): "~"}
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".tmp-", dir=out))
    try:
        scratch_root = work / "scratch"
        scratch = Paths(
            repo=scratch_root,
            raw=scratch_root / "raw",
            reports=scratch_root / "reports",
            lab=scratch_root / "lab",
        )
        scratch.reports.mkdir(parents=True)
        store = scratch.reports / STORE_FILE
        shutil.copyfile(root.reports / STORE_FILE, store)
        con = connect(store)
        try:
            scrub(con, replacements)
            con.execute("CHECKPOINT")
            export_snapshots(con, scratch.reports / "snapshots")
        finally:
            con.close()
        _write_fresh(store, work / STORE_FILE)
        log = root.repo / "lab" / "refresh" / refresh_log.LOG
        shutil.copyfile(log, work / REFRESH_LOG_FILE)
        (scratch.lab / "refresh").mkdir(parents=True)
        shutil.copyfile(log, scratch.lab / "refresh" / refresh_log.LOG)
        (scratch.lab / DEMO_MARKER).write_text("{}\n", encoding="utf-8")
        write_dashboard(
            scratch,
            out=work / "dashboard.html",
            lang="en",
            banner=(FRAME_SUMMARY, FRAME, CALIBRATION_NOTE),
        )
        reports = [
            changes.report(
                scratch.reports,
                scratch.lab,
                base=None,
                run=None,
                f=f,
                every=every,
                gene=None,
                marker=None,
                text=None,
            )
            for f, every in ((LOOSEST, True), (LOOSEST, False), (Filters(), False))
        ]
        with open_store(work / STORE_FILE) as shipped:
            runs = _runs(shipped)
            (work / "index.html").write_text(
                _landing(shipped, runs, reports), encoding="utf-8"
            )
            (work / "calibration.html").write_text(
                _calibration(lab, shipped, runs, replacements), encoding="utf-8"
            )
            manifest = {
                "version": genokeeper.__version__,
                "schema_version": SCHEMA_VERSION,
                "runs": runs,
                "date": datetime.now(UTC).date().isoformat(),
                "files": {name: sha256_file(work / name) for name in SHIPPED},
            }
            (work / MANIFEST).write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            files = [work / name for name in (*SHIPPED, MANIFEST)]
            scanned = [*files, *_texts(shipped, work / "text")]
        leaks = [hit for path in scanned for hit in scan_bytes(path, needles)]
        leaks += scan_private(scanned, patterns)
        if leaks:
            raise Leak("refusing to package: " + "; ".join(leaks))
        for path in files:
            os.replace(path, out / path.name)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return out
