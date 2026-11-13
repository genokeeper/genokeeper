"""`python -m genokeeper.demo package OUT`: what the HG002 demo publishes, written from the
demo root (`demo_paths(lab)`) into `OUT`.

The store is copied into a scratch root, every text value of every table (its VARCHAR,
JSON and VARCHAR[] columns) has its absolute paths rewritten (`scrub`: the demo root to
`<demo>`, the lab to `<lab>`, the checkout `package` runs from to `<genokeeper>`), and the shipped store is a fresh
database written from that copy (`COPY FROM DATABASE`): DuckDB keeps a replaced string in the
file's bytes after an UPDATE, and a fresh copy does not. From the scratch copy, with its
snapshots and a lab holding `DEMO_MARKER` and the demo's refresh-log record, come the
dashboard (English, the frame in its banner), the calibration page (each table's report
with the frame's calibration paragraph on top) and the landing page, its numbers from the
store. `MANIFEST` names each file's sha256.

The scan is the control: the adversary is the maintainer's own machine leaking into a public
file by mistake; the anchor is the outputs as written, read back. It covers every file
`package` ships, `MANIFEST` included, and every column of every table of the shipped store
(`_texts`). `scan_bytes` looks in each for the home directory, the lab, the demo root and the
checkout, as exact strings; `scan_paths` refuses any path left in a text output (`PATH`: one
starting at `/` or `~/`), which no rewrite above produced, and `scan_private` runs each pattern of the private-ID file the genome-data
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
import tomllib
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any
from urllib.parse import quote

import duckdb

import genokeeper
from genokeeper import changes, refresh_log
from genokeeper.calibrate.classes import read_class_table
from genokeeper.calibrate.table import CalibrationTable
from genokeeper.calibrate.tables import REPORT, TABLES, folder
from genokeeper.checksum import sha256_file
from genokeeper.config import _REPO, DEMO_MARKER, Paths, load_grading
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
from genokeeper.demo.install import FOLDER, QUESTIONS
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE
from genokeeper.explain import LOOSEST, Explanation, lookup, newest_medlineplus
from genokeeper.pipeline import _newest_dated
from genokeeper.query import RESEARCH_ASSOCIATION, Filters, open_store, query
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
# The header every page shows; the dashboard's links are views in its URL hash.
SITE_NAV = (
    ("Overview", "index.html"),
    ("Results", "dashboard.html#tab=feed"),
    ("Changes", "dashboard.html#tab=changes&rel=0&ev=0"),
    ("How often it is wrong", "calibration.html"),
    ("Sources", "sources.html"),
    ("GitHub", REPOSITORY),
)
# The dashboard's slider order: its `#rel=<i>&ev=<i>` index these.
REL = (Level.LOW, Level.MEDIUM, Level.HIGH)
EV = (Validity.CONFLICTING, Validity.LIMITED, Validity.MODERATE, Validity.STRONG)
_REL_NAMES = ("any genotype", "medium or better", "high only")
_EV_NAMES = ("any evidence", "limited or more", "moderate or more", "strong")
_SHORT = {"medium or better": "medium+", "high only": "high"}
_SHORT |= {"limited or more": "limited+", "moderate or more": "moderate+"}
_SHORT |= {"any genotype": "any", "any evidence": "any", "strong": "strong"}
# The recorded question and answer the landing shows, when the maintainer has recorded one.
TRANSCRIPT: Traversable = files("genokeeper").joinpath(
    "resources", "demo", "transcript.md"
)
_SAY = "genokeeper demo package"
# DeepVariant checkpoints by file name: what to call one, and Google's published figures.
CHECKPOINTS = {
    "weights-51-0.995354.ckpt": (
        "Google's Complete Genomics T7 checkpoint",
        (
            "SNP F1 0.997509 and INDEL F1 0.998207 on HG001 chr20 (DeepVariant r1.10,"
            " Complete Genomics T7 case study)"
        ),
    ),
}
_MOVES = ("left", "entered")  # the change kinds the landing shows one of each
REVERSED = "ClinVar reversed itself"
CLASS_TABLES = ("deletions", "cnv", "insertions")
_COUNTS = ("calls", "false", "found", "truth")
_WRONG_HEAD = (
    "Class",
    "Truth set; measured on",
    "Events in truth",
    "Found",
    "Share of calls false",
    "What the tool does with it",
    "For HG002",
)
# Per group of `citations.toml`: what genokeeper itself measured for it.
GROUP_MEASURED = {
    "Reads to calls": "Measured here: small-variant genotype error rates on HG003 and"
    " HG002, and the run-time checks on all four libraries.",
    "Truth sets and calibration": "Used here as the truth every error rate is measured"
    " against.",
    "Structural variants, copy number and repeats": "Measured here: each class's chance"
    " of being wrong on HG002's two libraries, and for repeats, copy number and paralog"
    " loci on HG003 and HG004 too.",
    "Pharmacogenomics": "No calibration table measures the star-allele calls.",
    "Clinical and population evidence": "These set the strength of evidence, which no"
    " truth set measures.",
    "Variant effect prediction": "Not measured here: a prediction carries no error rate"
    " of its own.",
    "Polygenic scores and ancestry": "Polygenic scores are not measured here; ancestry"
    " is, on HG002, HG003 and HG004.",
    "Haplogroups, HLA and other typing": "Measured here: HLA, chrM variants, runs of"
    " homozygosity and blood groups on HG002, HG003 and HG004, and the Y haplogroup on"
    " HG002's two libraries.",
}


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


# A path from a root of two characters or more, or from `~/`: never after a character that
# makes it part of a word, a URL (`https://…`), a placeholder (`<lab>/…`, `&gt;/…`) or a
# fraction; the demo folder's own default (`FOLDER`), which the landing tells users to open,
# is not one.
PATH = re.compile(rb"(?<![\w.:;/~<>%-])/[\w.@+-]{2,}/[\w.@+-]+|(?<![\w/])~/")


def scan_paths(path: Path) -> list[str]:
    """`<file name>: a path` when the text file holds one (`PATH`)."""
    text = path.read_bytes().replace(FOLDER.encode(), b"")
    return [f"{path.name}: a path"] if PATH.search(text) else []


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


def header() -> str:
    """The site's header: the wordmark, then `SITE_NAV`'s links."""
    links = "".join(
        f' <a href="{html.escape(href)}">{name}</a>' for name, href in SITE_NAV
    )
    mark = '<a class="mark" href="index.html">genokeeper</a>'
    return f'<nav class="site" aria-label="genokeeper">{mark}{links}</nav>'


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
        "header": header(),
        **parts,
    }
    page = folder_.joinpath(template).read_text(encoding="utf-8")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: every[m[1]], page)


def _f1(lab: Path) -> dict[str, float]:
    """HG003's SNP and INDEL F1 from the small-variant table; none without it."""
    try:
        return CalibrationTable.read(folder(lab, "small_variants")).f1
    except OSError, ValueError:
        return {}


def _class_values(path: Path) -> tuple[str, str, str]:
    """A class table's events in truth (a library's, over its bands), each library's found,
    and the share of calls false over every band and library."""
    table = read_class_table(path)
    sums = {
        library: Counter(
            {k: sum(int(bands[b][k]) for b in table.bands) for k in _COUNTS}
        )
        for library, bands in table.measured.items()
    }
    calls, false = (sum(c[k] for c in sums.values()) for k in ("calls", "false"))
    return (
        ", ".join(f"{t:,}" for t in sorted({c["truth"] for c in sums.values()})),
        ", ".join(f"{c['found']:,} ({library})" for library, c in sums.items()),
        f"{false / calls:.0%}" if calls else "no call",
    )


def _values(lab: Path, name: str) -> tuple[str, str, str]:
    """Events in truth, found and share of calls false, as the table `name` holds them."""
    none = ("—", "—", "—")
    if name == "small_variants" and (f1 := _f1(lab)):
        return "—", f"SNP F1 {f1['snp']:.5f}, INDEL F1 {f1['indel']:.5f}", "—"
    if name not in CLASS_TABLES:
        return none
    try:
        return _class_values(folder(lab, name))
    except OSError, ValueError, KeyError:
        return none


def _wrong(lab: Path, tools: Mapping[str, Any]) -> str:
    """ "How often it is wrong": one row per calibration table, from its values, and the
    run's DeepVariant model beside the small-variant F1 and Google's published one."""
    head = "".join(f"<th>{h}</th>" for h in _WRONG_HEAD)
    rows = []
    for name, m in TABLE_SAMPLES.items():
        cells = [m.measures, f"{m.truth}; {m.samples}", *_values(lab, name), m.does]
        same = ' class="same"' if m.for_hg002 == "same sample" else ""
        rows.append(
            "<tr>"
            + "".join(f"<td>{html.escape(c)}</td>" for c in cells)
            + f"<td{same}>{m.for_hg002}</td></tr>"
        )
    model, published = _checkpoint(tools)
    version = f" {v}" if (v := tools.get("deepvariant")) else ""
    said = f"Small variants are called by DeepVariant{version} with {model}."
    if f1 := _f1(lab):
        said += (
            f" Measured here on HG003, chr1-22: SNP F1 {f1['snp']:.5f}, INDEL F1"
            f" {f1['indel']:.5f}."
        )
    if published:
        said += f" Google's published figures for that checkpoint: {published}."
    return (
        f'<div class="wrap"><table class="cal">\n<thead><tr>{head}</tr></thead>\n'
        f"<tbody>\n{'\n'.join(rows)}\n</tbody>\n</table></div>\n"
        f'<p class="sub">{html.escape(said)}</p>'
    )


def _checkpoint(tools: Mapping[str, Any]) -> tuple[str, str | None]:
    """The run's DeepVariant checkpoint by its recorded file name (`CHECKPOINTS`'s, else the
    name) and Google's published figures for it."""
    found = re.search(r"weights-[\w.-]+?\.ckpt", str(tools.get("checkpoint_md5", "")))
    if found is None:
        return "a checkpoint the run did not record", None
    return CHECKPOINTS.get(found[0], (found[0], None))


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
    # A conflict has several causes (ClinVar's, BA1/BS1, a negative curation): name it.
    conflicting = side["validity"] == "conflicting" and side["reason"]
    why = f" ({side['reason']})" if conflicting else ""
    parts = [
        None if side["kind"] == "finding" else str(side["kind"]).replace("_", " "),
        None
        if side["source_level"] is None
        else str(side["source_level"]).replace("_", " ").replace("|", ", "),
        None if side["validity"] is None else f"{side['validity']} evidence{why}",
        None if side["level"] is None else f"genotype {side['level']}",
    ]
    return ", ".join(str(p) for p in parts if p is not None)


def _pathogenic(level: Any) -> bool:
    text = str(level).lower()
    return "pathogenic" in text and "conflicting" not in text


def _carried(side: Mapping[str, Any] | None) -> bool:
    return side is not None and side["dosage"] not in (None, "0", "unknown")


def _examples(
    every: Sequence[Mapping[str, Any]], listed: Sequence[Mapping[str, Any]]
) -> list[tuple[str, Mapping[str, Any]]]:
    """The changes the landing shows, each labelled and the first in its report's order: one
    whose ClinVar classification crossed between pathogenic and not at a variant carried on
    both sides (ClinVar reversed itself), then one that left and one that entered `listed`."""
    reversed_ = (
        r
        for r in every
        if r["id"].startswith("clinvar:")
        and _carried(r["before"])
        and _carried(r["after"])
        and _pathogenic(r["before"]["source_level"])
        != _pathogenic(r["after"]["source_level"])
    )
    picks = [next(reversed_, None)]
    picks += [next((r for r in listed if r["kind"] == k), None) for k in _MOVES]
    found: dict[str, tuple[str, Mapping[str, Any]]] = {}
    for what, r in zip((REVERSED, *_MOVES), picks, strict=True):
        if r is not None:
            found.setdefault(r["id"], (what, r))  # one row once, by its first label
    return list(found.values())


def _diff(
    every: Sequence[Mapping[str, Any]],
    listed: Sequence[Mapping[str, Any]],
    before: str,
    after: str,
) -> str:
    """`_examples` as a unified diff, one before and one after line each."""
    rows = _examples(every, listed)
    if not rows:
        return ""
    days = html.escape(f"ClinVar {_day(before)} → {_day(after)}")
    lines = [f'<div class="h">{days}</div>']
    for what, row in rows:
        marker = row["marker"]["rsid"] or row["marker"]["grch38"] or row["id"]
        label = f"{row['gene'] or '(no gene)'} {marker}"
        lines += [
            f'<div class="g h">{html.escape(label)} ({what})</div>',
            f'<div class="m">- {html.escape(_side(row["before"]))}</div>',
            f'<div class="p">+ {html.escape(_side(row["after"]))}</div>',
        ]
    aria = f"{_count(len(rows), 'change')} between two ClinVar releases"
    return (
        f'<div class="diff" role="img" aria-label="{aria}">\n'
        + "\n".join(lines)
        + "\n</div>"
    )


def _grid() -> list[list[Filters]]:
    """The matrix's thresholds, rows by `REL` and columns by `EV`."""
    return [[Filters(min_reliability=r, min_validity=v) for v in EV] for r in REL]


def _label(name: str) -> str:
    """A matrix header: its words, and the short form a phone shows (`_SHORT`)."""
    return f'<span class="long">{name}</span><span class="short">{_SHORT[name]}</span>'


def _matrix(counts: Mapping[Filters, int]) -> str:
    """The figure's table: each cell the results `query` keeps at its pair of thresholds,
    coloured by the pair (`.m<rel><ev>`: blue down the reliability axis, green across the
    evidence one), named by them and linking to the dashboard at that view; "high only" is
    the top row, and each link keeps its `REL` index."""
    head = "".join(f"<th>{_label(n)}</th>" for n in _EV_NAMES)
    rows = [f'<tr><th class="axis">results that remain</th>{head}</tr>']
    for r, row in reversed(list(enumerate(_grid()))):
        cells: list[str] = []
        for v, f in enumerate(row):
            default = f == Filters()
            small = "default" if default else "everything" if f == LOOSEST else ""
            name = f"{counts[f]:,} results: genotype {_REL_NAMES[r]}, evidence"
            name += f" {_EV_NAMES[v]}" + (", the default" if default else "")
            cells.append(
                f'<td><a class="m{r}{v}{" default" if default else ""}"'
                f' href="dashboard.html#rel={r}&amp;ev={v}" aria-label="{name}">'
                f"{counts[f]:,}{f'<small>{small}</small>' if small else ''}</a></td>"
            )
        rows.append(f"<tr><th>{_label(_REL_NAMES[r])}</th>{''.join(cells)}</tr>")
    return "\n".join(rows)


def _normal(finding: Mapping[str, Any]) -> bool:
    claim, call = finding["claim"] or {}, finding["gene_call"] or {}
    said = f"{claim.get('classification')} {call.get('phenotype')}".lower()
    return "normal" in said or "benign" in said


def _hero(
    findings: Sequence[dict[str, Any]],
    explained: Callable[[dict[str, Any]], Explanation | None],
) -> tuple[dict[str, Any], Explanation] | None:
    """The result the landing opens: carried, not normal, with a shipped explanation; the
    strictest first (genotype reliability, then evidence), then the lower measured chance
    of a wrong genotype (none measured last), then the id."""

    def rank(f: Mapping[str, Any]) -> tuple[int, int, float, str]:
        rel, p = f["reliability"], f["reliability"]["p_wrong"]
        level = EV.index(Validity(f["validity"]["level"]))
        return -REL.index(Level(rel["level"])), -level, 1.0 if p is None else p, f["id"]

    shown = (f for f in findings if _carried(f["genotype"]) and not _normal(f))
    for f in sorted(shown, key=rank):
        if (said := explained(f)) is not None:
            return f, said
    return None


def _card(finding: Mapping[str, Any], said: Explanation) -> str:
    """One result opened as the dashboard opens it: both bars, the explanation (its gist
    first, when it has one) and its model, and a link to its row."""
    rel, ev, e = finding["reliability"], finding["validity"], html.escape
    r, v = REL.index(Level(rel["level"])), EV.index(Validity(ev["level"]))
    copies = {"1": "one copy", "2": "two copies"}.get(finding["genotype"]["dosage"])
    wrong = "no measured error rate"
    if (p := rel["p_wrong"]) is not None:
        lo, hi = rel["p_wrong_ci"] or (None, None)
        wrong = f"wrong with a measured chance of {p:.2g}" + (
            "" if lo is None else f" [{lo:.2g}, {hi:.2g}]"
        )
    source = (
        finding["sources"][0] if finding["sources"] else {"source": "", "version": ""}
    )
    level = str(ev["source_level"] or "").replace("_", " ")
    marker = finding["marker"]["rsid"] or finding["marker"]["grch38"] or ""
    href = f"dashboard.html#rel={r}&amp;ev={v}&amp;q={quote(finding['id'], safe='')}"
    gist = f'<p class="gist">{e(said.gist)}</p>' if said.gist else ""
    return f"""    <aside class="card" aria-label="One result, opened">
      <p class="sub">{e(finding["topic"])} · {e(finding["gene"] or "")} {e(marker)}{f" · {copies} carried" if copies else ""}</p>
      <h3>{e(said.title)}</h3>
      <div class="pair"><span>genotype</span><span class="bar"><i class="r{r}"></i></span><b>{e(rel["level"])}</b></div>
      <p class="sub">{e(wrong)}</p>
      <div class="pair"><span>evidence</span><span class="bar"><i class="e{v}"></i></span><b>{e(ev["level"])}</b></div>
      <p class="sub">{e(f"{source['source']} {source['version']}: {level}")}</p>
      {gist}<p>{e(said.text)}</p>
      <p class="sub">Generated by {e(said.model)}: an interpretation of the sources ({e(", ".join(said.cites))}), not one of them.</p>
      <p class="sub"><a href="{href}">Open its row</a></p>
    </aside>"""


def _transcript() -> str:
    """The recorded question and answer (`TRANSCRIPT`, its front matter naming when and with
    which version), escaped; none without the file."""
    if not TRANSCRIPT.is_file():
        return ""
    _, head, body = TRANSCRIPT.read_text(encoding="utf-8").split("---\n", 2)
    meta = dict(line.split(": ", 1) for line in head.splitlines() if ": " in line)
    said = f"Recorded on {meta['recorded']} from genokeeper {meta['version']}."
    return (
        f'  <p class="sub">{html.escape(said)}</p>\n'
        f'  <pre class="transcript">{html.escape(body.strip())}</pre>'
    )


def _withheld(con: duckdb.DuckDBPyConnection, run_id: str) -> int:
    """The values a calibration table withheld in the run: each a claim's `WITHHELD` entry."""
    found = con.execute(
        "SELECT sum(json_array_length(claim, '$.withheld')) FROM statement"
        " WHERE run_id = ? AND json_exists(claim, '$.withheld')",
        [run_id],
    ).fetchone()
    return int(found[0] or 0) if found is not None else 0


def citations() -> list[dict[str, str]]:
    """`resources/demo/citations.toml`'s entries, in its order: the pipeline's."""
    text = files("genokeeper").joinpath("resources", "demo", "citations.toml")
    return tomllib.loads(text.read_text(encoding="utf-8"))["tool"]


def _built_on(tools: Mapping[str, Any]) -> str:
    """Each group's measured line and entries, each with the run's version of its
    `version_key`, where it names one."""
    groups: dict[str, list[str]] = {}
    for t in citations():
        version = tools.get(t["version_key"]) if t["version_key"] else None
        if version is not None and t["version_key"] == "deepvariant":
            version = f"{version} with {_checkpoint(tools)[0]}"
        groups.setdefault(t["group"], []).append(
            f'<li><a href="{html.escape(t["url"])}">{html.escape(t["name"])}</a>'
            + ("" if version is None else f" {html.escape(str(version))}")
            + f". {html.escape(t['citation'])}.</li>"
        )
    return "\n".join(
        f'<h3>{html.escape(g)}</h3>\n<p class="sub">{html.escape(GROUP_MEASURED[g])}'
        f' <a href="calibration.html">Calibration</a></p>\n'
        f'<ul class="plain methods">\n{"\n".join(items)}\n</ul>'
        for g, items in groups.items()
    )


def _tools(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, Any]:
    """The tool versions the run recorded."""
    found = con.execute("SELECT tools FROM run WHERE run_id = ?", [run_id]).fetchone()
    return json.loads(found[0]) if found is not None else {}


def _runs(con: duckdb.DuckDBPyConnection) -> list[dict[str, Any]]:
    """The store's runs, oldest first, each with its ClinVar release."""
    rows = con.execute(
        "SELECT run_id, tools->>'clinvar' FROM run ORDER BY started_at, run_id"
    ).fetchall()
    return [{"run_id": run_id, "clinvar": clinvar} for run_id, clinvar in rows]


def _landing(
    lab: Path,
    con: duckdb.DuckDBPyConnection,
    runs: Sequence[Mapping[str, Any]],
    reports: Sequence[Mapping[str, Any]],
    explained: Callable[[dict[str, Any]], Explanation | None],
) -> str:
    """The landing page; `reports` are every change, then those entering or leaving the
    list of everything called and the list at the default thresholds; `explained` gives a
    result's shipped explanation."""
    run, g = runs[-1], load_grading()
    grid = [f for row in _grid() for f in row]
    results = {f: query(con, run["run_id"], f) for f in grid}
    kept = {f: found.findings for f, found in results.items()}
    counts = {f: len(found) for f, found in kept.items()}
    research = results[LOOSEST].hidden.get(RESEARCH_ASSOCIATION, 0)
    hero = _hero(kept[LOOSEST], explained)
    tools = _tools(con, run["run_id"])
    loosest, default = counts[LOOSEST], counts[Filters()]
    said = sum(explained(f) is not None for f in kept[LOOSEST])
    print(
        f"{_SAY}: {said:,} of {loosest:,} results carry a shipped explanation",
        file=sys.stderr,
    )
    if hero is None:
        print(
            f"{_SAY}: warning: no carried, non-normal result has a shipped explanation;"
            " the landing's opened-result card is empty",
            file=sys.stderr,
        )
    topics = Counter(str(found["topic"]) for found in kept[Filters()])
    by_topic = ", ".join(
        f"{t} {n:,}" for t, n in sorted(topics.items(), key=lambda kv: (-kv[1], kv[0]))
    )
    clinvar = {r["run_id"]: r["clinvar"] for r in runs}
    base = clinvar.get(reports[0]["base"]) if reports[0]["base"] else None
    rows, everything, defaults = (r["changes"] for r in reports)
    changed = (
        f"The reads are the same. Between its releases of {_day(base)} and"
        f" {_day(run['clinvar'])}, ClinVar brought {_count(len(rows), 'change')};"
        f" {_count(len(everything), 'result')} entered or left the list of everything"
        f" called, and {len(defaults):,} entered or left the list at the default"
        " thresholds."
    )
    description = (
        f"{loosest:,} results for GIAB's HG002, {default:,} at the default"
        " thresholds, each graded by a measured error rate and the strength of its evidence."
    )
    version = genokeeper.__version__
    return _fill(
        "index.html",
        {
            "description": html.escape(description),
            "og_title": html.escape(f"One genome, {loosest:,} results"),
            "og_image": SOCIAL_IMAGE,
            "doi": DOI,
            "question": html.escape(QUESTIONS[0]),
            "loosest": f"{loosest:,}",
            "default": f"{default:,}",
            "matrix": _matrix(counts),
            "card": "" if hero is None else _card(*hero),
            "counted": html.escape(
                "A result is any graded statement, a normal one included; a finding is a"
                " result that is not normal. A research layer outside these counts holds"
                f" {_count(research, 'GWAS association')}."
            ),
            "transcript": _transcript(),
            "caption": html.escape(
                f"The {default:,} at the default thresholds by topic: {by_topic}."
                f" Run {run['run_id']}, ClinVar {_day(run['clinvar'])}, genokeeper {version}."
            ),
            "changed": html.escape(changed),
            "diff": _diff(rows, everything, base or "", run["clinvar"] or ""),
            "calibration": _wrong(lab, tools),
            "levels": html.escape(
                "Genotype reliability is high when the measured chance the genotype is"
                f" wrong is at most {g.p_high:g} and the site is in no difficult region,"
                f" low above {g.p_low:g} or at a no-call, and medium between. For a ClinVar"
                " claim, evidence is conflicting when its submitters conflict, gnomAD's"
                f" frequency is above BA1 ({g.ba1:g}) or the gene's BS1, or a ClinGen"
                " curation is disputed or refuted; strong at 3 or 4 review stars; below"
                " that, limited without a positive curation of the gene, else strong at 2"
                " stars, moderate at 1 and limited at none; a drug-response row quotes CPIC"
                " or ClinPGx, never the stars. Any evidence counts the conflicting claims"
                " too."
            ),
            "built": html.escape(
                f"Genokeeper runs {len(citations())} published tools, truth sets and"
                " knowledge sources."
            ),
            "withheld": html.escape(
                f"This run holds {_count(_withheld(con, run['run_id']), 'withheld value')}:"
                " a calibration table could not measure the item, so it is not graded and"
                " keeps the table's reason."
            ),
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
    tools = _tools(con, runs[-1]["run_id"])
    versions = ", ".join(f"{k} {v}" for k, v in sorted(tools.items()))
    releases = " and ".join(_day(r["clinvar"]) for r in runs)
    methods = [
        f"Run {runs[-1]['run_id']}, genokeeper {genokeeper.__version__}.",
        "Reference: GIAB's GRCh38 with the GRC false duplications masked.",
        f"ClinVar releases: {releases}.",
    ]
    return _fill(
        "calibration.html",
        {
            "note": _note("calibration.html"),
            "calibration": _wrong(lab, tools),
            "reports": "\n".join(sections) or "<p>No table's report was found.</p>",
            "methods": "\n".join(f"    <li>{html.escape(m)}</li>" for m in methods),
            "versions": html.escape(f"Tools and versions: {versions}."),
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


def package(
    lab: Path,
    out: Path,
    *,
    home: Path,
    private_ids: Path | None,
    no_private_ids: bool = False,
) -> Path:
    """Write the demo's files from `demo_paths(lab)` into `out` (see the module); return
    `out`. Raises `Leak`, with nothing under a final name, when an output holds the home
    directory, the lab, the demo root or a private-ID pattern, and `FileExistsError` when
    `out` holds anything but an earlier run's `.tmp-*` folder; `Leak` too without a
    pattern in `private_ids`, unless `no_private_ids` waives that check (and the run says so)."""
    root = demo_paths(lab)
    if out.is_dir() and (
        held := sorted(p.name for p in out.iterdir() if not p.name.startswith(".tmp-"))
    ):
        raise FileExistsError(
            f"{out} holds {', '.join(held)}: package writes only into an empty folder"
        )
    needles = {"the home directory": str(home), "the lab": str(lab)}
    needles |= {"the demo root": str(root.repo), "the checkout": str(_REPO)}
    patterns = {} if private_ids is None else _private_patterns(private_ids)
    if not patterns and not no_private_ids:
        raise Leak(
            "no private-ID file, or one with no pattern (genome-private-ids in the"
            " checkout's git folder): write it, or pass --no-private-ids to package"
            " without that check"
        )
    if not patterns:
        print(
            f"{_SAY}: --no-private-ids: the private-ID check is skipped",
            file=sys.stderr,
        )
    replacements = {
        str(root.repo): "<demo>",
        str(lab): "<lab>",
        str(_REPO): "<genokeeper>",
    }
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
        # The demo root's explanations, keyed on prompts grounded in the lab's MedlinePlus.
        if (explanations := root.reports / "explanations").is_dir():
            shutil.copytree(explanations, scratch.reports / explanations.name)
        if summaries := _newest_dated(lab / "medlineplus", [MEDLINEPLUS_FILE]):
            shutil.copytree(summaries, scratch.lab / "medlineplus" / summaries.name)
        write_dashboard(
            scratch,
            out=work / "dashboard.html",
            lang="en",
            banner=(FRAME_SUMMARY, FRAME, CALIBRATION_NOTE),
        )
        # The site's header in the dashboard's banner slot, which escapes what it is given.
        page = (work / "dashboard.html").read_text(encoding="utf-8")
        slot = '<main class="frame"><details class="frame-note">'
        assert page.count(slot) == 1, "the dashboard's banner slot moved"
        page = page.replace(slot, slot.replace("<details", f"{header()}<details"))
        (work / "dashboard.html").write_text(page, encoding="utf-8")
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
        medlineplus = newest_medlineplus(scratch.lab)

        def explained(finding: dict[str, Any]) -> Explanation | None:
            return lookup(finding, medlineplus, scratch.reports / "explanations")

        with open_store(work / STORE_FILE) as shipped:
            runs = _runs(shipped)
            (work / "index.html").write_text(
                _landing(lab, shipped, runs, reports, explained), encoding="utf-8"
            )
            (work / "calibration.html").write_text(
                _calibration(lab, shipped, runs, replacements), encoding="utf-8"
            )
            tools = _tools(shipped, runs[-1]["run_id"])
            (work / "sources.html").write_text(
                _fill("sources.html", {"built": _built_on(tools)}), encoding="utf-8"
            )
            manifest = {
                "version": genokeeper.__version__,
                "schema_version": SCHEMA_VERSION,
                "runs": runs,
                "files": {name: sha256_file(work / name) for name in SHIPPED},
            }
            (work / MANIFEST).write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            files = [work / name for name in (*SHIPPED, MANIFEST)]
            scanned = [*files, *_texts(shipped, work / "text")]
        leaks = [hit for path in scanned for hit in scan_bytes(path, needles)]
        leaks += [hit for p in scanned if p.name != STORE_FILE for hit in scan_paths(p)]
        leaks += scan_private(scanned, patterns)
        if leaks:
            raise Leak("refusing to package: " + "; ".join(leaks))
        for path in files:
            os.replace(path, out / path.name)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return out
