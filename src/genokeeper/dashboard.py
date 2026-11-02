"""The dashboard (spec item 10): one static page over a run's findings.

The page's data comes from the query alone. `payload` runs `query.query` once per slider and
topic combination and keeps each result as a view of finding ids beside the query's own
hidden and unassessed counts, so the page's script only looks a view up and never filters.

The Changes tab's data (`changes`) comes from the run snapshots (`genokeeper.history`): the
comparison with each of the last snapshots before the run, each change's kind at every
slider pair, so the tab's line between trusted and not is always the feed's. Each snapshot
is queried once at the loosest filters, and a slider pair's trusted set is the findings
of that result the pair's filters keep, as `query.query` would keep them.

`render` fills `resources/dashboard/page.html`. Its Content-Security-Policy allows nothing
but the inlined script and style, by the sha256 of their exact bytes, and the data sits in a
`<script type="application/json">` the browser never runs, with every `<` escaped, so no
source text can close it or add markup; the script builds nodes with `textContent` only.
Its one URL is PubMed's prefix, an anchor's `href` before a PMID of digits only; a
REViewer image is linked as the local file it is.

The page is in one language, `genokeeper.lang.language()`'s: every text it shows is a message
of `resources/dashboard/messages.json` (`messages`), embedded in the data in that language
and, for the page's own headings, filled into `page.html`. An English message stored as null
is the analyser's own English text (`_english`), kept once.

The carrier topic lists what no run looks at (`carrier.NOT_COVERED`) and each gene whose main
mechanism is not assessed (the packaged mechanisms file), in the page's language, less a gene
the run's own stage-3 statement assessed.

Under the feed, the "Excluded by screening" block lists the run's flagged non-findings
(`_excluded`), each with the test that would confirm it, and the stage-3 block the run's
other stage-3 watch items and non-findings (`_screened`), which the feed and its counts never
hold, a watch item apart from a non-finding.
"""

import argparse
import base64
import hashlib
import html
import json
import re
import sys
from dataclasses import asdict
from importlib.resources import files
from importlib.resources.abc import Traversable
from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper import changes, check, history, refresh_log
from genokeeper.carrier import (
    CFTR_REPEAT_ID,
    NOT_COVERED,
    STAGE3_SOURCES,
    load_mechanisms,
)
from genokeeper.config import Paths, load_paths
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE, MedlinePlus
from genokeeper.explain import LOOSEST, lookup, newest_medlineplus
from genokeeper.history import Change, ChangeKind, Dates, Snapshot
from genokeeper.lang import LANGUAGES, Language, language
from genokeeper.query import (
    HAPLOID_DOSAGE_NOTE,
    OWNS_GENE_CALL,
    Filters,
    RunNotFound,
    _hidden_by,
    open_store,
    query,
    text_fields,
)
from genokeeper.schedule import STALE_DAYS
from genokeeper.sources.research import is_research_layer
from genokeeper.store.db import Kind, Level, Topic, Validity

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

# The view the page opens on: the query's default filters, every topic.
_DEFAULT = Filters()
# The snapshots the Changes tab compares the run with: refreshes first, then plain runs.
BASES = 12
_PAIRS = [Filters(min_reliability=r, min_validity=v) for r in Level for v in Validity]
# The claim's and gene call's result fields a row shows, which its search string holds too.
_RESULTS = ("haplogroup", "class", "consequence", "hgvsc", "hgvsp")
# The messages `page.html` holds, filled in by `render`.
_PAGE = ("page_title", "tab_feed", "tab_changes", "search_label", "search_placeholder")


def _folder() -> Traversable:
    return files("genokeeper").joinpath("resources", "dashboard")


def _english() -> dict[str, Any]:
    """The catalog's English messages that are the analyser's own texts."""
    return {
        "haploid_note": HAPLOID_DOSAGE_NOTE,
        "carrier_not_covered": list(NOT_COVERED),
        "carrier_mechanisms": load_mechanisms(),
    }


def messages(lang: Language) -> dict[str, Any]:
    """Each message of `resources/dashboard/messages.json` in `lang`: a text with `{name}`
    placeholders, its forms by plural category, or a table; an English null is `_english`'s."""
    text = _folder().joinpath("messages.json").read_text(encoding="utf-8")
    catalog: dict[str, dict[str, Any]] = json.loads(text)
    own = _english() if lang == "en" else {}
    return {
        key: own[key] if entry[lang] is None else entry[lang]
        for key, entry in catalog.items()
    }


def _key(f: Filters) -> str:
    """A view's key: `<reliability>|<validity>|<topic>`, the topic empty for all of them."""
    return f"{f.min_reliability}|{f.min_validity}|{f.topic or ''}"


def _pair(f: Filters) -> str:
    """A Changes view's key: `<reliability>|<validity>`."""
    return f"{f.min_reliability}|{f.min_validity}"


def _day(snapshot: Snapshot) -> str:
    return history.local_day(snapshot).isoformat()


def _trusted(root: Path, run_id: str) -> dict[str, frozenset[str]]:
    """The snapshot's trusted ids at each slider pair, from one query at the loosest
    filters: a finding is trusted at a pair when the pair's filters keep it."""
    found = history.findings(root, run_id, LOOSEST)
    return {
        _pair(f): frozenset(x["id"] for x in found if _hidden_by(x, f) is None)
        for f in _PAIRS
    }


def _kinds(
    change: Change,
    before: dict[str, frozenset[str]],
    after: dict[str, frozenset[str]],
) -> dict[str, ChangeKind]:
    """The kind of `change` at each slider pair, given each snapshot's trusted ids there."""
    return {
        _pair(f): history.kind_of(
            change,
            f,
            change.before is not None
            and change.before.statement_id in before[_pair(f)],
            change.after is not None and change.statement_id in after[_pair(f)],
        )
        for f in _PAIRS
    }


def _timeline(
    root: Path, snaps: list[Snapshot], ids: dict[str, set[str]]
) -> dict[str, list[dict[str, Any]]]:
    """For each statement, its kind, validity and level in each snapshot of `snaps`,
    found under any of its `ids` (a ClinVar ID merge); absent as nulls."""
    wanted = set().union(*ids.values()) if ids else set()
    timeline: dict[str, list[dict[str, Any]]] = {sid: [] for sid in ids}
    for snap in snaps:
        held = history.states(root, snap.run_id, wanted)
        for sid, names in ids.items():
            state = next((held[n] for n in sorted(names) if n in held), None)
            timeline[sid].append(
                {
                    "run": snap.run_id,
                    "date": _day(snap),
                    "kind": None if state is None else state.kind,
                    "validity": None if state is None else state.validity,
                    "level": None if state is None else state.level,
                }
            )
    return timeline


def _changes(
    reports: Path | None,
    lab: Path | None,
    run_id: str,
    dates: dict[str, Dates],
) -> dict[str, Any]:
    """The Changes tab's data for `run_id`: the sources its refresh could not fetch, the
    daily update check's `ran_at`, `newer` and `failed` (no age), and the comparison with each of the last `BASES` snapshots before it (refreshes first)."""
    block: dict[str, Any] = {
        "run": run_id,
        "date": None,
        "oldest": None,
        "default_base": None,
        "note": changes.NO_BASE,
        "first_tier": [str(k) for k in ChangeKind if k in changes.FIRST_TIER],
        "failed": [],
        "check": {**check.summary(None), "stale_days": STALE_DAYS},
        "bases": {},
    }
    if lab is not None:
        # The age is the viewer's to compute from `ran_at`: the page outlives the check.
        block["check"] = {
            **check.summary(check.read_state(lab)),
            "stale_days": STALE_DAYS,
        }
    if reports is None:
        return block
    root = reports / "snapshots"
    records = refresh_log.read(lab) if lab is not None else []
    if lab is not None:
        block["failed"] = changes.failed(changes.refresh_of(lab, run_id))
    snaps = history.snapshots(root)
    order = {s.run_id: i for i, s in enumerate(snaps)}
    if run_id not in order:
        return block
    block["date"], block["oldest"] = _day(snaps[order[run_id]]), _day(snaps[0])
    older = snaps[: order[run_id]][::-1]
    refreshed = {r.run for r in records if r.outcome in ("complete", "partial")}
    bases = [s for s in older if s.run_id in refreshed]
    bases = (bases + [s for s in older if s.run_id not in refreshed])[:BASES]
    default = refresh_log.default_base(lab, root, run_id) if lab is not None else None
    if default == run_id:  # a refresh that made no run: only its publications are new
        bases.insert(0, snaps[order[run_id]])
    if not bases:
        return block
    block["note"] = None
    block["default_base"] = default or bases[0].run_id
    trusted = {
        s.run_id: _trusted(root, s.run_id) for s in [snaps[order[run_id]], *bases]
    }
    ids: dict[str, set[str]] = {}
    for base in bases:
        found, publications = changes.compared(reports, lab, base.run_id, run_id)
        rows: list[dict[str, Any]] = []
        views: dict[str, dict[str, str]] = {_pair(f): {} for f in _PAIRS}
        research: dict[tuple[ChangeKind, tuple[str, ...]], int] = {}
        for sid, change in found.items():
            state = change.after or change.before
            assert state is not None
            kinds = _kinds(change, trusted[base.run_id], trusted[run_id])
            default = kinds[_pair(_DEFAULT)]
            if is_research_layer(state.source):
                counted = (default, change.causes)
                research[counted] = research.get(counted, 0) + 1
                continue
            for key, kind in kinds.items():
                views[key][sid] = str(kind)
            listed = changes.row(sid, change, default, dates, publications)
            rows.append({**listed, "source": state.source})
            names = {sid} | ({change.before.statement_id} if change.before else set())
            ids.setdefault(sid, set()).update(names)
        rows.sort(key=changes.order)
        block["bases"][base.run_id] = {
            "date": _day(base),
            "refresh": base.run_id in refreshed,
            "rows": {r["id"]: r for r in rows},
            "views": views,
            "research": [
                {"kind": str(k), "causes": list(causes), "count": n}
                for (k, causes), n in sorted(
                    research.items(),
                    key=lambda kv: (list(ChangeKind).index(kv[0][0]), kv[0][1]),
                )
            ],
        }
    timeline = _timeline(root, snaps[: order[run_id] + 1], ids)
    for compared in block["bases"].values():
        for sid, listed in compared["rows"].items():
            listed["history"] = timeline[sid]
    return block


def _not_covered(
    con: duckdb.DuckDBPyConnection,
    run_id: str,
    listed: list[str],
    mechanisms: dict[str, str],
) -> dict[str, Any]:
    """The carrier topic's not-covered lists, `listed` and the genes of `mechanisms`: a
    gene whose stage-3 statement this run assessed is covered, as `carrier.genotype` reads
    it."""
    marks = ", ".join("?" * len(STAGE3_SOURCES))
    assessed = {
        str(gene).upper()
        for (gene,) in con.execute(
            "SELECT gene FROM statement WHERE run_id = ? AND vid IS NULL"
            f" AND source IN ({marks}) AND kind <> ?",
            [run_id, *STAGE3_SOURCES, str(Kind.UNASSESSED)],
        ).fetchall()
    }
    genes = {g: text for g, text in mechanisms.items() if g not in assessed}
    return {"all": listed, "genes": genes}


def _excluded(con: duckdb.DuckDBPyConnection, run_id: str) -> list[dict[str, str]]:
    """The run's exclusions (the "Excluded by screening" block): each non-finding flagged
    `exclusion`, by gene then id, with its stored reason and its claim's confirming test,
    which the writer refused to store a flagged statement without."""
    rows = con.execute(
        "SELECT statement_id, gene, reason, claim->>'confirm' FROM statement"
        " WHERE run_id = ? AND kind = ? AND exclusion ORDER BY gene, statement_id",
        [run_id, str(Kind.NON_FINDING)],
    ).fetchall()
    return [
        {"id": sid, "gene": gene, "wording": reason, "confirm": confirm}
        for sid, gene, reason, confirm in rows
    ]


def _screened(
    con: duckdb.DuckDBPyConnection, run_id: str
) -> dict[str, list[dict[str, str | None]]]:
    """The run's stage-3 gene-level watch items (`watch`) and non-findings (`normal`) not
    flagged `exclusion`, each by gene then id: its gene with its gene call's diplotype
    (`label`) and its stored reason (`wording`, or null). Not a gene call's phenotype, which
    can say more than `carrier.genotype` answers (SMN1's "not a carrier"), nor CFTR's
    poly-T/TG, which `carrier` keeps apart: its "no 5T allele" is no cystic-fibrosis test."""
    marks = ", ".join("?" * len(STAGE3_SOURCES))
    rows = con.execute(
        "SELECT s.statement_id, s.kind, s.gene, g.diplotype, s.reason"
        " FROM statement AS s LEFT JOIN gene_call AS g"
        " ON g.run_id = s.run_id AND g.gene = s.gene"
        f" AND {OWNS_GENE_CALL.format(s='s')}"
        " WHERE s.run_id = ? AND s.vid IS NULL AND s.kind IN (?, ?) AND NOT s.exclusion"
        f" AND s.source IN ({marks}) AND s.statement_id <> ?"
        " ORDER BY s.gene, s.statement_id",
        [
            run_id,
            str(Kind.WATCH),
            str(Kind.NON_FINDING),
            *STAGE3_SOURCES,
            CFTR_REPEAT_ID,
        ],
    ).fetchall()
    found: dict[str, list[dict[str, str | None]]] = {"watch": [], "normal": []}
    for sid, kind, gene, diplotype, wording in rows:
        # SMN's and a paralog's diplotype names its gene already ("SMN1 2, SMN2 2").
        label = (
            gene
            if diplotype is None
            else diplotype
            if diplotype.startswith(gene)
            else f"{gene} {diplotype}"
        )
        listed = {"id": sid, "label": label, "wording": wording}
        found["watch" if kind == Kind.WATCH else "normal"].append(listed)
    return found


def payload(
    con: duckdb.DuckDBPyConnection,
    run_id: str | None,
    cache_dir: Path,
    medlineplus: MedlinePlus | None = None,
    *,
    reports: Path | None = None,
    lab: Path | None = None,
    lang: Language = "en",
) -> dict[str, Any]:
    """The page's data in `lang` for `run_id` (None: the latest run): its messages
    (`messages` less the carrier lists), every finding once, by id (`items`), each with its
    explanation from `cache_dir` or null, its search string and its dates, and one view per
    reliability, validity and topic (or none) combination;
    and `changes`, the Changes tab's, from `reports/snapshots` and the refresh log in
    `lab` (without `reports`, no comparison); `excluded`, the run's exclusions; `screened`,
    its other stage-3 watch items and non-findings."""
    dates = history.dates(reports / "snapshots", run_id) if reports is not None else {}
    loosest = query(con, run_id, LOOSEST, dates=dates if reports is not None else None)

    def _item(finding: dict[str, Any]) -> dict[str, Any]:
        found = lookup(finding, medlineplus, cache_dir, lang=lang)
        claim, call = finding["claim"] or {}, finding["gene_call"] or {}
        results = [claim.get(k) for k in _RESULTS] + [
            call.get("diplotype"),
            call.get("phenotype"),
        ]
        fields = text_fields(finding) + [str(r) for r in results if r is not None]
        if found is not None:
            fields += [found.title, found.text]
        # One field per line, so that a search never matches across two of them; the
        # page lower-cases the typed text and reads `_` as a space the same way.
        search = "\n".join(fields).lower().replace("_", " ")
        explanation = None if found is None else asdict(found)
        return {**finding, "explanation": explanation, "search": search}

    views: dict[str, Any] = {}
    for reliability, validity, topic in product(Level, Validity, (None, *Topic)):
        f = Filters(topic=topic, min_reliability=reliability, min_validity=validity)
        result = query(con, loosest.run_id, f)
        views[_key(f)] = {
            "findings": [finding["id"] for finding in result.findings],
            "hidden": result.hidden,
            "unassessed": result.unassessed,
        }
    text = messages(lang)
    listed, mechanisms = text.pop("carrier_not_covered"), text.pop("carrier_mechanisms")
    return {
        "lang": lang,
        "messages": text,
        "run_id": loosest.run_id,
        "versions": loosest.versions,
        "items": {finding["id"]: _item(finding) for finding in loosest.findings},
        "views": views,
        "default": _key(_DEFAULT),
        "changes": _changes(reports, lab, loosest.run_id, dates),
        "not_covered": _not_covered(con, loosest.run_id, listed, mechanisms),
        "excluded": _excluded(con, loosest.run_id),
        "screened": _screened(con, loosest.run_id),
    }


def _sha256(text: str) -> str:
    return base64.b64encode(hashlib.sha256(text.encode()).digest()).decode()


def render(payload: dict[str, Any]) -> str:
    """The page for `payload`: `page.html` with its language, CSP, style, script, data and
    headings filled in."""
    folder = _folder()
    css = folder.joinpath("dashboard.css").read_text(encoding="utf-8")
    js = folder.joinpath("dashboard.js").read_text(encoding="utf-8")
    parts = {
        "csp": (
            f"default-src 'none'; script-src 'sha256-{_sha256(js)}';"
            f" style-src 'sha256-{_sha256(css)}'"
        ),
        "css": css,
        "js": js,
        "data": json.dumps(payload, ensure_ascii=False).replace("<", "\\u003c"),
        "lang": html.escape(payload["lang"]),
        **{key: html.escape(payload["messages"][key]) for key in _PAGE},
    }
    # One pass over the template: a filled-in part is never searched for placeholders.
    page = folder.joinpath("page.html").read_text(encoding="utf-8")
    return re.sub(r"\{\{(\w+)\}\}", lambda m: parts[m[1]], page)


def _write(
    db: Path,
    run_id: str | None,
    cache: Path,
    out: Path,
    reports: Path,
    lab: Path,
    lang: Language,
) -> Path:
    """Write the page in `lang` for `run_id` of the store `db` to `out`, its history from
    `reports/snapshots`; a missing MedlinePlus file is one warning on stderr. Raises
    `RunNotFound`."""
    try:
        medlineplus = newest_medlineplus(lab)
        missing = f"no {lab}/medlineplus/<YYYYMMDD>/{MEDLINEPLUS_FILE}"
        warning = missing if medlineplus is None else None
    except ValueError as e:
        medlineplus, warning = None, str(e)
    with open_store(db) as con:
        data = payload(
            con, run_id, cache, medlineplus, reports=reports, lab=lab, lang=lang
        )
    page = render(data)
    # After the store, so that a refusal stays one line.
    if warning is not None:
        print(
            f"genokeeper dashboard: {warning}: explanations grounded in MedlinePlus are"
            " not shown; run scripts/bootstrap.sh",
            file=sys.stderr,
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    tmp.write_text(page, encoding="utf-8")
    tmp.replace(out)
    return out


def write_dashboard(paths: Paths, run_id: str | None = None) -> Path:
    """Write `reports/dashboard.html` for `run_id` (None: the latest run) of
    `reports/genomelab.duckdb`, with both tabs, in `language()`'s language; return its
    path. Raises `RunNotFound`."""
    reports = paths.reports
    return _write(
        reports / "genomelab.duckdb",
        run_id,
        reports / "explanations",
        reports / "dashboard.html",
        reports,
        paths.lab,
        language(),
    )


def _dashboard_command(args: argparse.Namespace) -> int:
    try:
        lang = language(args.lang)
    except ValueError as e:  # an unknown $GENOKEEPER_LANG
        print(f"genokeeper dashboard: {e}", file=sys.stderr)
        return 2
    try:
        # The snapshots and refresh files sit beside the store, as `query` reads them.
        reports, lab = args.db.parent, load_paths().lab
        out = _write(args.db, args.run, args.cache, args.out, reports, lab, lang)
    except RunNotFound as e:
        print(f"genokeeper dashboard: {e}", file=sys.stderr)
        return 2
    print(f"genokeeper dashboard: wrote {out}")
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `dashboard` sub-command."""
    reports = load_paths().reports
    d = subparsers.add_parser("dashboard", help="write the findings page")
    d.add_argument(
        "--db",
        type=Path,
        default=reports / "genomelab.duckdb",
        help="the findings store (default: reports/genomelab.duckdb)",
    )
    d.add_argument("--run", metavar="RUN_ID", help="a run id (default: the latest run)")
    d.add_argument(
        "--cache",
        type=Path,
        default=reports / "explanations",
        help="the explanations' cache (default: reports/explanations/)",
    )
    d.add_argument(
        "--out",
        type=Path,
        default=reports / "dashboard.html",
        help="the page (default: reports/dashboard.html)",
    )
    d.add_argument(
        "--lang",
        choices=LANGUAGES,
        help="the page's language (default: $GENOKEEPER_LANG, else en)",
    )
    d.set_defaults(handler=_dashboard_command)
