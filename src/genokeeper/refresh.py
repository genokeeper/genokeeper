"""`genokeeper refresh`: take each refreshed source's newest release and grade the user's
calls on it again (spec, "The refresh").

In this order, under `<work>/.refresh.lock` and umask 077:

1. A pre-run on the releases already on disk, only when the analyser differs from the
   latest snapshot's (or that snapshot records none).
2. Each refreshed source fetched by `scripts/bootstrap.sh --refresh STEP`, recorded new,
   unchanged or failed.
3. A run when a grading input's newest file on disk, the gnomAD FAF manifest or a
   calibration table's `versions.json` differs (by sha256) from what the latest snapshot
   recorded, or the analyser still differs.
4. The publications of the trusted findings and the watch-list VUS, the explanations (only
   with `--explain`), the dashboard and the retention of dated folders.
5. One record in `<lab>/refresh/log.jsonl`.

A run that raises fails the refresh, as does any other error, which is still recorded: its
full text goes to `reports/refresh/<started>.error.txt`, and the log gets its class only.
No fetch names anything of the user's: `fetch` gets a name from `REFRESHED`. Explaining
sends each new finding's data to Anthropic (`genokeeper explain`), so it runs only when
the user passes `--explain`.
"""

import argparse
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from genokeeper import history, refresh_log
from genokeeper.calibrate.tables import versions_files
from genokeeper.checksum import sha256_file
from genokeeper.config import Paths, load_paths
from genokeeper.dashboard import write_dashboard
from genokeeper.evidence.cgd import newest_cgd
from genokeeper.evidence.clingen import DOSAGE_FILE
from genokeeper.evidence.frequency import MANIFEST as FAF_MANIFEST
from genokeeper.evidence.frequency import faf_folder
from genokeeper.evidence.gencc import newest_gencc
from genokeeper.evidence.gwas import DOWNLOAD as GWAS_DOWNLOAD
from genokeeper.evidence.pgx import CPIC_FILES
from genokeeper.evidence.publications import (
    CITATIONS_FILE,
    LITVAR_FILE,
    cited_pmids,
    litvar_pmids,
)
from genokeeper.explain import (
    ExplainError,
    claude_ask,
    explain,
    newest_medlineplus,
    to_explain,
)
from genokeeper.history import Snapshot
from genokeeper.identity import analyser_identity
from genokeeper.lang import language
from genokeeper.pipeline import (
    STORE_NAME,
    RunRefused,
    _newest_dated,
    fastq_pair,
    newest_clinvar38,
    newest_gene_validity,
    newest_pharmgkb,
    open_faf_table,
    run_genome,
)
from genokeeper.query import Filters, open_store
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus
from genokeeper.store.db import Kind

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

REFRESHED: tuple[str, ...] = (
    "clinvar",
    "clingen",
    "clingen_dosage",
    "cpic",
    "pharmgkb",
    "gwas",
    "medlineplus",
    "litvar",
    "clinvar_citations",
    "gencc",
)
# Each refreshed source's folder under the lab, holding one `<YYYYMMDD>/` per release.
FOLDERS: dict[str, str] = {
    "clinvar": "clinvar",
    "clingen": "clingen",
    "clingen_dosage": "clingen-dosage",
    "cpic": "cpic",
    "pharmgkb": "pharmgkb",
    "gwas": "gwas",
    "medlineplus": "medlineplus",
    "litvar": "litvar",
    "clinvar_citations": "clinvar-citations",
    "gencc": "gencc",
}
# The `run.inputs` keys a refresh compares with the newest files on disk ("pgx.cpic" is
# `inputs["pgx"]["cpic"]`). MedlinePlus and the publication sources are not among them.
GRADING_INPUTS: tuple[str, ...] = (
    "clinvar38",
    "gene_validity",
    "cgd",
    "pgx.cpic",
    "pgx.pharmgkb",
    "gwas",
    "clingen_dosage",
    "gencc",
)
# The `run.inputs` keys of what the lab builds itself, compared the same way: the gnomAD FAF
# table's manifest (`genokeeper build-faf`) and each calibration table's `versions.json`
# (`genokeeper calibrate`), so a refresh after either grades on it.
LAB_TABLES: tuple[str, ...] = ("faf", "tables")
LOCK_NAME = ".refresh.lock"  # under reports/work, beside the run's own lock
MANIFEST = "bootstrap-manifest.tsv"  # the bootstrap's step, path and sha256 per file
EXIT_CODES = {"complete": 0, "partial": 3, "failed": 1}
# Seconds per source: LitVar is about 2 GB; a stalled fetch is then a failed source.
FETCH_TIMEOUT = 3 * 3600

type Fetch = Callable[[str], subprocess.CompletedProcess[str]]
type Run = Callable[[], str]

_DATED = re.compile(r"[0-9]{8}")
_DIGITS = re.compile(r"[0-9]+")


class RefreshRefused(RuntimeError):
    """Another refresh holds `<work>/.refresh.lock`."""


@contextmanager
def refresh_lock(work: Path) -> Iterator[None]:
    """Hold `<work>/.refresh.lock` for the block; refuse (`RefreshRefused`) while another
    refresh holds it."""
    work.mkdir(parents=True, exist_ok=True)
    lock = work / LOCK_NAME
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RefreshRefused(f"another refresh holds {lock}") from None
        yield
    finally:
        os.close(fd)


# --- the sources ----------------------------------------------------------------------------


def _dated(folder: Path) -> list[Path]:
    """The `<YYYYMMDD>/` folders under `folder`, oldest first."""
    if not folder.is_dir():
        return []
    return sorted(
        d for d in folder.iterdir() if d.is_dir() and _DATED.fullmatch(d.name)
    )


def _manifest_rows(lab: Path, step: str) -> list[str]:
    path = lab / MANIFEST
    if not path.is_file():
        return []
    rows = path.read_text(encoding="utf-8").splitlines()
    return sorted(row for row in rows if row.split("\t", 1)[0] == step)


def _filed(
    done: subprocess.CompletedProcess[str], name: str, folder: Path
) -> str | None:
    """The dated folder the bootstrap's `refreshed: STEP FOLDER` line names, when it is
    one of `folder`'s."""
    for line in done.stdout.splitlines():
        found = line.removeprefix(f"refreshed: {name} ")
        if found == line:
            continue
        path = Path(found)
        if (
            _DATED.fullmatch(path.name)
            and path.is_dir()
            and path.parent.resolve() == folder.resolve()
        ):
            return path.name
    return None


def _last_line(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines[-1] if lines else ""


def fetch_sources(
    lab: Path, fetch: Fetch, details: dict[str, str] | None = None
) -> tuple[SourceOutcome, ...]:
    """Fetch each of `REFRESHED` through `fetch(name)` and say what came of it: `new` when
    the bootstrap filed another dated folder or other files (its manifest rows changed),
    `unchanged` otherwise, `failed` with an error class when it exited non-zero or named no
    dated folder. A failure's detail (stderr's last line, which may carry tool output) goes
    into `details`, never into the outcome."""
    outcomes: list[SourceOutcome] = []
    for name in REFRESHED:
        folder = lab / FOLDERS[name]
        dated = _dated(folder)
        before = dated[-1].name if dated else None
        rows = _manifest_rows(lab, name)
        after: str | None = None
        try:
            done = fetch(name)
        except Exception as e:  # noqa: BLE001 - recorded as the source's failure
            error, detail = type(e).__name__, " ".join(str(e).split())
        else:
            after = _filed(done, name, folder)
            detail = _last_line(done.stderr)
            if done.returncode != 0:
                error = f"bootstrap:exit:{done.returncode}"
            elif after is None:
                error = "bootstrap:no-folder"
                detail = f"no `refreshed: {name}` line naming a folder of {folder}"
            else:
                error = None
        if error is not None:
            outcomes.append(
                SourceOutcome(name, before, None, SourceStatus.FAILED, error)
            )
            if details is not None:
                details[name] = detail
            continue
        new = after != before or _manifest_rows(lab, name) != rows
        status = SourceStatus.NEW if new else SourceStatus.UNCHANGED
        outcomes.append(SourceOutcome(name, before, after, status))
    return tuple(outcomes)


# --- when a run is needed -------------------------------------------------------------------


def _refused_is_none(find: Callable[[Path], Path], lab: Path) -> Path | None:
    try:
        return find(lab)
    except RunRefused:
        return None


def _in(folder: Path | None, *names: str) -> dict[str, Path] | None:
    return None if folder is None else {name: folder / name for name in names}


def _one(path: Path | None) -> dict[str, Path] | None:
    return None if path is None else {"": path}


def _member(folder: Path, name: str) -> dict[str, Path] | None:
    """`name` in the newest dated folder of `folder` holding it, as one file input."""
    newest = _newest_dated(folder, [name])
    return _one(None if newest is None else newest / name)


# Each grading input's and lab table's newest file(s) on disk, found as `run_genome` finds
# them.
_NEWEST: dict[str, Callable[[Path], dict[str, Path] | None]] = {
    "clinvar38": lambda lab: _one(_refused_is_none(newest_clinvar38, lab)),
    "gene_validity": lambda lab: _one(_refused_is_none(newest_gene_validity, lab)),
    "cgd": lambda lab: _one(newest_cgd(lab)),
    "pgx.cpic": lambda lab: _in(
        _newest_dated(lab / "cpic", list(CPIC_FILES.values())), *CPIC_FILES.values()
    ),
    "pgx.pharmgkb": lambda lab: _one(newest_pharmgkb(lab)),
    "gwas": lambda lab: _member(lab / "gwas", GWAS_DOWNLOAD),
    "clingen_dosage": lambda lab: _member(lab / "clingen-dosage", DOSAGE_FILE),
    "gencc": lambda lab: _one(newest_gencc(lab)),
    "faf": lambda lab: _one(faf_folder(lab) / FAF_MANIFEST),
    "tables": versions_files,
}


def _on_disk(lab: Path, name: str) -> dict[str, str] | None:
    files = _NEWEST[name](lab)
    if files is None or not all(p.is_file() for p in files.values()):
        return None
    return {member: sha256_file(path) for member, path in files.items()}


def _recorded(inputs: Mapping[str, Any], name: str) -> dict[str, str] | None:
    """The sha256s `run.inputs` holds for `name`: `{"": sha}` for a file, `{member: sha}`
    for a map of files; None for an input the run recorded absent (GenCC's
    `{"present": false, ...}`), as for one it never recorded."""
    node: Any = inputs
    for part in name.split("."):
        node = node.get(part) if isinstance(node, dict) else None
    if not isinstance(node, dict):
        return None
    if "sha256" in node:
        return {"": str(node["sha256"])}
    if node.get("present") is False:
        return None
    return {
        str(k): str(v["sha256"])
        for k, v in node.items()
        if isinstance(v, dict) and "sha256" in v
    }


def run_needed(lab: Path, latest: Snapshot | None, identity: str) -> bool:
    """Whether a run on the newest releases is due: no snapshot yet, an analyser other than
    the latest snapshot's, or a grading input or lab table whose newest file on disk has
    another sha256 than the one that snapshot recorded (a missing one included)."""
    if latest is None or latest.tools.get("analyser") != identity:
        return True
    return any(
        _on_disk(lab, name) != _recorded(latest.inputs, name)
        for name in GRADING_INPUTS + LAB_TABLES
    )


# --- publications ---------------------------------------------------------------------------


def _newest_file(lab: Path, source: str, name: str) -> Path | None:
    found = [d / name for d in _dated(lab / FOLDERS[source]) if (d / name).is_file()]
    return found[-1] if found else None


def _watch_list(root: Path, run_id: str) -> frozenset[str]:
    path = history._parquet(root, run_id, "statement")
    with duckdb.connect() as con:
        rows = con.execute(
            "SELECT statement_id FROM read_parquet(?) WHERE kind = ?",
            [path, str(Kind.WATCH)],
        ).fetchall()
    return frozenset(row[0] for row in rows)


def _pmids(found: frozenset[str] | None) -> list[str]:
    return sorted(found or (), key=int)


def publications(lab: Path, root: Path, run_id: str) -> dict[str, dict[str, list[str]]]:
    """For each trusted finding (at the default filters) and each watch-list VUS of the
    snapshot `run_id`: the PMIDs the newest LitVar file lists for its rsID and the PMIDs
    the newest ClinVar citations attach to its VariationID. Both files are read locally."""
    ids = history.trusted_ids(root, run_id, Filters()) | _watch_list(root, run_id)
    states = history.states(root, run_id, ids)
    rsids = {sid: s.rsid for sid, s in states.items() if s.rsid}
    variations = {
        sid: sid.removeprefix("clinvar:")
        for sid in ids
        if sid.startswith("clinvar:")
        and _DIGITS.fullmatch(sid.removeprefix("clinvar:"))
    }
    litvar = _newest_file(lab, "litvar", LITVAR_FILE)
    cited = _newest_file(lab, "clinvar_citations", CITATIONS_FILE)
    by_rsid = litvar_pmids(litvar, set(rsids.values())) if litvar else {}
    by_variation = cited_pmids(cited, set(variations.values())) if cited else {}
    return {
        sid: {
            "litvar": _pmids(by_rsid.get(rsids.get(sid, ""))),
            "clinvar": _pmids(by_variation.get(variations.get(sid, ""))),
        }
        for sid in sorted(ids)
    }


def write_publications(
    reports: Path, started: datetime, pubs: Mapping[str, Mapping[str, list[str]]]
) -> None:
    """The refresh's own `publications.json` (`refresh_log.publications_path`), as
    `history.read_publications` reads it."""
    path = refresh_log.publications_path(reports, started)
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(pubs, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


# --- retention ------------------------------------------------------------------------------


def _paths(node: Any) -> Iterator[Path]:
    """Every `path` a `run.inputs` document records, at any depth."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "path" and isinstance(value, str):
                yield Path(value)
            else:
                yield from _paths(value)
    elif isinstance(node, list):
        for value in node:
            yield from _paths(value)


def retain(lab: Path, records: Sequence[RefreshRecord], root: Path) -> list[Path]:
    """Delete each refreshed source's dated folder that is not its newest and that neither
    of the last two complete or partial refreshes' runs read; return what was deleted.
    A run records no MedlinePlus, LitVar or citations file, so those keep their newest
    only."""
    runs = [r.run for r in records if r.outcome != "failed" and r.run is not None]
    kept = list(dict.fromkeys(reversed(runs)))[:2]
    by_id = {s.run_id: s for s in history.snapshots(root)}
    read = [p for rid in kept if rid in by_id for p in _paths(by_id[rid].inputs)]
    deleted: list[Path] = []
    for name in REFRESHED:
        for folder in _dated(lab / FOLDERS[name])[:-1]:
            if any(p.is_relative_to(folder) for p in read):
                continue
            shutil.rmtree(folder)
            deleted.append(folder)
    return deleted


# --- the refresh ----------------------------------------------------------------------------


type Failure = tuple[str, BaseException | str]


def _latest(root: Path) -> Snapshot | None:
    found = history.snapshots(root)
    return found[-1] if found else None


def _base(lab: Path, latest: Snapshot | None) -> str | None:
    """The comparison base, decided before anything runs: the run of the previous complete
    or partial refresh, else the latest snapshot."""
    runs = [r.run for r in refresh_log.read(lab) if r.outcome != "failed" and r.run]
    return runs[-1] if runs else (latest.run_id if latest is not None else None)


def error_path(reports: Path, started: datetime) -> Path:
    """`reports/refresh/<started>.error.txt`, the refresh's failures in full."""
    return reports / "refresh" / f"{refresh_log.stamp(started)}.error.txt"


def _write_errors(
    reports: Path, started: datetime, failures: Sequence[Failure]
) -> Path:
    """Every failure's full text, which can carry tool output naming the user's variants,
    into `reports/refresh/<started>.error.txt`."""
    path = error_path(reports, started)
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    parts = [
        f"{label}: {failure}\n"
        if isinstance(failure, str)
        else f"== {label}\n{''.join(traceback.format_exception(failure))}"
        for label, failure in failures
    ]
    path.write_text("".join(parts), encoding="utf-8")
    return path


def _steps(
    record: RefreshRecord,
    failures: list[Failure],
    steps: Mapping[str, Callable[[], Any]],
) -> RefreshRecord:
    """Run each step after a run; a step that raises leaves its text in `failures` and,
    for the first, `step:Class` as the record's error, and does not fail the refresh."""
    error = record.error
    for label, step in steps.items():
        try:
            step()
        except Exception as e:  # noqa: BLE001 - recorded; the refresh goes on
            failures.append((label, e))
            error = error or f"{label}:{type(e).__name__}"
    return replace(record, error=error)


def _refresh(
    paths: Paths,
    started: datetime,
    failures: list[Failure],
    *,
    fetch: Fetch,
    run: Run,
    explain: Callable[[str], None] | None,
    dashboard: Callable[[str], None],
    now: Callable[[], datetime],
) -> RefreshRecord:
    root = paths.reports / "snapshots"
    identity = analyser_identity()
    latest = _latest(root)
    base = _base(paths.lab, latest)
    pre_run: str | None = None
    made: str | None = None
    refused: str | None = None
    if latest is None or latest.tools.get("analyser") != identity:
        try:
            pre_run = run()
        except Exception as e:  # noqa: BLE001 - the refresh fails; its class is logged
            failures.append(("pre-run", e))
            refused = type(e).__name__
        latest = _latest(root)
    details: dict[str, str] = {}
    sources = fetch_sources(paths.lab, fetch, details)
    failures.extend(details.items())
    if run_needed(paths.lab, latest, identity):
        try:
            made = run()
        except Exception as e:  # noqa: BLE001 - the refresh fails; its class is logged
            failures.append(("run", e))
            refused = refused or type(e).__name__
    latest = _latest(root)
    current = made or pre_run or (latest.run_id if latest is not None else None)
    if refused is not None or current is None:
        # No Changes table for a failed refresh; the previous refresh stays the base.
        return RefreshRecord(
            started, now(), "failed", sources, None, pre_run, None, refused or "NoRun"
        )
    failed = any(s.status is SourceStatus.FAILED for s in sources)
    record = RefreshRecord(
        started=started,
        ended=started,
        outcome="partial" if failed else "complete",
        sources=sources,
        base_run=base,
        pre_run=pre_run,
        run=current,
    )
    steps: dict[str, Callable[[], Any]] = {
        "publications": lambda: write_publications(
            paths.reports, started, publications(paths.lab, root, current)
        )
    }
    if explain is not None:
        steps["explain"] = lambda: explain(current)
    steps["dashboard"] = lambda: dashboard(current)
    record = _steps(record, failures, steps)
    record = _steps(
        record,
        failures,
        {
            "retain": lambda: retain(
                paths.lab, [*refresh_log.read(paths.lab), record], root
            )
        },
    )
    return replace(record, ended=now())


def refresh(
    paths: Paths,
    *,
    fetch: Fetch,
    run: Run,
    explain: Callable[[str], None] | None,
    dashboard: Callable[[str], None],
    now: Callable[[], datetime],
) -> RefreshRecord:
    """Refresh the sources and grade on them (the module's order), append the record to
    the log and return it."""
    started = now()
    failures: list[Failure] = []
    # Everything a refresh writes (snapshots, publications, error texts) is the user's.
    previous = os.umask(0o077)
    try:
        try:
            with refresh_lock(paths.reports / "work"):
                try:
                    record = _refresh(
                        paths,
                        started,
                        failures,
                        fetch=fetch,
                        run=run,
                        explain=explain,
                        dashboard=dashboard,
                        now=now,
                    )
                # Recorded; a Ctrl-C is then raised again, after its log record.
                except (Exception, KeyboardInterrupt) as e:
                    failures.append(("refresh", e))
                    record = RefreshRecord(
                        started, now(), "failed", (), None, None, None, type(e).__name__
                    )
                    if isinstance(e, KeyboardInterrupt):
                        _finish(paths, started, failures, record)
                        raise
                _finish(paths, started, failures, record)
        except RefreshRefused as e:
            failures.append(("lock", e))
            record = RefreshRecord(
                started, now(), "failed", (), None, None, None, type(e).__name__
            )
            _finish(paths, started, failures, record)
    finally:
        os.umask(previous)
    return record


def _finish(
    paths: Paths, started: datetime, failures: Sequence[Failure], record: RefreshRecord
) -> None:
    if failures:
        _write_errors(paths.reports, started, failures)
    refresh_log.append(paths.lab, record)


# --- command line ---------------------------------------------------------------------------


def _command(args: argparse.Namespace) -> int:
    paths = load_paths()
    script = paths.repo / "scripts" / "bootstrap.sh"

    def fetch(name: str) -> subprocess.CompletedProcess[str]:
        # A fetch can take minutes (LitVar is about 2 GB), so each one is announced.
        print(f"genokeeper refresh: fetching {name}", file=sys.stderr, flush=True)
        return subprocess.run(
            ["bash", str(script), "--refresh", name],
            capture_output=True,
            text=True,
            check=False,
            timeout=FETCH_TIMEOUT,
        )

    def run() -> str:
        # As `genokeeper run` calls it, with every default.
        print(
            "genokeeper refresh: grading the calls (cached steps are reused)",
            file=sys.stderr,
            flush=True,
        )
        try:
            fastq1, fastq2 = fastq_pair(paths.raw)
            return run_genome(
                paths,
                fastq1=fastq1,
                fastq2=fastq2,
                clinvar38=newest_clinvar38(paths.lab),
                faf=open_faf_table(paths.lab),
            )
        except Exception as e:
            print(
                f"genokeeper refresh: the run was refused: {type(e).__name__}: {e}",
                file=sys.stderr,
                flush=True,
            )
            raise

    def explain_run(run_id: str) -> None:
        # As `genokeeper explain` does: only findings and stage-3 watch items without a
        # cached answer are asked, in the language `$GENOKEEPER_LANG` chooses.
        with open_store(paths.reports / STORE_NAME) as con:
            findings = to_explain(con, run_id)
        medlineplus = newest_medlineplus(paths.lab)
        if medlineplus is None:
            raise FileNotFoundError(
                f"no MedlinePlus summaries in {paths.lab}/medlineplus"
            )
        counts = explain(
            findings,
            medlineplus,
            paths.reports / "explanations",
            claude_ask,
            lang=language(),
        )
        if counts["failed"]:
            # As `genokeeper explain` exits 1: the step fails and the refresh records it.
            raise ExplainError(f"{counts['failed']} explain calls failed")

    record = refresh(
        paths,
        fetch=fetch,
        run=run,
        explain=explain_run if args.explain else None,
        dashboard=lambda run_id: print(
            f"genokeeper refresh: wrote {write_dashboard(paths, run_id)}"
        ),
        now=lambda: datetime.now(UTC),
    )
    for s in record.sources:
        version = s.after if s.status is not SourceStatus.FAILED else s.error
        print(f"{s.name}: {s.status} ({s.before} -> {version})")
    print(
        f"refresh {record.outcome}: run {record.run}, compared with {record.base_run}"
    )
    errors = error_path(paths.reports, record.started)
    if errors.is_file():
        print(f"genokeeper refresh: what failed is in {errors}", file=sys.stderr)
    return EXIT_CODES[record.outcome]


def register(subparsers: Subparsers) -> None:
    """Add the `refresh` sub-command."""
    r = subparsers.add_parser(
        "refresh",
        help="fetch the knowledge bases' newest releases and grade the calls on them",
    )
    r.add_argument(
        "--explain",
        action="store_true",
        help="also explain new findings: sends their data to Anthropic, as `genokeeper"
        " explain` does",
    )
    r.set_defaults(handler=_command)
