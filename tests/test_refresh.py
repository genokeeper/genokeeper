"""`genokeeper refresh`: the lock, the pre-run, the fetch, the run when a grading input
changed, publications, explanations, the dashboard, retention and the log record.

Every collaborator is a fake: `fetch` files dated folders in a tmp lab as the bootstrap's
`--refresh STEP` does, `run` exports a `snapshots_mini` snapshot recording the newest files
on disk as `run_genome` records them, and `explain`/`dashboard` log their calls.
"""

import contextlib
import fcntl
import gzip
import hashlib
import os
import re
import shlex
import shutil
import subprocess
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import pytest
from fixtures import snapshots_mini as sm
from fixtures.snapshots_mini import call, run_id, started, statement

from genokeeper import changes, cli, refresh
from genokeeper.calibrate import tables
from genokeeper.config import Paths
from genokeeper.evidence.cgd import FILE_NAME as CGD_FILE
from genokeeper.evidence.clingen import DOSAGE_FILE
from genokeeper.evidence.clingen import FILE_NAME as VALIDITY_FILE
from genokeeper.evidence.frequency import MANIFEST as FAF_MANIFEST
from genokeeper.evidence.frequency import faf_folder
from genokeeper.evidence.gencc import GENCC_FILE
from genokeeper.evidence.gwas import DOWNLOAD as GWAS_FILE
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE
from genokeeper.evidence.pgx import CPIC_FILES
from genokeeper.evidence.publications import CITATIONS_FILE, LITVAR_FILE
from genokeeper.explain import MODEL, ExplainError
from genokeeper.history import read_publications, snapshots
from genokeeper.pipeline import RunRefused
from genokeeper.query import Filters
from genokeeper.refresh import GRADING_INPUTS, REFRESHED
from genokeeper.refresh_log import RefreshRecord, publications_path, read

FIXTURES = Path(__file__).parent / "fixtures"
BOOTSTRAP = Path(__file__).parents[1] / "scripts" / "bootstrap.sh"
IDENTITY = "a" * 64
# The source -> its lab folder and the file each dated folder holds.
FOLDERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "clinvar": ("clinvar", ("clinvar_GRCh38.vcf.gz",)),
    "clingen": ("clingen", (VALIDITY_FILE,)),
    "clingen_dosage": ("clingen-dosage", (DOSAGE_FILE,)),
    "cpic": ("cpic", tuple(CPIC_FILES.values())),
    "pharmgkb": ("pharmgkb", ("summaryAnnotations.zip",)),
    "gwas": ("gwas", (GWAS_FILE,)),
    "medlineplus": ("medlineplus", (MEDLINEPLUS_FILE,)),
    "litvar": ("litvar", (LITVAR_FILE,)),
    "clinvar_citations": ("clinvar-citations", (CITATIONS_FILE,)),
    "gencc": ("gencc", (GENCC_FILE,)),
}
FIRST = "20261001"
EXTRA_PAPER = (
    '{"_id": "litvar@rs397507801#x", "rsid": "rs397507801", "pmids": [99999999]}\n'
)
# A trusted ClinVar finding at VariationID 2 (cited in the citations fixture) and a
# watch-list VUS whose rsID LitVar lists; an untrusted finding at a low call.
VID2, VID3, VID4 = (f"GRCh38:chr17:{p}:C:T" for p in (2002, 3003, 4004))
STATEMENTS = [
    statement("clinvar:2", vid=VID2, rsid="rs397704705"),
    statement("clinvar:3", vid=VID3, rsid="rs397507801", kind="watch", validity=None),
    statement("clinvar:4", vid=VID4, rsid="rs746985154"),
]
CALLS = [call(VID2), call(VID3), call(VID4, level="low")]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def put(lab: Path, source: str, date: str, content: str = "v1") -> Path:
    """A dated folder of `source` holding its files."""
    folder, names = FOLDERS[source]
    dest = lab / folder / date
    dest.mkdir(parents=True, exist_ok=True)
    for name in names:
        if source == "litvar":
            shutil.copy(FIXTURES / "litvar_head.jsonl.gz", dest / name)
            if date != FIRST:  # a later release lists one more paper for clinvar:3
                with gzip.open(dest / name, "at") as out:
                    out.write(EXTRA_PAPER)
        elif source == "clinvar_citations":
            shutil.copy(FIXTURES / "var_citations_head.txt", dest / name)
        else:
            (dest / name).write_text(f"{source} {name} {content}\n")
    return dest


def newest(lab: Path, source: str) -> Path:
    folder, names = FOLDERS[source]
    dated = sorted(
        d for d in (lab / folder).iterdir() if all((d / n).is_file() for n in names)
    )
    return dated[-1]


def recorded(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": sha(path)}


def inputs_of(lab: Path) -> dict[str, Any]:
    """`run.inputs` as `run_genome` records the newest grading inputs on disk."""

    def one(source: str) -> Path:
        return newest(lab, source) / FOLDERS[source][1][0]

    cpic = newest(lab, "cpic")
    return {
        "clinvar38": recorded(one("clinvar")),
        "gene_validity": {**recorded(one("clingen")), "date": FIRST, "present": True},
        "cgd": {**recorded(lab / "cgd" / "20250201" / CGD_FILE), "present": True},
        "pgx": {
            "cpic": {n: recorded(cpic / n) for n in CPIC_FILES.values()},
            "pharmgkb": {**recorded(one("pharmgkb")), "created": FIRST},
        },
        "gwas": {**recorded(one("gwas")), "release": FIRST},
        "clingen_dosage": {**recorded(one("clingen_dosage")), "date": FIRST},
        "gencc": {**recorded(one("gencc")), "date": FIRST},
        **(
            {"faf": recorded(faf)}
            if (faf := faf_folder(lab) / FAF_MANIFEST).is_file()
            else {}
        ),
        "tables": {
            name: recorded(tables.folder(lab, name) / tables.VERSIONS)
            for name in tables.TABLES
            if (tables.folder(lab, name) / tables.VERSIONS).is_file()
        },
        "steps": {},
    }


@dataclass
class Lab:
    """A tmp repository and lab, with fakes for every collaborator of `refresh`."""

    paths: Paths
    events: list[str] = field(default_factory=list)
    fetched: list[str] = field(default_factory=list)
    # source -> "same" | "rewrite" | "fail" | a new dated folder name
    plan: dict[str, str] = field(default_factory=dict)
    runs: list[str] = field(default_factory=list)
    refuse: list[BaseException] = field(default_factory=list)
    explain_error: BaseException | None = None
    umasks: list[int] = field(default_factory=list)
    clock: datetime = datetime(2026, 10, 27, 9, 0, tzinfo=UTC)
    day: int = 10
    statements: list[dict[str, Any]] = field(default_factory=lambda: list(STATEMENTS))

    @property
    def lab(self) -> Path:
        return self.paths.lab

    @property
    def root(self) -> Path:
        return self.paths.reports / "snapshots"

    def snapshot(self, analyser: str | None = IDENTITY) -> str:
        self.day += 1
        rid = run_id(self.day)
        tools = {"clinvar": FIRST} | (
            {} if analyser is None else {"analyser": analyser}
        )
        sm.write(
            self.root,
            rid,
            started(rid),
            tools,
            inputs_of(self.lab),
            self.statements,
            CALLS,
        )
        return rid

    def fetch(self, name: str) -> subprocess.CompletedProcess[str]:
        self.events.append(f"fetch:{name}")
        self.fetched.append(name)
        how = self.plan.get(name, "same")
        args = ["bash", "bootstrap.sh", "--refresh", name]
        if how == "fail":
            return subprocess.CompletedProcess(
                args, 1, "", "progress\nbootstrap: curl failed for BRCA1 chr17:2002\n"
            )
        if how == "rewrite":  # the same dated folder, filed again with other bytes
            dest = put(self.lab, name, newest(self.lab, name).name, "rewritten")
            with (self.lab / "bootstrap-manifest.tsv").open("a") as manifest:
                for path in sorted(dest.iterdir()):
                    manifest.write(f"{name}\t{path}\t{sha(path)}\n")
        elif how == "same":
            dest = newest(self.lab, name)
        else:
            dest = put(self.lab, name, how, how)
        return subprocess.CompletedProcess(args, 0, f"refreshed: {name} {dest}\n", "")

    def run(self) -> str:
        self.events.append("run")
        self.umasks.append(_umask())
        if self.refuse:
            raise self.refuse.pop(0)
        rid = self.snapshot()
        self.runs.append(rid)
        return rid

    def explain(self, rid: str) -> None:
        self.events.append(f"explain:{rid}")
        if self.explain_error is not None:
            raise self.explain_error

    def dashboard(self, rid: str) -> None:
        self.events.append(f"dashboard:{rid}")

    def now(self) -> datetime:
        self.clock += timedelta(minutes=1)
        return self.clock

    def refresh(self, *, explain: bool = True) -> RefreshRecord:
        return refresh.refresh(
            self.paths,
            fetch=self.fetch,
            run=self.run,
            explain=self.explain if explain else None,
            dashboard=self.dashboard,
            now=self.now,
        )

    def changes(self) -> dict[str, Any]:
        """`genokeeper changes --all --json` at its defaults."""
        return changes.report(
            self.paths.reports,
            self.lab,
            base=None,
            run=None,
            f=Filters(),
            every=True,
            gene=None,
            marker=None,
            text=None,
        )


def _umask() -> int:
    current = os.umask(0)
    os.umask(current)
    return current


@pytest.fixture
def lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Lab]:
    paths = Paths(
        repo=tmp_path,
        raw=tmp_path / "raw",
        reports=tmp_path / "reports",
        lab=tmp_path / "lab",
    )
    for source in REFRESHED:
        put(paths.lab, source, FIRST)
    (paths.lab / "cgd" / "20250201").mkdir(parents=True)
    (paths.lab / "cgd" / "20250201" / CGD_FILE).write_text("cgd\n")
    monkeypatch.setattr(refresh, "analyser_identity", lambda: IDENTITY)
    yield Lab(paths)


def statuses(record: RefreshRecord) -> dict[str, str]:
    return {s.name: str(s.status) for s in record.sources}


def lab_text(lab: Path) -> str:
    """Every file `refresh` writes under the lab (outside the fetched sources)."""
    return "".join(p.read_text() for p in (lab / "refresh").rglob("*") if p.is_file())


def test_first_refresh_ever(lab: Lab) -> None:
    # Review focus 2: the only snapshot was written before the analyser was recorded.
    first = lab.snapshot(analyser=None)
    record = lab.refresh()
    assert lab.events[0] == "run"  # the pre-run, before any fetch
    assert lab.events[1 : 1 + len(REFRESHED)] == [f"fetch:{name}" for name in REFRESHED]
    assert lab.events.count("run") == 1  # nothing new: no second run
    (pre,) = lab.runs
    assert record.outcome == "complete"
    assert (record.pre_run, record.run, record.base_run) == (pre, pre, first)
    assert record.error is None
    assert lab.events[-2:] == [f"explain:{pre}", f"dashboard:{pre}"]
    assert lab.umasks == [0o077]
    # Publications: the trusted finding and the watch-list VUS, not the untrusted one.
    path = publications_path(lab.paths.reports, record.started)
    pubs = read_publications(path)
    assert set(pubs) == {"clinvar:2", "clinvar:3"}
    assert pubs["clinvar:2"]["clinvar"] == [
        "20613862",
        "25741868",
        "26467025",
        "40081374",
    ]
    assert pubs["clinvar:3"]["litvar"] == [
        "12552570",
        "15172753",
        "32377563",
        "36922933",
    ]
    assert path.stat().st_mode & 0o077 == 0
    # A first list is a baseline: the change feed shows no publications change.
    report = lab.changes()
    assert (report["base"], report["run"]) == (first, pre)
    assert not [c for c in report["changes"] if c["publications"]["new"]]
    assert read(lab.lab) == [record]


def test_nothing_new_runs_nothing(lab: Lab) -> None:
    latest = lab.snapshot()
    record = lab.refresh()
    assert "run" not in lab.events
    assert record.outcome == "complete"
    assert set(statuses(record).values()) == {"unchanged"}
    assert (record.pre_run, record.run) == (None, latest)
    assert {s.before for s in record.sources} == {FIRST}
    assert lab.events[-2:] == [f"explain:{latest}", f"dashboard:{latest}"]


def test_medlineplus_or_litvar_alone_runs_nothing(lab: Lab) -> None:
    lab.snapshot()
    lab.plan = {"medlineplus": "20261020", "litvar": "20261026"}
    record = lab.refresh()
    assert "run" not in lab.events  # neither is a grading input
    assert statuses(record)["medlineplus"] == statuses(record)["litvar"] == "new"
    litvar = next(s for s in record.sources if s.name == "litvar")
    assert (litvar.before, litvar.after) == (FIRST, "20261026")
    assert record.outcome == "complete"


def test_new_clinvar_runs_once(lab: Lab) -> None:
    previous = lab.snapshot()
    lab.plan = {"clinvar": "20261018"}
    record = lab.refresh()
    assert lab.events.count("run") == 1
    assert lab.events.index("run") > lab.events.index("fetch:clinvar_citations")
    (rid,) = lab.runs
    assert (record.pre_run, record.run, record.base_run) == (None, rid, previous)
    assert statuses(record)["clinvar"] == "new"
    assert {n for n, s in statuses(record).items() if s == "unchanged"} == set(
        REFRESHED
    ) - {"clinvar"}
    latest = snapshots(lab.root)[-1]
    assert latest.inputs["clinvar38"]["path"].endswith("20261018/clinvar_GRCh38.vcf.gz")


def test_same_folder_filed_again_with_other_bytes_is_new(lab: Lab) -> None:
    # CPIC is dated by its fetch day: a second fetch that day overwrites the folder.
    lab.snapshot()
    lab.plan = {"cpic": "rewrite"}
    record = lab.refresh()
    cpic = next(s for s in record.sources if s.name == "cpic")
    assert (cpic.before, cpic.after, str(cpic.status)) == (FIRST, FIRST, "new")
    assert lab.events.count("run") == 1


@pytest.mark.parametrize("name", GRADING_INPUTS)
def test_each_grading_input_on_disk_decides_a_run(lab: Lab, name: str) -> None:
    lab.snapshot()
    latest = snapshots(lab.root)[-1]
    assert not refresh.run_needed(lab.lab, latest, IDENTITY)
    assert refresh.run_needed(lab.lab, latest, "b" * 64)
    assert refresh.run_needed(lab.lab, None, IDENTITY)
    source = {
        "clinvar38": "clinvar",
        "gene_validity": "clingen",
        "pgx.cpic": "cpic",
        "pgx.pharmgkb": "pharmgkb",
        "gwas": "gwas",
        "clingen_dosage": "clingen_dosage",
        "gencc": "gencc",
    }.get(name)
    if source is None:  # cgd: pinned, but a newer table on disk is still read by a run
        (lab.lab / "cgd" / "20260101").mkdir()
        (lab.lab / "cgd" / "20260101" / CGD_FILE).write_text("cgd 2\n")
    else:
        put(lab.lab, source, "20261101", "v2")
    assert refresh.run_needed(lab.lab, latest, IDENTITY)


@pytest.mark.parametrize("table", ["small_variants", "stage3", "slice3f"])
def test_a_new_calibration_table_decides_a_run(lab: Lab, table: str) -> None:
    """`genokeeper calibrate` writes a table with other versions: the next refresh runs,
    so the store stops reading "not calibrated"; a table the snapshot already graded with,
    and a stale one left as it was, run nothing."""
    tables.write_table(lab.lab, table, lambda _: None, {"bwa": "0.7.18"})
    lab.snapshot()
    latest = snapshots(lab.root)[-1]
    assert not refresh.run_needed(lab.lab, latest, IDENTITY)
    tables.write_table(lab.lab, table, lambda _: None, {"bwa": "0.7.19"})
    assert refresh.run_needed(lab.lab, latest, IDENTITY)


def test_a_rebuilt_faf_table_decides_a_run(lab: Lab) -> None:
    folder = faf_folder(lab.lab)
    folder.mkdir(parents=True)
    (folder / FAF_MANIFEST).write_text('{"floor": 0.001}')
    lab.snapshot()
    latest = snapshots(lab.root)[-1]
    assert not refresh.run_needed(lab.lab, latest, IDENTITY)
    (folder / FAF_MANIFEST).write_text('{"floor": 0.0005}')
    assert refresh.run_needed(lab.lab, latest, IDENTITY)


def test_gencc_refresh_registered(lab: Lab) -> None:
    """GenCC is a refreshed source and a grading input: a newer export on disk than the
    one the latest snapshot graded on asks for a run."""
    assert "gencc" in REFRESHED and "gencc" in GRADING_INPUTS
    assert refresh.FOLDERS["gencc"] == "gencc"
    lab.snapshot()
    latest = snapshots(lab.root)[-1]
    assert not refresh.run_needed(lab.lab, latest, IDENTITY)
    put(lab.lab, "gencc", "20261101", "v2")
    assert refresh.run_needed(lab.lab, latest, IDENTITY)


def test_analyser_changed_pre_runs_before_fetching(lab: Lab) -> None:
    previous = lab.snapshot(analyser="b" * 64)
    lab.plan = {"clinvar": "20261018"}
    record = lab.refresh()
    fetches = [f"fetch:{name}" for name in REFRESHED]
    pre, rid = lab.runs
    assert lab.events == [
        "run",
        *fetches,
        "run",
        f"explain:{rid}",
        f"dashboard:{rid}",
    ]
    # The base is what existed before the refresh: the pre-run's changes are reported.
    assert (record.pre_run, record.run, record.base_run) == (pre, rid, previous)
    assert record.outcome == "complete"


# clinvar:2 graded limited: it leaves the trusted tier.
LIMITED = [statement("clinvar:2", vid=VID2, rsid="rs397704705", validity="limited")]


def rows(document: dict[str, Any]) -> list[tuple[str, str, list[str]]]:
    return [(c["id"], c["kind"], c["causes"]) for c in document["changes"]]


@pytest.mark.parametrize("plan", [{}, {"clinvar": "fail"}], ids=["same", "failed"])
def test_a_refresh_that_made_no_run_shows_no_changes(
    lab: Lab, plan: dict[str, str]
) -> None:
    lab.snapshot()
    lab.plan = {"clinvar": "20261018"}
    lab.statements = LIMITED + STATEMENTS[1:]
    first = lab.refresh()
    assert [r[:2] for r in rows(lab.changes())] == [("clinvar:2", "left")]
    lab.plan = plan
    lab.events.clear()
    second = lab.refresh()
    assert "run" not in lab.events
    # The base is the previous refresh's run, which is this one's: nothing to compare.
    assert (second.base_run, second.run) == (rid := first.run, rid)
    report = lab.changes()
    assert (report["base"], report["run"], report["changes"]) == (rid, rid, [])
    assert report["failed"] == (
        [{"name": "clinvar", "before": "20261018"}] if plan else []
    )


def test_first_refresh_ever_with_a_new_source_keeps_the_analysers_changes(
    lab: Lab,
) -> None:
    first = lab.snapshot(analyser=None)
    lab.statements = LIMITED + STATEMENTS[1:]  # the current analyser grades it limited
    lab.plan = {"clinvar": "20261018"}
    record = lab.refresh()
    pre, rid = lab.runs
    assert (record.pre_run, record.run, record.base_run) == (pre, rid, first)
    assert rows(lab.changes()) == [("clinvar:2", "left", ["analyser"])]


def test_publications_are_kept_per_refresh(
    lab: Lab, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab.snapshot()
    first = lab.refresh()
    # A refresh whose publications step fails writes none ...
    real = refresh.publications
    monkeypatch.setattr(refresh, "publications", lambda *a: 1 / 0)
    failed = lab.refresh()
    assert failed.error == "publications:ZeroDivisionError"
    monkeypatch.setattr(refresh, "publications", real)
    # ... and a LitVar-only refresh, which makes no run, compares with the last list.
    lab.plan = {"litvar": "20261026"}
    third = lab.refresh()
    assert third.run == first.run and "run" not in lab.events
    (row,) = lab.changes()["changes"]
    assert (row["id"], row["kind"]) == ("clinvar:3", "publications")
    assert row["publications"] == {"new": ["99999999"], "attached_by_lab": []}
    # Each refresh keeps its own list: the first one's is not overwritten.
    kept = read_publications(publications_path(lab.paths.reports, first.started))
    assert "99999999" not in kept["clinvar:3"]["litvar"]


def test_every_source_failed(lab: Lab) -> None:
    # Review focus 3: the network is down for the whole refresh.
    latest = lab.snapshot()
    lab.plan = dict.fromkeys(REFRESHED, "fail")
    record = lab.refresh()
    assert "run" not in lab.events
    assert record.outcome == "partial"
    assert record.run == latest
    assert set(statuses(record).values()) == {"failed"}
    banner = changes.failed(changes.refresh_of(lab.lab, latest))
    assert banner == [{"name": name, "before": FIRST} for name in REFRESHED]
    assert {s.error for s in record.sources} == {"bootstrap:exit:1"}
    assert lab.events[-1] == f"dashboard:{latest}"


def test_every_source_failed_analyser_changed(lab: Lab) -> None:
    lab.snapshot(analyser="b" * 64)
    lab.plan = dict.fromkeys(REFRESHED, "fail")
    record = lab.refresh()
    assert lab.events.count("run") == 1  # only the pre-run
    assert lab.events[0] == "run"
    (pre,) = lab.runs
    assert (record.outcome, record.pre_run, record.run) == ("partial", pre, pre)
    banner = changes.failed(changes.refresh_of(lab.lab, pre))
    assert [b["name"] for b in banner] == list(REFRESHED)


def test_refused_run_then_next_refresh_runs(lab: Lab) -> None:
    # Review focus 4: a refused run fails the refresh, the previous refresh stays the base,
    # and the next refresh runs although the fetched release is already on disk.
    first = lab.snapshot()
    assert lab.refresh().outcome == "complete"
    lab.events.clear()
    lab.plan = {"clinvar": "20261018"}
    lab.refuse = [RunRefused("QC failed: mean depth 9.1 < 15")]
    failed = lab.refresh()
    assert failed.outcome == "failed"
    assert failed.error == "RunRefused"
    assert (failed.pre_run, failed.run) == (None, None)
    assert not [e for e in lab.events if e.startswith(("explain", "dashboard"))]
    lab.events.clear()
    lab.plan = {}
    second = lab.refresh()
    assert statuses(second)["clinvar"] == "unchanged"  # fetched by the refused refresh
    assert lab.events.count("run") == 1
    (rid,) = lab.runs
    assert (second.outcome, second.run, second.base_run) == ("complete", rid, first)
    assert [r.outcome for r in read(lab.lab)] == ["complete", "failed", "complete"]


def test_refused_pre_run_still_fetches_and_runs(lab: Lab) -> None:
    lab.snapshot(analyser="b" * 64)
    lab.refuse = [RunRefused("another run is under way")]
    record = lab.refresh()
    # A2: after the fetch the analyser still differs, so the run is tried again.
    assert lab.events[: 1 + len(REFRESHED)] == [
        "run",
        *(f"fetch:{name}" for name in REFRESHED),
    ]
    assert lab.events.count("run") == 2
    assert len(lab.runs) == 1  # the run after the fetch made a snapshot
    assert (record.outcome, record.error) == ("failed", "RunRefused")
    assert (record.pre_run, record.run) == (None, None)
    assert not [e for e in lab.events if e.startswith(("explain", "dashboard"))]


def test_second_refresh_while_locked_is_refused(lab: Lab) -> None:
    lab.snapshot()
    work = lab.paths.reports / "work"
    work.mkdir(parents=True)
    fd = os.open(work / ".refresh.lock", os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record = lab.refresh()
    finally:
        os.close(fd)
    assert lab.events == []  # nothing fetched, nothing run
    assert (record.outcome, record.error) == ("failed", "RefreshRefused")
    assert read(lab.lab) == [record]
    assert lab.refresh().outcome == "complete"  # the lock is free again


def test_without_explain_nothing_is_explained(lab: Lab) -> None:
    """Privacy audit I1: explaining sends findings to Anthropic, so a refresh does it only
    when asked; the dashboard is still written."""
    lab.snapshot(analyser=None)
    record = lab.refresh(explain=False)
    (pre,) = lab.runs
    assert (record.outcome, record.error) == ("complete", None)
    assert not [e for e in lab.events if e.startswith("explain")]
    assert lab.events[-1] == f"dashboard:{pre}"


def test_explain_failure_does_not_fail_the_refresh(lab: Lab) -> None:
    latest = lab.snapshot()
    lab.explain_error = RuntimeError("claude refused to explain BRCA1 chr17:2002")
    record = lab.refresh()
    assert record.outcome == "complete"
    assert record.error == "explain:RuntimeError"
    assert lab.events[-1] == f"dashboard:{latest}"
    assert "BRCA1" not in lab_text(lab.lab)


def test_error_text_stays_out_of_the_lab(lab: Lab) -> None:
    lab.snapshot()
    lab.plan = {"clinvar": "20261018", "gwas": "fail"}
    lab.refuse = [RunRefused("QC refused the sample: BRCA1 chr17:2002 depth 3")]
    record = lab.refresh()
    assert (record.outcome, record.error) == ("failed", "RunRefused")
    text = lab_text(lab.lab)
    assert "RunRefused" in text
    assert "BRCA1" not in text and "chr17" not in text and "depth" not in text
    errors = list((lab.paths.reports / "refresh").glob("*.error.txt"))
    assert [p.name for p in errors] == [
        record.started.strftime("%Y%m%dT%H%M%SZ") + ".error.txt"
    ]
    detail = errors[0].read_text()
    assert "QC refused the sample: BRCA1 chr17:2002 depth 3" in detail
    assert "gwas: bootstrap: curl failed for BRCA1 chr17:2002" in detail
    assert errors[0].stat().st_mode & 0o077 == 0


def test_an_unexpected_error_is_recorded_as_a_failed_refresh(
    lab: Lab, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab.snapshot()

    def unreadable(*args: Any) -> bool:
        raise OSError("cannot hash /lab/clinvar/BRCA1.vcf.gz")

    monkeypatch.setattr(refresh, "run_needed", unreadable)
    record = lab.refresh()
    assert (record.outcome, record.error) == ("failed", "OSError")
    assert read(lab.lab) == [record]
    assert "BRCA1" not in lab_text(lab.lab)
    assert "BRCA1" in refresh.error_path(lab.paths.reports, record.started).read_text()


def test_ctrl_c_still_leaves_a_log_record(
    lab: Lab, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab.snapshot()

    def interrupted(name: str) -> subprocess.CompletedProcess[str]:
        raise KeyboardInterrupt

    monkeypatch.setattr(lab, "fetch", interrupted)
    with pytest.raises(KeyboardInterrupt):
        lab.refresh()
    [record] = read(lab.lab)
    assert (record.outcome, record.error) == ("failed", "KeyboardInterrupt")


def test_a_torn_log_line_does_not_block_the_next_refresh(lab: Lab) -> None:
    lab.snapshot()
    first = lab.refresh()
    with (lab.lab / "refresh" / "log.jsonl").open("a") as log:
        log.write('{"started": "2026-10-')
    second = lab.refresh()
    assert second.outcome == "complete" and second.base_run == first.run
    assert read(lab.lab) == [first, second]  # the cut line does not swallow the next


def test_retain_keeps_what_the_last_two_refreshes_read(lab: Lab) -> None:
    lab_ = lab.lab
    old = lab.snapshot()  # reads clinvar FIRST
    put(lab_, "clinvar", "20261008")
    put(lab_, "litvar", "20261008")
    put(lab_, "medlineplus", "20261008")
    a = lab.snapshot()  # reads 20261008
    put(lab_, "clinvar", "20261015")
    b = lab.snapshot()  # reads 20261015
    put(lab_, "clinvar", "20261022")  # the newest, which no run read yet
    put(lab_, "litvar", "20261022")
    put(lab_, "medlineplus", "20261022")

    def done(rid: str, outcome: str = "complete") -> RefreshRecord:
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, outcome, (), None, None, rid)  # pyright: ignore[reportArgumentType]

    records = [done(old), done(a), done(b, "partial"), done(old, "failed")]
    deleted = refresh.retain(lab_, records, lab.root)
    assert sorted(p.relative_to(lab_).as_posix() for p in deleted) == [
        f"clinvar/{FIRST}",  # only `old` read it, and it is not among the last two
        f"litvar/{FIRST}",
        "litvar/20261008",  # LitVar keeps its newest only
        f"medlineplus/{FIRST}",
        "medlineplus/20261008",
    ]
    assert sorted(p.name for p in (lab_ / "clinvar").iterdir()) == [
        "20261008",
        "20261015",
        "20261022",
    ]
    # Every other source keeps its only folder, and CGD (pinned) is never touched.
    assert (lab_ / "cpic" / FIRST).is_dir() and (lab_ / "cgd" / "20250201").is_dir()


def test_retention_runs_after_a_complete_refresh(lab: Lab) -> None:
    lab.snapshot()
    lab.plan = {"clinvar": "20261018"}
    lab.refresh()
    lab.plan = {"clinvar": "20261025"}
    lab.refresh()
    # The last two refreshes read 20261018 and 20261025; FIRST went.
    assert sorted(p.name for p in (lab.lab / "clinvar").iterdir()) == [
        "20261018",
        "20261025",
    ]


# --- what reaches the network ---------------------------------------------------------------

_VAR = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*|[0-9])")
_VALUED = {
    *("-o", "-D", "-C", "--retry", "--retry-delay", "--range"),
    *("--connect-timeout", "--speed-limit", "--speed-time"),
}
# Variables a refreshed step's URL may hold besides the script's literal constants, each
# with what makes it the same for every user.
_UPSTREAM = {
    ("step_clinvar", "name"): "the file name ClinVar's own .md5 names",
    ("step_cpic", "entry"): "each entry of the constant $cpic_files",
    (
        "publication_step",
        "url",
    ): "its fourth argument, a constant (checked at each call)",
}


def _functions(script: str) -> dict[str, list[str]]:
    """Each function's body as logical lines (continuations joined, comments dropped)."""
    joined = re.sub(r"\\\n\s*", " ", script)
    bodies: dict[str, list[str]] = {}
    current: list[str] | None = None
    for line in joined.splitlines():
        if m := re.fullmatch(r"([a-z_0-9]+)\(\) \{", line):
            current = bodies.setdefault(m.group(1), [])
        elif line == "}":
            current = None
        elif current is not None and not line.strip().startswith("#"):
            current.append(line.strip())
    return bodies


def _constants(script: str) -> set[str]:
    """The top-level variables assigned a literal: no command substitution, and only other
    such constants expanded in the value."""
    found: set[str] = set()
    for line in script.splitlines():
        m = re.fullmatch(r"([a-z_0-9]+)=(.*)", line)
        if m is None or "$(" in m.group(2) or "`" in m.group(2):
            continue
        if all(v in found for v in _VAR.findall(m.group(2))):
            found.add(m.group(1))
    return found


def _urls(line: str) -> list[str]:
    """The URL arguments of a `curl`, `download` (curl's arguments), `fetch` or
    `publication_step` line."""
    if line.split(" ", 1)[0] not in ("curl", "download", "fetch", "publication_step"):
        return []
    words = shlex.split(line)
    if words[0] in ("curl", "download"):
        urls: list[str] = []
        skip = False
        for word in words[1:]:
            if skip:
                skip = False
            elif word in _VALUED:
                skip = True
            elif not word.startswith("-"):
                urls.append(word)
        return urls
    if words[0] == "fetch":
        return [words[1]]
    if words[0] == "publication_step":
        return [words[4]]
    return []


def url_violations(script: str) -> list[str]:
    """Each URL in a refreshed step (or a helper it calls) built from a variable that is
    neither a literal constant, a local built only from constants, nor a named upstream
    value."""
    constants = _constants(script)
    bodies = _functions(script)
    checked = [f"step_{name}" for name in REFRESHED] + ["publication_step"]
    bad: list[str] = []
    for function in checked:
        local = set(constants)
        for line in bodies[function]:
            if m := re.fullmatch(r"([a-z_0-9]+)=(\S*)", line):
                if "$(" not in m.group(2) and all(
                    v in local for v in _VAR.findall(m.group(2))
                ):
                    local.add(m.group(1))
                else:  # a constant's name given another value is no constant here
                    local.discard(m.group(1))
            for url in _urls(line):
                for var in _VAR.findall(url):
                    if var not in local and (function, var) not in _UPSTREAM:
                        bad.append(f"{function}: {var} in {url}")
    return bad


def _helper_urls(script: str) -> dict[str, list[str]]:
    bodies = _functions(script)
    return {
        f: [u for line in bodies[f] for u in _urls(line)] for f in ("download", "fetch")
    }


def test_no_request_depends_on_the_genome(lab: Lab) -> None:
    lab.snapshot()
    lab.refresh()
    assert lab.fetched == list(REFRESHED)  # names only, never anything of the user's
    script = BOOTSTRAP.read_text()
    assert url_violations(script) == []
    # `download` (the one curl) requests only the arguments it is given, and `fetch URL
    # DEST` its first argument only, which each call above passes.
    assert _helper_urls(script) == {"download": ["$@"], "fetch": ["$1"]}
    assert set(REFRESHED) <= {
        m.group(1) for m in re.finditer(r"^step_([a-z_]+)\(\)", script, re.MULTILINE)
    }


def test_a_url_built_from_another_variable_is_caught() -> None:
    script = BOOTSTRAP.read_text().replace(
        'publication_step litvar litvar "$litvar_file" "$litvar_url"',
        'curl -fL -o "$work/x" "https://www.ncbi.nlm.nih.gov/litvar/$rsid"\n'
        '  publication_step litvar litvar "$litvar_file" "$litvar_url"',
    )
    assert url_violations(script) == [
        "step_litvar: rsid in https://www.ncbi.nlm.nih.gov/litvar/$rsid"
    ]


# --- the command ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("outcome", "code"), [("complete", 0), ("partial", 3), ("failed", 1)]
)
def test_command_exit_code_follows_the_outcome(
    lab: Lab, monkeypatch: pytest.MonkeyPatch, outcome: str, code: int
) -> None:
    seen: dict[str, Any] = {}

    def fake(paths: Paths, **kw: Callable[..., Any]) -> RefreshRecord:
        seen["paths"] = paths
        seen["fetch"] = kw["fetch"]
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, outcome, (), None, None, None)  # pyright: ignore[reportArgumentType]

    monkeypatch.setattr(refresh, "refresh", fake)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    calls: list[list[str]] = []

    def fake_run(argv: list[str], **kw: Any) -> subprocess.CompletedProcess[str]:
        calls.append(argv)
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(refresh.subprocess, "run", fake_run)
    assert cli.main(["refresh"]) == code
    seen["fetch"]("clinvar")
    script = lab.paths.repo / "scripts" / "bootstrap.sh"
    assert calls == [["bash", str(script), "--refresh", "clinvar"]]
    assert seen["paths"] == lab.paths


@pytest.mark.parametrize(("argv", "asked"), [([], False), (["--explain"], True)])
def test_command_explains_only_with_the_flag(
    lab: Lab, monkeypatch: pytest.MonkeyPatch, argv: list[str], asked: bool
) -> None:
    seen: dict[str, Any] = {}

    def fake(paths: Paths, **kw: Callable[..., Any]) -> RefreshRecord:
        seen.update(kw)
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, "failed", (), None, None, None)

    monkeypatch.setattr(refresh, "refresh", fake)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    cli.main(["refresh", *argv])
    assert (seen["explain"] is not None) is asked


def test_a_stalled_fetch_is_killed_after_the_timeout(
    lab: Lab, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: dict[str, Any] = {}

    def fake(paths: Paths, **kw: Callable[..., Any]) -> RefreshRecord:
        seen.update(kw)
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, "failed", (), None, None, None)

    monkeypatch.setattr(refresh, "refresh", fake)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    script = lab.paths.repo / "scripts" / "bootstrap.sh"
    script.parent.mkdir()
    script.write_text("sleep 30\n")
    monkeypatch.setattr(refresh, "FETCH_TIMEOUT", 0.2)
    cli.main(["refresh"])
    with pytest.raises(subprocess.TimeoutExpired):
        seen["fetch"]("clinvar")


def test_a_failed_explain_call_fails_the_explain_step(
    lab: Lab, monkeypatch: pytest.MonkeyPatch
) -> None:
    # An expired login fails every call; the step must raise for `_steps` to record it,
    # as `genokeeper explain` exits 1.
    seen: dict[str, Any] = {}

    def fake(paths: Paths, **kw: Callable[..., Any]) -> RefreshRecord:
        seen.update(kw)
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, "complete", (), None, None, None)

    @contextlib.contextmanager
    def store(_: Path) -> Iterator[None]:
        yield None

    monkeypatch.setattr(refresh, "refresh", fake)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    monkeypatch.setattr(refresh, "open_store", store)
    # The step asks for what `genokeeper explain` asks for (`to_explain`: the findings and
    # the stage-3 watch items, BR-025).
    runs: list[str] = []
    monkeypatch.setattr(refresh, "to_explain", lambda _con, run: runs.append(run) or [])
    monkeypatch.setattr(refresh, "newest_medlineplus", lambda _: object())
    langs: list[tuple[str, str, str]] = []

    def explain(*_: Any, lang: str, model: str, tone: str) -> dict[str, int]:
        langs.append((lang, model, tone))
        return {"cached": 0, "written": 0, "failed": 2}

    monkeypatch.setattr(refresh, "explain", explain)
    # The language, model and tone are `genokeeper explain`'s: the environment's, else
    # the defaults; the model and tone as the command started, checked before the run.
    monkeypatch.setenv("GENOKEEPER_TONE", "engaging")
    cli.main(["refresh", "--explain"])
    monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    with pytest.raises(ExplainError, match="2 explain calls failed"):
        seen["explain"]("run-1")
    monkeypatch.setenv("GENOKEEPER_LANG", "ru")
    monkeypatch.delenv("GENOKEEPER_TONE")
    with pytest.raises(ExplainError):
        seen["explain"]("run-1")
    assert langs == [("en", MODEL, "engaging"), ("ru", MODEL, "engaging")]
    assert runs == ["run-1", "run-1"]


def test_command_announces_each_fetch_and_a_refused_run(
    lab: Lab, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """A fetch can take minutes and a refused run is otherwise silent until the end."""
    seen: dict[str, Any] = {}

    def fake(paths: Paths, **kw: Callable[..., Any]) -> RefreshRecord:
        seen.update(kw)
        at = datetime(2026, 10, 1, tzinfo=UTC)
        return RefreshRecord(at, at, "failed", (), None, None, None)

    monkeypatch.setattr(refresh, "refresh", fake)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    monkeypatch.setattr(
        refresh.subprocess,
        "run",
        lambda argv, **_: subprocess.CompletedProcess(argv, 0, "", ""),  # pyright: ignore[reportUnknownLambdaType, reportUnknownArgumentType]
    )
    cli.main(["refresh"])
    capsys.readouterr()
    seen["fetch"]("litvar")
    assert "fetching litvar" in capsys.readouterr().err
    with pytest.raises(Exception) as refused:  # the lab holds no FASTQ pair
        seen["run"]()
    err = capsys.readouterr().err
    assert "grading the calls" in err
    assert f"the run was refused: {type(refused.value).__name__}" in err


def test_an_unknown_tone_refuses_the_refresh_before_it_runs(
    lab: Lab, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def never(*_: Any, **__: Any) -> RefreshRecord:
        raise AssertionError("the refresh ran")

    monkeypatch.setattr(refresh, "refresh", never)
    monkeypatch.setattr(refresh, "load_paths", lambda: lab.paths)
    monkeypatch.setenv("GENOKEEPER_TONE", "lively")
    assert cli.main(["refresh"]) == 2
    assert "unknown tone 'lively'" in capsys.readouterr().err
