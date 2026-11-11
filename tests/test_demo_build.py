"""`python -m genokeeper.demo build`: two runs of one pair, a refresh record, and its
refusals before anything starts."""

import functools
import hashlib
import http.server
import threading
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from genokeeper import config, refresh_log
from genokeeper.calibrate import hg002
from genokeeper.demo import build as demo_build
from genokeeper.demo import demo_paths
from genokeeper.demo.build import BASE_CLINVAR, BuildRefused, build
from genokeeper.pipeline import RunRefused

BODY = b"##fileDate=2026-10-11\nbase clinvar bytes\n"


class _Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: Any) -> None:
        pass


def _publish(folder: Path, body: bytes, md5: str | None = None) -> None:
    name = f"clinvar_{BASE_CLINVAR}.vcf.gz"
    (folder / name).write_bytes(body)
    (folder / f"{name}.tbi").write_bytes(b"tbi")
    digest = md5 or hashlib.md5(body).hexdigest()
    (folder / f"{name}.md5").write_text(f"{digest}  {name}\n")


@pytest.fixture
def ncbi(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    site = tmp_path / "ncbi"
    site.mkdir()
    _publish(site, BODY)
    server = http.server.ThreadingHTTPServer(
        ("127.0.0.1", 0), functools.partial(_Quiet, directory=str(site))
    )
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setattr(
        demo_build, "CLINVAR_BASE_URL", f"http://127.0.0.1:{server.server_port}/"
    )
    yield site
    server.shutdown()
    server.server_close()


class Calls:
    def __init__(self) -> None:
        self.kwargs: list[dict[str, Any]] = []

    def __call__(self, paths: Any, **kwargs: Any) -> str:
        self.kwargs.append({"paths": paths, **kwargs})
        return f"run-{len(self.kwargs)}"


@pytest.fixture
def lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, ncbi: Path) -> Path:
    lab = tmp_path / "lab"
    for p in hg002.downsampled_pair(lab, "pcrfree"):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"reads")
    newest = lab / "clinvar" / "20261108" / "clinvar_GRCh38.vcf.gz"
    newest.parent.mkdir(parents=True)
    newest.write_bytes(b"new")
    monkeypatch.setattr(demo_build, "need_space", lambda *a: None)
    monkeypatch.setattr(demo_build.pipeline, "open_faf_table", lambda lab: "faf")
    return lab


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> Calls:
    calls = Calls()
    monkeypatch.setattr(demo_build.pipeline, "run_genome", calls)
    return calls


def test_build_runs_twice_on_hg002(lab: Path, calls: Calls) -> None:
    assert build(lab, threads=3) == ("run-1", "run-2")
    pair = hg002.downsampled_pair(lab, "pcrfree")
    first, second = calls.kwargs
    for call in calls.kwargs:
        assert call["paths"] == demo_paths(lab)
        assert (call["fastq1"], call["fastq2"]) == pair
        assert call["delivered"] is None
        assert call["threads"] == 3
    assert first["clinvar38"] == demo_build.base_clinvar_path(demo_paths(lab).repo)
    assert first["clinvar38"].read_bytes() == BODY
    assert second["clinvar38"] == lab / "clinvar" / "20261108" / "clinvar_GRCh38.vcf.gz"


def test_build_writes_refresh_record_in_demo_lab(lab: Path, calls: Calls) -> None:
    build(lab, threads=None)
    (record,) = refresh_log.read(demo_paths(lab).repo / "lab")
    assert (record.base_run, record.run, record.pre_run) == ("run-1", "run-2", None)
    clinvar = next(s for s in record.sources if s.name == "clinvar")
    assert (clinvar.before, clinvar.after) == (BASE_CLINVAR, "20261108")
    others = [s for s in record.sources if s is not clinvar]
    assert others and all(s.status.value == "unchanged" for s in others)
    assert not (lab / "refresh").exists()


def test_build_never_touches_repo_reports(lab: Path, calls: Calls) -> None:
    reports = demo_paths(lab).reports
    assert config._REPO not in (reports, *reports.parents)


def test_build_refuses_missing_pair(lab: Path, calls: Calls) -> None:
    hg002.downsampled_pair(lab, "pcrfree")[1].unlink()
    with pytest.raises(BuildRefused, match="pair is missing"):
        build(lab, threads=None)
    assert calls.kwargs == []
    assert not (demo_paths(lab).repo / "clinvar").exists()


def test_build_refuses_low_disk(
    lab: Path, calls: Calls, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[int] = []

    def short(folder: Path, need: int, what: str) -> None:
        seen.append(need)
        raise RunRefused("needs 150 GB free")

    monkeypatch.setattr(demo_build, "need_space", short)
    with pytest.raises(BuildRefused, match="150 GB"):
        build(lab, threads=None)
    assert seen == [130 * 10**9]
    assert calls.kwargs == []
    assert not (demo_paths(lab).repo / "clinvar").exists()


def test_base_clinvar_md5_refused(lab: Path, calls: Calls, ncbi: Path) -> None:
    _publish(ncbi, BODY, md5="0" * 32)
    with pytest.raises(BuildRefused, match="md5"):
        build(lab, threads=None)
    assert not demo_build.base_clinvar_path(demo_paths(lab).repo).exists()
    assert calls.kwargs == []


def test_build_reuses_base_clinvar(
    lab: Path, calls: Calls, monkeypatch: pytest.MonkeyPatch
) -> None:
    build(lab, threads=None)

    def refetch(*args: Any) -> str:
        pytest.fail("fetched again")

    monkeypatch.setattr(demo_build, "_fetch", refetch)
    build(lab, threads=None)
    assert len(calls.kwargs) == 4
