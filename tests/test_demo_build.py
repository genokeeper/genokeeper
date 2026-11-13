"""`python -m genokeeper.demo build`: two runs of one pair, a refresh record, and its
refusals before anything starts."""

import hashlib
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from fixtures.demo_site import serve

from genokeeper import refresh_log
from genokeeper.calibrate import hg002
from genokeeper.demo import STORE_FILE, demo_paths
from genokeeper.demo import build as demo_build
from genokeeper.demo.build import BASE_CLINVAR, BuildRefused, build
from genokeeper.pipeline import RunRefused

BODY = b"##fileDate=2026-05-03\nbase clinvar bytes\n"


def _publish(folder: Path, md5: str | None = None) -> None:
    name = f"clinvar_{BASE_CLINVAR}.vcf.gz"
    (folder / f"{name}.tbi").write_bytes(b"tbi")
    digest = md5 or hashlib.md5(BODY).hexdigest()
    (folder / f"{name}.md5").write_text(f"{digest}  {name}\n")


def _download(url: str, md5: str, out: Path) -> Path:
    """`hg002.download` without curl: BODY, refused (as it does) on another md5."""
    if hashlib.md5(BODY).hexdigest() != md5:
        raise ValueError(f"{url}: md5 mismatch")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(BODY)
    return out


@pytest.fixture
def ncbi(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    site = tmp_path / "ncbi" / BASE_CLINVAR[:4]
    site.mkdir(parents=True)
    _publish(site)
    monkeypatch.setattr(demo_build.hg002, "download", _download)
    with serve(site.parent) as base:
        monkeypatch.setattr(demo_build, "CLINVAR_BASE_URL", base)
        yield site


class Calls:
    def __init__(self) -> None:
        self.kwargs: list[dict[str, Any]] = []

    def __call__(self, paths: Any, **kwargs: Any) -> str:
        self.kwargs.append({"paths": paths, **kwargs})
        return f"run-{len(self.kwargs)}"


def _newest(lab: Path, file_date: str) -> Path:
    """The lab's newest ClinVar, in the folder 20261108, whose `##fileDate` is `file_date`."""
    newest = lab / "clinvar" / "20261108" / "clinvar_GRCh38.vcf.gz"
    newest.parent.mkdir(parents=True, exist_ok=True)
    newest.write_text(f"##fileDate={file_date}\n")
    return newest


@pytest.fixture
def lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, ncbi: Path) -> Path:
    lab = tmp_path / "lab"
    for p in hg002.downsampled_pair(lab, "pcrfree"):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"reads")
    _newest(lab, "2026-11-07")
    monkeypatch.setattr(demo_build, "need_space", lambda *a: None)
    monkeypatch.setattr(demo_build.pipeline, "open_faf_table", lambda lab: "faf")
    # The release as the run reads it: `##fileDate`, not the folder's name.
    monkeypatch.setattr(
        demo_build.pipeline,
        "_clinvar_release",
        lambda vcf: vcf.read_text().split("=")[1][:10].replace("-", ""),
    )
    return lab


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> Calls:
    calls = Calls()
    monkeypatch.setattr(demo_build.pipeline, "run_genome", calls)
    return calls


def _refused_before_anything(lab: Path, calls: Calls, match: str) -> None:
    with pytest.raises(BuildRefused, match=match):
        build(lab, threads=None)
    assert calls.kwargs == []
    assert not (demo_paths(lab).repo / "clinvar").exists()


def test_build_runs_twice_on_hg002(lab: Path, calls: Calls) -> None:
    assert build(lab, threads=3) == ("run-1", "run-2")
    pair = hg002.downsampled_pair(lab, "pcrfree")
    first, second = calls.kwargs
    for call in calls.kwargs:
        assert call["paths"] == demo_paths(lab)
        assert (call["fastq1"], call["fastq2"]) == pair
        assert call["delivered"] is None
        assert call["threads"] == 3
        assert call["skip_read_names"] is True
    assert first["clinvar38"] == demo_build.base_clinvar_path(demo_paths(lab).repo)
    assert first["clinvar38"].read_bytes() == BODY
    assert second["clinvar38"] == lab / "clinvar" / "20261108" / "clinvar_GRCh38.vcf.gz"


def test_build_writes_refresh_record_in_demo_lab(lab: Path, calls: Calls) -> None:
    build(lab, threads=None)
    (record,) = refresh_log.read(demo_paths(lab).repo / "lab")
    assert (record.base_run, record.run, record.pre_run) == ("run-1", "run-2", None)
    clinvar = next(s for s in record.sources if s.name == "clinvar")
    assert (clinvar.before, clinvar.after) == (BASE_CLINVAR, "20261107")
    others = [s for s in record.sources if s is not clinvar]
    assert others and all(s.status.value == "unchanged" for s in others)
    assert not (lab / "refresh").exists()


def test_build_refuses_missing_pair(lab: Path, calls: Calls) -> None:
    hg002.downsampled_pair(lab, "pcrfree")[1].unlink()
    _refused_before_anything(lab, calls, "pair is missing")


def test_build_refuses_existing_store(lab: Path, calls: Calls) -> None:
    store = demo_paths(lab).reports / STORE_FILE
    store.parent.mkdir(parents=True)
    store.write_bytes(b"runs 1 and 2")
    _refused_before_anything(lab, calls, "exists: to build again, delete it")


def test_build_refuses_a_lab_clinvar_not_newer(lab: Path, calls: Calls) -> None:
    _newest(lab, "2026-05-03")
    _refused_before_anything(lab, calls, f"not newer than {BASE_CLINVAR}")


def test_build_refuses_low_disk(
    lab: Path, calls: Calls, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[int] = []

    def short(folder: Path, need: int, what: str) -> None:
        seen.append(need)
        raise RunRefused("needs 150 GB free")

    monkeypatch.setattr(demo_build, "need_space", short)
    _refused_before_anything(lab, calls, "150 GB")
    assert seen == [130 * 10**9]


def test_resumed_build_asks_only_for_what_its_root_lacks(
    lab: Path, calls: Calls, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[int] = []

    def short(folder: Path, need: int, what: str) -> None:
        seen.append(need)
        raise RunRefused("short")

    monkeypatch.setattr(demo_build, "need_space", short)
    cache = demo_paths(lab).repo / "work" / "cache.bin"
    cache.parent.mkdir(parents=True)
    cache.write_bytes(b"x" * 10_000)
    os.link(cache, cache.parent / "hard")  # one inode, counted once
    (cache.parent / "soft").symlink_to(cache)  # not followed
    held = cache.stat().st_blocks * 512
    for need in (10**6, held - 1):
        monkeypatch.setattr(demo_build, "NEED_BYTES", need)
        with pytest.raises(BuildRefused, match="short"):
            build(lab, threads=None)
    assert seen == [10**6 - held, 0]


def test_base_clinvar_md5_refused(lab: Path, calls: Calls, ncbi: Path) -> None:
    _publish(ncbi, md5="0" * 32)
    with pytest.raises(BuildRefused, match="md5"):
        build(lab, threads=None)
    assert not demo_build.base_clinvar_path(demo_paths(lab).repo).exists()
    assert calls.kwargs == []


def test_base_clinvar_fetched_from_its_years_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[str] = []

    def offline(url: str, timeout: int) -> None:
        seen.append(url)
        raise OSError("offline")

    monkeypatch.setattr(demo_build, "BASE_CLINVAR", "20250105")
    monkeypatch.setattr(demo_build.urllib.request, "urlopen", offline)
    with pytest.raises(BuildRefused, match="offline"):
        demo_build.fetch_base_clinvar(tmp_path)
    archive = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/"
    assert seen == [f"{archive}2025/clinvar_20250105.vcf.gz.md5"]
