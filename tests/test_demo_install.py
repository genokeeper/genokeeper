"""`genokeeper demo`: the HG002 store downloaded from the demo's release, checked against its
manifest, and written into a folder Claude Code can work in."""

import html
import json
import shlex
import shutil
import tomllib
from pathlib import Path

import pytest
from fixtures.demo_site import RUN1, RUN2, packaged, serve

import genokeeper
from genokeeper import cli, refresh_log
from genokeeper.config import DEMO_MARKER, skill_resource
from genokeeper.demo import MANIFEST, REFRESH_LOG_FILE, STORE_FILE
from genokeeper.demo import install as inst
from genokeeper.demo.frame import FRAME, FRAME_SUMMARY
from genokeeper.demo.install import QUESTIONS, DemoRefused, install
from genokeeper.store.db import SCHEMA_VERSION


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return packaged(tmp_path_factory.mktemp("site"))


def _copy(site: Path, tmp_path: Path) -> Path:
    """A copy of the site a test may tamper with."""
    return Path(shutil.copytree(site, tmp_path / "site"))


def _files(folder: Path) -> list[str]:
    return sorted(str(p.relative_to(folder)) for p in folder.rglob("*"))


def test_install_writes_folder(site: Path, tmp_path: Path) -> None:
    with serve(site) as base:
        folder = install(tmp_path / "d", base=base)
    assert folder == tmp_path / "d"
    assert json.loads((folder / DEMO_MARKER).read_text()) == json.loads(
        (site / MANIFEST).read_text()
    )
    reports = folder / "reports"
    assert (reports / STORE_FILE).read_bytes() == (site / STORE_FILE).read_bytes()
    assert sorted(p.name for p in (reports / "snapshots").iterdir()) == [RUN1, RUN2]
    assert (folder / "lab" / DEMO_MARKER).is_file()
    log = folder / "lab" / "refresh" / refresh_log.LOG
    assert log.read_bytes() == (site / REFRESH_LOG_FILE).read_bytes()
    page = (folder / "dashboard.html").read_text(encoding="utf-8")
    assert html.escape(FRAME_SUMMARY) in page and html.escape(FRAME) in page
    project = tomllib.loads((folder / "pyproject.toml").read_text())
    assert project["project"]["dependencies"] == [
        f"genokeeper=={genokeeper.__version__}"
    ]
    assert project["project"]["requires-python"] == ">=3.14"
    assert "build-system" not in project
    skill = folder / ".claude" / "skills" / "genome-ask" / "SKILL.md"
    assert skill.read_bytes() == skill_resource().read_bytes()
    claude = (folder / "CLAUDE.md").read_text(encoding="utf-8")
    assert FRAME in claude and "explanations are not generated here" in claude
    assert not list(tmp_path.rglob("*.part"))


def test_install_refuses_sha_mismatch(site: Path, tmp_path: Path) -> None:
    served = _copy(site, tmp_path)
    with (served / STORE_FILE).open("ab") as f:
        f.write(b"\0")
    with serve(served) as base, pytest.raises(DemoRefused, match="sha256"):
        install(tmp_path / "d", base=base)
    assert not (tmp_path / "d" / "reports" / STORE_FILE).exists()
    assert not list(tmp_path.rglob("*.part"))


def test_install_refuses_schema_mismatch(site: Path, tmp_path: Path) -> None:
    served = _copy(site, tmp_path)
    manifest = json.loads((served / MANIFEST).read_text())
    manifest["schema_version"] = SCHEMA_VERSION + 1
    (served / MANIFEST).write_text(json.dumps(manifest))
    with serve(served) as base, pytest.raises(DemoRefused, match="upgrade genokeeper"):
        install(tmp_path / "d", base=base)
    assert not (tmp_path / "d").exists()


@pytest.mark.parametrize("marker", [None, '{"name": "a project of its own"}'])
def test_install_refuses_foreign_folder(
    site: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    marker: str | None,
) -> None:
    folder = tmp_path / "d"
    folder.mkdir()
    (folder / "notes.txt").write_text("mine")
    if marker is not None:  # a `demo.json` that is not a demo manifest
        (folder / DEMO_MARKER).write_text(marker)
    with serve(site) as base:
        with pytest.raises(DemoRefused, match="not a demo folder"):
            install(folder, base=base)
        monkeypatch.setattr(inst, "DEMO_BASE", base)
        assert cli.main(["demo", "--dir", str(folder)]) == 2
    assert _files(folder) == sorted(["notes.txt", *([DEMO_MARKER] if marker else [])])
    assert (folder / "notes.txt").read_text() == "mine"
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1 and "not a demo folder" in err[0]


def test_install_updates_demo_folder(site: Path, tmp_path: Path) -> None:
    folder = tmp_path / "d"
    with serve(site) as base:
        install(folder, base=base)
        (folder / ".venv").mkdir()
        (folder / ".venv" / "pyvenv.cfg").write_text("kept")
        (folder / "notes.md").write_text("mine")
        install(folder, base=base)
    assert (folder / ".venv" / "pyvenv.cfg").read_text() == "kept"
    assert (folder / "notes.md").read_text() == "mine"
    assert (folder / "reports" / STORE_FILE).read_bytes() == (
        site / STORE_FILE
    ).read_bytes()


def test_install_cut_download_leaves_no_store(site: Path, tmp_path: Path) -> None:
    folder = tmp_path / "d"
    with (
        serve(site, cut_store_at=4096) as base,
        pytest.raises(DemoRefused, match="sha256"),
    ):
        install(folder, base=base)
    assert not (folder / "reports" / STORE_FILE).exists()
    assert not list(tmp_path.rglob("*.part"))
    with serve(site) as base:
        install(folder, base=base)
    assert (folder / "reports" / STORE_FILE).read_bytes() == (
        site / STORE_FILE
    ).read_bytes()


def test_install_prints_next_steps(
    site: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The dashboard's path, `cd <folder> && claude` and three questions, last on stdout;
    the size line on stderr, no progress line without a terminal, no bootstrap to run."""
    folder = tmp_path / "d"
    with serve(site) as base:
        monkeypatch.setattr(inst, "DEMO_BASE", base)
        assert cli.main(["demo", "--dir", str(folder)]) == 0
    out, err = capsys.readouterr()
    lines = out.splitlines()
    assert [line.strip() for line in lines[-3:]] == [
        f"{n}. {q}" for n, q in enumerate(QUESTIONS, 1)
    ]
    cd = lines.index(f"  cd {shlex.quote(str(folder))} && claude")
    assert any(str(folder / "dashboard.html") in line for line in lines[:cd])
    claude = (folder / "CLAUDE.md").read_text(encoding="utf-8")
    assert FRAME_SUMMARY in claude and all(q in claude for q in QUESTIONS)
    size = (site / STORE_FILE).stat().st_size
    assert f"{size / 1e6:.1f} MB" in err
    assert "\r" not in err and "bootstrap" not in err
