"""The analyser's identity (`genokeeper.identity`): a hash of the package's code and
resources and of config/grading.toml, so a change to any of them reads as a different
analyser."""

import shutil
from collections.abc import Callable, Iterator
from contextlib import AbstractContextManager, contextmanager
from pathlib import Path

import pytest

from genokeeper import identity
from genokeeper.identity import analyser_identity

type Change = Callable[[Path, pytest.MonkeyPatch], AbstractContextManager[None]]


def _package_copy(tmp_path: Path) -> Path:
    """A copy of the package under `tmp_path`, which the identity then reads."""
    copy = tmp_path / "genokeeper"
    shutil.copytree(identity.PACKAGE_ROOT, copy, ignore=shutil.ignore_patterns("*.pyc"))
    return copy


def _append_byte(path: Path) -> None:
    path.write_bytes(path.read_bytes() + b"\n")


@contextmanager
def _edited_grading(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    grading = tmp_path / "grading.toml"
    shutil.copyfile(identity.GRADING_FILE, grading)
    _append_byte(grading)
    monkeypatch.setattr(identity, "GRADING_FILE", grading)
    yield


@contextmanager
def _edited_module(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    copy = _package_copy(tmp_path)
    _append_byte(copy / "identity.py")
    monkeypatch.setattr(identity, "PACKAGE_ROOT", copy)
    yield


@contextmanager
def _edited_resource(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    copy = _package_copy(tmp_path)
    resource = next(p for p in sorted((copy / "resources").rglob("*")) if p.is_file())
    _append_byte(resource)
    monkeypatch.setattr(identity, "PACKAGE_ROOT", copy)
    yield


@pytest.mark.parametrize("change", [_edited_grading, _edited_module, _edited_resource])
def test_identity_is_stable_and_follows_code_resources_and_grading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: Change
) -> None:
    first = analyser_identity()
    assert first == analyser_identity() and len(first) == 64
    with change(tmp_path, monkeypatch):
        assert analyser_identity() != first


def test_identity_ignores_bytecode_and_reads_a_copy_alike(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An unchanged copy of the package has the identity of the original, with or without
    compiled files beside it."""
    first = analyser_identity()
    copy = _package_copy(tmp_path)
    (copy / "__pycache__").mkdir(exist_ok=True)
    (copy / "__pycache__" / "x.cpython-314.pyc").write_bytes(b"junk")
    (copy / "stray.pyc").write_bytes(b"junk")
    monkeypatch.setattr(identity, "PACKAGE_ROOT", copy)
    assert analyser_identity() == first


def test_skill_edit_leaves_analyser_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The packaged skill is documentation for the reader, not part of the analyser."""
    first = analyser_identity()
    copy = _package_copy(tmp_path)
    _append_byte(copy / "resources" / "skills" / "genome-ask" / "SKILL.md")
    monkeypatch.setattr(identity, "PACKAGE_ROOT", copy)
    assert analyser_identity() == first


def test_demo_site_edit_leaves_analyser_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The demo site's pages and its recorded answer present the findings; they are not
    part of the analyser."""
    first = analyser_identity()
    copy = _package_copy(tmp_path)
    _append_byte(copy / "resources" / "demo" / "index.html")
    (copy / "resources" / "demo" / "transcript.md").write_text("---\n---\nx\n")
    monkeypatch.setattr(identity, "PACKAGE_ROOT", copy)
    assert analyser_identity() == first
