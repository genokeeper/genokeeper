"""Each environment variable is read under its name, else under its name before the rename
from genomelab (`genokeeper.env`)."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from genokeeper.config import load_paths
from genokeeper.env import getenv
from genokeeper.lang import language


def test_the_new_name_wins_and_an_empty_value_is_unset() -> None:
    both = {"GENOKEEPER_X": "new", "GENOMELAB_X": "old"}
    assert getenv("X", both) == "new"
    assert getenv("X", {"GENOKEEPER_X": "", "GENOMELAB_X": "old"}) == "old"
    assert getenv("X", {"GENOMELAB_X": ""}) is None
    assert getenv("X", {}) is None


def test_home_falls_back(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("GENOKEEPER_HOME", raising=False)
    monkeypatch.setenv("GENOMELAB_HOME", str(tmp_path))
    assert load_paths().lab == tmp_path


def test_lang_falls_back(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GENOKEEPER_LANG", raising=False)
    monkeypatch.setenv("GENOMELAB_LANG", "ru")
    assert language() == "ru"


def test_gate_falls_back() -> None:
    """`tools.GATE` is read at import, so in a fresh interpreter."""
    env = {k: v for k, v in os.environ.items() if k != "GENOKEEPER_GATE"}
    done = subprocess.run(
        [sys.executable, "-c", "from genokeeper import tools; print(tools.GATE)"],
        env=env | {"GENOMELAB_GATE": "/old-gate"},
        capture_output=True,
        text=True,
        check=True,
    )
    assert done.stdout == "/old-gate\n"
