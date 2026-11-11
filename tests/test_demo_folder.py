"""The genome-ask skill's commands answer in a demo folder from the installed wheel: a real
"outside a checkout" (no `config/grading.toml`, the repository is the working folder), where
the lab is the folder's own whatever `$GENOKEEPER_HOME` says."""

import json
import os
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from fixtures.demo_site import GENE, RUN1, RUN2, packaged, serve

from genokeeper import refresh_log
from genokeeper.demo.install import install
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus

REPO = Path(__file__).resolve().parents[1]
# The planted lab's refresh date: it must reach no output.
PLANTED = datetime(2031, 7, 17, 9, tzinfo=UTC)


def _installed(tmp_path: Path) -> Path:
    """The `genokeeper` of the built wheel in an isolated environment."""
    wheels = tmp_path / "wheel"
    subprocess.run(
        ["uv", "build", "--wheel", "--out-dir", str(wheels)],
        cwd=REPO,
        check=True,
        capture_output=True,
    )
    (wheel,) = wheels.glob("*.whl")
    env = tmp_path / "env"
    for command in (
        ["uv", "venv", "--python", sys.executable, str(env)],
        ["uv", "pip", "install", "--python", str(env / "bin" / "python"), str(wheel)],
    ):
        subprocess.run(command, check=True, capture_output=True)
    return env / "bin" / "genokeeper"


def test_skill_commands_in_demo_folder(tmp_path: Path) -> None:
    genokeeper = _installed(tmp_path)
    with serve(packaged(tmp_path / "site")) as base:
        folder = install(tmp_path / "d", base=base)
    planted = tmp_path / "home-lab"
    clinvar = SourceOutcome("clinvar", "20261011", "20261108", SourceStatus.NEW)
    record = RefreshRecord(PLANTED, PLANTED, "complete", (clinvar,), RUN1, None, RUN2)
    refresh_log.append(planted, record)
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV")
    }
    env["GENOKEEPER_HOME"] = str(planted)

    def run(*args: str) -> Any:
        done = subprocess.run(
            [str(genokeeper), *args],
            cwd=folder,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        assert done.returncode == 0, done.stderr
        assert PLANTED.date().isoformat() not in done.stdout + done.stderr
        return json.loads(done.stdout)

    found = run("query", "--json")
    assert found["run_id"] == RUN2
    first = found["findings"][0]["id"]
    assert run("show", first, "--json")["finding"]["id"] == first
    assert run("genotype", "--gene", GENE, "--json")["sites"]
    changes = run("changes", "--json")
    assert changes["run"] == RUN2 and changes["refresh"] is not None
    assert "check" not in changes
