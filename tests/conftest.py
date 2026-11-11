"""Shared test setup. Every VCF, BAM, CRAM and FASTA a test reads is built by `fixtures` at
test time, in the test's tmp_path; none is committed, except the public HG002 tools' outputs
under `fixtures/hg002/` (ExpansionHunter's VCF is there as `.vcf.txt`)."""

import os
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from genokeeper import env, tools


@pytest.fixture(autouse=True)
def no_old_variables(monkeypatch: pytest.MonkeyPatch) -> None:
    """A shell that still sets a variable under its name before the rename
    (`GENOMELAB_HOME`, ...) reaches no test: each falls back to it (`genokeeper.env`)."""
    for name in [n for n in os.environ if n.startswith(env.OLD_PREFIX)]:
        monkeypatch.delenv(name)


@pytest.fixture(autouse=True)
def ample_free_disk(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every test sees a disk with 1 TB free: `shutil.disk_usage` is the probe
    `calibrate.records.need_space` reads, and a test that downloads or re-aligns must not pass
    or fail with the free space of the machine running it. A test of the refusal stubs a
    low value itself."""
    usage = SimpleNamespace(total=2 * 10**12, used=10**12, free=10**12)
    monkeypatch.setattr(shutil, "disk_usage", lambda _path: usage)


@pytest.fixture(autouse=True)
def no_machine_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every subprocess a test starts inherits a gate path that does not exist, and the
    tests' own process uses it too (`genokeeper.tools.GATE`, read at import, before this
    runs), so nothing queues behind the machine's core budget."""
    monkeypatch.setenv("GENOKEEPER_GATE", "/nonexistent-gate")
    monkeypatch.setattr(tools, "GATE", Path("/nonexistent-gate"))


@pytest.fixture(scope="session")
def wheel(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The package's wheel, built once per session (per worker under xdist)."""
    out = tmp_path_factory.mktemp("wheel")
    repo = Path(__file__).resolve().parents[1]
    command = ["uv", "build", "--wheel", "--out-dir", str(out)]
    subprocess.run(command, cwd=repo, check=True, capture_output=True)
    (built,) = out.glob("*.whl")
    return built
