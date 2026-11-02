"""Shared test setup. Every VCF, BAM, CRAM and FASTA a test reads is built by `fixtures` at
test time, in the test's tmp_path; none is committed, except the public HG002 tools' outputs
under `fixtures/hg002/` (ExpansionHunter's VCF is there as `.vcf.txt`)."""

import os
from pathlib import Path

import pytest

from genokeeper import env, tools


@pytest.fixture(autouse=True)
def no_old_variables(monkeypatch: pytest.MonkeyPatch) -> None:
    """A shell that still sets a variable under its name before the rename
    (`GENOMELAB_HOME`, ...) reaches no test: each falls back to it (`genokeeper.env`)."""
    for name in [n for n in os.environ if n.startswith(env.OLD_PREFIX)]:
        monkeypatch.delenv(name)


@pytest.fixture(autouse=True)
def no_machine_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every subprocess a test starts inherits a gate path that does not exist, and the
    tests' own process uses it too (`genokeeper.tools.GATE`, read at import, before this
    runs), so nothing queues behind the machine's core budget."""
    monkeypatch.setenv("GENOKEEPER_GATE", "/nonexistent-gate")
    monkeypatch.setattr(tools, "GATE", Path("/nonexistent-gate"))
