import os
import subprocess
import sys
from pathlib import Path

import pytest

from genokeeper import tools
from genokeeper.tools import ToolFailed, ToolMissing, gated, run_tool


def test_failure_carries_the_tools_stderr() -> None:
    script = "import sys; sys.stderr.write('sequence chr2 not found'); sys.exit(3)"
    with pytest.raises(ToolFailed, match="sequence chr2 not found") as caught:
        run_tool([sys.executable, "-c", script])
    assert isinstance(caught.value, subprocess.CalledProcessError)
    assert caught.value.returncode == 3


def test_env_is_added_to_the_tools_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("GENOKEEPER_TOOLS_KEPT", "kept")
    script = (
        "import os; print(os.environ['GENOKEEPER_TOOLS_KEPT'], os.environ['ADDED'])"
    )
    out = run_tool([sys.executable, "-c", script], env={"ADDED": "added"})
    assert out.split() == ["kept", "added"]


def test_a_tool_reads_no_input_from_this_process() -> None:
    # A prompt (a tool's e-mail request) must never wait on, or read, our stdin: run
    # with stdin a pipe that holds an answer, the tool still sees an empty, closed input.
    script = "import sys; print(repr(sys.stdin.read()))"
    read, write = os.pipe()
    os.write(write, b"someone@example.org\n")
    os.close(write)
    saved = os.dup(0)
    os.dup2(read, 0)
    try:
        out = run_tool([sys.executable, "-c", script])
    finally:
        os.dup2(saved, 0)
        os.close(saved)
        os.close(read)
    assert out.strip() == "''"


def test_output_that_is_not_utf8_is_replaced_not_fatal() -> None:
    # Debian's samtools 1.19 prints a Latin-1 byte in `--version`'s compilation details.
    script = r"import sys; sys.stdout.buffer.write(b'samtools 1.19.2\n  flags \xab-O2\xbb\n')"
    out = run_tool([sys.executable, "-c", script])
    assert out.splitlines()[0] == "samtools 1.19.2"
    assert "�" in out


def test_missing_tool_is_named() -> None:
    with pytest.raises(ToolMissing, match="no-such-tool-xyz"):
        run_tool(["no-such-tool-xyz", "--version"])


def test_missing_tool_behind_the_gate_is_named(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gate = tmp_path / "coreslot-gate"
    gate.write_text('#!/usr/bin/env bash\nshift 2\nexec "$@"\n')
    gate.chmod(0o755)
    monkeypatch.setattr(tools, "GATE", gate)
    with pytest.raises(ToolMissing, match="no-such-tool-xyz"):
        run_tool(gated("dna:qc", ["no-such-tool-xyz", "--version"]))


def test_a_test_subprocess_never_reaches_the_machines_gate() -> None:
    """A subprocess a test starts (a bootstrap step's `python -m genokeeper...`) imports
    `genokeeper.tools` afresh: the conftest's `GENOKEEPER_GATE` must reach it."""
    printed = subprocess.run(
        [sys.executable, "-c", "import genokeeper.tools as t; print(t.GATE)"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert printed.stdout.strip() == "/nonexistent-gate"


def test_a_test_never_reaches_the_machines_gate_in_process() -> None:
    """`genokeeper.tools` read `GENOKEEPER_GATE` when pytest imported it, before the conftest
    set it: the conftest patches `tools.GATE` itself too."""
    assert gated("x", ["true"]) == ["true"]
