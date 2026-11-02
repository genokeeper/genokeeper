"""External command-line tools (bcftools, samtools) and the errors raised around them."""

import os
import shutil
import subprocess
from collections.abc import Mapping
from pathlib import Path

from genokeeper.env import getenv

# The machine's core budget gate: CPU-heavy steps run as `<GATE> dna:annotate -- <command>`
# when it exists. `$GENOKEEPER_GATE` (`$GENOMELAB_GATE` before the rename) names another
# path: tests set it (tests/conftest.py) so that the subprocesses they start never reach the
# machine's gate.
GATE = Path(getenv("GATE") or Path.home() / ".local" / "bin" / "coreslot-gate")


def gated(label: str, args: list[str]) -> list[str]:
    """`args` behind `<GATE> <label> --` when the machine's gate exists, else `args` as given.

    One gate per unit of work: a pipeline goes behind it whole (as one `bash -c`), never
    stage by stage, or a later stage waits for the slot an earlier one holds while that one
    blocks on the full pipe."""
    return [str(GATE), label, "--", *args] if GATE.exists() else args


class ToolMissing(RuntimeError):
    """A required command-line tool is not on PATH."""

    def __init__(self, tool: str) -> None:
        super().__init__(
            f"{tool} is not on PATH; install it (e.g. `brew install {tool}`)"
        )
        self.tool = tool


def require(*tools: str) -> None:
    """Raise `ToolMissing` for the first of `tools` not on PATH."""
    for tool in tools:
        if shutil.which(tool) is None:
            raise ToolMissing(tool)


class ToolFailed(subprocess.CalledProcessError):
    """A tool exited non-zero; the message carries what it wrote to stderr."""

    def __str__(self) -> str:
        stderr = (self.stderr or "").strip()
        return f"{super().__str__()}\n{stderr}" if stderr else super().__str__()


def run_tool(
    args: list[str], *, env: Mapping[str, str] | None = None, cwd: Path | None = None
) -> str:
    """Run `args`, whose first element names the tool, and return its stdout; `env` is
    added to this process's environment for the tool, which runs in `cwd` when given.

    The tool reads nothing from this process: its stdin is `/dev/null`, so a prompt
    (a tool asking for an e-mail on a terminal) can neither wait nor be answered. Bytes
    that are not UTF-8 come back replaced (Debian's samtools prints one in `--version`).
    Raises `ToolMissing` when the tool (behind the gate, the one after `--`) is absent and
    `ToolFailed` (a `CalledProcessError` whose message includes stderr) when it fails.
    """
    require(args[3] if args[0] == str(GATE) else args[0])
    full_env = None if env is None else {**os.environ, **env}
    try:
        return subprocess.run(
            args,
            check=True,
            capture_output=True,
            text=True,
            errors="replace",
            env=full_env,
            cwd=cwd,
            stdin=subprocess.DEVNULL,
        ).stdout
    except subprocess.CalledProcessError as failed:
        raise ToolFailed(
            failed.returncode, failed.cmd, failed.stdout, failed.stderr
        ) from failed
