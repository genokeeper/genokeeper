"""How long each built step and each run took on this machine: one JSON line each in
`<work>/timings.jsonl`, for `genokeeper refresh --plan`'s forecast.

A line holds a step name or a run id, the seconds and when it finished: nothing of the
user's genome. A reused step writes nothing; a line cut part-way is skipped.
"""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

FILE = "timings.jsonl"


def record(
    work: Path, seconds: float, *, step: str | None = None, run: str | None = None
) -> None:
    """Append one step's (`step`) or one run's (`run`) duration. Best effort: a file that
    cannot be written (full disk, no permission) never fails the step or run it times."""
    line = {
        "finished": datetime.now(UTC).isoformat(timespec="seconds"),
        "seconds": round(seconds),
    }
    line |= {"step": step} if step is not None else {"run": run}
    try:
        work.mkdir(parents=True, exist_ok=True)
        with (work / FILE).open("a", encoding="utf-8") as handle:
            handle.write("\n" + json.dumps(line) + "\n")  # after a line cut part-way
    except OSError:
        return


def read(work: Path) -> list[dict[str, Any]]:
    """Every line, oldest first; none without the file."""
    path = work / FILE
    if not path.is_file():
        return []
    lines: list[dict[str, Any]] = []
    for text in path.read_text(encoding="utf-8").splitlines():
        if not text:
            continue
        try:
            line = json.loads(text)
            if isinstance(line, dict) and isinstance(line.get("seconds"), int):
                lines.append(line)  # pyright: ignore[reportUnknownArgumentType]
        except ValueError:
            continue
    return lines


def runs(work: Path) -> list[tuple[str, int, list[str]]]:
    """Each run, oldest first: its id, its seconds and the steps it built (the step lines
    since the run before it)."""
    found: list[tuple[str, int, list[str]]] = []
    built: list[str] = []
    for line in read(work):
        if "step" in line:
            built.append(str(line["step"]))
        elif "run" in line:
            found.append((str(line["run"]), int(line["seconds"]), built))
            built = []
    return found


def steps(work: Path) -> dict[str, int]:
    """Each step's seconds at its latest build."""
    return {
        str(line["step"]): int(line["seconds"]) for line in read(work) if "step" in line
    }
