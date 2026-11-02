"""The calibration tables and the one command that rebuilds them (`genokeeper calibrate`).

Every table is a row of `TABLES`: its name, the tables it needs built first, the versions it
is measured under (`versions(lab)`, the run's own computation) and its builder
(`build(lab, threads)`, which raises `RunRefused` naming a missing input and the command
that supplies it). A table lives in `folder(lab, name)`, written whole or not at all by
`write_table` with a `versions.json` beside what its writer left (and the writer's optional
`report.txt`). `open_table` is the one staleness check: no folder or no `versions.json` is
`NOT_CALIBRATED`, a versions key that differs is "stale: <key>, ... changed; run `genokeeper
calibrate`". `calibrate` builds each stale table, every table before those that need it.

3c's rows are `small_variants` (2b's HG003 table), `hg002_cells` (needing it), `stage3`
(3a), `deletions` (3b) and `foundation` (`calibrate.foundation`); 3d's are `cnv`, `insertions`
and `paralogs`, none needing another. Their versions and builders are
`genokeeper.calibrate.build.ROWS`, reached here only when called, because `build` imports
the pipeline, which imports this module. A later slice adds its rows there, never a command
or a flag.

The check's anchor comes from the versions each row computes from what the lab runs, which
a table measured under others cannot move because `open_table` never reads the expected
versions from the table it judges.
"""

import json
import os
import secrets
import shutil
import sys
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

from genokeeper.calibrate import calibration_dir

VERSIONS = "versions.json"
REPORT = "report.txt"
RUN = "run `genokeeper calibrate`"
NOT_CALIBRATED = f"not calibrated: {RUN}"


@dataclass(frozen=True)
class Table:
    """A calibration table: its name, the tables it needs, its versions and its builder."""

    name: str
    needs: tuple[str, ...]
    versions: Callable[[Path], dict[str, str]]
    build: Callable[[Path, int], None]


@dataclass(frozen=True)
class Current:
    """A table measured under the versions asked for, in `folder`."""

    folder: Path
    versions: dict[str, str]


@dataclass(frozen=True)
class Stale:
    """A table that is missing or measured under other versions, and why."""

    reason: str


def _row(name: str, needs: tuple[str, ...] = ()) -> Table:
    def versions(lab: Path) -> dict[str, str]:
        from genokeeper.calibrate.build import ROWS  # build imports the pipeline

        return ROWS[name][0](lab)

    def build(lab: Path, threads: int) -> None:
        from genokeeper.calibrate.build import ROWS

        ROWS[name][1](lab, threads)

    return Table(name, needs, versions, build)


TABLES: dict[str, Table] = {
    table.name: table
    for table in (
        _row("small_variants"),
        _row("hg002_cells", ("small_variants",)),
        _row("stage3"),
        _row("deletions"),
        _row("cnv"),
        _row("insertions"),
        _row("paralogs"),
        _row("foundation"),
        _row("slice3f"),
    )
}


def folder(lab: Path, name: str) -> Path:
    """Where the table `name` lives."""
    return calibration_dir(lab) / name


def write_table(
    lab: Path, name: str, write: Callable[[Path], None], versions: Mapping[str, str]
) -> Path:
    """`write(staging)` into a hidden staging folder beside the table's folder, then
    `versions.json`, swapped in whole: a reader sees the old table or the new one, never
    half of one, and a writer that fails leaves the old one. Returns the folder."""
    target = folder(lab, name)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = target.with_name(f".{name}.{secrets.token_hex(4)}")
    staging.mkdir()
    try:
        write(staging)
        (staging / VERSIONS).write_text(
            json.dumps(dict(versions), indent=1, sort_keys=True)
        )
        old = target.with_name(f".{name}.old-{secrets.token_hex(4)}")
        if target.exists():
            os.replace(target, old)
        os.replace(staging, target)
        shutil.rmtree(old, ignore_errors=True)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return target


def _one_line(error: BaseException) -> str:
    return " ".join(str(error).split())


def open_table(
    lab: Path, name: str, *, versions: Mapping[str, str] | None = None
) -> Current | Stale:
    """The table `name` when it was measured under `versions` (default: its row's,
    `TABLES[name].versions(lab)`; a run passes the ones it computed itself); otherwise why
    not. Versions that cannot be computed (an input missing) make it stale with the cause."""
    target = folder(lab, name)
    try:
        recorded = json.loads((target / VERSIONS).read_text())
    except FileNotFoundError:
        return Stale(NOT_CALIBRATED)
    except (OSError, ValueError) as e:
        return Stale(f"unreadable {VERSIONS} ({_one_line(e)}); {RUN}")
    if versions is None:
        try:
            versions = TABLES[name].versions(lab)
        except (RuntimeError, OSError, ValueError) as e:  # RuntimeError: RunRefused
            return Stale(f"versions not computed ({_one_line(e)}); {RUN}")
    differ = sorted(
        key for key in {*recorded, *versions} if recorded.get(key) != versions.get(key)
    )
    if differ:
        return Stale(f"stale: {', '.join(differ)} changed; {RUN}")
    return Current(target, dict(recorded))


def _ordered(names: set[str]) -> list[str]:
    """`names` with every table before those that need it, registry order breaking ties."""
    ordered: list[str] = []
    left = [name for name in TABLES if name in names]
    while left:
        ready = next(
            (
                name
                for name in left
                if all(n in ordered or n not in names for n in TABLES[name].needs)
            ),
            None,
        )
        if ready is None:
            raise ValueError(f"calibration tables need each other: {', '.join(left)}")
        ordered.append(ready)
        left.remove(ready)
    return ordered


def stale_tables(lab: Path) -> list[str]:
    """The stale tables' names, every table before those that need it."""
    return _ordered(
        {name for name in TABLES if isinstance(open_table(lab, name), Stale)}
    )


def _needed(name: str) -> list[str]:
    """Every table `name` needs, directly or through another, in registry order."""
    found: set[str] = set()
    todo = list(TABLES[name].needs)
    while todo:
        need = todo.pop()
        if need not in found:
            found.add(need)
            todo += TABLES[need].needs
    return [n for n in TABLES if n in found]


def _report(lab: Path, name: str) -> None:
    path = folder(lab, name) / REPORT
    if path.is_file():
        print(path.read_text(), end="")


def calibrate(lab: Path, *, only: str | None, report: bool, threads: int) -> int:
    """Build each stale table, every table before those that need it (with `only`, that
    table alone, refused while a table it needs is stale); a table whose builder refuses
    is skipped with one line, and so is every table that needs it. A current table prints
    "<name>: current". With `report`, each built table's report is printed (with `only`,
    a current one's too). Returns 2 when any table was skipped or refused, else 0."""
    from genokeeper.pipeline import RunRefused  # the pipeline imports this module

    stale = stale_tables(lab)
    if only is not None:
        if only not in TABLES:
            print(
                f"genokeeper calibrate: no table {only}; one of {', '.join(TABLES)}",
                file=sys.stderr,
            )
            return 2
        behind = [n for n in _needed(only) if n in stale]
        if behind:
            print(
                f"genokeeper calibrate: {only} needs {', '.join(behind)}, which is not"
                f" current; {RUN}",
                file=sys.stderr,
            )
            return 2
        names = [only]
    else:
        names = list(_ordered(set(TABLES)))
    skipped: set[str] = set()
    code = 0
    for name in names:
        if name in skipped:
            continue
        if name not in stale:
            print(f"{name}: current")
            if report and only is not None:
                _report(lab, name)
            continue
        try:
            TABLES[name].build(lab, threads)
        except RunRefused as e:
            dependents = [n for n in TABLES if name in _needed(n)]
            skipped |= set(dependents)
            said = f"; skipped {', '.join(dependents)}" if dependents else ""
            print(
                f"genokeeper calibrate: {name}: {_one_line(e)}{said}", file=sys.stderr
            )
            code = 2
            continue
        if report:
            _report(lab, name)
    return code
