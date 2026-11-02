"""The one calibration command over a registry of tables (`genokeeper.calibrate.tables`):
staleness, dependency order, the whole-or-nothing swap, skipping what needs a refused
table, and the openers' one reason. The registry is a set of fake tables monkeypatched in,
each builder recording its call."""

import json
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from genokeeper import cli
from genokeeper.calibrate import tables
from genokeeper.calibrate.deletions import open_deletions
from genokeeper.calibrate.stage3 import open_stage3
from genokeeper.calibrate.tables import (
    NOT_CALIBRATED,
    Current,
    Stale,
    Table,
    calibrate,
    folder,
    open_table,
    stale_tables,
    write_table,
)
from genokeeper.pipeline import RunRefused, open_calibration

RUN = "run `genokeeper calibrate`"


class Fakes:
    """A registry of fake tables: each one's versions are `wanted[name]`, and its builder
    records its call and writes its report, or raises `RunRefused` when `refuse` names it."""

    def __init__(self, lab: Path, needs: dict[str, tuple[str, ...]]) -> None:
        self.lab = lab
        self.wanted = {name: {"v": "2"} for name in needs}
        self.built: list[str] = []
        self.refuse: set[str] = set()
        self.registry = {
            name: Table(name, deps, self._versions(name), self._builder(name))
            for name, deps in needs.items()
        }

    def _versions(self, name: str) -> Callable[[Path], dict[str, str]]:
        return lambda _lab: dict(self.wanted[name])

    def _builder(self, name: str) -> Callable[[Path, int], None]:
        def build(lab: Path, threads: int) -> None:
            assert threads == 3
            if name in self.refuse:
                raise RunRefused(f"missing input: {name}'s input; run make-{name}")
            self.built.append(name)

            def write(staging: Path) -> None:
                (staging / "report.txt").write_text(f"report of {name}\n")

            write_table(lab, name, write, self.wanted[name])

        return build

    def current(self, *names: str) -> None:
        for name in names:
            write_table(self.lab, name, lambda _s: None, self.wanted[name])


@pytest.fixture
def fakes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Fakes]:
    # `b` is registered before the `a` it needs: the order comes from `needs`.
    found = Fakes(tmp_path / "lab", {"b": ("a",), "a": (), "c": ()})
    monkeypatch.setattr(tables, "TABLES", found.registry)
    yield found


def test_stale_tables_in_dependency_order(fakes: Fakes) -> None:
    fakes.current("c")
    assert stale_tables(fakes.lab) == ["a", "b"]
    fakes.current("a")
    assert stale_tables(fakes.lab) == ["b"]


def test_open_table_reasons(fakes: Fakes) -> None:
    assert open_table(fakes.lab, "a") == Stale(NOT_CALIBRATED)
    assert NOT_CALIBRATED == f"not calibrated: {RUN}"
    folder(fakes.lab, "a").mkdir(parents=True)
    assert open_table(fakes.lab, "a") == Stale(NOT_CALIBRATED)  # no versions.json
    write_table(fakes.lab, "a", lambda _s: None, {"x": "1"})
    fakes.wanted["a"] = {"x": "2"}
    assert open_table(fakes.lab, "a") == Stale(f"stale: x changed; {RUN}")
    fakes.wanted["a"] = {"x": "1"}
    assert open_table(fakes.lab, "a") == Current(folder(fakes.lab, "a"), {"x": "1"})
    # A run passes the versions it computed itself.
    assert open_table(fakes.lab, "a", versions={"x": "3", "y": "1"}) == Stale(
        f"stale: x, y changed; {RUN}"
    )


def test_write_table_is_whole_or_nothing(fakes: Fakes) -> None:
    def old(staging: Path) -> None:
        (staging / "data.txt").write_text("old\n")

    target = write_table(fakes.lab, "a", old, {"v": "1"})
    assert target == folder(fakes.lab, "a")

    def broken(staging: Path) -> None:
        (staging / "data.txt").write_text("new\n")
        raise OSError("disk full")

    with pytest.raises(OSError, match="disk full"):
        write_table(fakes.lab, "a", broken, {"v": "2"})
    assert (target / "data.txt").read_text() == "old\n"
    assert json.loads((target / "versions.json").read_text()) == {"v": "1"}
    assert [p.name for p in target.parent.iterdir()] == ["a"]


def test_calibrate_builds_in_order_and_reports(
    fakes: Fakes, capsys: pytest.CaptureFixture[str]
) -> None:
    assert calibrate(fakes.lab, only=None, report=True, threads=3) == 0
    assert fakes.built == ["a", "b", "c"]
    out = capsys.readouterr().out
    assert "report of a\n" in out and "report of b\n" in out
    assert calibrate(fakes.lab, only=None, report=True, threads=3) == 0
    assert fakes.built == ["a", "b", "c"]  # all current: nothing built again
    out = capsys.readouterr().out
    assert out.splitlines() == ["a: current", "b: current", "c: current"]
    # With --only, a current table's report is printed too.
    assert calibrate(fakes.lab, only="a", report=True, threads=3) == 0
    assert capsys.readouterr().out == "a: current\nreport of a\n"


def test_calibrate_skips_what_needs_a_refused_table(
    fakes: Fakes, capsys: pytest.CaptureFixture[str]
) -> None:
    fakes.current("c")
    fakes.refuse = {"a"}
    assert calibrate(fakes.lab, only=None, report=False, threads=3) == 2
    assert fakes.built == []
    captured = capsys.readouterr()
    (line,) = captured.err.strip().splitlines()
    assert line.startswith("genokeeper calibrate: a: missing input: a's input")
    assert "run make-a" in line and "b" in line.split("run make-a")[1]
    assert "c: current" in captured.out.splitlines()
    assert not folder(fakes.lab, "b").exists()


def test_only_refuses_a_stale_need(
    fakes: Fakes, capsys: pytest.CaptureFixture[str]
) -> None:
    assert calibrate(fakes.lab, only="b", report=False, threads=3) == 2
    assert fakes.built == []
    (line,) = capsys.readouterr().err.strip().splitlines()
    assert line.startswith("genokeeper calibrate: b needs a") and line.endswith(RUN)
    fakes.current("a")
    assert calibrate(fakes.lab, only="b", report=False, threads=3) == 0
    assert fakes.built == ["b"]


def test_old_flags_are_gone() -> None:
    for flag in ("--stage3", "--deletions", "--fastq1=x", "--fastq2=x"):
        with pytest.raises(SystemExit) as raised:
            cli.main(["calibrate", flag])
        assert raised.value.code == 2


def test_openers_use_the_one_reason(tmp_path: Path) -> None:
    """On an empty lab the stage-3 and deletion openers give the one reason; the
    calibration table's opener returns None there (the run refuses after its QC step,
    which records the bases HG003 is downsampled to) and its table's state is that
    reason; a stale table refuses with the stale reason."""
    lab = tmp_path / "lab"
    for opener in (open_stage3, open_deletions):
        table, reason = opener(lab, {})
        assert table is None and reason == NOT_CALIBRATED
    assert NOT_CALIBRATED.endswith(RUN)
    assert open_calibration(lab, {}) is None
    assert open_table(lab, "small_variants", versions={}) == Stale(NOT_CALIBRATED)
    write_table(lab, "stage3", lambda _s: None, {"fasta": "old"})
    table, reason = open_stage3(lab, {"fasta": "new"})
    assert table is None and reason == f"stale: fasta changed; {RUN}"
    write_table(lab, "small_variants", lambda _s: None, {"bwa": "old"})
    with pytest.raises(
        RunRefused, match="stale: bwa changed; run `genokeeper calibrate`$"
    ):
        open_calibration(lab, {"bwa": "new"})
