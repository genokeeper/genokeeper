"""The Y haplogroup from Yleaf (`genokeeper.yhaplogroup`): its HG002 output
(tests/fixtures/hg002) as the grading reads it, and its container run, replaced by a stub that
checks the arguments and writes the fixture into the output directory."""

import shutil
from pathlib import Path

import pytest

from genokeeper import yhaplogroup
from genokeeper.yhaplogroup import (
    YHaplogroup,
    derived_markers,
    read_yleaf,
    run_yleaf,
)

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
PREDICTION = FIXTURES / "hg002.yleaf.hg_prediction.txt"
OUT_FILE = FIXTURES / "hg002.yleaf.out.txt"


def test_read_yleaf_hg002() -> None:
    assert read_yleaf(PREDICTION) == YHaplogroup("J-C134621", "C134621", 1.0, 170393)
    # The haplogroup HG002 is in is J-L816: its L816 marker is read derived.
    assert "L816" in derived_markers(OUT_FILE)


def test_derived_markers_are_the_derived_state_only() -> None:
    derived = derived_markers(OUT_FILE)
    ancestral = [
        row.split("\t")[2]
        for row in OUT_FILE.read_text().splitlines()[1:]
        if row.split("\t")[10] == "A"
    ]
    assert ancestral and not derived & set(ancestral)
    assert "L816" in derived and "FGC24756" not in derived


def test_read_yleaf_empty_none(tmp_path: Path) -> None:
    header = PREDICTION.read_text().splitlines()[0]
    path = tmp_path / "hg_prediction.hg"
    path.write_text(header + "\n")
    assert read_yleaf(path) is None
    path.write_text(header + "\nsample\tNA\tNA\t0\t0\tNA\tNA\tNA\tNA\n")
    assert read_yleaf(path) is None


def test_read_yleaf_two_samples_refuses(tmp_path: Path) -> None:
    text = PREDICTION.read_text()
    path = tmp_path / "hg_prediction.hg"
    path.write_text(text + text.splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="one sample"):
        read_yleaf(path)


def test_run_yleaf_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = (
        tmp_path / "in" / "x.cram",
        tmp_path / "ref" / "g.fa",
        tmp_path / "out",
    )
    for path in (cram, fasta):
        path.parent.mkdir()
        path.write_bytes(b"")
    out.mkdir()
    calls: list[tuple[str, list[str], list[tuple[Path, str, bool]], str | None]] = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        calls.append((name, list(args), list(mounts), label))
        # Yleaf 4.1.6 prompts (and dies on the closed stdin) when `-o` exists.
        target = args[args.index("-o") + 1]
        assert target.startswith("/out/"), f"-o {target} is the mount root"
        folder = out / target.removeprefix("/out/")
        assert not folder.exists(), f"-o {target} already exists"
        folder.mkdir()
        shutil.copy(PREDICTION, folder / "hg_prediction.hg")

    monkeypatch.setattr(yhaplogroup, "run_container", fake)
    result = run_yleaf(cram, fasta, out, threads=3)
    assert result == out / "yleaf" / "hg_prediction.hg"
    assert read_yleaf(result) is not None
    ((name, args, mounts, label),) = calls
    assert (name, label) == ("yleaf", "dna:yleaf")
    assert args[0] == "Yleaf"
    assert args[args.index("-cram") + 1] == "/cram/x.cram"
    assert args[args.index("-cr") + 1] == "/ref/g.fa"
    assert args[args.index("--ref-fasta") + 1] == "/ref/g.fa"
    assert args[args.index("-rg") + 1] == "hg38"
    assert args[args.index("-o") + 1] == "/out/yleaf"
    assert args[args.index("-t") + 1] == "3"
    assert "--no-update-check" in args
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }
    with pytest.raises(FileExistsError, match="yleaf"):
        run_yleaf(cram, fasta, out, threads=3)
    assert len(calls) == 1
