"""Truvari between the two callers and against the truth (`genokeeper.deletions.truvari`):
its HG002 chr20 output folders (tests/fixtures/hg002, written into a folder as Truvari names
its files) read, our own calls VCF written, and the container runs, whose container is
replaced by a stub that records the arguments. `write_calls` runs bcftools for real, so
these tests need bcftools.

The counts are the tests' own reading of the fixtures (`_count`), never the code under
test."""

import gzip
import subprocess
from pathlib import Path

import pytest

from genokeeper.deletions import truvari
from genokeeper.deletions.screen import Deletion
from genokeeper.deletions.truvari import (
    AGAINST_TRUTH,
    BETWEEN_CALLERS,
    against_truth,
    judged,
    match_callers,
    read_union,
    write_calls,
)

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
PARTS = ("tp-base", "tp-comp", "fn", "fp")


def _folder(tmp_path: Path, run: str) -> Path:
    """The `hg002.<run>.*` fixtures as Truvari writes them: `<part>.vcf.gz`."""
    folder = tmp_path / run
    folder.mkdir()
    for part in PARTS:
        text = (FIXTURES / f"hg002.{run}.{part}.vcf.txt").read_bytes()
        (folder / f"{part}.vcf.gz").write_bytes(gzip.compress(text))
    return folder


def _count(run: str, part: str) -> int:
    path = FIXTURES / f"hg002.{run}.{part}.vcf.txt"
    return sum(1 for line in path.open() if not line.startswith("#"))


def test_read_union_counts(tmp_path: Path) -> None:
    union = read_union(_folder(tmp_path, "callers"))
    both, dysgu, delly = (
        _count("callers", "tp-base"),
        _count("callers", "fn"),
        _count("callers", "fp"),
    )
    assert (both, dysgu, delly) == (57, 43, 30)
    assert len(union) == both + dysgu + delly
    keys = [frozenset(d.callers) for d in union]
    assert keys.count(frozenset({"dysgu", "delly"})) == both
    assert keys.count(frozenset({"dysgu"})) == dysgu
    assert keys.count(frozenset({"delly"})) == delly
    assert all(d.depth_ratio is None for d in union if set(d.callers) == {"dysgu"})


def test_read_union_joins_partners(tmp_path: Path) -> None:
    union = read_union(_folder(tmp_path, "callers"))
    both = [d for d in union if "dysgu" in d.callers and "delly" in d.callers]
    assert len(both) == _count("callers", "tp-base")
    # By hand: Dysgu's 6159026 (tp-base, MatchId=2.0,2.0) and DELLY's DEL00000002 (tp-comp,
    # MatchId=2.0,2.0; RCL 565, RC 94, RCR 475); Dysgu's coordinates are kept.
    (d,) = (d for d in both if d.callers["dysgu"].startswith("6159026 "))
    assert d.callers == {"dysgu": "6159026 QUAL 29", "delly": "DEL00000002 QUAL 960"}
    assert (d.chrom, d.pos, d.end, d.gt) == ("chr20", 613784, 613837, "1/1")
    assert d.depth_ratio == pytest.approx(94 / (565 + 475))


def _fasta(tmp_path: Path) -> Path:
    fasta = tmp_path / "g.fa"
    fasta.write_bytes(b"")
    Path(f"{fasta}.fai").write_text(
        "chr1\t248956422\t0\t60\t61\nchr2\t242193529\t0\t60\t61\n"
    )
    return fasta


def test_write_calls_carries_class(tmp_path: Path) -> None:
    later = Deletion("chr2", 5000, 7000, "G", "0/1", {"delly": "D1 QUAL 9"}, 0.4)
    earlier = Deletion("chr1", 100, 260, "A", "1/1", {"dysgu": "7 QUAL 3"}, None)
    out = write_calls(
        [(later, "a"), (earlier, "b")], _fasta(tmp_path), tmp_path / "calls.vcf.gz"
    )
    assert Path(f"{out}.tbi").is_file()
    shown = subprocess.run(
        ["bcftools", "view", "-H", str(out)], check=True, capture_output=True, text=True
    ).stdout
    rows = [line.split("\t") for line in shown.splitlines()]
    assert [(r[0], r[1], r[3], r[4], r[6]) for r in rows] == [
        ("chr1", "100", "A", "<DEL>", "PASS"),
        ("chr2", "5000", "G", "<DEL>", "PASS"),
    ]
    infos = [dict(item.split("=") for item in r[7].split(";")) for r in rows]
    assert infos == [
        {"END": "260", "SVTYPE": "DEL", "SVLEN": "-160", "CLASS": "b"},
        {"END": "7000", "SVTYPE": "DEL", "SVLEN": "-2000", "CLASS": "a"},
    ]
    assert [r[9] for r in rows] == ["1/1", "0/1"]
    header = subprocess.run(
        ["bcftools", "view", "-h", str(out)], check=True, capture_output=True, text=True
    ).stdout
    assert "##contig=<ID=chr2,length=242193529>" in header


def test_write_calls_refuses_a_class_that_breaks_info(tmp_path: Path) -> None:
    d = Deletion("chr1", 100, 260, "A", "1/1", {"dysgu": "7 QUAL 3"}, None)
    with pytest.raises(ValueError, match="CLASS"):
        write_calls([(d, "a;b")], _fasta(tmp_path), tmp_path / "calls.vcf.gz")


class _Fake:
    """Stands for `run_container`: records each call."""

    def __init__(self) -> None:
        self.calls: list[
            tuple[str, list[str], list[tuple[Path, str, bool]], str | None]
        ] = []

    def __call__(
        self,
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        self.calls.append((name, list(args), list(mounts), label))


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> _Fake:
    stub = _Fake()
    monkeypatch.setattr(truvari, "run_container", stub)
    return stub


def _contains(args: list[str], options: tuple[str, ...]) -> bool:
    return any(args[i : i + len(options)] == list(options) for i in range(len(args)))


def _mounted(mounts: list[tuple[Path, str, bool]], path: Path, inside: str) -> bool:
    """`inside` (a container path) is `path` under one of `mounts`."""
    return any(
        inside == f"{target}/{path.name}" and host == path.parent
        for host, target, _ in mounts
    )


def test_match_callers_args(tmp_path: Path, fake: _Fake) -> None:
    dysgu, delly, out = (
        tmp_path / "a" / "dysgu.vcf.gz",
        tmp_path / "b" / "delly.vcf.gz",
        tmp_path / "out",
    )
    assert match_callers(dysgu, delly, out, min_length=50) == out / "both"
    ((name, args, mounts, label),) = fake.calls
    assert (name, label) == ("truvari", "dna:truvari")
    assert args[:2] == ["truvari", "bench"]
    assert _mounted(mounts, dysgu, args[args.index("-b") + 1])
    assert _mounted(mounts, delly, args[args.index("-c") + 1])
    assert (out, "/out", False) in mounts
    assert args[args.index("-o") + 1] == "/out/both"
    assert _contains(args, BETWEEN_CALLERS)
    assert args[args.index("--sizemin") + 1] == "50"
    assert args[args.index("--sizefilt") + 1] == "25"
    assert "-f" not in args and "--reference" not in args
    assert all(ro for _, t, ro in mounts if t != "/out")


def test_against_truth_args(tmp_path: Path, fake: _Fake) -> None:
    calls, truth, bench, out = (
        tmp_path / "w" / "calls.vcf.gz",
        tmp_path / "giab" / "truth.vcf.gz",
        tmp_path / "giab" / "bench.bed",
        tmp_path / "out",
    )
    assert against_truth(calls, truth, bench, out) == out / "truth"
    ((name, args, mounts, label),) = fake.calls
    assert (name, label) == ("truvari", "dna:truvari")
    assert args[:2] == ["truvari", "bench"]
    assert _mounted(mounts, truth, args[args.index("-b") + 1])
    assert _mounted(mounts, calls, args[args.index("-c") + 1])
    assert _mounted(mounts, bench, args[args.index("--includebed") + 1])
    assert args[args.index("-o") + 1] == "/out/truth"
    assert _contains(args, AGAINST_TRUTH)
    assert not {"-f", "--reference", "--refine"} & set(args)
    assert all(ro for _, t, ro in mounts if t != "/out")


def test_options_as_decided() -> None:
    assert BETWEEN_CALLERS == (
        "--pctseq", "0", "--pctsize", "0.5", "--pctovl", "0.5", "--passonly",
        "--pick", "single", "--sizemax", "-1",
    )  # fmt: skip
    assert AGAINST_TRUTH == (
        "--pctseq", "0", "--pick", "ac", "--passonly", "-r", "2000", "-C", "5000",
        "--sizemin", "35", "--sizemax", "-1",
    )  # fmt: skip


def test_judged_hg002(tmp_path: Path) -> None:
    calls, truth = judged(_folder(tmp_path, "truth"))
    true_calls, false_calls = _count("truth", "tp-comp"), _count("truth", "fp")
    found, missed = _count("truth", "tp-base"), _count("truth", "fn")
    assert (true_calls, false_calls, found, missed) == (84, 24, 84, 293)
    assert [ok for _, _, ok in calls].count(True) == true_calls
    assert [ok for _, _, ok in calls].count(False) == false_calls
    assert {c for c, _, _ in calls} == {"x"}
    assert [ok for _, ok in truth].count(True) == found
    assert [ok for _, ok in truth].count(False) == missed
    # By hand: fp's first record, x4 at chr20:2883370, END=2883840; fn's first, truth74 at
    # chr20:393038, END=393161.
    assert ("x", 2883840 - 2883370, False) in calls
    assert (393161 - 393038, False) in truth
    assert all(length > 0 for _, length, _ in calls)


# --- any SV type (slice 3d) -----------------------------------------------------------------


def test_write_calls_ins_and_dup(tmp_path: Path) -> None:
    """An insertion at its POS (`END=POS`) and a duplication, each with a positive `SVLEN`;
    bcftools reads both."""
    ins = Deletion(
        "chr1", 500, 501, "C", "0/1", {"dysgu": "3 QUAL 9"}, None,
        svtype="INS", svlen=320,
    )  # fmt: skip
    dup = Deletion(
        "chr1", 1000, 4000, "G", "0/1", {"dysgu": "4 QUAL 9"}, None, svtype="DUP"
    )
    out = write_calls(
        [(dup, "d"), (ins, "i")], _fasta(tmp_path), tmp_path / "calls.vcf.gz"
    )
    shown = subprocess.run(
        ["bcftools", "view", "-H", str(out)], check=True, capture_output=True, text=True
    ).stdout
    rows = [line.split("\t") for line in shown.splitlines()]
    assert [(r[1], r[3], r[4]) for r in rows] == [
        ("500", "C", "<INS>"),
        ("1000", "G", "<DUP>"),
    ]
    infos = [dict(item.split("=") for item in r[7].split(";")) for r in rows]
    assert infos == [
        {"END": "500", "SVTYPE": "INS", "SVLEN": "320", "CLASS": "i"},
        {"END": "4000", "SVTYPE": "DUP", "SVLEN": "3000", "CLASS": "d"},
    ]


def test_read_union_insertions(tmp_path: Path) -> None:
    """Task 6's Truvari run of Dysgu's insertions (base) against INSurVeyor's (comparison)."""
    names = ("dysgu", "insurveyor")
    union = read_union(_folder(tmp_path, "ins"), names, svtype="INS")
    both, dysgu, insurveyor = (_count("ins", p) for p in ("tp-base", "fn", "fp"))
    assert (both, dysgu, insurveyor) == (34, 137, 49)
    keys = [frozenset(d.callers) for d in union]
    assert keys.count(frozenset(names)) == both
    assert keys.count(frozenset({"dysgu"})) == dysgu
    assert keys.count(frozenset({"insurveyor"})) == insurveyor
    assert all(d.svtype == "INS" and d.end == d.pos + 1 for d in union)
    assert all(d.svlen is not None and d.svlen > 0 for d in union)


def test_match_callers_names(tmp_path: Path, fake: _Fake) -> None:
    base, comp, out = (
        tmp_path / "a" / "dysgu.ins.vcf.gz",
        tmp_path / "b" / "insurveyor.vcf.gz",
        tmp_path / "out",
    )
    found = match_callers(base, comp, out, min_length=50, names=("dysgu", "insurveyor"))
    assert found == out / "both"
    ((_, args, mounts, _),) = fake.calls
    assert _mounted(mounts, base, args[args.index("-b") + 1])
    assert _mounted(mounts, comp, args[args.index("-c") + 1])
    # INSurVeyor's floor is the screen's: only DELLY's is halved.
    assert args[args.index("--sizefilt") + 1] == "50"
