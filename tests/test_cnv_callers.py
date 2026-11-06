"""The depth callers and Dysgu's DUP junctions (`genokeeper.cnv.callers`): their HG002 chr20
outputs (tests/fixtures/hg002, Task 6) read as events, and the container runs, whose
container is replaced by a stub that checks the arguments and writes the output.

The counts are the tests' own reading of the fixtures (`_records`), never the code under
test."""

from collections.abc import Iterator
from pathlib import Path

import pytest

from genokeeper.cnv import callers
from genokeeper.cnv.callers import (
    CNV_LABEL,
    CNVPYTOR_LABEL,
    expected_ploidy,
    read_cnvpytor,
    read_delly_cnv,
    read_junctions,
    run_cnvpytor,
    run_delly_cnv,
)
from genokeeper.config import load_grading
from genokeeper.deletions.exons import PRIMARY

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
DELLY_CNV = FIXTURES / "hg002.dellycnv.vcf.txt"
CNVPYTOR = FIXTURES / "hg002.cnvpytor.calls.txt"
DYSGU_DUP = FIXTURES / "hg002.dysgu.dup.vcf.txt"
MALE = expected_ploidy(load_grading(), male=True)
FEMALE = expected_ploidy(load_grading(), male=False)


def _records(path: Path) -> Iterator[tuple[list[str], dict[str, str], dict[str, str]]]:
    """Each record's columns, INFO and first sample's FORMAT values."""
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        pairs = (item.partition("=") for item in cols[7].split(";"))
        info = {name: value for name, _, value in pairs}
        values = dict(zip(cols[8].split(":"), cols[9].split(":"), strict=True))
        yield cols, info, values


def test_read_delly_cnv_hg002() -> None:
    """Every PASS record off two copies on chr20, a loss below and a gain above; DELLY
    refines a record's ends after its `-z` window, so a few are under 1,000 bp (the
    screen's floor cuts them)."""
    found = read_delly_cnv(DELLY_CNV, ploidy=lambda chrom, pos: 2)
    expected = [
        (cols, info, int(values["CN"]))
        for cols, info, values in _records(DELLY_CNV)
        if cols[6] == "PASS" and values["CN"] != "2"
    ]
    assert 0 < len(expected) == len(found)
    assert len(expected) < sum(1 for _ in _records(DELLY_CNV))
    for d, (cols, info, cn) in zip(found, expected, strict=True):
        assert (d.chrom, d.pos, d.end) == (cols[0], int(cols[1]) - 1, int(info["END"]))
        assert d.svtype == ("DEL" if cn < 2 else "DUP")
        assert d.callers == {"dellycnv": f"{cols[2]} CN {cn}"}
        assert d.gt == ("1/1" if cn == 0 else "0/1")
        assert d.length == int(info["END"]) - int(cols[1]) + 1
    assert {d.svtype for d in found} == {"DEL", "DUP"}
    assert sum(d.length >= 1000 for d in found) > len(found) * 0.9


def test_haploid_x_not_loss(tmp_path: Path) -> None:
    """A male's chrX outside the PARs has one copy: CN 1 there is normal, CN 2 a gain; in
    PAR1, and on a female's chrX, CN 1 is a loss."""
    header = "".join(line for line in DELLY_CNV.open() if line.startswith("#"))
    cols, _, _ = next(_records(DELLY_CNV))
    rows = [
        ("chrX", "5000001", "X1", "1"),
        ("chrX", "6000001", "X2", "2"),
        ("chrX", "100001", "PAR", "1"),
    ]
    lines = []
    for chrom, pos, vid, cn in rows:
        info = f"IMPRECISE;SVTYPE=CNV;END={int(pos) + 4999}"
        sample = cols[9].split(":")
        sample[cols[8].split(":").index("CN")] = cn
        lines.append(
            "\t".join([chrom, pos, vid, "N", "<CNV>", "50", "PASS", info, cols[8]])
            + "\t"
            + ":".join(sample)
        )
    path = tmp_path / "x.vcf"
    path.write_text(header + "\n".join(lines) + "\n")
    male = read_delly_cnv(path, ploidy=MALE)
    assert [(d.callers["dellycnv"], d.svtype) for d in male] == [
        ("X2 CN 2", "DUP"),
        ("PAR CN 1", "DEL"),
    ]
    female = read_delly_cnv(path, ploidy=FEMALE)
    assert [(d.callers["dellycnv"], d.svtype) for d in female] == [
        ("X1 CN 1", "DEL"),
        ("PAR CN 1", "DEL"),
    ]
    assert MALE("chr20", 5000001) == 2
    assert MALE("chrY", 3_000_000) == 1


def test_read_cnvpytor_hg002() -> None:
    """Its region is 1-based inclusive: `chr20:1-60000` is `[0, 60000)`."""
    found = read_cnvpytor(CNVPYTOR, min_length=1000)
    rows = [line.split("\t") for line in CNVPYTOR.read_text().splitlines()]
    kept = [r for r in rows if int(r[2]) >= 1000]
    assert 0 < len(kept) == len(found)
    assert all(d.length >= 1000 for d in found)
    for d, row in zip(found, kept, strict=True):
        assert d.svtype == {"deletion": "DEL", "duplication": "DUP"}[row[0]]
        assert d.length == int(row[2])
        assert d.callers == {"cnvpytor": f"{row[1]} e1 {row[4]}"}
        assert d.gt == (
            "1/1" if row[0] == "deletion" and float(row[3]) < 0.25 else "0/1"
        )
    first = found[0]
    assert (first.chrom, first.pos, first.end) == ("chr20", 0, 60000)
    assert {d.svtype for d in found} == {"DEL", "DUP"}
    longer = read_cnvpytor(CNVPYTOR, min_length=50_000)
    assert 0 < len(longer) < len(found)


def test_read_junctions_dysgu() -> None:
    found = read_junctions(DYSGU_DUP)
    passed = [
        (cols, info) for cols, info, _ in _records(DYSGU_DUP) if cols[6] == "PASS"
    ]
    assert 0 < len(passed) == len(found) < sum(1 for _ in _records(DYSGU_DUP))
    assert all(d.svtype == "DUP" for d in found)
    for d, (cols, info) in zip(found, passed, strict=True):
        assert (d.pos, d.end) == (int(cols[1]), int(info["END"]))
        assert set(d.callers) == {"dysgu"}


def _cram_and_fasta(tmp_path: Path, contigs: list[str]) -> tuple[Path, Path, Path]:
    cram, fasta, out = (
        tmp_path / "in" / "x.cram",
        tmp_path / "ref" / "g.fa",
        tmp_path / "out",
    )
    for path in (cram, fasta):
        path.parent.mkdir()
        path.write_bytes(b"")
    Path(f"{fasta}.fai").write_text(
        "".join(f"{c}\t{1000 * (i + 1)}\t0\t60\t61\n" for i, c in enumerate(contigs))
    )
    out.mkdir()
    return cram, fasta, out


TEMPLATE = "chr1\t0\t10000\ttelomere\nchrM\n"
CONTIGS = ["chr1", "chr20", "chrX", "chrY", "chrM", "chrEBV"]
MASKED = "chr1\t5000\t6000\nchr20\t100\t900\n"


class _Delly:
    """Stands for `run_container` and `run_tool`: copies the template and records each
    `delly cnv` container's arguments, mounts, label and exclude list."""

    def __init__(self) -> None:
        self.runs: list[
            tuple[list[str], list[tuple[Path, str, bool]], str | None, str]
        ] = []
        self.tools: list[list[str]] = []

    def container(
        self,
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert name == "delly"
        (out,) = (host for host, target, ro in mounts if target == "/out" and not ro)
        if args[0] == "cp":
            (out / "exclude.template.tsv").write_text(TEMPLATE)
            return
        self.runs.append(
            (list(args), list(mounts), label, (out / "exclude.tsv").read_text())
        )
        (out / "delly.bcf").write_bytes(b"")

    def tool(self, args: list[str], **_: object) -> str:
        self.tools.append(list(args))
        return ""


@pytest.mark.parametrize(
    ("haploid", "shards", "outside"),
    [
        (("chrX", "chrY"), PRIMARY, ["chrM", "chrEBV"]),
        ((), [c for c in PRIMARY if c != "chrY"], ["chrY", "chrM", "chrEBV"]),
    ],
)
def test_run_delly_cnv_args(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    haploid: tuple[str, ...],
    shards: list[str],
    outside: list[str],
) -> None:
    """One container over the genome, whatever `shards` holds: DELLY reads `-x` only on
    its finished calls, so its exclude list is what keeps them on `shards`."""
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    masked = tmp_path / "masked.bed"
    masked.write_text(MASKED)
    fake = _Delly()
    monkeypatch.setattr("genokeeper.deletions.callers.run_container", fake.container)
    monkeypatch.setattr(callers, "run_container", fake.container)
    monkeypatch.setattr(callers, "run_tool", fake.tool)
    found = run_delly_cnv(
        cram, fasta, out, shards=shards, haploid=haploid, masked=masked,
        min_length=1000,
    )  # fmt: skip
    assert found == out / "dellycnv.vcf.gz"
    ((args, mounts, label, exclude),) = fake.runs
    assert label == CNV_LABEL == "dna:cnv"
    assert args[:2] == ["delly", "cnv"]
    assert args[args.index("-z") + 1] == "1000"
    assert args[args.index("-g") + 1] == "/ref/g.fa"
    assert args[args.index("-x") + 1] == "/out/exclude.tsv"
    assert args[args.index("-o") + 1] == "/out/delly.bcf"
    assert args[-1] == "/cram/x.cram"
    if haploid:
        assert args[args.index("--haploid-chr") + 1] == "chrX,chrY"
    else:
        assert "--haploid-chr" not in args
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }
    assert (out, "/out", False) in mounts
    assert exclude == TEMPLATE + MASKED + "".join(f"{c}\n" for c in outside)
    view, index = fake.tools
    assert view[:2] == ["bcftools", "view"]
    assert view[view.index("-f") + 1] == "PASS"
    assert view[view.index("-o") + 1] == str(found)
    assert view[-1] == str(out / "delly.bcf")
    assert index[:2] == ["bcftools", "index"] and index[-1] == str(found)


def test_run_cnvpytor_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    seen: list[tuple[str, list[str], list[tuple[Path, str, bool]], str | None]] = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), list(mounts), label))

    monkeypatch.setattr(callers, "run_container", fake)
    assert run_cnvpytor(cram, fasta, out, threads=4, bin_size=500) == out / "calls.tsv"
    ((name, args, mounts, label),) = seen
    assert (name, label) == ("cnvpytor", CNVPYTOR_LABEL)
    assert CNVPYTOR_LABEL == "dna:cnvpytor"
    assert args[:2] == ["bash", "-c"]
    script = args[2]
    steps = [s.strip() for s in script.split("&&")]
    assert "cd /out" in steps
    (rd,) = [s for s in steps if " -rd " in s]
    assert rd.startswith("cnvpytor -root s.pytor -rd /cram/x.cram -T /ref/g.fa")
    assert "-j 4" in rd
    assert "cnvpytor -root s.pytor -his 500" in steps
    assert "cnvpytor -root s.pytor -partition 500" in steps
    (call,) = [s for s in steps if " -call " in s]
    assert call.startswith("cnvpytor -root s.pytor -call 500 > calls.tsv")
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }
    assert (out, "/out", False) in mounts
