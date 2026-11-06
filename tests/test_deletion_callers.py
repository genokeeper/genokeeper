"""DELLY's and Dysgu's deletions (`genokeeper.deletions.callers`): their HG002 chr20 outputs
(tests/fixtures/hg002, unfiltered) filtered and read, and the container runs, whose container
is replaced by a stub that checks the arguments and writes the output. `deletions_only` and
`run_dysgu` run bcftools for real, so these tests need bcftools.

The counts are the tests' own reading of the fixtures (`_records`, `_kept`), never the code
under test."""

import gzip
import threading
from collections.abc import Iterator
from pathlib import Path

import pytest

from genokeeper.deletions import callers
from genokeeper.deletions.callers import (
    caller_filter,
    caller_min_length,
    deletions_only,
    filter_expression,
    read_deletions,
    read_records,
    run_delly,
    run_dysgu,
)
from genokeeper.deletions.exons import PRIMARY

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
DELLY = FIXTURES / "hg002.delly.vcf.txt"
DYSGU = FIXTURES / "hg002.dysgu.vcf.txt"
DYSGU_INS = FIXTURES / "hg002.dysgu.ins.vcf.txt"
DYSGU_DUP = FIXTURES / "hg002.dysgu.dup.vcf.txt"
INSURVEYOR = FIXTURES / "hg002.insurveyor.vcf.txt"
MIN_LENGTH = {"delly": 25, "dysgu": 50}
CARRIED = {"0/1", "1/1", "0|1", "1|0", "1|1"}


def _records(path: Path) -> Iterator[tuple[list[str], dict[str, str], dict[str, str]]]:
    """Each record's columns, INFO and first sample's FORMAT values."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            info: dict[str, str] = {}
            for item in cols[7].split(";"):
                name, _, value = item.partition("=")
                info[name] = value
            values = dict(zip(cols[8].split(":"), cols[9].split(":"), strict=True))
            yield cols, info, values


def _length(cols: list[str], info: dict[str, str]) -> int:
    if cols[4].startswith("<"):
        return int(info["END"]) - int(cols[1])
    return len(cols[3]) - 1


def _kept(path: Path, caller: str) -> int:
    """The fixture's records `deletions_only` should keep, by the test's own reading."""
    return sum(
        1
        for cols, info, values in _records(path)
        if info.get("SVTYPE") == "DEL"
        and cols[6] == "PASS"
        and (caller != "delly" or values["FT"] == "PASS")
        and values["GT"] in CARRIED
        and _length(cols, info) >= MIN_LENGTH[caller]
    )


def _header(path: Path) -> str:
    return "".join(line for line in path.open() if line.startswith("#"))


def test_caller_floors_one_source() -> None:
    """Each caller's floor under the screen's 50 bp (`MIN_LENGTH`), and the filter the step
    key and the table's versions record at it; an unknown caller is refused."""
    for caller, floor in MIN_LENGTH.items():
        assert caller_min_length(caller, 50) == floor
        assert caller_filter(caller, 50) == filter_expression(caller, floor)
    with pytest.raises(ValueError, match="manta"):
        caller_min_length("manta", 50)


def test_delly_depth_ratio(tmp_path: Path) -> None:
    (cols, _, values) = next(r for r in _records(DELLY) if r[0][2] == "DEL00000002")
    assert (values["RC"], values["RCL"], values["RCR"]) == ("94", "565", "475")
    keys = cols[8].split(":")
    flat = [
        "0" if k in ("RCL", "RCR") else v
        for k, v in zip(keys, cols[9].split(":"), strict=True)
    ]
    left = [
        "0" if k == "RCL" else v for k, v in zip(keys, cols[9].split(":"), strict=True)
    ]
    zero = [*cols[:2], "ZERO", *cols[3:9], ":".join(flat)]
    one = [*cols[:2], "LEFT", *cols[3:9], ":".join(left)]
    path = tmp_path / "three.vcf"
    path.write_text(
        _header(DELLY) + "".join(f"{chr(9).join(r)}\n" for r in (cols, zero, one))
    )
    hand, flanks_zero, left_zero = read_deletions(path, caller="delly")
    assert hand.depth_ratio == pytest.approx(94 / (565 + 475))
    assert hand.callers == {"delly": "DEL00000002 QUAL 960"}
    assert flanks_zero.depth_ratio is None
    assert left_zero.depth_ratio is None


def test_deletions_only_measures_from_ref(tmp_path: Path) -> None:
    """Dysgu's SVLEN is one more than its REF says (plan Task 9): a 49 bp deletion with
    `SVLEN=50` is cut, a 50 bp one kept."""
    cols, info, _ = next(
        r
        for r in _records(DYSGU)
        if r[1].get("SVTYPE") == "DEL"
        and r[0][6] == "PASS"
        and r[2]["GT"] in CARRIED
        and not r[0][4].startswith("<")
    )
    rows = []
    for name, deleted in (("cut", 49), ("kept", 50)):
        ref = "A" * (deleted + 1)
        fields = {
            **info,
            "END": str(int(cols[1]) + deleted + 1),
            "SVLEN": str(deleted + 1),
        }
        text = ";".join(f"{k}={v}" if v else k for k, v in fields.items())
        rows.append(
            [cols[0], cols[1], name, ref, "A", cols[5], "PASS", text, *cols[8:]]
        )
    path = tmp_path / "two.vcf"
    path.write_text(_header(DYSGU) + "".join(f"{chr(9).join(r)}\n" for r in rows))
    out = deletions_only(path, tmp_path / "out.vcf.gz", caller="dysgu", min_length=50)
    assert [cols[2] for cols, _, _ in _records(out)] == ["kept"]


def test_read_delly_hg002(tmp_path: Path) -> None:
    out = deletions_only(
        DELLY, tmp_path / "delly.vcf.gz", caller="delly", min_length=25
    )
    found = read_deletions(out, caller="delly")
    expected = _kept(DELLY, "delly")
    assert 0 < expected < sum(1 for _ in _records(DELLY))
    assert len(found) == expected
    assert all(d.end > d.pos and d.gt in {"0/1", "1/1"} for d in found)
    assert all(set(d.callers) == {"delly"} for d in found)


def test_read_dysgu_hg002(tmp_path: Path) -> None:
    out = deletions_only(
        DYSGU, tmp_path / "dysgu.vcf.gz", caller="dysgu", min_length=50
    )
    found = read_deletions(out, caller="dysgu")
    expected = _kept(DYSGU, "dysgu")
    assert 0 < expected < sum(1 for _ in _records(DYSGU))
    assert len(found) == expected
    assert all(d.end > d.pos and d.gt in {"0/1", "1/1"} for d in found)
    assert all(d.depth_ratio is None for d in found)
    # A sequence-resolved record: its length from its REF, though Dysgu's END and SVLEN
    # are one more (plan Task 9).
    cols, info, _ = next(
        r for r in _records(out) if not r[0][4].startswith("<") and len(r[0][3]) > 50
    )
    (d,) = (d for d in found if d.chrom == cols[0] and d.pos == int(cols[1]))
    assert d.length == len(cols[3]) - 1
    assert int(info["END"]) == d.end + 1


def test_read_deletions_truth_shape(tmp_path: Path) -> None:
    deleted = "ACGT" * 300
    path = tmp_path / "truth.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tsnv\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        f"chr1\t200\tins\tA\tA{deleted}\t50\tPASS\t.\tGT\t0/1\n"
        f"chr1\t1000\tdel\tT{deleted}\tT\t50\tPASS\t.\tGT\t1/1\n"
    )
    (d,) = read_deletions(path, caller="truth")
    assert (d.chrom, d.pos, d.end, d.ref, d.gt) == (
        "chr1",
        1000,
        1000 + 1201 - 1,
        "T",
        "1/1",
    )
    assert d.length == 1200
    assert d.callers == {"truth": "del QUAL 50"}


@pytest.mark.parametrize(("caller", "fixture"), [("delly", DELLY), ("dysgu", DYSGU)])
def test_deletions_only_filters(tmp_path: Path, caller: str, fixture: Path) -> None:
    out = deletions_only(
        fixture, tmp_path / f"{caller}.vcf.gz", caller=caller,
        min_length=MIN_LENGTH[caller],
    )  # fmt: skip
    assert Path(f"{out}.tbi").is_file()
    rows = list(_records(out))
    assert rows
    for cols, info, values in rows:
        assert _length(cols, info) >= MIN_LENGTH[caller]
        assert info["SVTYPE"] == "DEL"
        assert cols[6] == "PASS"
        assert values.get("FT") != "LowQual"
        assert values["GT"] not in {"0/0", "./.", "."}
    positions = [(cols[0], int(cols[1])) for cols, _, _ in rows]
    assert positions == sorted(positions)


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


TEMPLATE = "chr1\t0\t10000\ttelomere\nchrM\nM\nchrMT\n"
# chrY absent: its shard is skipped. chr20 after chr2 by the `.fai`, not by name.
CONTIGS = ["chr1", "chr2", "chr20", "chrX", "chrM", "chrUn_KI270302v1", "chrEBV"]


class _Delly:
    """Stands for `run_container` and `run_tool`: copies the template, records each shard's
    arguments and exclude list, and fails the shard named `fail`."""

    def __init__(self, fail: str | None = None) -> None:
        self.fail = fail
        self.lock = threading.Lock()
        self.template: list[
            tuple[list[str], list[tuple[Path, str, bool]], str | None]
        ] = []
        self.shards: dict[
            str, tuple[list[str], list[tuple[Path, str, bool]], str | None, str]
        ] = {}
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
            self.template.append((list(args), list(mounts), label))
            (out / "exclude.template.tsv").write_text(TEMPLATE)
            return
        if out.name == self.fail:
            raise RuntimeError(f"shard {out.name} failed")
        with self.lock:
            self.shards[out.name] = (
                list(args),
                list(mounts),
                label,
                (out / "exclude.tsv").read_text(),
            )
        (out / "delly.bcf").write_bytes(b"")

    def tool(self, args: list[str], **_: object) -> str:
        self.tools.append(list(args))
        return ""


def test_run_delly_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    fake = _Delly()
    monkeypatch.setattr(callers, "run_container", fake.container)
    monkeypatch.setattr(callers, "run_tool", fake.tool)
    assert (
        run_delly(cram, fasta, out, threads=2, shards=PRIMARY) == out / "delly.vcf.gz"
    )
    ((args, mounts, label),) = fake.template
    assert args == ["cp", "/opt/delly/human.hg38.excl.tsv", "/out/exclude.template.tsv"]
    assert mounts == [(out, "/out", False)] and label is None
    shards = ["chr1", "chr2", "chr20", "chrX"]
    assert sorted(fake.shards) == sorted(shards)
    for chrom, (args, mounts, label, exclude) in fake.shards.items():
        assert label == "dna:delly"
        assert args[:2] == ["delly", "sr"]
        assert args[args.index("-t") + 1] == "DEL"
        assert args[args.index("-g") + 1] == "/ref/g.fa"
        assert args[args.index("-x") + 1] == "/out/exclude.tsv"
        assert args[args.index("-o") + 1] == "/out/delly.bcf"
        assert args[-1] == "/cram/x.cram"
        assert not {"-p", "-@", "--threads", "--procs"} & set(args)
        assert {(t, ro) for _, t, ro in mounts} == {
            ("/cram", True),
            ("/ref", True),
            ("/out", False),
        }
        assert (out / "delly-shards" / chrom, "/out", False) in mounts
        assert exclude.startswith(TEMPLATE)
        bare = exclude.removeprefix(TEMPLATE).splitlines()
        assert bare == [c for c in CONTIGS if c != chrom]
        assert chrom not in exclude.splitlines()
    concat, index = fake.tools
    bcfs = [str(out / "delly-shards" / c / "delly.bcf") for c in shards]
    assert concat[:2] == ["bcftools", "concat"]
    assert concat[concat.index("-o") + 1] == str(out / "delly.vcf.gz")
    assert concat[-len(bcfs) :] == bcfs
    assert index[:2] == ["bcftools", "index"] and index[-1] == str(out / "delly.vcf.gz")


def test_run_delly_shard_fails_the_step(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    fake = _Delly(fail="chr2")
    monkeypatch.setattr(callers, "run_container", fake.container)
    monkeypatch.setattr(callers, "run_tool", fake.tool)
    with pytest.raises(RuntimeError, match="shard chr2 failed"):
        run_delly(cram, fasta, out, threads=1, shards=PRIMARY)
    assert fake.tools == []
    assert sorted(fake.shards) == ["chr1"]  # chr20 and chrX, behind it, never started


def test_run_dysgu_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    calls: list[tuple[list[str], list[tuple[Path, str, bool]], str | None, str]] = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert name == "dysgu"
        (shard,) = (host for host, target, _ in mounts if target == "/out")
        calls.append(
            (list(args), list(mounts), label, (shard / "search.bed").read_text())
        )
        # The fixture's records on the shard's chromosome, then one on another (as a
        # translocation's may be), which `bcftools concat` would refuse as a second block.
        lines = DYSGU.read_text().splitlines(keepends=True)
        records = [line.split("\t", 1)[1] for line in lines if not line.startswith("#")]
        other = "chr2" if shard.name == "chr1" else "chr1"
        (shard / "dysgu.vcf").write_text(
            "".join(line for line in lines if line.startswith("#"))
            + "".join(f"{shard.name}\t{rest}" for rest in records)
            + f"{other}\t{records[0]}"
        )

    monkeypatch.setattr(callers, "run_container", fake)
    result = run_dysgu(cram, fasta, out, threads=3, shards=PRIMARY)
    shards = ["chr1", "chr2", "chr20", "chrX"]
    assert [m for _, mounts, _, _ in calls for m in mounts if m[1] == "/out"] == [
        (out / "dysgu-shards" / c, "/out", False) for c in shards
    ]
    lengths = dict(
        line.split("\t")[:2] for line in Path(f"{fasta}.fai").read_text().splitlines()
    )
    for chrom, (args, mounts, label, bed) in zip(shards, calls, strict=True):
        assert label == "dna:dysgu"
        assert args[:2] == ["dysgu", "run"]
        assert args[args.index("-p") + 1] == "3"
        assert args[args.index("--search") + 1] == "/out/search.bed"
        assert bed == f"{chrom}\t0\t{lengths[chrom]}\n"
        assert args[args.index("-o") + 1] == "/out/dysgu.vcf"
        reference, work, reads = args[-3:]
        assert (reference, reads) == ("/ref/g.fa", "/cram/x.cram")
        assert work.startswith("/out/") and work != "/out/"
        assert {(t, ro) for _, t, ro in mounts} == {
            ("/cram", True),
            ("/ref", True),
            ("/out", False),
        }
    assert result == out / "dysgu.vcf.gz" and Path(f"{result}.tbi").is_file()
    assert not (out / "dysgu.vcf").exists()
    chroms = [cols[0] for cols, _, _ in _records(result)]
    assert chroms == sorted(chroms, key=shards.index)
    assert len(chroms) == len(shards) * sum(1 for _ in _records(DYSGU))


def test_run_dysgu_shard_fails_the_step(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path, CONTIGS)
    tools: list[list[str]] = []

    def fake(name: str, args: list[str], **_: object) -> None:
        raise RuntimeError("shard failed")

    monkeypatch.setattr(callers, "run_container", fake)
    monkeypatch.setattr(callers, "run_tool", lambda args, **_: tools.append(args))
    with pytest.raises(RuntimeError, match="shard failed"):
        run_dysgu(cram, fasta, out, threads=1, shards=PRIMARY)
    assert tools == [] and not (out / "dysgu.vcf.gz").exists()


# --- insertions and duplications (slice 3d) -------------------------------------------------


def test_read_insertions_and_dups(tmp_path: Path) -> None:
    """A sequence-resolved insertion's length from its ALT, a symbolic one's from `SVLEN`
    (absolute); a duplication's end from `END`, sequence-resolved or not (Dysgu writes
    both); a lower bound from Dysgu's `LPREC=0` or INSurVeyor's `INCOMPLETE_ASSEMBLY`."""
    inserted = "ACGT" * 75
    path = tmp_path / "sv.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        f"chr1\t100\tseq\tA\tA{inserted}\t20\tPASS\tSVTYPE=INS;END=100;LPREC=1\tGT\t0/1\n"
        "chr1\t200\tsym\tC\t<INS>\t30\tPASS\tSVTYPE=INS;END=200;SVLEN=-1200;LPREC=0"
        "\tGT\t1/1\n"
        "chr1\t300\tpart\tC\t<INS>\t.\tPASS\tEND=300;SVTYPE=INS;SVINSSEQ="
        f"{inserted[:90]};INCOMPLETE_ASSEMBLY\tGT\t1\n"
        "chr1\t1000\tdup\tG\t<DUP>\t40\tPASS\tSVTYPE=DUP;END=6000;SVLEN=5000\tGT\t0/1\n"
        f"chr1\t7000\tseqdup\tT\tT{inserted}\t40\tPASS\tSVTYPE=DUP;END=7300\tGT\t0/1\n"
    )
    seq, sym, part = read_deletions(path, caller="dysgu", svtype="INS")
    assert (seq.svtype, seq.pos, seq.end, seq.svlen, seq.length) == (
        "INS",
        100,
        101,
        300,
        300,
    )
    assert (sym.end, sym.svlen, sym.gt, sym.lower_bound) == (201, 1200, "1/1", True)
    assert seq.lower_bound is False
    assert (part.svlen, part.lower_bound) == (90, True)
    assert seq.callers == {"dysgu": "seq QUAL 20"}
    dup, seqdup = read_deletions(path, caller="dysgu", svtype="DUP")
    assert (dup.svtype, dup.pos, dup.end, dup.svlen, dup.length) == (
        "DUP",
        1000,
        6000,
        None,
        5000,
    )
    assert (seqdup.end, seqdup.length) == (7300, 300)
    assert read_deletions(path, caller="dysgu") == []
    assert [d.svtype for d in read_deletions(path, caller="dysgu", svtype="DUP")] == [
        "DUP",
        "DUP",
    ]
    passed = read_records(DYSGU_DUP, caller="dysgu", svtype="DUP", pass_only=True)
    assert len(passed) == sum(
        1 for cols, _, _ in _records(DYSGU_DUP) if cols[6] == "PASS"
    )


def _ins_length(cols: list[str], info: dict[str, str]) -> int:
    """INSurVeyor's incomplete assemblies carry no `SVLEN`: their assembled part's length."""
    if not cols[4].startswith("<"):
        return len(cols[4]) - len(cols[3])
    if "SVLEN" in info:
        return abs(int(info["SVLEN"]))
    return len(info["SVINSSEQ"])


def test_insurveyor_insertion_of_unknown_length_never_kept(tmp_path: Path) -> None:
    """An incomplete INSurVeyor assembly is kept only with an assembled part of
    `min_length` or more, its length's lower bound: one with neither `SVLEN` nor a part is
    cut here, so the screen never sees an insertion of unknown length from a caller."""
    path = tmp_path / "insurveyor.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="l">\n'
        '##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="s">\n'
        '##INFO=<ID=INCOMPLETE_ASSEMBLY,Number=0,Type=Flag,Description="i">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + "".join(
            f"chr1\t{pos}\t{name}\tC\t<INS>\t.\tPASS\tEND={pos};SVTYPE=INS{extra}"
            ";INCOMPLETE_ASSEMBLY\tGT\t1\n"
            for pos, name, extra in (
                (100, "none", ""),
                (200, "dot", ";SVINSSEQ=."),
                (300, "part", ";SVINSSEQ=" + "ACGT" * 15),
            )
        )
    )
    out = deletions_only(
        path, tmp_path / "out.vcf.gz", caller="insurveyor", min_length=50, svtype="INS"
    )
    (kept,) = read_deletions(out, caller="insurveyor", svtype="INS")
    assert (kept.pos, kept.svlen, kept.lower_bound) == (300, 60, True)


@pytest.mark.parametrize(
    ("caller", "svtype", "fixture"),
    [
        ("dysgu", "INS", DYSGU_INS),
        ("insurveyor", "INS", INSURVEYOR),
        ("dysgu", "DUP", DYSGU_DUP),
    ],
)
def test_deletions_only_ins_and_dup(
    tmp_path: Path, caller: str, svtype: str, fixture: Path
) -> None:
    """`_KEEP` per (caller, type): FILTER PASS, a non-reference GT, the type's length as
    `read_records` measures it; counted by the test's own reading."""
    out = deletions_only(
        fixture, tmp_path / "out.vcf.gz", caller=caller, min_length=50, svtype=svtype
    )
    expected = sum(
        1
        for cols, info, values in _records(fixture)
        if info.get("SVTYPE") == svtype
        and cols[6] == "PASS"
        and values["GT"] in CARRIED | {"1"}
        and (
            _ins_length(cols, info)
            if svtype == "INS"
            else int(info["END"]) - int(cols[1])
        )
        >= 50
    )
    found = read_deletions(out, caller=caller, svtype=svtype)
    assert 0 < expected == len(found)
    assert all(d.svtype == svtype and d.length >= 50 for d in found)
    if caller == "insurveyor":
        assert sum(d.lower_bound for d in found) == sum(
            1 for _, info, _ in _records(fixture) if "INCOMPLETE_ASSEMBLY" in info
        )
    with pytest.raises(ValueError, match="no INS filter for caller 'delly'"):
        filter_expression("delly", 50, svtype="INS")
