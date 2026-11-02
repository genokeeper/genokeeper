"""The paralog loci's callers (`genokeeper.paralogs`): Parascopy's six loci, Gauchian, KILDA
and Aldy CYP2A6 read from their HG002 outputs (tests/fixtures/hg002, Task 6), and their
container runs, whose container (or Docker itself, for the offline flag) is a stub that
records the arguments and writes the output.

The values are the tests' own reading of the fixtures (`_rows`), never the code under
test."""

import gzip
import shutil
from collections.abc import Callable
from pathlib import Path

import pytest

from genokeeper import containers, paralogs, smn, tools
from genokeeper.deletions.exons import load_spans, locus_span
from genokeeper.paralogs import (
    GAUCHIAN_LABEL,
    KILDA_LABEL,
    LOCI,
    PARALOGS_LABEL,
    GauchianCall,
    PsCall,
    aldy_copies,
    depth_ratio,
    read_gauchian,
    read_kilda,
    read_locus,
    run_aldy_cyp2a6,
    run_gauchian,
    run_kilda,
    unique_part,
    wild_type,
)
from genokeeper.pgx import cyp2d6
from genokeeper.smn import parascopy_cn, parascopy_depth

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
LOCI_BED = FIXTURES / "hg002.parascopy.loci.bed.txt"
GAUCHIAN_TSV = FIXTURES / "hg002.gauchian.tsv.txt"
GAUCHIAN_JSON = FIXTURES / "hg002.gauchian.json.txt"
KILDA = FIXTURES / "hg002.kilda.txt"
ALDY = FIXTURES / "hg002.aldy.CYP2A6.txt"
MODELS = ("STRC", "RHCE", "AMY1C", "C4A", "CYP2A7", "GBA")
HEADER = (
    "#chrom\tstart\tend\tsample\tcn_filter\tcn\tcn_qual\tagcn_filter\tagcn\tagcn_qual"
    "\tmain_region\n"
)
type Mount = tuple[Path, str, bool]
type Seen = list[tuple[str, list[str], list[Mount], str | None]]


def _rows(chrom: str, start: int, end: int) -> list[list[str]]:
    """The fixture's rows overlapping 0-based `[start, end)`, read by hand."""
    lines = LOCI_BED.read_text().splitlines()
    rows = [line.split("\t") for line in lines if not line.startswith("#")]
    return [r for r in rows if r[0] == chrom and int(r[1]) < end and int(r[2]) > start]


def _span(locus: str) -> tuple[str, int, int]:
    return locus_span(load_spans(), LOCI[locus].genes)


def test_loci_six() -> None:
    assert set(LOCI) == {"STRC", "RHD", "AMY1", "C4", "CYP2A6", "GBA1"}
    assert tuple(locus.model for locus in LOCI.values()) == MODELS
    assert {k for k, v in LOCI.items() if v.role == "carrier"} == {"STRC", "GBA1"}
    assert {k for k, v in LOCI.items() if v.measure == "agcn"} == {"AMY1", "C4", "GBA1"}
    assert LOCI["STRC"].unique == "CATSPER2"
    assert all(v.unique is None for k, v in LOCI.items() if k != "STRC")


def test_read_loci_hg002() -> None:
    calls = {
        name: read_locus(LOCI_BED, locus, _span(name)) for name, locus in LOCI.items()
    }
    assert all(c.rows >= 1 for c in calls.values())
    # By hand: psCN (column 6) for STRC and RHD, agCN (column 9) for the C4 total.
    for name, column in (("STRC", 5), ("RHD", 5), ("C4", 8)):
        rows = _rows(*_span(name))
        (value,) = {r[column] for r in rows}
        filters = {r[column - 1] for r in rows}
        assert filters == {"PASS"}
        assert calls[name] == PsCall(
            cn=int(value),
            qual=min(float(r[column + 1]) for r in rows),
            passed=True,
            rows=len(rows),
        )
    assert (calls["STRC"].cn, calls["RHD"].cn, calls["C4"].cn) == (2, 2, 4)
    # HG002's truth (resources/paralogs/hg002_truth.tsv): AMY1 6, CYP2A6 2; AMY1's psCN row
    # fails its filter, its agCN row passes.
    assert (calls["AMY1"].cn, calls["CYP2A6"].cn) == (6, 2)
    # GBA1 over its exons 9-11 (`GBA1_ex9-11`), where Gauchian counts GBA+GBAP1: agCN 4.
    assert LOCI["GBA1"].genes == ("GBA1_ex9-11",)
    (gba,) = _rows(*load_spans()["GBA1_ex9-11"])
    assert (gba[7], gba[8]) == ("PASS", "4")
    assert calls["GBA1"] == PsCall(4, float(gba[9]), True, 1)
    (amy1,) = _rows(*_span("AMY1"))
    assert amy1[4] != "PASS" and amy1[7] == "PASS"


def test_gba1_whole_span_disagrees() -> None:
    """Over GBA1's whole MANE span the rows disagree, 4 over the GBAP1-homologous part
    and 2 over the rest, which has no GBAP1 copy: why the locus reads only exons 9-11."""
    whole = load_spans()["GBA1"]
    assert sorted(r[8] for r in _rows(*whole)) == ["2", "4"]
    assert read_locus(LOCI_BED, LOCI["GBA1"], whole) == PsCall(None, 40.0, True, 2)


def test_read_locus_gzip(tmp_path: Path) -> None:
    path = tmp_path / "res.paralog.bed.gz"
    path.write_bytes(gzip.compress(LOCI_BED.read_bytes()))
    span = _span("STRC")
    assert read_locus(path, LOCI["STRC"], span) == read_locus(
        LOCI_BED, LOCI["STRC"], span
    )


def test_locus_spans() -> None:
    spans = load_spans()
    assert _span("AMY1") == ("chr1", spans["AMY1A"][1], spans["AMY1C"][2])
    chrom, start, end = _span("C4")
    for gene in ("C4A", "C4B"):
        assert spans[gene][0] == chrom
        assert start <= spans[gene][1] and spans[gene][2] <= end
    with pytest.raises(ValueError, match="chromosomes"):
        locus_span(spans, ("STRC", "RHD"))


def _bed(tmp_path: Path, rows: list[str]) -> Path:
    path = tmp_path / "res.paralog.bed"
    path.write_text("## parascopy\n" + HEADER + "".join(f"{r}\n" for r in rows))
    return path


def test_rows_disagree_none(tmp_path: Path) -> None:
    strc = LOCI["STRC"]
    span = ("chr15", 1000, 3000)
    one, two = (
        "chr15\t900\t2000\tS\tPASS\t1\t40\tPASS\t3\t30\tchr15:1-2",
        "chr15\t2000\t3100\tS\tPASS\t2\t50\tPASS\t4\t30\tchr15:1-2",
    )
    assert read_locus(_bed(tmp_path, [one, two]), strc, span) == PsCall(
        None, 40.0, True, 2
    )
    failed = one.replace("\tPASS\t1\t", "\tLowQual\t1\t")
    assert read_locus(_bed(tmp_path, [failed]), strc, span) == PsCall(
        None, 40.0, False, 1
    )
    unknown = one.replace("\tPASS\t1\t40\t", "\tPASS\t?\t?\t")
    assert read_locus(_bed(tmp_path, [unknown]), strc, span) == PsCall(
        None, None, True, 1
    )
    far = "chr15\t5000\t6000\tS\tPASS\t2\t50\tPASS\t4\t30\tchr15:1-2"
    other = "chr1\t1000\t3000\tS\tPASS\t2\t50\tPASS\t4\t30\tchr1:1-2"
    assert read_locus(_bed(tmp_path, [far, other]), strc, span) == PsCall(
        None, None, False, 0
    )
    assert read_locus(_bed(tmp_path, [one]), strc, span) == PsCall(1, 40.0, True, 1)


def _model(tmp_path: Path) -> Path:
    """The STRC model's lines this test needs: its main region and the two duplications
    over CATSPER2 that the HG002 table's rows name (`main_region`)."""
    model = tmp_path / "STRC.gz"
    with gzip.open(model, "wt") as f:
        f.write("main\tchr15:43558001-43750000\t*\tname=STRC\n")
        f.write(
            "duplication\tchr15:43630563-43631371\tchr15:43727043-43727850:+\tix=0\n"
        )
        f.write(
            "duplication\tchr15:43644125-43645617\tchr15:43740914-43742400:+\tix=1\n"
        )
    return model


def test_unique_part_strc(tmp_path: Path) -> None:
    catsper2 = load_spans()["CATSPER2"]
    parts = unique_part(catsper2, _model(tmp_path))
    assert parts
    chrom, start, end = catsper2
    assert all(c == chrom and start <= s < e <= end for c, s, e in parts)
    assert sum(e - s for _, s, e in parts) == (end - start) - (43631371 - 43630562) - (
        43645617 - 43644124
    )
    assert ("chr15", 43630562, 43631371) not in parts


def test_depth_ratio(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: list[tuple[list[str], str]] = []

    def bedcov(args: list[str], **_: object) -> str:
        bed = Path(args[-2])
        seen.append((list(args), bed.read_text()))
        return "chr15\t100\t200\t2000\nchr15\t300\t400\t1000\n"

    monkeypatch.setattr(paralogs, "run_tool", bedcov)
    cram, fasta = tmp_path / "x.cram", tmp_path / "g.fa"
    parts = [("chr15", 100, 200), ("chr15", 300, 400)]
    assert depth_ratio(cram, fasta, parts, 30.0) == pytest.approx(3000 / 200 / 30)
    ((args, bed),) = seen
    assert args[:2] == ["samtools", "bedcov"]
    assert args[args.index("--reference") + 1] == str(fasta)
    assert args[-1] == str(cram)
    assert bed == "chr15\t100\t200\nchr15\t300\t400\n"
    assert depth_ratio(cram, fasta, [], 30.0) is None
    assert len(seen) == 1


def test_read_gauchian_hg002() -> None:
    call = read_gauchian(GAUCHIAN_TSV, GAUCHIAN_JSON)
    header, row = (line.split("\t") for line in GAUCHIAN_TSV.read_text().splitlines())
    cell = dict(zip(header, row, strict=True))
    assert call.total_cn == int(cell["CN(GBA+GBAP1)"]) == 4
    assert call.carrier is (cell["is_carrier(GBAP1-like_variant_exon9-11)"] == "True")
    assert call.biallelic is (
        cell["is_biallelic(GBAP1-like_variant_exon9-11)"] == "True"
    )
    assert call == GauchianCall(
        total_cn=4, carrier=False, biallelic=False, variants=(), other=(), reads=None
    )
    # Task 6 recorded HG002 as wild type.
    assert wild_type(call)


def test_read_gauchian_carrier(tmp_path: Path) -> None:
    """A called variant: its name in both files, its alt count from `variant_raw_count`;
    `None` and `?` cells read as None, never wild type."""
    tsv, js = tmp_path / "g.tsv", tmp_path / "g.json"
    header = GAUCHIAN_TSV.read_text().splitlines()[0]
    tsv.write_text(
        f"{header}\nsample\tFalse\tTrue\t4\tN/A\tc.1448T>G(p.Leu483Arg)\tNone\n"
    )
    text = GAUCHIAN_JSON.read_text().replace(
        '"GBAP1_like_variant_exon9_11": ""',
        '"GBAP1_like_variant_exon9_11": "c.1448T>G(p.Leu483Arg)"',
    )
    js.write_text(
        text.replace(
            '"c.1448T>G(p.Leu483Arg)": "0,67"', '"c.1448T>G(p.Leu483Arg)": "7,30"'
        )
    )
    call = read_gauchian(tsv, js)
    assert call == GauchianCall(4, True, False, ("c.1448T>G(p.Leu483Arg)",), (), 7)
    assert not wild_type(call)
    tsv.write_text(f"{header}\nsample\t?\tNone\tNone\tN/A\tNone\tc.1A>T\n")
    unknown = read_gauchian(tsv, GAUCHIAN_JSON)
    assert unknown == GauchianCall(None, None, None, (), ("c.1A>T",), None)
    assert not wild_type(unknown)
    # A reported breakpoint (a fusion) beside four copies is no wild type.
    tsv.write_text(f"{header}\nsample\tFalse\tFalse\t4\tintron9\tNone\tNone\n")
    fusion = read_gauchian(tsv, GAUCHIAN_JSON)
    assert fusion.breakpoint == "intron9" and not wild_type(fusion)


def test_read_kilda_hg002(tmp_path: Path) -> None:
    value = read_kilda(KILDA)
    assert isinstance(value, float)
    assert value == pytest.approx(34.47)
    na = tmp_path / "kilda_kiv2.CNs"
    na.write_text("ID\tKIV2_CN\tquantile\nsample\tNA\tNA\n")
    assert read_kilda(na) is None


def test_aldy_copies() -> None:
    assert aldy_copies("*1/*4") == 1
    assert aldy_copies("*1x2/*1") == 3
    assert aldy_copies("*1/*1") == 2
    assert aldy_copies("*4/*4") == 0
    assert aldy_copies("*1.009/*12+*1") == 3
    assert aldy_copies("garbage") is None
    assert aldy_copies("*1/*1 or *1/*4") is None


def test_run_aldy_cyp2a6(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: Seen = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[Mount],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), list(mounts), label))
        (out,) = (host for host, _, ro in mounts if not ro)
        shutil.copy(ALDY, out / "CYP2A6.aldy")

    monkeypatch.setattr(cyp2d6, "run_container", fake)
    cram, fasta, out = _inputs(tmp_path)
    call = run_aldy_cyp2a6(cram, fasta, out)
    ((name, args, _, label),) = seen
    assert (name, label) == ("aldy", PARALOGS_LABEL)
    assert PARALOGS_LABEL == "dna:paralogs"
    assert args[args.index("-g") + 1] == "CYP2A6"
    assert call.diplotype == "*1/*1"
    assert aldy_copies(call.diplotype) == 2


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    cram, fasta, out = (
        tmp_path / "in" / "x.cram",
        tmp_path / "ref" / "g.fa",
        tmp_path / "out",
    )
    for path in (cram, fasta):
        path.parent.mkdir()
        path.write_bytes(b"")
    out.mkdir()
    return cram, fasta, out


def _lab(tmp_path: Path, models: tuple[str, ...] = ("SMN1", *MODELS)) -> Path:
    lab = tmp_path / "lab"
    root = lab / "parascopy"
    paths = [
        root / "homology_table" / "GRCh38.bed.gz",
        *(root / "models_GRCh38_1KGP" / "EUR" / f"{m}.gz" for m in models),
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"")
    return lab


def _recorder(seen: Seen) -> Callable[..., None]:
    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[Mount],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), list(mounts), label))

    return fake


def test_parascopy_cn_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """One `cn-using` with the six models, on a depth file another step wrote once per
    CRAM (its folder read-only), `/parascopy` read-only, behind the label given."""
    seen: Seen = []
    monkeypatch.setattr(smn, "run_container", _recorder(seen))
    cram, fasta, out = _inputs(tmp_path)
    shared = tmp_path / "shared"
    shared.mkdir()
    depth = parascopy_depth(cram, fasta, shared, threads=2)
    assert depth == shared / "depth" / "depth.csv"
    lab = _lab(tmp_path)
    models = [locus.model for locus in LOCI.values()]
    result = parascopy_cn(
        cram, fasta, depth, out, lab=lab, threads=2, models=models, label=PARALOGS_LABEL
    )
    assert result == out / "cn" / "res.paralog.bed.gz"
    assert len(seen) == 2
    name, args, mounts, label = seen[1]
    assert (name, label) == ("parascopy", "dna:paralogs")
    assert args[:2] == ["parascopy", "cn-using"]
    assert args[2 : 2 + len(MODELS)] == [
        f"/parascopy/models_GRCh38_1KGP/EUR/{m}.gz" for m in MODELS
    ]
    assert args.count("cn-using") == 1
    assert args[args.index("-d") + 1] == "/depth/depth.csv"
    assert (depth.parent, "/depth", True) in mounts
    assert (lab / "parascopy", "/parascopy", True) in mounts
    assert all(ro for _, target, ro in mounts if target != "/out")
    missing = _lab(tmp_path / "other", ("SMN1",))
    with pytest.raises(FileNotFoundError, match="bootstrap"):
        parascopy_cn(
            cram, fasta, depth, out, lab=missing, threads=2, models=models,
            label=PARALOGS_LABEL,
        )  # fmt: skip


class _Docker:
    """Stands for Docker under `run_container`: records each full command, the gate's label
    included (a gate path that exists), and writes the tool's output into `/out`."""

    def __init__(self, write: Callable[[Path], None]) -> None:
        self.commands: list[list[str]] = []
        self.write = write

    def __call__(self, args: list[str], **_: object) -> str:
        self.commands.append(list(args))
        out = next(
            Path(a.split("source=")[1].split(",")[0])
            for a in args
            if a.endswith(",target=/out")
        )
        self.write(out)
        return ""


def _docker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, write: Callable[[Path], None]
) -> _Docker:
    gate = tmp_path / "gate"
    gate.write_text("")
    docker = _Docker(write)
    monkeypatch.setattr(tools, "GATE", gate)
    monkeypatch.setattr(containers, "_docker_info", lambda: None)
    monkeypatch.setattr(containers, "run_tool", docker)
    return docker


def _offline(command: list[str], label: str) -> list[str]:
    """The command after the gate and its label; `--network none` before the image."""
    assert command[1:3] == [label, "--"]
    run = command[3:]
    assert run[:2] == ["docker", "run"]
    assert run[run.index("--network") + 1] == "none"
    return run


def test_run_gauchian_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def write(out: Path) -> None:
        shutil.copy(GAUCHIAN_TSV, out / "gauchian.tsv")
        shutil.copy(GAUCHIAN_JSON, out / "gauchian.json")

    docker = _docker(tmp_path, monkeypatch, write)
    cram, fasta, out = _inputs(tmp_path)
    assert run_gauchian(cram, fasta, out, threads=3) == out / "gauchian.tsv"
    ((command,),) = [docker.commands]
    run = _offline(command, GAUCHIAN_LABEL)
    assert GAUCHIAN_LABEL == "dna:gauchian"
    image = run.index(containers.IMAGES["gauchian"])
    args = run[image + 1 :]
    assert args[0] == "gauchian"
    assert args[args.index("--manifest") + 1] == "/out/manifest.txt"
    assert (out / "manifest.txt").read_text() == "/cram/x.cram\n"
    assert args[args.index("--genome") + 1] == "38"
    assert args[args.index("--reference") + 1] == "/ref/g.fa"
    assert args[args.index("--prefix") + 1] == "gauchian"
    assert args[args.index("--outDir") + 1] == "/out"
    assert args[args.index("--threads") + 1] == "3"
    assert (
        f"type=bind,source={(tmp_path / 'in').resolve()},target=/cram,readonly" in run
    )


def test_run_kilda_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def write(out: Path) -> None:
        (out / "kilda-output").mkdir()
        shutil.copy(KILDA, out / "kilda-output" / "kilda_kiv2.CNs")

    docker = _docker(tmp_path, monkeypatch, write)
    cram, fasta, out = _inputs(tmp_path)
    lpa = load_spans()["LPA"]
    result = run_kilda(cram, fasta, out, threads=4, span=lpa)
    assert result == out / "kilda-output" / "kilda_kiv2.CNs"
    assert read_kilda(result) == pytest.approx(34.47)
    ((command,),) = [docker.commands]
    run = _offline(command, KILDA_LABEL)
    assert KILDA_LABEL == "dna:kilda"
    image = run.index(containers.IMAGES["kilda"])
    assert run[image + 1 : image + 3] == ["bash", "-c"]
    steps = [s.strip() for s in run[image + 3].split("&&")]
    region = f"chr6:{lpa[1] - 5000 + 1}-{lpa[2] + 5000}"
    assert (
        f"samtools view -b -T /ref/g.fa -o /out/lpa.bam /cram/x.cram {region}" in steps
    )
    assert (
        "samtools view -b -f 4 -T /ref/g.fa -o /out/unmapped.bam /cram/x.cram '*'"
        in (steps)
    )
    assert "samtools cat -o /out/both.bam /out/lpa.bam /out/unmapped.bam" in steps
    (count,) = [s for s in steps if s.startswith("jellyfish count")]
    assert "-m 31" in count and "-t 4" in count
    (kilda,) = [s for s in steps if "kilda.py" in s]
    assert "-c /out/counts_list.tsv" in kilda and "-o /out/kilda-output" in kilda
    assert "-k KIV2_hg38_kmers_6copies_specific.tsv" in kilda
    assert "cd /opt/kilda/data/kmers_GRCh38_LPA_k31" in steps
    assert (out / "counts_list.tsv").read_text() == "sample\t/out/s.counts.tsv\n"
    assert steps[-1].startswith("rm -f /out/lpa.bam")
    for mount in ("target=/cram,readonly", "target=/ref,readonly"):
        assert any(a.endswith(mount) for a in run)
