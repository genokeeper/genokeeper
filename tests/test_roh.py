"""Runs of homozygosity (slice 3f, "ROH"): our genotypes and the truth's at the AF file's
alleles (ruling R7), `bcftools roh`'s runs, their summary and the `roh:summary` statement.

`fixtures/hg002/hg002.roh.txt` is `bcftools roh -G30 --AF-file af.EUR -O r` on HG002's truth
genotypes at the EUR AF sites inside the benchmark (Task 5): 2 runs of 5 Mb or more, none of
10 Mb."""

import gzip
import json
import random
import subprocess
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path

import pysam
import pytest
from fixtures.gvcf_mini import block, ref, write_fasta, write_gvcf
from fixtures.run_context import run_context

from genokeeper import roh
from genokeeper.config import RohThresholds
from genokeeper.genotype.gvcf import NOCALL_REASON, Called
from genokeeper.prs.pgsc import write_target_vcf
from genokeeper.roh import Run, read_runs, summarise
from genokeeper.sources.roh import ROH_TOOL, RohSource
from genokeeper.store.db import Level, RowType

type Allele = tuple[str, int, str, str]

HG002 = Path(__file__).parent / "fixtures" / "hg002" / "hg002.roh.txt"
GRADING = RohThresholds()
MB = 1_000_000


def test_read_runs_hg002() -> None:
    runs = read_runs(HG002)
    assert len(runs) == 5811
    assert runs[0] == Run("chr1", 1050360, 1100425, 50066, 24.8)
    summary = summarise(runs, grading=GRADING)
    assert summary.n_long == 2
    assert summary.total_long == 5437484 + 6578824
    assert summary.longest == Run("chr11", 48025953, 54604776, 6578824, 72.2)
    assert summary.isodisomy == ()


def _run(chrom: str, length: int) -> Run:
    return Run(chrom, 1_000_000, 1_000_000 + length - 1, length, 60.0)


def test_summary_counts_and_isodisomy() -> None:
    runs = [_run("chr1", 4_900_000), _run("chr1", 5_000_000), _run("chr4", 12_000_000)]
    summary = summarise(runs, grading=GRADING)
    assert (summary.n_long, summary.total_long) == (2, 17_000_000)
    assert summary.longest == runs[2]
    assert summary.isodisomy == (runs[2],)
    assert summarise([], grading=GRADING).longest is None


def _stored_text(tmp_path: Path, source: RohSource) -> tuple[dict[str, object], str]:
    collected = source.collect(run_context(tmp_path))
    assert collected.tool == ROH_TOOL
    (row,) = collected.rows
    statement = row.statement
    assert statement["statement_id"] == "roh:summary"
    assert (statement["topic"], statement["actionable"]) == ("ancestry", False)
    return statement, json.dumps(statement, default=str)


def test_roh_has_no_relatedness_wording(tmp_path: Path) -> None:
    """The finding states the runs and names a long single run as possible isodisomy,
    never parental relatedness (R3 section 4)."""
    runs = [_run("chr1", 5_000_000), _run("chr4", 12_000_000)]
    summary = summarise(runs, grading=GRADING)
    statement, text = _stored_text(
        tmp_path, RohSource(summary=summary, not_run=None, withheld=None)
    )
    assert statement["kind"] == "finding"
    reason = str(statement["reason"])
    assert "17.0 Mb in 2 runs of 5 Mb or more" in reason
    assert (
        "a single run of 12.0 Mb on chr4 may suggest isodisomy;"
        " heterodisomy cannot be seen" in reason
    )
    claim = json.loads(str(statement["claim"]))
    assert (claim["total_long"], claim["n_long"]) == (17_000_000, 2)
    assert claim["longest"]["chrom"] == "chr4"
    for word in ("related", "consanguin"):
        assert word not in text.lower()
    for not_run, withheld in (("not calibrated", None), (None, "not covered: x")):
        statement, _ = _stored_text(
            tmp_path, RohSource(summary=summary, not_run=not_run, withheld=withheld)
        )
        assert (statement["kind"], statement["reason"]) == (
            "unassessed",
            not_run or withheld,
        )


# --- the input: our genotypes and the truth's at the AF file's alleles ---------------------

AF_ROWS = [
    "chr1\t100\tA,G\t0.2",
    "chr1\t200\tC,T\t0.1",
    "chr1\t300\tG,A\t0.3",
    "chr1\t400\tT,C\t0.4",
    "chr1\t500\tA,C\t0.1",  # two ALTs at one position: left out
    "chr1\t500\tA,T\t0.1",
    "chr2\t100\tG,C\t0.5",
]


def _af(tmp_path: Path) -> Path:
    path = tmp_path / "af.EUR.tsv.gz"
    path.write_bytes(gzip.compress(("\n".join(AF_ROWS) + "\n").encode()))
    return path


def _called(dosage: str, gt: str, row_type: RowType = RowType.PASS) -> Called:
    reasons = (NOCALL_REASON,) if dosage == "unknown" else ()
    return Called(gt, dosage, row_type, 40, 30, None, None, reasons)


def _rows(vcf: Path) -> list[list[str]]:
    with gzip.open(vcf, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line[0] != "#"]


def test_sites_vcf_rows_carry_af_alleles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every AF position with one ALT is a row with the AF file's REF and ALT; a reference
    call is `0/0`, so `bcftools roh` scores it (R7); an unknown or a low call is `./.`."""
    given = {
        ("chr1", 100, "A", "G"): _called("0", "0/0", RowType.REF_BLOCK),
        ("chr1", 200, "C", "T"): _called("1", "0/1"),
        ("chr1", 300, "G", "A"): _called("unknown", "./.", RowType.NOCALL),
        ("chr1", 400, "T", "C"): _called("2", "1/1"),
        ("chr2", 100, "G", "C"): _called("2", "1/1"),
    }
    asked: list[list[Allele]] = []

    def fake_calls_at(
        gvcf: Path,
        alleles: Iterable[Allele],
        *,
        fasta: Path,
        haploid: Callable[[str, int], bool],
        scratch: Path,
    ) -> dict[Allele, Called]:
        assert (gvcf.name, fasta.name, scratch) == ("s.g.vcf.gz", "ref.fa", tmp_path)
        wanted = list(alleles)
        asked.append(wanted)
        return {a: given[a] for a in wanted}

    monkeypatch.setattr(roh, "calls_at", fake_calls_at)

    def level(allele: Allele, called: Called) -> Level | None:
        return Level.LOW if allele == ("chr2", 100, "G", "C") else Level.HIGH

    out = roh.sites_vcf(
        tmp_path / "s.g.vcf.gz",
        _af(tmp_path),
        tmp_path / "ref.fa",
        tmp_path / "sites.vcf.gz",
        haploid=lambda chrom, pos: False,
        scratch=tmp_path,
        level=level,
    )
    assert [a for chunk in asked for a in chunk] == list(given)
    assert {a[0] for a in asked[0]} == {"chr1"}  # one chromosome at a time
    assert [(r[0], r[1], r[3], r[4], r[9]) for r in _rows(out)] == [
        ("chr1", "100", "A", "G", "0/0"),
        ("chr1", "200", "C", "T", "0/1"),
        ("chr1", "300", "G", "A", "./."),
        ("chr1", "400", "T", "C", "1/1"),
        ("chr2", "100", "G", "C", "./."),  # graded low
    ]


@pytest.mark.tools
def test_sites_vcf_makes_its_scratch_folder(tmp_path: Path) -> None:
    """The run's `scratch` does not exist yet when the real call layer writes under it."""
    base = ref("chr1", 2001)
    alt = "A" if base != "A" else "C"
    af = tmp_path / "af.EUR.tsv.gz"
    af.write_bytes(gzip.compress(f"chr1\t2001\t{base},{alt}\t0.2\n".encode()))
    out = roh.sites_vcf(
        write_gvcf(tmp_path / "s.g.vcf.gz", [block("chr1", 1900, 2100)]),
        af,
        write_fasta(tmp_path / "ref.fa"),
        tmp_path / "sites.vcf.gz",
        haploid=lambda chrom, pos: False,
        scratch=tmp_path / "work" / "scratch",
        level=lambda allele, called: Level.HIGH,
    )
    assert [r[9] for r in _rows(out)] == ["0/0"]


TRUTH = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0|1
chr1\t300\t.\tG\tC\t50\tPASS\t.\tGT\t1/1
chr1\t400\t.\tT\tC,G\t50\tPASS\t.\tGT\t1/2
chr2\t100\t.\tG\tC\t50\tPASS\t.\tGT\t1|1
"""


def test_truth_sites_fill_benchmark(tmp_path: Path) -> None:
    """The truth's genotypes at the AF alleles inside the benchmark, with the AF file's REF
    and ALT: an absent site is `0/0`, a site where the truth calls another allele `./.`
    (as our calls write it), and a site outside the benchmark is left out."""
    truth = tmp_path / "truth.vcf"
    truth.write_text(TRUTH)
    indexed = Path(pysam.tabix_index(str(truth), preset="vcf", force=True))
    bench = tmp_path / "bench.bed"
    bench.write_text("chr1\t50\t450\nchr2\t0\t99\n")  # chr2:100 is outside
    out = roh.truth_sites_vcf(indexed, bench, _af(tmp_path), tmp_path / "t.vcf.gz")
    assert [(r[0], r[1], r[3], r[4], r[9]) for r in _rows(out)] == [
        ("chr1", "100", "A", "G", "0/0"),  # absent from the truth
        ("chr1", "200", "C", "T", "0/1"),
        ("chr1", "300", "G", "A", "./."),  # the truth's ALT is C
        ("chr1", "400", "T", "C", "./."),  # 1/2: another allele beside C
    ]


def _plain(path: Path, rows: Sequence[str]) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    path.write_bytes(gzip.compress((header + "\n".join(rows) + "\n").encode()))
    return path


def test_spike_rewrites_only_inside_region(tmp_path: Path) -> None:
    window = ("chr2", 50 * MB, 62 * MB)
    rows = [
        f"chr2\t{49 * MB}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  # before the window
        f"chr2\t{50 * MB}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  # its first base
        f"chr2\t{55 * MB}\t.\tC\tT\t.\tPASS\t.\tGT\t0/0",
        f"chr2\t{56 * MB}\t.\tC\tT\t.\tPASS\t.\tGT\t./.",
        f"chr2\t{62 * MB}\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",  # its last base
        f"chr2\t{62 * MB + 1}\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",
        f"chr3\t{55 * MB}\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",
    ]
    vcf = _plain(tmp_path / "in.vcf.gz", rows)
    out = roh.spike(vcf, window, tmp_path / "spiked.vcf.gz")
    assert [r[9] for r in _rows(out)] == [
        "0/1",
        "1/1",
        "0/0",
        "./.",
        "1/1",
        "0/1",
        "0/1",
    ]
    assert [r[:9] for r in _rows(out)] == [r.split("\t")[:9] for r in rows]


def test_run_roh_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`bcftools roh -G30 --AF-file` in the `truvari` image, offline, gated `dna:roh`; the
    calibration's benchmark regions as a targets file."""
    seen: list[tuple[str, list[str], dict[str, tuple[Path, bool]], str | None]] = []

    def fake(
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        targets = {t: (h, ro) for h, t, ro in mounts}
        seen.append((name, list(args), targets, label))
        out = targets["/out"][0] / Path(args[args.index("-o") + 1]).name
        out.write_text(HG002.read_text())

    monkeypatch.setattr(roh, "run_container", fake)
    vcf, af = tmp_path / "in" / "sites.vcf.gz", tmp_path / "af" / "af.EUR.tsv.gz"
    bed = tmp_path / "bench" / "bench.bed"
    out = roh.run_roh(vcf, af, tmp_path / "out" / "roh.txt", regions=bed)
    assert out == tmp_path / "out" / "roh.txt"
    assert len(read_runs(out)) == 5811
    ((name, args, mounts, label),) = seen
    assert (name, label) == ("truvari", "dna:roh")
    assert args == [
        "bcftools", "roh", "-G30", "--AF-file", "/af/af.EUR.tsv.gz",
        "-T", "/regions/bench.bed", "-O", "r", "-o", "/out/roh.txt",
        "/in/sites.vcf.gz",
    ]  # fmt: skip
    assert mounts["/in"] == (vcf.parent, True)
    assert mounts["/af"] == (af.parent, True)
    assert mounts["/regions"] == (bed.parent, True)
    assert mounts["/out"] == (out.parent, False)


@pytest.mark.tools
def test_bcftools_roh_reads_the_target_writer(tmp_path: Path) -> None:
    """`bcftools roh` (the bcftools on PATH, standing in for the `truvari` image's) reads
    what `pgsc.write_target_vcf` writes, gzipped and without `##contig` lines, against a
    tabix-indexed AF file: every site is processed and a homozygous stretch between two
    heterozygous ones is one run."""
    rng = random.Random(1)
    sites, af = [], []
    for i in range(2000):
        pos = 1000 + i * 1000
        het = (i < 200 or i > 1800) and rng.random() < 0.4
        gt, dosage = ("0/1", "1") if het else ("0/0", "0")
        called = Called(gt, dosage, RowType.PASS, 50, 30, None, None, ())
        sites.append((("chr1", pos, "A", "G"), called, Level.HIGH))
        af.append(f"chr1\t{pos}\tA,G\t0.4\n")
    vcf = write_target_vcf(iter(sites), tmp_path / "target.vcf.gz")
    plain = tmp_path / "af.tsv"
    plain.write_text("".join(af))
    indexed = pysam.tabix_index(
        str(plain), seq_col=0, start_col=1, end_col=1, force=True
    )
    out = tmp_path / "roh.txt"
    done = subprocess.run(
        ["bcftools", "roh", "-G30", "--AF-file", indexed, "-O", "r", "-o", str(out)]
        + [str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "processed: 2000/2000" in done.stderr
    (run,) = read_runs(out)
    assert run.chrom == "chr1" and run.length > 1_500_000
