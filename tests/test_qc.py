"""The QC gate and sex inference, on mini CRAMs measured by the real mosdepth, samtools and
bcftools. VerifyBamID2's container is replaced by a stub that writes its `.selfSM`."""

from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import pysam
import pytest
from fixtures.grch38_mini import CONTIGS, write_cram, write_reference, write_regions

from genokeeper import qc, tools
from genokeeper.config import Grading, load_grading
from genokeeper.qc import run_qc
from genokeeper.regions import QC_AUTOSOMES, SEX_REGIONS, giab38_dir
from genokeeper.sex import SexRefused, copies, infer_sex, measure_ratios

# VerifyBamID2's output columns (its README, "Output"), with FREEMIX set per test.
SELFSM_HEADER = (
    "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\tFREE_RH"
    "\tFREE_RA\tCHIPMIX\tCHIPLK1\tCHIPLK0\tCHIP_RH\tCHIP_RA\tDPREF\tRDPHET\tRDPALT"
)


@pytest.fixture(autouse=True)
def no_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")


@pytest.fixture
def lab(tmp_path: Path) -> Path:
    """A lab whose derived region BEDs cover the mini reference: chr1 as the autosomes,
    chrY as the X-degenerate region, chrX as chrX's mappable non-PAR region."""
    home = tmp_path / "lab"
    giab = giab38_dir(home)
    autosomes = ("chr1", 0, CONTIGS["chr1"], "autosomes")
    write_regions(giab / QC_AUTOSOMES, [autosomes])
    write_regions(
        giab / SEX_REGIONS,
        [
            autosomes,
            ("chrY", 0, CONTIGS["chrY"], "xdr"),
            ("chrX", 0, CONTIGS["chrX"], "xnonpar"),
        ],
    )
    return home


class FakeVerifyBamID:
    """Stands in for `run_container`: writes `vb.selfSM` with FREEMIX `value` into the
    directory mounted at /out, and keeps each call's arguments."""

    def __init__(self) -> None:
        self.value = 0.001
        self.calls: list[list[str]] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, platform, label) == ("verifybamid2", "linux/amd64", "dna:qc")
        self.calls.append(list(args))
        out = next(host for host, target, _ in mounts if target == "/out")
        row = ["mini", "NA", "NA", "9000", "90000", "30", str(self.value), *["NA"] * 12]
        (out / "vb.selfSM").write_text(f"{SELFSM_HEADER}\n" + "\t".join(row) + "\n")


@pytest.fixture
def freemix(monkeypatch: pytest.MonkeyPatch) -> FakeVerifyBamID:
    fake = FakeVerifyBamID()
    monkeypatch.setattr(qc, "run_container", fake)
    return fake


@pytest.mark.tools
def test_qc_fails_low_depth(
    tmp_path: Path, lab: Path, freemix: FakeVerifyBamID
) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    cram = write_cram(tmp_path / "low.cram", fasta, {"chr1": 5})
    metrics = run_qc(cram, fasta, grading=load_grading(), lab=lab, threads=2)
    assert 4 < metrics["mean_depth"] < 6
    assert metrics["fraction_20x"] == 0
    assert metrics["mapped_fraction"] == 1
    assert metrics["duplicate_fraction"] == 0
    assert metrics["median_insert"] == 300
    assert metrics["freemix"] == 0.001
    failed = qc.judge(metrics, load_grading().qc)
    assert [f.split()[0] for f in failed] == ["mean_depth", "fraction_20x"]


@pytest.mark.tools
def test_contaminated_sample_refused(
    tmp_path: Path, lab: Path, freemix: FakeVerifyBamID
) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    cram = write_cram(tmp_path / "deep.cram", fasta, {"chr1": 30})
    grading = load_grading()
    # Random placement leaves part of chr1 below 20x: only contamination is under test.
    grading = replace(grading, qc=replace(grading.qc, min_fraction_20x=0.0))
    freemix.value = 0.08
    metrics = run_qc(cram, fasta, grading=grading, lab=lab, threads=2)
    assert qc.judge(metrics, grading.qc) == ["freemix 0.08 > 0.03"]
    args = freemix.calls[-1]
    assert args[args.index("--Reference") + 1] == f"/ref/{fasta.name}"
    assert args[args.index("--BamFile") + 1] == f"/reads/{cram.name}"
    assert args[args.index("--NumThread") + 1] == "2"


# A sample inside every shipped [qc] threshold.
PASSING = {
    "mean_depth": 30.0,
    "fraction_20x": 0.9,
    "mapped_fraction": 0.995,
    "duplicate_fraction": 0.05,
    "median_insert": 350.0,
    "freemix": 0.005,
    "titv": 2.05,
    "het_hom": 1.6,
}


def test_passing_sample_fails_nothing() -> None:
    assert qc.judge(PASSING, load_grading().qc) == []


@pytest.mark.parametrize(
    ("metric", "value", "line"),
    [
        ("mapped_fraction", 0.97, "mapped_fraction 0.97 < 0.98"),
        ("duplicate_fraction", 0.25, "duplicate_fraction 0.25 > 0.1"),
        ("median_insert", 180.0, "median_insert 180 < 250"),
        ("het_hom", 3.5, "het_hom 3.5 > 3"),
        ("titv", 1.5, "titv 1.5 < 1.9"),
    ],
)
def test_each_metric_outside_its_threshold_refused(
    metric: str, value: float, line: str
) -> None:
    """The metrics no measured fixture drives past its threshold: each one alone fails."""
    assert qc.judge({**PASSING, metric: value}, load_grading().qc) == [line]


def _vcf(path: Path, rows: list[tuple[int, str, str, str, str]]) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=RefCall,Description="Reference call">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CONTIGS.items())
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmini\n"
    )
    body = "".join(
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t30\t{flt}\t.\tGT\t{gt}\n"
        for pos, ref, alt, flt, gt in rows
    )
    plain = path.with_suffix("")
    plain.write_text(header + body)
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


@pytest.mark.tools
def test_titv_and_het_hom_from_pass_snvs(tmp_path: Path) -> None:
    vcf = _vcf(
        tmp_path / "calls.vcf.gz",
        [
            (100, "A", "G", "PASS", "0/1"),  # transition, het
            (200, "C", "T", "PASS", "1/1"),  # transition, hom
            (300, "A", "C", "PASS", "0/1"),  # transversion, het
            (400, "G", "A", "PASS", "0/1"),  # transition, het
            (500, "A", "T", "RefCall", "1/1"),  # not PASS: not counted
            (600, "AC", "A", "PASS", "1/1"),  # an indel: not an SNV
        ],
    )
    stats = qc.vcf_stats(vcf)
    assert stats == {"titv": 3.0, "het_hom": 3.0}
    assert qc.judge(PASSING | stats, load_grading().qc) == ["titv 3 > 2.3"]


@pytest.mark.tools
@pytest.mark.parametrize(
    ("chrx_depth", "chry_depth", "expected"),
    [(10, 10, "male"), (10, 2, "male"), (20, 0, "female"), (15, 0, None)],
)
def test_sex_from_x_copies_depth(
    tmp_path: Path, lab: Path, chrx_depth: int, chry_depth: int, expected: str | None
) -> None:
    # The X copies decide: a Y at 0.1 of the autosomes (mosaic loss) is still a male.
    fasta = write_reference(tmp_path / "ref.fa")
    depths = {"chr1": 20, "chrX": chrx_depth} | (
        {"chrY": chry_depth} if chry_depth else {}
    )
    cram = write_cram(tmp_path / "sample.cram", fasta, depths)
    grading: Grading = load_grading()
    if expected is None:
        # chrX at 0.75 of the autosomes: neither one copy nor two.
        with pytest.raises(SexRefused, match="--sex"):
            infer_sex(cram, grading=grading, fasta=fasta, lab=lab, threads=2)
    else:
        assert (
            infer_sex(cram, grading=grading, fasta=fasta, lab=lab, threads=2)
            == expected
        )


@pytest.mark.tools
def test_sex_ignores_reads_below_the_mapq(tmp_path: Path, lab: Path) -> None:
    # A female's reads from chrX homologues land on chrY at low MAPQ: they are not a Y.
    fasta = write_reference(tmp_path / "ref.fa")
    cram = tmp_path / "sample.cram"
    write_cram(tmp_path / "auto.cram", fasta, {"chr1": 20, "chrX": 20})
    write_cram(tmp_path / "y.cram", fasta, {"chrY": 10}, mapq=10, seed=3)
    pysam.merge(
        "-f",
        "--write-index",
        "-O",
        "cram",
        "--reference",
        str(fasta),
        str(cram),
        str(tmp_path / "auto.cram"),
        str(tmp_path / "y.cram"),
    )
    grading = load_grading()
    x_ratio, y_ratio = measure_ratios(
        cram, grading=grading, fasta=fasta, lab=lab, threads=2
    )
    assert y_ratio == 0.0 and copies(x_ratio, grading.dosage.tolerance) == 2
    assert infer_sex(cram, grading=grading, fasta=fasta, lab=lab, threads=2) == "female"
