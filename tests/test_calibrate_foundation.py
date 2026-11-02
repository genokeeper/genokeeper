"""The foundation table (slice 3c): the `[dosage]` and `[af_guard]` bounds proposed from the
trio's controls, the true heterozygotes' allele fractions read from vcfeval, and the
refusal of an HG004 that failed its own truth."""

from pathlib import Path

import pysam
import pytest

from genokeeper.calibrate import measure, tables
from genokeeper.calibrate.foundation import (
    CONTROLS,
    TRUE_HET_SHARE,
    build_foundation,
    propose_af_bound,
    propose_dosage,
    share_below,
)
from genokeeper.calibrate.measure import true_het_vafs
from genokeeper.calibrate.records import (
    HG004_MIN_SNP_F1,
    Record,
    reference_sha256,
    write_record,
)
from genokeeper.config import grch38_fasta
from genokeeper.pipeline import RunRefused
from genokeeper.sex import Dosage


def _measured(**changed: Dosage) -> dict[str, Dosage]:
    return {
        "pcrfree": Dosage(0.49, 0.51, 0.02, 2000),
        "pcr": Dosage(0.52, 0.47, None, 0),  # aligned only, never called
        "hg003": Dosage(0.50, 0.50, 0.03, 2000),
        "hg004": Dosage(0.98, 0.03, 0.58, 2000),
    } | changed


def test_controls_and_share() -> None:
    assert CONTROLS == {
        "pcrfree": "male",
        "pcr": "male",
        "hg003": "male",
        "hg004": "female",
    }
    assert TRUE_HET_SHARE == 0.006


def test_propose_dosage() -> None:
    proposed = propose_dosage(_measured(), 0.15)
    # (0.03 + 0.58) / 2 = 0.305, rounded up; 2 x 0.03 is 0.06 exactly, not 0.07.
    assert proposed == {"x_het_min": 0.31, "y_absent_max": 0.06}
    assert "tolerance" not in proposed
    # The floor: a female control's Y noise below 0.01 still proposes 0.02.
    quiet = _measured(hg004=Dosage(0.98, 0.004, 0.58, 2000))
    assert propose_dosage(quiet, 0.15)["y_absent_max"] == 0.02


def test_propose_dosage_refuses_a_wide_control() -> None:
    with pytest.raises(ValueError, match=r"pcr\b.*0\.08"):
        propose_dosage(_measured(pcr=Dosage(0.58, 0.47, None, 0)), 0.15)


def test_propose_dosage_refuses_overlap() -> None:
    with pytest.raises(ValueError, match="hg004"):
        propose_dosage(_measured(hg004=Dosage(0.98, 0.03, 0.08, 2000)), 0.15)


def test_propose_af_bound() -> None:
    vafs = [0.1, 0.12, 0.15, 0.2, 0.25, 0.28, 0.293] + [0.5] * 993
    assert len(vafs) == 1000
    # Six of 1,000 lie below 0.29 (at most 0.6%); seven lie below 0.30.
    assert propose_af_bound(vafs) == 0.29
    assert share_below(vafs, 0.30) == 7 / 1000


_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##INFO=<ID=CALL,Number=1,Type=String,Description="vcfeval">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
    '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction">\n'
    "##contig=<ID=chr1,length=100000>\n"
    "##contig=<ID=chrX,length=100000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


def _vcf(path: Path, rows: list[str]) -> Path:
    plain = path.with_suffix("")
    plain.write_text(_HEADER + "".join(f"{row}\n" for row in rows))
    return Path(pysam.tabix_index(str(plain), preset="vcf", force=True))


def test_true_het_vafs_reads_tp_hets(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The TP heterozygous SNVs on chr1-22, each its FORMAT/VAF (DeepVariant's AD[alt] /
    DP, every good read in DP), not AD[alt] / (AD[ref] + AD[alt])."""
    fmt = "GT:AD:DP:VAF"
    annotated = _vcf(
        tmp_path / "annotated.vcf.gz",
        [
            f"chr1\t100\t.\tA\tG\t50\tPASS\tCALL=TP\t{fmt}\t0/1:14,16:32:0.5",
            f"chr1\t200\t.\tA\tG\t50\tPASS\tCALL=FP\t{fmt}\t0/1:14,16:30:0.533",
            f"chr1\t300\t.\tA\tG\t50\tPASS\tCALL=TP\t{fmt}\t1/1:0,30:30:1",
            f"chr1\t400\t.\tAT\tA\t50\tPASS\tCALL=TP\t{fmt}\t0/1:15,15:30:0.5",
            f"chrX\t100\t.\tA\tG\t50\tPASS\tCALL=TP\t{fmt}\t0/1:15,15:30:0.5",
        ],
    )
    calls = _vcf(
        tmp_path / "calls.vcf.gz",
        [f"chr1\t100\t.\tA\tG\t50\tPASS\t.\t{fmt}\t0/1:14,16:32:0.5"],
    )
    seen: list[bool] = []

    def vcfeval(*args: object, squash: bool, threads: int) -> None:
        seen.append(squash)
        out = args[5]
        assert isinstance(out, Path)
        out.mkdir(parents=True)
        (out / "calls.vcf.gz").write_bytes(annotated.read_bytes())

    monkeypatch.setattr(measure, "_vcfeval", vcfeval)
    rtg = tmp_path / "rtg"
    rtg.write_text("placeholder\n")
    bench = tmp_path / "bench.bed"
    bench.write_text("chr1\t0\t100000\n")
    vafs = true_het_vafs(calls, calls, bench, sdf=tmp_path / "sdf", rtg=rtg, threads=2)
    assert vafs == [0.5]
    assert seen == [False]  # genotype-exact


def test_foundation_refuses_unaccepted_hg004(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    fasta = grch38_fasta(lab)
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\nACGT\n")
    cram = tmp_path / "hg004.cram"
    cram.write_bytes(b"")
    write_record(
        lab,
        "hg004",
        Record(
            cram=str(cram),
            key="k",
            metrics={},
            gvcf=None,
            reference_sha256=reference_sha256(lab),
            accepted=False,
            snp_f1=0.62,
            replaced=None,
        ),
    )
    with pytest.raises(RunRefused, match=str(HG004_MIN_SNP_F1)):
        build_foundation(lab, 2)
    assert HG004_MIN_SNP_F1 == 0.99
    assert not tables.folder(lab, "foundation").exists()


def test_foundation_is_registered_after_3c_tables() -> None:
    """3c's rows in order, then slice 3f's (`slice3f`, which needs none of them)."""
    assert list(tables.TABLES) == [
        "small_variants",
        "hg002_cells",
        "stage3",
        "deletions",
        "cnv",
        "insertions",
        "paralogs",
        "foundation",
        "slice3f",
    ]
    assert tables.TABLES["hg002_cells"].needs == ("small_variants",)
    assert tables.TABLES["foundation"].needs == ()
