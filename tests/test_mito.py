"""chrM from mutserve: its output read with each variant's level, and a chrM claim's call from
those levels. mutserve's container is replaced by a stub that writes what mutserve 2.0.3
wrote on a synthetic chrM BAM (run on 2026-10-20: a homoplasmic change, a 20% and a 3%
heteroplasmy at 200x; the BioContainers image `genokeeper.containers.IMAGES["mutserve"]`
wrote the same columns on such a BAM)."""

import random
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import pysam
import pytest
from fixtures.run_context import run_context

from genokeeper import mito, tools
from genokeeper.config import load_grading
from genokeeper.containers import PLATFORM
from genokeeper.genotype.gvcf import NO_COVERAGE_REASON, OTHER_ALLELE_WARNING, Called
from genokeeper.genotype.reliability import OUTSIDE_CALIBRATED, StrataIndex
from genokeeper.haplogroup import write_chrm_vcf
from genokeeper.mito import (
    HETEROPLASMY_WARNING,
    INDEL_REASON,
    MATERNAL_NOTE,
    TISSUE_NOTE,
    MitoCall,
    call_chrm,
    grade_mito,
    mito_called,
)
from genokeeper.panel.build import Claim, ClaimSet
from genokeeper.pgx.mtrnr1 import mtrnr1_outside
from genokeeper.query import GENOTYPE_UNKNOWN
from genokeeper.sources import Row
from genokeeper.sources.panel import ClaimsSource
from genokeeper.store.db import Level, RowType
from genokeeper.variants.key import make_vid

# mutserve 2.0.3's text output, verbatim but for the sample ID.
MUTSERVE_TXT = (
    "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMajorBase\tMajorLevel\tMinorBase"
    "\tMinorLevel\tCoverage\tCoverageFWD\tCoverageREV\tType\tMeanBaseQuality\n"
    "chrM.bam\tPASS\t3000\tA\tG\t1.00\tG\t\t-\t\t200\t94\t106\t1\t40\n"
    "chrM.bam\tPASS\t5000\tA\tG\t0.198\tA\t0.802\tG\t0.198\t202\t112\t90\t2\t40\n"
    "chrM.bam\tPASS\t8000\tG\tA\t0.029\tG\t0.971\tA\t0.029\t171\t91\t80\t2\t40\n"
)
# The first columns of its raw per-position table (`--write-raw`), as it wrote them.
RAW_HEADER = (
    "SAMPLE\tPOS\tREF\tTOP-FWD\tMINOR-FWD\tTOP-REV\tMINOR-REV\tCOV-FWD\tCOV-REV"
    "\tCOV-TOTAL\tTYPE\tLEVEL"
)
RAW_ROWS = {4: 2, 100: 150, 3000: 200, 6000: 3}


def _raw() -> str:
    rows = [
        f"chrM.bam\t{pos}\tC\tC\t-\t-\t-\t{cov}\t0\t{cov}\t0\t0.0"
        for pos, cov in RAW_ROWS.items()
    ]
    return RAW_HEADER + "\n" + "\n".join(rows) + "\n"


class FakeMutserve:
    """Stands in for `run_container` in `genokeeper.mito`: checks the chrM inputs it is
    handed and writes mutserve's outputs, its VCF among them (`--output` names it; the
    text table is named from the part before the first dot), into the work directory."""

    def __init__(self) -> None:
        self.args: list[str] = []
        self.chrm_reads = -1

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, platform, label) == ("mutserve", PLATFORM, "dna:mito")
        self.args = list(args)
        ((work, target, read_only),) = mounts
        assert (target, read_only) == ("/work", False)
        with pysam.AlignmentFile(str(work / "chrM.bam")) as bam:
            self.chrm_reads = sum(1 for _ in bam)
            assert {r.reference_name for r in bam.fetch()} <= {"chrM"}
        with pysam.FastaFile(str(work / "chrM.fa")) as fa:
            assert list(fa.references) == ["chrM"]
        assert args[list(args).index("--output") + 1] == "/work/mutserve.vcf.gz"
        (work / "mutserve.txt").write_text(MUTSERVE_TXT)
        (work / "mutserve_raw.txt").write_text(_raw())
        (work / "mutserve.vcf.gz").write_bytes(b"\x1f\x8b vcf")


def _reads(tmp_path: Path) -> tuple[Path, Path]:
    """A FASTA with chr1 and chrM, and a BAM with 5 reads on each."""
    rng = random.Random(3)
    seqs = {"chr1": 2000, "chrM": 16569}
    fasta = tmp_path / "ref.fa"
    text = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in seqs.items()}
    fasta.write_text("".join(f">{c}\n{s}\n" for c, s in text.items()))
    pysam.faidx(str(fasta))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}}
    header["SQ"] = [{"SN": c, "LN": n} for c, n in seqs.items()]  # pyright: ignore[reportArgumentType]
    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for contig in seqs:
            for i in range(5):
                read = pysam.AlignedSegment(out.header)
                read.query_name, read.flag = f"{contig}-{i}", 0
                read.reference_name, read.reference_start = contig, 100 * i
                read.mapping_quality, read.cigarstring = 60, "50M"
                read.query_sequence = text[contig][100 * i : 100 * i + 50]
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(read)
    pysam.index(str(bam))
    return fasta, bam


@pytest.mark.tools
def test_mutserve_output_parsed_with_levels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")
    fake = FakeMutserve()
    monkeypatch.setattr(mito, "run_container", fake)
    fasta, bam = _reads(tmp_path)
    out = tmp_path / "new" / "mito"  # mutserve exits 255 without its output folder

    calls = call_chrm(bam, fasta, out)

    assert calls == [
        MitoCall(pos=3000, ref="A", alt="G", level=1.0, dp=200),
        MitoCall(pos=5000, ref="A", alt="G", level=0.198, dp=202),
        MitoCall(pos=8000, ref="G", alt="A", level=0.029, dp=171),
    ]
    # Only the chrM reads and the chrM sequence are handed to the container.
    assert fake.chrm_reads == 5
    assert fake.args[:2] == ["mutserve", "call"]
    assert "--write-raw" in fake.args and "/work/chrM.bam" in fake.args
    assert fake.args[fake.args.index("--level") + 1] == "0.01"
    assert fake.args[fake.args.index("--reference") + 1] == "/work/chrM.fa"
    assert mito.chrm_coverage(out) == RAW_ROWS
    assert mito.read_calls(out) == calls
    # The extracted reads are not kept beside the calls.
    assert sorted(p.name for p in out.iterdir()) == [
        "mutserve.txt",
        "mutserve.vcf.gz",
        "mutserve_raw.txt",
    ]


def test_mito_threshold_in_config() -> None:
    assert load_grading().mito.min_homoplasmy == 0.9


@pytest.mark.parametrize(
    ("allele", "dosage", "row_type", "warning", "level"),
    [
        # A homoplasmic call: one copy of the chrM genome, never high without a
        # calibrated error estimate.
        (("chrM", 3000, "A", "G"), "1", RowType.PASS, None, Level.MEDIUM),
        # A heteroplasmy: carried at its level, low.
        (("chrM", 5000, "A", "G"), "1", RowType.PASS, HETEROPLASMY_WARNING, Level.LOW),
        # Another base called at the claim's position.
        (("chrM", 3000, "A", "T"), "0", RowType.PASS, OTHER_ALLELE_WARNING, None),
        # No call, and reads cover the position: the reference.
        (("chrM", 100, "C", "T"), "0", RowType.REF_BLOCK, None, Level.MEDIUM),
        # Too few reads, or none.
        (
            ("chrM", 6000, "C", "T"),
            "unknown",
            RowType.NONE,
            NO_COVERAGE_REASON,
            Level.LOW,
        ),
        (
            ("chrM", 7000, "C", "T"),
            "unknown",
            RowType.NONE,
            NO_COVERAGE_REASON,
            Level.LOW,
        ),
        # mutserve calls no indels here.
        (("chrM", 100, "CA", "C"), "unknown", RowType.NONE, INDEL_REASON, Level.LOW),
    ],
)
def test_chrm_claims_called_from_levels(
    allele: tuple[str, int, str, str],
    dosage: str,
    row_type: RowType,
    warning: str | None,
    level: Level | None,
) -> None:
    calls = [
        MitoCall(pos=3000, ref="A", alt="G", level=0.99, dp=200),
        MitoCall(pos=5000, ref="A", alt="G", level=0.2, dp=202),
    ]
    grading = load_grading()
    called = mito_called(
        allele, calls, RAW_ROWS, min_dp=grading.min_dp, min_homoplasmy=0.9
    )
    assert (called.dosage, called.row_type) == (dosage, row_type)
    assert called.gq is None
    if warning is not None:
        assert any(w.startswith(warning) for w in called.warnings)
    if dosage == "1":
        assert called.vaf == (0.99 if allele[1] == 3000 else 0.2)
    graded = grade_mito(
        called, allele[0], allele[1], strata=StrataIndex({}), min_homoplasmy=0.9
    )
    if level is not None:
        assert graded.level == level
    assert graded.level != Level.HIGH and graded.p_wrong is None
    assert mito.MUTSERVE_WARNING in graded.warnings or dosage == "unknown"
    # HG003's truth covers chr1-22: a known chrM genotype is outside calibrated regions.
    assert (OUTSIDE_CALIBRATED in graded.warnings) == (dosage != "unknown")


def test_a_variant_on_too_few_reads_is_not_carried() -> None:
    # The haplogroup's VCF (`write_chrm_vcf`) leaves it out as no evidence; so does the claim.
    calls = [MitoCall(pos=1555, ref="A", alt="G", level=1.0, dp=3)]
    called = mito_called(
        ("chrM", 1555, "A", "G"), calls, {1555: 3}, min_dp=10, min_homoplasmy=0.9
    )
    assert (called.dosage, called.dp) == ("unknown", 3)
    assert NO_COVERAGE_REASON in called.warnings


# --- slice 3f: the heteroplasmy floor, contamination and a withheld heteroplasmy ---

LHON = ("chrM", 11778, "G", "A")


def _lhon_call(level: float, **keywords: Any) -> Called:
    calls = [MitoCall(pos=11778, ref="G", alt="A", level=level, dp=2000)]
    return mito_called(
        LHON, calls, {11778: 2000}, min_dp=10, min_homoplasmy=0.9, **keywords
    )


def _lhon_claim(tmp_path: Path, level: float, **keywords: Any) -> Row:
    """A ClinVar pathogenic chrM claim at m.11778G>A (LHON) on mutserve's call at `level`."""
    claim = Claim(
        variation_id=9681,
        gene="MT-ND4",
        clnsig="Pathogenic",
        clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
        disease_ids="MedGen:C0917796",
        condition="Leber optic atrophy",
        chrom="chrM",
        pos=11778,
        ref="G",
        alt="A",
        vid=make_vid(*LHON),
        rsid="rs199476112",
    )
    called = _lhon_call(level, **keywords)
    graded = grade_mito(
        called, "chrM", 11778, strata=StrataIndex({}), min_homoplasmy=0.9
    )
    ctx = replace(run_context(tmp_path), haploid=lambda chrom, pos: chrom == "chrM")
    source = ClaimsSource(
        ClaimSet(claims=[claim], input_rows=1, dropped={}), [(claim, called, graded)]
    )
    (row,) = source.collect(ctx).rows
    return row


def test_heteroplasmy_below_floor_unassessed(tmp_path: Path) -> None:
    """Review Focus 1: a LHON variant on a ClinVar pathogenic claim at 8% heteroplasmy,
    below the floor, is unassessed with the level and the floor named, never a finding
    and never normal."""
    low = _lhon_claim(tmp_path, 0.08, floor=0.10)
    assert low.statement["kind"] == "unassessed"
    assert low.call is not None
    assert (low.call["dosage"], low.call["vaf"]) == ("unknown", 0.08)
    warnings = list(low.call["warnings"])  # pyright: ignore[reportArgumentType]
    (floor_warning,) = [w for w in warnings if "below the floor" in w]
    assert floor_warning.startswith(
        GENOTYPE_UNKNOWN
    )  # what `query` reads as the reason
    assert "heteroplasmy at level 0.08, below the floor 0.1" in floor_warning
    assert MATERNAL_NOTE not in warnings

    mid = _lhon_claim(tmp_path, 0.12, floor=0.10)
    assert mid.statement["kind"] == "finding"
    assert mid.call is not None
    assert TISSUE_NOTE in mid.call["warnings"]  # pyright: ignore[reportOperatorIssue]
    assert MATERNAL_NOTE in mid.call["warnings"]  # pyright: ignore[reportOperatorIssue]

    high = _lhon_claim(tmp_path, 0.95, floor=0.10)
    assert high.statement["kind"] == "finding"
    assert high.call is not None
    assert high.call["level"] == str(Level.MEDIUM)
    warnings = list(high.call["warnings"])  # pyright: ignore[reportArgumentType]
    assert not any(w.startswith(HETEROPLASMY_WARNING) for w in warnings)
    assert TISSUE_NOTE not in warnings and MATERNAL_NOTE in warnings


def test_mito_called_defaults_are_todays() -> None:
    """Without the new keywords a 0.08 call is what `main` gives: carried, at its level,
    warned a heteroplasmy, no tissue note (the existing caller's contract)."""
    assert _lhon_call(0.08) == Called(
        gt="1",
        dosage="1",
        row_type=RowType.PASS,
        gq=None,
        dp=2000,
        ad=None,
        vaf=0.08,
        warnings=("heteroplasmy at level 0.08",),
    )


def test_withheld_heteroplasmy_unknown() -> None:
    reason = "not covered: false heteroplasmy at 0.25 in calibration"
    for floor in (None, 0.1):
        called = _lhon_call(0.40, floor=floor, heteroplasmy_withheld=reason)
        assert (called.dosage, called.vaf) == ("unknown", 0.40)
        assert any(
            w.startswith(GENOTYPE_UNKNOWN) and reason in w for w in called.warnings
        )
        homoplasmic = _lhon_call(0.95, floor=floor, heteroplasmy_withheld=reason)
        assert (homoplasmic.dosage, homoplasmic.warnings) == ("1", ())


def test_mt_rnr1_and_haplogroup_unchanged(tmp_path: Path) -> None:
    """MT-RNR1 and the haplogroup keep [mito] min_homoplasmy: a call above the floor but
    below it is neither MT-RNR1's allele nor in the haplogroup's VCF."""
    calls = [
        MitoCall(pos=1555, ref="A", alt="G", level=0.12, dp=500),
        MitoCall(pos=263, ref="A", alt="G", level=0.5, dp=500),
        MitoCall(pos=750, ref="A", alt="G", level=0.95, dp=500),
    ]
    outside = mtrnr1_outside(calls, ["m.1555A>G"], min_homoplasmy=0.9)
    assert outside.diplotype is None
    assert outside.unresolved_reason == "m.1555A>G at heteroplasmy 0.12, below 0.9"
    vcf = write_chrm_vcf(
        calls,
        {1555: 500, 263: 500, 750: 500},
        tmp_path / "chrM.vcf",
        min_homoplasmy=0.9,
        min_dp=10,
    )
    rows = [line for line in vcf.read_text().splitlines() if not line.startswith("#")]
    assert [row.split("\t")[1] for row in rows] == ["750"]


def test_floor_in_config() -> None:
    assert load_grading().mito.heteroplasmy_floor == 0.10
