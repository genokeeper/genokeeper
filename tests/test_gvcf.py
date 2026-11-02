"""The gVCF call layer: a genotype, or unknown with its reason, at every claimed allele.

Records are DeepVariant 1.10's real shapes (`fixtures.gvcf_mini`); `bcftools norm` runs for
real, so these tests need bcftools.
"""

import tempfile
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import pytest
from fixtures.gvcf_mini import (
    HOMOPOLYMER,
    PAR_END,
    block,
    other_base,
    ref,
    star_record,
    variant,
    write_fasta,
    write_gvcf,
)

from genokeeper.calibrate.table import VERSION_KEYS, CalibrationTable
from genokeeper.config import load_grading
from genokeeper.genotype import gvcf as gvcf_module
from genokeeper.genotype.gvcf import (
    INDEL_BLOCK_WARNING,
    INSIDE_DELETION_REASON,
    MULTI_BASE_REASON,
    NO_COVERAGE_REASON,
    NOCALL_REASON,
    OTHER_ALLELE_WARNING,
    Called,
    calls_at,
)
from genokeeper.genotype.reliability import (
    DEPTH_WARNING,
    HAPLOID_HET,
    OUTSIDE_CALIBRATED,
    Reliability,
    StrataIndex,
    grade_called,
)
from genokeeper.store.db import Level, RowType
from genokeeper.tools import run_tool

pytestmark = pytest.mark.tools

GRADING = load_grading()
NO_STRATA = StrataIndex.from_beds({})
# No benchmark and no cells: every call is graded from its GQ, as outside HG003's benchmark
# (the measured grades: tests/test_reliability.py).
UNMEASURED = CalibrationTable(
    cells={}, versions={key: f"{key}-value" for key in VERSION_KEYS}
)
MEDIANS = {"chr1": 30.0, "chrX": 15.0}

type Claim = tuple[str, int, str, str]


def diploid(chrom: str, pos: int) -> bool:
    del chrom, pos
    return False


def male(chrom: str, pos: int) -> bool:
    """One copy: a male's chrX outside the mini PAR."""
    return chrom == "chrX" and pos > PAR_END


def _calls(
    tmp_path: Path,
    records: list[str],
    claims: list[Claim],
    haploid: Callable[[str, int], bool] = diploid,
) -> dict[Claim, Called]:
    fasta = write_fasta(tmp_path / "ref.fa")
    gvcf = write_gvcf(tmp_path / "mini.g.vcf.gz", records)
    return calls_at(gvcf, claims, fasta=fasta, haploid=haploid, scratch=tmp_path)


def _snv(chrom: str, pos: int, k: int = 0) -> Claim:
    base = ref(chrom, pos)
    return (chrom, pos, base, other_base(base, k))


def _grade(
    called: Called, claim: Claim, *, indel: bool = False, haploid: bool = False
) -> Reliability:
    return grade_called(
        called,
        claim[0],
        claim[1],
        strata=NO_STRATA,
        medians=MEDIANS,
        grading=GRADING,
        indel=indel,
        haploid=haploid,
        table=UNMEASURED,
        benchmark=NO_STRATA,
    )


def test_variant_record_dosage(tmp_path: Path) -> None:
    het, hom = _snv("chr1", 2001), _snv("chr1", 2101)
    records = [
        block("chr1", 1900, 2000),
        variant("chr1", 2001, 1, [het[3]], gt="0/1", gq=37, ad=[14, 16]),
        block("chr1", 2002, 2100),
        variant("chr1", 2101, 1, [hom[3]], gt="1/1", gq=40, ad=[0, 26]),
        block("chr1", 2102, 2200),
    ]
    calls = _calls(tmp_path, records, [het, hom])
    assert replace(calls[het], vaf=None) == Called(
        gt="0/1",
        dosage="1",
        row_type=RowType.PASS,
        gq=37,
        dp=30,
        ad=(14, 16, 0),
        vaf=None,
        warnings=(),
    )
    assert calls[het].vaf == pytest.approx(16 / 30, abs=1e-5)
    assert (calls[hom].gt, calls[hom].dosage, calls[hom].gq) == ("1/1", "2", 40)


def test_called_alleles_written_only_under_the_runs_work_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The alleles `bcftools norm` reads (the user's calls) go under the run's work
    directory, which the deny rules cover, and are gone afterwards; never the system temp."""
    system_tmp = tmp_path / "system-tmp"
    system_tmp.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(system_tmp))
    scratch = tmp_path / "reports" / "work"
    scratch.mkdir(parents=True)
    read: list[Path] = []

    def spy(args: list[str]) -> str:
        read.append(Path(args[-1]))
        return run_tool(args)

    monkeypatch.setattr(gvcf_module, "run_tool", spy)
    het = _snv("chr1", 2001)
    records = [variant("chr1", 2001, 1, [het[3]], gt="0/1", gq=37, ad=[14, 16])]
    fasta = write_fasta(tmp_path / "ref.fa")
    gvcf = write_gvcf(tmp_path / "mini.g.vcf.gz", records)
    calls = calls_at(gvcf, [het], fasta=fasta, haploid=diploid, scratch=scratch)
    assert calls[het].dosage == "1"
    assert len(read) == 2  # the claims, then the gVCF's records near them
    assert all(p.is_relative_to(scratch) and not p.exists() for p in read)
    assert (list(scratch.iterdir()), list(system_tmp.iterdir())) == ([], [])


def test_refcall_zero_with_network_gq(tmp_path: Path) -> None:
    claim = _snv("chr1", 3001)
    # A RefCall deletion: its REF span has no block (as after a called deletion).
    spanned = _snv("chr1", 3202)
    records = [
        block("chr1", 2900, 3000),
        variant(
            "chr1",
            3001,
            1,
            [claim[3]],
            gt="0/0",
            gq=32,
            ad=[10, 2],
            filt="RefCall",
            qual=0,
        ),
        block("chr1", 3002, 3200),
        variant(
            "chr1",
            3201,
            3,
            [ref("chr1", 3201)],
            gt="0/0",
            gq=21,
            ad=[12, 3],
            filt="RefCall",
            qual=0,
        ),
        block("chr1", 3204, 3300),
    ]
    calls = _calls(tmp_path, records, [claim, spanned])
    assert (calls[spanned].dosage, calls[spanned].row_type, calls[spanned].gq) == (
        "0",
        RowType.REFCALL,
        21,
    )
    called = calls[claim]
    assert (called.gt, called.dosage, called.row_type, called.gq) == (
        "0/0",
        "0",
        RowType.REFCALL,
        32,
    )
    assert called.warnings == ()


def test_block_zero_with_block_gq(tmp_path: Path) -> None:
    inside, capped, shallow = _snv("chr1", 1550), _snv("chr1", 1750), _snv("chr1", 1950)
    records = [
        block("chr1", 1500, 1600, gq=37, min_dp=15),
        block("chr1", 1700, 1800, gq=60, min_dp=25),  # never seen from DeepVariant
        block("chr1", 1900, 2000, gq=12, min_dp=4),
    ]
    calls = _calls(tmp_path, records, [inside, capped, shallow])
    assert calls[inside] == Called(
        gt="0/0",
        dosage="0",
        row_type=RowType.REF_BLOCK,
        gq=37,
        dp=15,
        ad=None,
        vaf=None,
        warnings=(),
    )
    assert calls[capped].gq == 50  # binomial, capped at 50 (spec, 2a)
    assert _grade(calls[inside], inside).level == Level.HIGH
    # A block's depth is its MIN_DP.
    assert calls[shallow].dp == 4
    assert DEPTH_WARNING in _grade(calls[shallow], shallow).warnings
    assert DEPTH_WARNING not in _grade(calls[inside], inside).warnings


def test_block_at_indel_capped_medium(tmp_path: Path) -> None:
    deletion = ("chr1", 1550, ref("chr1", 1550, 3), ref("chr1", 1550))
    snv = _snv("chr1", 1550)
    calls = _calls(tmp_path, [block("chr1", 1500, 1600)], [deletion, snv])
    assert (calls[deletion].dosage, calls[deletion].row_type) == (
        "0",
        RowType.REF_BLOCK,
    )
    graded = _grade(calls[deletion], deletion, indel=True)
    assert graded.level == Level.MEDIUM
    assert INDEL_BLOCK_WARNING in graded.warnings
    assert _grade(calls[snv], snv).level == Level.HIGH  # the cap is the indel's alone


def test_nocall_unknown(tmp_path: Path) -> None:
    record, in_block = _snv("chr1", 2001), _snv("chr1", 2201)
    records = [
        variant(
            "chr1",
            2001,
            1,
            [record[3]],
            gt="./.",
            gq=14,
            ad=[8, 8],
            filt="NoCall",
            qual=0.2,
        ),
        block("chr1", 2201, 2201, gq=0, min_dp=34, gt="./."),
    ]
    calls = _calls(tmp_path, records, [record, in_block])
    for claim in (record, in_block):
        called = calls[claim]
        assert (called.gt, called.dosage, called.row_type) == (
            "./.",
            "unknown",
            RowType.NOCALL,
        )
        assert called.warnings == (NOCALL_REASON,)
        assert _grade(called, claim).level == Level.LOW
    assert calls[record].gq == 14


def test_nocall_record_on_a_deleted_base_unknown(tmp_path: Path) -> None:
    # No record at the deletion's position; a NoCall record on a base it deletes.
    deletion = ("chr1", 3001, ref("chr1", 3001, 4), ref("chr1", 3001))
    alt = other_base(ref("chr1", 3003))
    records = [
        block("chr1", 2900, 3002),
        variant(
            "chr1", 3003, 1, [alt], gt="./.", gq=14, ad=[8, 8], filt="NoCall", qual=0.2
        ),
        block("chr1", 3004, 3400),
    ]
    called = _calls(tmp_path, records, [deletion])[deletion]
    assert (called.dosage, called.row_type) == ("unknown", RowType.NOCALL)
    assert called.warnings == (NOCALL_REASON,)


def test_no_coverage_unknown(tmp_path: Path) -> None:
    uncovered, absent, other_contig = (
        _snv("chr1", 2050),
        _snv("chr1", 3500),
        _snv("chrX", 3500),
    )
    records = [block("chr1", 2001, 2100, gq=1, min_dp=0)]
    calls = _calls(tmp_path, records, [uncovered, absent, other_contig])
    assert calls[uncovered] == Called(
        gt=None,
        dosage="unknown",
        row_type=RowType.NONE,
        gq=None,
        dp=0,
        ad=None,
        vaf=None,
        warnings=(NO_COVERAGE_REASON,),
    )
    for claim in (absent, other_contig):
        assert calls[claim].dosage == "unknown"
        assert calls[claim].row_type == RowType.NONE
        assert calls[claim].warnings == (NO_COVERAGE_REASON,)


def test_other_alt_dosage_zero_with_warning(tmp_path: Path) -> None:
    # `A C,G,<*>` at 0/1: G's reads did not make a genotype (observed at 0/1, GQ 13).
    pos = 2001
    base = ref("chr1", pos)
    c, g, t = (other_base(base, k) for k in range(3))
    in_record, not_in_record, called_alt = (
        ("chr1", pos, base, g),
        ("chr1", pos, base, t),
        ("chr1", pos, base, c),
    )
    deletion = ("chr1", pos, ref("chr1", pos, 2), base)
    records = [
        block("chr1", 1900, 2000),
        variant("chr1", pos, 1, [c, g], gt="0/1", gq=43, ad=[10, 10, 6]),
        block("chr1", 2002, 2100),
    ]
    calls = _calls(tmp_path, records, [in_record, not_in_record, called_alt, deletion])
    for claim in (in_record, not_in_record, deletion):
        assert (calls[claim].dosage, calls[claim].row_type) == ("0", RowType.PASS)
        assert calls[claim].warnings == (OTHER_ALLELE_WARNING,)
    assert (calls[called_alt].dosage, calls[called_alt].warnings) == ("1", ())
    # GQ 43 is high for an SNV; for an indel the inference is capped at medium.
    assert _grade(calls[in_record], in_record).level == Level.HIGH
    assert _grade(calls[deletion], deletion, indel=True).level == Level.MEDIUM


def test_multiallelic_record_counts_claim_allele(tmp_path: Path) -> None:
    pos = 2001
    base = ref("chr1", pos)
    first, second = other_base(base, 0), other_base(base, 1)
    snv_pos = 2501  # `TGG AGG,T,<*>`: an SNV and a deletion of one position, 1/2
    snv_ref = ref("chr1", snv_pos, 3)
    snv_alt = other_base(snv_ref[0])
    claims = [
        ("chr1", pos, base, second),
        ("chr1", pos, base, first),
        ("chr1", snv_pos, snv_ref[0], snv_alt),
        ("chr1", snv_pos, snv_ref, snv_ref[0]),
    ]
    records = [
        variant("chr1", pos, 1, [first, second], gt="1/2", gq=76, ad=[0, 20, 19]),
        variant(
            "chr1",
            snv_pos,
            3,
            [snv_alt + snv_ref[1:], snv_ref[0]],
            gt="1/2",
            gq=73,
            ad=[0, 21, 16],
        ),
    ]
    calls = _calls(tmp_path, records, claims)
    for claim in claims:
        assert (calls[claim].gt, calls[claim].dosage, calls[claim].warnings) == (
            "1/2",
            "1",
            (),
        )
    assert calls[claims[0]].vaf == pytest.approx(19 / 39, abs=1e-5)


def test_inside_deletion_unknown(tmp_path: Path) -> None:
    # `CTGCGT C,<*>` at 0/1, an SNV on the other haplotype at its third deleted base: a
    # record of its own; the next block starts after the deletion's REF.
    pos = 3000
    inside, own_record = _snv("chr1", pos + 2), _snv("chr1", pos + 3)
    star_pos = 3501
    star = _snv("chr1", star_pos, 2)
    star_alt = other_base(ref("chr1", star_pos), 0)
    # A called multi-base substitution (`TGG AGG`, alone): its second base is no deletion,
    # and is not worked out either.
    mnv_pos = 3701
    in_mnv = _snv("chr1", mnv_pos + 1)
    mnv_ref = ref("chr1", mnv_pos, 3)
    records = [
        block("chr1", 2900, pos - 1),
        variant("chr1", pos, 6, [ref("chr1", pos)], gt="0/1", gq=30, ad=[10, 20]),
        variant("chr1", pos + 3, 1, [own_record[3]], gt="0/1", gq=14, ad=[0, 9]),
        block("chr1", pos + 6, 3400),
        star_record("chr1", star_pos, star_alt, gt="0/2"),
        variant(
            "chr1",
            mnv_pos,
            3,
            [other_base(mnv_ref[0]) + mnv_ref[1:]],
            gt="0/1",
            gq=40,
            ad=[12, 12],
        ),
    ]
    calls = _calls(tmp_path, records, [inside, own_record, star, in_mnv])
    assert (calls[in_mnv].dosage, calls[in_mnv].warnings) == (
        "unknown",
        (MULTI_BASE_REASON,),
    )
    assert calls[inside].dosage == "unknown"
    assert calls[inside].warnings == (INSIDE_DELETION_REASON,)
    assert calls[inside].row_type == RowType.NONE
    assert (calls[own_record].dosage, calls[own_record].row_type) == ("1", RowType.PASS)
    assert calls[star].dosage == "unknown"
    assert calls[star].warnings == (INSIDE_DELETION_REASON,)


def test_haploid_one_copy(tmp_path: Path) -> None:
    hom, in_block, refcall = _snv("chrX", 2001), _snv("chrX", 2500), _snv("chrX", 3001)
    in_par = _snv("chrX", 501)
    records = [
        block("chrX", 400, 500),
        variant("chrX", 501, 1, [in_par[3]], gt="1/1", gq=45, ad=[0, 14]),
        variant("chrX", 2001, 1, [hom[3]], gt="1/1", gq=45, ad=[0, 14]),
        block("chrX", 2002, 3000, min_dp=8),
        variant(
            "chrX",
            3001,
            1,
            [refcall[3]],
            gt="0/0",
            gq=32,
            ad=[10, 2],
            filt="RefCall",
            qual=0,
        ),
    ]
    calls = _calls(tmp_path, records, [hom, in_block, refcall, in_par], haploid=male)
    assert (calls[hom].gt, calls[hom].dosage) == ("1", "1")
    # DeepVariant writes a male's non-PAR reference blocks 0/0: one copy of the reference.
    assert (calls[in_block].gt, calls[in_block].dosage) == ("0", "0")
    assert calls[in_block].row_type == RowType.REF_BLOCK
    assert (calls[refcall].gt, calls[refcall].dosage) == ("0", "0")
    assert (calls[in_par].gt, calls[in_par].dosage) == ("1/1", "2")


def test_haploid_het_low(tmp_path: Path) -> None:
    claim = _snv("chrX", 2001)
    records = [variant("chrX", 2001, 1, [claim[3]], gt="0/1", gq=40, ad=[7, 7])]
    called = _calls(tmp_path, records, [claim], haploid=male)[claim]
    assert (called.gt, called.dosage) == ("0/1", "1")
    graded = _grade(called, claim, haploid=True)
    assert graded.level == Level.LOW
    assert HAPLOID_HET in graded.warnings
    assert OUTSIDE_CALIBRATED in graded.warnings  # chrX: outside HG003's benchmark
    assert _grade(called, claim, haploid=False).level == Level.HIGH


def test_alleles_compared_after_normalisation(tmp_path: Path) -> None:
    # `CAA C,CA,<*>` at 1/2 in a homopolymer (observed): its CA allele is CA>C.
    contig, pos, run = HOMOPOLYMER
    one_a = (contig, pos, run[:2], run[0])
    two_a = (contig, pos, run[:3], run[0])
    one_a_shifted = (
        contig,
        pos + 7,
        run[7:9],
        run[7],
    )  # the same deletion, written right
    records = [
        block(contig, 900, pos - 1),
        variant(contig, pos, 3, [run[0], run[:2]], gt="1/2", gq=78, ad=[2, 14, 8]),
        block(contig, pos + 3, 1100),
    ]
    calls = _calls(tmp_path, records, [one_a, two_a, one_a_shifted])
    for claim in (one_a, two_a, one_a_shifted):
        assert (calls[claim].dosage, calls[claim].warnings) == ("1", ()), claim
    assert calls[one_a].vaf == pytest.approx(8 / 24, abs=1e-5)
