"""Sex-chromosome dosage (slice 3c): the X and Y copy counts, the chrX allele balance, the
karyotype they suggest, and the run's sex from the X copies."""

import json
from dataclasses import replace
from pathlib import Path

import pysam
import pytest
from fixtures.run_context import run_context

from genokeeper.config import load_grading
from genokeeper.sex import (
    Dosage,
    SexRefused,
    copies,
    karyotype,
    sex_from_ratios,
    x_het_share,
)
from genokeeper.sources.dosage import DOSAGE_TOOL, DosageSource
from genokeeper.store.db import Kind, Topic, Validity

XNONPAR = [("chrX", 3_000_000, 4_000_000), ("chrX", 5_000_000, 6_000_000)]
VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chrX,length=156040895>\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


def _vcf(path: Path, rows: list[tuple[int, str, tuple[int, int], int | None]]) -> Path:
    """A haploid chrX VCF of (pos, filter, AD, DP) records (DP defaults to AD's sum)."""
    body = "".join(
        f"chrX\t{pos}\t.\tA\tG\t30\t{flt}\t.\tGT:DP:AD\t1:{sum(ad) if dp is None else dp}"
        f":{ad[0]},{ad[1]}\n"
        for pos, flt, ad, dp in rows
    )
    plain = path.with_suffix("")
    plain.write_text(VCF_HEADER + body)
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def test_x_het_share_reads_ad_not_gt(tmp_path: Path) -> None:
    vcf = _vcf(
        tmp_path / "calls.vcf.gz",
        [
            (100_000, "PASS", (15, 15), None),  # PAR1: out
            (3_000_100, "PASS", (15, 16), None),  # balanced: a heterozygous share
            (3_000_200, "PASS", (0, 30), None),
            (3_000_300, "PASS", (2, 28), None),
            (3_000_400, "RefCall", (14, 13), None),  # RefCall counts, and is balanced
            (3_000_500, "PASS", (3, 2), 5),  # DP 5: below [reliability] min_dp
            (4_500_000, "PASS", (15, 15), None),  # between the intervals (XTR): out
        ],
    )
    shipped = load_grading()

    def at(min_x_sites: int) -> tuple[float | None, int]:
        dosage = replace(shipped.dosage, min_x_sites=min_x_sites)
        grading = replace(shipped, dosage=dosage)
        return x_het_share(vcf, xnonpar=XNONPAR, grading=grading)

    assert at(4) == (2 / 4, 4)
    assert at(5) == (None, 4)


@pytest.mark.parametrize(
    ("ratio", "expected"), [(0.52, 1), (0.40, 1), (0.74, None), (1.03, 2), (0.01, 0)]
)
def test_copies(ratio: float, expected: int | None) -> None:
    # 0.40: a ratio below 0.5 but within the tolerance stays one copy.
    assert copies(ratio, 0.15) == expected


def test_sex_from_x_copies() -> None:
    t = load_grading().dosage
    assert sex_from_ratios(0.50, 0.25, t) == "male"  # a high Y ratio does not change it
    assert sex_from_ratios(0.50, 0.10, t) == "male"  # nor a low Y above y_absent_max
    assert sex_from_ratios(0.50, 0.0, t) == "male"
    assert sex_from_ratios(1.0, 0.0, t) == "female"
    with pytest.raises(SexRefused, match="--sex") as refused:
        sex_from_ratios(0.75, 0.0, t)
    assert "0.750" in str(refused.value) and "dosage statement" in str(refused.value)


@pytest.mark.parametrize(
    ("x", "y", "het", "label", "kind", "topic"),
    [
        (1.0, 0.49, 0.6, "47,XXY", Kind.FINDING, Topic.RISKS),
        (0.5, 1.02, 0.0, "47,XYY", Kind.FINDING, Topic.RISKS),
        (0.5, 0.01, 0.0, "45,X", Kind.FINDING, Topic.RISKS),
        (0.5, 0.10, 0.0, "possible mosaic loss of Y", Kind.FINDING, Topic.RESEARCH),
        (0.5, 0.25, 0.0, "possible mosaic loss of Y", Kind.FINDING, Topic.RESEARCH),
        (
            0.75,
            0.0,
            0.3,
            "possible mosaic sex-chromosome pattern",
            Kind.FINDING,
            Topic.RISKS,
        ),
        (
            1.0,
            0.75,
            0.3,
            "possible mosaic sex-chromosome pattern",
            Kind.FINDING,
            Topic.RISKS,
        ),
        (0.51, 0.5, 0.4, "possible 46,XY/47,XXY mosaic", Kind.FINDING, Topic.RISKS),
        (0.5, 0.5, 0.02, "46,XY", Kind.NON_FINDING, Topic.RISKS),
        (1.0, 0.0, 0.5, "46,XX", Kind.NON_FINDING, Topic.RISKS),
        (1.5, 0.0, 0.5, "47,XXX", Kind.FINDING, Topic.RISKS),
    ],
)
def test_karyotype_table(
    x: float, y: float, het: float, label: str, kind: Kind, topic: Topic
) -> None:
    t = load_grading().dosage
    k = karyotype(Dosage(x_ratio=x, y_ratio=y, x_het_share=het, x_sites=2000), t)
    assert (k.label, k.kind, k.topic) == (label, kind, topic)
    if label == "possible mosaic loss of Y":
        assert k.validity == Validity.LIMITED
        share = "80%" if y == 0.10 else "50%"
        assert k.reason == (
            f"estimated share of cells without Y: {share} (information; not actionable)"
        )
    elif kind == Kind.FINDING:
        # A doubtful pattern, worded "possible", is limited; a whole count is moderate.
        possible = label.startswith("possible")
        assert k.validity == (Validity.LIMITED if possible else Validity.MODERATE)
        assert k.reason.startswith("screening result, needs a karyotype")
    if label == "46,XY":
        assert k.reason == "XY, no sign of an extra X or Y"
    if label == "46,XX":
        assert k.reason == "XX, no sign of a missing or extra X"


def _statement(source: DosageSource, tmp_path: Path) -> dict[str, object]:
    collected = source.collect(run_context(tmp_path))
    assert (collected.tool, collected.input_rows, collected.dropped) == (
        DOSAGE_TOOL,
        1,
        {},
    )
    (row,) = collected.rows
    assert row.variant is None and row.call is None
    return row.statement


def test_dosage_source_not_run_unassessed(tmp_path: Path) -> None:
    """`not_run` (a failed step, or the foundation table's reason) leaves it unassessed,
    whatever `calibrated_on` says."""
    t = replace(load_grading().dosage, calibrated_on="20261030")
    dosage = Dosage(x_ratio=0.5, y_ratio=0.5, x_het_share=0.01, x_sites=2000)
    for measured in (dosage, None):
        statement = _statement(DosageSource(measured, "stale: x", t), tmp_path)
        assert (statement["kind"], statement["reason"]) == (Kind.UNASSESSED, "stale: x")
        assert statement["actionable"] is False


def test_dosage_source_finding_wording(tmp_path: Path) -> None:
    t = replace(load_grading().dosage, calibrated_on="20261030")
    dosage = Dosage(x_ratio=1.0, y_ratio=0.49, x_het_share=0.6, x_sites=2000)
    statement = _statement(DosageSource(dosage, None, t), tmp_path)
    assert statement["statement_id"] == "dosage:sex-chromosomes"
    assert (statement["kind"], statement["topic"]) == (Kind.FINDING, Topic.RISKS)
    assert str(statement["reason"]).startswith("screening result, needs a karyotype")
    assert statement["actionable"] is False
    assert statement["gene"] is None and statement["vid"] is None
    assert (statement["source"], statement["source_version"]) == ("dosage", "20261030")
    assert statement["source_level"] == "mosdepth depth ratios and chrX allele balance"
    assert json.loads(str(statement["claim"])) == {
        "schema": "genomelab.claim/1",
        "type": "sex_chromosome_dosage",
        "label": "47,XXY",
        "x_ratio": 1.0,
        "y_ratio": 0.49,
        "x_het_share": 0.6,
        "x_sites": 2000,
    }
    normal = _statement(DosageSource(Dosage(0.5, 0.5, 0.01, 2000), None, t), tmp_path)
    assert (normal["kind"], normal["reason"]) == (
        Kind.NON_FINDING,
        "XY, no sign of an extra X or Y",
    )
