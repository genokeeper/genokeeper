"""Per-gene coverage and callability over MANE CDS (slice 3c)."""

from importlib.resources import files
from pathlib import Path

from genokeeper.config import load_grading
from genokeeper.coverage import (
    REFERENCE_LIMITED,
    CoverageRow,
    cds_bed,
    coverage_fields,
    coverage_genes,
    coverage_rows,
    expected_copies,
    load_cds,
    share_inside,
)
from genokeeper.deletions.exons import MANE_RELEASE, MANE_SHA256, Exon, cut_exons
from genokeeper.genelists import load_gene_lists
from genokeeper.regions import read_bed


def test_resource_covers_the_genes() -> None:
    genes = coverage_genes(load_gene_lists(), load_grading().repeats.loci)
    assert {"SMN1", "SMN2", "HTT", "FXN", *REFERENCE_LIMITED} <= genes
    cds = load_cds()
    # SLC37A4's only MANE transcript is on chr11_KN196481v1_fix, off the primary assembly.
    assert set(cds) == genes - {"SLC37A4"}
    assert {"KCNE1", "CRYAA", "U2AF1", "CBS"} <= set(cds)
    assert all(cds[g] for g in cds)
    header = files("genokeeper").joinpath("resources", "coverage", "cds.GRCh38.tsv")
    lines = header.read_text().splitlines()
    assert f"# release: {MANE_RELEASE}" in lines
    assert f"# sha256: {MANE_SHA256}" in lines


def _line(feature: str, start: int, end: int, number: int) -> str:
    attributes = (
        'gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_name "AAA"; '
        f'exon_number {number}; tag "MANE_Select"; '
    )
    return f"chr1\thavana\t{feature}\t{start}\t{end}\t.\t+\t.\t{attributes}\n"


def test_cut_reads_cds_rows_only(tmp_path: Path) -> None:
    gtf = tmp_path / "mane.gtf"
    gtf.write_text(_line("exon", 901, 1100, 1) + _line("CDS", 1001, 1100, 1))
    assert cut_exons(gtf, {"AAA"}, feature="CDS") == [
        (Exon("AAA", "ENST1.1", 1, "chr1", 1000, 1100), "MANE_Select")
    ]
    # 3b's exon cut is unchanged.
    assert cut_exons(gtf, {"AAA"}) == [
        (Exon("AAA", "ENST1.1", 1, "chr1", 900, 1100), "MANE_Select")
    ]


def _exon(gene: str, chrom: str, start: int, end: int) -> Exon:
    return Exon(gene, "ENST1.1", 1, chrom, start, end)


def test_cds_bed_margin_merge_clip(tmp_path: Path) -> None:
    cds = {"AAA": [_exon("AAA", "chr1", 1103, 1200), _exon("AAA", "chr1", 1000, 1100)]}
    out = tmp_path / "cds.bed"
    genes = cds_bed(cds, margin=2, contigs={"chr1": 5000}, out=out)
    assert genes == {"AAA": [("chr1", 998, 1202)]}
    assert out.read_text() == "chr1\t998\t1202\tAAA\n"
    clipped = cds_bed(
        {**cds, "BBB": [_exon("BBB", "chr9", 10, 20)]},
        margin=2,
        contigs={"chr1": 1150},
        out=out,
    )
    assert clipped == {"AAA": [("chr1", 998, 1150)]}  # BBB's contig is absent: dropped
    assert read_bed(out) == [("chr1", 998, 1150)]


def test_share_inside() -> None:
    assert share_inside([], [("chr1", 0, 10)]) == 0.0
    intervals = [("chr1", 0, 10), ("chr2", 0, 10)]
    assert share_inside(intervals, [("chr1", 5, 20), ("chr1", 0, 2)]) == 0.35


def test_coverage_rows_values() -> None:
    genes = {
        "CBS": [("chr21", 1000, 2000)],
        "PMS2": [("chr7", 0, 1000)],
        "SMN1": [("chr5", 0, 100)],
    }
    rows = coverage_rows(
        {"CBS": (34.0, 0.99), "PMS2": (21.0, 0.80), "SMN1": (30.0, 0.97)},
        genes,
        difficult=[("chr7", 500, 1000)],
        benchmark=[("chr21", 0, 5000)],
        owners={"SMN1": "smn"},
        min_callable=0.95,
    )
    by_gene = {r.gene: r for r in rows}
    assert by_gene["CBS"] == CoverageRow(
        gene="CBS",
        bases=1000,
        callable=0.99,
        mean_depth=34.0,
        difficult=0.0,
        outside_benchmark=0.0,
        owner=None,
        reliable=True,
    )
    pms2 = by_gene["PMS2"]
    assert (pms2.reliable, pms2.difficult, pms2.outside_benchmark) == (False, 0.5, 1.0)
    assert {g: r.owner for g, r in by_gene.items()} == {
        "CBS": None,
        "PMS2": None,
        "SMN1": "smn",
    }
    fields = coverage_fields("r1", by_gene["SMN1"])
    assert fields == {
        "run_id": "r1",
        "gene": "SMN1",
        "bases": 100,
        "callable": 0.97,
        "mean_depth": 30.0,
        "difficult": 0.0,
        "outside_benchmark": 1.0,
        "owner": "smn",
        "reliable": True,
    }


def test_callable_per_expected_copy() -> None:
    par = load_grading().par
    genes = {
        "DMD": [("chrX", 31_000_000, 33_000_000)],
        "CBS": [("chr21", 43_000_000, 43_100_000)],
        "SHOX": [("chrX", 600_000, 700_000)],  # inside PAR1
        "SRY": [("chrY", 2_786_000, 2_787_000)],
    }
    assert expected_copies(genes, male=True, par=par) == {
        "DMD": 1,
        "CBS": 2,
        "SHOX": 2,
        "SRY": 1,
    }
    assert set(expected_copies(genes, male=False, par=par).values()) == {2}
    # A male's DMD at a mean of 18x with 99% of its bases at 10x (its one copy's 10x).
    (dmd,) = coverage_rows(
        {"DMD": (18.0, 0.99)},
        {"DMD": genes["DMD"]},
        difficult=[],
        benchmark=[],
        owners={},
        min_callable=load_grading().coverage.min_callable,
    )
    assert dmd.reliable
