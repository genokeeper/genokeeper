import gzip
from pathlib import Path

from genokeeper.config import load_grading
from genokeeper.regions import (
    AUTOSOMES,
    CHRX_AMPLICONIC,
    CHRX_NONPAR,
    CHRX_XTR,
    CHRY_AMPLICONIC,
    CHRY_NONPAR,
    CHRY_XTR,
    GAPS,
    LOW_MAPPABILITY,
    QC_AUTOSOMES,
    SEX_REGIONS,
    derive,
    giab38_downloads,
    subtract,
)


def test_subtract_cuts_overlaps_and_keeps_the_rest() -> None:
    keep = [("chr1", 0, 100), ("chr1", 200, 300), ("chr2", 0, 50)]
    remove = [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 90, 210), ("chr2", 50, 60)]
    assert subtract(keep, remove) == [
        ("chr1", 0, 10),
        ("chr1", 30, 90),
        ("chr1", 210, 300),
        ("chr2", 0, 50),
    ]


def _bed(path: Path, rows: list[tuple[str, int, int]]) -> None:
    path.write_bytes(
        gzip.compress("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows).encode())
    )


def _chrx(folder: Path, nonpar: tuple[int, int] = (0, 1000)) -> None:
    _bed(folder / CHRX_NONPAR, [("chrX", *nonpar)])
    _bed(folder / CHRX_XTR, [])
    _bed(folder / CHRX_AMPLICONIC, [])


def test_xdr_is_nonpar_without_xtr_ampliconic_gaps_or_low_mappability(
    tmp_path: Path,
) -> None:
    _bed(tmp_path / CHRY_NONPAR, [("chrY", 1000, 9000)])
    _bed(tmp_path / CHRY_XTR, [("chrY", 2000, 3000)])
    _bed(tmp_path / CHRY_AMPLICONIC, [("chrY", 5000, 6000)])
    _bed(tmp_path / GAPS, [("chrY", 8500, 9500), ("chr1", 0, 100)])
    _bed(tmp_path / LOW_MAPPABILITY, [("chrY", 1000, 1100), ("chr1", 400, 500)])
    _bed(tmp_path / AUTOSOMES, [("chr1", 0, 1000)])
    _chrx(tmp_path)
    assert derive(tmp_path) == [tmp_path / SEX_REGIONS, tmp_path / QC_AUTOSOMES]
    assert (tmp_path / SEX_REGIONS).read_text().splitlines() == [
        "chr1\t100\t400\tautosomes",
        "chr1\t500\t1000\tautosomes",
        "chrY\t1100\t2000\txdr",
        "chrY\t3000\t5000\txdr",
        "chrY\t6000\t8500\txdr",
        "chrX\t0\t1000\txnonpar",
    ]
    # The QC gate's autosomes keep low-mappability sequence: only gaps are not genome.
    assert (tmp_path / QC_AUTOSOMES).read_text().splitlines() == [
        "chr1\t100\t1000\tautosomes"
    ]


def test_downloads_cover_the_strata_and_the_sources() -> None:
    names = giab38_downloads(load_grading())
    assert set(load_grading().strata.values()) <= set(names)
    assert {CHRY_NONPAR, CHRY_XTR, CHRY_AMPLICONIC, AUTOSOMES, GAPS} <= set(names)
    assert {CHRX_NONPAR, CHRX_XTR, CHRX_AMPLICONIC} <= set(names)
    assert len(names) == len(set(names))


def test_derive_writes_xnonpar(tmp_path: Path) -> None:
    _bed(tmp_path / CHRY_NONPAR, [("chrY", 1000, 9000)])
    _bed(tmp_path / CHRY_XTR, [("chrY", 2000, 3000)])
    _bed(tmp_path / CHRY_AMPLICONIC, [])
    _bed(tmp_path / AUTOSOMES, [("chr1", 0, 1000)])
    _bed(tmp_path / CHRX_NONPAR, [("chrX", 0, 1000)])
    _bed(tmp_path / CHRX_XTR, [("chrX", 100, 200)])
    _bed(tmp_path / CHRX_AMPLICONIC, [("chrX", 300, 400)])
    _bed(tmp_path / GAPS, [("chrX", 500, 550)])
    _bed(tmp_path / LOW_MAPPABILITY, [("chrX", 700, 750)])
    derive(tmp_path)
    rows = [r.split("\t") for r in (tmp_path / SEX_REGIONS).read_text().splitlines()]
    assert [(c, int(s), int(e)) for c, s, e, name in rows if name == "xnonpar"] == [
        ("chrX", 0, 100),
        ("chrX", 200, 300),
        ("chrX", 400, 500),
        ("chrX", 550, 700),
        ("chrX", 750, 1000),
    ]
    # xdr is chrY's alone: chrX's XTR, ampliconic, gaps and low mappability leave it.
    assert [(c, int(s), int(e)) for c, s, e, name in rows if name == "xdr"] == [
        ("chrY", 1000, 2000),
        ("chrY", 3000, 9000),
    ]
