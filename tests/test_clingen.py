"""ClinGen's gene-disease validity export, read from a small synthetic CSV built here."""

from pathlib import Path

import pytest
from fixtures.clingen_mini import (
    write_dosage,
    write_gene_validity,
    write_recurrent,
    write_regions,
)

from genokeeper.evidence.clingen import (
    DOSAGE_FILE,
    POSITIVE_VALIDITY,
    RECURRENT_FILE,
    REGIONS_FILE,
    Curation,
    DosageList,
    gene_validity_date,
    load_curations,
    load_dosage,
    load_gene_validity,
    load_recurrent,
    load_regions,
    regions_date,
)


def test_dosage_list_names_its_haploinsufficient_genes(tmp_path: Path) -> None:
    """Score 3 is haploinsufficient; 30 (autosomal recessive), 40 (unlikely), 1 and 2 are
    not. The list's date files it; a page without the header is refused."""
    path = write_dosage(
        tmp_path / DOSAGE_FILE,
        [("BRCA1", "3"), ("CFTR", "30"), ("GENEW", "1"), ("GENEZ", "40")],
    )
    assert load_dosage(path) == DosageList(
        "20261023",
        frozenset({"BRCA1"}),
        frozenset(),
        {g: ("chr1", 0, 2) for g in ("BRCA1", "CFTR", "GENEW", "GENEZ")},
    )
    page = tmp_path / "page.tsv"
    page.write_text("#ClinGen Gene Curation Results\n#23 Oct,2026\n<html>gone</html>\n")
    with pytest.raises(ValueError, match="header"):
        load_dosage(page)


def test_dosage_triplosensitive_and_spans(tmp_path: Path) -> None:
    """Triplosensitivity score 3 is triplosensitive (slice 3d); each gene's location,
    1-based inclusive in the file, is read 0-based half-open, and `tbd` gives no span."""
    path = write_dosage(
        tmp_path / DOSAGE_FILE,
        [("MECP2", "30"), ("PMP22", "3"), ("GENEX", "0")],
        triplo={"MECP2": "3", "PMP22": "3"},
        locations={"MECP2": "chrX:154021573-154097717", "GENEX": "tbd"},
    )
    dosage = load_dosage(path)
    assert dosage.triplosensitive == {"MECP2", "PMP22"}
    assert dosage.haploinsufficient == {"PMP22"}
    assert dosage.spans == {
        "MECP2": ("chrX", 154021572, 154097717),
        "PMP22": ("chr1", 0, 2),
    }


def test_regions_keep_scored_and_recurrent(tmp_path: Path) -> None:
    """Region curations with a haploinsufficiency or triplosensitivity score of 3 and a
    location; the recurrent BED's orange rows, not its black breakpoint clusters. A page
    that is neither is refused."""
    tsv = write_regions(
        tmp_path / REGIONS_FILE,
        [
            ("22q11.2 recurrent (DGS) region (proximal, A-D)", "chr22:18924719-21111383", "3", "1"),
            ("21q21.3 region (includes APP)", "chr21:25734570-26573286", "0", "3"),
            ("16p12.3 recurrent region", "chr16:16762939-18073544", "0", "0"),
            ("19p13.3 Critical Region 1", "tbd", "3", "0"),
            ("22q11.21 population region", "chr22:21512587-21554998", "40", "0"),
        ],
    )  # fmt: skip
    assert load_regions(tsv) == [
        ("chr22", 18924718, 21111383, "22q11.2 recurrent (DGS) region (proximal, A-D)"),
        ("chr21", 25734569, 26573286, "21q21.3 region (includes APP)"),
    ]
    assert regions_date(tsv) == "20260428"
    bed = write_recurrent(
        tmp_path / RECURRENT_FILE,
        [
            ("chr22", 18912231, 18924718, "BP A", "0,0,0"),
            ("chr22", 18924718, 21111383, "22q11.2 recurrent (DGS) region (A-D)", "255,111,0"),
        ],
    )  # fmt: skip
    assert load_recurrent(bed) == [
        ("chr22", 18924718, 21111383, "22q11.2 recurrent (DGS) region (A-D)")
    ]
    page = tmp_path / "page.txt"
    page.write_text("<html>maintenance</html>\n")
    with pytest.raises(ValueError):
        load_regions(page)
    with pytest.raises(ValueError, match="BED"):
        load_recurrent(page)
    black = write_recurrent(tmp_path / "black.bed", [("chr1", 1, 2, "BP1", "0,0,0")])
    with pytest.raises(ValueError, match="recurrent"):
        load_recurrent(black)


def test_gene_validity_takes_the_strongest_positive(tmp_path: Path) -> None:
    csv = write_gene_validity(
        tmp_path / "gene-validity.csv",
        [
            ("BRCA1", "Moderate"),
            ("BRCA1", "Definitive"),
            ("BRCA1", "Disputed"),
            ("GENEL", "Refuted"),
            ("GENEL", "Limited"),
            ("GENEP", "Limited"),
            ("GENEP", "Definitive"),
            ("GENED", "Disputed"),
            ("GENED", "No Known Disease Relationship"),
            ("GENER", "Disputed"),
            ("GENER", "Refuted"),
            ("GENEN", "No Known Disease Relationship"),
        ],
    )
    assert load_gene_validity(csv) == {
        "BRCA1": "Definitive",
        "GENEL": "Refuted",
        "GENEP": "Definitive",
        "GENED": "Disputed",
        "GENER": "Refuted",
        "GENEN": "No Known Disease Relationship",
    }


def test_gene_validity_date_is_the_files_own(tmp_path: Path) -> None:
    csv = write_gene_validity(tmp_path / "g.csv", [("BRCA1", "Definitive")])
    assert gene_validity_date(csv) == "20261018"


def test_gene_validity_reads_the_columns_by_name(tmp_path: Path) -> None:
    # A column added before CLASSIFICATION moves it; the header still names it.
    csv = write_gene_validity(
        tmp_path / "g.csv", [("BRCA1", "Strong")], extra_column=True
    )
    assert load_gene_validity(csv) == {"BRCA1": "Strong"}


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ('"GENE ID (HGNC)","CLASSIFICATION"\n', "GENE SYMBOL"),
        ('"FILE CREATED: 2026-10-18"\n"GENE SYMBOL","MOI"\n', "CLASSIFICATION"),
        ('"FILE CREATED: 2026-10-18"\n"GENE SYMBOL","CLASSIFICATION"\n', "no curation"),
    ],
)
def test_gene_validity_refuses_an_unexpected_file(
    tmp_path: Path, text: str, match: str
) -> None:
    csv = tmp_path / "g.csv"
    csv.write_text(text)
    with pytest.raises(ValueError, match=match):
        load_gene_validity(csv)


def test_gene_validity_date_refuses_a_file_without_one(tmp_path: Path) -> None:
    csv = tmp_path / "g.csv"
    csv.write_text('"GENE SYMBOL","CLASSIFICATION"\n"BRCA1","Definitive"\n')
    with pytest.raises(ValueError, match="FILE CREATED"):
        gene_validity_date(csv)


def test_curations_keep_each_disease(tmp_path: Path) -> None:
    csv = write_gene_validity(
        tmp_path / "g.csv",
        [
            ("BRCA1", "Definitive", "MONDO:0011450"),
            ("BRCA1", "Disputed Evidence", "MONDO:0007254"),
            ("GENEL", "Limited", "MONDO:0000777"),
        ],
    )
    assert load_curations(csv) == {
        "BRCA1": [
            Curation("MONDO:0011450", "Definitive"),
            Curation("MONDO:0007254", "Disputed"),
        ],
        "GENEL": [Curation("MONDO:0000777", "Limited")],
    }


def test_curations_without_a_mondo_column_match_nothing(tmp_path: Path) -> None:
    csv = write_gene_validity(
        tmp_path / "g.csv", [("BRCA1", "Definitive")], mondo_column=False
    )
    assert load_curations(csv) == {"BRCA1": [Curation("", "Definitive")]}
    assert load_gene_validity(csv) == {"BRCA1": "Definitive"}


def test_positive_validity_leaves_limited_out() -> None:
    assert POSITIVE_VALIDITY == {"definitive", "strong", "moderate"}
