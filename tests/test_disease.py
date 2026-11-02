"""ClinVar's diseases (CLNDISDB) matched to ClinGen's curations by MONDO ID."""

from genokeeper.evidence.clingen import Curation
from genokeeper.evidence.disease import CurationMatch, match_curation, parse_clndisdb


def test_parse_clndisdb_keeps_one_set_of_mondo_ids_per_disease() -> None:
    assert parse_clndisdb("MedGen:C1|MONDO:MONDO:0011450,OMIM:604370") == [
        frozenset(),
        frozenset({"MONDO:0011450"}),
    ]


def test_parse_clndisdb_of_no_value() -> None:
    assert parse_clndisdb(".") == [frozenset()]


def _diseases(*mondo: str) -> list[frozenset[str]]:
    return [frozenset({m}) for m in mondo]


CURATIONS = {
    "GENE": [
        Curation("MONDO:0000001", "Moderate"),
        Curation("MONDO:0000002", "Definitive"),
        Curation("MONDO:0000003", "Disputed"),
        Curation("MONDO:0000004", "Limited"),
        Curation("MONDO:0000005", "Refuted"),
    ]
}


def test_strongest_matched_positive_wins() -> None:
    found = match_curation(
        "GENE", _diseases("MONDO:0000001", "MONDO:0000002"), CURATIONS, "Definitive"
    )
    assert found == CurationMatch("Definitive", matched=True)
    # A matched positive beats a matched negative.
    found = match_curation(
        "GENE", _diseases("MONDO:0000003", "MONDO:0000001"), CURATIONS, "Definitive"
    )
    assert found == CurationMatch("Moderate", matched=True)


def test_matched_negative_only_without_positive() -> None:
    found = match_curation("GENE", _diseases("MONDO:0000003"), CURATIONS, "Definitive")
    assert found == CurationMatch("Disputed", matched=True)
    found = match_curation(
        "GENE", _diseases("MONDO:0000003", "MONDO:0000005"), CURATIONS, "Definitive"
    )
    assert found == CurationMatch("Refuted", matched=True)


def test_unmatched_diseases_ignored_when_one_matches() -> None:
    # MONDO:0000009 has no curation; the gene's Refuted curation is for another disease.
    found = match_curation(
        "GENE", _diseases("MONDO:0000009", "MONDO:0000001"), CURATIONS, "Definitive"
    )
    assert found == CurationMatch("Moderate", matched=True)


def test_no_match_falls_back_to_gene_level() -> None:
    found = match_curation("GENE", _diseases("MONDO:0000009"), CURATIONS, "Definitive")
    assert found == CurationMatch("Definitive", matched=False)
    assert match_curation("GENE", [frozenset()], CURATIONS, "Definitive") == (
        CurationMatch("Definitive", matched=False)
    )
    assert match_curation("OTHER", _diseases("MONDO:0000001"), CURATIONS, None) == (
        CurationMatch(None, matched=False)
    )


def test_a_curation_without_a_mondo_id_matches_nothing() -> None:
    curations = {"GENE": [Curation("", "Definitive")]}
    assert match_curation("GENE", [frozenset()], curations, "Definitive") == (
        CurationMatch("Definitive", matched=False)
    )


def test_limited_is_not_positive() -> None:
    found = match_curation(
        "GENE", _diseases("MONDO:0000004", "MONDO:0000003"), CURATIONS, "Definitive"
    )
    assert found == CurationMatch("Disputed", matched=True)
    found = match_curation("GENE", _diseases("MONDO:0000004"), CURATIONS, "Definitive")
    assert found == CurationMatch("Limited", matched=True)
