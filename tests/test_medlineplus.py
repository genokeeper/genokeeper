"""MedlinePlus Genetics summaries: a gene's function by its exact symbol, a condition's
description by its name or a synonym after normalisation, and nothing fuzzier."""

from pathlib import Path

import pytest
from fixtures.medlineplus_mini import condition, gene, write_medlineplus

from genokeeper.evidence.medlineplus import Summary, load_medlineplus

HEMOCHROMATOSIS = condition(
    "hereditary-hemochromatosis",
    "Hereditary hemochromatosis",
    "Hereditary hemochromatosis is a disorder.",
    "There are four types.",
    synonyms=("Iron overload disorder",),
)


def test_gene_function_text_by_symbol(tmp_path: Path) -> None:
    path = write_medlineplus(
        tmp_path / "ghr.xml",
        gene("HFE", "The <html:i>HFE</html:i> gene makes a protein.", "It binds."),
        HEMOCHROMATOSIS,
    )
    found = load_medlineplus(path)
    assert found.genes == {
        "HFE": Summary(
            id="medlineplus:gene:HFE",
            name="hfe protein",
            text="The HFE gene makes a protein.\n\nIt binds.",
            published="2025-08-25",
        )
    }
    # A gene named only inside a condition's related-gene list is no summary.
    assert "RELATED1" not in found.genes


def test_condition_by_name_and_synonym(tmp_path: Path) -> None:
    found = load_medlineplus(
        write_medlineplus(tmp_path / "ghr.xml", gene("HFE", "x"), HEMOCHROMATOSIS)
    )
    summary = Summary(
        id="medlineplus:condition:hereditary-hemochromatosis",
        name="Hereditary hemochromatosis",
        text="Hereditary hemochromatosis is a disorder.\n\nThere are four types.",
        published="2020-09-08",
    )
    assert found.condition("Hereditary_hemochromatosis") == summary
    assert found.condition("  HEREDITARY   hemochromatosis ") == summary
    assert found.condition("iron_overload_disorder") == summary
    # A gene summary's related-condition names and synonyms are not conditions.
    assert found.condition("HFE related disorder") is None
    assert found.condition("HFE disorder") is None


def test_no_fuzzy_match(tmp_path: Path) -> None:
    found = load_medlineplus(
        write_medlineplus(
            tmp_path / "ghr.xml",
            gene("HFE", "x"),
            HEMOCHROMATOSIS,
            # Two conditions share a synonym: it names neither.
            condition("a-syndrome", "A syndrome", "a", synonyms=("AS",)),
            condition("b-syndrome", "B syndrome", "b", synonyms=("AS",)),
        )
    )
    assert found.condition("hemochromatosis type 1") is None
    assert found.condition("hemochromatosis") is None
    assert found.condition("AS") is None
    assert found.condition("A syndrome") is not None


@pytest.mark.parametrize(
    "text",
    [
        "<html><body>maintenance</body></html>",
        # The right root, but no gene summary.
        '<summaries xmlns="https://medlineplus.gov/download/x.xsd"></summaries>',
        "not xml at all",
    ],
)
def test_not_a_summaries_file_refused(tmp_path: Path, text: str) -> None:
    path = tmp_path / "ghr.xml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=str(path)):
        load_medlineplus(path)
