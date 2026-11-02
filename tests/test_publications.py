"""LitVar2's variants file and ClinVar's var_citations.txt, read from small cuts of the real
public files (tests/fixtures/litvar_head.jsonl.gz, var_citations_head.txt: the first lines of
each as served on 2026-10-27, LitVar2's records picked by line)."""

import gzip
import json
from pathlib import Path

import pytest

from genokeeper.evidence.publications import (
    CITATIONS_FILE,
    LITVAR_FILE,
    check_citations,
    check_litvar,
    cited_pmids,
    litvar_pmids,
)

FIXTURES = Path(__file__).parent / "fixtures"
LITVAR = FIXTURES / "litvar_head.jsonl.gz"
CITATIONS = FIXTURES / "var_citations_head.txt"


def _litvar(path: Path, *records: dict[str, object]) -> Path:
    path.write_bytes(
        gzip.compress(b"".join(json.dumps(r).encode() + b"\n" for r in records))
    )
    return path


def test_file_names_are_the_bootstraps() -> None:
    assert (LITVAR_FILE, CITATIONS_FILE) == (
        "litvar2_variants.json.gz",
        "var_citations.txt",
    )


def test_litvar_merges_records_sharing_an_rsid(tmp_path: Path) -> None:
    # Real records: the dbSNP record and the ClinGen-allele record of each rsID.
    got = litvar_pmids(LITVAR, {"rs397507801", "rs746985154"})
    assert got == {
        "rs397507801": frozenset({"15172753", "12552570", "32377563", "36922933"}),
        "rs746985154": frozenset({"27480553", "27397505"}),
    }
    # The union, not the first or the largest record: each record has a PMID the other lacks.
    path = _litvar(
        tmp_path / "two.json.gz",
        {"_id": "litvar@rs1##", "rsid": "rs1", "pmids": [1, 2]},
        {"_id": "litvar@CA9#rs1##", "rsid": "rs1", "pmids": [2, 3]},
    )
    assert litvar_pmids(path, {"rs1"}) == {"rs1": frozenset({"1", "2", "3"})}


def test_litvar_keeps_only_asked_rsids() -> None:
    assert litvar_pmids(LITVAR, {"rs746985154", "rs999"}) == {
        "rs746985154": frozenset({"27480553", "27397505"})
    }
    assert litvar_pmids(LITVAR, set()) == {}
    # Records that name no rsID (a protein change, "All NR2E3 variants") answer to none.
    assert litvar_pmids(LITVAR, {"p.S44X", "rs10002"}) == {}


def test_cited_pmids_by_variation_id_pubmed_only(tmp_path: Path) -> None:
    got = cited_pmids(CITATIONS, {"2", "5", "404"})
    assert got["2"] == frozenset({"20613862", "25741868", "26467025", "40081374"})
    assert "4544753" not in got["2"]  # a PubMedCentral number, not a PMID
    assert set(got) == {"2", "5"}  # asked and cited; VariationID 3 and 404 not asked
    assert len(got["5"]) == 6
    # Variation 9 has a bookshelf row (NBK1440) among its PubMed rows.
    assert not {p for p in cited_pmids(CITATIONS, {"9"})["9"] if not p.isdigit()}
    assert "NBK1440" not in cited_pmids(CITATIONS, {"9"})["9"]
    # A VariationID cited only by other sources is absent, not an empty set.
    only_pmc = tmp_path / CITATIONS_FILE
    only_pmc.write_text(
        "#AlleleID\tVariationID\trs\tnsv\tcitation_source\tcitation_id\torganization_ids\n"
        "1\t7\t\t\tPubMedCentral\t123\t1\n"
        "1\t7\t\t\tNCBIBookShelf\tNBK1\t1\n"
    )
    assert cited_pmids(only_pmc, {"7"}) == {}
    assert cited_pmids(CITATIONS, set()) == {}


def test_checks_accept_the_real_cuts() -> None:
    check_litvar(LITVAR)
    check_citations(CITATIONS)


def test_checks_refuse_a_file_without_the_fields(tmp_path: Path) -> None:
    good = {"_id": "litvar@rs1##", "rsid": "rs1", "pmids": [1]}
    refused = {
        "html": tmp_path / "html.gz",
        "no pmids": _litvar(tmp_path / "np.gz", {"_id": "x", "rsid": "rs1"}),
        "rsid not an rsID": _litvar(tmp_path / "rs.gz", {**good, "rsid": "1234"}),
        "no rsID in the file": _litvar(tmp_path / "nr.gz", {**good, "rsid": None}),
        "a PMID that is not digits": _litvar(
            tmp_path / "pm.gz", {**good, "pmids": ["x"]}
        ),
        "empty": _litvar(tmp_path / "empty.gz"),
    }
    refused["html"].write_bytes(gzip.compress(b"<html>maintenance</html>\n"))
    for path in refused.values():
        with pytest.raises(ValueError):
            check_litvar(path)
    # Not gzip at all (an error page served as the file), and truncated gzip.
    plain = tmp_path / "plain.gz"
    plain.write_text("<html>maintenance</html>")
    with pytest.raises(ValueError, match="not a readable gzip"):
        check_litvar(plain)
    cut = tmp_path / "cut.gz"
    cut.write_bytes(LITVAR.read_bytes()[:40])
    with pytest.raises(ValueError):
        check_litvar(cut)
    # A line whose spacing the reader's fast pattern would miss is refused here.
    spaced = tmp_path / "spaced.gz"
    spaced.write_bytes(gzip.compress(b'{"_id":"x","rsid" :"rs1","pmids":[1]}\n'))
    with pytest.raises(ValueError, match="would not find"):
        check_litvar(spaced)

    header = "#AlleleID\tVariationID\trs\tnsv\tcitation_source\tcitation_id\torganization_ids\n"
    row = "1\t7\t\t\tPubMed\t123\t1\n"
    texts = {
        "html": "<html>maintenance</html>\n",
        "no citation_id column": header.replace("citation_id", "x") + row,
        "no rows": header,
        "non-numeric VariationID": header + row.replace("\t7\t", "\tabc\t"),
        "a PubMed id that is not digits": header + row.replace("123", "NBK1"),
        "short row": header + "1\t7\n",
        "empty": "",
    }
    for name, text in texts.items():
        path = tmp_path / f"{name.replace(' ', '_')}.txt"
        path.write_text(text)
        with pytest.raises(ValueError):
            check_citations(path)
