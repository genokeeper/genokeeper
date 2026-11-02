"""GenCC beside ClinGen for gene-disease validity (slice 3f, "GenCC"): the new-format export
read by column name, its classification titles in ClinGen's words, its submissions
aggregated per gene and disease (a dispute caps the strongest), and only for genes ClinGen
does not curate; a reason naming a GenCC validity says so.

`fixtures/gencc_head.tsv` is the export's header and its first 50 rows (2026-10-30).
"""

from pathlib import Path

from fixtures.faf_mini import write_faf_table

from genokeeper.config import load_grading
from genokeeper.evidence.clingen import Curation, gene_level
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.gencc import (
    GENCC_FILE,
    aggregate,
    gencc_date,
    load_gencc,
    merge,
    newest_gencc,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.panel.build import Claim
from genokeeper.sources import RunContext
from genokeeper.sources.panel import grade_claim
from genokeeper.store.db import RowType, Validity

HEAD = Path(__file__).parent / "fixtures" / "gencc_head.tsv"
TWO_STARS = "criteria_provided,_multiple_submitters,_no_conflicts"


def gencc(*classifications: str, mondo: str = "MONDO:1") -> list[Curation]:
    return [Curation(mondo, c) for c in classifications]


def test_load_gencc_head() -> None:
    """Every row is one submission, keyed by gene; titles read as ClinGen's words
    (Supportive is Limited, neither is positive)."""
    loaded = load_gencc(HEAD)
    assert sum(map(len, loaded.values())) == 50
    rims1 = loaded["RIMS1"]
    assert [(c.mondo, c.classification) for c in rims1[:5]] == [
        ("MONDO:0011355", "Limited"),
        ("MONDO:0005258", "Limited"),
        ("MONDO:0011355", "Disputed"),
        ("MONDO:0011355", "Strong"),
        ("MONDO:0015993", "Limited"),  # Orphanet's Supportive
    ]
    assert {c.classification for c in loaded["MYH8"]} == {
        "Limited",
        "Strong",
        "Refuted",
    }
    assert [c.classification for c in loaded["TREH"]] == ["Limited"]
    words = {c.classification for found in loaded.values() for c in found}
    assert words <= {
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed",
        "Refuted",
    }
    assert gencc_date(HEAD) == "20261025"  # its newest submitted_run_date (KIF1B)


def test_gencc_disputed_caps() -> None:
    """Per disease the strongest classification, unless any submitter disputes or refutes
    it: one lab's Definitive never outweighs another's dispute."""
    assert aggregate(gencc("Definitive", "Disputed")) == gencc("Disputed")
    assert aggregate(gencc("Strong", "Limited", "Strong")) == gencc("Strong")
    assert aggregate(gencc("Refuted")) == gencc("Disputed")
    assert aggregate(
        gencc("Limited", mondo="MONDO:2") + gencc("Moderate", "Refuted")
    ) == gencc("Limited", mondo="MONDO:2") + gencc("Disputed")
    loaded = load_gencc(HEAD)
    rims1 = {c.mondo: c.classification for c in aggregate(loaded["RIMS1"])}
    assert rims1 == {
        "MONDO:0011355": "Disputed",  # Strong (Labcorp) against Disputed (Illumina)
        "MONDO:0005258": "Limited",
        "MONDO:0015993": "Limited",
        "MONDO:0005260": "Strong",
    }
    myh8 = {c.mondo: c.classification for c in aggregate(loaded["MYH8"])}
    assert myh8 == {"MONDO:0008016": "Strong", "MONDO:0012137": "Disputed"}


def test_clingen_wins_over_gencc() -> None:
    """A gene ClinGen curates keeps ClinGen's curations only; the others take GenCC's,
    aggregated, and are named as GenCC's."""
    clingen = {
        "SKI": [Curation("MONDO:0008426", "Moderate")],
        "KIF1B": [Curation("MONDO:0007308", "Limited")],
    }
    curations, sources = merge(clingen, load_gencc(HEAD))
    assert curations["SKI"] == clingen["SKI"]
    assert curations["KIF1B"] == clingen["KIF1B"]
    assert curations["SLC45A1"] == [
        Curation("MONDO:0044322", "Strong"),
        Curation("MONDO:0019502", "Limited"),
    ]
    assert sources["SLC45A1"] == sources["TREH"] == "GenCC"
    assert "SKI" not in sources and "KIF1B" not in sources


def _context(
    tmp_path: Path, curations: dict[str, list[Curation]], sources: dict[str, str]
) -> RunContext:
    return RunContext(
        run_id="run-1",
        lists=load_gene_lists(),
        grading=load_grading(),
        curations=curations,
        gene_level={gene: gene_level(found) for gene, found in curations.items()},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=False,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
        validity_source=sources,
    )


def _claim(gene: str, mondo: str) -> Claim:
    return Claim(
        variation_id=1,
        gene=gene,
        clnsig="Pathogenic",
        clnrevstat=TWO_STARS,
        disease_ids=f"MONDO:{mondo}",
        condition="some disease",
        chrom="chr1",
        pos=100,
        ref="C",
        alt="T",
        vid="GRCh38:chr1:100:C:T",
        rsid=None,
    )


def test_gencc_fills_uncurated_gene(tmp_path: Path) -> None:
    """A 2-star pathogenic claim in a GenCC-Strong gene ClinGen lacks is strong, as under a
    ClinGen Strong; in a GenCC-Limited one it is limited, its reason naming GenCC."""
    curations, sources = merge({}, load_gencc(HEAD))
    ctx = _context(tmp_path, curations, sources)
    het = Called("0/1", "1", RowType.PASS, 40, 30, None, None, ())
    strong = grade_claim(_claim("SLC45A1", "MONDO:0044322"), het, ctx, None)
    assert strong.level is Validity.STRONG
    limited = grade_claim(_claim("TREH", "MONDO:0012803"), het, ctx, None)
    assert limited.level is Validity.LIMITED
    assert limited.reason is not None
    assert limited.reason.split("; ")[0] == "gene validity Limited (GenCC)"
    # No GenCC source given: the reason reads as before.
    plain = _context(tmp_path, curations, {})
    reason = grade_claim(_claim("TREH", "MONDO:0012803"), het, plain, None).reason
    assert reason is not None and reason.split("; ")[0] == "gene validity Limited"


def test_newest_gencc(tmp_path: Path) -> None:
    assert newest_gencc(tmp_path) is None
    for day in ("20261025", "20261101"):
        (tmp_path / "gencc" / day).mkdir(parents=True)
    (tmp_path / "gencc" / "20261025" / GENCC_FILE).write_text("x")
    assert newest_gencc(tmp_path) == tmp_path / "gencc" / "20261025" / GENCC_FILE
