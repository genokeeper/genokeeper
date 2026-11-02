"""CGD inheritance (the Clinical Genomic Database's own download) for genes on neither list."""

from pathlib import Path

import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.genome_mini import write_cgd

from genokeeper.config import load_grading
from genokeeper.evidence.cgd import load_cgd
from genokeeper.evidence.clinvar import INHERITANCE_REASON
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim, ClaimSet
from genokeeper.sources import RunContext
from genokeeper.sources.panel import ClaimsSource
from genokeeper.store.db import Kind, Level, RowType


def test_load_cgd_maps_single_modes_only(tmp_path: Path) -> None:
    table = write_cgd(
        tmp_path / "CGD.txt.gz",
        [
            ("GENEX", "AR"),
            # Several modes: left out.
            ("GENEY", "AD/AR"),
            ("GENEZ", "AD"),
            ("GENEZ", "AR"),  # rows that disagree: left out
            ("GENEW", " XL "),
            ("GENEV", "Mitochondrial"),
            ("GENEU", "AD"),
            ("GENEU", "AD"),
        ],
    )
    assert load_cgd(table) == {"GENEX": "AR", "GENEW": "XL", "GENEU": "AD"}


def test_load_cgd_refuses_another_table(tmp_path: Path) -> None:
    other = tmp_path / "other.txt.gz"
    other.write_bytes(b"\x1f\x8b" + b"not gzip")
    with pytest.raises(ValueError, match="CGD"):
        load_cgd(other)


def _claim(variation_id: int, gene: str, pos: int) -> Claim:
    return Claim(
        variation_id=variation_id,
        gene=gene,
        clnsig="Pathogenic",
        clnrevstat="criteria_provided,_multiple_submitters,_no_conflicts",
        disease_ids="",
        condition="Some disease",
        chrom="chr1",
        pos=pos,
        ref="A",
        alt="G",
        vid=f"GRCh38:chr1:{pos}:A:G",
        rsid=None,
    )


def test_cgd_inheritance_for_unlisted_gene(tmp_path: Path) -> None:
    """A gene on neither list takes its inheritance from CGD: one pathogenic copy in a
    CGD-recessive gene is no finding, while a gene CGD lacks stays "inheritance unknown"."""
    inheritance = load_cgd(write_cgd(tmp_path / "CGD.txt.gz", [("GENEX", "AR")]))
    claims = [_claim(1, "GENEX", 100), _claim(2, "GENEY", 200)]
    het = Called(
        gt="0/1",
        dosage="1",
        row_type=RowType.PASS,
        gq=50,
        dp=30,
        ad=(15, 15, 0),
        vaf=0.5,
        warnings=(),
    )
    graded = Reliability(p_wrong=1e-5, level=Level.HIGH, flags=(), warnings=())
    ctx = RunContext(
        run_id="20261020T000000Z-abcdef",
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance=inheritance,
        male=True,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
    )
    source = ClaimsSource(
        ClaimSet(claims=claims, input_rows=2, dropped={}),
        [(c, het, graded) for c in claims],
    )
    rows = {r.statement["statement_id"]: r for r in source.collect(ctx).rows}
    recessive, unknown = rows["clinvar:1"].statement, rows["clinvar:2"].statement
    assert recessive["kind"] == Kind.NON_FINDING
    assert unknown["kind"] == Kind.FINDING
    assert INHERITANCE_REASON in str(unknown["reason"])
    # Schema 3's rows: the variant keyed by its GRCh38 vid alone, the call as called.
    variant, call = rows["clinvar:1"].variant, rows["clinvar:1"].call
    assert variant is not None and call is not None
    assert set(variant) == {"run_id", "vid", "rsid"}
    assert call["row_type"] == "PASS" and "mapping" not in call
