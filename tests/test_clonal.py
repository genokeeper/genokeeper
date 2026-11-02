"""The clonal genes' spans: the packaged resource the allele-fraction guard reads."""

from genokeeper.config import load_grading
from genokeeper.genotype.clonal import ClonalIndex, load_spans

NINE = ("TP53", "CHEK2", "DNMT3A", "TET2", "ASXL1", "PPM1D", "JAK2", "SF3B1", "SRSF2")


def test_spans_resource_holds_the_nine() -> None:
    spans = load_spans()
    assert set(spans) == set(NINE)
    assert set(load_grading().af_guard.clonal_genes) == set(NINE)
    assert spans["TP53"].chrom == "chr17"
    index = ClonalIndex.load(NINE)
    assert index.at("chr17", 7_675_000)
    assert index.at("17", 7_675_000)  # contig names are normalised
    assert not index.at("chr17", 1_000)
    # Only the genes asked for: TP53's span is not clonal for an index without it.
    assert not ClonalIndex.load(("JAK2",)).at("chr17", 7_675_000)
