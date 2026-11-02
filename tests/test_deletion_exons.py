"""The MANE exons of the listed genes: the cut, the packaged resource, and the screen's
targets."""

import csv
from importlib.resources import files
from pathlib import Path

import pytest

from genokeeper.deletions.exons import (
    MANE_RELEASE,
    MANE_SHA256,
    OWNED_CYP2D,
    PARALOG_SPAN_GENES,
    PRIMARY,
    SCREEN_GENES_NOTE,
    Exon,
    Target,
    cut_exons,
    cut_spans,
    listed_genes,
    load_exons,
    load_spans,
    locus_span,
    targets,
)
from genokeeper.genelists import load_gene_lists, load_screen_genes

RESOURCES = files("genokeeper").joinpath("resources", "deletions")


def _exon_line(
    chrom: str,
    start: int,
    end: int,
    gene: str,
    transcript: str,
    number: int,
    tags: tuple[str, ...],
) -> str:
    attributes = (
        f'gene_id "ENSG{gene}.1"; transcript_id "{transcript}.1"; gene_name "{gene}"; '
        f"exon_number {number}; " + "".join(f'tag "{t}"; ' for t in tags)
    )
    return f"{chrom}\thavana\texon\t{start}\t{end}\t.\t+\t.\t{attributes}\n"


def test_resource_covers_lists_and_screen_genes() -> None:
    """ACMG SF, carrier Tier 3 and the screen genes (DPYD, CHEK2): 192 genes, 191 in the
    resource. DPYD's MANE Select (MANE v1.5: NM_000110.4, ENST00000370192.8) has 23 exons,
    read at the cut; CHEK2's holds exons 9 and 10, the del5395 founder deletion's."""
    lists = load_gene_lists()
    listed = (
        set(lists["acmg_sf"].genes)
        | set(lists["carrier_tier3"].genes)
        | set(load_screen_genes())
    )
    assert len(listed) == 192
    assert listed_genes() == listed
    exons = load_exons()
    assert len({e.number for e in exons["DPYD"]}) == 23
    assert {e.transcript for e in exons["DPYD"]} == {"ENST00000370192.8"}
    assert {9, 10} <= {e.number for e in exons["CHEK2"]}
    # SLC37A4's only MANE transcript is on chr11_KN196481v1_fix, off the primary assembly.
    assert set(exons) == listed - {"SLC37A4"}
    assert all(exons[g] for g in exons)
    assert {e.chrom for found in exons.values() for e in found} <= set(PRIMARY)
    # The tag column of the resource: every gene has a MANE Select transcript's exons.
    rows = [
        r
        for r in csv.DictReader(
            (
                line
                for line in RESOURCES.joinpath("exons.GRCh38.tsv")
                .read_text()
                .splitlines()
                if not line.startswith("#")
            ),
            delimiter="\t",
        )
    ]
    assert {r["gene"] for r in rows if r["tag"] == "MANE_Select"} == set(exons)


def test_cut_keeps_mane_only(tmp_path: Path) -> None:
    gtf = tmp_path / "mane.gtf"
    gtf.write_text(
        "#!genome-build GRCh38\n"
        + _exon_line("chr1", 1001, 1100, "AAA", "ENST1", 1, ("MANE_Select", "basic"))
        + _exon_line("chr1", 2001, 2100, "AAA", "ENST1", 2, ("MANE_Select",))
        + _exon_line(
            "chr1", 1501, 1600, "AAA", "ENST2", 1, ("basic", "MANE_Plus_Clinical")
        )
        + _exon_line("chr1", 1701, 1800, "AAA", "ENST3", 1, ("basic",))
        + _exon_line("chr2", 501, 600, "BBB", "ENST4", 1, ("MANE_Select",))
        + 'chr1\thavana\tCDS\t1001\t1100\t.\t+\t0\tgene_name "AAA"; tag "MANE_Select";\n'
    )
    assert cut_exons(gtf, {"AAA"}) == [
        (Exon("AAA", "ENST1.1", 1, "chr1", 1000, 1100), "MANE_Select"),
        (Exon("AAA", "ENST2.1", 1, "chr1", 1500, 1600), "MANE_Plus_Clinical"),
        (Exon("AAA", "ENST1.1", 2, "chr1", 2000, 2100), "MANE_Select"),
    ]


def test_cut_keeps_primary_contigs_only(tmp_path: Path) -> None:
    """An exon on a fix or alt contig is dropped: no caller shard runs there."""
    gtf = tmp_path / "mane.gtf"
    gtf.write_text(
        _exon_line("chr11", 1001, 1100, "AAA", "ENST1", 1, ("MANE_Select",))
        + _exon_line(
            "chr11_KN196481v1_fix", 2001, 2100, "BBB", "ENST2", 1, ("MANE_Select",)
        )
        + _exon_line("chrY", 501, 600, "CCC", "ENST3", 1, ("MANE_Select",))
    )
    assert cut_exons(gtf, {"AAA", "BBB", "CCC"}) == [
        (Exon("AAA", "ENST1.1", 1, "chr11", 1000, 1100), "MANE_Select"),
        (Exon("CCC", "ENST3.1", 1, "chrY", 500, 600), "MANE_Select"),
    ]


def test_targets_apply_margin() -> None:
    exon = Exon("AAA", "ENST1", 1, "chr1", 1000, 1100)
    assert targets({"AAA": [exon]}, 20) == {"AAA": [Target("chr1", 980, 1120, (1,))]}


def test_targets_merge_close_exons() -> None:
    def exon(number: int, start: int, end: int) -> Exon:
        return Exon("AAA", "ENST1", number, "chr1", start, end)

    first = exon(1, 1000, 1100)
    close = targets({"AAA": [first, exon(2, 1130, 1200), exon(3, 1230, 1300)]}, 20)
    assert close == {"AAA": [Target("chr1", 980, 1320, (1, 2, 3))]}
    apart = targets({"AAA": [first, exon(2, 1150, 1200)]}, 20)
    assert apart == {
        "AAA": [Target("chr1", 980, 1120, (1,)), Target("chr1", 1130, 1220, (2,))]
    }


def test_source_names_release_and_sha() -> None:
    assert MANE_RELEASE == "MANE v1.5"
    assert MANE_SHA256.startswith("71d4b3c89c9d948683bf0db5d81b00d9")
    header = [
        line
        for line in RESOURCES.joinpath("exons.GRCh38.tsv").read_text().splitlines()
        if line.startswith("#")
    ]
    assert f"# release: {MANE_RELEASE}" in header
    assert f"# sha256: {MANE_SHA256}" in header
    assert any(
        line.startswith("# source: https://ftp.ncbi.nlm.nih.gov/") for line in header
    )
    source = RESOURCES.joinpath("SOURCE.md").read_text()
    assert MANE_RELEASE in source and MANE_SHA256 in source
    assert any(
        SCREEN_GENES_NOTE in line for line in header if line.startswith("# cut:")
    )


def _transcript_line(
    chrom: str, start: int, end: int, gene: str, tags: tuple[str, ...]
) -> str:
    attributes = f'gene_name "{gene}"; ' + "".join(f'tag "{t}"; ' for t in tags)
    return f"{chrom}\thavana\ttranscript\t{start}\t{end}\t.\t+\t.\t{attributes}\n"


def test_cut_spans_reads_mane_select_on_primary_contigs(tmp_path: Path) -> None:
    """A MANE Select transcript's span, 0-based half-open; a Plus Clinical or untagged
    transcript is not read, nor an alt contig's (MANE v1.5 places a C4A on
    chr6_GL000254v2_alt, which must not replace chr6's)."""
    gtf = tmp_path / "mane.gtf"
    gtf.write_text(
        _transcript_line("chr6", 1001, 2000, "C4A", ("MANE_Select",))
        + _transcript_line("chr6_GL000254v2_alt", 5001, 6000, "C4A", ("MANE_Select",))
        + _transcript_line("chr6", 901, 2500, "C4A", ("MANE_Plus_Clinical",))
        + _transcript_line("chr1", 101, 200, "AAA", ("basic",))
        + _exon_line("chr1", 101, 200, "BBB", "ENST2", 1, ("MANE_Select",))
    )
    assert cut_spans(gtf, {"C4A", "AAA", "BBB"}) == {"C4A": ("chr6", 1000, 2000)}


def test_spans_cover_paralog_genes() -> None:
    """Every paralog gene has its MANE Select span on its chromosome; GSTM1's is pinned as
    the MANE v1.5 GTF gives it (chr1:109,687,817-109,693,745), and GBA1_ex9-11 is GBA1's
    exons 9 to 11, the part GBAP1 copies."""
    spans = load_spans()
    expected = {
        "STRC": "chr15", "CATSPER2": "chr15", "RHD": "chr1", "RHCE": "chr1",
        "UGT2B17": "chr4", "UGT2B15": "chr4", "AMY1A": "chr1", "AMY1B": "chr1",
        "AMY1C": "chr1", "GSTM1": "chr1", "C4A": "chr6", "C4B": "chr6",
        "CYP2A6": "chr19", "CYP2A7": "chr19", "GBA1": "chr1", "LPA": "chr6",
        "CYP2D6": "chr22", "SMN1": "chr5", "SMN2": "chr5",
    }  # fmt: skip
    assert set(PARALOG_SPAN_GENES) == set(expected)
    assert {g: spans[g][0] for g in PARALOG_SPAN_GENES} == expected
    assert spans["GSTM1"] == ("chr1", 109687816, 109693745)
    assert spans["C4A"] == ("chr6", 31982056, 32002681)
    assert spans["GBA1_ex9-11"] == ("chr1", 155234451, 155235844)
    gba1 = spans["GBA1"]
    assert gba1[1] == spans["GBA1_ex9-11"][1] and gba1[2] > spans["GBA1_ex9-11"][2]


def test_locus_span_unites_one_chromosome() -> None:
    spans = load_spans()
    amy1 = locus_span(spans, ("AMY1A", "AMY1B", "AMY1C"))
    assert amy1 == ("chr1", spans["AMY1A"][1], spans["AMY1C"][2])
    with pytest.raises(ValueError, match="chromosome"):
        locus_span(spans, ("SMN1", "CYP2D6"))


def test_owned_cyp2d_covers_pseudogenes() -> None:
    """CYP2D6's span, 5 kb before and 25 kb after it, holds CYP2D7 and CYP2D8P (NCBI Gene
    1564 and 1568, GRCh38: chr22:42,139,576-42,144,483 and chr22:42,149,883-42,155,016,
    1-based, read 2026-10-30)."""
    chrom, start, end = OWNED_CYP2D
    cyp2d6 = load_spans()["CYP2D6"]
    assert (chrom, start, end) == ("chr22", cyp2d6[1] - 5_000, cyp2d6[2] + 25_000)
    for first, last in ((42_139_576, 42_144_483), (42_149_883, 42_155_016)):
        assert start <= first - 1 and last <= end
