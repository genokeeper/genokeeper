"""The genome-wide claims of a ClinVar GRCh38 release (`panel.build.build_claims`)."""

from pathlib import Path

import pysam
import pytest
from fixtures import clinvar_mini as cv
from fixtures import gvcf_mini as gm

from genokeeper.genelists import GeneList
from genokeeper.panel.build import (
    CLAIM_DROP_REASONS,
    Claim,
    build_claims,
    check_clinvar38,
)

REVIEWED = "criteria_provided,_multiple_submitters,_no_conflicts"
CLAIM_LISTS = {
    "acmg_sf": GeneList(name="acmg_sf", version="0", source="t", genes={"GENEA": "AD"}),
    "carrier_tier3": GeneList(
        name="carrier_tier3", version="0", source="t", genes={"GENEB": "AR"}
    ),
}


def _cv(gene: str, clnsig: str | None, **extra: str) -> dict[str, str]:
    info = {"CLNREVSTAT": REVIEWED, "GENEINFO": gene, **extra}
    if clnsig is not None:
        info["CLNSIG"] = clnsig
    return info


def _snv(pos: int, contig: str = "chr1") -> tuple[str, str]:
    base = gm.ref(contig, pos)
    return base, gm.other_base(base)


def test_claims_genome_wide_classes(tmp_path: Path) -> None:
    """Every class the spec evaluates is a claim in any gene; benign records, VUS outside
    the lists, other assertions and records the reference refutes are counted, not kept."""
    fasta = gm.write_fasta(tmp_path / "ref.fa")
    homopolymer = gm.HOMOPOLYMER[1]  # C, eight A, G
    wrong_ref = gm.other_base(gm.ref("chr1", 1100))
    pathogenic = _cv("GENEA:3", "Pathogenic")
    insertion = gm.ref("chr1", 1800) + "ACGT" * 12 + "ACG"  # 51 bases after the anchor
    records: list[cv.Record] = [
        ("1", 100, 7001, *_snv(100), _cv("UNLISTED:1", "Pathogenic", RS="701")),
        (
            "1",
            200,
            7002,
            *_snv(200),
            _cv("UNLISTED:1", "Conflicting_classifications_of_pathogenicity"),
        ),
        ("1", 300, 7003, *_snv(300), _cv("UNLISTED:1", "drug_response")),
        ("1", 400, 7004, *_snv(400), _cv("UNLISTED:1", "risk_factor")),
        ("1", 500, 7005, *_snv(500), _cv("UNLISTED:1", "Pathogenic|Benign")),
        ("1", 600, 7006, *_snv(600), _cv("GENEB:2", "Uncertain_significance")),
        ("1", 700, 7007, *_snv(700), _cv("UNLISTED:1", "Uncertain_significance")),
        ("1", 800, 7008, *_snv(800), _cv("GENEA:3", "Likely_benign")),
        ("1", 900, 7009, *_snv(900), _cv("UNLISTED:1", "not_provided")),
        ("1", 950, 7010, *_snv(950), _cv("UNLISTED:1", None)),
        ("1", 1100, 7011, wrong_ref, "T", _cv("GENEA:3", "Pathogenic")),
        ("1", 1200, 7012, gm.ref("chr1", 1200), "A,C", _cv("GENEA:3", "Pathogenic")),
        ("1", 1300, 7013, gm.ref("chr1", 1300), ".", _cv("GENEA:3", "Pathogenic")),
        ("NT_187361.1", 500, 7014, "A", "G", _cv("GENEA:3", "Pathogenic")),
        ("MT", 3243, 7015, "A", "G", _cv("MT-TL1:4567", "Pathogenic")),
        # A homopolymer deletion written at its right end: CA>C once left-aligned.
        ("1", homopolymer + 7, 7016, "AA", "A", _cv("GENEA:3", "Pathogenic")),
        ("X", 2000, 7017, *_snv(2000, "chrX"), _cv("UNLISTED:1|GENEA:3", "Pathogenic")),
        # Deletions of 50 and 51 bases and an insertion of 51: only the 50 is a
        # claim, the others are structural variants SNV and small-indel calling cannot see.
        ("1", 1500, 7018, gm.ref("chr1", 1500, 51), gm.ref("chr1", 1500), pathogenic),
        ("1", 1600, 7019, gm.ref("chr1", 1600, 52), gm.ref("chr1", 1600), pathogenic),
        ("1", 1800, 7020, gm.ref("chr1", 1800), insertion, pathogenic),
    ]
    clinvar = cv.write_clinvar(tmp_path / "clinvar38.vcf.gz", "GRCh38", records)

    found = build_claims(clinvar, CLAIM_LISTS, fasta=fasta)

    kept = {c.variation_id: c for c in found.claims}
    assert sorted(kept) == [7001, 7002, 7003, 7004, 7005, 7006, 7016, 7017, 7018]
    assert {k: v for k, v in found.dropped.items() if v} == {
        "vus_outside_lists": 1,
        "benign": 1,
        "not_in_scope": 2,
        "ref_disagrees_with_reference": 1,
        "multiallelic": 1,
        "no_alt": 1,
        "non_primary_contig": 1,
        "contig_not_in_reference": 1,
        "structural_variant": 2,
    }
    assert set(found.dropped) == set(CLAIM_DROP_REASONS)
    assert found.input_rows == len(records)
    assert found.input_rows == len(found.claims) + sum(found.dropped.values())

    ref, alt = _snv(100)
    assert kept[7001] == Claim(
        variation_id=7001,
        gene="UNLISTED",
        clnsig="Pathogenic",
        clnrevstat=REVIEWED,
        disease_ids="",
        condition="",
        chrom="chr1",
        pos=100,
        ref=ref,
        alt=alt,
        vid=f"GRCh38:chr1:100:{ref}:{alt}",
        rsid="rs701",
    )
    deletion = kept[7016]
    assert (deletion.chrom, deletion.pos, deletion.ref, deletion.alt) == (
        "chr1",
        homopolymer,
        "CA",
        "C",
    )
    assert deletion.vid == f"GRCh38:chr1:{homopolymer}:CA:C"
    # The listed gene wins among several symbols.
    assert kept[7017].gene == "GENEA" and kept[7017].chrom == "chrX"


def test_check_clinvar38_takes_what_the_run_reads_and_refuses_the_rest(
    tmp_path: Path,
) -> None:
    """The bootstrap's pre-filing check: a release the reader reads and holding a claim
    passes; one with no classified record, or a record with no numeric ID, does not."""
    good = cv.write_clinvar(
        tmp_path / "good.vcf.gz",
        "GRCh38",
        [("1", 100, 1, "A", "G", _cv("GENEA:3", "Pathogenic"))],
    )
    check_clinvar38(good)
    unclassified = cv.write_clinvar(
        tmp_path / "none.vcf.gz",
        "GRCh38",
        [("1", 100, 1, "A", "G", {"GENEINFO": "GENEA:3"})],
    )
    with pytest.raises(ValueError, match="no record with a germline classification"):
        check_clinvar38(unclassified)
    # A record without an ID: the reader stops on it, so the check does.
    plain = tmp_path / "noid.vcf"
    text = str(pysam.VariantFile(str(good)).header)
    plain.write_text(
        text + "1\t100\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENEA:3\n"
    )
    noid = tmp_path / "noid.vcf.gz"
    pysam.tabix_compress(str(plain), str(noid), force=True)
    with pytest.raises(ValueError):
        check_clinvar38(noid)
