"""`actionable` on a ClinVar statement in an ACMG SF gene: true only for a finding ACMG SF
v3.3 Table 1 would report. Every claim here is made up; the genes are the list's."""

import pytest

from genokeeper.config import load_grading
from genokeeper.evidence.clinvar import grade_clinvar
from genokeeper.genelists import load_gene_lists, parse_gene_list
from genokeeper.sources.panel import clinvar_statement

LISTS = load_gene_lists()
HFE_C282Y = "GRCh38:chr6:26092913:G:A"
OTHER_VID = "GRCh38:chr1:1000:A:G"
NONSENSE = "SO:0001587|nonsense"
MISSENSE = "SO:0001583|missense_variant"


def _actionable(
    gene: str,
    dosage: str,
    clnsig: str = "Pathogenic",
    *,
    vid: str = OTHER_VID,
    consequence: str = "",
    male: bool = True,
    haploid: bool = False,
) -> bool:
    revstat = "criteria_provided,_multiple_submitters,_no_conflicts"
    grade = grade_clinvar(
        clnsig=clnsig,
        clnrevstat=revstat,
        gene=gene,
        gene_validity="Definitive",
        lists=LISTS,
        faf_grpmax=None,
        dosage=dosage,
        grading=load_grading(),
        male=male,
        haploid=haploid,
        vid=vid,
    )
    row = clinvar_statement(
        "r",
        variation_id=1,
        vid=vid,
        gene=gene,
        clnsig=clnsig,
        clnrevstat=revstat,
        condition="not provided",
        consequence=consequence,
        grade=grade,
        source_version="20261011",
        lists=LISTS,
    )
    return bool(row["actionable"])


def test_dominant_finding_actionable() -> None:
    assert _actionable("BRCA1", "1")
    assert _actionable("BRCA1", "2", "Likely_pathogenic")


@pytest.mark.parametrize("gene", ["ATP7B", "MUTYH", "GAA", "RPE65"])
def test_recessive_het_not_actionable(gene: str) -> None:
    """One copy in an AR gene: a carrier finding (Tier 3) or a non-finding, never reported."""
    assert not _actionable(gene, "1")


@pytest.mark.parametrize("gene", ["ATP7B", "MUTYH", "GAA", "RPE65"])
def test_recessive_hom_actionable(gene: str) -> None:
    assert _actionable(gene, "2")


def test_ttn_truncating_only() -> None:
    """Table 1 footnote f: only TTN truncating variants."""
    assert not _actionable("TTN", "1", consequence=MISSENSE)
    assert not _actionable(
        "TTN", "1"
    )  # no molecular consequence: not known to truncate
    assert _actionable("TTN", "1", consequence=NONSENSE)
    assert _actionable(
        "TTN", "1", consequence=f"{MISSENSE},SO:0001589|frameshift_variant"
    )


def test_hfe_c282y_homozygotes_only() -> None:
    assert not _actionable("HFE", "1", vid=HFE_C282Y)
    assert _actionable("HFE", "2", vid=HFE_C282Y)
    assert not _actionable("HFE", "2", vid="GRCh38:chr6:26093141:C:G")


def test_non_finding_or_watch_not_actionable() -> None:
    assert not _actionable("BRCA1", "1", "Benign")
    assert not _actionable("BRCA1", "1", "Uncertain_significance")  # a watch item
    assert not _actionable("BRCA1", "0")
    assert not _actionable("BRCA1", "unknown")  # unassessed
    assert not _actionable("BRCA1", "1", "Conflicting_classifications_of_pathogenicity")


def test_x_linked() -> None:
    """GLA and OTC: hemizygous, heterozygous and homozygous; ABCD1: hemizygous, homozygous
    or two heterozygous variants, never one (footnote k)."""
    assert _actionable("GLA", "1", male=False)
    assert _actionable("OTC", "1", male=True, haploid=True)
    assert not _actionable("ABCD1", "1", male=False)
    assert _actionable("ABCD1", "2", male=False)
    assert _actionable("ABCD1", "1", male=True, haploid=True)


def test_every_sf_gene_has_a_report_rule() -> None:
    sf = LISTS["acmg_sf"]
    assert set(sf.report) == set(sf.genes)
    assert {g for g, r in sf.report.items() if r == "biallelic"} == {
        "MUTYH", "CASQ2", "TRDN", "BTD", "CYP27A1", "GAA", "RPE65", "ATP7B", "ABCD1"
    }  # fmt: skip
    assert sf.report["TTN"] == "truncating"
    assert sf.report["HFE"] == f"homozygous:{HFE_C282Y}"


def test_x_linked_biallelic_genes_are_carrier_genes() -> None:
    """Placement makes a woman's single XL copy a carrier finding only in a Tier 3 gene; an
    XL biallelic SF gene outside Tier 3 would make one copy a risk finding, and actionable."""
    sf, tier3 = LISTS["acmg_sf"], LISTS["carrier_tier3"]
    for gene, rule in sf.report.items():
        if rule == "biallelic" and sf.genes[gene] == "XL":
            assert gene in tier3.genes, gene


def test_unknown_report_rule_refused() -> None:
    head = "# source: https://x\n# version: 1\ngene\tinheritance\treport\n"
    with pytest.raises(ValueError, match="report rule"):
        parse_gene_list("x", head + "APC\tAD\tsome\n")
    with pytest.raises(ValueError, match="report rule"):
        parse_gene_list("x", head + "HFE\tAR\thomozygous:chr6\n")
    assert parse_gene_list("x", head + "APC\tAD\tall\n").report == {"APC": "all"}
