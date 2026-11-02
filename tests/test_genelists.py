import subprocess
import sys
from collections import Counter

import pytest

from genokeeper.genelists import (
    DPYD_NOTE,
    carrier_genes,
    load_gene_lists,
    load_screen_genes,
    parse_screen_genes,
    screen_overrides,
)


def test_acmg_sf_has_84_genes() -> None:
    sf = load_gene_lists()["acmg_sf"]
    assert len(sf.genes) == 84
    assert {"PLN", "ABCD1", "CYP27A1"} <= set(sf.genes)
    assert sf.version == "3.3"
    assert sf.source.startswith("https://")


def test_carrier_tier3_has_113_genes() -> None:
    tier3 = load_gene_lists()["carrier_tier3"]
    assert len(tier3.genes) == 113
    assert Counter(tier3.genes.values()) == {"AR": 97, "XL": 16}


def test_every_inheritance_is_known() -> None:
    for gene_list in load_gene_lists().values():
        assert set(gene_list.genes.values()) <= {"AD", "AR", "XL"}, gene_list.name


def test_screen_genes_loaded() -> None:
    """DPYD and CHEK2, the structural screen's genes off both ACMG lists (slice 3d, R4 and
    R7); DPYD's placement overrides the gene gate (A4). The module stays a leaf: it imports
    no source and no evidence module, which import it (C23)."""
    screen = load_screen_genes()
    assert set(screen) == {"DPYD", "CHEK2"}
    assert (screen["DPYD"].inheritance, screen["DPYD"].placement) == (
        "AR",
        "drug_response",
    )
    assert (screen["CHEK2"].inheritance, screen["CHEK2"].placement) == ("AD", "lof")
    assert all(g.reason for g in screen.values())
    assert screen_overrides() == {"DPYD": ("drug_response", DPYD_NOTE)}
    code = (
        "import sys, genokeeper.genelists;"
        "print([m for m in sys.modules"
        " if m.startswith(('genokeeper.sources', 'genokeeper.evidence'))])"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert out.stdout.strip() == "[]"


def test_screen_genes_refuse_an_unknown_placement() -> None:
    head = "# source: x\n# version: 1\ngene\tinheritance\tplacement\treason\n"
    with pytest.raises(ValueError, match="placement"):
        parse_screen_genes(head + "DPYD\tAR\tdrugs\twhy\n")
    with pytest.raises(ValueError, match="inheritance"):
        parse_screen_genes(head + "DPYD\tAD/AR\tlof\twhy\n")
    with pytest.raises(ValueError, match="reason"):
        parse_screen_genes(head + "DPYD\tAR\tlof\t\n")


def test_tier4_list_loaded() -> None:
    """ClinGen's AR and XL curations at Definitive or Strong with no AD or SD curation at
    Moderate or above (2026-10-30 export): 1,309 genes; `carrier_genes` is Tier 3 and
    Tier 4 together, Tier 3's inheritance winning."""
    lists = load_gene_lists()
    tier4 = lists["carrier_tier4"]
    assert len(tier4.genes) == 1309
    assert Counter(tier4.genes.values()) == {"AR": 1160, "XL": 149}
    assert not {"ATM", "PKD1", "RYR1", "F5"} & set(tier4.genes)
    assert tier4.version == "2026-10-30"
    assert tier4.source.startswith("https://search.clinicalgenome.org/")
    carrier = carrier_genes(lists)
    tier3 = lists["carrier_tier3"].genes
    assert all(carrier[g] == mode for g, mode in tier3.items())
    assert set(tier4.genes) <= set(carrier)
    assert set(carrier) == set(tier3) | set(tier4.genes)
    assert carrier["ABCA4"] == "AR" and "ABCA4" not in tier3
