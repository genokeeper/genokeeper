"""Versioned gene lists (ACMG SF, ACMG carrier Tiers 3 and 4, the structural screen's own
genes), read from the package's TSV files.

Each file starts with `# source: <URL>` and `# version: <x>`; other `# key: value` lines are
notes for a reader. Then comes a `gene<TAB>inheritance` header and one row per gene, or a
`gene<TAB>inheritance<TAB>report` header and each gene's report rule (`REPORT_RULES`, or
`homozygous:<vid>`): ACMG SF's "variants to report", which `evidence.clinvar.sf_actionable`
reads.

Carrier status reads Tier 3 and Tier 4 together (`carrier_genes`): ACMG's Tier 4 is every
gene ClinGen curates AR or XL at Definitive or Strong, cut from its gene-validity export.
`screen_genes.tsv` names the genes the structural screens read beyond the two ACMG lists,
each with how a loss there is placed and why it is read; its four columns are its own
loader's (`load_screen_genes`), not `parse_gene_list`'s. This module is a leaf: it imports
no source and no evidence module, which import it.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib.resources import files

INHERITANCE = frozenset({"AD", "AR", "XL"})
REPORT_ALL = "all"
REPORT_BIALLELIC = "biallelic"
REPORT_TRUNCATING = "truncating"
REPORT_HOMOZYGOUS = "homozygous"  # homozygous:<vid>, that one allele's homozygotes only
REPORT_RULES = frozenset({REPORT_ALL, REPORT_BIALLELIC, REPORT_TRUNCATING})

# The lists whose genes the VUS watch rule, the panel's ClinVar claims, SpliceAI and the
# coverage table read; Tier 4 enters only the carrier gate (`carrier_genes`).
LISTED = ("acmg_sf", "carrier_tier3")
_FILES = {
    "acmg_sf": "acmg_sf_v3.3.tsv",
    "carrier_tier3": "acmg_carrier_tier3_2021.tsv",
    "carrier_tier4": "acmg_carrier_tier4.tsv",
}
SCREEN_GENES = "screen_genes.tsv"
SCREEN_COLUMNS = ("gene", "inheritance", "placement", "reason")
# How a loss in a screen gene is placed: through the gene gate, as any listed gene's
# (`lof`), or as a drug-response finding instead (`drug_response`, A4).
PLACEMENTS = frozenset({"lof", "drug_response"})
DPYD_NOTE = "may reduce DPD activity: CPIC advises testing before fluoropyrimidines"


@dataclass(frozen=True)
class GeneList:
    """A named gene list: gene symbol -> inheritance (AD, AR or XL), and gene symbol ->
    report rule where the file has the column."""

    name: str
    version: str
    source: str
    genes: dict[str, str]
    report: dict[str, str] = field(default_factory=dict)


def _report_rule(rule: str) -> bool:
    kind, sep, vid = rule.partition(":")
    if not sep:
        return rule in REPORT_RULES
    return kind == REPORT_HOMOZYGOUS and len(vid.split(":")) == 5


def parse_gene_list(name: str, text: str) -> GeneList:
    """Parse one gene-list TSV; refuse an unknown inheritance or report rule, a duplicate or
    missing metadata."""
    meta: dict[str, str] = {}
    genes: dict[str, str] = {}
    report: dict[str, str] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        if line.startswith("#"):
            key, sep, value = line.removeprefix("#").partition(":")
            if sep:
                meta.setdefault(key.strip(), value.strip())
            continue
        gene, inheritance, *rule = line.split("\t")
        if gene == "gene":
            continue
        if inheritance not in INHERITANCE:
            raise ValueError(f"{name}:{number}: unknown inheritance {inheritance!r}")
        if gene in genes:
            raise ValueError(f"{name}:{number}: duplicate gene {gene}")
        genes[gene] = inheritance
        if rule:
            if len(rule) > 1 or not _report_rule(rule[0]):
                raise ValueError(f"{name}:{number}: unknown report rule {rule!r}")
            report[gene] = rule[0]
    for key in ("source", "version"):
        if key not in meta:
            raise ValueError(f"{name}: no '# {key}:' line")
    return GeneList(
        name=name,
        version=meta["version"],
        source=meta["source"],
        genes=genes,
        report=report,
    )


def load_gene_lists() -> dict[str, GeneList]:
    """The packaged gene lists, keyed `acmg_sf`, `carrier_tier3` and `carrier_tier4`."""
    folder = files("genokeeper").joinpath("resources", "genelists")
    return {
        name: parse_gene_list(name, folder.joinpath(filename).read_text())
        for name, filename in _FILES.items()
    }


def carrier_genes(lists: Mapping[str, GeneList]) -> dict[str, str]:
    """The carrier genes, Tier 3 and Tier 4, with their inheritance; Tier 3's wins on a
    gene in both."""
    return lists["carrier_tier4"].genes | lists["carrier_tier3"].genes


def gene_mode(
    gene: str, lists: Mapping[str, GeneList], inheritance: str | None
) -> str | None:
    """The gene's inheritance: ACMG SF's, else Tier 3's, else the caller's `inheritance`
    (CGD's, which states one only where its rows agree), else Tier 4's. Tier 4 comes last
    because CGD disagrees with it on 11 of its genes (CGD 2025-02-21): ADK, FASLG, HMX1,
    KNG1, LDB3, MSH3, MYPN, PSAP, RARS1 and SYNE1 are AD in CGD, CSF2RA XL; one copy there
    keeps CGD's handling, not carrier status. Its anchor comes from ACMG SF, Tier 3 and
    CGD, which Tier 4's row cannot move because they are read first."""
    return (
        lists["acmg_sf"].genes.get(gene)
        or lists["carrier_tier3"].genes.get(gene)
        or inheritance
        or lists["carrier_tier4"].genes.get(gene)
    )


@dataclass(frozen=True)
class ScreenGene:
    """A gene the structural screens read beyond the ACMG lists: its inheritance, how a
    loss there is placed (`PLACEMENTS`) and why the screen reads it."""

    gene: str
    inheritance: str
    placement: str
    reason: str


def parse_screen_genes(text: str) -> dict[str, ScreenGene]:
    """Parse `screen_genes.tsv`; refuse an unknown inheritance or placement, an empty
    reason, a duplicate, another column set or missing metadata."""
    meta: dict[str, str] = {}
    genes: dict[str, ScreenGene] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        if line.startswith("#"):
            key, sep, value = line.removeprefix("#").partition(":")
            if sep:
                meta.setdefault(key.strip(), value.strip())
            continue
        fields = tuple(line.split("\t"))
        if fields == SCREEN_COLUMNS:
            continue
        if len(fields) != len(SCREEN_COLUMNS):
            raise ValueError(f"{SCREEN_GENES}:{number}: expected {SCREEN_COLUMNS}")
        gene, inheritance, placement, reason = fields
        if inheritance not in INHERITANCE:
            raise ValueError(
                f"{SCREEN_GENES}:{number}: unknown inheritance {inheritance!r}"
            )
        if placement not in PLACEMENTS:
            raise ValueError(
                f"{SCREEN_GENES}:{number}: unknown placement {placement!r}"
            )
        if not reason.strip():
            raise ValueError(f"{SCREEN_GENES}:{number}: {gene} has no reason")
        if gene in genes:
            raise ValueError(f"{SCREEN_GENES}:{number}: duplicate gene {gene}")
        genes[gene] = ScreenGene(gene, inheritance, placement, reason)
    for key in ("source", "version"):
        if key not in meta:
            raise ValueError(f"{SCREEN_GENES}: no '# {key}:' line")
    return genes


def load_screen_genes() -> dict[str, ScreenGene]:
    """The packaged screen genes (DPYD, CHEK2)."""
    folder = files("genokeeper").joinpath("resources", "genelists")
    return parse_screen_genes(folder.joinpath(SCREEN_GENES).read_text())


def screen_overrides() -> dict[str, tuple[str, str]]:
    """gene -> (the class it is placed as, the note its reason carries) for each screen
    gene whose placement is not the gene gate's (A4: a DPYD loss is a drugs finding)."""
    return {
        g.gene: ("drug_response", DPYD_NOTE)
        for g in load_screen_genes().values()
        if g.placement == "drug_response"
    }
