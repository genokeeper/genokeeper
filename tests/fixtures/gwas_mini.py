"""A miniature GWAS Catalog associations download: the zip the bootstrap fetches, its one
TSV member under the name and with the 38 columns of the release of 2026-10-05 (read from
the real file then), each row an association the test describes."""

import io
import zipfile
from dataclasses import dataclass
from pathlib import Path

MEMBER = "gwas-catalog-download-associations-alt-full.tsv"
# The member's date in the zip, which names the release (YYYYMMDD).
BUILT = (2026, 10, 3, 7, 15, 0)
COLUMNS = (
    "DATE ADDED TO CATALOG", "PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL", "LINK", "STUDY",
    "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "REPLICATION SAMPLE SIZE", "REGION", "CHR_ID",
    "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID",
    "SNP_GENE_IDS", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE",
    "STRONGEST SNP-RISK ALLELE", "SNPS", "MERGED", "SNP_ID_CURRENT", "CONTEXT",
    "INTERGENIC", "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG", "P-VALUE (TEXT)",
    "OR or BETA", "95% CI (TEXT)", "PLATFORM [SNPS PASSING QC]", "CNV", "MAPPED_TRAIT",
    "MAPPED_TRAIT_URI", "STUDY ACCESSION", "GENOTYPING TECHNOLOGY",
)  # fmt: skip
EFO = "http://www.ebi.ac.uk/efo/"
OBO = "http://purl.obolibrary.org/obo/"


@dataclass(frozen=True)
class Assoc:
    """One row: its SNP, GRCh38 position, risk allele, p-value, effect, traits (URIs
    joined by `, `), study, PubMed ID and mapped gene; `snps` and `risk` as the catalog
    writes them."""

    pos: int
    risk: str
    p: str = "1E-10"
    effect: str = "1.14"
    ci: str = "[1.08-1.2]"
    trait: str = "type 2 diabetes mellitus"
    uris: str = f"{OBO}MONDO_0005148"
    study: str = "GCST000001"
    pubmed: str = "17463246"
    chrom: str = "1"
    snps: str | None = None
    gene: str = "PPARG"

    def row(self) -> dict[str, str]:
        snps = self.snps if self.snps is not None else f"rs{self.pos}"
        return dict.fromkeys(COLUMNS, "") | {
            "PUBMEDID": self.pubmed,
            "CHR_ID": self.chrom,
            "CHR_POS": str(self.pos),
            "STRONGEST SNP-RISK ALLELE": f"{snps}-{self.risk}",
            "SNPS": snps,
            "P-VALUE": self.p,
            "OR or BETA": self.effect,
            "95% CI (TEXT)": self.ci,
            "MAPPED_TRAIT": self.trait,
            "MAPPED_TRAIT_URI": self.uris,
            "STUDY ACCESSION": self.study,
            "MAPPED_GENE": self.gene,
        }


def write_gwas(path: Path, rows: list[Assoc]) -> Path:
    """The download at `path`."""
    text = "\t".join(COLUMNS) + "\n"
    for assoc in rows:
        text += "\t".join(assoc.row()[c] for c in COLUMNS) + "\n"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(zipfile.ZipInfo(MEMBER, date_time=BUILT), text)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(buffer.getvalue())
    return path
