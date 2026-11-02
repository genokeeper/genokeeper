"""A miniature ClinPGx annotations download, built in a test's tmp_path from the rows a
test gives, in either shape the real download has: the live `summaryAnnotations.zip`
(`summary_annotations.tsv` and `summary_ann_alleles.tsv`, keyed by `Summary Annotation
ID`, read 2026-10-27, its CREATED member of 2026-10-05) or the frozen
`clinicalAnnotations.zip` (`clinical_annotations.tsv` and `clinical_ann_alleles.tsv`,
keyed by `Clinical Annotation ID`). Every other column has the same name in both. Each
holds a `CREATED_<date>.txt` member and a licence. No ClinPGx data is committed: every row
is a test's own.

Each genotype row's annotation text speaks of a dose, as many real ones do, so a test can
show no reader copies it.
"""

import io
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

Layout = Literal["summary", "clinical"]
# Each layout: its download's name, its two tables and the name of their ID column.
LAYOUTS: dict[Layout, tuple[str, str, str, str]] = {
    "summary": (
        "summaryAnnotations.zip",
        "summary_annotations.tsv",
        "summary_ann_alleles.tsv",
        "Summary Annotation ID",
    ),
    "clinical": (
        "clinicalAnnotations.zip",
        "clinical_annotations.tsv",
        "clinical_ann_alleles.tsv",
        "Clinical Annotation ID",
    ),
}
ANNOTATION_COLUMNS = (
    "Clinical Annotation ID",
    "Variant/Haplotypes",
    "Gene",
    "Level of Evidence",
    "Level Override",
    "Level Modifiers",
    "Score",
    "Phenotype Category",
    "PMID Count",
    "Evidence Count",
    "Drug(s)",
    "Phenotype(s)",
    "Latest History Date (YYYY-MM-DD)",
    "URL",
    "Specialty Population",
)
ALLELE_COLUMNS = (
    "Clinical Annotation ID",
    "Genotype/Allele",
    "Annotation Text",
    "Allele Function",
)
CREATED = "2025-07-05"
DOSING_TEXT = "Patients with this genotype may require a decreased dose of {drug}."


@dataclass(frozen=True)
class Ann:
    """One clinical annotation: its id, variant (an rsID or haplotypes), gene, level,
    category, drugs and phenotypes (`;`-joined as the download writes them), and its
    genotype rows' labels."""

    id: str
    variant: str
    gene: str
    level: str
    genotypes: tuple[str, ...] = ("CC", "CT", "TT")
    category: str = "Efficacy"
    drugs: str = "somedrug"
    phenotypes: str = "Some disease"


def zip_bytes(
    annotations: list[Ann], *, created: str = CREATED, layout: Layout = "summary"
) -> bytes:
    """The download's bytes holding `annotations`, in the live (`summary`) or the frozen
    (`clinical`) layout."""
    _, annotations_tsv, alleles_tsv, id_column = LAYOUTS[layout]
    rows = ["\t".join((id_column, *ANNOTATION_COLUMNS[1:]))]
    alleles = ["\t".join((id_column, *ALLELE_COLUMNS[1:]))]
    for a in annotations:
        rows.append(
            "\t".join(
                [a.id, a.variant, a.gene, a.level, "", "", "1.0", a.category, "1", "1"]
                + [a.drugs, a.phenotypes, "2024-01-01"]
                + [f"https://www.clinpgx.org/clinicalAnnotation/{a.id}", ""]
            )
        )
        for label in a.genotypes:
            text = DOSING_TEXT.format(drug=a.drugs.split(";")[0])
            alleles.append(f"{a.id}\t{label}\t{text}\t")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        z.writestr("LICENSE.txt", "Creative Commons Attribution-ShareAlike 4.0\n")
        z.writestr(f"CREATED_{created}.txt", f"Created on {created}.\n")
        z.writestr(annotations_tsv, "\n".join(rows) + "\n")
        z.writestr(alleles_tsv, "\n".join(alleles) + "\n")
    return buffer.getvalue()


def write_annotations(
    path: Path,
    annotations: list[Ann],
    *,
    created: str = CREATED,
    layout: Layout = "summary",
) -> Path:
    """The download at `path`, in the live (`summary`) or the frozen (`clinical`) layout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(zip_bytes(annotations, created=created, layout=layout))
    return path
