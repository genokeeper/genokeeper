"""PharmGKB's clinical annotations, from the ClinPGx download (track design, slice 1b,
"PharmGKB clinical annotations").

The download (`DOWNLOADS`, fetched whole by the bootstrap's step_pharmgkb into
`<lab>/pharmgkb/<its CREATED date>/`) holds one row per annotation (its variant, gene,
level of evidence, category, drugs and phenotypes) and one row per genotype of it. ClinPGx
serves it in two layouts of the same columns, read alike: the live `summaryAnnotations.zip`
(`summary_annotations.tsv`, `summary_ann_alleles.tsv`, keyed by `Summary Annotation ID`)
and the frozen `clinicalAnnotations.zip`, which an older bootstrap fetched
(`clinical_annotations.tsv`, `clinical_ann_alleles.tsv`, keyed by `Clinical Annotation
ID`); an id has the same grades in both. An annotation is kept
(`read_clinical_annotations`):

- at levels 1B, 2A and 2B only: 1A restates a guideline, which PharmCAT reports, and 3 and 4
  fall below the bar. 1B is moderate validity, 2A and 2B limited (`validity`). Level 3 (one
  study or an unconfirmed association) is kept apart (`ClinicalAnnotations.research`), read
  and filtered alike, for the research layer (slice 3f), which grades it limited; the
  `level_3` count still says where those rows went from the 1B–2B statements' tally;
- in a gene PharmCAT calls (`pharmcat.PHARMCAT_GENES`), only when PharmCAT's report for the
  run gives that gene no drug at all (`left_by_pharmcat`: F2, F5 and IFNL3 on 3.4.0);
  otherwise PharmCAT's gene result answers. No drug name is matched.

Every other row is counted by why (`PHARMGKB_DROP_REASONS`); `pharmcat_gene` once the
report is read.

The download names no position, so an annotation is located (`locate`) at its rsID's allele
in the run's ClinVar GRCh38 release (`clinvar_alleles`), normalised as every claim's allele
is (`panel.build.normalise_records`) (builder's call). Single-base variants only: the
annotation's genotypes (`AG`) must hold both the allele's REF and ALT, which, PharmGKB
writing them on GRCh38's plus strand, also refuses an allele read on the other strand; and
exactly one of the rsID's alleles must match. An annotation of haplotypes (`CYP2A6*1,
CYP2A6*2`, an HLA allele) has no caller here: it is counted, not located.

No text of an annotation is read, not its genotypes' `Annotation Text` (which often speaks
of a dose) nor its level override's note. This control, against dosing text in the store:
its anchor comes from the columns this module names, which an annotation's text cannot add
to, because every value is read from a named column and no free-text column is named; and a
category, drug or phenotype that speaks of a dose (`pharmcat.DOSING`, the `Dosage` category
among them) is left out of what is kept.
"""

import csv
import io
import re
import zipfile
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.genotype.gvcf import Allele
from genokeeper.panel.build import is_primary, normalise_records
from genokeeper.pgx.pharmcat import DOSING, PHARMCAT_GENES
from genokeeper.store.db import Validity
from genokeeper.variants.key import ucsc_contig

# The live download first, then the frozen one an older bootstrap fetched.
DOWNLOADS = ("summaryAnnotations.zip", "clinicalAnnotations.zip")
# The member naming the day ClinPGx built the download: `CREATED_2025-07-05.txt`.
_CREATED = re.compile(r"CREATED_(\d{4})-(\d{2})-(\d{2})\.txt")
# The columns read; no other is.
_VARIANT, _GENE, _LEVEL = "Variant/Haplotypes", "Gene", "Level of Evidence"
_CATEGORY, _DRUGS, _PHENOTYPES = "Phenotype Category", "Drug(s)", "Phenotype(s)"
_GENOTYPE = "Genotype/Allele"


@dataclass(frozen=True)
class _Layout:
    """One download layout: its two tables and the column keying them."""

    annotations: str
    alleles: str
    id: str


# The two layouts, the live one first; the members present choose.
_LAYOUTS = (
    _Layout(
        "summary_annotations.tsv", "summary_ann_alleles.tsv", "Summary Annotation ID"
    ),
    _Layout(
        "clinical_annotations.tsv", "clinical_ann_alleles.tsv", "Clinical Annotation ID"
    ),
)
LEVELS = {"1B": Validity.MODERATE, "2A": Validity.LIMITED, "2B": Validity.LIMITED}
# The level the research layer reads (slice 3f); limited, as every research association is.
RESEARCH_LEVEL = "3"
# Where every annotation that is no statement went.
PHARMGKB_DROP_REASONS = (
    "level_1a",
    "level_3",
    "level_4",
    "other_level",
    "pharmcat_gene",
    "not_an_rsid",
    "rsid_not_in_clinvar",
    "not_a_snv",
    "alleles_differ",
    "multiallelic",
)
_LEVEL_DROPS = {"1A": "level_1a", "3": "level_3", "4": "level_4"}
_RSID = re.compile(r"rs\d+")
_SNV_GENOTYPE = re.compile(r"[ACGT]{2}")
_BASES = re.compile(r"[ACGT]+")


@dataclass(frozen=True)
class Annotation:
    """One clinical annotation: its id, variant (an rsID, or haplotypes), first gene, its
    genes PharmCAT calls, level of evidence, its categories, drugs and phenotypes (none
    speaking of a dose), and the labels of its genotypes (`AG`)."""

    annotation_id: str
    variant: str
    gene: str
    pharmcat_genes: frozenset[str]
    level: str
    categories: tuple[str, ...]
    drugs: tuple[str, ...]
    phenotypes: tuple[str, ...]
    genotypes: tuple[str, ...]


@dataclass(frozen=True)
class ClinicalAnnotations:
    """The annotations kept, how many rows the download holds, where the others went, the
    date ClinPGx built it (YYYYMMDD), the statements' source version, and the level-3
    annotations kept for the research layer (also counted `level_3` in `dropped`)."""

    annotations: list[Annotation]
    input_rows: int
    dropped: dict[str, int]
    created: str
    research: list[Annotation]


def _listed(text: str) -> tuple[str, ...]:
    """A `;`-joined field's items, each that speaks of a dose left out."""
    items = (part.strip() for part in text.split(";"))
    return tuple(item for item in items if item and not DOSING.search(item))


def _rows(
    z: zipfile.ZipFile, name: str, columns: Iterable[str]
) -> list[dict[str, str]]:
    text = z.read(name).decode("utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE)
    missing = [c for c in columns if c not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{name} has no column {', '.join(missing)}")
    return list(reader)


def read_clinical_annotations(path: Path) -> ClinicalAnnotations:
    """The download at `path`, in either layout, read as the module docstring says. A file
    that is not the download (not a zip, or without the two tables of one layout and a
    CREATED member) is refused (`ValueError`)."""
    try:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            created = [m for n in names if (m := _CREATED.fullmatch(n))]
            layout = next(
                (x for x in _LAYOUTS if {x.annotations, x.alleles} <= names), None
            )
            if len(created) != 1 or layout is None:
                tables = " or ".join(f"{x.annotations} + {x.alleles}" for x in _LAYOUTS)
                raise ValueError(f"no tables ({tables}) or CREATED date")
            id_column = layout.id
            rows = _rows(
                z,
                layout.annotations,
                (id_column, _VARIANT, _GENE, _LEVEL, _CATEGORY, _DRUGS, _PHENOTYPES),
            )
            genotypes: dict[str, list[str]] = {}
            for row in _rows(z, layout.alleles, (id_column, _GENOTYPE)):
                genotypes.setdefault(row[id_column], []).append(row[_GENOTYPE].strip())
    except (zipfile.BadZipFile, KeyError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(
            f"{path}: not the ClinPGx clinical annotations download ({e})"
        ) from e
    dropped: dict[str, int] = dict.fromkeys(PHARMGKB_DROP_REASONS, 0)
    kept: list[Annotation] = []
    research: list[Annotation] = []
    for row in rows:
        level = row[_LEVEL].strip()
        genes = [g.strip() for g in row[_GENE].split(";") if g.strip()]
        if level not in LEVELS:
            dropped[_LEVEL_DROPS.get(level, "other_level")] += 1
        if level in LEVELS or level == RESEARCH_LEVEL:
            (kept if level in LEVELS else research).append(
                Annotation(
                    annotation_id=row[id_column].strip(),
                    variant=row[_VARIANT].strip(),
                    gene=genes[0] if genes else "",
                    pharmcat_genes=PHARMCAT_GENES.intersection(genes),
                    level=level,
                    categories=_listed(row[_CATEGORY]),
                    drugs=_listed(row[_DRUGS]),
                    phenotypes=_listed(row[_PHENOTYPES]),
                    genotypes=tuple(genotypes.get(row[id_column], ())),
                )
            )
    year, month, day = created[0].groups()
    return ClinicalAnnotations(
        kept, len(rows), dropped, f"{year}{month}{day}", research
    )


def left_by_pharmcat(
    annotations: Sequence[Annotation], drugged: Collection[str]
) -> tuple[list[Annotation], int]:
    """The `annotations` PharmCAT's report leaves to PharmGKB, those in none of the genes it
    gives a drug in this run (`drugged`), and how many it answers (`pharmcat_gene`)."""
    left = [a for a in annotations if not a.pharmcat_genes & set(drugged)]
    return left, len(annotations) - len(left)


def validity(level: str) -> Validity:
    """An annotation's validity: moderate for 1B, limited for 2A, 2B and 3."""
    return Validity.LIMITED if level == RESEARCH_LEVEL else LEVELS[level]


def _rs_numbers(value: object) -> set[str]:
    if value is None:
        return set()
    parts = value if isinstance(value, tuple) else (value,)
    return {n for p in parts for n in str(p).split(",") if n}


def clinvar_alleles(
    clinvar38: Path, rsids: Collection[str], *, fasta: Path
) -> dict[str, list[Allele]]:
    """Each of `rsids`' alleles in the ClinVar GRCh38 release `clinvar38` (its INFO `RS`),
    on the primary contigs of `fasta`, normalised against it (an allele whose REF the FASTA
    does not have is left out), in the release's order; an rsID ClinVar does not name has
    none."""
    wanted = {r.removeprefix("rs") for r in rsids if _RSID.fullmatch(r)}
    with pysam.FastaFile(str(fasta)) as f:
        contigs = set(f.references)
    raw: list[tuple[str, Allele]] = []
    with pysam.VariantFile(str(clinvar38)) as vcf:
        for rec in vcf:
            hits = _rs_numbers(rec.info.get("RS")) & wanted
            chrom = ucsc_contig(rec.chrom)
            if not hits or rec.ref is None or not is_primary(rec.chrom):
                continue
            if chrom not in contigs:
                continue
            ref = rec.ref.upper()
            for alt in rec.alts or ():
                if _BASES.fullmatch(ref) and _BASES.fullmatch(alt.upper()):
                    raw += [
                        (f"rs{n}", (chrom, rec.pos, ref, alt.upper())) for n in hits
                    ]
    normalised = normalise_records([allele for _, allele in raw], fasta=fasta)
    found: dict[str, list[Allele]] = {}
    for i, (rsid, _) in enumerate(raw):
        allele = normalised.get(i)
        if allele is not None and allele not in found.setdefault(rsid, []):
            found[rsid].append(allele)
    return found


def _snv(allele: Allele) -> bool:
    return len(allele[2]) == 1 and len(allele[3]) == 1


def locate(
    annotations: Sequence[Annotation], alleles: Mapping[str, Sequence[Allele]]
) -> tuple[list[tuple[Annotation, Allele]], dict[str, int]]:
    """Each annotation at its rsID's one matching allele of `alleles`
    (`clinvar_alleles`), and where the others went (see the module docstring)."""
    located: list[tuple[Annotation, Allele]] = []
    dropped: dict[str, int] = dict.fromkeys(PHARMGKB_DROP_REASONS, 0)
    for a in annotations:
        found = alleles.get(a.variant, ())
        snvs = [allele for allele in found if _snv(allele)]
        labels = [g for g in a.genotypes if g]
        bases = {b for g in labels for b in g}
        if not _RSID.fullmatch(a.variant):
            reason = "not_an_rsid"
        elif not found:
            reason = "rsid_not_in_clinvar"
        elif not snvs or not labels or not all(map(_SNV_GENOTYPE.fullmatch, labels)):
            reason = "not_a_snv"
        else:
            matching = [s for s in snvs if {s[2], s[3]} <= bases]
            if len(matching) == 1:
                located.append((a, matching[0]))
                continue
            reason = "multiallelic" if matching else "alleles_differ"
        dropped[reason] += 1
    return located, dropped
