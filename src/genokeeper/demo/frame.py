"""The demo's frame, in one place: every page, and the demo folder's `CLAUDE.md`, render it
from here. `FRAME_SUMMARY` is the line always shown, `FRAME` and `CALIBRATION_NOTE` the two
paragraphs under it.

`TABLE_SAMPLES` says, for each calibration table (`calibrate.tables.TABLES`), what it
measures and on which libraries, so the calibration page can say which rates are in-sample
for HG002. It is kept by hand from the table builders (`calibrate.build`, `slice3d.TRIO`,
`slice3f.LIBRARIES`, `foundation.CONTROLS`); a test fails when a table is missing or extra.
"""

from dataclasses import dataclass
from typing import Literal

FRAME_SUMMARY = (
    "This is HG002, a public reference genome, not anyone's health record. It is also one"
    " of the samples the tool was calibrated on."
)
FRAME = (
    "HG002 comes from a Personal Genome Project participant under PGP's open consent. The"
    " results illustrate the tool. They are not a statement about anyone's health, and"
    " nobody should try to identify the person. The dashboard is written as its owner would"
    ' read it, so it says "you".'
)
CALIBRATION_PAGE = "https://genokeeper.github.io/genokeeper/calibration.html"
# The note's last clause; each site page shows it its own way (`package._note`).
POINTER = f"; the calibration page ({CALIBRATION_PAGE}) says"
CALIBRATION_NOTE = (
    "Most small-variant error rates were measured on HG003, HG002's father: a different"
    " library of a related sample. The rates for chrX, chrY and the regions only GIAB's"
    " v5.0q truth reaches, and those for repeats, deletions, copy number, insertions,"
    " paralog loci and the run-time checks, were measured on HG002's own libraries, this"
    " one among them, often beside HG003 and HG004, its mother: they are shown on data they"
    " were measured on. HG002, HG003 and HG004 are one Ashkenazi Jewish family, so every"
    " rate comes from one family of one ancestry. Inside GIAB's benchmark regions a call"
    " can be checked against a truth set; a call outside them is marked as outside the"
    f" calibrated regions{POINTER} which table was measured on which samples."
)


@dataclass(frozen=True)
class Measured:
    """A calibration table: what it measures, the libraries it was measured on, and whether
    HG002 is among them."""

    measures: str
    samples: str
    for_hg002: Literal["related sample", "same sample"]


TABLE_SAMPLES: dict[str, Measured] = {
    "small_variants": Measured(
        "Small-variant genotypes, chr1-22", "HG003", "related sample"
    ),
    "hg002_cells": Measured(
        "Small-variant genotypes on chrX, chrY and v5.0q-only regions",
        "HG002 PCR-free",
        "same sample",
    ),
    "stage3": Measured(
        "SMN1 and SMN2 copies, repeat lengths, the Y haplogroup",
        "HG002 PCR-free and PCR; HG003 and HG004 for repeats",
        "same sample",
    ),
    "deletions": Measured(
        "Deletions of 50 bp or more in ACMG SF and carrier genes",
        "HG002 PCR-free and PCR",
        "same sample",
    ),
    "cnv": Measured(
        "Copy number from read depth",
        "HG002 PCR-free and PCR, HG003, HG004",
        "same sample",
    ),
    "insertions": Measured(
        "Insertions and mobile elements", "HG002 PCR-free and PCR", "same sample"
    ),
    "paralogs": Measured(
        "Paralog loci (GBA1, Lp(a) and others)",
        "HG002 PCR-free and PCR, HG003, HG004",
        "same sample",
    ),
    "foundation": Measured(
        "Run-time checks: sex dosage, allele fraction",
        "HG002 PCR-free and PCR, HG003, HG004",
        "same sample",
    ),
    "slice3f": Measured(
        "HLA, chrM variants, HHV-6, runs of homozygosity, blood groups, ancestry",
        "HG002 PCR-free, HG003, HG004",
        "same sample",
    ),
}
