"""Gene inheritance from the Clinical Genomic Database (CGD), for genes on neither list.

The bootstrap fetches CGD's own download whole (`research.nhgri.nih.gov/CGD/download/txt/
CGD.txt.gz`) into `<lab>/cgd/<its Last-Modified date>/CGD.txt.gz`. The gene lists (ACMG SF,
carrier Tier 3) still decide a gene's inheritance where they list it (`evidence.clinvar`);
CGD fills the others.

Its INHERITANCE column holds AD, AR or XL for most genes, and otherwise several modes
(`AD/AR`, `AR/Digenic`), blood groups, maternal or Y-linked inheritance. Only AD, AR and XL
are kept, and a gene whose rows disagree is left out: a claim in such a gene grades as
"inheritance unknown", never under a mode the table does not state alone.
"""

import csv
import gzip
import re
from pathlib import Path

from genokeeper.genelists import INHERITANCE

FILE_NAME = "CGD.txt.gz"
_DATE = re.compile(r"\d{8}")
_GENE, _MODE = "#GENE", "INHERITANCE"


def newest_cgd(lab: Path) -> Path | None:
    """The CGD table of the newest `<lab>/cgd/<YYYYMMDD>/`, or None."""
    dated = sorted(
        (d for d in (lab / "cgd").glob("*") if _DATE.fullmatch(d.name)), reverse=True
    )
    return next((d / FILE_NAME for d in dated if (d / FILE_NAME).is_file()), None)


def load_cgd(path: Path) -> dict[str, str]:
    """Gene -> AD, AR or XL from the CGD table at `path` (gzipped TSV); raises
    `ValueError` for a file whose header is not CGD's."""
    modes: dict[str, set[str]] = {}
    try:
        with gzip.open(path, "rt", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, [])
            if _GENE not in header or _MODE not in header:
                raise ValueError(f"{path}: not the CGD table (header {header[:5]})")
            gene_at, mode_at = header.index(_GENE), header.index(_MODE)
            for row in reader:
                if len(row) <= max(gene_at, mode_at) or not row[gene_at].strip():
                    continue
                gene = row[gene_at].strip()
                modes.setdefault(gene, set()).add(row[mode_at].strip())
    except (OSError, EOFError, UnicodeDecodeError, csv.Error) as e:
        raise ValueError(f"{path}: not a readable CGD table ({e})") from e
    return {
        gene: next(iter(values))
        for gene, values in modes.items()
        if len(values) == 1 and values <= INHERITANCE
    }
