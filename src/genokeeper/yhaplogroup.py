"""The Y haplogroup from Yleaf (slice 3a).

Yleaf runs offline in its image (`containers.IMAGES["yleaf"]`, built from docker/, with its
tree and position files) behind the core gate as `dna:yleaf`, on the CRAM, its index and the
FASTA (with its `.fai`) mounted read-only and one output directory mounted read-write.
`--ref-fasta` keeps it from fetching a reference (it has no network) and `--no-update-check`
from checking for new data.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

from genokeeper.containers import run_container

YLEAF_VERSION = "4.1.6"
LABEL = "dna:yleaf"
PREDICTION = "hg_prediction.hg"
FOLDER = "yleaf"


@dataclass(frozen=True)
class YHaplogroup:
    """Yleaf's haplogroup, its defining marker, its QC score (the share of the path's markers
    in a consistent state) and the markers it read."""

    name: str
    marker: str
    qc_score: float
    valid_markers: int


def run_yleaf(cram: Path, fasta: Path, out_dir: Path, *, threads: int) -> Path:
    """Yleaf on `cram` (GRCh38) into `out_dir/yleaf`, which must not exist yet: Yleaf asks
    on stdin before reusing an output folder, and the container has none. Returns its
    prediction `out_dir/yleaf/hg_prediction.hg`; the sample's markers are in
    `out_dir/yleaf/<sample>/<sample>.out`."""
    folder = out_dir / FOLDER
    if folder.exists():
        raise FileExistsError(f"{folder}: Yleaf needs a new output folder")
    run_container(
        "yleaf",
        [
            "Yleaf", "-cram", f"/cram/{cram.name}", "-cr", f"/ref/{fasta.name}",
            "--ref-fasta", f"/ref/{fasta.name}", "-rg", "hg38", "-o", f"/out/{FOLDER}",
            "-t", str(threads), "--no-update-check",
        ],
        mounts=[
            (cram.parent, "/cram", True),
            (fasta.parent, "/ref", True),
            (out_dir, "/out", False),
        ],
        label=LABEL,
    )  # fmt: skip
    return folder / PREDICTION


def read_yleaf(path: Path) -> YHaplogroup | None:
    """The one sample's haplogroup in Yleaf's prediction at `path`; None when Yleaf predicted
    nothing (no row, or `NA`); `ValueError` for more than one sample."""
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) > 1:
        raise ValueError(f"{path}: {len(rows)} samples, expected one sample")
    if not rows or rows[0]["Hg"] in ("", "NA"):
        return None
    row = rows[0]
    return YHaplogroup(
        row["Hg"], row["Hg_marker"], float(row["QC-score"]), int(row["Valid_markers"])
    )


def derived_markers(out_file: Path) -> set[str]:
    """The markers a sample's Yleaf `.out` file calls derived (`state` is `D`)."""
    with out_file.open(newline="") as f:
        rows = csv.DictReader(f, delimiter="\t")
        return {row["marker_name"] for row in rows if row["state"] == "D"}
