"""Inherited chromosomally integrated HHV-6 (slice 3f, "HHV-6 band").

An integrated HHV-6 genome sits in a chromosome's telomere of every cell, so its reads come
at about the human depth of one copy; reads of HHV-6 are unmapped on GRCh38. The run takes
the CRAM's pairs with either read unmapped (`extract_unmapped`, one samtools pipeline behind
the core gate), aligns them with bwa mem to HHV-6A (`NC_001664.4`) and HHV-6B
(`NC_000898.1`) concatenated (`typing_files.TypingFiles.hhv6_fasta`; `align_viral`), and
reads `samtools coverage`'s mean depth over each genome's unique region, between its
terminal direct repeats (`typing_files.HHV6_UNIQUE`; `viral_coverage`, `read_coverage`).

`classify` divides the higher of the two depths by the haploid human depth (the run's mean
depth / 2): one integrated copy per cell is a ratio of 1. A ratio in the band
(`Hhv6Thresholds.min_ratio` to `max_ratio`) is a finding, one at or below the calibration's
background ceiling a non-finding, and anything else (between the two, above the band, or no
human depth to divide by) unassessed.
"""

import math
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from genokeeper.config import Hhv6Thresholds
from genokeeper.store.db import Kind
from genokeeper.tools import gated, run_tool
from genokeeper.typing_files import HHV6_UNIQUE

# The core gate's label for the extraction and the alignment.
LABEL = "dna:hhv6"
FASTQ = ("unmapped_1.fq", "unmapped_2.fq")
SPECIES = {"NC_001664.4": "HHV-6A", "NC_000898.1": "HHV-6B"}
NO_DEPTH = "not assessed: no human depth to compare the viral depth with"


@dataclass(frozen=True)
class Hhv6Call:
    """The species of the deeper genome (None when neither has a read), the ratio of its
    depth to the haploid human depth (NaN without a human depth), and what it is."""

    species: str | None
    ratio: float
    kind: Kind


def extract_unmapped(cram: Path, fasta: Path, out_dir: Path) -> tuple[Path, Path]:
    """The pairs of `cram` (indexed, written against `fasta`) with either read unmapped,
    without secondary, supplementary or duplicate records, as `out_dir/unmapped_1.fq` and
    `unmapped_2.fq`. One samtools pipeline behind one gate; the files are renamed into place
    only once it has finished."""
    out_dir.mkdir(parents=True, exist_ok=True)
    finals = [out_dir / name for name in FASTQ]
    parts = [path.with_name(path.name + ".part") for path in finals]
    script = (
        "set -euo pipefail; "
        'samtools view -u -h -F 0xD00 -e "flag.unmap || flag.munmap" -T "$1" "$2" '
        '| samtools collate -u -O - "$3/collate" '
        '| samtools fastq -n -0 /dev/null -s /dev/null -1 "$4" -2 "$5" -'
    )
    work = out_dir / ".collate"
    work.mkdir(exist_ok=True)
    try:
        run_tool(
            gated(
                LABEL,
                [
                    "bash", "-c", script, "extract-unmapped",
                    str(fasta), str(cram), str(work), str(parts[0]), str(parts[1]),
                ],
            )
        )  # fmt: skip
        for part, final in zip(parts, finals, strict=True):
            part.replace(final)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        for part in parts:
            part.unlink(missing_ok=True)
    return finals[0], finals[1]


def align_viral(fq1: Path, fq2: Path, fasta: Path, out: Path, *, threads: int) -> Path:
    """`fq1`/`fq2` aligned natively by bwa mem to `fasta` (its bwa index beside it), sorted
    and indexed as the BAM `out`, behind the gate as `LABEL`; renamed into place once the
    pipeline has finished."""
    out.parent.mkdir(parents=True, exist_ok=True)
    part = out.with_name(f".{out.name}.part.bam")
    script = (
        "set -euo pipefail; "
        f'bwa mem -t {threads} "$1" "$2" "$3" | samtools sort -o "$4" - '
        '&& samtools index "$4"'
    )
    try:
        run_tool(
            gated(
                LABEL,
                [
                    "bash", "-c", script, "align-hhv6",
                    str(fasta), str(fq1), str(fq2), str(part),
                ],
            )
        )  # fmt: skip
        Path(f"{part}.bai").replace(f"{out}.bai")
        part.replace(out)
    finally:
        part.unlink(missing_ok=True)
        Path(f"{part}.bai").unlink(missing_ok=True)
    return out


def viral_coverage(bam: Path, out: Path) -> Path:
    """`samtools coverage -r <contig>:<unique region>` of `bam` for each genome
    (`HHV6_UNIQUE`), concatenated under one header into `out`."""
    lines: list[str] = []
    for contig, (start, end) in HHV6_UNIQUE.items():
        text = run_tool(
            ["samtools", "coverage", "-r", f"{contig}:{start}-{end}", str(bam)]
        )
        rows = text.splitlines()
        lines += rows if not lines else rows[1:]
    out.write_text("\n".join(lines) + "\n")
    return out


def read_coverage(path: Path) -> dict[str, float]:
    """`samtools coverage`'s mean depth per contig."""
    lines = path.read_text().splitlines()
    header = lines[0].lstrip("#").split("\t")
    name, depth = header.index("rname"), header.index("meandepth")
    return {
        fields[name]: float(fields[depth])
        for fields in (line.split("\t") for line in lines[1:] if line.strip())
    }


def classify(
    depths: Mapping[str, float],
    haploid_depth: float,
    *,
    grading: Hhv6Thresholds,
    background: float,
) -> Hhv6Call:
    """The call from each genome's depth (`read_coverage`), the haploid human depth and the
    calibration's background ceiling `background`."""
    contig, depth = max(depths.items(), key=lambda item: item[1])
    species = SPECIES.get(contig, contig) if depth > 0 else None
    if math.isnan(haploid_depth) or haploid_depth <= 0:
        return Hhv6Call(species=species, ratio=math.nan, kind=Kind.UNASSESSED)
    ratio = depth / haploid_depth
    if ratio <= background:
        kind = Kind.NON_FINDING
    elif grading.min_ratio <= ratio <= grading.max_ratio:
        kind = Kind.FINDING
    else:
        kind = Kind.UNASSESSED
    return Hhv6Call(species=species, ratio=ratio, kind=kind)
