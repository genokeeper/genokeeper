"""The GIAB v3.6 GRCh38 BEDs stage 2 reads, and the two region BEDs derived from them.

`scripts/bootstrap.sh` (step_giab38) downloads `giab38_downloads(grading)` into
`giab38_dir(lab)`, md5-checked, then writes the derived BEDs there:

    python -m genokeeper.regions downloads      # the GIAB file names, one per line
    python -m genokeeper.regions derived        # the derived file names
    python -m genokeeper.regions derive DIR     # write the derived BEDs into DIR

GIAB v3.6 has no X-degenerate BED, so the sex regions are derived (builder's call): chrY's
non-PAR minus its X-transposed and ampliconic regions is the X-degenerate region (Skaletsky
et al. 2003); minus GIAB's gaps (+15 kb) and low-mappability regions it keeps about 8.5 Mb,
the single-copy sequence where a read at MAPQ >= 30 can only come from a Y. Slice 3c adds
chrX's mappable non-PAR region the same way (`xnonpar`): chrX's non-PAR minus its
X-transposed and ampliconic regions, gaps and low mappability, where a normal male's reads
are X's alone (Y-derived reads in the PARs, the XTR and the ampliconic blocks would balance
the allele depths), for the X copy count and the chrX allele balance (`genokeeper.sex`).
"""

import argparse
import gzip
import os
import sys
from collections.abc import Iterable
from pathlib import Path

from genokeeper.config import Grading, load_grading

type Interval = tuple[str, int, int]

CHRY_NONPAR = "GRCh38_chrY_nonPAR.bed.gz"
CHRY_XTR = "GRCh38_chrY_XTR.bed.gz"
CHRY_AMPLICONIC = "GRCh38_chrY_ampliconic.bed.gz"
AUTOSOMES = "GRCh38_AllAutosomes.bed.gz"
GAPS = "GRCh38_gaps_slop15kb.bed.gz"
LOW_MAPPABILITY = "GRCh38_lowmappabilityall.bed.gz"
CHRX_NONPAR = "GRCh38_chrX_nonPAR.bed.gz"
CHRX_XTR = "GRCh38_chrX_XTR.bed.gz"
CHRX_AMPLICONIC = "GRCh38_chrX_ampliconic.bed.gz"
_SOURCES = (
    CHRY_NONPAR,
    CHRY_XTR,
    CHRY_AMPLICONIC,
    AUTOSOMES,
    GAPS,
    LOW_MAPPABILITY,
    CHRX_NONPAR,
    CHRX_XTR,
    CHRX_AMPLICONIC,
)

# Derived: name `xdr` on chrY's mappable X-degenerate region, `xnonpar` on chrX's mappable
# non-PAR region and `autosomes` on the mappable autosomes (genokeeper.sex); and the
# autosomes outside gaps, where the QC gate measures depth.
SEX_REGIONS = "GRCh38_sex_regions.derived.bed"
QC_AUTOSOMES = "GRCh38_qc_autosomes.derived.bed"
DERIVED = (SEX_REGIONS, QC_AUTOSOMES)


def giab38_dir(lab: Path) -> Path:
    """Where the lab keeps the GIAB v3.6 GRCh38 BEDs and the derived ones."""
    return lab / "giab" / "v3.6-GRCh38"


def giab38_downloads(grading: Grading) -> list[str]:
    """Every GIAB file name stage 2 needs: the [strata] BEDs and the derivations' sources."""
    return sorted(set(grading.strata.values()) | set(_SOURCES))


def read_bed(path: Path, chrom: str | None = None) -> list[Interval]:
    """The (chrom, start, end) intervals of a BED, gzipped or not, optionally of one chrom."""
    opener = gzip.open if path.suffix == ".gz" else open
    intervals: list[Interval] = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            name, start, end = line.split("\t")[:3]
            if chrom is None or name == chrom:
                intervals.append((name, int(start), int(end)))
    return intervals


def read_named(path: Path, name: str) -> list[Interval]:
    """The intervals of a named BED (`SEX_REGIONS`, `QC_AUTOSOMES`) whose name is `name`."""
    found: list[Interval] = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 4 and fields[3] == name:
            found.append((fields[0], int(fields[1]), int(fields[2])))
    return found


def merge(intervals: Iterable[Interval]) -> list[Interval]:
    """Sorted, with overlapping and touching intervals joined."""
    out: list[Interval] = []
    for chrom, start, end in sorted(intervals):
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            out[-1] = (chrom, out[-1][1], max(out[-1][2], end))
        else:
            out.append((chrom, start, end))
    return out


def subtract(keep: Iterable[Interval], remove: Iterable[Interval]) -> list[Interval]:
    """The parts of `keep` outside `remove`, merged and sorted."""
    cuts: dict[str, list[Interval]] = {}
    for interval in merge(remove):
        cuts.setdefault(interval[0], []).append(interval)
    out: list[Interval] = []
    index: dict[str, int] = {}
    for chrom, start, end in merge(keep):
        chrom_cuts = cuts.get(chrom, [])
        i = index.get(chrom, 0)
        while i < len(chrom_cuts) and chrom_cuts[i][2] <= start:
            i += 1
        index[chrom] = i
        cursor = start
        j = i
        while j < len(chrom_cuts) and chrom_cuts[j][1] < end:
            if chrom_cuts[j][1] > cursor:
                out.append((chrom, cursor, chrom_cuts[j][1]))
            cursor = max(cursor, chrom_cuts[j][2])
            j += 1
        if cursor < end:
            out.append((chrom, cursor, end))
    return out


def derive(folder: Path) -> list[Path]:
    """Write the derived BEDs from the GIAB files in `folder` (each through a temp name)."""
    gaps = read_bed(folder / GAPS)
    low_map = read_bed(folder / LOW_MAPPABILITY)
    xdr = subtract(
        read_bed(folder / CHRY_NONPAR, "chrY"),
        [
            *read_bed(folder / CHRY_XTR),
            *read_bed(folder / CHRY_AMPLICONIC),
            *gaps,
            *low_map,
        ],
    )
    xnonpar = subtract(
        read_bed(folder / CHRX_NONPAR, "chrX"),
        [
            *read_bed(folder / CHRX_XTR),
            *read_bed(folder / CHRX_AMPLICONIC),
            *gaps,
            *low_map,
        ],
    )
    autosomes = subtract(read_bed(folder / AUTOSOMES), gaps)
    mappable = subtract(autosomes, low_map)
    if not xdr or not xnonpar or not mappable:
        raise ValueError(f"{folder}: a derived region is empty")
    return [
        _write(
            folder / SEX_REGIONS,
            [(i, "autosomes") for i in mappable]
            + [(i, "xdr") for i in xdr]
            + [(i, "xnonpar") for i in xnonpar],
        ),
        _write(folder / QC_AUTOSOMES, [(i, "autosomes") for i in autosomes]),
    ]


def _write(path: Path, rows: list[tuple[Interval, str]]) -> Path:
    tmp = path.with_name(f"{path.name}.tmp")
    with tmp.open("w") as f:
        for (chrom, start, end), name in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    os.replace(tmp, path)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.regions")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("downloads")
    sub.add_parser("derived")
    derive_parser = sub.add_parser("derive")
    derive_parser.add_argument("folder", type=Path)
    args = parser.parse_args(argv)
    if args.command == "downloads":
        print(*giab38_downloads(load_grading()), sep="\n")
    elif args.command == "derived":
        print(*DERIVED, sep="\n")
    else:
        for path in derive(args.folder):
            print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
