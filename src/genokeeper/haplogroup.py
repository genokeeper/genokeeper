"""The mitochondrial haplogroup: Haplogrep3 on mutserve's homoplasmic chrM calls (track
design, slice 1d, "Haplogroup").

`write_chrm_vcf` writes the calls Haplogrep3 classifies: each mutserve SNV at a level of at
least [mito] min_homoplasmy, at a position at least [reliability] min_dp reads cover.
Haplogrep3 reads every ALT in a VCF as present, so a heteroplasmy is left out; and it reads
every position the VCF does not name as rCRS, so a position with fewer reads is listed
(`uncovered`) and named beside the haplogroup. When more than `MAX_UNCOVERED` of chrM is
uncovered, no haplogroup is stated (`sources.mito`). This control, against an uncovered
chrM read as rCRS (an all-uncovered one is classified H2a2a1): its anchor comes from
mutserve's own count of the reads at every position (`mito.chrm_coverage`), which
Haplogrep3 cannot move because the count is read before its VCF is written, and a position
the count lacks has no reads.

`run_haplogrep` runs Haplogrep3 in its image (`containers.IMAGES["haplogrep3"]`, built from
docker/haplogrep3.Dockerfile), offline (`--network none`) behind the core gate as
`dna:haplogroup`. Haplogrep3 fetches its trees from the network on first use, which the
run forbids, so the tree is chosen when the image is built and baked into it (`TREE`,
recorded in `run.tools`), never by `haplogrep3 trees` at run time.
"""

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from genokeeper.containers import run_container
from genokeeper.mito import CHRM, MitoCall

HAPLOGREP_VERSION = "3.3.2"
# The rule (track design): phylotree-rcrs@17.3, or the forensic-update tree if Haplogrep3
# lists a newer one. Haplogrep3 3.3.2's list (its haplogrep3.yaml and the trees.json of
# github.com/genepi/haplogrep-trees, read 2026-10-21) has phylotree-fu-rcrs@1.3, released
# 2026-06-08 16:08 UTC, before phylotree-rcrs@17.3 (16:55 the same day): not newer.
TREE = "phylotree-rcrs@17.3"
LABEL = "dna:haplogroup"
# rCRS, GRCh38's chrM.
RCRS_LENGTH = 16569
# Above this share of chrM below min_dp no haplogroup is stated (builder's call): at 30x
# WGS chrM is covered hundreds of times over, so more than 1% uncovered (166 positions)
# means the reads were lost, not that a few positions were hard.
MAX_UNCOVERED = 0.01
VCF_NAME = "chrM.vcf"
OUTPUT = "haplogroups.txt"
_JAR = "/opt/haplogrep3/haplogrep3.jar"


@dataclass(frozen=True)
class Haplogroup:
    """Haplogrep3's haplogroup and its quality (0-1: how well the sample's mutations fit
    the haplogroup's)."""

    name: str
    quality: float


def uncovered(coverage: Mapping[int, int], *, min_dp: int) -> list[int]:
    """The chrM positions fewer than `min_dp` reads cover (`mito.chrm_coverage`), in
    order; a position the coverage lacks has none."""
    return [p for p in range(1, RCRS_LENGTH + 1) if coverage.get(p, 0) < min_dp]


def too_uncovered(positions: Sequence[int]) -> bool:
    """Whether `uncovered` positions are too many for a haplogroup (`MAX_UNCOVERED`)."""
    return len(positions) > MAX_UNCOVERED * RCRS_LENGTH


def write_chrm_vcf(
    calls: Sequence[MitoCall],
    coverage: Mapping[int, int],
    out: Path,
    *,
    min_homoplasmy: float,
    min_dp: int,
) -> Path:
    """Write the homoplasmic SNVs of `calls` at covered positions (see the module
    docstring) as a one-sample VCF on chrM, each one copy (`1`). `out` holds the user's
    chrM calls: the caller puts it under the run's work directory."""
    rows = [
        f"{CHRM}\t{c.pos}\t.\t{c.ref}\t{c.alt}\t.\tPASS\t.\tGT\t1\n"
        for c in sorted(calls, key=lambda c: c.pos)
        if c.level >= min_homoplasmy
        and coverage.get(c.pos, 0) >= min_dp
        and len(c.ref) == len(c.alt) == 1
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CHRM},length={RCRS_LENGTH}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(rows)
    )
    return out


def read_haplogroup(path: Path) -> Haplogroup:
    """The haplogroup and quality of Haplogrep3's one sample in its output at `path`;
    `ValueError` for an output that holds not exactly one sample."""
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: {len(rows)} samples, expected one sample")
    return Haplogroup(rows[0]["Haplogroup"], float(rows[0]["Quality"]))


def run_haplogrep(vcf: Path, out_dir: Path, *, tree: str) -> Haplogroup:
    """Haplogrep3's haplogroup of `vcf` on `tree` (one the image holds), its output written
    into `out_dir`."""
    run_container(
        "haplogrep3",
        [
            "java", "-jar", _JAR, "classify", "--in", f"/in/{vcf.name}",
            "--tree", tree, "--out", f"/out/{OUTPUT}",
        ],
        mounts=[(vcf.parent, "/in", True), (out_dir, "/out", False)],
        label=LABEL,
    )  # fmt: skip
    return read_haplogroup(out_dir / OUTPUT)
