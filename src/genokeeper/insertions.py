"""Insertions of 50 bp or more in the listed genes' exons (slice 3d): Dysgu 1.9.0 and
INSurVeyor 1.1.3, matched by Truvari, with Scramble 1.0.2's mobile elements.

INSurVeyor and Scramble run offline in their images (`containers.IMAGES`) behind the core
gate, on the CRAM, its index and the FASTA mounted read-only and one output directory
read-write, as 3b's callers do; bcftools and samtools run natively (`run_tool`). Dysgu's
insertions come from 3b's own Dysgu run. The two callers' union is 3b's
(`deletions_only`, `match_callers`, `read_union`, of type `INS`), screened by 3b's
`screen_all` with `INSERTION_CALLERS` and `[insertions]` (`Grading.insertions`) and stated
by 3b's `DeletionSource` in `INSERTION_WORDING`: no screen or source of its own.

Scramble's image has no samtools, so `run_scramble` cuts the target exons' reads natively
first (`samtools view -M`, the index's multi-region iterator: without it samtools streams the
whole CRAM). A Scramble call is an insertion at a point whose length is unknown (a lower
bound, C20) and whose genotype it does not give (one copy, `0/1`, as a depth caller's call
off a homozygous loss). Within `[insertions] max_distance` of a union call it labels that call
with its family, which keeps its class (`screen.LABELS`); any other is an event of its own,
band `mei`.
"""

import shutil
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path

from genokeeper.containers import run_container
from genokeeper.deletions.callers import _inputs, deletions_only, sort_index
from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import Deletion
from genokeeper.deletions.truvari import match_callers, read_union
from genokeeper.regions import merge
from genokeeper.sources.deletions import DYSGU_VERSION, TRUVARI_VERSION, Wording
from genokeeper.tools import run_tool

INSURVEYOR_VERSION = "1.1.3"
SCRAMBLE_VERSION = "1.0.2"
INSURVEYOR_LABEL = "dna:insurveyor"
SCRAMBLE_LABEL = "dna:scramble"
INSERTION_CALLERS = ("dysgu", "insurveyor")
INSERTIONS_TOOL = "insertions"
# bp each side of a screen gene's exon whose reads Scramble reads (the plan's Task 10).
SCRAMBLE_MARGIN = 500
# Inside Scramble's image (docker/scramble.Dockerfile): SCRAMble.R's folder and the MEI
# consensus sequences, as its upstream tag lays them out.
SCRAMBLE_BIN = "/opt/scramble/cluster_analysis/bin"
MEI_REFS = "/opt/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa"
_PREFIX = "scramble"
# Scramble's `MEI_Family` values, as the reason names them.
_FAMILIES = {"alu": "Alu", "l1": "L1", "sva": "SVA"}
MEI_SCOPE = (
    "mobile elements labelled by Scramble on the target exons only; MELT and xTea are not"
    " open to a private user"
)


def mei_family(d: Deletion) -> str | None:
    """The mobile-element family Scramble named for the call (`Alu`, `L1`, `SVA`), if any."""
    label = d.callers.get("scramble")
    return None if label is None else label.split(" ")[0]


INSERTION_WORDING = Wording(
    tool=INSERTIONS_TOOL,
    noun="insertion",
    claim_type="insertion",
    # R3: MLPA misreads a mobile-element insertion as a deletion.
    confirm_by="needs confirmation by long reads or PCR across the insertion",
    alt="<INS>",
    source_level=(
        f"Dysgu {DYSGU_VERSION} and INSurVeyor {INSURVEYOR_VERSION}, matched by Truvari"
        f" {TRUVARI_VERSION}; Scramble {SCRAMBLE_VERSION} for mobile elements"
    ),
    source_version=(
        f"Dysgu {DYSGU_VERSION}; INSurVeyor {INSURVEYOR_VERSION};"
        f" Truvari {TRUVARI_VERSION}; Scramble {SCRAMBLE_VERSION}"
    ),
    callers=INSERTION_CALLERS,
    qualifier=mei_family,
)


def run_insurveyor(cram: Path, fasta: Path, out_dir: Path, *, threads: int) -> Path:
    """INSurVeyor on `cram` (GRCh38) with `threads`, its work directory `out_dir/work`,
    which the tool does not create; its PASS calls (`out.pass.vcf.gz`, every record
    `SVTYPE=INS`) sorted, bgzipped and indexed natively as `out_dir/insurveyor.vcf.gz`."""
    (out_dir / "work").mkdir(exist_ok=True)
    run_container(
        "insurveyor",
        [
            "insurveyor.py", "--threads", str(threads), f"/cram/{cram.name}", "/out/work",
            f"/ref/{fasta.name}",
        ],
        mounts=_inputs(cram, fasta, out_dir),
        label=INSURVEYOR_LABEL,
    )  # fmt: skip
    return sort_index(
        out_dir / "work" / "out.pass.vcf.gz", out_dir / "insurveyor.vcf.gz"
    )


def scramble_regions(targets: Mapping[str, Sequence[Target]], out: Path) -> Path:
    """The BED `run_scramble` reads: `targets` (the screen genes' exons, built with margin
    `SCRAMBLE_MARGIN`) merged and sorted."""
    rows = merge((t.chrom, t.start, t.end) for ts in targets.values() for t in ts)
    out.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))
    return out


def run_scramble(
    cram: Path, fasta: Path, out_dir: Path, *, regions: Path, threads: int
) -> Path:
    """Scramble's mobile-element insertions in `cram` over `regions` (`scramble_regions`):
    the reads there cut natively into `out_dir/exons.bam` and indexed, then
    `cluster_identifier` and `SCRAMble.R --eval-meis` (no deletions) in one container behind
    `dna:scramble`; its MEI table moved to `out_dir/meis.txt`."""
    bam = out_dir / "exons.bam"
    run_tool(
        [
            "samtools", "view", "-M", "-b", "-@", str(threads), "-L", str(regions),
            "-T", str(fasta), "-o", str(bam), str(cram),
        ]
    )  # fmt: skip
    run_tool(["samtools", "index", str(bam)])
    rscript = " ".join(
        [
            "Rscript --vanilla SCRAMble.R",
            f"--out-name /out/{_PREFIX}",
            "--cluster-file /out/clusters.txt",
            f"--install-dir {SCRAMBLE_BIN}",
            f"--mei-refs {MEI_REFS}",
            f"--ref /ref/{fasta.name}",
            "--eval-meis",
        ]
    )
    script = " && ".join(
        [
            f"cluster_identifier /out/{bam.name} > /out/clusters.txt",
            f"cd {SCRAMBLE_BIN}",
            rscript,
        ]
    )
    run_container(
        "scramble",
        ["bash", "-c", script],
        mounts=_inputs(cram, fasta, out_dir),
        label=SCRAMBLE_LABEL,
    )
    meis = out_dir / "meis.txt"
    shutil.move(out_dir / f"{_PREFIX}_MEIs.txt", meis)
    return meis


def read_scramble(path: Path) -> list[Deletion]:
    """One `INS` per row of Scramble's MEI table (a header-only file holds none): its
    `Insertion` `chr:p` (0-based p, so the VCF base before it is 1-based p) as `pos`, `end =
    pos + 1`, length unknown (`svlen` None, `lower_bound`), `callers = {"scramble":
    "<family> <strand> <clipped reads>"}`. A family other than Alu, L1 or SVA is refused:
    its anchor comes from Scramble's MEI consensus file, which names only those three, so a
    table from another tool or version cannot pass unread."""
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    column = {name: i for i, name in enumerate(header)}
    found: list[Deletion] = []
    for line in lines[1:]:
        fields = line.split("\t")
        chrom, _, at = fields[column["Insertion"]].rpartition(":")
        family = fields[column["MEI_Family"]]
        if family not in _FAMILIES:
            raise ValueError(f"{path}: unknown MEI family {family!r}")
        strand = fields[column["Insertion_Direction"]]
        reads = fields[column["Clipped_Reads_In_Cluster"]]
        pos = int(at)
        found.append(
            Deletion(
                chrom=chrom,
                pos=pos,
                end=pos + 1,
                ref="N",
                gt="0/1",
                callers={"scramble": f"{_FAMILIES[family]} {strand} {reads}"},
                depth_ratio=None,
                svtype="INS",
                svlen=None,
                lower_bound=True,
            )
        )
    return found


def insertion_union(
    dysgu_ins: Path, insurveyor: Path, out_dir: Path, *, min_length: int
) -> list[Deletion]:
    """The Dysgu+INSurVeyor union of insertions of `min_length` bp or more: INSurVeyor's
    calls kept as 3b keeps a caller's (`deletions_only`, `out_dir/insurveyor.ins.vcf.gz`),
    matched against Dysgu's (`dysgu_ins`, already kept the same way) by Truvari in
    `out_dir/both`, read back as 3b's union."""
    kept = deletions_only(
        insurveyor,
        out_dir / "insurveyor.ins.vcf.gz",
        caller="insurveyor",
        min_length=min_length,
        svtype="INS",
    )
    folder = match_callers(
        dysgu_ins, kept, out_dir, min_length=min_length, names=INSERTION_CALLERS
    )
    return read_union(folder, INSERTION_CALLERS, svtype="INS")


def add_scramble(
    union: Sequence[Deletion], meis: Sequence[Deletion], *, max_distance: int
) -> list[Deletion]:
    """`union` with each call that a Scramble MEI lies within `max_distance` of labelled
    (its `callers` gain the nearest such MEI's `scramble`; its class is unchanged), then
    every MEI within `max_distance` of no union call, as an event of its own."""
    nearest: dict[int, tuple[int, Deletion]] = {}
    alone: list[Deletion] = []
    for mei in meis:
        near = [
            (abs(d.pos - mei.pos), i)
            for i, d in enumerate(union)
            if d.chrom == mei.chrom and abs(d.pos - mei.pos) <= max_distance
        ]
        if not near:
            alone.append(mei)
            continue
        distance, i = min(near)
        if i not in nearest or distance < nearest[i][0]:
            nearest[i] = (distance, mei)
    labelled = [
        d
        if i not in nearest
        else replace(d, callers={**d.callers, **nearest[i][1].callers})
        for i, d in enumerate(union)
    ]
    return [*labelled, *alone]
