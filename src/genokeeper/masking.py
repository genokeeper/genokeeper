"""The GRC false duplications GIAB masks in its GRCh38 no-alt analysis set
(`GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2`), read from the FASTA
itself.

GIAB's `GRC_exclusions_T2Tv2.bed` (2023-02) post-dates the v2 FASTA (2021-07), so the mask is
not taken from a BED: `masked_intervals` compares the masked FASTA with NCBI's unmasked one
base by base, and every position the masked one sets to N is masked; any other difference
refuses. `derive` wrote the committed copy, `MASK_RESOURCE`; the bootstrap's `verify` derives
the mask again from the lab's two FASTAs and refuses unless it equals that copy, writing it
beside the FASTA (`masked_bed`), where the benchmark subtractions read it. So the mask every
check reads is the one in git.

`region_resources` and `overlaps` are a functional check, not a control: a region file under
`resources/` (whatever slice adds it) or a stage-3 window (`WINDOWS`) inside a masked block
would read zero depth, so the CI test refuses one.

    python -m genokeeper.masking verify UNMASKED MASKED OUT   # the BED, refused unless committed
    python -m genokeeper.masking derive UNMASKED MASKED       # writes MASK_RESOURCE
    python -m genokeeper.masking index FASTA                  # bwa index beside FASTA
"""

import argparse
import json
import os
import re
import secrets
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import pysam

from genokeeper.config import grch38_fasta
from genokeeper.pgx.hla import MHC_REGION
from genokeeper.regions import Interval, merge, read_bed, subtract
from genokeeper.repeats.expansionhunter import REGION_EXTENSION
from genokeeper.tools import ToolFailed, ToolMissing, gated, require, run_tool

# The masked FASTA as the bootstrap pins it (scripts/bootstrap.sh, masked_name and
# masked_sha256: the sha256 of the .gz fetched on 2026-10-30).
MASKED_NAME = (
    "GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"
)
MASKED_SHA256 = "74482c74f2a47e8908404028d5305167636c9942fba4129ea0d70c5f04cab5ad"
MASKED_BED = "GRCh38_maskedGRC.intervals.bed"
_PACKAGE = Path(__file__).resolve().parent
_REPO = _PACKAGE.parents[1]
RESOURCES = _PACKAGE / "resources"
MASK_RESOURCE = RESOURCES / "reference" / "masked_intervals.GRCh38.tsv"
WINDOW = 1_000_000  # bases compared at a time
BWA_EXTS = ("amb", "ann", "bwt", "pac", "sa")
# Stage 3's windows: SMN1 and SMN2 with Parascopy's model and SMNCopyNumberCaller's sites,
# CYP2D6 with CYP2D7/8, CYP2B6, and the MHC (HLA).
WINDOWS: dict[str, Interval] = {
    "smn": ("chr5", 69_000_000, 72_000_000),
    "cyp2d6": ("chr22", 42_100_000, 42_200_000),
    "cyp2b6": ("chr19", 40_980_000, 41_030_000),
    "hla": MHC_REGION,
}
_REGION_SUFFIXES = (".GRCh38.json", ".GRCh38.tsv")


def masked_bed(lab: Path) -> Path:
    """The mask the bootstrap derived from the lab's FASTAs, beside `grch38_fasta`."""
    return grch38_fasta(lab).with_name(MASKED_BED)


def without_mask(lab: Path, bed: Path, out: Path) -> Path:
    """`bed`'s intervals minus the lab's mask (`masked_bed`), written whole to `out`: a
    benchmark read on the masked reference, where a masked block has no reads."""
    _write(out, _bed(subtract(read_bed(bed), read_bed(masked_bed(lab)))))
    return out


def _first_difference(contig: str, start: int, a: bytes, b: bytes) -> ValueError:
    at = next(i for i, (x, y) in enumerate(zip(a, b, strict=True)) if x != y)
    return ValueError(
        f"{contig}:{start + at} (0-based): the masked FASTA has {chr(b[at])!r} where the"
        f" unmasked one has {chr(a[at])!r}; only a base set to N may differ"
    )


def masked_intervals(unmasked: Path, masked: Path) -> list[Interval]:
    """The positions `masked` sets to N where `unmasked` has a base (case ignored: soft
    masking is case, not sequence), merged into 0-based half-open runs. Both FASTAs'
    `.fai` must list the same contigs with the same lengths in the same order; they are
    read contig by contig in `WINDOW`-sized windows, so their line widths may differ.
    Raises `ValueError` on other contigs, or naming the first position where they differ
    other than by N."""
    found: list[Interval] = []
    with pysam.FastaFile(str(unmasked)) as u, pysam.FastaFile(str(masked)) as m:
        ours = list(zip(u.references, u.lengths, strict=True))
        theirs = list(zip(m.references, m.lengths, strict=True))
        if ours != theirs:
            first = next(
                (f"{a} against {b}" for a, b in zip(ours, theirs) if a != b),
                f"{len(ours)} contigs against {len(theirs)}",
            )
            raise ValueError(
                f"the FASTAs' contigs differ ({unmasked} against {masked}): {first}"
            )
        for contig, length in ours:
            for start in range(0, length, WINDOW):
                end = min(start + WINDOW, length)
                a = u.fetch(contig, start, end).upper().encode()
                b = m.fetch(contig, start, end).upper().encode()
                if a == b:
                    continue
                at = 0
                for run in re.finditer(rb"N+", b):
                    s, e = run.span()
                    if a[at:s] != b[at:s]:
                        raise _first_difference(contig, start + at, a[at:s], b[at:s])
                    for base in re.finditer(rb"[^N]+", a[s:e]):
                        found.append(
                            (contig, start + s + base.start(), start + s + base.end())
                        )
                    at = e
                if a[at:] != b[at:]:
                    raise _first_difference(contig, start + at, a[at:], b[at:])
    return merge(found)


def _read_intervals(text: str) -> list[Interval]:
    out: list[Interval] = []
    for line in text.splitlines():
        if not line or line.startswith(("#", "chrom\t")):
            continue
        chrom, start, end = line.split("\t")[:3]
        out.append((chrom, int(start), int(end)))
    return out


def load_mask() -> list[Interval]:
    """The committed mask, `MASK_RESOURCE`: 0-based half-open intervals."""
    return _read_intervals(MASK_RESOURCE.read_text())


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{secrets.token_hex(4)}")
    tmp.write_text(text)
    os.replace(tmp, path)


def _bed(intervals: Sequence[Interval]) -> str:
    return "".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)


def _catalog_regions(path: Path) -> list[Interval]:
    data = json.loads(path.read_text())
    if not isinstance(data, list) or not data:
        raise ValueError(f"{path}: not an ExpansionHunter catalog (a list of loci)")
    out: list[Interval] = []
    for entry in data:
        found = entry.get("ReferenceRegion") if isinstance(entry, dict) else None
        if found is None:
            raise ValueError(f"{path}: a locus without a ReferenceRegion: {entry!r}")
        for span in [found] if isinstance(found, str) else found:
            chrom, bounds = span.split(":")
            start, end = bounds.split("-")
            # ExpansionHunter reads REGION_EXTENSION bases on each side of the repeat.
            out.append(
                (chrom, int(start) - REGION_EXTENSION, int(end) + REGION_EXTENSION)
            )
    return out


def _table_regions(path: Path) -> list[Interval]:
    rows = [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    header = rows[0] if rows else []
    if {"chrom", "start", "end"} <= set(header):
        c, s, e = (header.index(k) for k in ("chrom", "start", "end"))
        return [(r[c], int(r[s]), int(r[e])) for r in rows[1:]]
    if {"chrom", "pos"} <= set(header):
        c, p = header.index("chrom"), header.index("pos")
        return [(r[c], int(r[p]) - 1, int(r[p])) for r in rows[1:]]
    raise ValueError(
        f"{path}: a GRCh38 region file needs chrom, start and end (0-based) or chrom and"
        f" pos (1-based) columns; it has {', '.join(header) or 'none'}"
    )


def region_resources(root: Path = RESOURCES) -> dict[str, list[Interval]]:
    """Every file under `root` whose name declares GRCh38 regions (`*.GRCh38.json`, an
    ExpansionHunter catalog, each region widened by `REGION_EXTENSION`; `*.GRCh38.tsv`
    with `chrom`, `start`, `end` or `chrom`, `pos` columns), the mask resource itself
    excluded, keyed by its path (from the repository root when under it). A declared file of
    neither shape raises `ValueError`: no region file is skipped silently."""
    found: dict[str, list[Interval]] = {}
    for path in sorted(root.rglob("*")):
        if not path.name.endswith(_REGION_SUFFIXES) or not path.is_file():
            continue
        if path.resolve() == MASK_RESOURCE.resolve():
            continue
        read = _catalog_regions if path.name.endswith(".json") else _table_regions
        resolved = path.resolve()
        key = (
            str(resolved.relative_to(_REPO))
            if resolved.is_relative_to(_REPO)
            else str(path)
        )
        found[key] = read(path)
    return found


def overlaps(
    masked: Sequence[Interval], regions: Mapping[str, Sequence[Interval]]
) -> dict[str, list[Interval]]:
    """Per region name, the masked intervals any of its intervals overlaps (empty: none)."""
    by_chrom: dict[str, list[Interval]] = {}
    for interval in masked:
        by_chrom.setdefault(interval[0], []).append(interval)
    out: dict[str, list[Interval]] = {}
    for name, intervals in regions.items():
        hit = {
            m
            for chrom, start, end in intervals
            for m in by_chrom.get(chrom, [])
            if m[1] < end and start < m[2]
        }
        out[name] = sorted(hit)
    return out


def _summary(intervals: Sequence[Interval]) -> str:
    return (
        f"{len(intervals)} masked intervals, {sum(e - s for _, s, e in intervals)} bp"
    )


def _verify(args: argparse.Namespace) -> int:
    found = masked_intervals(args.unmasked, args.masked)
    committed = load_mask()
    if found != committed:
        first = next(
            (f"{a} against {b}" for a, b in zip(found, committed) if a != b),
            "one is a prefix of the other",
        )
        raise ValueError(
            f"the mask derived from {args.masked} differs from the committed"
            f" {MASK_RESOURCE.name} ({_summary(found)} against {_summary(committed)};"
            f" first: {first})"
        )
    _write(args.out, _bed(found))
    print(f"{_summary(found)} -> {args.out}")
    return 0


def _derive(args: argparse.Namespace) -> int:
    found = masked_intervals(args.unmasked, args.masked)
    _write(
        MASK_RESOURCE,
        f"# source: {MASKED_NAME}\n# sha256: {MASKED_SHA256}\n"
        "# derived: the positions the masked FASTA sets to N where NCBI's no-alt analysis"
        " set has a base, merged; 0-based half-open\n"
        "chrom\tstart\tend\n" + _bed(found),
    )
    print(f"{_summary(found)} -> {MASK_RESOURCE}")
    return 0


def _index(args: argparse.Namespace) -> int:
    """`bwa index` of FASTA (gated `dna:align`) into hidden temporaries beside it, each
    renamed to `<FASTA>.<ext>` once all five are written."""
    require("bwa")
    fasta: Path = args.fasta
    for old in fasta.parent.glob(f".{fasta.name}.index-*"):
        old.unlink()
    prefix = fasta.with_name(f".{fasta.name}.index-{secrets.token_hex(4)}")
    try:
        run_tool(gated("dna:align", ["bwa", "index", "-p", str(prefix), str(fasta)]))
        for ext in BWA_EXTS:
            os.replace(f"{prefix}.{ext}", f"{fasta}.{ext}")
    finally:
        for ext in BWA_EXTS:
            Path(f"{prefix}.{ext}").unlink(missing_ok=True)
    print(f"bwa index -> {fasta}.{{{','.join(BWA_EXTS)}}}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.masking")
    sub = parser.add_subparsers(dest="command", required=True)
    verify = sub.add_parser("verify")
    derive = sub.add_parser("derive")
    for p in (verify, derive):
        p.add_argument("unmasked", type=Path)
        p.add_argument("masked", type=Path)
    verify.add_argument("out", type=Path)
    verify.set_defaults(handler=_verify)
    derive.set_defaults(handler=_derive)
    index = sub.add_parser("index")
    index.add_argument("fasta", type=Path)
    index.set_defaults(handler=_index)
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (ValueError, OSError, ToolFailed, ToolMissing) as e:
        print(f"masking {args.command}: {' '.join(str(e).split())}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
