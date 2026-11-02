"""ExpansionHunter (5.0.0 built at Illumina's commit 9e93322f4a, the last Apache-2.0 one,
which has the RFC1 motif analysis) and REViewer 0.2.7 on the repeat loci (slices 3a, 3e).

Both run offline in their images (`containers.IMAGES`, built from docker/) behind the core
gate as `dna:repeats`, on the CRAM, its index and the FASTA (with its `.fai`) mounted
read-only and one output directory mounted read-write. The variant catalog in
`resources/repeats/` is the cut of Illumina's ExpansionHunter v5.0.0 GRCh38 catalog merged
with the entries it lacks (`added.GRCh38.json`), plus RFC1 and CFTR (`cut_catalog`,
`SOURCE.md`); the thresholds come from STRchive (`[repeat_loci]` in config/grading.toml).
"""

import hashlib
import json
import re
import sys
import tempfile
from collections.abc import Iterable, Sequence
from importlib.resources import as_file, files
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any

from genokeeper.containers import run_container
from genokeeper.repeats.loci import Allele, RepeatCall, RepeatLocus
from genokeeper.sex import Sex
from genokeeper.tools import run_tool

EH_VERSION = "5.0.0 at 9e93322f4a"
REVIEWER_VERSION = "0.2.7"
LABEL = "dna:repeats"
CATALOG_NAME = "stage3a.GRCh38.json"
CATALOG: Traversable = files("genokeeper").joinpath(
    "resources", "repeats", CATALOG_NAME
)
ADDED: Traversable = files("genokeeper").joinpath(
    "resources", "repeats", "added.GRCh38.json"
)
# Illumina's entry these replace with the added file's (PABPN1: Illumina counts the pure GCG
# run, the disease counts every GCN).
REPLACED = frozenset({"PABPN1"})
# Entries graded apart from `[repeat_loci]`, closing the catalog: RFC1 (Illumina's, with
# the motif analysis on) and CFTR (the added file's).
EXTRA = ("RFC1", "CFTR")
# ExpansionHunter's own default, which REViewer must be given to draw the reads
# ExpansionHunter realigned.
REGION_EXTENSION = 1000
_REALIGNED = "eh_realigned.bam"
SORTED_BAM = "eh_realigned.sorted.bam"
# C9orf72's flanks are read for the GC-flank ratio, this far each side of the repeat.
FLANK = 500


def _chr(region: str) -> str:
    return region if region.startswith("chr") else f"chr{region}"


def _entry(e: dict[str, Any]) -> dict[str, Any]:
    """An Illumina entry without `OfftargetRegions`, its contigs `chr` prefixed."""
    found = e["ReferenceRegion"]
    regions = _chr(found) if isinstance(found, str) else [_chr(r) for r in found]
    return {k: v for k, v in e.items() if k != "OfftargetRegions"} | {
        "ReferenceRegion": regions
    }


def cut_catalog(
    source: Path, loci: Iterable[RepeatLocus], *, added: Path
) -> list[dict[str, Any]]:
    """The catalog: the entries of Illumina's ExpansionHunter catalog at `source` whose
    `LocusId` is a locus's, in file order and unchanged but for `OfftargetRegions` (left
    out) and each `ReferenceRegion`'s contig (`chr` prefixed), a `REPLACED` one swapped in
    place for `added`'s entry; then `added`'s other entries of the loci, in its order; then
    Illumina's RFC1 entry with its motif analysis on, and `added`'s CFTR entry.

    A locus in neither file, a `LocusId` in both but not in `REPLACED`, a `REPLACED` locus
    missing from `added` and a missing `EXTRA` entry are errors."""
    # The merge never lets one file silently shadow the other. Its anchor comes from
    # `REPLACED` and `EXTRA`, this module's constants, which an edited catalog file cannot
    # move because they are code, not read from either file.
    wanted = {locus.locus_id for locus in loci}
    illumina: list[dict[str, Any]] = json.loads(source.read_text())
    extra: list[dict[str, Any]] = json.loads(added.read_text())
    ours = {e["LocusId"]: e for e in extra}
    theirs = {e["LocusId"]: e for e in illumina}
    rfc1, cftr = EXTRA
    for names, why in (
        ((REPLACED & wanted) - set(ours), f"in REPLACED but not in {added}"),
        (
            (set(ours) & set(theirs)) - REPLACED,
            f"in {source} and {added}, not REPLACED",
        ),
        (wanted - set(ours) - set(theirs), f"a locus in neither {source} nor {added}"),
        ({rfc1} - set(theirs), f"not in {source}"),
        ({cftr} - set(ours), f"not in {added}"),
    ):
        if names:
            raise ValueError(f"{', '.join(sorted(names))}: {why}")
    cut = [
        ours[e["LocusId"]] if e["LocusId"] in REPLACED else _entry(e)
        for e in illumina
        if e["LocusId"] in wanted
    ]
    cut += [e for e in extra if e["LocusId"] in wanted and e["LocusId"] not in REPLACED]
    return [*cut, _entry(theirs[rfc1]) | {"RFC1MotifAnalysis": True}, ours[cftr]]


def catalog_sha256() -> str:
    """The sha256 of the catalog resource ExpansionHunter and REViewer are given."""
    return hashlib.sha256(CATALOG.read_bytes()).hexdigest()


def region(locus: RepeatLocus) -> tuple[str, int, int]:
    """The locus's repeat in GRCh38 (1-based, inclusive), from its catalog entry's
    `ReferenceRegion` (0-based start), the one its `VariantId` (or, for a one-repeat entry
    without it, its `LocusId`) names `locus.variant_id`."""
    entry = next(
        e for e in json.loads(CATALOG.read_text()) if e["LocusId"] == locus.locus_id
    )
    found = entry["ReferenceRegion"]
    regions = [found] if isinstance(found, str) else found
    ids = entry.get("VariantId") or [locus.locus_id]
    (span,) = (r for r, v in zip(regions, ids, strict=True) if v == locus.variant_id)
    chrom, bounds = span.split(":")
    start, end = bounds.split("-")
    return chrom, int(start) + 1, int(end)


def _inputs(cram: Path, fasta: Path, out_dir: Path) -> list[tuple[Path, str, bool]]:
    return [
        (cram.parent, "/cram", True),
        (fasta.parent, "/ref", True),
        (out_dir, "/out", False),
    ]


def run_expansionhunter(
    cram: Path, fasta: Path, out_dir: Path, *, sex: Sex, threads: int
) -> Path:
    """ExpansionHunter's VCF (`out_dir/eh.vcf`) of the catalog's loci on `cram`; its JSON
    and its realigned reads, sorted and indexed natively (`SORTED_BAM`), land beside it."""
    with as_file(CATALOG) as catalog:
        run_container(
            "expansionhunter",
            [
                "ExpansionHunter", "--reads", f"/cram/{cram.name}",
                "--reference", f"/ref/{fasta.name}",
                "--variant-catalog", "/catalog.json",
                "--output-prefix", "/out/eh", "--sex", sex,
                "--threads", str(threads), "--analysis-mode", "seeking",
            ],
            mounts=[*_inputs(cram, fasta, out_dir), (catalog, "/catalog.json", True)],
            label=LABEL,
        )  # fmt: skip
    bam, sorted_bam = out_dir / _REALIGNED, out_dir / SORTED_BAM
    run_tool(["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), str(bam)])
    run_tool(["samtools", "index", "-@", str(threads), str(sorted_bam)])
    return out_dir / "eh.vcf"


def _fields(text: str, count: int, where: str) -> list[str]:
    parts = text.split("/")
    if len(parts) != count:
        raise ValueError(f"{where}: {text!r} has {len(parts)} values, expected {count}")
    return parts


def _pairs(counts: str) -> list[tuple[int, int]]:
    """The (units, reads) pairs of a `CountsOf…Reads` string such as
    `"(1, 1), (6, 1), (25, 1)"`; none for `"()"`."""
    return [(int(u), int(r)) for u, r in re.findall(r"\((\d+), (\d+)\)", counts)]


def _longest(counts: str) -> int:
    """The largest unit count with a nonzero read count in a `CountsOf…Reads` string."""
    return max((units for units, reads in _pairs(counts) if reads > 0), default=0)


def read_calls(vcf: Path, json_path: Path) -> dict[str, RepeatCall]:
    """ExpansionHunter's calls in `vcf`, keyed by `VARID`; `json_path` (its JSON) gives the
    read length, the longest flanking read and the in-repeat reads (the VCF's `ADIR` counts
    only those consistent with an allele at least a read long, so it is 0 for every allele
    shorter than the read). A variant whose `REPCN` is `.` is left out."""
    results = json.loads(json_path.read_text())["LocusResults"].values()
    read_length = {
        variant: locus["ReadLength"]
        for locus in results
        for variant in locus["Variants"]
    }
    flanking = {
        variant: _longest(found.get("CountsOfFlankingReads", "()"))
        for locus in results
        for variant, found in locus["Variants"].items()
    }
    inrepeat = {
        variant: sum(r for _, r in _pairs(found.get("CountsOfInrepeatReads", "()")))
        for locus in results
        for variant, found in locus["Variants"].items()
    }
    calls: dict[str, RepeatCall] = {}
    for line in vcf.read_text().splitlines():
        if line.startswith("#"):
            continue
        row = line.split("\t")
        info = dict(kv.split("=", 1) for kv in row[7].split(";") if "=" in kv)
        varid = info["VARID"]
        sample = dict(zip(row[8].split(":"), row[9].split(":"), strict=True))
        if "." in sample["REPCN"]:
            continue
        count = len(sample["GT"].split("/"))
        counts = [
            _fields(sample[key], count, f"{varid} {key}")
            for key in ("REPCN", "REPCI", "ADSP", "ADFL", "ADIR")
        ]
        alleles = tuple(
            Allele(
                int(n),
                (int(ci.split("-")[0]), int(ci.split("-")[1])),
                int(sp),
                int(fl),
                int(ir),
            )
            for n, ci, sp, fl, ir in zip(*counts, strict=True)
        )
        calls[varid] = RepeatCall(
            variant_id=varid,
            alleles=alleles,
            coverage=float(sample["LC"]),
            low_depth="LowDepth" in row[6].split(";"),
            read_length=read_length[varid],
            longest_flanking=flanking[varid],
            inrepeat=inrepeat[varid],
        )
    return calls


def review(
    out_dir: Path, *, vcf: Path, fasta: Path, locus_ids: Sequence[str]
) -> dict[str, Path]:
    """One REViewer run over `locus_ids` on ExpansionHunter's `vcf` (inside `out_dir`, with
    its sorted realigned reads); the image of each locus REViewer drew, keyed by locus."""
    with as_file(CATALOG) as catalog:
        run_container(
            "reviewer",
            [
                "REViewer", "--reads", f"/out/{SORTED_BAM}",
                "--vcf", f"/out/{vcf.relative_to(out_dir)}",
                "--reference", f"/ref/{fasta.name}",
                "--catalog", "/catalog.json", "--locus", ",".join(locus_ids),
                "--region-extension-length", str(REGION_EXTENSION),
                "--output-prefix", "/out/rev",
            ],
            mounts=[
                (fasta.parent, "/ref", True),
                (out_dir, "/out", False),
                (catalog, "/catalog.json", True),
            ],
            label=LABEL,
        )  # fmt: skip
    images = {locus_id: out_dir / f"rev.{locus_id}.svg" for locus_id in locus_ids}
    return {locus_id: path for locus_id, path in images.items() if path.is_file()}


def flank_ratio(
    cram: Path, fasta: Path, *, region: tuple[str, int, int], mean_depth: float
) -> float:
    """The depth over the `FLANK` bases each side of `region` (1-based, inclusive: the
    repeat) per base, over `mean_depth`."""
    chrom, start, end = region
    bed = [(chrom, start - FLANK - 1, start - 1), (chrom, end, end + FLANK)]
    with tempfile.TemporaryDirectory() as scratch:
        regions = Path(scratch) / "flanks.bed"
        regions.write_text("".join(f"{c}\t{a}\t{b}\n" for c, a, b in bed))
        out = run_tool(
            ["samtools", "bedcov", "--reference", str(fasta), str(regions), str(cram)]
        )
    rows = [line.split("\t") for line in out.splitlines()]
    bases = sum(int(r[2]) - int(r[1]) for r in rows)
    return sum(int(r[3]) for r in rows) / bases / mean_depth


def main(argv: Sequence[str] | None = None) -> int:
    """`python -m genokeeper.repeats.expansionhunter cut FILE`: write the catalog resource
    from Illumina's ExpansionHunter catalog at FILE and the added entries (`ADDED`)."""
    from genokeeper.config import load_grading

    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2 or args[0] != "cut":
        print(
            "usage: python -m genokeeper.repeats.expansionhunter cut FILE",
            file=sys.stderr,
        )
        return 2
    with as_file(ADDED) as added:
        cut = cut_catalog(
            Path(args[1]), load_grading().repeats.loci.values(), added=added
        )
    target = (
        Path(__file__).resolve().parents[1] / "resources" / "repeats" / CATALOG_NAME
    )
    target.write_text(json.dumps(cut, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
