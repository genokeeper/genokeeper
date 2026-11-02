"""DELLY 2.7.0's `cnv` mode and CNVpytor 1.3.2 on the CRAM, and Dysgu's duplication junctions,
read as `Deletion` events (slice 3d).

Both depth callers run offline in their images (`containers.IMAGES`) behind the core gate, on
the CRAM, its index and the FASTA mounted read-only and one output directory read-write, as
3b's callers do. `delly cnv` shards per chromosome exactly as 3b's `delly sr`
(`deletions.callers.delly_shards`): the chr20 shard found 46 of 46 of the whole-genome run's
chr20 PASS calls at reciprocal overlap 0.5 on HG002 (plan Task 6, 2026-10-30). Its exclude
list adds 3c's masked intervals to DELLY's template. CNVpytor's four commands run in one
container; its call table is its stdout.

A depth event is a `Deletion` of type `DEL` (a loss) or `DUP` (a gain), its deleted or
copied bases `[pos, end)` 0-based as 3b's deletions: DELLY's `<CNV>` record covers POS..END
and CNVpytor's region start..end (both 1-based inclusive), so `pos` is one less than the
first. DELLY writes only a copy number (`FORMAT/CN`), so a loss or a gain is known only
against the copy number expected there (`expected_ploidy`): with `--haploid-chr` a male's
normal chrX is written at CN 1 (260 such records on HG002), which that rule skips. CNVpytor
fits the sex chromosomes on their own, so its type column is already relative to the
sample's ploidy. Dysgu's junctions are the PASS `SVTYPE=DUP` records of 3b's raw
`dysgu.vcf.gz` (the plan's "Junctions": no second `delly call -t DUP`).
"""

from collections.abc import Callable, Collection, Sequence
from pathlib import Path

from genokeeper.config import Grading
from genokeeper.containers import run_container
from genokeeper.deletions.callers import (
    _inputs,
    delly_shards,
    read_records,
    vcf_records,
)
from genokeeper.deletions.screen import Deletion
from genokeeper.genotype.reliability import on_haploid_contig
from genokeeper.regions import read_bed
from genokeeper.tools import run_tool

CNV_LABEL = "dna:cnv"
CNVPYTOR_LABEL = "dna:cnvpytor"
# CNVpytor's rd ratio under which a deletion is read as both copies lost.
HOMOZYGOUS_BELOW = 0.25
_CNVPYTOR_TYPES = {"deletion": "DEL", "duplication": "DUP"}


def expected_ploidy(grading: Grading, *, male: bool) -> Callable[[str, int], int]:
    """The copy number expected at a 1-based GRCh38 position: 1 on a male's chrX or chrY
    outside the PARs (`[par]` in grading.toml, GRC's PAR1 and PAR2), 2 elsewhere."""

    def ploidy(chrom: str, pos: int) -> int:
        return 1 if on_haploid_contig(chrom, pos, grading=grading, male=male) else 2

    return ploidy


def run_delly_cnv(
    cram: Path,
    fasta: Path,
    out_dir: Path,
    *,
    threads: int,
    shards: Collection[str],
    haploid: Sequence[str],
    masked: Path,
    min_length: int,
) -> Path:
    """`delly cnv` on `cram` (GRCh38), its events from `min_length` bp (`-z`; `[cnv]`),
    `haploid` (a male's `("chrX", "chrY")`, a female's `()`) passed as `--haploid-chr`,
    sharded per chromosome in `shards` (`delly_shards`, `out_dir/dellycnv-shards/`, label
    `dna:cnv`) with every interval of the `masked` BED (3c's masked intervals) excluded;
    their BCFs concatenated in the `.fai`'s order and kept to FILTER PASS as
    `out_dir/dellycnv.vcf.gz`, bgzipped and indexed."""
    command = ["delly", "cnv", "-g", f"/ref/{fasta.name}", "-z", str(min_length)]
    if haploid:
        command += ["--haploid-chr", ",".join(haploid)]
    bcfs = delly_shards(
        cram,
        fasta,
        out_dir,
        threads=threads,
        shards=shards,
        folder="dellycnv-shards",
        command=command,
        label=CNV_LABEL,
        extra=[f"{c}\t{s}\t{e}" for c, s, e in read_bed(masked)],
    )
    merged, out = out_dir / "dellycnv.all.bcf", out_dir / "dellycnv.vcf.gz"
    run_tool(
        [
            "bcftools", "concat", "--no-version", "-Ob", "-o", str(merged),
            *map(str, bcfs),
        ]
    )  # fmt: skip
    run_tool(
        [
            "bcftools", "view", "--no-version", "-f", "PASS", "-Oz", "-o", str(out),
            str(merged),
        ]
    )  # fmt: skip
    run_tool(["bcftools", "index", "-f", "-t", str(out)])
    merged.unlink(missing_ok=True)
    return out


def read_delly_cnv(vcf: Path, *, ploidy: Callable[[str, int], int]) -> list[Deletion]:
    """Each FILTER PASS record of `vcf` whose copy number (`FORMAT/CN`) is off the one
    `ploidy` expects at its POS: below it a `DEL`, above it a `DUP`; `callers =
    {"dellycnv": "<ID> CN <cn>"}`, `gt` `"1/1"` at copy number 0, else `"0/1"`."""
    found: list[Deletion] = []
    for cols, info, values in vcf_records(vcf):
        chrom, pos, vid, ref, _, _, filt, _ = cols
        cn = values.get("CN", ".")
        if filt != "PASS" or not cn.isdigit():
            continue
        copies, expected = int(cn), ploidy(chrom, int(pos))
        if copies == expected:
            continue
        found.append(
            Deletion(
                chrom=chrom,
                pos=int(pos) - 1,
                end=int(info["END"]),
                ref=ref[0],
                gt="1/1" if copies == 0 else "0/1",
                callers={"dellycnv": f"{vid} CN {copies}"},
                depth_ratio=None,
                svtype="DEL" if copies < expected else "DUP",
            )
        )
    return found


def run_cnvpytor(
    cram: Path, fasta: Path, out_dir: Path, *, threads: int, bin_size: int
) -> Path:
    """CNVpytor on `cram` (GRCh38; it detects hg38 from the CRAM) in one container behind
    `dna:cnvpytor`, in `out_dir`: read depth (`-rd`, `threads` jobs), its histograms,
    partition and calls at `bin_size` (`[cnv]`); returns its call table, `calls.tsv` (the
    `-call` step's stdout; its stderr in `call.log`)."""
    b = str(bin_size)
    root = "cnvpytor -root s.pytor"
    script = " && ".join(
        [
            "cd /out",
            "rm -f s.pytor calls.tsv",
            f"{root} -rd /cram/{cram.name} -T /ref/{fasta.name} -j {threads}",
            f"{root} -his {b}",
            f"{root} -partition {b}",
            f"{root} -call {b} > calls.tsv 2> call.log",
        ]
    )
    run_container(
        "cnvpytor",
        ["bash", "-c", script],
        mounts=_inputs(cram, fasta, out_dir),
        label=CNVPYTOR_LABEL,
    )
    return out_dir / "calls.tsv"


def read_cnvpytor(path: Path, *, min_length: int) -> list[Deletion]:
    """CNVpytor's `deletion` and `duplication` rows of `path` (no header; type, region
    `chr:start-end` 1-based inclusive, length, rd ratio, e-val1, ...) at least `min_length`
    long; `callers = {"cnvpytor": "<region> e1 <e-val1>"}`, `gt` `"1/1"` for a deletion
    whose rd ratio is under `HOMOZYGOUS_BELOW`, else `"0/1"`."""
    found: list[Deletion] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        kind, region, _, ratio, e1 = line.split("\t")[:5]
        chrom, _, span = region.rpartition(":")
        start, end = (int(x) for x in span.split("-"))
        svtype = _CNVPYTOR_TYPES[kind]
        if end - start + 1 < min_length:
            continue
        found.append(
            Deletion(
                chrom=chrom,
                pos=start - 1,
                end=end,
                ref="N",
                gt="1/1"
                if svtype == "DEL" and float(ratio) < HOMOZYGOUS_BELOW
                else "0/1",
                callers={"cnvpytor": f"{region} e1 {e1}"},
                depth_ratio=None,
                svtype=svtype,
            )
        )
    return found


def read_junctions(vcf: Path) -> list[Deletion]:
    """Dysgu's FILTER PASS duplication records in `vcf` (3b's raw `dysgu.vcf.gz`), as
    `read_records` reads them: the junctions a depth gain is matched with."""
    found = read_records(vcf, caller="dysgu", svtype="DUP", pass_only=True)
    return [d for d, _ in found]
