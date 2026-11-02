"""HG003's gVCF scored against GIAB's truth: errors per cell (`genokeeper.calibrate.table`).

What counts as an error (track design, stage 2 "2b. Calibration"):

- **Variant records** (PASS): `rtg vcfeval --output-mode annotate` inside the benchmark BED,
  run twice. Genotype-exact, a call annotated `FP` (no such allele in the truth) or `FP_CA`
  (the allele, another genotype) is a genotype error; with `--squash-ploidy`, which ignores
  zygosity, an `FP` is an allele error. `TP` is right; records vcfeval skips (`IGN`: not
  PASS, or no ALT in the genotype) or finds outside the BED (`OUT`) are not counted. The
  allele class is the called alleles': an indel if any changes length.
- **RefCalls and reference blocks**: vcfeval skips records whose genotype has no ALT, so
  these are counted here by interval arithmetic. The units are their reference positions
  inside the benchmark BED (a block from POS to END, a RefCall over its REF), and an error
  is a truth variant there: an SNV at its first changed base, a deletion at its first
  deleted base, an insertion at the base before it. A position can hold a claim of either
  class, so each position counts once in the SNV cell and once in the indel cell, and a
  truth SNV is an error of the SNV cell, a truth indel of the indel cell. The two error
  kinds coincide here. A no-call block or RefCall, and a block over no reads (MIN_DP 0),
  are no calls (`genokeeper.genotype.gvcf`) and are not counted; a block's GQ is capped at
  `BLOCK_GQ_CAP` as the call layer caps it.

The region is "difficult" where any of the grading's GIAB v3.6 strata (`StrataIndex`) holds
the position, "easy" elsewhere. Counts are kept per chromosome (`count_cells`), so the
held-out checks can build a table from any set of chromosomes without scoring again.

rtg-tools runs natively, not in a container: it reads HG003's calls, a public sample, never
the user's; the bootstrap's rtg.cfg turns off its crash reports and usage logging.
"""

import functools
import gzip
import os
import shutil
import subprocess
import tempfile
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Iterator, Mapping
from pathlib import Path

from genokeeper.calibrate import rtg_path
from genokeeper.calibrate.table import (
    ALLELE,
    ALLELE_CLASSES,
    GENOTYPE,
    KINDS,
    CalibrationTable,
    CellKey,
    Counts,
    gq_bin,
)
from genokeeper.config import load_paths
from genokeeper.genotype.gvcf import BLOCK_GQ_CAP
from genokeeper.genotype.reliability import StrataIndex
from genokeeper.tools import ToolFailed, ToolMissing, gated, run_tool
from genokeeper.variants.key import ucsc_contig

type Cells = dict[tuple[str, CellKey], Counts]
# Merged, sorted 0-based half-open interval starts and ends of one contig.
type _Intervals = tuple[list[int], list[int]]

# vcfeval's annotations (annotate mode) of a scored call, and which are errors per kind.
_SCORED = {"TP", "FP", "FP_CA"}
_ERRORS = {GENOTYPE: {"FP", "FP_CA"}, ALLELE: {"FP"}}
_SYMBOLIC = "<*>"


def measure(
    gvcf: Path,
    truth_vcf: Path,
    bench_bed: Path,
    *,
    sdf: Path,
    strata: StrataIndex,
    versions: Mapping[str, str],
    rtg: Path | None = None,
    threads: int | None = None,
) -> CalibrationTable:
    """The calibration table of `gvcf` (DeepVariant's, on the reference `sdf` holds)
    against `truth_vcf` inside `bench_bed`, with the pipeline's `versions`. `rtg` defaults
    to the lab's (`rtg_path`), `threads` to this machine's CPUs."""
    total: Cells = {}
    by_chrom, f1 = count_cells(
        gvcf, truth_vcf, bench_bed, sdf=sdf, strata=strata, rtg=rtg, threads=threads
    )
    for cells in by_chrom.values():
        for key, found in cells.items():
            total.setdefault(key, Counts()).merge(found)
    return CalibrationTable.from_counts(total, versions, f1=f1, by_chrom=by_chrom)


def count_cells(
    gvcf: Path,
    truth_vcf: Path,
    bench_bed: Path,
    *,
    sdf: Path,
    strata: StrataIndex,
    rtg: Path | None = None,
    threads: int | None = None,
) -> tuple[dict[str, Cells], dict[str, float]]:
    """Errors and units per (error kind, cell), per chromosome (see the module), and the
    genotype-exact F1 of the SNPs and of the INDELs (every other call) over all PASS
    calls: the last row, at the lowest GQ, of vcfeval's snp and non-snp ROCs."""
    f1: dict[str, float] = {}
    rtg = rtg if rtg is not None else rtg_path(load_paths().lab)
    if not rtg.is_file():
        raise FileNotFoundError(f"no rtg at {rtg}: run scripts/bootstrap.sh")
    threads = threads if threads is not None else (os.cpu_count() or 4)
    difficult = _union(strata)
    bench = _read_bed(bench_bed)
    by_chrom: dict[str, Cells] = {}
    work = Path(tempfile.mkdtemp(prefix=".vcfeval-", dir=gvcf.parent))
    try:
        # vcfeval scores PASS records only: it is handed those alone, not the blocks.
        variants = work / "variants.vcf.gz"
        run_tool(
            gated(
                "dna:calibrate",
                ["bcftools", "view", "--no-version", "-f", "PASS", "-Oz"]
                + ["-o", str(variants), "--write-index=tbi", str(gvcf)],
            )
        )
        for kind in KINDS:
            out = work / kind
            _vcfeval(
                rtg,
                truth_vcf,
                variants,
                bench_bed,
                sdf,
                out,
                squash=kind == ALLELE,
                threads=threads,
            )
            _count_variants(out / "calls.vcf.gz", kind, difficult, by_chrom)
        for allele_type, roc in (("snp", "snp"), ("indel", "non_snp")):
            with gzip.open(work / GENOTYPE / f"{roc}_roc.tsv.gz", "rt") as f:
                rows = [line.split("\t") for line in f if not line.startswith("#")]
            # Its f_measure column is rounded to 4 places: F1 from the counts instead.
            tp_base, fp, tp_call, fn = map(float, rows[-1][1:5]) if rows else (0,) * 4
            f1[allele_type] = (
                2
                * tp_call
                * tp_base
                / (tp_call * (tp_base + fn) + tp_base * (tp_call + fp))
                if tp_call and tp_base
                else 0.0
            )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    _count_reference(gvcf, _truth_positions(truth_vcf), bench, difficult, by_chrom)
    return by_chrom, f1


_AUTOSOMES = frozenset(f"chr{n}" for n in range(1, 23))


def true_het_vafs(
    gvcf: Path,
    truth_vcf: Path,
    bench_bed: Path,
    *,
    sdf: Path,
    rtg: Path,
    threads: int,
) -> list[float]:
    """The FORMAT/VAF of each true heterozygous SNV call of `gvcf` on chr1-22: its PASS
    records scored by genotype-exact vcfeval (as `count_cells`) inside `bench_bed`, the TP
    records whose genotype holds the reference and one ALT, an SNV. VAF is DeepVariant's
    AD[alt] / DP, the quantity `reliability.guard_allele_fraction` compares; vcfeval's
    annotated output keeps the call's FORMAT."""
    if not rtg.is_file():
        raise FileNotFoundError(f"no rtg at {rtg}: run scripts/bootstrap.sh")
    vafs: list[float] = []
    work = Path(tempfile.mkdtemp(prefix=".vcfeval-", dir=gvcf.parent))
    try:
        variants = work / "variants.vcf.gz"
        run_tool(
            gated(
                "dna:calibrate",
                ["bcftools", "view", "--no-version", "-f", "PASS", "-Oz"]
                + ["-o", str(variants), "--write-index=tbi", str(gvcf)],
            )
        )
        out = work / GENOTYPE
        _vcfeval(
            rtg, truth_vcf, variants, bench_bed, sdf, out, squash=False, threads=threads
        )
        for chrom, ref, alt, call, gt, vaf in _query(
            out / "calls.vcf.gz", "%CHROM\t%REF\t%ALT\t%INFO/CALL\t[%GT]\t[%VAF]\n"
        ):
            indices = gt.replace("|", "/").split("/")
            called = _gt_alleles(gt)
            if call != "TP" or ucsc_contig(chrom) not in _AUTOSOMES:
                continue
            if len(indices) != 2 or "0" not in indices or len(called) != 1:
                continue
            (k,) = called
            if len(ref) == 1 and len(alt.split(",")[k - 1]) == 1:
                vafs.append(float(vaf.split(",")[k - 1]))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return vafs


def _vcfeval(
    rtg: Path,
    truth: Path,
    calls: Path,
    bench: Path,
    sdf: Path,
    out: Path,
    *,
    squash: bool,
    threads: int,
) -> None:
    args = [
        str(rtg), "vcfeval", "-b", str(truth), "-c", str(calls), "-e", str(bench),
        "-t", str(sdf), "-o", str(out), "--output-mode", "annotate", "-T", str(threads),
    ]  # fmt: skip
    roc = ["--squash-ploidy", "--no-roc"] if squash else ["--roc-subset", "snp,non-snp"]
    run_tool(gated("dna:calibrate", args + roc))


def _query(path: Path, fmt: str) -> Iterator[list[str]]:
    """`bcftools query -f fmt` over `path`, one row of fields at a time (a gVCF has hundreds
    of millions of rows: streamed, never held). Raises `ToolFailed` when bcftools fails."""
    if shutil.which("bcftools") is None:
        raise ToolMissing("bcftools")
    args = ["bcftools", "query", "-f", fmt, str(path)]
    with (
        tempfile.TemporaryFile(mode="w+") as err,
        subprocess.Popen(
            args,
            stdout=subprocess.PIPE,
            stderr=err,
            stdin=subprocess.DEVNULL,
            text=True,
        ) as proc,
    ):
        assert proc.stdout is not None
        for line in proc.stdout:
            yield line.rstrip("\n").split("\t")
        proc.wait()
        if proc.returncode != 0:
            err.seek(0)
            raise ToolFailed(proc.returncode, args, None, err.read())


def _gt_alleles(gt: str) -> set[int]:
    """The ALT numbers a genotype calls."""
    return {int(i) for i in gt.replace("|", "/").split("/") if i not in (".", "0")}


def _allele_class(ref: str, alts: list[str]) -> str:
    return "indel" if any(alt == "*" or len(alt) != len(ref) for alt in alts) else "snv"


def _count_variants(
    calls: Path, kind: str, difficult: dict[str, _Intervals], by_chrom: dict[str, Cells]
) -> None:
    """Each scored PASS record of vcfeval's annotated calls into its variant cell."""
    for chrom, pos, ref, alt, call, gt, gq, dp in _query(
        calls, "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CALL\t[%GT]\t[%GQ]\t[%DP]\n"
    ):
        if call not in _SCORED or gq == ".":
            continue
        alts = alt.split(",")
        called = [alts[i - 1] for i in _gt_alleles(gt)]
        if not called:
            continue
        contig = ucsc_contig(chrom)
        region = "difficult" if _inside(difficult.get(contig), int(pos) - 1) else "easy"
        key = CellKey("variant", _allele_class(ref, called), region, gq_bin(int(gq)))
        cells = by_chrom.setdefault(contig, {})
        cells.setdefault((kind, key), Counts()).add(
            int(call in _ERRORS[kind]), 1, 0.0 if dp == "." else float(dp)
        )


def _affected(pos: int, ref: str, alt: str) -> tuple[int, str]:
    """A truth allele's 0-based position (see the module) and its class."""
    if len(ref) == len(alt):
        first = next(
            (i for i, (r, a) in enumerate(zip(ref, alt, strict=True)) if r != a), 0
        )
        return pos - 1 + first, "snv"
    shared = len(os.path.commonprefix([ref, alt]))
    if len(ref) > len(alt):
        return pos - 1 + shared, "indel"
    return pos - 1 + max(shared - 1, 0), "indel"


def _truth_positions(truth: Path) -> dict[str, tuple[list[int], list[str]]]:
    """Each called truth allele's position and class, per contig, sorted by position; two
    alleles of one class at one position (a `1|2` site) are one error of its unit."""
    found: dict[str, list[tuple[int, str]]] = {}
    for chrom, pos, ref, alt, gt in _query(truth, "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n"):
        alts = alt.split(",")
        for i in _gt_alleles(gt):
            found.setdefault(ucsc_contig(chrom), []).append(
                _affected(int(pos), ref.upper(), alts[i - 1].upper())
            )
    out: dict[str, tuple[list[int], list[str]]] = {}
    for contig, found_alleles in found.items():
        alleles = sorted(set(found_alleles))
        out[contig] = ([p for p, _ in alleles], [c for _, c in alleles])
    return out


def _count_reference(
    gvcf: Path,
    truth: dict[str, tuple[list[int], list[str]]],
    bench: dict[str, _Intervals],
    difficult: dict[str, _Intervals],
    by_chrom: dict[str, Cells],
) -> None:
    """Each RefCall's and reference block's positions inside the benchmark into their
    cells, with the truth alleles there as errors."""
    for chrom, pos, end, filt, alt, gt, gq, min_dp, dp in _query(
        gvcf, "%CHROM\t%POS\t%END\t%FILTER\t%ALT\t[%GT]\t[%GQ]\t[%MIN_DP]\t[%DP]\n"
    ):
        if alt == _SYMBOLIC:
            row_type, depth = "ref_block", min_dp
        elif filt == "RefCall":
            row_type, depth = "RefCall", dp
        else:
            continue
        if "." in gt or gq == "." or depth in (".", "0"):
            continue  # a no-call, or no reads: not a call
        contig = ucsc_contig(chrom)
        within = bench.get(contig)
        if within is None:
            continue
        quality = min(int(gq), BLOCK_GQ_CAP) if row_type == "ref_block" else int(gq)
        bin_ = gq_bin(quality)
        positions, classes = truth.get(contig, ([], []))
        cells = by_chrom.setdefault(contig, {})
        for a, b in _pieces(within, int(pos) - 1, int(end)):
            for start, stop, region in _split(difficult.get(contig), a, b):
                first, last = (
                    bisect_left(positions, start),
                    bisect_left(positions, stop),
                )
                hits = Counter(classes[first:last]) if first < last else Counter()
                for allele_class in ALLELE_CLASSES:
                    key = _cell_key(row_type, allele_class, region, bin_)
                    for kind in KINDS:
                        cells.setdefault((kind, key), Counts()).add(
                            hits[allele_class], stop - start, float(depth)
                        )


@functools.cache
def _cell_key(row_type: str, allele_class: str, region: str, gq: int) -> CellKey:
    """One `CellKey` per cell, made once: a gVCF has hundreds of millions of rows."""
    return CellKey(row_type, allele_class, region, gq)


# --- intervals ---------------------------------------------------------------------------


def _merged(raw: dict[str, list[tuple[int, int]]]) -> dict[str, _Intervals]:
    out: dict[str, _Intervals] = {}
    for contig, intervals in raw.items():
        starts: list[int] = []
        ends: list[int] = []
        for start, end in sorted(intervals):
            if starts and start <= ends[-1]:
                ends[-1] = max(ends[-1], end)
            else:
                starts.append(start)
                ends.append(end)
        out[contig] = (starts, ends)
    return out


def _read_bed(path: Path) -> dict[str, _Intervals]:
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            raw.setdefault(ucsc_contig(chrom), []).append((int(start), int(end)))
    return _merged(raw)


def _union(strata: StrataIndex) -> dict[str, _Intervals]:
    """Every stratum's intervals together: where any flag applies."""
    raw: dict[str, list[tuple[int, int]]] = {}
    for contigs in strata.intervals.values():
        for contig, (starts, ends) in contigs.items():
            raw.setdefault(contig, []).extend(zip(starts, ends, strict=True))
    return _merged(raw)


def _inside(intervals: _Intervals | None, offset: int) -> bool:
    if intervals is None:
        return False
    starts, ends = intervals
    i = bisect_right(starts, offset) - 1
    return i >= 0 and offset < ends[i]


def _pieces(intervals: _Intervals, start: int, end: int) -> Iterator[tuple[int, int]]:
    """The parts of [start, end) inside `intervals`."""
    starts, ends = intervals
    i = bisect_right(ends, start)
    while i < len(starts) and starts[i] < end:
        a, b = max(start, starts[i]), min(end, ends[i])
        if a < b:
            yield a, b
        i += 1


def _split(
    difficult: _Intervals | None, start: int, end: int
) -> Iterator[tuple[int, int, str]]:
    """[start, end) cut into its difficult and easy parts."""
    cursor = start
    for a, b in _pieces(difficult, start, end) if difficult is not None else ():
        if cursor < a:
            yield cursor, a, "easy"
        yield a, b, "difficult"
        cursor = b
    if cursor < end:
        yield cursor, end, "easy"
