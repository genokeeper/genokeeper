"""The stage 2a QC gate: measures the aligned sample and refuses it before calling when a
metric falls outside config/grading.toml [qc].

- mosdepth: mean depth and the share of bases at 20x or more over the non-gap autosomes
  (`genokeeper.regions.QC_AUTOSOMES`), reads at MAPQ >= [qc] depth_mapq, duplicates excluded;
- samtools stats: mapped and duplicate fractions, the median insert size;
- VerifyBamID2 (container, offline): FREEMIX, the contamination estimate;
- with a VCF, bcftools stats over PASS SNVs: Ti/Tv and the het/hom ratio.

The gate's anchor comes from the thresholds in config/grading.toml, each with its published
source, which the sample under test cannot move because every metric is measured here from
the CRAM and VCF rather than read from a delivery report.
"""

import gzip
import shutil
import tempfile
from pathlib import Path

from genokeeper.config import Grading, QcThresholds, load_paths
from genokeeper.containers import PLATFORM, run_container
from genokeeper.regions import QC_AUTOSOMES, giab38_dir
from genokeeper.tools import gated, run_tool

# Where the BioContainers image keeps VerifyBamID2's 1000 Genomes SVD resources (GRCh38,
# 100k markers, chr names).
SVD_PREFIX = (
    "/usr/local/share/verifybamid2-2.0.3-0/resource/1000g.phase3.100k.b38.vcf.gz.dat"
)


def mosdepth(
    cram: Path,
    fasta: Path,
    regions: Path,
    prefix: Path,
    *,
    mapq: int,
    threads: int,
    thresholds: int | None = None,
) -> dict[str, tuple[float, float]]:
    """Run mosdepth over the named BED `regions` (behind the gate as `dna:qc`) and return,
    per region name, (mean depth, share of bases at `thresholds`x or more; 0 without it),
    each weighted by interval length."""
    args = ["mosdepth", "-t", str(threads), "-n", "-Q", str(mapq), "-f", str(fasta)]
    args += ["--by", str(regions)]
    if thresholds is not None:
        args += ["--thresholds", str(thresholds)]
    run_tool(gated("dna:qc", [*args, str(prefix), str(cram)]))
    length: dict[str, int] = {}
    bases: dict[str, float] = {}
    covered: dict[str, int] = {}
    with gzip.open(f"{prefix}.regions.bed.gz", "rt") as f:
        for line in f:
            _, start, end, name, mean = line.rstrip("\n").split("\t")
            n = int(end) - int(start)
            length[name] = length.get(name, 0) + n
            bases[name] = bases.get(name, 0.0) + float(mean) * n
    if thresholds is not None:
        with gzip.open(f"{prefix}.thresholds.bed.gz", "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                _, _, _, name, count = line.rstrip("\n").split("\t")
                covered[name] = covered.get(name, 0) + int(count)
    return {
        name: (bases[name] / n, covered.get(name, 0) / n)
        for name, n in length.items()
        if n
    }


def _samtools_stats(cram: Path, fasta: Path, threads: int) -> dict[str, float]:
    out = run_tool(
        gated(
            "dna:qc",
            [
                "samtools",
                "stats",
                "-@",
                str(threads),
                "--reference",
                str(fasta),
                str(cram),
            ],
        )
    )
    summary: dict[str, float] = {}
    inserts: list[tuple[int, int]] = []
    for line in out.splitlines():
        fields = line.split("\t")
        if fields[0] == "SN":
            summary[fields[1].rstrip(":")] = float(fields[2])
        elif fields[0] == "IS":
            inserts.append((int(fields[1]), int(fields[2])))
    total = summary["raw total sequences"]
    mapped = summary["reads mapped"]
    return {
        "mapped_fraction": mapped / total if total else 0.0,
        "duplicate_fraction": summary["reads duplicated"] / mapped if mapped else 0.0,
        "median_insert": _median(inserts),
        # Every primary read's bases: the FASTQ's, which `genokeeper calibrate` matches.
        "total_bases": summary["total length"],
    }


def _median(histogram: list[tuple[int, int]]) -> float:
    """The median of a (value, count) histogram; 0 when it is empty."""
    total = sum(count for _, count in histogram)
    seen = 0
    for value, count in sorted(histogram):
        seen += count
        if count and seen * 2 >= total:
            return float(value)
    return 0.0


def _freemix(cram: Path, fasta: Path, work: Path, threads: int) -> float:
    run_container(
        "verifybamid2",
        [
            "verifybamid2",
            "--SVDPrefix",
            SVD_PREFIX,
            "--Reference",
            f"/ref/{fasta.name}",
            "--BamFile",
            f"/reads/{cram.name}",
            "--Output",
            "/out/vb",
            "--NumThread",
            str(threads),
        ],
        mounts=[
            (fasta.parent, "/ref", True),
            (cram.parent, "/reads", True),
            (work, "/out", False),
        ],
        platform=PLATFORM,
        label="dna:qc",
    )
    header, values = (work / "vb.selfSM").read_text().splitlines()[:2]
    return float(
        dict(zip(header.split("\t"), values.split("\t"), strict=True))["FREEMIX"]
    )


def vcf_stats(vcf: Path) -> dict[str, float]:
    """Ti/Tv and het/hom over the PASS SNVs of a one-sample VCF (bcftools stats)."""
    out = run_tool(
        ["bcftools", "stats", "-f", "PASS", "-i", 'TYPE="snp"', "-s", "-", str(vcf)]
    )
    metrics: dict[str, float] = {}
    for line in out.splitlines():
        fields = line.split("\t")
        if fields[0] == "TSTV":
            metrics["titv"] = float(fields[4])
        elif fields[0] == "PSC":
            hom, het = int(fields[4]), int(fields[5])
            metrics["het_hom"] = het / hom if hom else float("inf")
    return metrics


def judge(metrics: dict[str, float], qc: QcThresholds) -> list[str]:
    """One line per metric outside its threshold, in the order the gate measures them."""
    lower = {
        "mean_depth": qc.min_mean_depth,
        "fraction_20x": qc.min_fraction_20x,
        "mapped_fraction": qc.min_mapped_fraction,
        "median_insert": qc.min_median_insert,
        "titv": qc.min_titv,
    }
    upper = {
        "duplicate_fraction": qc.max_duplicate_fraction,
        "freemix": qc.max_freemix,
        "titv": qc.max_titv,
        "het_hom": qc.max_het_hom,
    }
    failed: list[str] = []
    for name, value in metrics.items():
        if name in lower and value < lower[name]:
            failed.append(f"{name} {value:.4g} < {lower[name]:.4g}")
        elif name in upper and value > upper[name]:
            failed.append(f"{name} {value:.4g} > {upper[name]:.4g}")
    return failed


def run_qc(
    cram: Path,
    fasta: Path,
    *,
    grading: Grading,
    lab: Path | None = None,
    threads: int = 4,
) -> dict[str, float]:
    """The metrics of `cram` (indexed, written against `fasta`); the derived QC region BED
    comes from `lab` (default: the configured lab). Raises what the tools raise; `judge`
    names the metrics outside [qc]."""
    lab = lab if lab is not None else load_paths().lab
    regions = giab38_dir(lab) / QC_AUTOSOMES
    if not regions.is_file():
        raise FileNotFoundError(f"{regions}: run scripts/bootstrap.sh")
    work = Path(tempfile.mkdtemp(prefix=".qc-", dir=cram.parent))
    try:
        depth = mosdepth(
            cram,
            fasta,
            regions,
            work / "depth",
            mapq=grading.qc.depth_mapq,
            threads=threads,
            thresholds=20,
        )
        mean, at_20x = depth["autosomes"]
        metrics = {"mean_depth": mean, "fraction_20x": at_20x}
        metrics |= _samtools_stats(cram, fasta, threads)
        metrics["freemix"] = _freemix(cram, fasta, work, threads)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return metrics
