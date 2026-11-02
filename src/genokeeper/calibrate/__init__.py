"""Stage 2b, calibration: genotype error rates measured on GIAB's HG003 (NA24149) sequenced
on the original DNBSEQ-T7, through exactly the run's pipeline, in place of DeepVariant's own
estimate (track design, stage 2 "2b. Calibration").

`scripts/bootstrap.sh --with-calibration` (step_hg003) fetches HG003's FASTQ and its GIAB
v4.2.1 GRCh38 benchmark, `--with-masked-reference` (step_masked_reference) builds rtg's SDF of
the lab's GRCh38 FASTA, and `--with-hg004` (step_hg004) fetches HG004's benchmark, into the
places `layout` names; `genokeeper.calibrate.sample` downsamples the FASTQ and runs it through the
run's align, QC and call steps, and `genokeeper.calibrate.measure` scores the gVCF against the
benchmark into a `genokeeper.calibrate.table.CalibrationTable`. `genokeeper calibrate`
(`genokeeper.calibrate.tables`, its `small_variants` row) does both and writes the table to
`table_dir`, which every `genokeeper run` requires.

    python -m genokeeper.calibrate layout     # key<TAB>path, one per line
"""

import argparse
from pathlib import Path

from genokeeper.config import grch38_fasta, load_paths

# GIAB's HG003 v4.2.1 GRCh38 benchmark (chr1-22), as the bootstrap pins and names them.
TRUTH_VCF = "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BENCH_BED = "HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
# GIAB's HG004 v4.2.1 GRCh38 benchmark (chr1-22), which `hg002 accept hg004` scores against.
HG004_TRUTH_VCF = "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_BENCH_BED = "HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
# The rtg-tools release the bootstrap installs (step_rtg).
RTG_VERSION = "3.13"


def calibration_dir(lab: Path) -> Path:
    """The lab's calibration samples, their runs and the tables (deny-ruled for reads)."""
    return lab / "calibration"


def table_dir(lab: Path) -> Path:
    """The calibration table `genokeeper calibrate` writes and `genokeeper run` requires
    (`tables.folder(lab, "small_variants")`)."""
    from genokeeper.calibrate.tables import folder  # tables imports this package

    return folder(lab, "small_variants")


def hg003_dir(lab: Path) -> Path:
    """HG003's FASTQ, downsampled FASTQ and run, under the calibration folder."""
    return calibration_dir(lab) / "hg003"


def hg004_truth(lab: Path) -> tuple[Path, Path]:
    """GIAB's HG004 v4.2.1 truth VCF (its `.tbi` beside it) and benchmark BED in the lab."""
    folder = lab / "giab" / "HG004-v4.2.1-GRCh38"
    return folder / HG004_TRUTH_VCF, folder / HG004_BENCH_BED


def layout(lab: Path) -> dict[str, Path]:
    """Where the bootstrap puts the calibration inputs: `fastq` the full HG003 FASTQ,
    `downsampled` the pair `downsample` writes (once there, `genokeeper calibrate` has
    deleted the full FASTQ, and it is not fetched again), `truth` GIAB's HG003 benchmark,
    `hg004_truth` HG004's, `sdf` rtg's SDF of `fasta`, the lab's GRCh38 FASTA."""
    fasta = grch38_fasta(lab)
    return {
        "fastq": hg003_dir(lab) / "fastq",
        "downsampled": hg003_dir(lab) / "downsampled",
        "truth": lab / "giab" / "HG003-v4.2.1-GRCh38",
        "hg004_truth": hg004_truth(lab)[0].parent,
        "fasta": fasta,
        "sdf": fasta.with_suffix(".sdf"),
    }


def rtg_path(lab: Path) -> Path:
    """The `rtg` launcher the bootstrap installs (with its rtg.cfg turning off talkback)."""
    return lab / "tools" / f"rtg-tools-{RTG_VERSION}" / "rtg"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.calibrate")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("layout")
    parser.parse_args(argv)
    for key, path in layout(load_paths().lab).items():
        print(f"{key}\t{path}")
    return 0
