"""DeepVariant 1.10.0 with the Complete Genomics T7 checkpoint: a VCF and a gVCF.

The checkpoint is `complete-t7/weights-51-0.995354.ckpt` from Google's
complete-case-study-testdata bucket (the bootstrap's step_dv_model pins its four files'
md5s). Its `example_info.json` reads
`{"version": "1.6.1", "shape": [100, 221, 7], "channels": [1, 2, 3, 4, 5, 6, 19]}`: model
version 1.6.1, a 100 x 221 pileup with 7 channels, the six defaults (read base, base quality,
mapping quality, strand, read supports variant, base differs from reference) and channel 19,
insert size. `model.example_info.json` adds the pileup height and width as calling flags.

The T7 case study (docs/deepvariant-complete-t7-case-study.md, branch r1.10) runs
`run_deepvariant --model_type WGS --customized_model .../weights-51-0.995354.ckpt` and passes
no make_examples flags: run_deepvariant takes the channels from the example_info.json beside
the checkpoint, so the checkpoint's directory is mounted whole. The only make_examples
argument added here is `gvcf_gq_binsize=1`, which keeps every GQ value in the reference
blocks (spec, stage 2a "Calling"). A male sample gets chrX and chrY as haploid contigs outside
the GRCh38 PARs of config/grading.toml [par].

The container runs offline from the pinned image (genokeeper.containers) behind the core gate
as `dna:call`. DeepVariant writes into a work directory; the VCF, the gVCF and their indexes
reach the output directory only after the container exits 0 with all four written, so a
killed or failed run promotes nothing.
"""

import os
import secrets
import shutil
from pathlib import Path

from genokeeper.config import Grading
from genokeeper.containers import PLATFORM, docker_memory, run_container

RUN_DEEPVARIANT = "/opt/deepvariant/bin/run_deepvariant"
# Memory allowed per make_examples shard (`--num_shards`), so the shards fit Docker's VM.
# Measured on HG003 downsampled to a 30x run's bases, under Rosetta (2026-10-22, 15 to
# 20 minutes into make_examples): 1.62 GiB peak RSS per shard, 15 GiB for 11 shards; 2 GiB
# is that peak with a quarter's margin. The earlier 4 GiB (Google's n2-standard-96 run, 384
# GiB over 96 shards) left a quarter of the cores idle where Docker's memory set the
# shard count.
SHARD_MEMORY = 2 * 2**30
CHECKPOINT = "weights-51-0.995354.ckpt"
# The files run_deepvariant reads beside `--customized_model`.
CHECKPOINT_FILES = (
    "example_info.json",
    "model.example_info.json",
    f"{CHECKPOINT}.data-00000-of-00001",
    f"{CHECKPOINT}.index",
)
PAR_BED = "par38.bed"
# The arguments that decide what DeepVariant calls on the autosomes on any machine: the model
# type and make_examples' own extra arguments. The calibration table records them (its
# `make_examples_flags`). The shard count (this machine's CPUs and Docker's memory) is not
# among them, nor is a male's haploid chrX and chrY, which HG003's benchmark does not cover.
CALLING_FLAGS = ("--model_type=WGS", "--make_examples_extra_args=gvcf_gq_binsize=1")


def checkpoint_path(lab: Path) -> Path:
    """The T7 checkpoint prefix in the lab, as the bootstrap fetches it."""
    return lab / "deepvariant" / "complete-t7" / CHECKPOINT


def shard_count(threads: int, memory: int) -> int:
    """The shards `memory` bytes hold at `SHARD_MEMORY` each, at most `threads`, at least
    one."""
    return max(1, min(threads, memory // SHARD_MEMORY))


def _outputs(cram: Path) -> tuple[str, str]:
    return f"{cram.stem}.vcf.gz", f"{cram.stem}.g.vcf.gz"


def deepvariant_args(
    cram: Path,
    fasta: Path,
    out_dir: Path,
    *,
    checkpoint: Path,
    male: bool,
    par_bed: Path,
    threads: int,
) -> list[str]:
    """The container's command line, in container paths (see `_mounts`): `out_dir` is /out,
    and for a male `par_bed` is read from /par."""
    del out_dir  # mounted at /out
    vcf, gvcf = _outputs(cram)
    model_type, extra = CALLING_FLAGS
    args = [
        RUN_DEEPVARIANT,
        model_type,
        f"--customized_model=/model/{checkpoint.name}",
        f"--ref=/ref/{fasta.name}",
        f"--reads=/reads/{cram.name}",
        f"--output_vcf=/out/{vcf}",
        f"--output_gvcf=/out/{gvcf}",
        f"--num_shards={threads}",
        extra,
        "--intermediate_results_dir=/out/intermediate",
    ]
    if male:
        args += [
            "--haploid_contigs=chrX,chrY",
            f"--par_regions_bed=/par/{par_bed.name}",
        ]
    return args


def _mounts(
    cram: Path, fasta: Path, work: Path, checkpoint: Path, par_bed: Path | None
) -> list[tuple[Path, str, bool]]:
    mounts = [
        (fasta.parent, "/ref", True),
        (cram.parent, "/reads", True),
        (checkpoint.parent, "/model", True),
        (work, "/out", False),
    ]
    if par_bed is not None:
        mounts.append((par_bed.parent, "/par", True))
    return mounts


def _write_par_bed(path: Path, grading: Grading) -> Path:
    """[par] (1-based, inclusive) as a BED (0-based, half-open)."""
    rows = [
        f"{contig}\t{start - 1}\t{end}\n"
        for contig, intervals in sorted(grading.par.items())
        for start, end in intervals
    ]
    if not rows:
        raise ValueError("config/grading.toml has no [par] intervals")
    path.write_text("".join(rows))
    return path


def run_deepvariant(
    cram: Path,
    fasta: Path,
    out_dir: Path,
    *,
    checkpoint: Path,
    male: bool,
    grading: Grading,
    threads: int,
) -> tuple[Path, Path]:
    """Call `cram` (indexed, on `fasta`) into `out_dir`; return the VCF and the gVCF. Runs
    `threads` shards, fewer when Docker's VM lacks `SHARD_MEMORY` for each (`shard_count`).

    Raises `FileNotFoundError` when a checkpoint file is missing, `DockerUnavailable` when
    Docker is not running, `ToolFailed` when the container fails, and `RuntimeError` when it
    exits 0 without an output; in every case nothing is written to `out_dir`."""
    absent = [n for n in CHECKPOINT_FILES if not (checkpoint.parent / n).is_file()]
    if absent:
        raise FileNotFoundError(
            f"T7 checkpoint files missing in {checkpoint.parent}: {', '.join(absent)};"
            " run scripts/bootstrap.sh"
        )
    shards = shard_count(threads, docker_memory())
    work = out_dir / f".deepvariant-{secrets.token_hex(4)}"
    work.mkdir()
    try:
        par_bed = _write_par_bed(work / PAR_BED, grading) if male else None
        run_container(
            "deepvariant",
            deepvariant_args(
                cram,
                fasta,
                work,
                checkpoint=checkpoint,
                male=male,
                par_bed=par_bed if par_bed is not None else work / PAR_BED,
                threads=shards,
            ),
            mounts=_mounts(cram, fasta, work, checkpoint, par_bed),
            platform=PLATFORM,
            label="dna:call",
        )
        names = [n for output in _outputs(cram) for n in (output, f"{output}.tbi")]
        missing = [n for n in names if not (work / n).is_file()]
        if missing:
            raise RuntimeError(f"DeepVariant exited 0 without {', '.join(missing)}")
        # Indexes first: a VCF in out_dir always has its index beside it.
        for name in sorted(names, key=lambda n: not n.endswith(".tbi")):
            os.replace(work / name, out_dir / name)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    vcf, gvcf = _outputs(cram)
    return out_dir / vcf, out_dir / gvcf
