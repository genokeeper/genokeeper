"""HG003 through the run's own pipeline: its FASTQ downsampled by bases to the user's, then
aligned and called by the steps `genokeeper run` aligns and calls with.

`downsample` keeps a fraction of the read pairs with seqtk, one seed for both mates: seqtk
draws the same sequence of random numbers for each file, so the reads it keeps in one file
are the mates of those it keeps in the other (the files list the pairs in one order). The
fraction is the target over the bases of both files, as `seqtk size` counts them.

`run_sample` aligns, gates and calls through `pipeline.align_step`, `pipeline.qc_step`,
`pipeline.call_step` and `pipeline.call_qc_step` with `pipeline.qc_gate`, the functions the
run itself uses: the calibration measures the pipeline the user's calls come from, on a
sample that passes the user's QC gate, with the same bwa, samtools, image, checkpoint and
DeepVariant arguments for the sample's sex (the shard count follows this machine's CPUs and Docker's memory, as a run's
does; it is not one of the calibration's versions). Its work lives under
`<lab>/calibration/<name>/`, cached as a run's steps are, behind its own run lock.
"""

import json
import os
import secrets
import shlex
import subprocess
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from genokeeper.align import PlatformRefused, check_platform
from genokeeper.calibrate import calibration_dir
from genokeeper.checksum import sha256_cached
from genokeeper.config import Grading, Paths, load_grading
from genokeeper.containers import DockerUnavailable
from genokeeper.pipeline import (
    IDENTITIES,
    RunRefused,
    _Step,
    align_step,
    call_qc_step,
    call_step,
    calling_inputs,
    genome_tool_versions,
    qc_gate,
    qc_step,
    run_lock,
    sex_step,
)
from genokeeper.regions import QC_AUTOSOMES, SEX_REGIONS, giab38_dir
from genokeeper.sex import Sex
from genokeeper.tools import ToolFailed, gated, require, run_tool

SEED = 11  # seqtk's own default seed


def sample_fraction(total_bases: int, target_bases: int) -> float:
    """The fraction of reads that keeps `target_bases` of `total_bases`; a `ValueError` for
    a target that is not below the total (downsampling cannot add bases) or not positive."""
    if total_bases <= 0 or not 0 < target_bases < total_bases:
        raise ValueError(
            f"cannot downsample {total_bases} bases to {target_bases}: the target must be"
            " positive and below the total"
        )
    return target_bases / total_bases


def fastq_bases(fastq: Path) -> int:
    """The bases in `fastq`, as `seqtk size` counts them (it prints reads, then bases)."""
    out = run_tool(gated("dna:downsample", ["seqtk", "size", str(fastq)]))
    return int(out.split()[1])


def downsample(
    fastq1: Path, fastq2: Path, out_dir: Path, *, target_bases: int, seed: int = SEED
) -> tuple[Path, Path]:
    """Keep `target_bases` of the pair's bases: the same fraction of each file sampled by
    `seqtk sample` with `seed`, bgzipped as `out_dir/<each file's name>`. Raises
    `ValueError` for a target the pair cannot reach (before writing anything), `ToolMissing`
    without seqtk or bgzip, and `ToolFailed` when a stage fails; both files land only when
    both were written."""
    require("seqtk", "bgzip")
    outs = (out_dir / fastq1.name, out_dir / fastq2.name)
    if fastq1.name == fastq2.name or any(
        o.resolve() == f.resolve() for o, f in zip(outs, (fastq1, fastq2), strict=True)
    ):
        raise ValueError(f"downsample {fastq1} and {fastq2} into another folder")
    fraction = sample_fraction(fastq_bases(fastq1) + fastq_bases(fastq2), target_bases)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmps = [o.with_name(f".{o.name}.{secrets.token_hex(4)}") for o in outs]
    try:
        for fastq, tmp in zip((fastq1, fastq2), tmps, strict=True):
            sample = ["seqtk", "sample", "-s", str(seed), str(fastq), repr(fraction)]
            script = (
                f"set -o pipefail; {shlex.join(sample)} | bgzip -c"
                f" > {shlex.quote(str(tmp))}"
            )
            # One gate per pipeline: seqtk and bgzip hold one slot together.
            command = gated("dna:downsample", ["bash", "-c", script])
            done = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                errors="replace",
                stdin=subprocess.DEVNULL,
            )
            if done.returncode != 0:
                raise ToolFailed(done.returncode, command, None, done.stderr)
        for tmp, out in zip(tmps, outs, strict=True):
            os.replace(tmp, out)
    finally:
        for tmp in tmps:
            tmp.unlink(missing_ok=True)
    return outs


@dataclass(frozen=True)
class _Aligned:
    """What the align, QC and sex steps leave for calling: the alignment, its metrics, the
    inputs' hashes and the settings the steps ran with."""

    step: _Step
    metrics: dict[str, float]
    fasta: Path
    fasta_sha: str
    checkpoint: Path
    identities: Path
    bcftools: str
    grading: Grading
    threads: int
    work: Path


@contextmanager
def _aligned(
    paths: Paths,
    *,
    name: str,
    fastq1: Path,
    fastq2: Path,
    sex: Sex,
    grading: Grading | None,
    threads: int | None,
    from_cram: bool,
) -> Iterator[_Aligned]:
    """The run's align, QC gate and sex check of the pair, under `<lab>/calibration/<name>/`
    and its run lock, which stays held for the block (see `align_sample`). A pair written
    back from a calibration CRAM (`from_cram`) skips the platform gate's read-name test."""
    threads = threads if threads is not None else (os.cpu_count() or 4)
    versions = genome_tool_versions()
    for what, path in (("FASTQ 1", fastq1), ("FASTQ 2", fastq2)):
        if not path.is_file():
            raise RunRefused(f"missing input {what}: {path}")
    fasta, checkpoint = calling_inputs(paths.lab)
    for regions in (QC_AUTOSOMES, SEX_REGIONS):
        if not (giab38_dir(paths.lab) / regions).is_file():
            raise RunRefused(
                f"missing input: the QC regions {regions}; run scripts/bootstrap.sh"
            )
    grading = grading if grading is not None else load_grading()
    work = calibration_dir(paths.lab) / name
    with run_lock(work):
        try:
            check_platform(fastq1, names=not from_cram)
        except PlatformRefused as e:
            raise RunRefused(" ".join(str(e).split())) from e
        identities = work / IDENTITIES
        shas = {
            "fastq1": sha256_cached(fastq1, identities),
            "fastq2": sha256_cached(fastq2, identities),
            "fasta": sha256_cached(fasta, identities),
        }
        aligned = align_step(
            work,
            fastq1=fastq1,
            fastq2=fastq2,
            fasta=fasta,
            shas=shas,
            versions=versions,
            threads=threads,
        )
        measured = qc_step(
            work,
            aligned=aligned,
            fasta=fasta,
            lab=paths.lab,
            versions=versions,
            grading=grading,
            threads=threads,
        )
        metrics: dict[str, float] = json.loads(measured.path.read_text())
        qc_gate(metrics, grading)
        inferred = sex_step(
            work,
            aligned=aligned,
            fasta=fasta,
            lab=paths.lab,
            versions=versions,
            grading=grading,
            threads=threads,
        )
        found = json.loads(inferred.path.read_text())["sex"]
        if found != sex:
            raise RunRefused(f"{name}'s inferred sex is {found}, not {sex}")
        yield _Aligned(
            aligned,
            metrics,
            fasta,
            shas["fasta"],
            checkpoint,
            identities,
            versions["bcftools"],
            grading,
            threads,
            work,
        )


def align_sample(
    paths: Paths,
    *,
    name: str,
    fastq1: Path,
    fastq2: Path,
    sex: Sex,
    grading: Grading | None = None,
    threads: int | None = None,
    from_cram: bool = False,
) -> tuple[Path, str, dict[str, float]]:
    """Align the pair as `genokeeper run` does and gate it, under
    `<lab>/calibration/<name>/`; return the CRAM, its align step's key and the QC metrics.
    The sample passes the run's QC gate on the alignment, and the run's sex step must infer
    `sex` from it. Each step is cached by its inputs, so a rerun reuses what is done.
    `grading` defaults to config/grading.toml, `threads` to this machine's CPUs.
    `from_cram` says the pair was written back from a calibration CRAM
    (`genokeeper.calibrate.realign`): its reads carry ENA's accession, not the flowcell, and
    the platform gate tests their qualities only.

    Raises `RunRefused` on a missing tool or input (naming the bootstrap), FASTQ from
    another platform, a failed QC gate, another inferred sex or another calibration run
    holding the lock."""
    with _aligned(
        paths,
        name=name,
        fastq1=fastq1,
        fastq2=fastq2,
        sex=sex,
        grading=grading,
        threads=threads,
        from_cram=from_cram,
    ) as a:
        return a.step.path, a.step.key, a.metrics


@dataclass(frozen=True)
class SampleRun:
    """A called calibration sample: the CRAM, its align step's key, the QC metrics with the
    call's Ti/Tv and het/hom, and the gVCF (the VCF lies beside it)."""

    cram: Path
    key: str
    metrics: dict[str, float]
    gvcf: Path


def run_sample(
    paths: Paths,
    *,
    name: str,
    fastq1: Path,
    fastq2: Path,
    sex: Sex,
    grading: Grading | None = None,
    threads: int | None = None,
    from_cram: bool = False,
) -> SampleRun:
    """`align_sample`, then call the pair as `genokeeper run` does; return the `SampleRun`.
    The calls pass the run's QC gate too, against the alignment's metrics, or the sample is
    refused as a run's sample is.

    Raises what `align_sample` raises, and `RunRefused` when Docker is not running."""
    with _aligned(
        paths,
        name=name,
        fastq1=fastq1,
        fastq2=fastq2,
        sex=sex,
        grading=grading,
        threads=threads,
        from_cram=from_cram,
    ) as a:
        try:
            called = call_step(
                a.work,
                aligned=a.step,
                fasta=a.fasta,
                fasta_sha=a.fasta_sha,
                checkpoint=a.checkpoint,
                identities=a.identities,
                sex=sex,
                grading=a.grading,
                threads=a.threads,
            )
        except DockerUnavailable as e:
            raise RunRefused(" ".join(str(e).split())) from e
        stats = call_qc_step(a.work, called=called, bcftools=a.bcftools)
        metrics = a.metrics | json.loads(stats.path.read_text())
        qc_gate(metrics, a.grading)
    return SampleRun(
        cram=a.step.path, key=a.step.key, metrics=metrics, gvcf=called.path
    )


def recall_sample(
    paths: Paths,
    *,
    name: str,
    cram: Path,
    aligned_key: str,
    sex: Sex,
    threads: int | None = None,
) -> Path:
    """Call a recorded alignment again with this pipeline's caller, as `run_sample` calls a
    pair, under `<lab>/calibration/<name>/` and its run lock; the calls pass the run's QC gate
    on Ti/Tv and het/hom. Returns the new gVCF. Nothing is aligned again (BR-043)."""
    threads = threads if threads is not None else (os.cpu_count() or 4)
    versions = genome_tool_versions()
    fasta, checkpoint = calling_inputs(paths.lab)
    grading = load_grading()
    work = calibration_dir(paths.lab) / name
    with run_lock(work):
        identities = work / IDENTITIES
        try:
            called = call_step(
                work,
                aligned=_Step(key=aligned_key, path=cram, reused=True),
                fasta=fasta,
                fasta_sha=sha256_cached(fasta, identities),
                checkpoint=checkpoint,
                identities=identities,
                sex=sex,
                grading=grading,
                threads=threads,
            )
        except DockerUnavailable as e:
            raise RunRefused(" ".join(str(e).split())) from e
        stats = call_qc_step(work, called=called, bcftools=versions["bcftools"])
        qc_gate(json.loads(stats.path.read_text()), grading)
    return called.path
