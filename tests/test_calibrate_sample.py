"""HG003's inputs through the run's own pipeline: the FASTQ downsampled by bases with one
seed for both mates, then aligned, QC-gated and called exactly as `genokeeper run` aligns,
gates and calls."""

import gzip
import stat
from dataclasses import replace
from pathlib import Path

import pytest
from fixtures.genome_mini import (
    MINI_QC,
    FakeDeepVariant,
    FakeVerifyBamID,
    GenomeLab,
    bases,
    variant,
    write_lab,
)
from fixtures.grch38_mini import CONTIGS as GRCH38_CONTIGS
from fixtures.grch38_mini import READ, write_fastq_pair, write_reference

from genokeeper import pipeline, qc, tools
from genokeeper.calibrate import calibration_dir
from genokeeper.calibrate.sample import (
    align_sample,
    downsample,
    run_sample,
    sample_fraction,
)
from genokeeper.call import deepvariant
from genokeeper.call.deepvariant import (
    PAR_BED,
    checkpoint_path,
    deepvariant_args,
    shard_count,
)
from genokeeper.config import QcThresholds, grch38_fasta, load_grading
from genokeeper.pipeline import RunRefused
from genokeeper.sex import Sex

GIB = 2**30
GRADING = replace(load_grading(), qc=MINI_QC)


def _logging_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A core gate that logs each command it is handed, then runs it."""
    log = tmp_path / "gate.log"
    gate = tmp_path / "coreslot-gate"
    gate.write_text(f'#!/usr/bin/env bash\necho "$*" >> "{log}"\nshift 2\nexec "$@"\n')
    gate.chmod(gate.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setattr(tools, "GATE", gate)
    return log


def _faked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[GenomeLab, FakeDeepVariant, FakeVerifyBamID]:
    """The mini lab with DeepVariant (a het and a hom SNV, so the calls' QC has a het/hom
    ratio) and VerifyBamID2 faked."""
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")
    lab = write_lab(tmp_path)
    records = []
    for pos, gt, ad in ((1000, "0/1", [15, 15]), (2000, "1/1", [0, 30])):
        alt = "A" if bases(lab.fasta, "chr1", pos) != "A" else "G"
        records.append(variant(lab.fasta, "chr1", pos, [alt], gt=gt, gq=50, ad=ad))
    fake = FakeDeepVariant(records, dict(GRCH38_CONTIGS))
    verifybamid = FakeVerifyBamID()
    monkeypatch.setattr(deepvariant, "run_container", fake)
    monkeypatch.setattr(deepvariant, "docker_memory", lambda: 64 * GIB)
    monkeypatch.setattr(qc, "run_container", verifybamid)
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: "male")
    return lab, fake, verifybamid


def _names(fastq: Path, mate: int) -> list[str]:
    with gzip.open(fastq, "rt") as f:
        lines = f.read().splitlines()
    names = lines[0::4]
    assert all(n.endswith(f"/{mate}") for n in names)
    return [n.removesuffix(f"/{mate}") for n in names]


def test_sample_fraction_refuses_what_it_cannot_reach() -> None:
    assert sample_fraction(1000, 250) == 0.25
    for total, target in ((1000, 1000), (1000, 2000), (1000, 0), (0, 10)):
        with pytest.raises(ValueError, match="downsample"):
            sample_fraction(total, target)


@pytest.mark.tools
def test_downsample_fraction_from_bases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    log = _logging_gate(tmp_path, monkeypatch)
    fasta = write_reference(tmp_path / "ref.fa")
    (tmp_path / "in").mkdir()
    pairs = 400
    fq1, fq2 = write_fastq_pair(tmp_path / "in", fasta, pairs=pairs)
    total = 2 * pairs * READ
    out = tmp_path / "out"
    got = downsample(fq1, fq2, out, target_bases=total // 4)
    assert got == (out / fq1.name, out / fq2.name)
    # One seed for both mates keeps every sampled read's mate.
    kept = _names(got[0], 1)
    assert kept == _names(got[1], 2)
    # A quarter of 400 pairs: 100, sd 8.7.
    assert 60 <= len(kept) <= 140
    gated = log.read_text().splitlines()
    assert all(line.startswith("dna:downsample -- ") for line in gated)
    # The fraction is the target over both files' bases (seqtk size), the seed 11.
    samples = [line for line in gated if "seqtk sample" in line]
    assert [
        ("-s 11" in s, f"{fq} 0.25 |" in s)
        for s, fq in zip(samples, (fq1, fq2), strict=True)
    ] == [
        (True, True),
        (True, True),
    ]
    assert sorted(p.name for p in out.iterdir()) == sorted([fq1.name, fq2.name])
    with pytest.raises(ValueError, match="downsample"):
        downsample(fq1, fq2, tmp_path / "again", target_bases=total)
    assert not (tmp_path / "again").exists()  # refused before writing anything


@pytest.mark.tools
@pytest.mark.parametrize("sex", ["male", "female"])
def test_run_sample_uses_run_args(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sex: Sex
) -> None:
    """HG003's CRAM and gVCF come from the run's own align and call steps: DeepVariant gets
    the arguments `genokeeper run` gives a sample of that sex, once the run's sex step
    infers that sex."""
    lab, fake, _ = _faked(tmp_path, monkeypatch)
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: sex)
    run = run_sample(
        lab.paths,
        name="hg003",
        fastq1=lab.fastq1,
        fastq2=lab.fastq2,
        sex=sex,
        grading=GRADING,
        threads=2,
    )
    cram, gvcf = run.cram, run.gvcf
    assert cram.is_file() and Path(f"{cram}.crai").is_file()
    assert gvcf.is_file() and gvcf.name.endswith(".g.vcf.gz")
    home = calibration_dir(lab.paths.lab) / "hg003"
    assert home in cram.parents and home in gvcf.parents
    expected = deepvariant_args(
        cram,
        grch38_fasta(lab.paths.lab),
        Path("/unused"),
        checkpoint=checkpoint_path(lab.paths.lab),
        male=sex == "male",
        par_bed=Path(PAR_BED),
        threads=shard_count(2, 64 * GIB),
    )
    assert fake.args == [expected]
    assert ("--haploid_contigs=chrX,chrY" in expected) == (sex == "male")
    # A rerun reuses both steps: nothing is aligned or called again.
    again = run_sample(
        lab.paths,
        name="hg003",
        fastq1=lab.fastq1,
        fastq2=lab.fastq2,
        sex=sex,
        grading=GRADING,
        threads=2,
    )
    assert again == run
    assert fake.calls == 1


@pytest.mark.tools
@pytest.mark.parametrize(
    ("thresholds", "refusal", "called"),
    [
        (replace(MINI_QC, min_mean_depth=1000.0), r"QC gate: mean_depth .* < 1000", 0),
        (replace(MINI_QC, min_titv=1e9), r"QC gate: titv .* < 1e\+09", 1),
    ],
    ids=["before calling", "after calling"],
)
def test_run_sample_applies_the_runs_qc_gate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    thresholds: QcThresholds,
    refusal: str,
    called: int,
) -> None:
    """HG003 goes through the run's QC gate: the alignment's metrics before calling, the
    calls' after, against the same [qc] thresholds a run is judged by."""
    lab, fake, verifybamid = _faked(tmp_path, monkeypatch)
    with pytest.raises(RunRefused, match=refusal):
        run_sample(
            lab.paths,
            name="hg003",
            fastq1=lab.fastq1,
            fastq2=lab.fastq2,
            sex="male",
            grading=replace(GRADING, qc=thresholds),
            threads=2,
        )
    assert (verifybamid.calls, fake.calls) == (1, called)


@pytest.mark.tools
def test_align_sample_matches_run_sample(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`align_sample` returns the CRAM `run_sample` aligns, its align step's key and the
    QC metrics, without calling."""
    lab, fake, _ = _faked(tmp_path, monkeypatch)
    kwargs = {
        "name": "hg003",
        "fastq1": lab.fastq1,
        "fastq2": lab.fastq2,
        "sex": "male",
        "grading": GRADING,
        "threads": 2,
    }
    cram, key, metrics = align_sample(lab.paths, **kwargs)  # type: ignore[arg-type]
    assert fake.calls == 0
    assert cram.parent.name == key and cram.is_file()
    assert metrics["mean_depth"] > 0
    run = run_sample(lab.paths, **kwargs)  # type: ignore[arg-type]
    assert (run.cram, run.key) == (cram, key)


@pytest.mark.tools
def test_run_sample_returns_key_and_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`run_sample` returns what a record keeps: the CRAM, its align step's key, the QC
    metrics with the call's Ti/Tv and het/hom, and the gVCF."""
    lab, _, _ = _faked(tmp_path, monkeypatch)
    run = run_sample(
        lab.paths,
        name="hg003",
        fastq1=lab.fastq1,
        fastq2=lab.fastq2,
        sex="male",
        grading=GRADING,
        threads=2,
    )
    assert run.key == run.cram.parent.name
    assert run.cram.parent.parent == calibration_dir(lab.paths.lab) / "hg003" / "align"
    assert run.metrics["total_bases"] == 2 * 400 * READ
    assert "titv" in run.metrics and "het_hom" in run.metrics
    assert run.gvcf.name.endswith(".g.vcf.gz") and run.gvcf.is_file()
