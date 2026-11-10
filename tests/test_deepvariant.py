from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pytest

from genokeeper.call import deepvariant
from genokeeper.call.deepvariant import (
    CALLING_FLAGS,
    CHECKPOINT_FILES,
    deepvariant_args,
    run_deepvariant,
)
from genokeeper.config import load_grading
from genokeeper.containers import IMAGES, PLATFORM, DockerUnavailable
from genokeeper.pipeline import CALL_STEP_VERSION
from genokeeper.tools import ToolFailed

CHECKPOINT = "weights-51-0.995354.ckpt"
GIB = 2**30


@pytest.fixture(autouse=True)
def docker_vm(monkeypatch: pytest.MonkeyPatch) -> None:
    """A Docker VM with room for every shard, unless a test says otherwise."""
    monkeypatch.setattr(deepvariant, "docker_memory", lambda: 64 * GIB)


@pytest.mark.parametrize(
    ("threads", "memory", "shards"),
    [
        (14, int(17.54 * GIB), 8),  # a small Docker VM: memory caps the 14 CPUs
        (14, int(46.96 * GIB), 14),  # a larger Docker VM: capped by the CPUs
        (14, 1 * GIB, 1),  # never fewer than one
        (2, 64 * GIB, 2),
    ],
)
def test_shards_capped_by_docker_memory(threads: int, memory: int, shards: int) -> None:
    assert deepvariant.shard_count(threads, memory) == shards


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Path]:
    """Stand-in inputs: a CRAM, a FASTA and the checkpoint's four files, each a placeholder."""
    paths = {
        "cram": tmp_path / "work" / "sample.cram",
        "fasta": tmp_path / "reference" / "GRCh38_no_alt.fa",
        "checkpoint": tmp_path / "deepvariant" / "complete-t7" / CHECKPOINT,
        "out": tmp_path / "calls",
    }
    for key in ("cram", "fasta"):
        paths[key].parent.mkdir(parents=True)
        paths[key].write_text("placeholder\n")
    paths["checkpoint"].parent.mkdir(parents=True)
    for name in CHECKPOINT_FILES:
        (paths["checkpoint"].parent / name).write_text("placeholder\n")
    paths["out"].mkdir()
    return paths


def _args(inputs: dict[str, Path], *, male: bool) -> list[str]:
    return deepvariant_args(
        inputs["cram"],
        inputs["fasta"],
        inputs["out"],
        checkpoint=inputs["checkpoint"],
        male=male,
        par_bed=inputs["out"] / "par38.bed",
        threads=6,
    )


def test_args_pinned_and_gvcf(inputs: dict[str, Path]) -> None:
    args = _args(inputs, male=False)
    assert args[0] == "/opt/deepvariant/bin/run_deepvariant"
    assert "--model_type=WGS" in args
    assert f"--customized_model=/model/{CHECKPOINT}" in args
    assert "--ref=/ref/GRCh38_no_alt.fa" in args
    assert "--reads=/reads/sample.cram" in args
    assert "--output_vcf=/out/sample.vcf.gz" in args
    assert "--output_gvcf=/out/sample.g.vcf.gz" in args
    # Every GQ value kept in the reference blocks.
    assert "--make_examples_extra_args=gvcf_gq_binsize=1" in args
    assert "--num_shards=6" in args


def test_male_gets_haploid_contigs(inputs: dict[str, Path]) -> None:
    args = _args(inputs, male=True)
    assert "--haploid_contigs=chrX,chrY" in args
    assert "--par_regions_bed=/par/par38.bed" in args


def test_female_has_none(inputs: dict[str, Path]) -> None:
    args = _args(inputs, male=False)
    assert not [a for a in args if a.startswith(("--haploid_contigs", "--par_regions"))]


@dataclass(frozen=True)
class Call:
    name: str
    args: list[str]
    mounts: dict[str, tuple[Path, bool]]  # container path -> (host path, read-only)
    platform: str | None
    label: str | None
    par_bed: str | None  # the PAR BED's text when the container ran


class FakeContainer:
    """Stands in for `run_container`: writes the outputs DeepVariant would into the /out
    mount (all four, or only `partial`), then raises `error` if one is set."""

    def __init__(self, partial: Sequence[str] | None = None) -> None:
        self.partial = partial
        self.error: Exception | None = None
        self.calls: list[Call] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        by_target = {target: (host, ro) for host, target, ro in mounts}
        par = by_target["/par"][0] / "par38.bed" if "/par" in by_target else None
        text = par.read_text() if par is not None else None
        self.calls.append(Call(name, list(args), by_target, platform, label, text))
        out = by_target["/out"][0]
        names = ["sample.vcf.gz", "sample.vcf.gz.tbi"]
        names += ["sample.g.vcf.gz", "sample.g.vcf.gz.tbi"]
        for written in self.partial if self.partial is not None else names:
            (out / written).write_text("calls\n")
        (out / "intermediate").mkdir()
        if self.error is not None:
            raise self.error


def _run(inputs: dict[str, Path], *, male: bool = True) -> tuple[Path, Path]:
    return run_deepvariant(
        inputs["cram"],
        inputs["fasta"],
        inputs["out"],
        checkpoint=inputs["checkpoint"],
        male=male,
        grading=load_grading(),
        threads=6,
    )


def test_run_passes_the_capped_shard_count(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeContainer()
    monkeypatch.setattr(deepvariant, "run_container", fake)
    monkeypatch.setattr(deepvariant, "docker_memory", lambda: int(5.5 * GIB))
    _run(inputs)  # six threads, room for two shards
    assert "--num_shards=2" in fake.calls[0].args


def test_calls_promoted_after_a_clean_run(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeContainer()
    monkeypatch.setattr(deepvariant, "run_container", fake)
    vcf, gvcf = _run(inputs)
    assert (vcf, gvcf) == (
        inputs["out"] / "sample.vcf.gz",
        inputs["out"] / "sample.g.vcf.gz",
    )
    assert sorted(p.name for p in inputs["out"].iterdir()) == [
        "sample.g.vcf.gz",
        "sample.g.vcf.gz.tbi",
        "sample.vcf.gz",
        "sample.vcf.gz.tbi",
    ]
    (call,) = fake.calls
    assert (call.name, call.platform, call.label) == (
        "deepvariant",
        PLATFORM,
        "dna:call",
    )
    assert call.mounts["/model"] == (inputs["checkpoint"].parent, True)
    assert call.mounts["/ref"] == (inputs["fasta"].parent, True)
    assert call.mounts["/reads"] == (inputs["cram"].parent, True)
    assert call.mounts["/out"][1] is False
    assert call.mounts["/par"][1] is True
    assert "--par_regions_bed=/par/par38.bed" in call.args
    # The male PAR BED is [par], 0-based half-open.
    assert call.par_bed == (
        "chrX\t10000\t2781479\nchrX\t155701382\t156030895\n"
        "chrY\t10000\t2781479\nchrY\t56887902\t57217415\n"
    )


def test_female_run_mounts_no_par(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeContainer()
    monkeypatch.setattr(deepvariant, "run_container", fake)
    _run(inputs, male=False)
    (call,) = fake.calls
    assert "/par" not in call.mounts
    assert not [a for a in call.args if a.startswith("--haploid_contigs")]


@pytest.mark.parametrize(
    "error",
    [
        ToolFailed(137, ["docker", "run"], None, "Killed"),
        DockerUnavailable("Docker is not running: start Docker Desktop and rerun"),
    ],
)
def test_failed_container_promotes_nothing(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch, error: Exception
) -> None:
    # DeepVariant killed after writing the VCF but before the gVCF: nothing reaches the
    # output directory, and no work directory stays behind to be taken for a result.
    fake = FakeContainer(partial=["sample.vcf.gz", "sample.vcf.gz.tbi"])
    fake.error = error
    monkeypatch.setattr(deepvariant, "run_container", fake)
    with pytest.raises(type(error)):
        _run(inputs)
    assert list(inputs["out"].iterdir()) == []


def test_missing_output_promotes_nothing(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        deepvariant, "run_container", FakeContainer(partial=["sample.vcf.gz"])
    )
    with pytest.raises(RuntimeError, match=r"sample\.g\.vcf\.gz"):
        _run(inputs)
    assert list(inputs["out"].iterdir()) == []


def test_missing_checkpoint_refused(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeContainer()
    monkeypatch.setattr(deepvariant, "run_container", fake)
    (inputs["checkpoint"].parent / f"{CHECKPOINT}.index").unlink()
    with pytest.raises(FileNotFoundError, match="bootstrap"):
        _run(inputs)
    assert fake.calls == []


def test_image_pinned_by_digest() -> None:
    assert IMAGES["deepvariant"] == (
        "google/deepvariant@sha256:"
        "962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6"
    )


# The calling flags decide the calls but are not in the call step's key, so each change of
# them must come with a new CALL_STEP_VERSION, or every run and calibration reuses calls
# made with the old flags (BR-043). Add the new pair here when you bump it.
FLAGS_BY_CALL_STEP = {
    "1": ("--model_type=WGS", "--make_examples_extra_args=gvcf_gq_binsize=1"),
}


def test_calling_flags_change_only_with_the_call_step_version() -> None:
    assert FLAGS_BY_CALL_STEP[CALL_STEP_VERSION] == CALLING_FLAGS
