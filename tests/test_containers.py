import os
import re
import stat
from pathlib import Path

import pytest

from genokeeper import containers, tools
from genokeeper.containers import (
    BUILT,
    DOCKER_DIR,
    IMAGES,
    PLATFORM,
    DockerUnavailable,
    build_command,
    build_inputs_digest,
    docker_memory,
    run_container,
)
from genokeeper.tools import ToolFailed

# Logs each call's arguments, one per line, then `--`; `info` fails when $DOCKER_DOWN is set,
# hangs for 3 s when $DOCKER_HANG is set; `run` exits with $DOCKER_RUN_EXIT (default 0) after
# writing to stderr.
FAKE_DOCKER = """#!/usr/bin/env bash
printf '%s\\n' "$@" -- >> "$DOCKER_LOG"
case "$1" in
  info) [ -z "${DOCKER_HANG:-}" ] || exec sleep 3
    [ -z "${DOCKER_DOWN:-}" ] || exit 1; [ "${2:-}" != --format ] || echo 18833149952 ;;
  run) echo "container said no" >&2; exit "${DOCKER_RUN_EXIT:-0}" ;;
  *) exit 1 ;;
esac
"""


def _executable(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def docker(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A stub `docker` first on PATH; returns its log. No gate."""
    _executable(tmp_path / "bin" / "docker", FAKE_DOCKER)
    log = tmp_path / "docker.log"
    monkeypatch.setenv("PATH", f"{tmp_path / 'bin'}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("DOCKER_LOG", str(log))
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")
    return log


def _calls(log: Path) -> list[list[str]]:
    calls: list[list[str]] = [[]]
    for line in log.read_text().splitlines():
        if line == "--":
            calls.append([])
        else:
            calls[-1].append(line)
    return calls[:-1]


def test_container_always_offline_and_pinned(docker: Path, tmp_path: Path) -> None:
    ref, out = tmp_path / "ref", tmp_path / "out"
    # Arguments after the image are the tool's: a `--network host` there changes nothing.
    run_container(
        "verifybamid2",
        ["verifybamid2", "--network", "host"],
        mounts=[(ref, "/ref", True), (out, "/out", False)],
        platform=PLATFORM,
    )
    info, run = _calls(docker)
    assert info == ["info"]
    image = run.index(IMAGES["verifybamid2"])
    assert run[:image] == [
        "run",
        "--rm",
        "--network",
        "none",
        "--pull",
        "never",
        "--platform",
        "linux/amd64",
        "--mount",
        f"type=bind,source={ref.resolve()},target=/ref,readonly",
        "--mount",
        f"type=bind,source={out.resolve()},target=/out",
    ]
    assert run[image + 1 :] == ["verifybamid2", "--network", "host"]
    assert "@sha256:" in run[image]


def test_docker_not_running_names_docker(
    docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("DOCKER_DOWN", "1")
    with pytest.raises(DockerUnavailable, match="Docker"):
        run_container("verifybamid2", ["true"], mounts=[])
    assert _calls(docker) == [["info"]]  # nothing run


def test_image_missing_asks_docker(
    docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An image is missing only when Docker runs and lacks it (the stub knows none)."""
    assert containers.image_missing("haplocheck")
    assert _calls(docker)[-1] == ["image", "inspect", IMAGES["haplocheck"]]
    monkeypatch.setenv("DOCKER_DOWN", "1")
    assert not containers.image_missing("haplocheck")


def test_a_hung_docker_is_unavailable(
    docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("DOCKER_HANG", "1")
    monkeypatch.setattr(containers, "DOCKER_INFO_TIMEOUT", 0.5)
    with pytest.raises(DockerUnavailable, match="did not answer"):
        run_container("verifybamid2", ["true"], mounts=[])
    assert _calls(docker) == [["info"]]  # nothing run


def test_docker_memory_is_the_vms_memtotal(
    docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert docker_memory() == 18833149952
    assert _calls(docker) == [["info", "--format", "{{.MemTotal}}"]]
    monkeypatch.setenv("DOCKER_DOWN", "1")
    with pytest.raises(DockerUnavailable, match="Docker"):
        docker_memory()


def test_failed_container_raises_with_its_stderr(
    docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("DOCKER_RUN_EXIT", "137")
    with pytest.raises(ToolFailed, match="container said no") as failed:
        run_container("verifybamid2", ["true"], mounts=[])
    assert failed.value.returncode == 137
    # Docker's own failure (125: an image missing, among others): nothing ran.
    monkeypatch.setenv("DOCKER_RUN_EXIT", "125")
    with pytest.raises(DockerUnavailable, match="could not start.*container said no"):
        run_container("verifybamid2", ["true"], mounts=[])


def test_one_gate_wraps_the_container(
    docker: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gate_log = tmp_path / "gate.log"
    gate = _executable(
        tmp_path / "coreslot-gate",
        f'#!/usr/bin/env bash\necho "$1 $2" >> "{gate_log}"\nshift 2\nexec "$@"\n',
    )
    monkeypatch.setattr(tools, "GATE", gate)
    run_container("verifybamid2", ["true"], mounts=[], label="dna:qc")
    assert gate_log.read_text().splitlines() == ["dna:qc --"]
    assert [call[0] for call in _calls(docker)] == ["info", "run"]


def test_a_mount_source_docker_would_misread_is_refused(docker: Path) -> None:
    with pytest.raises(ValueError, match="comma"):
        run_container("verifybamid2", ["true"], mounts=[(Path("/a,b"), "/x", True)])


def test_built_images_named_by_their_inputs() -> None:
    """An image built here (no BioContainers image for its pinned version) is named by the
    sha256 of its Dockerfile and the lists it copies: an edited input names another image,
    which the bootstrap builds. Every file under docker/ is such an input, and every base
    image is pinned by digest."""
    inputs: set[str] = set()
    for name, files in BUILT.items():
        assert files[0] == f"{name}.Dockerfile"
        repository, tag = IMAGES[name].split(":", 1)
        # The namespace from before the rename: a new one would rebuild every image.
        assert repository == f"genomelab/{name}"
        assert tag.endswith(f"-{build_inputs_digest(name)}")
        assert len(tag) <= 128  # Docker's limit
        inputs.update(files)
        froms = [
            line.split()[1]
            for line in (DOCKER_DIR / files[0]).read_text().splitlines()
            if line.startswith("FROM ")
        ]
        assert froms and all("@sha256:" in image for image in froms)
    assert {p.name for p in DOCKER_DIR.iterdir()} == inputs
    assert build_command("aldy") == [
        "docker",
        "build",
        "--tag",
        IMAGES["aldy"],
        "--file",
        str(DOCKER_DIR / "aldy.Dockerfile"),
        str(DOCKER_DIR),
    ]


def test_inputs_digest_moves_with_any_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for name in BUILT["bcyrius"]:
        (tmp_path / name).write_bytes((DOCKER_DIR / name).read_bytes())
    monkeypatch.setattr(containers, "DOCKER_DIR", tmp_path)
    before = build_inputs_digest("bcyrius")
    assert before == IMAGES["bcyrius"].rsplit("-", 1)[1]
    (tmp_path / "bcyrius.requirements.txt").write_text("numpy==1\n")
    assert build_inputs_digest("bcyrius") != before


def test_a_build_from_other_inputs_is_refused(
    docker: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An edited input would build another image under the pinned name, which
    `image_pinned` then trusts after the file is reverted."""
    inputs = tmp_path / "docker"
    for name in BUILT["bcyrius"]:
        _executable(inputs / name, (DOCKER_DIR / name).read_text())
    monkeypatch.setattr(containers, "DOCKER_DIR", inputs)
    containers.main(["build", "bcyrius"])
    assert [call[0] for call in _calls(docker)] == ["build"]
    (inputs / "bcyrius.requirements.txt").write_text("numpy==1\n")
    assert containers.main(["build", "bcyrius"]) == 1
    assert [call[0] for call in _calls(docker)] == ["build"]  # no second build


def _lock_lines(lock: str, package: str) -> list[str]:
    """The lines of docker/`lock` for the conda package named `package` (name-version-build)."""
    return [
        line
        for line in (DOCKER_DIR / lock).read_text().splitlines()
        if line.rsplit("/", 1)[-1].startswith(f"{package}.")
    ]


def _assert_conda_image_pinned(name: str, platform: str, packages: list[str]) -> None:
    assert name in BUILT and name in IMAGES
    assert IMAGES[name].endswith(f"-{build_inputs_digest(name)}")
    assert containers.PLATFORMS[name] == platform
    assert BUILT[name] == (f"{name}.Dockerfile", f"{name}.lock")
    lock = (DOCKER_DIR / f"{name}.lock").read_text().splitlines()
    assert lock[0] == "@EXPLICIT"
    for package in packages:
        (line,) = _lock_lines(f"{name}.lock", package)
        assert re.fullmatch(r"https://\S+#[0-9a-f]{32}", line)


def test_smncaller_image_pinned() -> None:
    _assert_conda_image_pinned(
        "smncaller", "linux/arm64", ["smncopynumbercaller-1.1.2-py312h7e72e81_1"]
    )


def test_parascopy_image_pinned() -> None:
    _assert_conda_image_pinned(
        "parascopy", "linux/arm64", ["parascopy-1.19.0-py312haa49b00_0"]
    )


def test_yleaf_image_pinned() -> None:
    _assert_conda_image_pinned(
        "yleaf",
        "linux/arm64",
        ["samtools-1.24-h391949c_1", "python-3.12.15-h94ad73b_2_cpython"],
    )
    dockerfile = (DOCKER_DIR / "yleaf.Dockerfile").read_text()
    assert re.search(r"genid/Yleaf\.git#[0-9a-f]{40}\b", dockerfile)
    assert "ENV YLEAF_NO_UPDATE_CHECK=1" in dockerfile


def test_expansionhunter_image_pinned() -> None:
    """ExpansionHunter is built here from Illumina's source at the last Apache-2.0 commit,
    the one with the RFC1 motif analysis, not taken from bioconda's 5.0.0."""
    name = "expansionhunter"
    assert name in BUILT and name in IMAGES
    assert IMAGES[name].startswith("genomelab/expansionhunter:9e93322f4a-")
    assert IMAGES[name].endswith(f"-{build_inputs_digest(name)}")
    assert containers.PLATFORMS[name] == "linux/arm64"
    assert BUILT[name] == ("expansionhunter.Dockerfile", "expansionhunter.lock")
    dockerfile = (DOCKER_DIR / "expansionhunter.Dockerfile").read_text()
    lines = dockerfile.replace("\\\n", " ").splitlines()  # a backslash joins its lines
    assert (
        "github.com/Illumina/ExpansionHunter.git"
        "#9e93322f4ae72d450a277c43fc6dd54e5b610a95" in dockerfile
    )
    adds = [
        line
        for line in lines
        if re.match(r"ADD --checksum=sha256:[0-9a-f]{64}\b", line)
    ]
    assert len(adds) == 5
    assert sum("archives.boost.io" in line for line in adds) == 1
    assert not any("jfrog" in line for line in adds)
    offline = [line for line in lines if line.startswith("RUN --network=none")]
    assert any("cmake" in line and "make" in line for line in offline)
    assert any("RFC1MotifAnalysis" in line for line in offline)
    # The lab's CRAMs are 3.1, which htslib reads from 1.12 on: the commit's own 1.10.2 aborts.
    assert "htslib-1.22.1" in dockerfile
    # The commit's 1.10.2 URL stays only as the sed pattern that finds it in CMakeLists.txt, and
    # is replaced by the local 1.22.1 file; nothing fetches or names it otherwise.
    code = [x for x in lines if not x.startswith("#")]
    pattern = "https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2"
    assert not any("1.10.2" in x.replace(pattern + "|", "") for x in code)
    assert f"{pattern}|/opt/eh-deps/htslib-1.22.1.tar.bz2|" in " ".join(code)
    cram31 = [line for line in offline if "test_view" in line and "VERSION=3.1" in line]
    assert len(cram31) == 1 and "-C" in cram31[0] and "-D" in cram31[0]
    assert any("grep" in line and "1.22.1" in line for line in offline)
    # Two stages on one pinned base: the final one carries the binary and its conda
    # libraries, and checks offline that it runs.
    froms = [line for line in lines if line.startswith("FROM ")]
    assert len(froms) == 2 and froms[0].endswith(" AS build")
    assert froms[0].removesuffix(" AS build") == froms[1]
    assert re.fullmatch(r"FROM \S+@sha256:[0-9a-f]{64}", froms[1])
    final = lines[lines.index(froms[1]) :]
    assert "COPY --from=build /tmp/eh/ /opt/eh/" in final
    assert not any(line.startswith(("ADD", "COPY expansionhunter")) for line in final)
    checks = [line for line in final if line.startswith("RUN --network=none")]
    assert len(checks) == 1 and "ExpansionHunter --help" in checks[0]
    assert all(
        line.startswith("RUN --network=none") for line in final if "RUN" in line[:3]
    )
    lock = (DOCKER_DIR / "expansionhunter.lock").read_text().splitlines()
    assert lock[0] == "@EXPLICIT"
    assert all(re.fullmatch(r"https://\S+#[0-9a-f]{32}", line) for line in lock[1:])
    names = [line.rsplit("/", 1)[-1] for line in lock[1:]]
    assert sum(n.startswith("gxx_linux-aarch64-10.") for n in names) == 1
    assert sum(n.startswith("cmake-3.31.") for n in names) == 1
    assert not any(n.startswith("expansionhunter-") for n in names)


def test_reviewer_image_pinned() -> None:
    _assert_conda_image_pinned("reviewer", "linux/arm64", ["reviewer-0.2.7-h739925f_0"])


def test_delly_image_pinned() -> None:
    _assert_conda_image_pinned("delly", "linux/arm64", ["delly-2.7.0-ha658904_0"])
    dockerfile = (DOCKER_DIR / "delly.Dockerfile").read_text()
    assert re.search(r"dellytools/delly/[0-9a-f]{40}/excludeTemplates/", dockerfile)
    assert re.search(r"--checksum=sha256:[0-9a-f]{64}\b", dockerfile)


def test_dysgu_image_pinned() -> None:
    _assert_conda_image_pinned("dysgu", "linux/arm64", ["dysgu-1.9.0-py312hed307a4_0"])


def test_truvari_image_pinned() -> None:
    _assert_conda_image_pinned(
        "truvari",
        "linux/arm64",
        [
            "truvari-5.5.0-pyhdfd78af_1",
            "bcftools-1.25-he84ed4f_0",
            "htslib-1.24-h190289e_0",
        ],
    )


def _first_comment(name: str) -> str:
    return (DOCKER_DIR / f"{name}.Dockerfile").read_text().split("\nFROM ", 1)[0]


def test_cnvpytor_image_pinned() -> None:
    _assert_conda_image_pinned(
        "cnvpytor", "linux/arm64", ["cnvpytor-1.3.2-pyhdfd78af_0"]
    )
    dockerfile = (DOCKER_DIR / "cnvpytor.Dockerfile").read_text()
    assert "3adc56fc6433ab8d76db1f151720cc4a4e487f20" in dockerfile
    assert len(re.findall(r"--checksum=sha256:[0-9a-f]{64}\b", dockerfile)) == 2
    assert "MIT" in _first_comment("cnvpytor")


def test_insurveyor_image_pinned() -> None:
    _assert_conda_image_pinned(
        "insurveyor", "linux/arm64", ["insurveyor-1.1.3-h1e84f2d_2"]
    )
    comment = _first_comment("insurveyor")
    assert "GPL-3.0" in comment and "run, never redistributed" in comment


def test_scramble_image_pinned() -> None:
    _assert_conda_image_pinned(
        "scramble",
        "linux/arm64",
        [
            "htslib-1.25-h190289e_0",
            "r-base-4.4.3-h269e6ae_11",
            "bioconductor-biostrings-2.74.0-r44h59b44dc_1",
            "bioconductor-pwalign-1.2.0-r44h59b44dc_0",
            "bioconductor-rblast-1.2.0-r44hdfd78af_0",
            "blast-2.16.0-h6a93c2d_5",
        ],
    )
    dockerfile = (DOCKER_DIR / "scramble.Dockerfile").read_text()
    comment = _first_comment("scramble")
    assert "CC BY-NC-SA 4.0" in comment and "run, never redistributed" in comment
    assert re.search(r"GeneDx/scramble\.git#[0-9a-f]{40}\b", dockerfile)


def test_gauchian_image_pinned() -> None:
    _assert_conda_image_pinned(
        "gauchian",
        "linux/arm64",
        ["python-3.12.16-h94ad73b_0_cpython", "pysam-0.24.1-py312hfd8969c_0"],
    )
    dockerfile = (DOCKER_DIR / "gauchian.Dockerfile").read_text()
    assert re.search(r"Illumina/Gauchian\.git#[0-9a-f]{40}\b", dockerfile)
    assert "GPLv3" in _first_comment("gauchian")


def test_kilda_image_pinned() -> None:
    _assert_conda_image_pinned(
        "kilda",
        "linux/arm64",
        ["kmer-jellyfish-2.3.1-py312pl5321h42fd94b_6", "samtools-1.25-h391949c_0"],
    )
    dockerfile = (DOCKER_DIR / "kilda.Dockerfile").read_text()
    assert re.search(r"HCL-HUBL/KILDA\.git#[0-9a-f]{40}\b", dockerfile)
    assert "AGPL-3.0" in _first_comment("kilda")
    assert "kmer-jellyfish-2.3.1" in (DOCKER_DIR / "kilda.lock").read_text()


def _assert_hashed_pip_image_pinned(name: str, platform: str) -> None:
    assert name in BUILT and name in IMAGES
    assert IMAGES[name].endswith(f"-{build_inputs_digest(name)}")
    assert containers.PLATFORMS[name] == platform
    assert BUILT[name][:2] == (f"{name}.Dockerfile", f"{name}.requirements.txt")


def _pinned_lines(requirements: str) -> list[str]:
    return [
        line
        for line in (DOCKER_DIR / requirements).read_text().splitlines()
        if line and not line.startswith("#")
    ]


def test_haplocheck_image_pinned() -> None:
    _assert_conda_image_pinned("haplocheck", "linux/arm64", [])
    (line,) = [
        line
        for line in (DOCKER_DIR / "haplocheck.lock").read_text().splitlines()
        if "/haplocheck-1.3.3-" in line
    ]
    assert re.fullmatch(r"https://\S+/haplocheck-1\.3\.3-\S+#[0-9a-f]{32}", line)


def test_rbceq2_image_pinned() -> None:
    _assert_hashed_pip_image_pinned("rbceq2", "linux/arm64")
    lines = _pinned_lines("rbceq2.requirements.txt")
    assert all("--hash=sha256:" in line for line in lines)
    assert (
        "rbceq2==2.4.4 --hash=sha256:"
        "9091276c8d954081d673e6fd254936f26e7d4ff59edb581e3d8c1a089ded576c"
    ) in lines
    dockerfile = (DOCKER_DIR / "rbceq2.Dockerfile").read_text()
    assert re.search(
        r"^LABEL rbceq2\.db=\S+/rbceq2/resources/db\.tsv$", dockerfile, re.MULTILINE
    )


def test_fastmixture_image_pinned() -> None:
    _assert_hashed_pip_image_pinned("fastmixture", "linux/arm64")
    dockerfile = (DOCKER_DIR / "fastmixture.Dockerfile").read_text()
    assert "GPL-3.0: run, never\n# redistributed" in dockerfile
    assert re.search(
        r"ADD --checksum=sha256:[0-9a-f]{64} \\\s+\S+/fastmixture-1\.3\.0\.tar\.gz",
        dockerfile,
    )
    assert all(
        "--hash=sha256:" in line
        for line in _pinned_lines("fastmixture.requirements.txt")
    )


def test_hmmix_image_pinned() -> None:
    assert "hmmix" in BUILT and "hmmix" in IMAGES
    assert IMAGES["hmmix"].endswith(f"-{build_inputs_digest('hmmix')}")
    assert containers.PLATFORMS["hmmix"] == "linux/arm64"
    assert BUILT["hmmix"] == (
        "hmmix.Dockerfile",
        "hmmix.lock",
        "hmmix.requirements.txt",
    )
    for package in ("python-3.10.", "bcftools-1.25-", "vcftools-0.1.17-"):
        (line,) = [
            line
            for line in (DOCKER_DIR / "hmmix.lock").read_text().splitlines()
            if f"/{package}" in line
        ]
        assert re.fullmatch(r"https://\S+#[0-9a-f]{32}", line)
    lines = _pinned_lines("hmmix.requirements.txt")
    assert all("--hash=sha256:" in line for line in lines)
    assert any(line.startswith("hmmix==0.9.2 --hash=sha256:") for line in lines)
    dockerfile = (DOCKER_DIR / "hmmix.Dockerfile").read_text()
    assert "run, never\n# redistributed" in dockerfile
    assert "vcftools --version" in dockerfile
