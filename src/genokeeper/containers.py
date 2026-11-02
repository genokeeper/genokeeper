"""The container images the genome tools run in, each pinned to a digest, and the one way to
run them.

Every third-party tool that reads the user's reads or calls runs in one of these images
(spec, "Containment of genome tools"). `scripts/bootstrap.sh` (step_images) pulls each
reference below by its digest, or builds it (`BUILT`, `build_command`) where no image of the
pinned version exists for the platform it runs on (each `IMAGES` comment says why).

Four controls, each with its anchor:
- The pin, against an image swapped upstream: its anchor comes from the digest written here,
  read from the image's registry on the day it was pinned, which a re-pushed tag cannot move
  because Docker resolves `repo@sha256:<digest>` by content and refuses other content.
- The pins of a built image, against a changed download: their anchor comes from the base
  image digests, the source's sha256 or git commit and the wheels' sha256s written in its
  Dockerfile and lists under docker/, which a re-uploaded file cannot pass because the build
  checks each before it uses it. The image's own digest is not pinned: a build is not
  byte-reproducible across Docker versions, and its digest would stop no adversary the input
  pins do not (the 2a review's ruling).
- The name of a built image, against running an image built from other inputs: its anchor
  comes from the sha256 of those inputs written in its tag here (a lock file's without its
  comment lines and blank lines, which install nothing), which an edited Dockerfile or pin
  cannot keep because `main` recomputes it from the files before every build and refuses a
  mismatch (and `tests/test_containers.py` from the committed files).
- `--network none`, against a tool's own network call carrying genotype-derived data: its
  anchor comes from our `docker run` command line, which the tool inside the container cannot
  change because the network namespace is set before its process starts; the tool's
  arguments follow the image, where Docker reads no option of its own.
"""

import hashlib
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

from genokeeper.tools import ToolFailed, gated, run_tool

# The images published for linux/amd64 only; on Apple silicon each runs under Rosetta.
PLATFORM = "linux/amd64"

# Image name -> `repo@sha256:<digest>`. The built images keep the `genomelab/` namespace and
# docker/ its `genomelab.containers` comments from before the rename to genokeeper: the
# namespace is in the step keys and calibration versions that name an image, and the files'
# digest in each tag (`build_inputs_digest`), so renaming either would rebuild every image
# and re-run its steps.
IMAGES: dict[str, str] = {
    # DeepVariant 1.10.0, Google's own image: the digest Docker Hub's API gives for tag
    # 1.10.0 (hub.docker.com/v2/repositories/google/deepvariant/tags/1.10.0, read
    # 2026-10-19; one linux/amd64 manifest, pushed 2026-03-02).
    "deepvariant": "google/deepvariant@sha256:"
    "962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6",
    # VerifyBamID 2.0.3. Upstream (github.com/Griffan/VerifyBamID/releases) publishes
    # binaries, no image; BioContainers builds this one by bioconda from the v2.0.3 source
    # tarball (sha256 ac0189793e85eb1982f7ce7e801957b540d85caa3d7de07dcbca3f9f819e0ea6) and
    # ships the 1000 Genomes SVD resources under /usr/local/share/verifybamid2-2.0.3-0/.
    # Digest of tag 2.0.3--hc004090_0 from the quay.io API, 2026-10-19.
    "verifybamid2": "quay.io/biocontainers/verifybamid2@sha256:"
    "7f8b271e644c690fc870521e49bf542b3df1c392eb715f1e7f123e27e9d2da3c",
    # mutserve 2.0.3 and WhatsHap 2.8. Upstream publishes no image (mutserve a release zip,
    # WhatsHap PyPI wheels); BioContainers builds each by bioconda from that release. Digests
    # of tags 2.0.3--hdfd78af_0 and 2.8--py312hf731ba3_0 from the quay.io API, 2026-10-20.
    "mutserve": "quay.io/biocontainers/mutserve@sha256:"
    "d92233d837b02c1426b4472d9f33d61cd99f357acb8b61acae9bfa52f06c1756",
    "whatshap": "quay.io/biocontainers/whatshap@sha256:"
    "e50ebc2cc57c796b41248139c384d8c9d6885aa7713e3685058ddeaa9c856d78",
    # PharmCAT 3.4.0, PharmGKB's own image (github.com/PharmGKB/PharmCAT release v3.4.0,
    # 2026-07-14). Published for linux/amd64 and linux/arm64 under one OCI index: this is the
    # index's digest, the one Docker Hub's API gives for tag 3.4.0 and the registry's
    # Docker-Content-Digest (both read 2026-10-20; the index body hashes to it). The index
    # lists each platform's manifest by its own digest (amd64 sha256:5f742b22...,
    # arm64 sha256:db0d939c...), so pinning it pins both, and Docker runs the host's: native
    # on Apple silicon, no Rosetta for the JVM.
    "pharmcat": "pgkb/pharmcat@sha256:"
    "bc498d55c33094002bd7e14221a6b3fcd7243a587675369bc1a5ecf6ac0e5319",
    # OptiType 1.5.0 (GitHub tag v1.5.0; with RazerS3, YARA, GLPK and CBC) and T1K 1.0.10:
    # BioContainers builds each by bioconda from that release's GitHub tarball (sha256s
    # a59e917d... and cea98e29... in its recipes). Digests of tags 1.5.0--pyhdfd78af_1 and
    # 1.0.10--h5814d7d_0 from the quay.io API, 2026-10-20 (linux/amd64 only).
    "optitype": "quay.io/biocontainers/optitype@sha256:"
    "ce53ef12449bb8492971299044c8352ad46a43b703a86ebee155d5b99159420e",
    "t1k": "quay.io/biocontainers/t1k@sha256:"
    "250f92b6eb37f0df2abc066c4c9ddcdfe4084d9c565475d99a022713952f9a95",
    # Aldy 4.8.3 and BCyrius 2.2.0: no BioContainers or upstream image, so built here from
    # docker/ (`BUILT`) for the host's platform, each named `<version>-<inputs' sha256>`.
    "aldy": "genomelab/aldy:4.8.3-"
    "61121bdf3ef0411eaf135a9de48b7129035dd73036cb6c8e54641e85b7572aab",
    "bcyrius": "genomelab/bcyrius:2.2.0-"
    "e7d4c07589c8e635c7724c68574a08b8949de7b72a26c7e9cc007862e42a0b3e",
    # Haplogrep3 3.3.2: upstream publishes a release tarball, no image, and BioContainers
    # builds 3.2.2 at newest (quay.io, read 2026-10-21); built here, with its tree baked in
    # (docker/haplogrep3.Dockerfile).
    "haplogrep3": "genomelab/haplogrep3:3.3.2-"
    "3733ddceb7ab3bff017e288b2e5a6586084e6dd29ea85142c2c069d544e2a2b7",
    # pgsc_calc v2.3.0's images, one per container label of its conf/modules.config at the
    # commit the tag points to (72ee54f2): pygscatalog pgscatalog-utils-2.0.4, plink2
    # 2.00a5.10, zstd 2-beta, report 2-beta, pyyaml 6.0, fraposa_pgsc v1.0.2. Each pinned
    # to its linux/arm64 manifest, the digest ghcr.io's registry API lists for that
    # platform in the tag's index (read 2026-10-21; the manifest bodies hash to them), so
    # pgsc_calc's arm profile runs them natively. Nextflow starts them
    # (`genokeeper.prs.pgsc.nextflow_config`); `run_container` only plink2's, for the
    # panel's ancestry positions.
    "pgsc_utils": "ghcr.io/pgscatalog/pygscatalog@sha256:"
    "91e96c0543c26b9e1bccbcc5b39fb5fe1d74ba3aa90d54278fa1075f7fed939b",
    "pgsc_plink2": "ghcr.io/pgscatalog/plink2@sha256:"
    "fa92fba77492257e9096ac9d8e27c66b6108c80a2fd33c1580642b8b39d28498",
    "pgsc_zstd": "ghcr.io/pgscatalog/zstd@sha256:"
    "bd3dd3e8c7d30f9e3d4b161f9f87108ae68fc50d357b02205364f461ef1eccfb",
    "pgsc_report": "ghcr.io/pgscatalog/report@sha256:"
    "53e9094b6c3cb5f7757e63a16c098f0a4be359f55de55054fef50dfd105e88fb",
    "pgsc_pyyaml": "ghcr.io/pgscatalog/pyyaml@sha256:"
    "67287d2d3c87b4de25040f3f0268d77c71332efc24b89140bf96770244e1cc2e",
    "pgsc_fraposa": "ghcr.io/pgscatalog/fraposa_pgsc@sha256:"
    "8fca7722c94ebd8e7b4a8b9f6eb603c12641a941a966dbf927cbe16a36710c1d",
    # Ensembl VEP 116.3 with its plugins (VEP_plugins release/116, baked in at /plugins),
    # Ensembl's own image: the OCI index digest Docker Hub's API gives for tag
    # release_116.3 and the registry's Docker-Content-Digest (both read 2026-10-23; pushed
    # 2026-10-09). The index lists linux/amd64 and linux/arm64 manifests, so it runs
    # natively on Apple silicon.
    "vep": "ensemblorg/ensembl-vep@sha256:"
    "2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3",
    # UCSC's bigBedToBed 482 (bioconda's ucsc-bigbedtobed, built from UCSC's kent source
    # by BioContainers), for the bootstrap's NMD Escape track. Digest of tag
    # 482--h0b57e2e_0 from the quay.io API, 2026-10-23 (linux/amd64 only).
    "ucsc_bigbedtobed": "quay.io/biocontainers/ucsc-bigbedtobed@sha256:"
    "aea2a8e7a31957eebece9a309b7aa1e487f101240eb71a4925e99a7e4c580bf5",
    # SpliceAI 1.3.1: Illumina publishes PyPI files and no image; built here from docker/
    # (`BUILT`) for linux/arm64 (`PLATFORMS`).
    "spliceai": "genomelab/spliceai:1.3.1-"
    "ea4699ff4fe9c239bd2fc38150926403d9d3871c6a4f76917f5159de5d0fc676",
    # SMNCopyNumberCaller 1.1.2 and Parascopy 1.19.0: BioContainers builds both for
    # linux/amd64 only, so built here from docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`),
    # each from an explicit conda lockfile; two images because their htslib pins differ.
    "smncaller": "genomelab/smncaller:1.1.2-"
    "7a27cc8b40616910ee1a38001a86535035f0da1ab93b52cdfa9631a34b6c7ccb",
    "parascopy": "genomelab/parascopy:1.19.0-"
    "3254559caa1352ef43e4b7f5c0c2447b2fb606b96de0e9b0d1b3976104e006b9",
    # ExpansionHunter (Illumina's source at commit 9e93322f4a, the last Apache-2.0 one, which
    # has the RFC1 motif analysis; it reports v5.0.0): built here from docker/ (`BUILT`) for
    # linux/arm64 (`PLATFORMS`), on a toolchain from an explicit conda lockfile.
    "expansionhunter": "genomelab/expansionhunter:9e93322f4a-"
    "c672ce04a77ef09c526737194c4fc169769e1e7a55c9c86a91f412b3dfb61332",
    # REViewer 0.2.7: BioContainers builds it for linux/amd64 only, so built here from docker/
    # (`BUILT`) for linux/arm64 (`PLATFORMS`), from an explicit conda lockfile.
    "reviewer": "genomelab/reviewer:0.2.7-"
    "39ba9362bcf3b27aa7996c4bb3fc98ffa4b4239b3633ae06c5e33dfcd03bf799",
    # Yleaf 4.1.6 (github.com/genid/Yleaf tag v4.1.6): no image upstream, so built here from
    # docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`), with its tree and position files.
    "yleaf": "genomelab/yleaf:4.1.6-"
    "84bce812b5fae12df9120850a55bde2fdc32c758a3f37b097a04c32d01d0b147",
    # DELLY 2.7.0 (github.com/dellytools/delly tag v2.7.0): BioContainers' image is not used,
    # so built here from docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`), from an explicit
    # conda lockfile; the image also carries the GRCh38 exclude list at
    # /opt/delly/human.hg38.excl.tsv.
    "delly": "genomelab/delly:2.7.0-"
    "2eb0c82c51c084ff964e7dcb70c24fb40537afc995ecb66aa5307348be9cfb96",
    # Dysgu 1.9.0 (github.com/kcleal/dysgu tag v1.9.0): built here from docker/ (`BUILT`) for
    # linux/arm64 (`PLATFORMS`), from an explicit conda lockfile.
    "dysgu": "genomelab/dysgu:1.9.0-"
    "6e14f4246365e22b84d2259e62bb1f486aba81e41781bbcd69b4972e4ba74568",
    # Truvari 5.5.0 (github.com/ACEnglish/truvari tag v5.5.0): built here from docker/
    # (`BUILT`) for linux/arm64 (`PLATFORMS`), from an explicit conda lockfile; it carries
    # bcftools and tabix.
    "truvari": "genomelab/truvari:5.5.0-"
    "3eae51b735f80cb2c2fd61a881f0ddb1776574f0cec444b9f955cdb1a1b4c801",
    # CNVpytor 1.3.2 (github.com/abyzovlab/CNVpytor tag v1.3.2), INSurVeyor 1.1.3
    # (github.com/kensung-lab/INSurVeyor tag 1.1.3) and Scramble 1.0.2
    # (github.com/GeneDx/scramble tag 1.0.2): BioContainers builds the first two for
    # linux/amd64 only and bioconda has no linux-aarch64 Scramble, so built here from docker/
    # (`BUILT`) for linux/arm64 (`PLATFORMS`), each from an explicit conda lockfile.
    # CNVpytor's image also carries its GRCh38 GC and mask files inside the package, where it
    # looks for them; Scramble's compiles cluster_identifier from the tag's C sources and
    # keeps its R scripts and MEI consensus at /opt/scramble/cluster_analysis/.
    "cnvpytor": "genomelab/cnvpytor:1.3.2-"
    "55ae847dc92f14854c822fe48d2a9c42be608969f55d7b00dc85c10cd4213dc0",
    "insurveyor": "genomelab/insurveyor:1.1.3-"
    "77d3b37f54cb51a565402d36de7fc246ca6b3df7b3e696841baafc97f28ccc8c",
    "scramble": "genomelab/scramble:1.0.2-"
    "c2f5ebacca43a7179433e599414899af3303ccb35180e996d4a386780caa25a4",
    # Gauchian 1.0.2 (github.com/Illumina/Gauchian tag v1.0.2) and KILDA 1.2.3
    # (github.com/HCL-HUBL/KILDA tag v1.2.3): neither publishes an image, so built here from
    # docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`), each from an explicit conda lockfile;
    # Gauchian is pip-installed from its commit with its GRCh38 data inside the package, KILDA's
    # script and curated k-mer lists are kept at /opt/kilda/ beside jellyfish and samtools.
    "gauchian": "genomelab/gauchian:1.0.2-"
    "088d0619615e9dd9fbf23a9b430f71f91fda938fb1c826d619d6a2d815b0aced",
    "kilda": "genomelab/kilda:1.2.3-"
    "e4ca2966abbcc00c160d55d6bf7be5ea035e28184f8a17a28782c7fb30cffa31",
    # haplocheck 1.3.3 (github.com/genepi/haplocheck tag v1.3.3): BioContainers' image is
    # linux/amd64 only, so built here from docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`),
    # from an explicit conda lockfile.
    "haplocheck": "genomelab/haplocheck:1.3.3-"
    "1eb4f00b4d73c5e2927a078892fc0ecdce9bfd378f9edf67495f34adb3067d08",
    # RBCeq2 2.4.4 (github.com/limcintyre/RBCeq2 tag v2.4.4): no image upstream or at
    # BioContainers, so built here from docker/ (`BUILT`) for linux/arm64 (`PLATFORMS`), from
    # a hashed wheel list; the image carries its ISBT allele table (LABEL rbceq2.db).
    "rbceq2": "genomelab/rbceq2:2.4.4-"
    "44e37c5bbfb52fb91e0f3745c071c3cb9b8dbe2fc6f29c2caf07e12448427885",
    # fastmixture 1.3.0 (github.com/Rosemeis/fastmixture tag v1.3.0; GPL-3.0: run, never
    # redistributed): no image upstream or at BioContainers, so built here from docker/
    # (`BUILT`) for linux/arm64 (`PLATFORMS`), from its PyPI sdist by sha256.
    "fastmixture": "genomelab/fastmixture:1.3.0-"
    "9a20b1dc54ca0cfd604d4bf1da63b917573da7f7a5e124b128b450429e8b8f02",
    # hmmix 0.9.2 (github.com/LauritsSkov/Introgression-detection; no tags): no image
    # upstream or at BioContainers, so built here from docker/ (`BUILT`) for linux/arm64
    # (`PLATFORMS`), from an explicit conda lockfile (Python 3.10, bcftools, vcftools) and a
    # hashed wheel list.
    "hmmix": "genomelab/hmmix:0.9.2-"
    "ffe9a70f04f2f1a3d400ae2686dd8911c15cf2d6ccec98bdb162f407dc85b363",
}
# Images built here, not pulled: name -> the files under docker/ the build reads, its
# Dockerfile first.
BUILT: dict[str, tuple[str, ...]] = {
    "aldy": (
        "aldy.Dockerfile",
        "aldy.build-requirements.txt",
        "aldy.requirements.txt",
    ),
    "bcyrius": ("bcyrius.Dockerfile", "bcyrius.requirements.txt"),
    "haplogrep3": ("haplogrep3.Dockerfile", "haplogrep3.yaml"),
    "spliceai": ("spliceai.Dockerfile", "spliceai.requirements.txt"),
    "smncaller": ("smncaller.Dockerfile", "smncaller.lock"),
    "parascopy": ("parascopy.Dockerfile", "parascopy.lock"),
    "yleaf": ("yleaf.Dockerfile", "yleaf.lock"),
    "expansionhunter": ("expansionhunter.Dockerfile", "expansionhunter.lock"),
    "reviewer": ("reviewer.Dockerfile", "reviewer.lock"),
    "delly": ("delly.Dockerfile", "delly.lock"),
    "dysgu": ("dysgu.Dockerfile", "dysgu.lock"),
    "truvari": ("truvari.Dockerfile", "truvari.lock"),
    "cnvpytor": ("cnvpytor.Dockerfile", "cnvpytor.lock"),
    "insurveyor": ("insurveyor.Dockerfile", "insurveyor.lock"),
    "scramble": ("scramble.Dockerfile", "scramble.lock"),
    "gauchian": ("gauchian.Dockerfile", "gauchian.lock"),
    "kilda": ("kilda.Dockerfile", "kilda.lock"),
    "haplocheck": ("haplocheck.Dockerfile", "haplocheck.lock"),
    "rbceq2": ("rbceq2.Dockerfile", "rbceq2.requirements.txt"),
    "fastmixture": ("fastmixture.Dockerfile", "fastmixture.requirements.txt"),
    "hmmix": ("hmmix.Dockerfile", "hmmix.lock", "hmmix.requirements.txt"),
}
DOCKER_DIR = Path(__file__).resolve().parents[2] / "docker"
# Image name -> the platform it is pulled and run for. An image absent here is pinned by a
# multi-platform index digest and is pulled and run for the host's platform.
PLATFORMS: dict[str, str] = {
    name: PLATFORM
    for name in (
        "deepvariant",
        "verifybamid2",
        "mutserve",
        "whatshap",
        "optitype",
        "t1k",
        "ucsc_bigbedtobed",
    )
} | {name: "linux/arm64" for name in IMAGES if name.startswith("pgsc_")}
# SpliceAI's TensorFlow is published for linux/arm64 as tensorflow-cpu-aws, which its list
# pins; TensorFlow's linux/amd64 wheels need AVX, which Rosetta does not emulate.
PLATFORMS["spliceai"] = "linux/arm64"
# The conda images are built from linux-aarch64 lockfiles.
PLATFORMS |= dict.fromkeys(
    (
        "smncaller",
        "parascopy",
        "yleaf",
        "expansionhunter",
        "reviewer",
        "delly",
        "dysgu",
        "truvari",
        "cnvpytor",
        "insurveyor",
        "scramble",
        "gauchian",
        "kilda",
        "haplocheck",
        "hmmix",
    ),
    "linux/arm64",
)
# RBCeq2's and fastmixture's wheels (numpy, pandas, pyarrow, polars) are pinned for
# linux/arm64, and fastmixture's extensions are compiled for the build machine's CPU.
PLATFORMS["rbceq2"] = PLATFORMS["fastmixture"] = "linux/arm64"


class DockerUnavailable(RuntimeError):
    """Docker is not installed or not running, or could not start the container (an image
    missing, among others); nothing was run."""


# `docker run`'s exit status when Docker itself failed and the container never ran (an image
# missing, among others). It refuses the run rather than leaving a gene unassessed; there is
# no adversary: a run without the bootstrap's images is incomplete, and the message names the
# bootstrap. A tool that itself exits 125 is refused the same way.
DOCKER_RUN_FAILED = 125
# Seconds `docker info` may take: a hung Docker Desktop answers never.
DOCKER_INFO_TIMEOUT = 60


def _docker_info(*args: str) -> str:
    """`docker info` with `args`, its stdout; `DockerUnavailable` when it fails."""
    docker = shutil.which("docker")
    if docker is None:
        raise DockerUnavailable("docker is not on PATH: install Docker Desktop")
    try:
        info = subprocess.run(
            [docker, "info", *args],
            check=False,
            capture_output=True,
            text=True,
            stdin=subprocess.DEVNULL,
            timeout=DOCKER_INFO_TIMEOUT,
        )
    except subprocess.TimeoutExpired as e:
        raise DockerUnavailable(
            "Docker did not answer: restart Docker Desktop and rerun"
        ) from e
    if info.returncode != 0:
        said = (info.stderr.strip().splitlines() or ["docker info failed"])[-1]
        raise DockerUnavailable(
            f"Docker is not running: start Docker Desktop and rerun ({said})"
        )
    return info.stdout


def image_missing(name: str) -> bool:
    """Whether Docker runs but has no pinned image `name` (a lab bootstrapped before the
    image was added); False when Docker itself is not running."""
    try:
        _docker_info()
    except DockerUnavailable:
        return False
    inspect = ["docker", "image", "inspect", IMAGES[name]]
    return subprocess.run(inspect, check=False, capture_output=True).returncode != 0


def docker_memory() -> int:
    """The bytes of memory Docker's VM gives its containers (`docker info` MemTotal)."""
    return int(_docker_info("--format", "{{.MemTotal}}"))


def run_container(
    name: str,
    args: Sequence[str],
    *,
    mounts: Sequence[tuple[Path, str, bool]],
    platform: str | None = None,
    label: str | None = None,
) -> None:
    """Run `args` in the pinned image `name`, offline, with `mounts` as (host path, container
    path, read-only) bind mounts; behind the core gate as `label` when given.

    The image is never pulled here (`--pull never`): the bootstrap pulls and checks it.
    Raises `DockerUnavailable` when `docker info` fails or Docker cannot start the
    container (`DOCKER_RUN_FAILED`), `ToolFailed` (with the container's stderr) when the
    tool fails, and `ValueError` for a mount path Docker would misread."""
    _docker_info()
    command = ["docker", "run", "--rm", "--network", "none", "--pull", "never"]
    if platform is not None:
        command += ["--platform", platform]
    for host, target, read_only in mounts:
        source = str(host.resolve())
        if "," in source or "," in target:
            raise ValueError(f"a mount path with a comma: {source} -> {target}")
        option = f"type=bind,source={source},target={target}"
        command += ["--mount", f"{option},readonly" if read_only else option]
    command += [IMAGES[name], *args]
    try:
        run_tool(gated(label, command) if label is not None else command)
    except ToolFailed as e:
        if e.returncode != DOCKER_RUN_FAILED:
            raise
        # Docker's own error line, which a gate's notice may precede; else the first line.
        lines = [s for s in (e.stderr or "").splitlines() if s.strip()]
        said = next(
            (s for s in lines if s.startswith("docker:")), lines[0] if lines else ""
        )
        raise DockerUnavailable(
            f"Docker could not start {IMAGES[name]} ({said.strip()}):"
            " run scripts/bootstrap.sh"
        ) from e


def _installed(file: str, data: bytes) -> bytes:
    """`data` as it decides what the build installs: a lock file (`.txt`, `.lock`) without
    its comment lines and blank lines, which install nothing; any other file whole (a
    Dockerfile's `# syntax=` line matters, a copied configuration file is installed as it
    is). A `#` inside a line stays: a conda lock carries each package's md5 after one."""
    if not file.endswith((".txt", ".lock")):
        return data
    lines = data.splitlines(keepends=True)
    return b"".join(x for x in lines if x.strip() and not x.lstrip().startswith(b"#"))


def build_inputs_digest(name: str) -> str:
    """The sha256 over a `BUILT` image's input files, each file's name and its bytes in
    order (a lock file's without comment lines, `_installed`), so an edit to what the build
    installs names another image and an edit to a lock file's comments does not."""
    digest = hashlib.sha256()
    for file in BUILT[name]:
        digest.update(file.encode() + b"\0")
        digest.update(_installed(file, (DOCKER_DIR / file).read_bytes()) + b"\0")
    return digest.hexdigest()


def build_command(name: str) -> list[str]:
    """The `docker build` of a `BUILT` image, tagged with its `IMAGES` name, for its
    `PLATFORMS` platform where it has one, else the host's (every base image is a
    multi-platform index pinned by digest)."""
    platform = ["--platform", PLATFORMS[name]] if name in PLATFORMS else []
    return [
        "docker", "build", *platform, "--tag", IMAGES[name],
        "--file", str(DOCKER_DIR / BUILT[name][0]), str(DOCKER_DIR),
    ]  # fmt: skip


def main(argv: Sequence[str] | None = None) -> int:
    """`python -m genokeeper.containers build NAME`: build a `BUILT` image (the bootstrap's
    step_images runs it for an image it lacks)."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2 or args[0] != "build" or args[1] not in BUILT:
        print(
            f"usage: python -m genokeeper.containers build {{{','.join(BUILT)}}}",
            file=sys.stderr,
        )
        return 2
    if build_inputs_digest(args[1]) != IMAGES[args[1]].rsplit("-", 1)[1]:
        print(
            f"docker/ holds other inputs than {IMAGES[args[1]]} was named for:"
            " restore them, or rename the image in genokeeper.containers.IMAGES",
            file=sys.stderr,
        )
        return 1
    return subprocess.run(build_command(args[1]), check=False).returncode


if __name__ == "__main__":
    sys.exit(main())
