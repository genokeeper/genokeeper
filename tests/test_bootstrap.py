"""scripts/bootstrap.sh, step by step, with stub `curl`, `brew` and `docker` and nothing
downloaded.

Sourcing the script defines its functions without running a step; each test calls the steps it
covers in a fresh bash.
"""

import gzip
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
import tarfile
import warnings
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pysam
import pytest
from fixtures.clingen_mini import (
    write_dosage,
    write_gene_validity,
    write_recurrent,
    write_regions,
)
from fixtures.clinvar_mini import write_clinvar
from fixtures.faf_mini import write_faf_table
from fixtures.gwas_mini import Assoc, write_gwas
from fixtures.medlineplus_mini import gene as medlineplus_gene
from fixtures.medlineplus_mini import write_medlineplus
from fixtures.pharmgkb_mini import Ann, zip_bytes

from genokeeper.annotate.vep import vep_files
from genokeeper.calibrate import BENCH_BED, TRUTH_VCF, layout, rtg_path
from genokeeper.calibrate.hg002 import NOT_YET
from genokeeper.config import grch38_fasta, load_grading, unmasked_fasta
from genokeeper.containers import BUILT, IMAGES, PLATFORM, PLATFORMS, build_command
from genokeeper.deletions.files import sv_files
from genokeeper.evidence.clingen import DOSAGE_FILE, RECURRENT_FILE, REGIONS_FILE
from genokeeper.evidence.frequency import GNOMAD_CHROMS, faf_folder
from genokeeper.evidence.gwas import DOWNLOAD as GWAS_DOWNLOAD
from genokeeper.evidence.medlineplus import MEDLINEPLUS_FILE, MEDLINEPLUS_URL
from genokeeper.evidence.pgx import CPIC_FILES as CPIC_FILE_NAMES
from genokeeper.masking import MASKED_SHA256, masked_bed
from genokeeper.pgx.hla import IMGT_HLA_VERSION, t1k_index
from genokeeper.pgx.pharmcat import PHARMCAT_VERSION, positions_path
from genokeeper.pgx.stellarpgx import STELLARPGX_COMMIT, stellarpgx_dir
from genokeeper.prs.pgsc import (
    NEXTFLOW_VERSION,
    NF_SCHEMA_VERSION,
    PANEL_NAME,
    PANEL_PVAR_NAME,
    PGSC_CALC_COMMIT,
    PGSC_CALC_VERSION,
    ancestry_pvar,
    java_home,
    nextflow_path,
    nxf_home,
    panel_path,
    panel_pvar,
    pgsc_calc_dir,
    scoring_file,
)
from genokeeper.prs.scores import load_research_scores, load_scores
from genokeeper.regions import (
    AUTOSOMES,
    CHRX_AMPLICONIC,
    CHRX_NONPAR,
    CHRX_XTR,
    CHRY_NONPAR,
    DERIVED,
    SEX_REGIONS,
    giab38_dir,
    giab38_downloads,
)
from genokeeper.sources.panelapp import PANEL_VERSION, panelapp_path
from genokeeper.typing_files import ARCHAIC, HMMIX_HELPERS, TypingFiles, typing_files

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "bootstrap.sh"
LITVAR_FIXTURE = REPO / "tests" / "fixtures" / "litvar_head.jsonl.gz"
CITATIONS_FIXTURE = REPO / "tests" / "fixtures" / "var_citations_head.txt"
PANELAPP_URL = "https://panelapp.genomicsengland.co.uk/api/v1/panels/285/?format=json"
# Copies $CLINGEN_CSV to the -o file and logs the URL (the last argument).
FAKE_CLINGEN_CURL = """#!/usr/bin/env bash
out=
for a; do [ "${prev:-}" = -o ] && out=$a; prev=$a; done
echo "$prev" >> "$CURL_LOG"
cp "$CLINGEN_CSV" "$out"
"""
CLINGEN_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"


def _executable(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


@pytest.fixture
def lab(tmp_path: Path) -> dict[str, str]:
    """Environment for a lab under tmp_path."""
    home = tmp_path / "lab"
    home.mkdir()
    # No test may download or install: curl and brew fail unless a test stubs them.
    offline = tmp_path / "offline-bin"
    for tool in ("curl", "brew"):
        _executable(offline / tool, "#!/bin/sh\necho offline >&2\nexit 22\n")
    return {
        **os.environ,
        "PATH": f"{offline}:{os.environ['PATH']}",
        "GENOKEEPER_HOME": str(home),
    }


def _run(env: dict[str, str], steps: str) -> subprocess.CompletedProcess[str]:
    # No stdin, as under CI: a tool that asks (unzip over a duplicate entry) fails here.
    return subprocess.run(
        ["bash", "-c", f'source "$0"; {steps}', str(SCRIPT)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
        stdin=subprocess.DEVNULL,
    )


def test_the_lab_falls_back_to_its_variable_before_the_rename(
    lab: dict[str, str], tmp_path: Path
) -> None:
    env = {k: v for k, v in lab.items() if k != "GENOKEEPER_HOME"}
    done = _run(env | {"GENOMELAB_HOME": str(tmp_path / "old")}, 'echo "$lab"')
    assert done.stdout == f"{tmp_path / 'old'}\n", done.stderr
    both = env | {"GENOMELAB_HOME": "/old", "GENOKEEPER_HOME": str(tmp_path / "new")}
    assert _run(both, 'echo "$lab"').stdout == f"{tmp_path / 'new'}\n"


def test_check_without_a_venv_creates_none(lab: dict[str, str], tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "config").mkdir()
    shutil.copy(SCRIPT, repo / "scripts" / SCRIPT.name)
    shutil.copy(REPO / "config" / "grading.toml", repo / "config" / "grading.toml")
    (repo / "pyproject.toml").write_text('[project]\nname = "x"\nversion = "0"\n')
    result = subprocess.run(
        ["bash", str(repo / "scripts" / SCRIPT.name), "--check"],
        env=lab,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert "missing: project environment for genokeeper.regions (uv sync)" in lines
    assert "missing: project environment for the GRCh38 FASTA path (uv sync)" in lines
    assert not (repo / ".venv").exists()


def test_check_lists_missing_gene_validity(lab: dict[str, str]) -> None:
    result = _run(lab, "check=1; step_clingen")
    assert (result.returncode, result.stdout) == (
        0,
        f"missing: ClinGen gene-disease validity export ({CLINGEN_URL})\n",
    )


def _clingen_env(lab: dict[str, str], tmp_path: Path, csv: Path) -> dict[str, str]:
    bin_dir = tmp_path / "bin"
    _executable(bin_dir / "curl", FAKE_CLINGEN_CURL)
    return {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "CLINGEN_CSV": str(csv),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }


def test_gene_validity_is_filed_under_its_own_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    csv = write_gene_validity(tmp_path / "served.csv", [("BRCA1", "Definitive")])
    env = _clingen_env(lab, tmp_path, csv)
    home = Path(lab["GENOKEEPER_HOME"])
    # A partial download of an earlier day's file is never resumed into today's.
    (home / "clingen" / ".download").mkdir(parents=True)
    (home / "clingen" / ".download" / "gene-validity.csv.part").write_text("stale")
    done = _run(env, "step_clingen")
    assert done.returncode == 0, done.stderr
    dest = home / "clingen" / "20261018" / "gene-validity.csv"
    assert dest.read_bytes() == csv.read_bytes()
    assert not (home / "clingen" / ".download").exists()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert [row.split("\t")[:2] for row in manifest] == [["clingen", str(dest)]]
    assert Path(env["CURL_LOG"]).read_text().splitlines() == [CLINGEN_URL]
    again = _run(env, "step_clingen; check=1; step_clingen")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 1


def test_gene_validity_that_does_not_parse_is_not_recorded(
    lab: dict[str, str], tmp_path: Path
) -> None:
    served = tmp_path / "served.csv"
    served.write_text("<html>maintenance</html>\n")
    result = _run(_clingen_env(lab, tmp_path, served), "step_clingen")
    assert result.returncode == 1
    assert "ClinGen gene validity" in result.stderr
    home = Path(lab["GENOKEEPER_HOME"])
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not list((home / "clingen").glob("2*"))


def test_a_rebuilt_clingen_file_is_filed_only_when_its_curations_change(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """ClinGen writes the request day into the validity export and the build day into the
    dosage list: a refresh files either again only when a line other than that date
    differs, and otherwise names the newest folder, so it reads "unchanged"."""
    home = Path(lab["GENOKEEPER_HOME"])
    csv = tmp_path / "served.csv"
    env = _clingen_env(lab, tmp_path, csv)
    dosage_env, serve = _serving(lab, tmp_path)
    cases = (
        (
            "clingen",
            env,
            lambda day, gene: write_gene_validity(
                csv, [(gene, "Definitive")], created=f"2026-10-{day}"
            ),
        ),
        (
            "clingen_dosage",
            dosage_env,
            lambda day, gene: write_dosage(
                serve / DOSAGE_FILE, [(gene, "3")], date=f"{day} Oct,2026"
            ),
        ),
    )
    for step, step_env, serve_file in cases:
        folder = home / step.replace("_", "-")
        serve_file("18", "BRCA1")
        assert _refresh(step_env, step).returncode == 0
        serve_file("27", "BRCA1")  # rebuilt on a later day, nothing else changed
        same = _refresh(step_env, step)
        assert (same.returncode, same.stdout) == (
            0,
            f"refreshed: {step} {folder / '20261018'}\n",
        )
        assert _folders(folder) == ["20261018"]
        serve_file("27", "BRCA2")
        new = _refresh(step_env, step)
        assert new.stdout == f"refreshed: {step} {folder / '20261027'}\n"
        assert _folders(folder) == ["20261018", "20261027"]


def test_brew_installs_without_analytics(lab: dict[str, str], tmp_path: Path) -> None:
    stubs = tmp_path / "brew-bin"
    _executable(
        stubs / "brew",
        '#!/bin/sh\necho "${HOMEBREW_NO_ANALYTICS:-unset} $*" >> "$BREW_LOG"\n',
    )
    # No bcftools, samtools or htslib on this PATH: brew is asked for them. The system
    # directories are linked in without those tools, which a CI runner may have installed.
    system = tmp_path / "system-bin"
    system.mkdir()
    for directory in (Path("/usr/bin"), Path("/bin")):
        for tool in directory.iterdir():
            link = system / tool.name
            if tool.name not in {"bcftools", "samtools", "bgzip", "tabix"} and not (
                link.exists() or link.is_symlink()
            ):
                link.symlink_to(tool)
    env = {**lab, "PATH": f"{stubs}:{system}", "BREW_LOG": str(tmp_path / "b.log")}
    _run(env, "step_tools")
    assert Path(env["BREW_LOG"]).read_text().splitlines() == [
        "1 install bcftools samtools htslib"
    ]


# --- flags, the GRCh38 FASTA, the gnomAD FAF table ----------------------------------------


def test_flags_combine(lab: dict[str, str]) -> None:
    for flags, expected in (
        ("", "0 0 0"),
        ("--check", "1 0 0"),
        ("--with-gnomad", "0 1 0"),
        ("--check --with-gnomad", "1 1 0"),
        ("--with-gnomad --check", "1 1 0"),
        ("--with-calibration", "0 0 1"),
        ("--check --with-calibration --with-gnomad", "1 1 1"),
    ):
        result = _run(
            lab, f'parse_args {flags}; echo "$check $with_gnomad $with_calibration"'
        )
        assert (result.returncode, result.stdout) == (0, f"{expected}\n"), flags


def test_unknown_flag_exits_64(lab: dict[str, str]) -> None:
    for flags in (["--bogus"], ["--check", "--with-gnomed"]):
        result = subprocess.run(
            ["bash", str(SCRIPT), *flags],
            env=lab,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 64, flags
        assert "[--check] [--with-gnomad] [--with-calibration]" in result.stderr
        # Refused before any step ran.
        assert not (Path(lab["GENOKEEPER_HOME"]) / "bootstrap-manifest.tsv").exists()


# `genokeeper build-faf` without --check is logged, never run: it would stream 0.88 TB.
# Everything else (the FASTA path, `build-faf --check`) goes to the real uv.
FAKE_GNOMAD_UV = """#!/usr/bin/env bash
case " $* " in
  *" build-faf --check "*) ;;
  *" build-faf "*) echo "$*" >> "$UV_LOG"; exit 0 ;;
esac
exec "$REAL_UV" "$@"
"""
GNOMAD_MESSAGE = (
    "gnomAD FAF table missing: run with --with-gnomad (streams 0.88 TB, keeps a few GB)"
)


def _uv_calls(env: dict[str, str]) -> list[str]:
    log = Path(env["UV_LOG"])
    return log.read_text().splitlines() if log.exists() else []


@pytest.fixture
def gnomad(lab: dict[str, str], tmp_path: Path) -> dict[str, str]:
    """`lab` with a `uv` that logs `genokeeper build-faf` instead of running it."""
    real = shutil.which("uv")
    assert real is not None
    stubs = tmp_path / "gnomad-bin"
    _executable(stubs / "uv", FAKE_GNOMAD_UV)
    return {
        **lab,
        "PATH": f"{stubs}:{lab['PATH']}",
        "REAL_UV": real,
        "UV_LOG": str(tmp_path / "uv.log"),
    }


def test_gnomad_step_needs_flag(gnomad: dict[str, str]) -> None:
    without = _run(gnomad, "step_gnomad_faf")
    assert without.returncode == 0, without.stderr
    assert f"bootstrap: {GNOMAD_MESSAGE}" in without.stderr.splitlines()
    assert _uv_calls(gnomad) == []
    built = _run(gnomad, "with_gnomad=1; step_gnomad_faf")
    assert built.returncode == 0, built.stderr
    assert _uv_calls(gnomad) == [
        f"run --project {REPO} --frozen --quiet genokeeper build-faf"
    ]
    # A complete table: nothing to say without the flag.
    write_faf_table(faf_folder(Path(gnomad["GENOKEEPER_HOME"])))
    complete = _run(gnomad, "step_gnomad_faf")
    assert complete.returncode == 0, complete.stderr
    assert GNOMAD_MESSAGE not in complete.stderr


def test_check_lists_missing_faf(gnomad: dict[str, str]) -> None:
    folder = faf_folder(Path(gnomad["GENOKEEPER_HOME"]))
    write_faf_table(folder, chroms=GNOMAD_CHROMS[:-2])
    for flags in ("check=1", "check=1; with_gnomad=1"):
        result = _run(gnomad, f"{flags}; step_gnomad_faf")
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == [
            f"missing: gnomAD FAF table {c} (run with --with-gnomad)"
            for c in ("chrX", "chrY")
        ]
    assert _uv_calls(gnomad) == []  # --check streams nothing, even with the flag


# Serves $SERVE/<the URL's file name> to the -o file and logs the URL (the last argument).
FAKE_SERVING_CURL = """#!/usr/bin/env bash
out=
for a; do [ "${prev:-}" = -o ] && out=$a; prev=$a; done
echo "$prev" >> "$CURL_LOG"
[ -f "$SERVE/${prev##*/}" ] || exit 22
cp "$SERVE/${prev##*/}" "$out"
"""
GRCH38_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
    "seqs_for_alignment_pipelines.ucsc_ids"
)
GRCH38_NAME = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"


@pytest.mark.tools
def test_grch38_fasta_md5_checked(lab: dict[str, str], tmp_path: Path) -> None:
    serve = tmp_path / "serve"
    serve.mkdir()
    fasta_text = ">chr1\nACGTACGTAC\n>chr2\nGGGGCCCCAA\n"
    packed = gzip.compress(fasta_text.encode())
    (serve / GRCH38_NAME).write_bytes(packed)
    sums = serve / "md5checksums.txt"
    sums.write_text(f"{'0' * 32}  ./{GRCH38_NAME}\n")
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_SERVING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    home = Path(lab["GENOKEEPER_HOME"])
    dest = unmasked_fasta(home)
    manifest = home / "bootstrap-manifest.tsv"
    refused = _run(env, "step_grch38_fasta")
    assert refused.returncode == 1
    assert f"md5 mismatch for {GRCH38_BASE}/{GRCH38_NAME}" in refused.stderr
    assert not dest.exists() and not manifest.exists()
    sums.write_text(
        f"{'1' * 32}  ./other.fna.gz\n"
        f"{hashlib.md5(packed).hexdigest()}  ./{GRCH38_NAME}\n"
    )
    done = _run(env, "step_grch38_fasta")
    assert done.returncode == 0, done.stderr
    assert dest.read_text() == fasta_text
    assert Path(f"{dest}.fai").is_file()
    assert [row.split("\t")[:2] for row in manifest.read_text().splitlines()] == [
        ["grch38_fasta", str(dest)],
        ["grch38_fasta", f"{dest}.fai"],
    ]
    assert not (dest.parent / ".grch38").exists()
    fetched = Path(env["CURL_LOG"]).read_text().splitlines()
    assert set(fetched) == {
        f"{GRCH38_BASE}/md5checksums.txt",
        f"{GRCH38_BASE}/{GRCH38_NAME}",
    }
    again = _run(env, "step_grch38_fasta; check=1; step_grch38_fasta")
    assert (again.returncode, again.stdout) == (0, "")
    assert Path(env["CURL_LOG"]).read_text().splitlines() == fetched


CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly"


@pytest.mark.tools
def test_clinvar_grch38_release_filed_under_its_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The weekly GRCh38 VCF the .md5 names, md5-checked, filed under its ##fileDate; no
    GRCh37 file is asked for."""
    serve = tmp_path / "serve"
    serve.mkdir()
    name = "clinvar_20261011.vcf.gz"
    record = ("1", 100, 1, "A", "G", {"CLNSIG": "Pathogenic"})
    vcf = write_clinvar(serve / name, "GRCh38", [record])
    (serve / f"{name}.tbi").write_bytes(Path(f"{vcf}.tbi").read_bytes())
    md5 = serve / "clinvar.vcf.gz.md5"
    md5.write_text(f"{'0' * 32}  /netmnt/vcf_GRCh38/{name}\n")
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_SERVING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    home = Path(lab["GENOKEEPER_HOME"])
    manifest = home / "bootstrap-manifest.tsv"
    refused = _run(env, "step_clinvar")
    assert refused.returncode == 1
    assert f"md5 mismatch for {CLINVAR_BASE}/{name}" in refused.stderr
    assert not manifest.exists() and not list((home / "clinvar").glob("2*"))
    md5.write_text(f"{hashlib.md5(vcf.read_bytes()).hexdigest()}  /x/{name}\n")
    done = _run(env, "step_clinvar")
    assert done.returncode == 0, done.stderr
    dest = home / "clinvar" / "20261012"  # the ##fileDate write_clinvar writes
    assert sorted(p.name for p in dest.iterdir()) == [
        "clinvar_GRCh38.vcf.gz",
        "clinvar_GRCh38.vcf.gz.tbi",
    ]
    assert [row.split("\t")[:2] for row in manifest.read_text().splitlines()] == [
        ["clinvar", str(dest / "clinvar_GRCh38.vcf.gz")],
        ["clinvar", str(dest / "clinvar_GRCh38.vcf.gz.tbi")],
    ]
    fetched = Path(env["CURL_LOG"]).read_text().splitlines()
    assert all(url.startswith(f"{CLINVAR_BASE}/") for url in fetched)
    again = _run(env, "step_clinvar; check=1; step_clinvar")
    assert (again.returncode, again.stdout) == (0, "")
    assert Path(env["CURL_LOG"]).read_text().splitlines() == fetched


@pytest.mark.tools
def test_clinvar_part_of_an_earlier_week_is_not_resumed(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The weekly link names another file each week: a part cut from last week's release
    is never resumed onto this week's."""
    serve = tmp_path / "serve"
    serve.mkdir()
    name = "clinvar_20261011.vcf.gz"
    record = ("1", 100, 1, "A", "G", {"CLNSIG": "Pathogenic"})
    vcf = write_clinvar(serve / name, "GRCh38", [record])
    (serve / f"{name}.tbi").write_bytes(Path(f"{vcf}.tbi").read_bytes())
    md5 = hashlib.md5(vcf.read_bytes()).hexdigest()
    (serve / "clinvar.vcf.gz.md5").write_text(f"{md5}  /x/{name}\n")
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    work = Path(lab["GENOKEEPER_HOME"]) / "clinvar" / ".download"
    work.mkdir(parents=True)
    (work / "GRCh38.md5.part").write_bytes(b"0123")
    (work / "clinvar_20261004.vcf.gz.part").write_bytes(b"last week's first bytes")
    done = _run(env, "step_clinvar")
    assert done.returncode == 0, done.stderr


# --- stage 2 inputs: tools, GIAB GRCh38, the checkpoint, CGD, images -----------------------

DV_BASE = "https://storage.googleapis.com/deepvariant/complete-case-study-testdata/complete-t7"
DV_FILES = (
    "example_info.json",
    "model.example_info.json",
    "weights-51-0.995354.ckpt.data-00000-of-00001",
    "weights-51-0.995354.ckpt.index",
)
CGD_URL = "https://research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz"
RTG_URL = (
    "https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.13/"
    "rtg-tools-3.13-nojre.zip"
)


def _path_without(tmp_path: Path, tools: set[str], path: str) -> str:
    """A PATH of one directory linking every executable on `path`, the first one of each
    name, except `tools`: those are absent whatever this machine has installed."""
    linked = tmp_path / "path-without"
    linked.mkdir(exist_ok=True)
    for entry in path.split(os.pathsep):
        directory = Path(entry)
        if not directory.is_dir():
            continue
        for tool in directory.iterdir():
            link = linked / tool.name
            if tool.name in tools or link.exists() or link.is_symlink():
                continue
            if tool.is_file() and os.access(tool, os.X_OK):
                link.symlink_to(tool)
    return str(linked)


# Logs its arguments; answers `info` unless $DOCKER_DOWN is set, prints $DOCKER_DIGESTS for
# `image inspect` (failing when it is empty) and succeeds on `pull` and `build`.
FAKE_DOCKER = """#!/usr/bin/env bash
echo "$*" >> "$DOCKER_LOG"
case "$1" in
  info) [ -z "${DOCKER_DOWN:-}" ] ;;
  image) [ -n "${DOCKER_DIGESTS:-}" ] && printf '%s\\n' $DOCKER_DIGESTS ;;
  pull) exit 0 ;;
  build) exit 0 ;;
  *) exit 1 ;;
esac
"""


@pytest.fixture
def stage2(lab: dict[str, str], tmp_path: Path) -> dict[str, str]:
    """`lab` without bwa, mosdepth or rtg on PATH, and with a stub `docker`."""
    stubs = tmp_path / "docker-bin"
    _executable(stubs / "docker", FAKE_DOCKER)
    path = _path_without(tmp_path, {"bwa", "mosdepth", "rtg", "docker"}, lab["PATH"])
    return {
        **lab,
        "PATH": f"{stubs}:{path}",
        "DOCKER_LOG": str(tmp_path / "docker.log"),
        "DOCKER_DIGESTS": "",
    }


def test_check_lists_stage2_inputs(stage2: dict[str, str]) -> None:
    home = Path(stage2["GENOKEEPER_HOME"])
    steps = "step_tools38; step_giab38; step_dv_model; step_cgd"
    result = _run(stage2, f"check=1; {steps}; step_images")
    assert result.returncode == 0, result.stderr
    giab = giab38_dir(home)
    assert result.stdout.splitlines() == [
        "missing: tool bwa 0.7.19 (brew install bwa)",
        "missing: tool mosdepth 0.3.14 (brew install mosdepth)",
        (
            f"missing: rtg-tools 3.13 in {home / 'tools' / 'rtg-tools-3.13'}"
            f" (Homebrew's openjdk and {RTG_URL})"
        ),
        *(
            f"missing: GIAB v3.6 GRCh38 stratification {name}"
            for name in giab38_downloads(load_grading())
        ),
        *(f"missing: region BED {giab / name} (derived)" for name in DERIVED),
        *(
            f"missing: DeepVariant T7 checkpoint file {name} ({DV_BASE}/{name})"
            for name in DV_FILES
        ),
        f"missing: CGD gene file ({CGD_URL})",
        *(f"missing: container image {name} ({ref})" for name, ref in IMAGES.items()),
    ]
    assert "pull" not in Path(stage2["DOCKER_LOG"]).read_text()
    assert not (home / "bootstrap-manifest.tsv").exists()


GIAB38_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "genome-stratifications/v3.6/GRCh38@all"
)
GIAB38_MD5S = "GRCh38-genome-stratifications-md5s.txt"
# Each source's intervals: enough for every derived region to be non-empty.
GIAB38_ROWS = {
    CHRY_NONPAR: [("chrY", 1000, 9000)],
    CHRX_NONPAR: [("chrX", 1000, 9000)],
    AUTOSOMES: [("chr1", 0, 9000)],
}


def _bgzipped_bed(path: Path, rows: list[tuple[str, int, int]]) -> Path:
    plain = path.with_name(f"{path.name}.plain")
    plain.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    plain.unlink()
    return path


@pytest.mark.tools
def test_giab38_fetches_only_missing(lab: dict[str, str], tmp_path: Path) -> None:
    """A lab whose GIAB files predate slice 3c's three chrX BEDs fetches those alone."""
    names = giab38_downloads(load_grading())
    added = {CHRX_NONPAR, CHRX_XTR, CHRX_AMPLICONIC}
    serve = tmp_path / "serve"
    serve.mkdir()
    md5s = []
    for name in names:
        bed = _bgzipped_bed(serve / name, GIAB38_ROWS.get(name, [("chr1", 10, 20)]))
        md5s.append(f"{hashlib.md5(bed.read_bytes()).hexdigest()}  ./XY/{name}\n")
    (serve / GIAB38_MD5S).write_text("".join(md5s))
    giab = giab38_dir(Path(lab["GENOKEEPER_HOME"]))
    giab.mkdir(parents=True)
    recorded = []
    for name in sorted(set(names) - added):
        shutil.copy(serve / name, giab / name)
        pysam.tabix_index(str(giab / name), preset="bed", force=True)
        recorded += [str(giab / name), f"{giab / name}.tbi"]
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_SERVING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    old = _run(env, f"record giab38 {' '.join(recorded)}")
    assert old.returncode == 0, old.stderr
    done = _run(env, "step_giab38")
    assert done.returncode == 0, done.stderr
    fetched = Path(env["CURL_LOG"]).read_text().splitlines()
    assert sorted(url.rsplit("/", 1)[1] for url in fetched) == sorted(
        {GIAB38_MD5S, *added}
    )
    assert "xnonpar" in (giab / SEX_REGIONS).read_text()
    again = _run(env, "check=1; step_giab38")
    assert (again.returncode, again.stdout) == (0, "")


def test_check_names_docker_when_it_is_not_running(stage2: dict[str, str]) -> None:
    down = {**stage2, "DOCKER_DOWN": "1"}
    listed = _run(down, "check=1; step_images")
    assert listed.returncode == 0, listed.stderr
    assert listed.stdout.startswith("missing: Docker (start Docker Desktop")
    refused = _run(down, "step_images")
    assert refused.returncode == 1
    assert "Docker is not running" in refused.stderr


@pytest.mark.parametrize(
    ("bwa", "mosdepth", "warned"),
    [("0.7.17", "0.3.14", "bwa 0.7.17"), ("0.7.19", "0.3.15", "mosdepth 0.3.15")],
)
def test_bwa_version_pinned(
    stage2: dict[str, str], tmp_path: Path, bwa: str, mosdepth: str, warned: str
) -> None:
    """Another bwa or mosdepth version (the user after `brew upgrade`) is a warning, not a
    refusal: the run records the version and keys its caches on it."""
    stubs = tmp_path / "bwa-bin"
    # bwa prints its version to stderr and exits 1 when run without a command.
    _executable(
        stubs / "bwa",
        f'#!/bin/sh\necho "Program: bwa" >&2\necho "Version: {bwa}-r1188" >&2\nexit 1\n',
    )
    _executable(stubs / "mosdepth", f'#!/bin/sh\necho "mosdepth {mosdepth}"\n')
    env = {**stage2, "PATH": f"{stubs}:{stage2['PATH']}"}
    for steps in ("step_tools38", "check=1; step_tools38"):
        # rtg is stubbed out: only the two tools' lines are under test.
        done = _run(env, f"step_rtg() {{ :; }}; {steps}")
        assert done.returncode == 0, done.stderr
        assert f"warning: {warned} is not the tested " in done.stderr
        assert "missing: tool" not in done.stdout
        assert "offline" not in done.stderr  # nothing installed over it


def _grch38_fai(fasta: Path, contigs: list[tuple[str, int]]) -> Path:
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text("".join(f">{name}\nN\n" for name, _ in contigs))
    Path(f"{fasta}.fai").write_text(
        "".join(f"{name}\t{length}\t0\t1\t2\n" for name, length in contigs)
    )
    return fasta


CONTIGS = [("chr1", 248956422), ("chrY", 57227415)]


def test_checkpoint_on_disk_is_recognised_and_md5_checked(
    lab: dict[str, str], tmp_path: Path
) -> None:
    home = Path(lab["GENOKEEPER_HOME"])
    folder = home / "deepvariant" / "complete-t7"
    folder.mkdir(parents=True)
    contents = {name: f"{name}\n".encode() for name in DV_FILES}
    for name, data in contents.items():
        (folder / name).write_bytes(data)
    md5s = " ".join(f"{n}:{hashlib.md5(d).hexdigest()}" for n, d in contents.items())
    pins = f"dv_model_files='{md5s}'"
    done = _run(lab, f"{pins}; step_dv_model")  # curl is offline: nothing is fetched
    assert done.returncode == 0, done.stderr
    manifest = home / "bootstrap-manifest.tsv"
    assert [row.split("\t")[:2] for row in manifest.read_text().splitlines()] == [
        ["dv_model", str(folder / name)] for name in DV_FILES
    ]
    # A file whose md5 differs is fetched again; a served copy that still differs is
    # refused and not kept.
    (folder / DV_FILES[0]).write_bytes(b"changed\n")
    manifest.unlink()
    serve = tmp_path / "serve"
    serve.mkdir()
    (serve / DV_FILES[0]).write_bytes(b"still not it\n")
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_SERVING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    refused = _run(env, f"{pins}; step_dv_model")
    assert refused.returncode == 1
    assert f"md5 mismatch for {DV_BASE}/{DV_FILES[0]}" in refused.stderr
    fetched = Path(env["CURL_LOG"]).read_text().splitlines()
    assert fetched == [f"{DV_BASE}/{DV_FILES[0]}"]
    assert not manifest.exists()
    assert not (folder / ".download").exists()


# Copies $SERVE/CGD.txt.gz to the -o file, dated as the server's Last-Modified when given
# -R, and logs its arguments.
FAKE_CGD_CURL = """#!/usr/bin/env bash
out=
for a; do [ "${prev:-}" = -o ] && out=$a; prev=$a; done
echo "$*" >> "$CURL_LOG"
cp "$SERVE/CGD.txt.gz" "$out"
case " $* " in *" -R "*) touch -t 202502211704 "$out" ;; esac
"""
CGD_HEADER = "#GENE\tHGNC ID\tENTREZ GENE ID\tCONDITION\tINHERITANCE\n"


def test_cgd_is_filed_under_its_server_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    serve = tmp_path / "serve"
    serve.mkdir()
    served = gzip.compress((CGD_HEADER + "A2M\t7\t2\tx\tAD\n").encode())
    (serve / "CGD.txt.gz").write_bytes(served)
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_CGD_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    done = _run(env, "step_cgd")
    assert done.returncode == 0, done.stderr
    home = Path(lab["GENOKEEPER_HOME"])
    dest = home / "cgd" / "20250221" / "CGD.txt.gz"
    assert dest.read_bytes() == served
    manifest = home / "bootstrap-manifest.tsv"
    sha = hashlib.sha256(served).hexdigest()
    assert manifest.read_text().splitlines() == [f"cgd\t{dest}\t{sha}"]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert len(calls) == 1 and calls[0].endswith(CGD_URL)
    assert " -C " not in f" {calls[0]} "  # fetched whole, never resumed
    again = _run(env, "step_cgd; check=1; step_cgd")
    assert (again.returncode, again.stdout) == (0, "")
    # A page that is not the CGD table is not recorded.
    (serve / "CGD.txt.gz").write_bytes(gzip.compress(b"<html>maintenance</html>\n"))
    manifest.unlink()
    refused = _run(env, "step_cgd")
    assert refused.returncode == 1
    assert "not the CGD table" in refused.stderr
    assert not manifest.exists()


def test_missing_image_pulled_by_its_digest(stage2: dict[str, str]) -> None:
    """Only the image whose pinned digest the store lacks is pulled, by that digest."""
    refs = list(IMAGES.values())
    present = {**stage2, "DOCKER_DIGESTS": " ".join(refs[1:])}
    done = _run(present, "step_images")
    assert done.returncode == 0, done.stderr
    calls = Path(stage2["DOCKER_LOG"]).read_text().splitlines()
    assert [line for line in calls if line.startswith("pull")] == [
        f"pull --platform {PLATFORM} {refs[0]}"
    ]
    # Every pinned digest present: nothing is pulled.
    Path(stage2["DOCKER_LOG"]).unlink()
    pinned = {**stage2, "DOCKER_DIGESTS": " ".join(IMAGES.values())}
    done = _run(pinned, "step_images; check=1; step_images")
    assert (done.returncode, done.stdout) == (0, "")
    calls = Path(stage2["DOCKER_LOG"]).read_text().splitlines()
    assert not [line for line in calls if line.startswith("pull")]


# --- the calibration inputs (--with-calibration) -------------------------------------------

HG003_FASTQ_BASE = (
    "https://ftp.cngb.org/pub/CNSA/data7/CNP0008319/CNS1414628/CNX1267209/CNR1441268"
)
HG003_TRUTH_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
    "HG003_NA24149_father/NISTv4.2.1/GRCh38"
)
# Serves $SERVE/<the URL's file name> to the -o file and logs every argument. With -C - it
# appends what the file lacks, as a resumed transfer does. While $CUT_MARK does not exist,
# a FASTQ transfer writes $CUT_BYTES bytes, creates the mark and fails as a cut connection.
FAKE_RESUMING_CURL = """#!/usr/bin/env bash
out=; resume=0
for a; do
  [ "${prev:-}" = -o ] && out=$a
  [ "${prev:-}" = -C ] && [ "$a" = - ] && resume=1
  prev=$a
done
echo "$*" >> "$CURL_LOG"
src=$SERVE/${prev##*/}
[ -f "$src" ] || exit 22
have=0
if [ "$resume" = 1 ] && [ -f "$out" ]; then have=$(wc -c < "$out" | tr -d ' '); fi
if [ -n "${CUT_MARK:-}" ] && [ ! -e "$CUT_MARK" ] && [[ "$src" == *.fq.gz ]]; then
  touch "$CUT_MARK"
  head -c "$CUT_BYTES" "$src" > "$out"
  exit 18
fi
if [ "$have" -gt 0 ]; then tail -c +$((have + 1)) "$src" >> "$out"; else cp "$src" "$out"; fi
"""
# rtg as the bootstrap installs it: `version` names 3.13; `format -o DIR FASTA` writes an
# SDF's `done` and `mainIndex` and logs the FASTA.
FAKE_RTG = """#!/usr/bin/env bash
case "$1" in
  version) echo "Product: RTG Tools 3.13" ;;
  format)
    shift
    while [ $# -gt 0 ]; do
      case "$1" in -o) out=$2; shift 2 ;; *) fasta=$1; shift ;; esac
    done
    mkdir -p "$out"
    echo "$fasta" > "$out/mainIndex"
    echo "SDF created" > "$out/done"
    echo "format $fasta" >> "$RTG_LOG"
    ;;
  *) exit 1 ;;
esac
"""


@dataclass(frozen=True)
class Calibration:
    env: dict[str, str]
    pins: (
        str  # shell assignments pinning the served files, and the FASTA's manifest row
    )
    serve: Path
    fastq: tuple[str, str]
    truth: tuple[str, str, str]

    @property
    def home(self) -> Path:
        return Path(self.env["GENOKEEPER_HOME"])

    def curl_calls(self) -> list[str]:
        log = Path(self.env["CURL_LOG"])
        return log.read_text().splitlines() if log.exists() else []


@pytest.fixture
def calibration(lab: dict[str, str], tmp_path: Path) -> Calibration:
    """A lab with the GRCh38 FASTA and rtg (stubs), a curl serving two small FASTQ and
    three truth files, and the shell pins naming their md5s."""
    home = Path(lab["GENOKEEPER_HOME"])
    _grch38_fai(grch38_fasta(home), CONTIGS)
    _executable(home / "tools" / "rtg-tools-3.13" / "rtg", FAKE_RTG)
    serve = tmp_path / "serve"
    serve.mkdir()
    fastq = ("reads_1.fq.gz", "reads_2.fq.gz")
    truth = ("truth.vcf.gz", "truth.vcf.gz.tbi", "bench.bed")
    pins_fastq: list[str] = []
    pins_truth: list[str] = []
    for i, name in enumerate(fastq):
        data = gzip.compress(f"@E1L1C001R00100000{i}/{i + 1}\nACGT\n+\nIIII\n".encode())
        (serve / name).write_bytes(data)
        pins_fastq.append(f"{name}:{hashlib.md5(data).hexdigest()}:{len(data)}")
    for name in truth:
        data = f"{name}\n".encode()
        (serve / name).write_bytes(data)
        pins_truth.append(f"{name}:{hashlib.md5(data).hexdigest()}")
    bin_dir = tmp_path / "hg003-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
        "RTG_LOG": str(tmp_path / "rtg.log"),
    }
    pins = (
        f"hg003_fastq_files='{' '.join(pins_fastq)}';"
        f" hg003_truth_files='{' '.join(pins_truth)}';"
        f' record grch38_fasta "{grch38_fasta(home)}"'
    )
    return Calibration(env, pins, serve, fastq, truth)


def _hg003_paths(home: Path) -> dict[str, Path]:
    fastq = home / "calibration" / "hg003" / "fastq"
    return {
        "fastq": fastq,
        "downsampled": fastq.parent / "downsampled",
        "truth": home / "giab" / "HG003-v4.2.1-GRCh38",
    }


def test_calibration_needs_flag(calibration: Calibration) -> None:
    c = calibration
    where = _hg003_paths(c.home)
    without = _run(c.env, f"{c.pins}; step_hg003")
    assert without.returncode == 0, without.stderr
    assert "run with --with-calibration" in without.stderr
    assert c.curl_calls() == []
    assert not where["fastq"].exists()

    done = _run(c.env, f"{c.pins}; with_calibration=1; step_hg003")
    assert done.returncode == 0, done.stderr
    for name in c.fastq:
        assert (where["fastq"] / name).read_bytes() == (c.serve / name).read_bytes()
        assert not (where["fastq"] / f"{name}.part").exists()
    for name in c.truth:
        assert (where["truth"] / name).read_bytes() == (c.serve / name).read_bytes()
    assert not Path(c.env["RTG_LOG"]).exists()  # the SDF is the masked reference's
    calls = c.curl_calls()
    assert sorted(call.split()[-1] for call in calls) == sorted(
        [f"{HG003_FASTQ_BASE}/{n}" for n in c.fastq]
        + [f"{HG003_TRUTH_BASE}/{n}" for n in c.truth]
    )
    for call in calls:  # every transfer resumes a partial file
        assert " -C - " in f" {call} "
    manifest = (c.home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert sorted({row.split("\t")[0] for row in manifest}) == [
        "grch38_fasta",
        "hg003_truth",
    ]

    again = _run(
        c.env, f"{c.pins}; with_calibration=1; step_hg003; check=1; step_hg003"
    )
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert c.curl_calls() == calls


def test_check_lists_calibration_inputs(calibration: Calibration) -> None:
    c = calibration
    where = _hg003_paths(c.home)
    for flags in ("check=1", "check=1; with_calibration=1"):
        listed = _run(c.env, f"{c.pins}; {flags}; step_hg003")
        assert listed.returncode == 0, listed.stderr
        lines = listed.stdout.splitlines()
        assert lines[:3] == [
            f"missing: GIAB HG003 v4.2.1 truth file {n} ({HG003_TRUTH_BASE}/{n})"
            for n in c.truth
        ]
        assert [line.split(" (")[0] for line in lines[3:]] == [
            f"missing: HG003 FASTQ {where['fastq'] / n}" for n in c.fastq
        ]
    assert c.curl_calls() == []
    # Once downsampled, the full FASTQ is not wanted again: it is deleted, never refetched.
    where["downsampled"].mkdir(parents=True)
    for name in c.fastq:
        (where["downsampled"] / name).write_bytes(b"")
    listed = _run(c.env, f"{c.pins}; check=1; step_hg003")
    assert "FASTQ" not in listed.stdout
    assert len(listed.stdout.splitlines()) == 3


@pytest.mark.parametrize("held", ["aligned.json", "align/k/sample.cram"])
def test_hg003_not_wanted_once_recorded(calibration: Calibration, held: str) -> None:
    """HG003 with a record, or with a CRAM under `align/` and no record yet, is re-aligned
    from its CRAM (`hg002 realign hg003`): its FASTQ is neither listed nor fetched, under
    `--with-calibration` too; its truth still is."""
    c = calibration
    where = _hg003_paths(c.home)
    record = where["fastq"].parent / held
    record.parent.mkdir(parents=True)
    record.write_text("{}")
    listed = _run(c.env, f"{c.pins}; check=1; step_hg003")
    assert "FASTQ" not in listed.stdout and len(listed.stdout.splitlines()) == 3
    done = _run(c.env, f"{c.pins}; with_calibration=1; step_hg003")
    assert done.returncode == 0, done.stderr
    assert "FASTQ" not in done.stderr
    assert not where["fastq"].exists()
    assert sorted(call.split()[-1] for call in c.curl_calls()) == sorted(
        f"{HG003_TRUTH_BASE}/{n}" for n in c.truth
    )


def test_hg003_cut_transfer_resumes_and_md5_refuses(calibration: Calibration) -> None:
    c = calibration
    where = _hg003_paths(c.home)
    cut = {**c.env, "CUT_MARK": str(c.serve.parent / "cut"), "CUT_BYTES": "10"}
    first = _run(cut, f"{c.pins}; with_calibration=1; step_hg003")
    assert first.returncode != 0
    part = where["fastq"] / f"{c.fastq[0]}.part"
    assert part.stat().st_size == 10  # kept, to be resumed
    landed = where["fastq"] / c.fastq[0]
    assert not landed.exists()
    resumed = _run(cut, f"{c.pins}; with_calibration=1; step_hg003")
    assert resumed.returncode == 0, resumed.stderr
    assert landed.read_bytes() == (c.serve / c.fastq[0]).read_bytes()
    # A served file whose md5 is not the pinned one is refused and removed.
    landed.unlink()
    (c.serve / c.fastq[0]).write_bytes(b"other bytes")
    refused = _run(c.env, f"{c.pins}; with_calibration=1; step_hg003")
    assert refused.returncode == 1
    assert f"md5 mismatch for {HG003_FASTQ_BASE}/{c.fastq[0]}" in refused.stderr
    assert not landed.exists() and not part.exists()


def test_hg003_pins_match_genokeeper(lab: dict[str, str]) -> None:
    """The truth files the script pins are those `genokeeper.calibrate` reads, and the rtg
    it installs is the one `genokeeper.calibrate.rtg_path` names."""
    shown = _run(
        lab, 'echo "$hg003_truth_files"; echo "$hg003_fastq_files"; echo "$rtg_version"'
    )
    truth, fastq, rtg = shown.stdout.splitlines()
    assert [e.split(":")[0] for e in truth.split()] == [
        TRUTH_VCF,
        f"{TRUTH_VCF}.tbi",
        BENCH_BED,
    ]
    assert all(len(e.split(":")[1]) == 32 for e in truth.split())
    assert [e.split(":")[0] for e in fastq.split()] == [
        "E200029822_L01_UDB-451_1.fq.gz",
        "E200029822_L01_UDB-451_2.fq.gz",
    ]
    home = Path(lab["GENOKEEPER_HOME"])
    assert rtg_path(home) == home / "tools" / f"rtg-tools-{rtg}" / "rtg"


# --- slice 1b: PharmCAT's positions and CPIC's tables --------------------------------------

PHARMCAT_BASE = "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0"
PHARMCAT_FILES = (
    "pharmcat_positions_3.4.0.vcf.bgz",
    "pharmcat_positions_3.4.0.vcf.bgz.csi",
)
CPIC_API = "https://api.cpicpgx.org/v1"
# Each file's query and the rows CPIC's API answers it with.
CPIC_QUERIES = {
    "cpic_pairs.json": (
        "pair_view?select=genesymbol,drugname,cpiclevel",
        [{"genesymbol": "CYP2C19", "drugname": "clopidogrel", "cpiclevel": "A"}],
    ),
    "MT-RNR1_alleles.json": (
        "allele?genesymbol=eq.MT-RNR1&select=name",
        [{"name": "Reference"}, {"name": "m.1555A>G"}],
    ),
}
# FAKE_SERVING_CURL that logs every argument, not only the URL.
FAKE_LOGGING_CURL = """#!/usr/bin/env bash
out=
for a; do [ "${prev:-}" = -o ] && out=$a; prev=$a; done
echo "$*" >> "$CURL_LOG"
src=$SERVE/${prev##*/}
[ -f "$src" ] || exit 22
cp "$src" "$out"
"""


def _serving(lab: dict[str, str], tmp_path: Path) -> tuple[dict[str, str], Path]:
    serve = tmp_path / "serve"
    serve.mkdir(exist_ok=True)
    bin_dir = tmp_path / "serving-bin"
    _executable(bin_dir / "curl", FAKE_LOGGING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    return env, serve


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_pharmcat_positions_pinned_by_sha256(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The release assets are kept only when their sha256 is the one GitHub publishes,
    pinned in the script: other bytes are refused and nothing is recorded."""
    env, serve = _serving(lab, tmp_path)
    for name in PHARMCAT_FILES:
        (serve / name).write_bytes(f"served {name}".encode())
    home = Path(lab["GENOKEEPER_HOME"])
    dest = home / "pharmcat" / "3.4.0"
    refused = _run(env, "step_pharmcat")
    assert refused.returncode == 1
    assert f"sha256 mismatch for {PHARMCAT_BASE}/{PHARMCAT_FILES[0]}" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not list(dest.glob("**/*.bgz*"))
    listed = _run(env, "check=1; step_pharmcat")
    assert listed.stdout.splitlines() == [
        f"missing: PharmCAT 3.4.0 positions file {name} ({PHARMCAT_BASE}/{name})"
        for name in PHARMCAT_FILES
    ]
    # Pins that match the served bytes (a test cannot serve the real ones).
    pins = " ".join(f"{name}:{_sha256(serve / name)}" for name in PHARMCAT_FILES)
    done = _run(env, f'pharmcat_files="{pins}"; step_pharmcat')
    assert done.returncode == 0, done.stderr
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [
        f"pharmcat\t{dest / name}\t{_sha256(serve / name)}" for name in PHARMCAT_FILES
    ]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls[-2:]] == [
        f"{PHARMCAT_BASE}/{name}" for name in PHARMCAT_FILES
    ]
    again = _run(env, f'pharmcat_files="{pins}"; step_pharmcat; check=1; step_pharmcat')
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == len(calls)


def test_pharmcat_pins_match_genokeeper(lab: dict[str, str]) -> None:
    """The positions file the script pins is the one `genokeeper.pgx.pharmcat` reads, and
    each pin is a sha256."""
    shown = _run(lab, 'echo "$pharmcat_version"; echo "$pharmcat_files"')
    version, files = shown.stdout.splitlines()
    assert version == PHARMCAT_VERSION
    assert [e.split(":")[0] for e in files.split()] == list(PHARMCAT_FILES)
    assert all(len(e.split(":")[1]) == 64 for e in files.split())
    home = Path(lab["GENOKEEPER_HOME"])
    assert positions_path(home) == home / "pharmcat" / version / PHARMCAT_FILES[0]


IMGT_HLA_URL = (
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/v3.66.0-alpha/hla.dat.zip"
)
# Stands in for `python -m genokeeper.pgx.hla build-index DAT OUT_DIR`: writes the index
# files T1K's build writes, the first naming what it was built from.
T1K_BUILD_STUB = (
    't1k_build() { mkdir -p "$2"; cat "$1" > "$2/hla_dna_seq.fa";'
    ' printf ">rna\\n" > "$2/hla_rna_seq.fa"; }'
)


def _zip(members: dict[str, str]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        for name, text in members.items():
            z.writestr(name, text)
    return buffer.getvalue()


def test_t1k_index_built_from_the_pinned_release(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """IPD-IMGT/HLA's hla.dat.zip is fetched whole and kept only with its pinned sha256;
    T1K's index is built from its hla.dat and recorded with it."""
    env, serve = _serving(lab, tmp_path)
    (serve / "hla.dat.zip").write_bytes(_zip({"hla.dat": "ID   HLA00001\n"}))
    home = Path(lab["GENOKEEPER_HOME"])
    index = t1k_index(home)
    refused = _run(env, f"{T1K_BUILD_STUB}; step_t1k")
    assert refused.returncode == 1
    assert f"sha256 mismatch for {IMGT_HLA_URL}" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not list(index.parent.glob("*"))
    listed = _run(env, "check=1; step_t1k")
    assert listed.stdout.splitlines() == [
        f"missing: T1K HLA index from IPD-IMGT/HLA 3.66.0 ({IMGT_HLA_URL})"
    ]
    pin = _sha256(serve / "hla.dat.zip")
    done = _run(env, f"{T1K_BUILD_STUB}; imgt_hla_sha256={pin}; step_t1k")
    assert done.returncode == 0, done.stderr
    assert index.read_text() == "ID   HLA00001\n"
    zipped = index.parent / "hla.dat.zip"
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [
        f"t1k\t{zipped}\t{pin}",
        f"t1k\t{index}\t{_sha256(index)}",
    ]
    assert sorted(p.name for p in index.parent.iterdir()) == [
        "hla.dat.zip",
        "hla_dna_seq.fa",
        "hla_rna_seq.fa",
    ]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert calls[-1].split()[-1] == IMGT_HLA_URL and " -C " not in calls[-1]
    again = _run(env, f"imgt_hla_sha256={pin}; step_t1k; check=1; step_t1k")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == len(calls)
    # A zip without hla.dat is refused before any index is built.
    (serve / "hla.dat.zip").write_bytes(_zip({"other.txt": "x"}))
    pin = _sha256(serve / "hla.dat.zip")
    shutil.rmtree(index.parent)
    empty = _run(env, f"{T1K_BUILD_STUB}; imgt_hla_sha256={pin}; step_t1k")
    assert empty.returncode == 1 and "holds no hla.dat" in empty.stderr


def test_imgt_hla_pin_matches_genokeeper(lab: dict[str, str]) -> None:
    shown = _run(lab, 'echo "$imgt_hla_version"; echo "$imgt_hla_sha256"')
    version, pin = shown.stdout.splitlines()
    assert version == IMGT_HLA_VERSION and len(pin) == 64
    home = Path(lab["GENOKEEPER_HOME"])
    assert t1k_index(home).parent.name == f"IPD-IMGT-HLA-{version}"


def test_cpic_data_fetched_whole_under_its_fetch_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    env, serve = _serving(lab, tmp_path)
    for query, rows in CPIC_QUERIES.values():
        (serve / query).write_text(json.dumps(rows))
    home = Path(lab["GENOKEEPER_HOME"])
    # A folder of CPIC's old workbooks, recorded, is not the data the run reads.
    old = home / "cpic" / "20261019" / "cpic_gene-drug_pairs.xlsx"
    old.parent.mkdir(parents=True)
    old.write_bytes(b"PK")
    _run(env, f"record cpic {old}")
    done = _run(env, "step_cpic")
    assert done.returncode == 0, done.stderr
    dest = home / "cpic" / datetime.now(UTC).strftime("%Y%m%d")
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [
        f"cpic\t{dest / n}\t{_sha256(serve / q)}" for n, (q, _) in CPIC_QUERIES.items()
    ]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls] == [
        f"{CPIC_API}/{q}" for q, _ in CPIC_QUERIES.values()
    ]
    assert all(" -C " not in f" {c} " for c in calls)  # whole, never resumed
    assert all(" --speed-limit " in f" {c} " for c in calls)  # a stall fails
    again = _run(env, "step_cpic; check=1; step_cpic")
    assert (again.returncode, again.stdout) == (0, "")
    # An answer the reader does not take is refused, and nothing is recorded.
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "cpic")
    (serve / CPIC_QUERIES["MT-RNR1_alleles.json"][0]).write_text("[]")
    refused = _run(env, "step_cpic")
    assert refused.returncode == 1
    assert f"{CPIC_API}/ answered with no CPIC rows genokeeper reads" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not [p for p in (home / "cpic").iterdir() if not p.name.startswith(".")]
    listed = _run(env, "check=1; step_cpic")
    assert listed.stdout == (
        f"missing: CPIC gene-drug pairs and MT-RNR1 alleles ({CPIC_API}/)\n"
    )


def test_every_download_is_bounded() -> None:
    """One curl command line, with a connect timeout and a stall bound: a fetch that stops
    moving fails (and is retried) instead of blocking the bootstrap or a refresh."""
    calls = [
        line.split()
        for line in SCRIPT.read_text().splitlines()
        if "curl " in line and not line.lstrip().startswith("#")
    ]
    assert len(calls) == 1
    assert {"--connect-timeout", "--speed-limit", "--speed-time"} <= set(calls[0])


def test_cpic_rows_unchanged_since_an_earlier_day_are_not_filed_again(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """CPIC dates nothing: a fetch on a new day of the same rows is no new release."""
    env, serve = _serving(lab, tmp_path)
    for query, rows in CPIC_QUERIES.values():
        (serve / query).write_text(json.dumps(rows))
    home = Path(lab["GENOKEEPER_HOME"])
    old = home / "cpic" / "20200101"
    old.mkdir(parents=True)
    for name, (query, _) in CPIC_QUERIES.items():
        shutil.copy(serve / query, old / name)
    done = _refresh(env, "cpic")
    assert (done.returncode, done.stdout) == (0, f"refreshed: cpic {old}\n")
    assert _folders(home / "cpic") == ["20200101"]
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"cpic\t{old / n}\t{_sha256(old / n)}" for n in CPIC_QUERIES]
    # One row changed: filed under the fetch date.
    pairs = [{"genesymbol": "CYP2D6", "drugname": "codeine", "cpiclevel": "A"}]
    (serve / CPIC_QUERIES["cpic_pairs.json"][0]).write_text(json.dumps(pairs))
    new = home / "cpic" / datetime.now(UTC).strftime("%Y%m%d")
    assert _refresh(env, "cpic").stdout == f"refreshed: cpic {new}\n"
    assert _folders(home / "cpic") == sorted(["20200101", new.name])


def test_image_without_upstream_built_not_pulled(stage2: dict[str, str]) -> None:
    """Aldy and BCyrius have no BioContainers image: one whose input-named tag is absent
    is built from its Dockerfile, never pulled; one present is left alone."""
    absent = "aldy"
    others = " ".join(ref for name, ref in IMAGES.items() if name != absent)
    done = _run({**stage2, "DOCKER_DIGESTS": others}, "step_images")
    assert done.returncode == 0, done.stderr
    calls = Path(stage2["DOCKER_LOG"]).read_text().splitlines()
    assert [line for line in calls if line.startswith(("pull", "build"))] == [
        " ".join(build_command(absent)[1:])
    ]
    assert {"aldy", "bcyrius"} <= set(BUILT)
    # --check lists it as it lists an image to pull.
    listed = _run({**stage2, "DOCKER_DIGESTS": others}, "check=1; step_images")
    assert listed.stdout == f"missing: container image {absent} ({IMAGES[absent]})\n"


def test_multi_arch_image_pulled_for_the_host(stage2: dict[str, str]) -> None:
    """An image pinned by its multi-platform index digest is pulled without a platform:
    Docker takes the host's manifest from the pinned index."""
    others = [ref for name, ref in IMAGES.items() if name != "pharmcat"]
    present = {**stage2, "DOCKER_DIGESTS": " ".join(others)}
    done = _run(present, "step_images")
    assert done.returncode == 0, done.stderr
    calls = Path(stage2["DOCKER_LOG"]).read_text().splitlines()
    assert [line for line in calls if line.startswith("pull")] == [
        f"pull {IMAGES['pharmcat']}"
    ]
    assert "pharmcat" not in PLATFORMS


PHARMGKB_URL = "https://api.clinpgx.org/v1/download/file/data/summaryAnnotations.zip"


def test_pharmgkb_annotations_fetched_whole_under_their_created_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """ClinPGx's live summary annotations download is fetched whole (never resumed: the next
    request may serve a newer file), checked by genokeeper's reader, and filed under the
    date its own CREATED member names, with its sha256 in the manifest."""
    env, serve = _serving(lab, tmp_path)
    served = serve / "summaryAnnotations.zip"
    served.write_bytes(zip_bytes([Ann("1", "rs1", "TNF", "2B")], created="2026-09-30"))
    home = Path(lab["GENOKEEPER_HOME"])
    listed = _run(env, "check=1; step_pharmgkb")
    assert listed.stdout == f"missing: ClinPGx clinical annotations ({PHARMGKB_URL})\n"
    done = _run(env, "step_pharmgkb")
    assert done.returncode == 0, done.stderr
    dest = home / "pharmgkb" / "20260930" / "summaryAnnotations.zip"
    assert dest.read_bytes() == served.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"pharmgkb\t{dest}\t{_sha256(served)}"]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls] == [PHARMGKB_URL]
    assert " -C " not in f" {calls[0]} "
    again = _run(env, "step_pharmgkb; check=1; step_pharmgkb")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 1
    # A page that is not the download is refused, and nothing is recorded.
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "pharmgkb")
    served.write_bytes(b"<html>maintenance</html>")
    refused = _run(env, "step_pharmgkb")
    assert refused.returncode == 1
    assert f"{PHARMGKB_URL} is not the ClinPGx clinical annotations" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not [p for p in (home / "pharmgkb").iterdir() if not p.name.startswith(".")]


GWAS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip"
)


def test_gwas_associations_fetched_whole_under_their_build_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The GWAS Catalog's associations are fetched whole (never resumed), checked by
    genokeeper's reader, and filed under the build date their member carries, with their
    sha256 in the manifest; a file that is not the download is refused."""
    env, serve = _serving(lab, tmp_path)
    served = write_gwas(serve / GWAS_DOWNLOAD, [Assoc(102, "A")])
    home = Path(lab["GENOKEEPER_HOME"])
    listed = _run(env, "check=1; step_gwas")
    assert listed.stdout == f"missing: GWAS Catalog associations ({GWAS_URL})\n"
    done = _run(env, "step_gwas")
    assert done.returncode == 0, done.stderr
    dest = home / "gwas" / "20261003" / GWAS_DOWNLOAD
    assert dest.read_bytes() == served.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"gwas\t{dest}\t{_sha256(served)}"]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls] == [GWAS_URL]
    assert " -C " not in f" {calls[0]} "
    again = _run(env, "step_gwas; check=1; step_gwas")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 1
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "gwas")
    served.write_bytes(b"<html>maintenance</html>")
    refused = _run(env, "step_gwas")
    assert refused.returncode == 1
    assert f"{GWAS_URL} is not the GWAS Catalog associations download" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not [p for p in (home / "gwas").iterdir() if not p.name.startswith(".")]


def test_gwas_file_matches_genokeeper(lab: dict[str, str]) -> None:
    """The file the script files is the one the run looks for."""
    shown = _run(lab, 'echo "$gwas_file"; echo "$gwas_url"')
    name, url = shown.stdout.splitlines()
    assert name == GWAS_DOWNLOAD and url.endswith(f"/{GWAS_DOWNLOAD}")


def test_cpic_files_match_genokeeper(lab: dict[str, str]) -> None:
    """The CPIC data the script fetches is what the run reads."""
    shown = _run(lab, 'echo "$cpic_api"; echo "$cpic_files"')
    api, files = shown.stdout.splitlines()
    assert api == CPIC_API
    assert files.split() == [f"{n}:{q}" for n, (q, _) in CPIC_QUERIES.items()]
    assert list(CPIC_QUERIES) == list(CPIC_FILE_NAMES.values())


# --- slice 1d: polygenic scores ------------------------------------------------------------

PGS_SCORES_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"
PGS_PANEL_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst"
)
# Stands in for `python -m genokeeper.prs.pgsc ancestry-positions PANEL OUT_DIR`: writes the
# ancestry positions and the panel's variants as the name of the panel they came from, and
# logs it.
PGS_ANCESTRY_STUB = (
    'pgs_ancestry_build() { mkdir -p "$2"; echo "$1" >> "$BUILD_LOG";'
    ' basename "$1" > "$2/GRCh38_HGDP+1kGP_ancestry.pvar.zst";'
    ' basename "$1" > "$2/GRCh38_HGDP+1kGP_panel.pvar.zst"; }'
)
# `brew --prefix FORMULA` names a folder whose bin/java exists once `brew install` ran.
FAKE_JAVA_BREW = """#!/usr/bin/env bash
echo "$*" >> "$BREW_LOG"
case "$1" in
  --prefix) echo "$BREW_ROOT/opt/$2" ;;
  install) mkdir -p "$BREW_ROOT/opt/$2/bin"; touch "$BREW_ROOT/opt/$2/bin/java";
           chmod +x "$BREW_ROOT/opt/$2/bin/java" ;;
esac
"""


def test_pgs_scoring_files_checked_against_their_published_md5(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """Each harmonised GRCh38 scoring file is fetched whole with the md5 the PGS Catalog
    publishes beside it, and kept only when they agree, under the path the run reads; a
    re-harmonised file with its new md5 is kept with no edit to the script."""
    env, serve = _serving(lab, tmp_path)
    ids = [s.pgs_id for s in load_scores()]
    names = [f"{i}_hmPOS_GRCh38.txt.gz" for i in ids]
    for name in names:
        (serve / name).write_bytes(gzip.compress(f"#pgs_id={name[:9]}\n".encode()))
        (serve / f"{name}.md5").write_text(f"{'0' * 32}  {name}\n")
    home = Path(lab["GENOKEEPER_HOME"])
    refused = _run(env, "step_pgs_scores")
    assert refused.returncode == 1
    assert f"md5 mismatch for {PGS_SCORES_BASE}/{ids[0]}/" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    listed = _run(env, "check=1; step_pgs_scores")
    assert listed.stdout == (
        f"missing: {len(ids)} PGS Catalog scoring files ({PGS_SCORES_BASE}/)\n"
    )
    for name in names:
        md5 = hashlib.md5((serve / name).read_bytes()).hexdigest()
        (serve / f"{name}.md5").write_text(f"{md5}  {name}\n")
    done = _run(env, "step_pgs_scores")
    assert done.returncode == 0, done.stderr
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert [row.split("\t")[1] for row in manifest] == [
        str(scoring_file(home, i)) for i in ids
    ]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls[-2 * len(ids) :]] == [
        f"{PGS_SCORES_BASE}/{i}/ScoringFiles/Harmonized/{n}{suffix}"
        for i, n in zip(ids, names, strict=True)
        for suffix in ("", ".md5")
    ]
    assert all(" -C " not in f" {c} " for c in calls)
    again = _run(env, "step_pgs_scores; check=1; step_pgs_scores")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == len(calls)


def test_research_scoring_files_are_fetched_without_touching_the_chosen(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """A lab that holds the chosen scoring files gets only the research ones, each md5
    checked: a chosen file fetched again could come back re-harmonised and move the
    chosen run's cache key."""
    env, serve = _serving(lab, tmp_path)
    ids = [s.pgs_id for s in load_research_scores()]
    chosen = [s.pgs_id for s in load_scores()]
    for i in [*chosen, *ids]:
        name = f"{i}_hmPOS_GRCh38.txt.gz"
        (serve / name).write_bytes(gzip.compress(f"#pgs_id={i}\n".encode()))
        md5 = hashlib.md5((serve / name).read_bytes()).hexdigest()
        (serve / f"{name}.md5").write_text(f"{md5}  {name}\n")
    home = Path(lab["GENOKEEPER_HOME"])
    assert _run(env, "step_pgs_scores").returncode == 0
    before = {i: scoring_file(home, i).read_bytes() for i in chosen}
    calls = len(Path(env["CURL_LOG"]).read_text().splitlines())
    name = f"{ids[0]}_hmPOS_GRCh38.txt.gz"
    (serve / f"{name}.md5").write_text(f"{'0' * 32}  {name}\n")
    refused = _run(env, "step_pgs_research_scores")
    assert f"md5 mismatch for {PGS_SCORES_BASE}/{ids[0]}/" in refused.stderr
    assert not scoring_file(home, ids[0]).exists()
    assert (
        _run(env, "check=1; step_pgs_research_scores").stdout
        == f"missing: {len(ids)} PGS Catalog research scoring files"
        f" ({PGS_SCORES_BASE}/)\n"
    )
    md5 = hashlib.md5((serve / name).read_bytes()).hexdigest()
    (serve / f"{name}.md5").write_text(f"{md5}  {name}\n")
    done = _run(env, "step_pgs_research_scores")
    assert done.returncode == 0, done.stderr
    assert all(scoring_file(home, i).is_file() for i in ids)
    assert {i: scoring_file(home, i).read_bytes() for i in chosen} == before
    fetched = Path(env["CURL_LOG"]).read_text().splitlines()[calls:]
    assert not [c for c in fetched if any(f"/{i}/" in c for i in chosen)]
    again = _run(env, "step_pgs_scores; check=1; step_pgs_scores")
    assert (again.returncode, again.stdout) == (0, "")


def test_pgs_panel_moved_from_prefetch_when_its_md5_matches(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The HGDP+1kGP panel downloaded before under prefetch/ is moved into place when its
    md5 is the published one, and the ancestry positions are built from it; a copy with
    another md5 is left where it is and the panel is downloaded, md5-checked, instead."""
    serve = tmp_path / "serve"
    serve.mkdir()
    bin_dir = tmp_path / "resuming-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
        "BUILD_LOG": str(tmp_path / "build.log"),
    }
    home = Path(lab["GENOKEEPER_HOME"])
    panel = panel_path(home)
    prefetch = home / "prefetch" / "1d" / "pgsc_calc" / PANEL_NAME
    prefetch.parent.mkdir(parents=True)
    prefetch.write_bytes(b"prefetched panel")
    listed = _run(env, "check=1; step_pgs_panel")
    assert listed.stdout.splitlines() == [
        f"missing: HGDP+1kGP reference panel {panel} ({PGS_PANEL_URL}, 16 GB)",
        (
            f"missing: ancestry positions {ancestry_pvar(home)} and the panel's"
            f" variants {panel_pvar(home)} (built from the panel)"
        ),
    ]
    pin = hashlib.md5(b"prefetched panel").hexdigest()
    moved = _run(env, f"{PGS_ANCESTRY_STUB}; pgs_panel_md5={pin}; step_pgs_panel")
    assert moved.returncode == 0, moved.stderr
    assert panel.read_bytes() == b"prefetched panel" and not prefetch.exists()
    assert not Path(env["CURL_LOG"]).exists()
    assert Path(env["BUILD_LOG"]).read_text().splitlines() == [str(panel)]
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    pvar, full = ancestry_pvar(home), panel_pvar(home)
    assert manifest == [
        f"pgs_ancestry\t{pvar}\t{_sha256(pvar)}",
        f"pgs_ancestry\t{full}\t{_sha256(full)}",
    ]
    again = _run(env, f"pgs_panel_md5={pin}; step_pgs_panel; check=1; step_pgs_panel")
    assert (again.returncode, again.stdout) == (0, "")
    # Another md5 under prefetch/: left there, and the published panel downloaded.
    panel.unlink()
    (home / "bootstrap-manifest.tsv").unlink()
    prefetch.write_bytes(b"a truncated panel")
    (serve / PANEL_NAME).write_bytes(b"the published panel")
    pin = hashlib.md5(b"the published panel").hexdigest()
    fetched = _run(env, f"{PGS_ANCESTRY_STUB}; pgs_panel_md5={pin}; step_pgs_panel")
    assert fetched.returncode == 0, fetched.stderr
    assert panel.read_bytes() == b"the published panel"
    assert prefetch.read_bytes() == b"a truncated panel"
    (call,) = Path(env["CURL_LOG"]).read_text().splitlines()
    assert call.split()[-1] == PGS_PANEL_URL and " -C - " in f" {call} "
    # A download with another md5 is removed and refused; nothing is built from it.
    panel.unlink()
    refused = _run(
        env, f"{PGS_ANCESTRY_STUB}; pgs_panel_md5={'0' * 32}; step_pgs_panel"
    )
    assert refused.returncode == 1 and "md5 mismatch" in refused.stderr
    assert not panel.exists() and not panel.with_name(f"{PANEL_NAME}.part").exists()
    assert len(Path(env["BUILD_LOG"]).read_text().splitlines()) == 2


def _tar_gz(top: str, members: dict[str, str]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(f"{top}/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def _duplicated_zip(name: str) -> bytes:
    """A zip listing `name` twice, as nf-schema 2.0.0's release zip lists its
    classes/META-INF/ entries."""
    buffer = io.BytesIO()
    with warnings.catch_warnings(), zipfile.ZipFile(buffer, "w") as z:
        warnings.simplefilter("ignore")  # zipfile warns of the duplicate name
        z.writestr(name, "m\n")
        z.writestr(name, "m\n")
    return buffer.getvalue()


def test_pgs_tools_pinned_by_sha256(lab: dict[str, str], tmp_path: Path) -> None:
    """Nextflow, its nf-schema plugin and pgsc_calc at the commit tag v2.3.0 points to are
    each kept only with the sha256 pinned in the script, unpacked where the run reads them
    (the plugin under NXF_HOME, so Nextflow runs offline), with Java 21 from Homebrew."""
    env, serve = _serving(lab, tmp_path)
    env |= {"BREW_LOG": str(tmp_path / "brew.log"), "BREW_ROOT": str(tmp_path / "brew")}
    _executable(Path(env["PATH"].split(":")[0]) / "brew", FAKE_JAVA_BREW)
    served = {
        f"nextflow-{NEXTFLOW_VERSION}-dist": b"#!/bin/sh\necho nextflow\n",
        f"nf-schema-{NF_SCHEMA_VERSION}.zip": _duplicated_zip("META-INF/MANIFEST.MF"),
        f"{PGSC_CALC_COMMIT}.tar.gz": _tar_gz(
            f"pgsc_calc-{PGSC_CALC_COMMIT}", {"main.nf": "workflow {}\n"}
        ),
    }
    for name, data in served.items():
        (serve / name).write_bytes(data)
    home = Path(lab["GENOKEEPER_HOME"])
    refused = _run(env, "step_pgs_tools")
    assert refused.returncode == 1
    assert "sha256 mismatch for https://github.com/nextflow-io/" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()
    listed = _run(env, "check=1; step_pgs_tools")
    assert listed.stdout == (
        f"missing: Nextflow {NEXTFLOW_VERSION}, nf-schema {NF_SCHEMA_VERSION},"
        f" pgsc_calc {PGSC_CALC_VERSION} ({PGSC_CALC_COMMIT}) and openjdk@21\n"
    )
    pins = "; ".join(
        f"{var}={_sha256(serve / name)}"
        for var, name in zip(
            ("nextflow_sha256", "nf_schema_sha256", "pgsc_calc_sha256"),
            served,
            strict=True,
        )
    )
    done = _run(env, f"{pins}; step_pgs_tools")
    assert done.returncode == 0, done.stderr
    assert os.access(nextflow_path(home), os.X_OK)
    assert (pgsc_calc_dir(home) / "main.nf").read_text() == "workflow {}\n"
    plugin = nxf_home(home) / "plugins" / f"nf-schema-{NF_SCHEMA_VERSION}"
    assert (plugin / "META-INF" / "MANIFEST.MF").is_file()
    assert os.access(java_home(home) / "bin" / "java", os.X_OK)
    assert "install openjdk@21" in Path(env["BREW_LOG"]).read_text().splitlines()
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    again = _run(env, f"{pins}; step_pgs_tools; check=1; step_pgs_tools")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == len(calls)


def test_stellarpgx_is_unpacked_only_with_its_pin(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """StellarPGx's archive is kept only with the sha256 pinned in the script, and the commit
    pinned there is the one `genokeeper.pgx.stellarpgx` reads."""
    env, serve = _serving(lab, tmp_path)
    (serve / f"{STELLARPGX_COMMIT}.tar.gz").write_bytes(
        _tar_gz(f"StellarPGx-{STELLARPGX_COMMIT}", {"main.nf": "workflow {}\n"})
    )
    home = Path(lab["GENOKEEPER_HOME"])
    refused = _run(env, "step_stellarpgx")
    assert refused.returncode == 1
    assert "sha256 mismatch for https://github.com/SBIMB/" in refused.stderr
    assert _run(env, "check=1; step_stellarpgx").stdout == (
        f"missing: StellarPGx ({STELLARPGX_COMMIT})\n"
    )
    pin = f"stellarpgx_sha256={_sha256(serve / f'{STELLARPGX_COMMIT}.tar.gz')}"
    assert _run(env, f"{pin}; step_stellarpgx").returncode == 0
    assert (stellarpgx_dir(home) / "main.nf").read_text() == "workflow {}\n"
    assert _run(env, f"{pin}; check=1; step_stellarpgx").stdout == ""


def test_pgs_pins_match_genokeeper(lab: dict[str, str]) -> None:
    """The versions, commit and files the script pins are the ones `genokeeper.prs` reads:
    every chosen score, then every research score, in order."""
    shown = _run(
        lab,
        'echo "$nextflow_version"; echo "$nf_schema_version"; echo "$pgsc_calc_commit";'
        ' echo "$pgs_panel_name"; echo "$pgs_panel_pvar"; echo $pgs_ids;'
        " echo $pgs_research_ids",
    )
    nextflow, nf_schema, commit, panel, full, ids, research = shown.stdout.splitlines()
    assert (nextflow, nf_schema, commit, panel, full) == (
        NEXTFLOW_VERSION,
        NF_SCHEMA_VERSION,
        PGSC_CALC_COMMIT,
        PANEL_NAME,
        PANEL_PVAR_NAME,
    )
    assert ids.split() == [s.pgs_id for s in load_scores()]
    assert research.split() == [s.pgs_id for s in load_research_scores()]


# --- slice 1c: MedlinePlus Genetics summaries -----------------------------------------------


def test_medlineplus_fetched_whole_under_its_fetch_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """MedlinePlus rebuilds its summaries daily under one URL: the file is fetched whole
    (a resumed part could splice two builds), checked by genokeeper's loader, and filed
    under the fetch date with its sha256 in the manifest."""
    env, serve = _serving(lab, tmp_path)
    served = write_medlineplus(serve / MEDLINEPLUS_FILE, medlineplus_gene("HFE", "x"))
    home = Path(lab["GENOKEEPER_HOME"])
    listed = _run(env, "check=1; step_medlineplus")
    assert listed.stdout == (
        f"missing: MedlinePlus Genetics summaries ({MEDLINEPLUS_URL})\n"
    )
    done = _run(env, "step_medlineplus")
    assert done.returncode == 0, done.stderr
    dest = (
        home / "medlineplus" / datetime.now(UTC).strftime("%Y%m%d") / MEDLINEPLUS_FILE
    )
    assert dest.read_bytes() == served.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"medlineplus\t{dest}\t{_sha256(served)}"]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls] == [MEDLINEPLUS_URL]
    assert " -C " not in f" {calls[0]} "
    again = _run(env, "step_medlineplus; check=1; step_medlineplus")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 1
    # A page that is not the summaries file is refused, and nothing is recorded.
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "medlineplus")
    served.write_text("<html>maintenance</html>")
    refused = _run(env, "step_medlineplus")
    assert refused.returncode == 1
    assert (
        f"{MEDLINEPLUS_URL} is not the MedlinePlus Genetics summaries" in refused.stderr
    )
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not [
        p for p in (home / "medlineplus").iterdir() if not p.name.startswith(".")
    ]


# --- slice 1e: ClinGen's dosage list, VEP's cache and score files ---------------------------


DOSAGE_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"


def test_dosage_list_fetched_whole_under_its_own_date(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """ClinGen's dosage sensitivity list is fetched whole, checked by genokeeper's reader and
    filed under the date it names; a page the reader refuses is not recorded."""
    env, serve = _serving(lab, tmp_path)
    served = write_dosage(serve / DOSAGE_FILE, [("BRCA1", "3")])
    home = Path(lab["GENOKEEPER_HOME"])
    listed = _run(env, "check=1; step_clingen_dosage")
    assert listed.stdout == (
        f"missing: ClinGen dosage sensitivity curations ({DOSAGE_URL})\n"
    )
    done = _run(env, "step_clingen_dosage")
    assert done.returncode == 0, done.stderr
    dest = home / "clingen-dosage" / "20261023" / DOSAGE_FILE
    assert dest.read_bytes() == served.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"clingen_dosage\t{dest}\t{_sha256(served)}"]
    (call,) = Path(env["CURL_LOG"]).read_text().splitlines()
    assert call.split()[-1] == DOSAGE_URL and " -C " not in f" {call} "
    again = _run(env, "step_clingen_dosage; check=1; step_clingen_dosage")
    assert (again.returncode, again.stdout) == (0, "")
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "clingen-dosage")
    served.write_text("<html>maintenance</html>\n")
    refused = _run(env, "step_clingen_dosage")
    assert refused.returncode == 1
    assert "is not ClinGen's gene curation list" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()


REGIONS_URL = f"https://ftp.clinicalgenome.org/{REGIONS_FILE}"
RECURRENT_URL = f"https://ftp.clinicalgenome.org/{RECURRENT_FILE}"


def test_clingen_regions_step(lab: dict[str, str], tmp_path: Path) -> None:
    """ClinGen's region curations and recurrent CNV regions (slice 3d, A5) are fetched whole,
    checked by genokeeper's readers and filed together under the region list's date; a rerun
    fetches nothing; unchanged files are not filed again; a page either reader refuses is
    not recorded."""
    env, serve = _serving(lab, tmp_path)
    regions = write_regions(
        serve / REGIONS_FILE, [("22q11.2 recurrent (DGS)", "chr22:18924719-21111383", "3", "0")]
    )  # fmt: skip
    recurrent = write_recurrent(
        serve / RECURRENT_FILE, [("chr22", 18924718, 21111383, "22q11.2", "255,111,0")]
    )
    home = Path(lab["GENOKEEPER_HOME"])
    listed = _run(env, "check=1; step_clingen_regions")
    assert listed.stdout == (
        "missing: ClinGen region curations and recurrent CNV regions"
        f" ({REGIONS_URL}, {RECURRENT_URL})\n"
    )
    done = _run(env, "step_clingen_regions")
    assert done.returncode == 0, done.stderr
    folder = home / "clingen-regions" / "20260428"
    assert (folder / REGIONS_FILE).read_bytes() == regions.read_bytes()
    assert (folder / RECURRENT_FILE).read_bytes() == recurrent.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [
        f"clingen_regions\t{folder / REGIONS_FILE}\t{_sha256(regions)}",
        f"clingen_regions\t{folder / RECURRENT_FILE}\t{_sha256(recurrent)}",
    ]
    calls = Path(env["CURL_LOG"]).read_text().splitlines()
    assert [c.split()[-1] for c in calls] == [REGIONS_URL, RECURRENT_URL]
    assert not any(" -C " in f" {c} " for c in calls)
    again = _run(env, "step_clingen_regions; check=1; step_clingen_regions")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 2
    # Refreshed with the same curations under a new build date: nothing filed again.
    write_regions(
        serve / REGIONS_FILE,
        [("22q11.2 recurrent (DGS)", "chr22:18924719-21111383", "3", "0")],
        date="29 Apr,2026",
    )
    same = _run(env, "refresh=clingen_regions; step_clingen_regions")
    assert same.returncode == 0, same.stderr
    assert same.stdout == f"refreshed: clingen_regions {folder}\n"
    assert sorted(p.name for p in (home / "clingen-regions").iterdir()) == ["20260428"]
    # A page either reader refuses is not recorded.
    (home / "bootstrap-manifest.tsv").unlink()
    shutil.rmtree(home / "clingen-regions")
    recurrent.write_text("<html>maintenance</html>\n")
    refused = _run(env, "step_clingen_regions")
    assert refused.returncode == 1
    assert (
        "is not ClinGen's region curation list and recurrent CNV BED" in refused.stderr
    )
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert not [
        p for p in (home / "clingen-regions").iterdir() if not p.name.startswith(".")
    ]


GENCC_URL = "https://thegencc.org/download/action/submissions-export-tsv?format=new"
GENCC_HEAD = REPO / "tests" / "fixtures" / "gencc_head.tsv"


def _gencc_row(gene: str, title: str, run_date: str) -> str:
    """One submission row of the new format, the columns `load_gencc` reads filled."""
    header = GENCC_HEAD.read_text().split("\n", 1)[0].split("\t")
    values = {
        "gene_symbol": gene,
        "disease_curie": "MONDO:0000001",
        "classification_title": title,
        "submitted_run_date": run_date,
    }
    return "\t".join(values.get(column, "") for column in header) + "\n"


PANELAPP_NAME = "?format=json"  # the URL's last segment, the fake curl's file


def test_panelapp_panel_pinned_by_sha256(lab: dict[str, str], tmp_path: Path) -> None:
    """One file, kept only when its sha256 is the pin, taken with its `evidence` lists
    sorted (BR-041); another version of the panel, or a failed fetch, is not filed and
    leaves the bootstrap running, and a done step fetches nothing again."""
    env, serve = _serving(lab, tmp_path)
    served = serve / PANELAPP_NAME
    served.write_text('{"version": "other"}')
    home = Path(lab["GENOKEEPER_HOME"])
    dest = panelapp_path(home)
    listed = _run(env, "check=1; step_panelapp")
    assert listed.stdout == (
        "missing: PanelApp's Intellectual disability panel"
        f" {PANEL_VERSION} ({PANELAPP_URL})\n"
    )
    refused = _run(env, "step_panelapp")
    assert refused.returncode == 0, refused.stderr
    assert f"sha256 mismatch for {PANELAPP_URL}" in refused.stderr
    assert not dest.exists()
    assert not (home / "bootstrap-manifest.tsv").exists()
    assert (home / "panelapp" / "285" / PANEL_VERSION).is_dir()
    assert not [p.name for p in dest.parent.iterdir()]
    down = _run(env, "panelapp_url=https://no.where/x; step_panelapp")
    assert down.returncode == 0, down.stderr
    assert "could not fetch https://no.where/x" in down.stderr
    served.write_text('{"version": "x", "genes": [{"evidence": ["b", "a"]}]}')
    sorted_panel = '{"genes": [{"evidence": ["a", "b"]}], "version": "x"}'
    digest = hashlib.sha256(sorted_panel.encode()).hexdigest()
    pin = f"panelapp_sha256={digest}"
    done = _run(env, f"{pin}; step_panelapp")
    assert done.returncode == 0, done.stderr
    assert dest.read_text() == sorted_panel
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"panelapp\t{dest}\t{digest}"]
    calls = len(Path(env["CURL_LOG"]).read_text().splitlines())
    again = _run(env, f"{pin}; step_panelapp; check=1; step_panelapp")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == calls


def test_panelapp_pin_matches_genokeeper(lab: dict[str, str]) -> None:
    """The version the script pins is the one `genokeeper.sources.panelapp` reads, its
    sha256 is a whole digest, and the step is one of the bootstrap's."""
    shown = _run(lab, 'echo "$panelapp_version"; echo "$panelapp_sha256"')
    version, sha = shown.stdout.splitlines()
    assert version == PANEL_VERSION
    assert len(sha) == 64
    assert "  step_panelapp\n" in SCRIPT.read_text()


def test_gencc_filed_only_when_changed(lab: dict[str, str], tmp_path: Path) -> None:
    """GenCC's export is fetched whole (never resumed: it is built per request), checked
    by genokeeper's reader and filed under its newest submission run date only when it
    differs from the newest copy; a prefetched copy with its pin is moved, not
    downloaded; a page the reader refuses is not recorded."""
    env, serve = _serving(lab, tmp_path)
    home = Path(lab["GENOKEEPER_HOME"])
    served = serve / GENCC_URL.rsplit("/", 1)[1]
    head = GENCC_HEAD.read_text()
    served.write_text(head)
    listed = _run(env, "check=1; step_gencc")
    assert listed.stdout == f"missing: GenCC gene-disease submissions ({GENCC_URL})\n"
    first = _refresh(env, "gencc")
    folder = home / "gencc"
    assert (first.returncode, first.stdout) == (
        0,
        f"refreshed: gencc {folder / '20261025'}\n",
    ), first.stderr
    dest = folder / "20261025" / "gencc-submissions.tsv"
    assert dest.read_bytes() == served.read_bytes()
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [f"gencc\t{dest}\t{_sha256(served)}"]
    (call,) = Path(env["CURL_LOG"]).read_text().splitlines()
    assert call.split()[-1] == GENCC_URL and " -C " not in f" {call} "
    same = _refresh(env, "gencc")
    assert (same.returncode, same.stdout) == (0, first.stdout)
    assert _folders(folder) == ["20261025"]
    served.write_text(head + _gencc_row("GENEX", "Strong", "2026-10-28"))
    new = _refresh(env, "gencc")
    assert new.stdout == f"refreshed: gencc {folder / '20261028'}\n"
    assert _folders(folder) == ["20261025", "20261028"]
    # A prefetched copy with the pinned size and sha256 is moved into place.
    shutil.rmtree(folder)
    (home / "bootstrap-manifest.tsv").unlink()
    prefetch = home / "prefetch" / "3f" / "gencc-submissions.tsv"
    prefetch.parent.mkdir(parents=True)
    prefetch.write_text(head)
    pin = f"gencc_size={prefetch.stat().st_size}; gencc_sha256={_sha256(prefetch)}"
    moved = _run(env, f"{pin}; step_gencc")
    assert moved.returncode == 0, moved.stderr
    assert dest.read_text() == head and not prefetch.exists()
    assert len(Path(env["CURL_LOG"]).read_text().splitlines()) == 3
    shutil.rmtree(folder)
    (home / "bootstrap-manifest.tsv").unlink()
    served.write_text("<html>maintenance</html>\n")
    refused = _run(env, "step_gencc")
    assert refused.returncode == 1
    assert "is not GenCC's submissions export" in refused.stderr
    assert not (home / "bootstrap-manifest.tsv").exists()


def _bgzipped(path: Path, text: str) -> Path:
    plain = path.with_name(f"{path.name}.plain")
    plain.write_text(text)
    pysam.tabix_compress(str(plain), str(path), force=True)
    plain.unlink()
    return path


# `genokeeper_python -m genokeeper.annotate.vep nmd-bed BB OUT` converts nothing: it names the
# bigBed it was given and writes a BED line; every other genokeeper_python call runs.
NMD_STUB = (
    'eval "real_$(declare -f genokeeper_python)"; genokeeper_python() {'
    ' if [ "$3" = nmd-bed ]; then echo "$4" >> "$BUILD_LOG";'
    ' printf \'chr1\\t10\\t20\\tGENEA\\n\' > "$5"; else real_genokeeper_python "$@"; fi; }'
)


@dataclass(frozen=True)
class VepDownloads:
    env: dict[str, str]
    pins: str  # the shell assignment of `vep_downloads`, naming the fakes' sha256s
    serve: Path
    prefetch: Path
    made: Path  # every stand-in, with its pinned bytes


@pytest.fixture
def vep_downloads(lab: dict[str, str], tmp_path: Path) -> VepDownloads:
    """Small stand-ins for slice 1e's seven downloads, as their publishers shape them: the
    VEP cache and pext under prefetch/1e/ with the pinned bytes, REVEL there with other
    bytes, the rest served by a curl that resumes."""
    home = Path(lab["GENOKEEPER_HOME"])
    made = tmp_path / "made"
    made.mkdir()
    cache = io.BytesIO()
    with tarfile.open(fileobj=cache, mode="w:gz") as tar:
        data = b"species\thomo_sapiens\n"
        info = tarfile.TarInfo("homo_sapiens/116_GRCh38/info.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    (made / "homo_sapiens_vep_116_GRCh38.tar.gz").write_bytes(cache.getvalue())
    revel = io.BytesIO()
    with zipfile.ZipFile(revel, "w") as z:
        z.writestr(
            "revel_with_transcript_ids",
            "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid\n"
            "1,900,950,G,A,T,M,0.9,ENST2\n"
            "1,35142,.,G,C,T,R,0.035,ENST1\n"
            "1,800,850,G,C,T,R,0.5,ENST2\n",
        )
    (made / "revel-v1.3_all_chromosomes.zip").write_bytes(revel.getvalue())
    _bgzipped(
        made / "AlphaMissense_hg38.tsv.gz",
        "# Copyright 2023 DeepMind Technologies Limited\n#\n# CC BY-NC-SA 4.0\n"
        "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant"
        "\tam_pathogenicity\tam_class\n"
        "chr1\t69094\tG\tT\thg38\tQ8NH21\tENST00000335137.4\tV2L\t0.2937\tlikely_benign\n",
    )
    sv = _bgzipped(
        made / "SpliceVault_data_GRCh38.tsv.gz",
        "#chrom\tpos\tref\talt\ttranscript_id\tsplicevault_site_type\n"
        "chr1\t925770\tG\tGT\tENST00000342066\tDonor_loss\n",
    )
    pysam.tabix_index(str(sv), seq_col=0, start_col=1, end_col=1, force=True)
    (made / "gnomad.pext.gtex_v10.base_level.tsv.gz").write_bytes(
        gzip.compress(
            b"locus\tgene_id\tgene_symbol\texp_prop_mean\tAdipose\n"
            b"chr1:100\tENSG1\tGENEA\t2.0000e-01\t0\n"
            b"chr1:101\tENSG1\tGENEA\t8.0000e-01\t0\n"
            b"chr1:102\tENSG1\tGENEA\tNaN\t0\n"
            b"chr1:102\tENSG2\tGENEB\t1.0000e-01\t0\n"
        )
    )
    (made / "nmdEscMane.bb").write_bytes(b"bigBed bytes")
    subfolders = {
        "homo_sapiens_vep_116_GRCh38.tar.gz": "vep",
        "revel-v1.3_all_chromosomes.zip": "revel",
        "AlphaMissense_hg38.tsv.gz": "alphamissense",
        "SpliceVault_data_GRCh38.tsv.gz": "splicevault",
        "SpliceVault_data_GRCh38.tsv.gz.tbi": "splicevault",
        "gnomad.pext.gtex_v10.base_level.tsv.gz": "pext",
        "nmdEscMane.bb": "nmd",
    }
    serve = tmp_path / "serve"
    serve.mkdir()
    prefetch = home / "prefetch" / "1e"
    lines: list[str] = []
    for name, sub in subfolders.items():
        data = (made / name).read_bytes()
        pin = hashlib.sha256(data).hexdigest()
        lines.append(f"{name} {sub} {len(data)} {pin} https://example.test/{name}")
        if sub in ("vep", "pext", "revel"):
            (prefetch / sub).mkdir(parents=True)
            # REVEL's prefetched copy was cut short: it is downloaded instead.
            (prefetch / sub / name).write_bytes(data[:10] if sub == "revel" else data)
        if sub != "vep" and sub != "pext":
            (serve / name).write_bytes(data)
    bin_dir = tmp_path / "vep-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
        "BUILD_LOG": str(tmp_path / "build.log"),
    }
    pins = "vep_downloads='" + "\n".join(lines) + "'"
    return VepDownloads(env, pins, serve, prefetch, made)


def test_vep_step_needs_flag(vep_downloads: VepDownloads) -> None:
    """Without --with-vep nothing is fetched and the step says what it would download;
    --check lists every file the run reads and changes nothing."""
    v = vep_downloads
    home = Path(v.env["GENOKEEPER_HOME"])
    files = vep_files(home)
    without = _run(v.env, f"{v.pins}; step_vep")
    assert without.returncode == 0, without.stderr
    assert (
        "VEP 116 cache and score files missing: run with --with-vep" in without.stderr
    )
    assert not Path(v.env["CURL_LOG"]).exists()
    listed = _run(v.env, f"{v.pins}; check=1; with_vep=1; step_vep")
    assert listed.returncode == 0, listed.stderr
    assert listed.stdout.splitlines() == [
        f"missing: VEP 116 file {p} (run with --with-vep)" for p in files.required()
    ]
    assert not Path(v.env["CURL_LOG"]).exists()
    assert not files.folder.exists()


def test_vep_files_moved_from_prefetch_or_fetched_and_pinned(
    vep_downloads: VepDownloads,
) -> None:
    """A prefetched copy with the pinned sha256 is moved into place, one with other bytes
    left and the file downloaded; then the cache is unpacked, REVEL sorted on GRCh38 and
    indexed, AlphaMissense indexed, pext written as BED with each gene's largest value, and
    the NMD track converted. A rerun does nothing; a served file with other bytes is
    refused and nothing is recorded."""
    v = vep_downloads
    home = Path(v.env["GENOKEEPER_HOME"])
    files = vep_files(home)
    done = _run(v.env, f"{v.pins}; {NMD_STUB}; with_vep=1; step_vep")
    assert done.returncode == 0, done.stderr
    assert (files.folder / "cache/homo_sapiens/116_GRCh38/info.txt").is_file()
    assert not (v.prefetch / "vep").exists() or not any((v.prefetch / "vep").iterdir())
    assert (
        v.prefetch / "revel" / "revel-v1.3_all_chromosomes.zip"
    ).stat().st_size == 10
    fetched = [c.split()[-1] for c in Path(v.env["CURL_LOG"]).read_text().splitlines()]
    assert sorted(fetched) == sorted(
        f"https://example.test/{p.name}" for p in v.serve.iterdir()
    )
    with pysam.TabixFile(str(files.revel)) as revel:
        assert next(iter(revel.header)).startswith("#chr\thg19_pos\tgrch38_pos")
        rows = [r.split("\t")[2] for r in revel.fetch("1", 0, 10_000)]
    assert rows == ["850", "950"]  # sorted on GRCh38; no row without a GRCh38 position
    with pysam.TabixFile(str(files.alphamissense)) as am:
        assert len(list(am.fetch("chr1", 69000, 69100))) == 1
    with pysam.TabixFile(str(files.pext)) as pext:
        assert list(pext.fetch("chr1", 0, 1000)) == [
            "chr1\t99\t100\tENSG1\t2.0000e-01",
            "chr1\t100\t101\tENSG1\t8.0000e-01",
            "chr1\t101\t102\tENSG2\t1.0000e-01",
        ]
    assert files.pext_max.read_text().splitlines() == ["ENSG1\t0.8", "ENSG2\t0.1"]
    assert files.nmd.read_text() == "chr1\t10\t20\tGENEA\n"
    assert Path(v.env["BUILD_LOG"]).read_text().splitlines() == [
        str(files.folder / ".download" / "nmdEscMane.bb")
    ]
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert [row.split("\t")[1] for row in manifest] == [
        str(p) for p in files.required()
    ]
    assert not (files.folder / ".download").exists()
    calls = len(fetched)
    again = _run(v.env, f"{v.pins}; with_vep=1; step_vep; check=1; step_vep")
    assert (again.returncode, again.stdout) == (0, "")
    assert len(Path(v.env["CURL_LOG"]).read_text().splitlines()) == calls
    # A served file with other bytes: removed and refused, nothing recorded.
    shutil.rmtree(files.folder)
    (home / "bootstrap-manifest.tsv").unlink()
    for made in v.made.iterdir():
        shutil.copy(made, v.serve / made.name)
    (v.serve / "nmdEscMane.bb").write_bytes(b"other bytes")
    refused = _run(v.env, f"{v.pins}; {NMD_STUB}; with_vep=1; step_vep")
    assert refused.returncode == 1
    assert "sha256 mismatch for https://example.test/nmdEscMane.bb" in refused.stderr
    assert not (files.folder / ".download" / "nmdEscMane.bb.part").exists()
    assert not (home / "bootstrap-manifest.tsv").exists()


# --- slice 3a: Parascopy's data and HG002 ---------------------------------------------------

PARASCOPY_BASE = "https://zenodo.org/records/15019940/files"
HG002_TRUTH_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
    "HG002_NA24385_son/v5.0q"
)
HG002_TRUTH = (
    "HG002_GRCh38_v5.0q_smvar.vcf.gz",
    "HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi",
    "HG002_GRCh38_v5.0q_smvar.benchmark.bed",
)
# Serves $SERVE/<name> for a URL ending in <name> or <name>?download=1 (Zenodo's form).
FAKE_ZENODO_CURL = """#!/usr/bin/env bash
out=
for a; do [ "${prev:-}" = -o ] && out=$a; prev=$a; done
echo "$prev" >> "$CURL_LOG"
name=${prev##*/}; name=${name%%\\?*}
[ -f "$SERVE/$name" ] || exit 22
cp "$SERVE/$name" "$out"
"""


@dataclass(frozen=True)
class Slice3a:
    env: dict[str, str]
    pins: str
    serve: Path
    tarballs: tuple[str, str]
    fetched: (
        Path  # where the stubbed `genokeeper.calibrate.hg002 fetch` logs its libraries
    )

    @property
    def home(self) -> Path:
        return Path(self.env["GENOKEEPER_HOME"])

    def curl_calls(self) -> list[str]:
        log = Path(self.env["CURL_LOG"])
        return log.read_text().splitlines() if log.exists() else []

    def fetch_calls(self) -> list[str]:
        return self.fetched.read_text().split() if self.fetched.exists() else []


@pytest.fixture
def slice3a(lab: dict[str, str], tmp_path: Path) -> Slice3a:
    """A lab with a curl serving two small Parascopy tarballs and three HG002 truth files,
    the shell pins naming their md5s, and `fetch` of the HG002 module logged, not run."""
    serve = tmp_path / "serve"
    serve.mkdir()
    tarballs = ("GRCh38_v1.7.tar.gz", "models_GRCh38_1KGP_v1.7.tar.gz")
    pins_tar: list[str] = []
    for name in tarballs:
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            data = f"{name}\n".encode()
            info = tarfile.TarInfo(f"{name.removesuffix('.tar.gz')}/data.txt")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        (serve / name).write_bytes(buf.getvalue())
        pins_tar.append(f"{name}:{hashlib.md5(buf.getvalue()).hexdigest()}")
    pins_truth: list[str] = []
    for name in HG002_TRUTH:
        data = f"{name}\n".encode()
        (serve / name).write_bytes(data)
        pins_truth.append(f"{name}:{hashlib.md5(data).hexdigest()}")
    bin_dir = tmp_path / "slice3a-bin"
    _executable(bin_dir / "curl", FAKE_ZENODO_CURL)
    fetched = tmp_path / "fetched.log"
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    pins = (
        f"parascopy_files='{' '.join(pins_tar)}';"
        f" hg002_truth_files='{' '.join(pins_truth)}';"
        # Only `fetch` is stubbed: `layout` is the module's own answer.
        "genokeeper_python() { case \"$*\" in *' fetch '*) "
        f'echo "${{@: -1}}" >> "{fetched}" ;; *) genokeeper_run python "$@" ;; esac; }}'
    )
    return Slice3a(env, pins, serve, tarballs, fetched)


def test_parascopy_step_fetches_and_unpacks(slice3a: Slice3a) -> None:
    s = slice3a
    listed = _run(s.env, f"{s.pins}; check=1; step_parascopy")
    assert listed.returncode == 0, listed.stderr
    assert listed.stdout.splitlines() == [
        f"missing: Parascopy data {n} ({PARASCOPY_BASE}/{n}?download=1)"
        for n in s.tarballs
    ]
    assert s.curl_calls() == []
    done = _run(s.env, f"{s.pins}; step_parascopy")
    assert done.returncode == 0, done.stderr
    folder = s.home / "parascopy"
    for name in s.tarballs:
        assert (folder / name).read_bytes() == (s.serve / name).read_bytes()
        member = folder / name.removesuffix(".tar.gz") / "data.txt"
        assert member.read_text() == f"{name}\n"
        assert not (folder / f"{name}.part").exists()
    assert s.curl_calls() == [f"{PARASCOPY_BASE}/{n}?download=1" for n in s.tarballs]
    again = _run(s.env, f"{s.pins}; step_parascopy; check=1; step_parascopy")
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert len(s.curl_calls()) == 2


def test_parascopy_md5_mismatch_refuses(slice3a: Slice3a) -> None:
    s = slice3a
    (s.serve / s.tarballs[1]).write_bytes(b"other bytes")
    refused = _run(s.env, f"{s.pins}; step_parascopy")
    assert refused.returncode == 1
    assert f"md5 mismatch for {PARASCOPY_BASE}/{s.tarballs[1]}?download=1" in (
        refused.stderr
    )
    folder = s.home / "parascopy"
    assert not (folder / s.tarballs[1]).exists()
    assert not (folder / f"{s.tarballs[1]}.part").exists()
    assert not (s.home / "bootstrap-manifest.tsv").exists()


def test_hg002_without_flag_says_size(slice3a: Slice3a) -> None:
    s = slice3a
    without = _run(s.env, f"{s.pins}; step_hg002")
    assert without.returncode == 0, without.stderr
    lines = without.stderr.splitlines()
    assert len(lines) == 1
    assert "--with-hg002" in lines[0] and "300 GB" in lines[0]
    assert s.curl_calls() == [] and s.fetch_calls() == []
    assert not (s.home / "giab").exists()
    listed = _run(s.env, f"{s.pins}; check=1; step_hg002")
    assert listed.returncode == 0, listed.stderr
    out = listed.stdout.splitlines()
    assert out[:3] == [
        f"missing: GIAB HG002 v5.0q truth file {n} ({HG002_TRUTH_BASE}/{n})"
        for n in HG002_TRUTH
    ]
    assert [line.split(":")[1].strip() for line in out[3:]] == [
        "HG002 library pcrfree",
        "HG002 library pcr",
    ]
    assert s.curl_calls() == []


def test_hg002_skips_downsampled_library(slice3a: Slice3a) -> None:
    s = slice3a
    calibration = s.home / "calibration" / "hg002"
    (calibration / "pcrfree").mkdir(parents=True)
    (calibration / "pcrfree" / "aligned.json").write_text("{}")  # aligned, pair deleted
    pair = calibration / "pcr" / "downsampled"
    pair.mkdir(parents=True)
    for mate in (1, 2):
        (pair / f"SRR38517578_{mate}.fastq.gz").write_bytes(b"")
    done = _run(s.env, f"{s.pins}; with_hg002=1; step_hg002")
    assert done.returncode == 0, done.stderr
    assert s.fetch_calls() == []  # both libraries done: nothing streamed
    truth = s.home / "giab" / "HG002-v5.0q-GRCh38"
    for name in HG002_TRUTH:
        assert (truth / name).read_bytes() == (s.serve / name).read_bytes()
    assert sorted(c.rsplit("/", 1)[1] for c in s.curl_calls()) == sorted(HG002_TRUTH)
    # A library missing its pair is streamed; the finished truth is not fetched again.
    (calibration / "pcrfree" / "aligned.json").unlink()
    again = _run(s.env, f"{s.pins}; with_hg002=1; step_hg002")
    assert again.returncode == 0, again.stderr
    assert s.fetch_calls() == ["pcrfree"]
    assert len(s.curl_calls()) == 3


def test_hg002_step_reads_only_hg002_rows(slice3a: Slice3a, tmp_path: Path) -> None:
    """The layout lists HG003's and HG004's records too: `--check` and `fetch` see only
    HG002's two libraries."""
    s = slice3a
    gone = tmp_path / "nothing"
    rows = [
        f"pcrfree\tHG002\t{gone}/a.json\t{gone}/p1\t{gone}/p2",
        f"pcr\tHG002\t{gone}/b.json\t{gone}/q1\t{gone}/q2",
        f"hg003\tHG003\t{gone}/c.json\t-\t-",
        f"hg004\tHG004\t{gone}/d.json\t{gone}/r1\t{gone}/r2",
    ]
    env = {**s.env, "LAYOUT": "\n".join(rows)}
    fake = (
        "genokeeper_python() { case \"$*\" in *' fetch '*) "
        f'echo "${{@: -1}}" >> "{s.fetched}" ;;'
        ' *" layout"*) printf \'%s\\n\' "$LAYOUT" ;; *) exit 9 ;; esac; }'
    )
    listed = _run(env, f"{s.pins}; {fake}; check=1; step_hg002")
    assert listed.returncode == 0, listed.stderr
    assert [line.split(":")[1].strip() for line in listed.stdout.splitlines()[3:]] == [
        "HG002 library pcrfree",
        "HG002 library pcr",
    ]
    done = _run(env, f"{s.pins}; {fake}; with_hg002=1; step_hg002")
    assert done.returncode == 0, done.stderr
    assert s.fetch_calls() == ["pcrfree", "pcr"]


@pytest.mark.parametrize("step", ["step_hg002", "step_hg004"])
def test_a_library_fetch_before_a_run_lets_the_bootstrap_go_on(
    lab: dict[str, str], tmp_path: Path, step: str
) -> None:
    """Before `genokeeper run` recorded the bases, `fetch` exits `NOT_YET`: the step says the
    libraries wait for that run and returns, so the steps after it (the FAF table, which the
    run needs) still run; any other failure still stops the bootstrap."""
    gone = tmp_path / "nothing"
    rows = "\n".join(
        [
            f"pcrfree\tHG002\t{gone}/a.json\t{gone}/p1\t{gone}/p2",
            f"hg004\tHG004\t{gone}/d.json\t{gone}/r1\t{gone}/r2",
            f"hg004_truth\t{gone}/truth",
        ]
    )
    env = {**lab, "LAYOUT": rows}
    # The truth files count as recorded: only the library fetch is under test.
    fake = (
        "step_done() { return 0; }; "
        "genokeeper_python() { case \"$*\" in *' fetch '*) echo refused >&2; return $CODE ;;"
        ' *" layout"*) printf \'%s\\n\' "$LAYOUT" ;; *) exit 9 ;; esac; }'
    )
    flags = "with_hg002=1; with_hg004=1"
    waited = _run(env | {"CODE": str(NOT_YET)}, f"{fake}; {flags}; {step}; echo after")
    assert waited.returncode == 0, waited.stderr
    assert waited.stdout.splitlines() == ["after"]
    assert "genokeeper run" in waited.stderr.splitlines()[-1]
    failed = _run(env | {"CODE": "2"}, f"{fake}; {flags}; {step}; echo after")
    assert failed.returncode != 0 and "after" not in failed.stdout


# --- stage 4: --refresh STEP, LitVar2 and ClinVar citations ---------------------------------

REFRESH_STEPS = (
    "clinvar clingen clingen_dosage clingen_regions cpic pharmgkb gwas medlineplus litvar"
    " clinvar_citations gencc"
)
# FAKE_LOGGING_CURL that answers -D with a header block and, given -R, dates the file as
# the server's Last-Modified ($LAST_MODIFIED, UTC), as curl does.
FAKE_DATED_CURL = """#!/usr/bin/env bash
out=; hdr=
for a; do
  [ "${prev:-}" = -o ] && out=$a
  [ "${prev:-}" = -D ] && hdr=$a
  prev=$a
done
echo "$*" >> "$CURL_LOG"
src=$SERVE/${prev##*/}
[ -f "$src" ] || exit 22
cp "$src" "$out"
if [ -n "$hdr" ]; then
  printf 'HTTP/1.1 200 OK\\r\\n' > "$hdr"
  [ -z "${LAST_MODIFIED:-}" ] || printf 'Last-Modified: Sun, 25 Oct 2026 10:38:43 GMT\\r\\n' >> "$hdr"
  printf '\\r\\n' >> "$hdr"
fi
case " $* " in *" -R "*) [ -z "${LAST_MODIFIED:-}" ] || TZ=UTC touch -t "$LAST_MODIFIED" "$out" ;; esac
"""
LITVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz"
CITATIONS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"
)


def _refresh(env: dict[str, str], *args: str) -> subprocess.CompletedProcess[str]:
    """`bootstrap.sh --refresh ARGS`, run as `genokeeper refresh` runs it."""
    return subprocess.run(
        ["bash", str(SCRIPT), "--refresh", *args],
        env=env,
        capture_output=True,
        text=True,
        check=False,
        stdin=subprocess.DEVNULL,
    )


def _folders(parent: Path) -> list[str]:
    return sorted(p.name for p in parent.iterdir() if not p.name.startswith("."))


def test_refresh_flag_runs_one_step_past_step_done(
    lab: dict[str, str], tmp_path: Path
) -> None:
    env, serve = _serving(lab, tmp_path)
    served = serve / "summaryAnnotations.zip"
    served.write_bytes(zip_bytes([Ann("1", "rs1", "TNF", "2B")], created="2026-09-30"))
    home = Path(lab["GENOKEEPER_HOME"])
    log = Path(env["CURL_LOG"])
    assert _run(env, "step_pharmgkb").returncode == 0
    # The step is done: a plain run fetches nothing.
    assert _run(env, "step_pharmgkb").returncode == 0
    assert len(log.read_text().splitlines()) == 1
    served.write_bytes(zip_bytes([Ann("1", "rs1", "TNF", "2B")], created="2026-10-30"))
    done = _refresh(env, "pharmgkb")
    assert done.returncode == 0, done.stderr
    new = home / "pharmgkb" / "20261030"
    assert done.stdout == f"refreshed: pharmgkb {new}\n"
    assert _folders(home / "pharmgkb") == ["20260930", "20261030"]
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert manifest == [
        f"pharmgkb\t{new / 'summaryAnnotations.zip'}\t{_sha256(served)}"
    ]
    # Again, with the manifest now matching the disk: still fetched, still reported.
    again = _refresh(env, "pharmgkb")
    assert (again.returncode, again.stdout) == (0, done.stdout)
    assert len(log.read_text().splitlines()) == 3
    # Only the step asked for ran: no other source was fetched, no tool installed.
    assert {c.split()[-1] for c in log.read_text().splitlines()} == {PHARMGKB_URL}


def test_refresh_flag_failure_is_nonzero_and_keeps_the_previous_folder(
    lab: dict[str, str], tmp_path: Path
) -> None:
    env, serve = _serving(lab, tmp_path)
    served = serve / "summaryAnnotations.zip"
    served.write_bytes(zip_bytes([Ann("1", "rs1", "TNF", "2B")], created="2026-09-30"))
    home = Path(lab["GENOKEEPER_HOME"])
    assert _run(env, "step_pharmgkb").returncode == 0
    manifest = (home / "bootstrap-manifest.tsv").read_text()
    # The network is down (the offline curl of `lab`), then a page that is not the download.
    down = _refresh(lab, "pharmgkb")
    served.write_bytes(b"<html>maintenance</html>")
    page = _refresh(env, "pharmgkb")
    for failed in (down, page):
        assert failed.returncode != 0
        assert "refreshed:" not in failed.stdout
        assert _folders(home / "pharmgkb") == ["20260930"]
        assert (home / "bootstrap-manifest.tsv").read_text() == manifest
    assert f"{PHARMGKB_URL} is not the ClinPGx clinical annotations" in page.stderr


@pytest.mark.tools
def test_clinvar_rejected_by_the_reader_never_becomes_newest(
    lab: dict[str, str], tmp_path: Path
) -> None:
    """The md5 matches, the index reads, the ##fileDate is there, and still a release the
    run's own reader has no claim in is not filed: the previous one stays the newest."""
    serve = tmp_path / "serve"
    serve.mkdir()
    name = "clinvar_20261011.vcf.gz"
    bin_dir = tmp_path / "curl-bin"
    _executable(bin_dir / "curl", FAKE_SERVING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }

    def serve_release(info: dict[str, str]) -> Path:
        vcf = write_clinvar(serve / name, "GRCh38", [("1", 100, 1, "A", "G", info)])
        Path(f"{vcf}.tbi").replace(serve / f"{name}.tbi")
        md5 = hashlib.md5(vcf.read_bytes()).hexdigest()
        (serve / "clinvar.vcf.gz.md5").write_text(f"{md5}  /x/{name}\n")
        return vcf

    home = Path(lab["GENOKEEPER_HOME"])
    manifest = home / "bootstrap-manifest.tsv"
    unclassified = {"GENEINFO": "GENEA:3"}  # no CLNSIG: the reader keeps no claim of it
    serve_release(unclassified)
    refused = _run(env, "step_clinvar")
    assert refused.returncode == 1
    assert f"the ClinVar release {name} is not one genokeeper's reader takes" in (
        refused.stderr
    )
    assert not manifest.exists() and not list((home / "clinvar").glob("2*"))
    # A good release is filed; a later release the reader rejects does not replace it.
    good = serve_release({"CLNSIG": "Pathogenic", "GENEINFO": "GENEA:3"})
    assert _run(env, "step_clinvar").returncode == 0
    dest = home / "clinvar" / "20261012" / "clinvar_GRCh38.vcf.gz"
    kept = dest.read_bytes()
    assert kept == good.read_bytes()
    before = manifest.read_text()
    serve_release(unclassified)
    again = _refresh(env, "clinvar")
    assert again.returncode != 0 and "refreshed:" not in again.stdout
    assert (dest.read_bytes(), manifest.read_text()) == (kept, before)
    serve_release({"CLNSIG": "Pathogenic", "GENEINFO": "GENEA:3"})
    ok = _refresh(env, "clinvar")
    assert ok.returncode == 0, ok.stderr
    assert ok.stdout == f"refreshed: clinvar {dest.parent}\n"


def test_litvar_is_not_in_the_default_steps(lab: dict[str, str]) -> None:
    stub = (
        "for f in $(declare -F | awk '{print $3}' | grep '^step_' | grep -v '^step_done$');"
        ' do eval "$f() { echo $f; }"; done; '
    )

    def steps(*args: str) -> list[str]:
        result = _run(lab, f"{stub}main {' '.join(args)}")
        assert result.returncode == 0, result.stderr
        return result.stdout.split()

    default = steps()
    assert "step_litvar" not in default and "step_clinvar_citations" not in default
    assert "step_clinvar" in default and "step_gwas" in default
    both = steps("--with-publications")
    assert [s for s in both if "litvar" in s or "citations" in s] == [
        "step_litvar",
        "step_clinvar_citations",
    ]
    assert [s for s in both if s not in ("step_litvar", "step_clinvar_citations")] == (
        default
    )
    # --refresh runs that one step and no other, with or without the flag.
    assert steps("--refresh", "litvar") == ["step_litvar"]
    assert steps("--refresh", "cpic", "--with-publications") == ["step_cpic"]


def test_usage_lists_refresh(lab: dict[str, str]) -> None:
    refused = subprocess.run(
        ["bash", str(SCRIPT), "--bogus"],
        env=lab,
        capture_output=True,
        text=True,
        check=False,
    )
    assert refused.returncode == 64
    assert "[--with-publications]" in refused.stderr
    assert f"--refresh STEP   (STEP: {REFRESH_STEPS})" in refused.stderr
    # Every named step is a function of the script, and no other STEP is taken.
    for step in REFRESH_STEPS.split():
        assert _run(lab, f"declare -F step_{step}").returncode == 0, step
    for args in (
        ["nosuch"],
        ["clinvar clingen"],  # two names in one word of the list
        ["gnomad_faf"],
        ["tools"],
        [],
        ["cpic", "--check"],
        ["--check", "--refresh", "cpic"],
    ):
        flags = ["--refresh", *args] if args[:1] != ["--check"] else args
        result = subprocess.run(
            ["bash", str(SCRIPT), *flags],
            env=lab,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 64, args
        assert result.stdout == ""
        assert not (Path(lab["GENOKEEPER_HOME"]) / "bootstrap-manifest.tsv").exists()


@pytest.mark.parametrize(
    ("step", "folder", "file", "url", "fixture"),
    [
        ("litvar", "litvar", "litvar2_variants.json.gz", LITVAR_URL, LITVAR_FIXTURE),
        (
            "clinvar_citations",
            "clinvar-citations",
            "var_citations.txt",
            CITATIONS_URL,
            CITATIONS_FIXTURE,
        ),
    ],
)
def test_publication_sources_filed_under_their_last_modified_date(
    lab: dict[str, str],
    tmp_path: Path,
    step: str,
    folder: str,
    file: str,
    url: str,
    fixture: Path,
) -> None:
    """Fetched whole from one fixed URL (never resumed: LitVar2 is rebuilt in place), read
    by genokeeper.evidence.publications, filed under the Last-Modified date (UTC), and
    fetched only under --refresh or --with-publications."""
    env, serve = _serving(lab, tmp_path)
    _executable(tmp_path / "serving-bin" / "curl", FAKE_DATED_CURL)
    env["LAST_MODIFIED"] = "202610251038"
    served = serve / file
    served.write_bytes(fixture.read_bytes())
    home = Path(lab["GENOKEEPER_HOME"])
    log = Path(env["CURL_LOG"])
    listed = _run(env, f"check=1; step_{step}")
    assert listed.stdout.startswith("missing: ") and url in listed.stdout
    assert not log.exists()
    done = _refresh(env, step)
    assert done.returncode == 0, done.stderr
    dest = home / folder / "20261025"
    assert done.stdout == f"refreshed: {step} {dest}\n"
    assert (dest / file).read_bytes() == served.read_bytes()
    assert (home / "bootstrap-manifest.tsv").read_text().splitlines() == [
        f"{step}\t{dest / file}\t{_sha256(served)}"
    ]
    calls = log.read_text().splitlines()
    assert len(calls) == 1 and calls[0].split()[-1] == url
    assert " -C " not in f" {calls[0]} " and " -R " in f" {calls[0]} "
    assert not (home / folder / ".download").exists()
    # Done: a plain step fetches nothing; --refresh fetches again.
    plain = _run(env, f"step_{step}")
    assert (plain.returncode, plain.stdout) == (0, "")
    assert len(log.read_text().splitlines()) == 1
    assert _refresh(env, step).returncode == 0
    assert len(log.read_text().splitlines()) == 2
    # An error page, and an answer without a Last-Modified date, are filed under nothing.
    shutil.rmtree(home / folder)
    (home / "bootstrap-manifest.tsv").unlink()
    served.write_bytes(gzip.compress(b"<html>maintenance</html>\n"))
    if file.endswith(".txt"):
        served.write_text("<html>maintenance</html>\n")
    page = _refresh(env, step)
    assert page.returncode != 0 and page.stdout == ""
    served.write_bytes(fixture.read_bytes())
    env["LAST_MODIFIED"] = ""
    undated = _refresh(env, step)
    assert undated.returncode != 0 and "no Last-Modified date" in undated.stderr
    assert _folders(home / folder) == []
    assert not (home / "bootstrap-manifest.tsv").exists()


# --- slice 3b: gnomAD-SV and HG002's v5.0q SV truth ----------------------------------------

SV_NAMES = ("gnomad.v4.1.sv.sites.vcf.gz", "gnomad.v4.1.sv.sites.vcf.gz.tbi")
HG002_SV_NAMES = (
    "HG002_GRCh38_v5.0q_stvar.vcf.gz",
    "HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi",
    "HG002_GRCh38_v5.0q_stvar.benchmark.bed",
)


@dataclass(frozen=True)
class Slice3b:
    env: dict[str, str]
    pins: str  # the shell assignments of `sv_downloads` and `hg002_sv_downloads`
    serve: Path
    prefetch: Path  # <lab>/prefetch/3b, holding every stand-in with its pinned bytes

    @property
    def home(self) -> Path:
        return Path(self.env["GENOKEEPER_HOME"])

    def curl_calls(self) -> list[str]:
        log = Path(self.env["CURL_LOG"])
        return log.read_text().splitlines() if log.exists() else []


@pytest.fixture
def slice3b(lab: dict[str, str], tmp_path: Path) -> Slice3b:
    """Small stand-ins for the five downloads, all in prefetch/3b/ with their pinned bytes
    and also served by a curl that resumes."""
    home = Path(lab["GENOKEEPER_HOME"])
    prefetch = home / "prefetch" / "3b"
    prefetch.mkdir(parents=True)
    serve = tmp_path / "serve"
    serve.mkdir()
    rows: dict[str, list[str]] = {"sv": [], "hg002": []}
    for kind, names in (("sv", SV_NAMES), ("hg002", HG002_SV_NAMES)):
        for name in names:
            data = f"{name} stand-in\n".encode()
            (prefetch / name).write_bytes(data)
            (serve / name).write_bytes(data)
            pin = hashlib.sha256(data).hexdigest()
            rows[kind].append(f"{name} {len(data)} {pin} https://example.test/{name}")
    bin_dir = tmp_path / "slice3b-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    pins = (
        f"sv_downloads='{chr(10).join(rows['sv'])}';"
        f" hg002_sv_downloads='{chr(10).join(rows['hg002'])}'"
    )
    return Slice3b(env, pins, serve, prefetch)


def test_sv_moved_from_prefetch_without_flag(slice3b: Slice3b) -> None:
    """Both files in prefetch/3b/ with their pins: moved with no flag and no download,
    recorded under `sv`; a rerun and --check do nothing."""
    s = slice3b
    files = sv_files(s.home)
    stand_ins = [(s.prefetch / n).read_bytes() for n in SV_NAMES]
    done = _run(s.env, f"{s.pins}; step_sv")
    assert done.returncode == 0, done.stderr
    assert [p.read_bytes() for p in files.required()] == stand_ins
    assert not any((s.prefetch / n).exists() for n in SV_NAMES)
    assert s.curl_calls() == []
    manifest = (s.home / "bootstrap-manifest.tsv").read_text().splitlines()
    assert [row.split("\t")[:2] for row in manifest] == [
        ["sv", str(p)] for p in files.required()
    ]
    again = _run(s.env, f"{s.pins}; step_sv; check=1; step_sv")
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert s.curl_calls() == []


def test_sv_prefetch_wrong_sha_downloads_or_says(slice3b: Slice3b) -> None:
    """A prefetched file with other bytes: without --with-sv one `say` line names the flag
    and nothing is moved or fetched (--check lists both files); with it that file is
    downloaded and the good one still moved."""
    s = slice3b
    files = sv_files(s.home)
    (s.prefetch / SV_NAMES[0]).write_bytes(b"cut short")
    without = _run(s.env, f"{s.pins}; step_sv")
    assert without.returncode == 0, without.stderr
    (line,) = without.stderr.splitlines()
    assert "--with-sv" in line
    assert not files.gnomad.parent.exists() and s.curl_calls() == []
    assert (s.prefetch / SV_NAMES[0]).read_bytes() == b"cut short"
    assert (s.prefetch / SV_NAMES[1]).exists()
    listed = _run(s.env, f"{s.pins}; check=1; step_sv")
    assert listed.stdout.splitlines() == [
        f"missing: gnomAD-SV file {p} (run with --with-sv)" for p in files.required()
    ]
    done = _run(s.env, f"{s.pins}; with_sv=1; step_sv")
    assert done.returncode == 0, done.stderr
    assert files.gnomad.read_bytes() == (s.serve / SV_NAMES[0]).read_bytes()
    assert [c.split()[-1] for c in s.curl_calls()] == [
        f"https://example.test/{SV_NAMES[0]}"
    ]
    assert (s.prefetch / SV_NAMES[0]).read_bytes() == b"cut short"
    assert step_rows(s.home, "sv") == [str(p) for p in files.required()]


def step_rows(home: Path, step: str) -> list[str]:
    manifest = (home / "bootstrap-manifest.tsv").read_text().splitlines()
    return [r.split("\t")[1] for r in manifest if r.split("\t")[0] == step]


def test_sv_download_mismatch_refuses(slice3b: Slice3b) -> None:
    """--with-sv with nothing prefetched and a server serving other bytes: exit 1, the
    `.part` removed, nothing under the final name, nothing recorded."""
    s = slice3b
    files = sv_files(s.home)
    for n in SV_NAMES:
        (s.prefetch / n).unlink()
    (s.serve / SV_NAMES[0]).write_bytes(b"other bytes")
    refused = _run(s.env, f"{s.pins}; with_sv=1; step_sv")
    assert refused.returncode == 1
    assert f"sha256 mismatch for https://example.test/{SV_NAMES[0]}" in refused.stderr
    assert not files.gnomad.exists()
    assert not files.gnomad.with_name(f"{SV_NAMES[0]}.part").exists()
    assert not (s.home / "bootstrap-manifest.tsv").exists()


def test_hg002_sv_truth_from_prefetch(slice3b: Slice3b) -> None:
    """--with-hg002 moves the three v5.0q SV files into 3a's truth folder under
    `hg002_sv_truth`; a lab whose `hg002_truth` was recorded keeps that record; without
    the flag the step says so in one line, and --check lists the three."""
    s = slice3b
    truth = s.home / "giab" / "HG002-v5.0q-GRCh38"
    truth.mkdir(parents=True)
    small = truth / "HG002_GRCh38_v5.0q_smvar.vcf.gz"
    small.write_bytes(b"3a's truth")
    prior = _run(s.env, f'record hg002_truth "{small}"')
    assert prior.returncode == 0, prior.stderr
    before = (s.home / "bootstrap-manifest.tsv").read_text()
    stand_ins = [(s.prefetch / n).read_bytes() for n in HG002_SV_NAMES]
    without = _run(s.env, f"{s.pins}; step_hg002_sv")
    assert without.returncode == 0, without.stderr
    (line,) = without.stderr.splitlines()
    assert "--with-hg002" in line
    assert (s.prefetch / HG002_SV_NAMES[0]).exists() and s.curl_calls() == []
    listed = _run(s.env, f"{s.pins}; check=1; step_hg002_sv")
    assert [x.split(" (")[0] for x in listed.stdout.splitlines()] == [
        f"missing: GIAB HG002 v5.0q SV truth file {n}" for n in HG002_SV_NAMES
    ]
    done = _run(s.env, f"{s.pins}; with_hg002=1; step_hg002_sv")
    assert done.returncode == 0, done.stderr
    assert [(truth / n).read_bytes() for n in HG002_SV_NAMES] == stand_ins
    assert s.curl_calls() == []
    manifest = (s.home / "bootstrap-manifest.tsv").read_text()
    assert manifest.startswith(before)
    assert step_rows(s.home, "hg002_sv_truth") == [
        str(truth / n) for n in HG002_SV_NAMES
    ]
    again = _run(
        s.env, f"{s.pins}; with_hg002=1; step_hg002_sv; check=1; step_hg002_sv"
    )
    assert (again.returncode, again.stdout) == (0, ""), again.stderr


def test_sv_say_lines_give_the_pinned_sizes(lab: dict[str, str]) -> None:
    """With the script's own pins and nothing prefetched, the lines name the flags and the
    download sizes (1.74 GB of gnomAD-SV, 50 MB of HG002 SV truth)."""
    said = _run(lab, "step_sv; step_hg002_sv")
    assert said.returncode == 0, said.stderr
    sv, hg002 = said.stderr.splitlines()
    assert "--with-sv" in sv and "1.74 GB" in sv
    assert "--with-hg002" in hg002 and "50 MB" in hg002


# --- slice 3c: the masked GRCh38 and HG004's truth -----------------------------------------

MASKED_NAME = (
    "GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"
)
HG004_TRUTH = (
    "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
    "HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
)
# bwa as `masking index` runs it: `index -p PREFIX FASTA` writes PREFIX's five files.
FAKE_BWA = """#!/usr/bin/env bash
[ "$1" = index ] || exit 1
shift
while [ $# -gt 0 ]; do
  case "$1" in -p) prefix=$2; shift 2 ;; *) fasta=$1; shift ;; esac
done
for ext in amb ann bwt pac sa; do echo "$ext" > "$prefix.$ext"; done
echo "index $fasta" >> "$BWA_LOG"
"""


@dataclass(frozen=True)
class Slice3c:
    env: dict[str, str]
    pins: str  # the masked FASTA's pins, and `masking verify` stubbed to $VERIFY_EXIT
    prefetch: Path
    text: str  # the masked FASTA's text

    @property
    def home(self) -> Path:
        return Path(self.env["GENOKEEPER_HOME"])

    def logged(self, name: str) -> list[str]:
        log = Path(self.env[name])
        return log.read_text().splitlines() if log.exists() else []


@pytest.fixture
def slice3c(lab: dict[str, str], tmp_path: Path) -> Slice3c:
    """The unmasked FASTA recorded, the masked one gzipped in prefetch/3c/ with its pins,
    fake `bwa` and `rtg`, a logging curl serving nothing, and `masking verify` stubbed: it
    writes its BED and exits $VERIFY_EXIT (the base-level comparison is test_masking's)."""
    home = Path(lab["GENOKEEPER_HOME"])
    unmasked = unmasked_fasta(home)
    unmasked.parent.mkdir(parents=True)
    unmasked.write_text(">chr1\nACGTACGTAC\n")
    pysam.faidx(str(unmasked))
    text = ">chr1\nACNNNNGTAC\n"
    prefetch = home / "prefetch" / "3c"
    prefetch.mkdir(parents=True)
    packed = gzip.compress(text.encode())
    (prefetch / MASKED_NAME).write_bytes(packed)
    _executable(home / "tools" / "rtg-tools-3.13" / "rtg", FAKE_RTG)
    bin_dir = tmp_path / "slice3c-bin"
    _executable(bin_dir / "bwa", FAKE_BWA)
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    serve = tmp_path / "serve"
    serve.mkdir()
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
        "RTG_LOG": str(tmp_path / "rtg.log"),
        "BWA_LOG": str(tmp_path / "bwa.log"),
        "VERIFY_EXIT": "0",
    }
    pins = (
        f"masked_size={len(packed)}; masked_sha256={hashlib.sha256(packed).hexdigest()};"
        f' record grch38_fasta "{unmasked}";'
        " genokeeper_python() { case \"$*\" in *'masking verify '*)"
        ' printf "chr1\\t2\\t6\\n" > "${@: -1}"; return "$VERIFY_EXIT" ;;'
        ' *) genokeeper_run python "$@" ;; esac; }'
    )
    return Slice3c(env, pins, prefetch, text)


def _masked_outputs(home: Path) -> list[Path]:
    fasta = grch38_fasta(home)
    sdf = layout(home)["sdf"]
    return [
        fasta,
        Path(f"{fasta}.fai"),
        masked_bed(home),
        *(Path(f"{fasta}.{ext}") for ext in ("amb", "ann", "bwt", "pac", "sa")),
        sdf / "done",
        sdf / "fasta.sha256",
    ]


@pytest.mark.tools
def test_masked_reference_from_prefetch(slice3c: Slice3c) -> None:
    """Prefetched with its pins: moved, unpacked to `grch38_fasta`, indexed, its mask
    verified into `masked_bed`, bwa-indexed and SDF-built, every output recorded under
    `masked_reference`, with no flag and no download; a rerun and --check do nothing."""
    s = slice3c
    done = _run(s.env, f"{s.pins}; step_masked_reference")
    assert done.returncode == 0, done.stderr
    fasta = grch38_fasta(s.home)
    assert fasta.read_text() == s.text
    assert not (s.prefetch / MASKED_NAME).exists()
    assert masked_bed(s.home).read_text() == "chr1\t2\t6\n"
    assert s.logged("BWA_LOG") == [f"index {fasta}"]
    assert s.logged("RTG_LOG") == [f"format {fasta}"]
    sdf = layout(s.home)["sdf"]
    sha = hashlib.sha256(s.text.encode()).hexdigest()
    assert (sdf / "fasta.sha256").read_text() == f"{sha}\n"
    assert step_rows(s.home, "masked_reference") == [
        str(p) for p in _masked_outputs(s.home)
    ]
    assert s.logged("CURL_LOG") == []
    assert not list(fasta.parent.glob(".*"))  # no temporary left
    again = _run(
        s.env, f"{s.pins}; step_masked_reference; check=1; step_masked_reference"
    )
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert len(s.logged("RTG_LOG")) == 1


@pytest.mark.tools
def test_masked_reference_mask_differs_dies(slice3c: Slice3c) -> None:
    """`masking verify` refusing (the derived BED differs from the committed mask): exit 1,
    nothing recorded under `masked_reference`, no index or SDF built, and no FASTA or .fai
    at the path every consumer reads."""
    s = slice3c
    refused = _run({**s.env, "VERIFY_EXIT": "2"}, f"{s.pins}; step_masked_reference")
    assert refused.returncode == 1
    assert "mask" in refused.stderr
    assert step_rows(s.home, "masked_reference") == []
    assert s.logged("BWA_LOG") == [] and s.logged("RTG_LOG") == []
    fasta = grch38_fasta(s.home)
    assert not fasta.exists() and not Path(f"{fasta}.fai").exists()


def test_masked_reference_without_flag_says_size(lab: dict[str, str]) -> None:
    """Nothing prefetched, the script's own pins: one `say` line naming the flag and the
    886 MB download; --check lists the masked reference; nothing fetched."""
    home = Path(lab["GENOKEEPER_HOME"])
    said = _run(lab, "step_masked_reference")
    assert said.returncode == 0, said.stderr
    (line,) = said.stderr.splitlines()
    assert "--with-masked-reference" in line and "886 MB" in line
    listed = _run(lab, "check=1; step_masked_reference")
    assert listed.returncode == 0, listed.stderr
    (missing,) = listed.stdout.splitlines()
    assert str(grch38_fasta(home)) in missing
    assert not (home / "bootstrap-manifest.tsv").exists()
    pinned = _run(lab, 'echo "$masked_name $masked_size $masked_sha256"')
    assert pinned.stdout.split() == [MASKED_NAME, "885955190", MASKED_SHA256]


def test_hg004_without_flag_says_size(lab: dict[str, str]) -> None:
    """Nothing there, the script's own pins: one `say` line naming --with-hg004, the 163 MB
    of truth and the 144 GB from CNGB; --check lists the truth and the pair."""
    home = Path(lab["GENOKEEPER_HOME"])
    said = _run(lab, "step_hg004")
    assert said.returncode == 0, said.stderr
    (line,) = said.stderr.splitlines()
    assert "--with-hg004" in line and "163 MB" in line and "144 GB" in line
    listed = _run(lab, "check=1; step_hg004")
    assert listed.returncode == 0, listed.stderr
    out = listed.stdout.splitlines()
    assert [x.split(" (")[0] for x in out[:3]] == [
        f"missing: GIAB HG004 v4.2.1 truth file {n}" for n in HG004_TRUTH
    ]
    assert out[3].startswith("missing: HG004 library hg004")
    assert len(out) == 4
    assert not (home / "giab").exists()


def test_hg004_truth_from_its_prefetch(lab: dict[str, str], tmp_path: Path) -> None:
    """The three truth files in prefetch/3c/ with their pins: --with-hg004 moves them into
    giab/HG004-v4.2.1-GRCh38/ under `hg004_truth`, downloads nothing, then fetches HG004's
    pair; with its record there, the pair is not fetched again."""
    home = Path(lab["GENOKEEPER_HOME"])
    prefetch = home / "prefetch" / "3c"
    prefetch.mkdir(parents=True)
    rows = []
    for name in HG004_TRUTH:
        data = f"{name} stand-in\n".encode()
        (prefetch / name).write_bytes(data)
        pin = hashlib.sha256(data).hexdigest()
        rows.append(f"{name} {len(data)} {pin} https://example.test/{name}")
    fetched = tmp_path / "fetched.log"
    bin_dir = tmp_path / "hg004-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(tmp_path),
        "CURL_LOG": str(tmp_path / "curl.log"),
    }
    pins = (
        f"hg004_truth_downloads='{chr(10).join(rows)}';"
        " genokeeper_python() { case \"$*\" in *' fetch '*) "
        f'echo "${{@: -1}}" >> "{fetched}" ;; *) genokeeper_run python "$@" ;; esac; }}'
    )
    done = _run(env, f"{pins}; with_hg004=1; step_hg004")
    assert done.returncode == 0, done.stderr
    truth = home / "giab" / "HG004-v4.2.1-GRCh38"
    assert [(truth / n).read_text() for n in HG004_TRUTH] == [
        f"{n} stand-in\n" for n in HG004_TRUTH
    ]
    assert not any((prefetch / n).exists() for n in HG004_TRUTH)
    assert not Path(env["CURL_LOG"]).exists()
    assert step_rows(home, "hg004_truth") == [str(truth / n) for n in HG004_TRUTH]
    assert fetched.read_text().split() == ["hg004"]
    aligned = home / "calibration" / "hg004" / "pcrfree" / "aligned.json"
    aligned.parent.mkdir(parents=True)
    aligned.write_text("{}")
    again = _run(env, f"{pins}; with_hg004=1; step_hg004; check=1; step_hg004")
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert fetched.read_text().split() == ["hg004"]


# --- slice 3f: HHV-6, hmmix and the panel-derived typing inputs ------------------------------

HHV6_PARTS = ("NC_001664.4", "NC_000898.1")


@dataclass(frozen=True)
class Slice3f:
    env: dict[str, str]
    pins: (
        str  # shell: `typing_downloads` and stand-ins for bwa and the two panel builds
    )
    serve: Path
    prefetch: Path  # <lab>/prefetch/3f, in the layout the script reads
    files: TypingFiles
    stand_ins: dict[Path, bytes]  # destination -> bytes, for every downloaded file

    @property
    def home(self) -> Path:
        return Path(self.env["GENOKEEPER_HOME"])

    def curl_calls(self) -> list[str]:
        log = Path(self.env["CURL_LOG"])
        return log.read_text().splitlines() if log.exists() else []

    def builds(self) -> list[str]:
        log = Path(self.env["BUILD_LOG"])
        return log.read_text().splitlines() if log.exists() else []

    def prefetched(self, destination: Path) -> Path:
        return _prefetched_3f(self.prefetch, destination)


def _prefetched_3f(prefetch: Path, destination: Path) -> Path:
    """Where prefetch/3f/ keeps the stand-in for `destination`: hmmix's helpers in hmmix/,
    the archaic BCFs in archaic/, the genomes at the top."""
    if destination.parent.name == "archaic":
        return prefetch / "archaic" / destination.name
    if destination.parent.name == "0.9.2":
        return prefetch / "hmmix" / destination.name
    return prefetch / destination.name


@pytest.fixture
def slice3f(lab: dict[str, str], tmp_path: Path) -> Slice3f:
    """Small stand-ins for the 50 downloads, all in prefetch/3f/ with their pinned bytes and
    also served by a curl that resumes; the bwa index and the two panel builds are stubs that
    write the files `genokeeper.typing_files` names and log a line."""
    home = Path(lab["GENOKEEPER_HOME"])
    files = typing_files(home)
    prefetch = home / "prefetch" / "3f"
    serve = tmp_path / "serve"
    serve.mkdir()
    destinations = [
        *(files.hhv6_fasta.with_name(f"{name}.fa") for name in HHV6_PARTS),
        *(files.hmmix / name for name in HMMIX_HELPERS),
        *(files.hmmix / "archaic" / name for name in ARCHAIC),
    ]
    stand_ins: dict[Path, bytes] = {}
    rows = []
    for dest in destinations:
        data = f">{dest.name}\nACGT\n".encode()
        stand_ins[dest] = data
        (serve / dest.name).write_bytes(data)
        source = _prefetched_3f(prefetch, dest)
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes(data)
        pin = hashlib.sha256(data).hexdigest()
        rows.append(
            f"{dest.relative_to(home)} {len(data)} {pin} https://example.test/{dest.name}"
        )
    required = files.required()
    index = required[1:6]
    panel = [p for p in required if p.parent == files.panel_ref]
    af = [p for p in required if p.parent == files.roh_af]
    (home / "pgs").mkdir()
    (home / "pgs" / PANEL_NAME).write_bytes(b"panel")
    bin_dir = tmp_path / "slice3f-bin"
    _executable(bin_dir / "curl", FAKE_RESUMING_CURL)

    def stub(name: str, made: list[Path]) -> str:
        parents = " ".join(f"'{p.parent}'" for p in made[:1])
        touch = " ".join(f"'{p}'" for p in made)
        return f'{name}() {{ echo {name} >> "$BUILD_LOG"; mkdir -p {parents}; touch {touch}; }}'

    pins = (
        f"typing_downloads='{chr(10).join(rows)}'; "
        + stub("typing_bwa_index", index)
        + "; "
        + stub("typing_build_panel", panel)
        + "; "
        + stub("typing_build_af", af)
    )
    env = {
        **lab,
        "PATH": f"{bin_dir}:{lab['PATH']}",
        "SERVE": str(serve),
        "CURL_LOG": str(tmp_path / "curl.log"),
        "BUILD_LOG": str(tmp_path / "build.log"),
    }
    return Slice3f(env, pins, serve, prefetch, files, stand_ins)


def test_typing_moved_from_prefetch_without_flag(slice3f: Slice3f) -> None:
    """Every file in prefetch/3f/ with its pin: moved with no flag and no download, the two
    genomes concatenated into hhv6.fa (indexed), the panel builds run, all recorded under
    `typing`; a rerun and --check do nothing."""
    s = slice3f
    done = _run(s.env, f"{s.pins}; step_typing")
    assert done.returncode == 0, done.stderr
    for dest, data in s.stand_ins.items():
        assert dest.read_bytes() == data
        assert not s.prefetched(dest).exists()
    assert s.curl_calls() == []
    assert s.files.hhv6_fasta.read_bytes() == b"".join(
        s.stand_ins[s.files.hhv6_fasta.with_name(f"{name}.fa")] for name in HHV6_PARTS
    )
    assert s.builds() == [
        "typing_bwa_index",
        "typing_build_panel",
        "typing_build_af",
    ]
    parts = [s.files.hhv6_fasta.with_name(f"{name}.fa") for name in HHV6_PARTS]
    assert step_rows(s.home, "typing") == [
        str(p) for p in [*s.files.required(), *parts]
    ]
    again = _run(s.env, f"{s.pins}; step_typing; check=1; step_typing")
    assert (again.returncode, again.stdout) == (0, ""), again.stderr
    assert s.curl_calls() == [] and len(s.builds()) == 3


def test_typing_wrong_sha_says_or_downloads(slice3f: Slice3f) -> None:
    """A prefetched file with other bytes: without --with-typing one `say` line names the
    flag and nothing is moved, fetched or built (--check lists every typing file); with it
    that file is downloaded and the others still moved."""
    s = slice3f
    cut = s.files.hmmix / "archaic" / "individuals_highcov.7.bcf"
    (s.prefetched(cut)).write_bytes(b"cut short")
    without = _run(s.env, f"{s.pins}; step_typing")
    assert without.returncode == 0, without.stderr
    (line,) = without.stderr.splitlines()
    assert "--with-typing" in line and "GB" in line
    assert not s.files.hhv6_fasta.parent.exists() and s.curl_calls() == []
    assert s.builds() == []
    assert s.prefetched(cut).read_bytes() == b"cut short"
    assert s.prefetched(s.files.hmmix / HMMIX_HELPERS[0]).exists()
    listed = _run(s.env, f"{s.pins}; check=1; step_typing")
    assert listed.stdout.splitlines() == [
        f"missing: typing file {p} (run with --with-typing)" for p in s.files.required()
    ]
    done = _run(s.env, f"{s.pins}; with_typing=1; step_typing")
    assert done.returncode == 0, done.stderr
    assert cut.read_bytes() == s.stand_ins[cut]
    assert [c.split()[-1] for c in s.curl_calls()] == [
        "https://example.test/individuals_highcov.7.bcf"
    ]
    assert s.prefetched(cut).read_bytes() == b"cut short"
    assert set(step_rows(s.home, "typing")) >= {str(p) for p in s.files.required()}


def test_typing_download_mismatch_refuses(slice3f: Slice3f) -> None:
    """--with-typing with nothing prefetched and a server serving other bytes: exit 1, the
    `.part` removed, nothing under the final name, nothing built or recorded."""
    s = slice3f
    for dest in s.stand_ins:
        s.prefetched(dest).unlink()
    first = s.files.hhv6_fasta.with_name(f"{HHV6_PARTS[0]}.fa")
    (s.serve / first.name).write_bytes(b"other bytes")
    refused = _run(s.env, f"{s.pins}; with_typing=1; step_typing")
    assert refused.returncode == 1
    assert f"sha256 mismatch for https://example.test/{first.name}" in refused.stderr
    assert not first.exists() and not first.with_name(f"{first.name}.part").exists()
    assert not s.files.hhv6_fasta.exists() and s.builds() == []
    assert not (s.home / "bootstrap-manifest.tsv").exists()


def test_typing_pins_name_the_layouts_downloads(lab: dict[str, str]) -> None:
    """The script's own 50 pins cover the two genomes and every hmmix file the layout
    requires, nothing else; with nothing prefetched the `say` line names the flag and the
    3.5 GB it downloads."""
    home = Path(lab["GENOKEEPER_HOME"])
    files = typing_files(home)
    listed = _run(lab, 'echo "$typing_downloads"; echo "$hhv6_parts"')
    assert listed.returncode == 0, listed.stderr
    *rows, parts = listed.stdout.splitlines()
    assert parts.split() == list(HHV6_PARTS)
    destinations = {home / row.split()[0] for row in rows}
    assert destinations == {
        *(files.hhv6_fasta.with_name(f"{name}.fa") for name in HHV6_PARTS),
        *(p for p in files.required() if p.is_relative_to(files.hmmix)),
    }
    for row in rows:
        _, size, pin, url = row.split()
        assert size.isdigit() and len(pin) == 64 and url.startswith("https://")
    assert sum(int(row.split()[1]) for row in rows) == 3_496_325_138
    said = _run(lab, "step_typing")
    assert said.returncode == 0, said.stderr
    (line,) = said.stderr.splitlines()
    assert "--with-typing" in line and "3.5 GB" in line
