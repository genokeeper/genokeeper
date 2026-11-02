"""scripts/check_no_genome_data.sh in a throwaway git repository: what it refuses to commit;
and the project settings that deny the agent's reads of genome data."""

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_no_genome_data.sh"
# An invented kit ID, and the private-ID pattern each test repository's guard reads from its
# git directory's info/genome-private-ids (the user's real one is never committed).
KIT_ID = "TESTKIT" + "0001"
PRIVATE_IDS = "# a comment line\nTESTKIT[0-9]+\n"
# A VCF's column header up to its first sample column.
VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
# A result folder below the top level, and the genotype formats the pipeline handles.
NESTED = (
    "sub/reports/owner.tsv",
    "a/data/raw/x.txt",
    "x/results/y.txt",
    "x/work/y.txt",
)
FORMATS = (
    "panel.pgen",
    "panel.psam",
    "panel.pvar.zst",
    "owner.bed",
    "owner.bim",
    "owner.fam",
    "owner.vcf.bgz",
    "owner.vcf.zst",
    "reads.fastq.bz2",
    "reads.fq.bz2",
    "store.sqlite3",
    "store.duckdb.wal",
)


def _git(repo: Path, *args: str) -> None:
    env = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, env=env
    )


def _init(repo: Path, private_ids: str | None = PRIVATE_IDS) -> None:
    _git(repo, "init", "-q")
    if private_ids is not None:
        (repo / ".git" / "info" / "genome-private-ids").write_text(private_ids)


def _guard(tmp_path: Path, name: str, text: str = "x\n") -> dict[str, int]:
    """Stage one file `name` in a new repository and run the guard both ways."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo)
    path = repo / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    _git(repo, "add", "-f", "--", name)
    return {
        mode: subprocess.run(
            ["bash", str(SCRIPT), mode], cwd=repo, capture_output=True, check=False
        ).returncode
        for mode in ("--staged", "--tracked")
    }


@pytest.mark.parametrize(
    "name",
    [
        "sample.vcf.gz",
        "X.VCF.GZ",  # macOS writes and git tracks any case
        "calls/Sample.Vcf",
        "sample.vcf.gz.tbi",
        "sample.vcf.gz.csi",
        "panel.tsv.gz.TBI",
        "regions.bed.gz",
        "reads.CRAM",
        "store.DuckDB",
        "Data/raw/notes.txt",
        f"notes/{KIT_ID}.txt",  # the kit ID in a path, not only in text
        *NESTED,
        *FORMATS,
    ],
)
def test_a_genomic_or_personal_path_is_refused(tmp_path: Path, name: str) -> None:
    assert _guard(tmp_path, name) == {"--staged": 1, "--tracked": 1}


def test_gitignore_ignores_the_genotype_formats(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo)
    shutil.copy(SCRIPT.parents[1] / ".gitignore", repo / ".gitignore")
    shown = subprocess.run(
        ["git", "-C", str(repo), "check-ignore", "--no-index", *NESTED, *FORMATS],
        capture_output=True,
        text=True,
        check=False,
    )
    assert shown.stdout.splitlines() == [*NESTED, *FORMATS]


def test_the_kit_id_in_text_is_refused(tmp_path: Path) -> None:
    assert _guard(tmp_path, "tests/fixtures/notes.txt", f"kit {KIT_ID}\n") == {
        "--staged": 1,
        "--tracked": 1,
    }


def test_without_the_private_id_file_an_id_passes(tmp_path: Path) -> None:
    """A fresh clone or CI has no private-ID file: the ID rule is off, the others hold."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo, private_ids=None)
    (repo / "notes.md").write_text(f"kit {KIT_ID}\n")
    _git(repo, "add", "notes.md")
    done = subprocess.run(
        ["bash", str(SCRIPT), "--staged"], cwd=repo, capture_output=True, check=False
    )
    assert done.returncode == 0


def test_the_kit_id_in_a_binary_file_is_refused(tmp_path: Path) -> None:
    """A NUL byte makes git's text diff print "Binary files differ" and no line of it."""
    assert _guard(tmp_path, "tests/fixtures/report.dat", f"kit {KIT_ID}\0binary\n") == {
        "--staged": 1,
        "--tracked": 1,
    }


def test_a_symlink_replaced_by_a_file_is_checked(tmp_path: Path) -> None:
    """A tracked symlink turned into a regular file is a type change (`T`), not a modify."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo)
    link = repo / "tests" / "fixtures" / "notes.txt"
    link.parent.mkdir(parents=True)
    link.symlink_to("elsewhere")
    _git(repo, "add", "-f", "--", "tests/fixtures/notes.txt")
    _git(repo, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "link")
    link.unlink()
    link.write_text(f"kit {KIT_ID}\n")
    _git(repo, "add", "-f", "--", "tests/fixtures/notes.txt")
    staged = subprocess.run(
        ["bash", str(SCRIPT), "--staged"], cwd=repo, capture_output=True, check=False
    )
    assert staged.returncode == 1
    _git(repo, "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "file")
    assert _range(repo, "HEAD~1..HEAD") == 1


def _range(repo: Path, spec: str) -> int:
    return subprocess.run(
        ["bash", str(SCRIPT), "--range", spec],
        cwd=repo,
        capture_output=True,
        check=False,
    ).returncode


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("data/raw/sample.vcf.gz", "x\n"),
        ("tests/fixtures/report.dat", f"kit {KIT_ID}\0binary\n"),
        ("tests/fixtures/eh.vcf.txt", f"{VCF_HEADER}\tsample\n"),
    ],
)
def test_a_range_refuses_what_a_later_commit_deleted(
    tmp_path: Path, name: str, text: str
) -> None:
    """Each commit is checked, not only the tip, whose tree no longer holds the file."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo)
    commit = ("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q")
    _git(repo, *commit, "--allow-empty", "-m", "base")
    _git(repo, "branch", "base")
    (repo / name).parent.mkdir(parents=True, exist_ok=True)
    (repo / name).write_text(text)
    _git(repo, "add", "-f", "--", name)
    _git(repo, *commit, "-m", "add")
    _git(repo, "rm", "-q", "--", name)
    _git(repo, *commit, "-m", "delete")
    tracked = subprocess.run(
        ["bash", str(SCRIPT), "--tracked"], cwd=repo, capture_output=True, check=False
    )
    assert (tracked.returncode, _range(repo, "base..HEAD")) == (0, 1)
    assert (_range(repo, "HEAD~1..HEAD"), _range(repo, "base")) == (0, 0)


@pytest.mark.parametrize(
    "name",
    [
        "README.md",
        "LICENSE",
        "NOTICE",
        "CITATION.cff",
        "src/genokeeper/vcf.py",
        "docs/bed.md",
        "ci/runner/Dockerfile",
        ".github/local-ci.json",
        "src/genokeeper/resources/repeats/loci.json",
        "tests/fixtures/hg002/hg002.eh3e.json",
        "docs/demo/dashboard.html",
    ],
)
def test_an_ordinary_file_passes(tmp_path: Path, name: str) -> None:
    assert _guard(tmp_path, name) == {"--staged": 0, "--tracked": 0}


@pytest.mark.parametrize(
    "name",
    [
        "dashboard.html",  # `genokeeper dashboard --out dashboard.html`
        "findings.json",  # `genokeeper query --json > findings.json`
        "reports.bak/explanations/x.json",  # `cp -R reports reports.bak`
        "docs/rev.ATXN1.svg",  # a REViewer image
        "calls.txt",
        "owner.vcf.bak",
        "owner.vcf.xz",
        "owner.ped",
        "owner.csv",
        "lab.zip",
        "reads.bam.part",
        "findings",
    ],
)
def test_a_data_file_outside_the_data_folders_is_refused(
    tmp_path: Path, name: str
) -> None:
    """Privacy audit I2: a derived output carries no kit ID, so where it may sit is the rule."""
    assert _guard(tmp_path, name) == {"--staged": 1, "--tracked": 1}


@pytest.mark.parametrize(
    ("samples", "code"),
    [("sample", 1), ("HG002\towner", 1), ("HG002", 0), ("HG002\tNA12878", 0)],
)
def test_a_vcf_fixture_names_only_public_samples(
    tmp_path: Path, samples: str, code: int
) -> None:
    """Privacy audit I3: the pipeline names every sample `sample`, the user's too, so a
    fixture from a public sample's run says which one."""
    text = f"##fileformat=VCFv4.2\n{VCF_HEADER}\t{samples}\nchr1\t1\t.\tA\tG\n"
    assert _guard(tmp_path, "tests/fixtures/x/eh.vcf.txt", text) == {
        "--staged": code,
        "--tracked": code,
    }


def test_a_sites_only_vcf_passes(tmp_path: Path) -> None:
    text = VCF_HEADER.removesuffix("\tFORMAT") + "\nchr1\t1\t.\tA\tG\n"
    assert _guard(tmp_path, "tests/fixtures/sites.vcf.txt", text) == {
        "--staged": 0,
        "--tracked": 0,
    }


def test_a_file_over_5_mb_is_refused(tmp_path: Path) -> None:
    """The size cap, the backstop for a genotype file under a name no rule lists."""
    assert _guard(tmp_path, "tests/fixtures/notes.txt", "x" * 5_242_881) == {
        "--staged": 1,
        "--tracked": 1,
    }
    (tmp_path / "at").mkdir()
    assert _guard(tmp_path / "at", "tests/fixtures/notes.txt", "x" * 5_242_880) == {
        "--staged": 0,
        "--tracked": 0,
    }


def test_a_file_whose_size_cannot_be_read_is_refused(tmp_path: Path) -> None:
    """A staged gitlink has no blob to measure, a tracked file gone from the disk no bytes."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _init(repo)
    _git(
        repo,
        "update-index",
        "--add",
        "--cacheinfo",
        f"160000,{'1' * 40},tests/fixtures/sub",
    )
    gone = repo / "tests/fixtures/gone.txt"
    gone.parent.mkdir(parents=True)
    gone.write_text("x\n")
    _git(repo, "add", "tests/fixtures/gone.txt")
    gone.unlink()
    runs = {
        mode: subprocess.run(
            ["bash", str(SCRIPT), mode],
            cwd=repo,
            capture_output=True,
            text=True,
            check=False,
        )
        for mode in ("--staged", "--tracked")
    }
    assert "cannot read the size of: tests/fixtures/sub" in runs["--staged"].stderr
    assert (
        "cannot read the size of: tests/fixtures/gone.txt" in runs["--tracked"].stderr
    )
    assert {mode: run.returncode for mode, run in runs.items()} == {
        "--staged": 1,
        "--tracked": 1,
    }


SETTINGS = Path(__file__).resolve().parents[1] / ".claude" / "settings.json"
# The raw delivery and the run outputs (the store, `reports/work/`, the snapshots) relative to
# a session at the repository root, the same anywhere on disk (a worktree's checkout holds
# neither: the files sit in the main checkout), and the calibration outputs.
GENOME_READS = {
    "Read(data/raw/**)",
    "Read(reports/**)",
    "Read(//**/data/raw/**)",
    "Read(//**/reports/**)",
    "Read(~/genome-lab/calibration/**)",
}


def _denies_genome_reads(settings: dict[str, Any]) -> bool:
    return GENOME_READS <= set(settings.get("permissions", {}).get("deny", []))


def test_settings_deny_genome_paths() -> None:
    """The project settings deny the agent's reads of genome data.

    A `Read` deny rule also covers Grep, Glob and LSP; a `Grep(path)` or `Glob(path)` rule
    is never matched (the harness warns and ignores it), so the rules are `Read` rules.
    """
    assert _denies_genome_reads(json.loads(SETTINGS.read_text()))


def test_other_settings_keep_the_deny_check() -> None:
    """A hook, an allow rule or another deny rule is the developer's; a dropped rule fails."""
    deny = [*GENOME_READS, "Read(.env)"]
    settings = {"permissions": {"deny": deny, "allow": ["Bash(ls)"]}, "hooks": {}}
    assert _denies_genome_reads(settings)
    deny.remove("Read(reports/**)")
    assert not _denies_genome_reads(settings)


PRE_PUSH = SCRIPT.parents[1] / ".githooks" / "pre-push"


@pytest.mark.parametrize(
    ("name", "base", "code"),
    [
        ("data/raw/x.txt", "new", 1),  # a new branch: every commit no remote ref holds
        ("data/raw/x.txt", "old", 1),  # an update: the commits after the remote's tip
        ("src/ok.py", "new", 0),
        ("src/ok.py", "old", 0),
    ],
)
def test_pre_push_checks_the_outgoing_commits(
    tmp_path: Path, name: str, base: str, code: int
) -> None:
    """Privacy audit I4: a commit no pre-commit hook saw (a rebase, `--no-verify`) is
    checked before git sends it."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    shutil.copy(SCRIPT, repo / "scripts" / SCRIPT.name)
    _init(repo)
    commit = ("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q")
    _git(repo, *commit, "--allow-empty", "-m", "base")
    old = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    (repo / name).parent.mkdir(parents=True, exist_ok=True)
    (repo / name).write_text("x\n")
    _git(repo, "add", "-f", "--", name)
    _git(repo, *commit, "--no-verify", "-m", "add")
    tip = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    remote = "0" * 40 if base == "new" else old
    pushed = subprocess.run(
        ["sh", str(PRE_PUSH), "origin", "git@example:x.git"],
        cwd=repo,
        input=f"refs/heads/b {tip} refs/heads/b {remote}\n",
        capture_output=True,
        text=True,
        check=False,
    )
    assert pushed.returncode == code, pushed.stderr
    deleted = subprocess.run(
        ["sh", str(PRE_PUSH), "origin", "git@example:x.git"],
        cwd=repo,
        input=f"(delete) {'0' * 40} refs/heads/b {tip}\n",
        capture_output=True,
        text=True,
        check=False,
    )
    assert deleted.returncode == 0
