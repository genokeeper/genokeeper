"""`python -m genokeeper.demo build`: two runs of GIAB's HG002 into the demo root, the first
on ClinVar `BASE_CLINVAR` and the second on the lab's newest, and the refresh-log record that
labels the diff between them a refresh.

`run_genome` branches on no sample's name: its identity check reads "not checked" without a
delivered file (`delivered=None`), the sex is inferred from the reads, and the GIAB files it
opens are the references every run uses. HG002 therefore runs through the same path as the
owner's genome, but for the platform gate's read-name test (see `build`).
"""

import os
import re
import shutil
import stat
import subprocess
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from genokeeper import pipeline, refresh_log
from genokeeper.calibrate import hg002
from genokeeper.calibrate.records import need_space
from genokeeper.demo import STORE_FILE, demo_paths
from genokeeper.refresh import FOLDERS, REFRESHED, _dated
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus

BASE_CLINVAR = "20260503"
# NCBI files each archived release under its year: `{CLINVAR_BASE_URL}{year}/`.
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/"
# From an empty root; a full build left 101 GB there (2026-11-13). With `MARGIN_BYTES`, 125 GB.
NEED_BYTES = 105 * 10**9


class BuildRefused(Exception):
    """The build cannot start: an input is missing, the disk is short or a download is wrong."""


def base_clinvar_path(root: Path) -> Path:
    """Where the demo keeps its base ClinVar: under its own root, never the lab's `clinvar/`."""
    return root / "clinvar" / BASE_CLINVAR / "clinvar_GRCh38.vcf.gz"


def fetch_base_clinvar(root: Path) -> Path:
    """ClinVar `BASE_CLINVAR` from NCBI's archive into `root` (`hg002.download`), checked against
    the `.md5` NCBI publishes beside it; its `.tbi` is fetched again with every new VCF."""
    out = base_clinvar_path(root)
    url = f"{CLINVAR_BASE_URL}{BASE_CLINVAR[:4]}/clinvar_{BASE_CLINVAR}.vcf.gz"
    try:
        with urllib.request.urlopen(url + ".md5", timeout=60) as response:
            words = response.read().decode("ascii", "replace").split()
    except OSError as e:
        raise BuildRefused(f"cannot read {url}.md5: {e}") from e
    expected = words[0].lower() if words else ""
    if not re.fullmatch(r"[0-9a-f]{32}", expected):
        raise BuildRefused(f"{url}.md5 holds no md5")
    fetched = not out.exists()
    tbi = out.with_name(out.name + ".tbi")
    part = tbi.with_name(tbi.name + ".part")
    try:
        hg002.download(url, expected, out)
        if fetched or not tbi.is_file():
            with (
                urllib.request.urlopen(url + ".tbi", timeout=60) as response,
                part.open("wb") as handle,
            ):
                shutil.copyfileobj(response, handle)
            part.replace(tbi)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as e:
        raise BuildRefused(f"cannot fetch ClinVar {BASE_CLINVAR}: {e}") from e
    return out


def _sources(lab: Path, clinvar: tuple[str, str]) -> tuple[SourceOutcome, ...]:
    """ClinVar moved from the base release; every other source is the lab's newest, both runs."""
    outcomes: list[SourceOutcome] = []
    for name in REFRESHED:
        if name == "clinvar":
            outcomes.append(SourceOutcome(name, *clinvar, SourceStatus.NEW))
        else:
            dated = _dated(lab / FOLDERS[name])
            release = dated[-1].name if dated else None
            outcomes.append(
                SourceOutcome(name, release, release, SourceStatus.UNCHANGED)
            )
    return tuple(outcomes)


def _held(root: Path) -> int:
    """The allocated bytes of the regular files under `root`, each inode once, no symlink
    followed; 0 when `root` does not exist."""
    sizes: dict[tuple[int, int], int] = {}
    for folder, _, names in os.walk(root):
        for name in names:
            st = os.lstat(os.path.join(folder, name))
            if stat.S_ISREG(st.st_mode):
                sizes[st.st_dev, st.st_ino] = st.st_blocks * 512
    return sum(sizes.values())


def build(lab: Path, *, threads: int | None) -> tuple[str, str]:
    """Run HG002's PCR-free pair twice into `demo_paths(lab)`, then log the second as a
    refresh of the first, each ClinVar release as its run records it. Returns the two run
    ids, oldest first. Refuses (`BuildRefused`) before any download or run when the demo
    store exists, the pair is missing, the disk is short of what the root does not hold yet
    or the lab's newest ClinVar is not newer than `BASE_CLINVAR`."""
    paths = demo_paths(lab)
    root = paths.repo
    if (paths.reports / STORE_FILE).exists():
        raise BuildRefused(
            f"{paths.reports / STORE_FILE} exists: to build again, delete it, its"
            f" snapshots ({paths.reports / 'snapshots'}) and {root / 'lab' / 'refresh'};"
            " the work caches stay and are reused"
        )
    pair = hg002.downsampled_pair(lab, "pcrfree")
    missing = [p for p in pair if not p.is_file()]
    if missing:
        raise BuildRefused(
            f"HG002's PCR-free pair is missing ({', '.join(p.name for p in missing)}"
            f" under {pair[0].parent}); fetch it with `genokeeper calibrate`"
        )
    nearest = next(p for p in (root, *root.parents) if p.exists())
    try:
        need_space(nearest, max(NEED_BYTES - _held(root), 0), "the HG002 demo build")
        faf = pipeline.open_faf_table(lab)
        newest = pipeline.newest_clinvar38(lab)
        after = pipeline._clinvar_release(newest)
        if after <= BASE_CLINVAR:
            raise BuildRefused(
                f"the lab's newest ClinVar ({after}) is not newer than {BASE_CLINVAR}:"
                " run `genokeeper refresh` first"
            )
        root.mkdir(parents=True, exist_ok=True)
        started = datetime.now(UTC)
        base = fetch_base_clinvar(root)
        before = pipeline._clinvar_release(base)
    except pipeline.RunRefused as e:
        raise BuildRefused(str(e)) from e
    first, second = (
        pipeline.run_genome(
            paths,
            fastq1=pair[0],
            fastq2=pair[1],
            clinvar38=clinvar,
            faf=faf,
            threads=threads,
            delivered=None,
            # Written back from a calibration CRAM (`calibrate.realign`): the names carry the
            # accession alone, and these reads passed the full gate when it was aligned.
            skip_read_names=True,
        )
        for clinvar in (base, newest)
    )
    refresh_log.append(
        root / "lab",
        RefreshRecord(
            started=started,
            ended=datetime.now(UTC),
            outcome="complete",
            sources=_sources(lab, (before, after)),
            base_run=first,
            pre_run=None,
            run=second,
        ),
    )
    return first, second
