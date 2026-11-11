"""`python -m genokeeper.demo build`: two runs of GIAB's HG002 into the demo root, the first
on ClinVar 20261011 and the second on the lab's newest, and the refresh-log record that
labels the diff between them a refresh.

`run_genome` branches on no sample's name: its identity check reads "not checked" without a
delivered file (`delivered=None`), the sex is inferred from the reads, and the GIAB files it
opens are the references every run uses. HG002 therefore runs through the same path as the
owner's genome.
"""

import hashlib
import os
import re
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from genokeeper import pipeline, refresh_log
from genokeeper.calibrate import hg002
from genokeeper.calibrate.records import need_space
from genokeeper.demo import demo_paths
from genokeeper.refresh import FOLDERS, REFRESHED
from genokeeper.refresh_log import RefreshRecord, SourceOutcome, SourceStatus

BASE_CLINVAR = "20261011"
CLINVAR_BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/"
)
NEED_BYTES = 130 * 10**9  # with `MARGIN_BYTES`, the spec's 150 GB
_CHUNK = 1 << 20
_DATE = re.compile(r"\d{8}")


class BuildRefused(Exception):
    """The build cannot start: an input is missing, the disk is short or a download is wrong."""


def base_clinvar_path(root: Path) -> Path:
    """Where the demo keeps ClinVar 20261011: under its own root, never the lab's `clinvar/`."""
    return root / "clinvar" / BASE_CLINVAR / "clinvar_GRCh38.vcf.gz"


def _fetch(url: str, out: Path) -> str:
    """Write `url` to `<out>.part` and return its md5; the caller renames it."""
    part = out.with_name(out.name + ".part")
    digest = hashlib.md5(usedforsecurity=False)
    with urllib.request.urlopen(url, timeout=60) as response, part.open("wb") as handle:
        while chunk := response.read(_CHUNK):
            digest.update(chunk)
            handle.write(chunk)
    return digest.hexdigest()


def _md5_of(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_base_clinvar(root: Path) -> Path:
    """ClinVar 20261011 from NCBI's archive into `root`, checked against the `.md5` NCBI
    publishes beside it; a file already there with that md5 is kept. Nothing appears under
    the final name unless the md5 matches."""
    out = base_clinvar_path(root)
    url = f"{CLINVAR_BASE_URL}clinvar_{BASE_CLINVAR}.vcf.gz"
    try:
        with urllib.request.urlopen(url + ".md5", timeout=60) as response:
            words = response.read().decode("ascii", "replace").split()
    except OSError as e:
        raise BuildRefused(f"cannot read {url}.md5: {e}") from e
    expected = words[0].lower() if words else ""
    if not re.fullmatch(r"[0-9a-f]{32}", expected):
        raise BuildRefused(f"{url}.md5 holds no md5")
    out.parent.mkdir(parents=True, exist_ok=True)
    if not (out.is_file() and _md5_of(out) == expected):
        out.unlink(missing_ok=True)
        part = out.with_name(out.name + ".part")
        try:
            got = _fetch(url, out)
        except OSError as e:
            raise BuildRefused(f"cannot download {url}: {e}") from e
        if got != expected:
            part.replace(out.with_name(f"{out.name}.mismatch-{got}"))
            raise BuildRefused(f"{url} has md5 {got}, NCBI's .md5 says {expected}")
        part.replace(out)
    tbi = out.with_name(out.name + ".tbi")
    if not tbi.is_file():
        try:
            _fetch(url + ".tbi", tbi)
        except OSError as e:
            raise BuildRefused(f"cannot download {url}.tbi: {e}") from e
        tbi.with_name(tbi.name + ".part").replace(tbi)
    return out


def _newest_release(folder: Path) -> str | None:
    dated = sorted(d.name for d in folder.glob("*") if _DATE.fullmatch(d.name))
    return dated[-1] if dated else None


def _sources(lab: Path, clinvar_after: str) -> tuple[SourceOutcome, ...]:
    """ClinVar moved from the base release; every other source is the lab's newest, both runs."""
    outcomes: list[SourceOutcome] = []
    for name in REFRESHED:
        if name == "clinvar":
            outcomes.append(
                SourceOutcome(name, BASE_CLINVAR, clinvar_after, SourceStatus.NEW)
            )
        else:
            release = _newest_release(lab / FOLDERS[name])
            outcomes.append(
                SourceOutcome(name, release, release, SourceStatus.UNCHANGED)
            )
    return tuple(outcomes)


def build(lab: Path, *, threads: int | None) -> tuple[str, str]:
    """Run HG002's PCR-free pair twice into `demo_paths(lab)`, then log the second as a
    refresh of the first. Returns the two run ids, oldest first. Refuses (`BuildRefused`)
    before any download or run when the pair is missing or the disk is short."""
    paths = demo_paths(lab)
    pair = hg002.downsampled_pair(lab, "pcrfree")
    missing = [p for p in pair if not p.is_file()]
    if missing:
        raise BuildRefused(
            f"HG002's PCR-free pair is missing ({', '.join(p.name for p in missing)}"
            f" under {pair[0].parent}); fetch it with `genokeeper calibrate`"
        )
    root = paths.repo
    nearest = next(p for p in (root, *root.parents) if p.exists())
    try:
        need_space(nearest, NEED_BYTES, "the HG002 demo build")
        faf = pipeline.open_faf_table(lab)
        newest = pipeline.newest_clinvar38(lab)
    except pipeline.RunRefused as e:
        raise BuildRefused(str(e)) from e
    root.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)
    base = fetch_base_clinvar(root)
    previous = os.umask(0o077)  # as `genokeeper run`: the store is its maker's alone
    try:
        first, second = (
            pipeline.run_genome(
                paths,
                fastq1=pair[0],
                fastq2=pair[1],
                clinvar38=clinvar,
                faf=faf,
                threads=threads,
                delivered=None,
            )
            for clinvar in (base, newest)
        )
    finally:
        os.umask(previous)
    refresh_log.append(
        root / "lab",
        RefreshRecord(
            started=started,
            ended=datetime.now(UTC),
            outcome="complete",
            sources=_sources(lab, newest.parent.name),
            base_run=first,
            pre_run=None,
            run=second,
        ),
    )
    return first, second
