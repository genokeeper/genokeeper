"""Stage 2a alignment: the platform gate on the delivered FASTQ, and bwa mem to a
coordinate-sorted, duplicate-marked CRAM on the GRCh38 no-alt FASTA.

The platform gate stops reads the T7 checkpoint was not trained on: a DNBSEQ-G400 run
(flowcell `V...`) or a DNBSEQ-T7+ run (flowcell `ML...`, base qualities binned to four
levels). Its anchor comes from the FASTQ itself, the flowcell ID the sequencer writes into
every read name and the quality characters of every base, which a mislabelled delivery or a
file mix-up cannot move because no one downstream of the instrument rewrites them. A
sequence archive (SRA, ENA, DDBJ) puts its run accession first and keeps the instrument's
name as the second field; the gate reads that one.
"""

import gzip
import os
import re
import secrets
import shlex
import shutil
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from genokeeper.tools import ToolFailed, gated, require

# The original DNBSEQ-T7 names its flowcells E...: the platform the T7 checkpoint wants.
T7_FLOWCELL = "E"
# T7+ bins base qualities to four levels; an unbinned T7 run writes dozens.
MIN_QUALITY_LEVELS = 5
_FLOWCELL = re.compile(r"^@([A-Z]+)\d")
_ARCHIVED = re.compile(r"^@[SED]RR\d+\.\d+\s+(?=[A-Z]+\d)")
_KNOWN = {"E": "DNBSEQ-T7", "V": "DNBSEQ-G400", "ML": "DNBSEQ-T7+"}


@dataclass(frozen=True)
class Platform:
    """What the FASTQ's first reads show: the flowcell's letter prefix and how many distinct
    base quality values occur."""

    flowcell_prefix: str
    quality_levels: int


class PlatformRefused(RuntimeError):
    """The FASTQ is not from the platform the T7 checkpoint was trained on."""


def _records(fastq: Path) -> Iterator[tuple[str, str]]:
    """(name line, quality line) of each record, gzipped or plain."""
    with fastq.open("rb") as raw:
        gzipped = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    with opener(fastq, "rt") as f:
        while name := f.readline():
            f.readline()
            f.readline()
            yield name.rstrip("\n"), f.readline().rstrip("\n")


def check_platform(
    fastq1: Path, *, reads: int = 100_000, names: bool = True
) -> Platform:
    """The platform of the first `reads` records of `fastq1`; raises `PlatformRefused` unless
    every read name carries an E... flowcell and the qualities take at least
    `MIN_QUALITY_LEVELS` values. Only aggregates leave this function.

    `names=False` skips the flowcell tests and reports the first alignment's `E`: a FASTQ
    written back from a calibration CRAM carries ENA's accession alone (bwa kept the
    name's first word), and its reads passed the full gate when that CRAM was aligned. The
    quality test runs as before."""
    prefixes: set[str] = set()
    qualities: set[str] = set()
    seen = 0
    for name, quality in _records(fastq1):
        match = _FLOWCELL.match(_ARCHIVED.sub("@", name, count=1))
        prefixes.add(match.group(1) if match else "?")
        qualities.update(quality)
        seen += 1
        if seen == reads:
            break
    if not seen:
        raise PlatformRefused(f"{fastq1}: no reads")
    if not names:
        prefixes = {T7_FLOWCELL}
    if len(prefixes) > 1:
        raise PlatformRefused(
            f"{fastq1}: reads from flowcells {' and '.join(sorted(prefixes))}; one run's"
            " FASTQ comes from one kind of flowcell"
        )
    prefix = prefixes.pop()
    if prefix != T7_FLOWCELL:
        kind = _KNOWN.get(prefix, "an unknown platform")
        raise PlatformRefused(
            f"{fastq1}: flowcell prefix {prefix} ({kind}); the T7 checkpoint is for the"
            f" original DNBSEQ-T7 ({T7_FLOWCELL}...)"
        )
    if len(qualities) < MIN_QUALITY_LEVELS:
        raise PlatformRefused(
            f"{fastq1}: {len(qualities)} base quality levels (binned, as DNBSEQ-T7+ writes"
            f" them); the T7 checkpoint was trained on unbinned qualities"
        )
    return Platform(flowcell_prefix=prefix, quality_levels=len(qualities))


def align(
    fastq1: Path, fastq2: Path, fasta: Path, out: Path, *, sample: str, threads: int
) -> Path:
    """`bwa mem | samtools fixmate -m | samtools sort | samtools markdup -O cram` into `out`,
    indexed, as one unit behind the core gate (`dna:align`). bwa finds its index beside
    `fasta` (`<fasta>.bwt` and the rest, as the bootstrap unpacks it).

    Every stage's exit status is checked: the first that failed raises `ToolFailed` with its
    stderr. The CRAM and its index are written under temp names and renamed only when every
    stage succeeded, so a failed or killed run leaves no `out`."""
    require("bwa", "samtools")
    work = out.parent / f".align-{secrets.token_hex(4)}"
    work.mkdir()
    tmp = out.with_name(f".{out.name}.{secrets.token_hex(4)}.cram")
    t = str(threads)
    group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:DNBSEQ"
    stages = [
        # -K: a fixed batch size makes the output independent of the thread count.
        ["bwa", "mem", "-t", t, "-K", "100000000", "-Y", "-R", group, str(fasta)]
        + [str(fastq1), str(fastq2)],
        ["samtools", "fixmate", "-m", "-u", "-@", t, "-", "-"],
        ["samtools", "sort", "-u", "-@", t, "-T", str(work / "sort"), "-"],
        ["samtools", "markdup", "-@", t, "-T", str(work / "markdup")]
        + ["-O", "cram", "--reference", str(fasta), "-", str(tmp)],
        ["samtools", "index", "-@", t, str(tmp)],
    ]
    status = work / "status"
    pipe = " | ".join(
        f"{shlex.join(stage)} 2> {shlex.quote(str(work / f'err{i}'))}"
        for i, stage in enumerate(stages[:-1])
    )
    index = f"{shlex.join(stages[-1])} 2> {shlex.quote(str(work / 'err4'))}"
    script = (
        f'{pipe}\ns=("${{PIPESTATUS[@]}}")\n'
        f'if [ "${{s[*]}}" = "0 0 0 0" ]; then {index}; s+=($?); fi\n'
        f'echo "${{s[*]}}" > {shlex.quote(str(status))}\n'
    )
    try:
        ran = subprocess.run(
            gated("dna:align", ["bash", "-c", script]),
            check=False,  # each stage's status is read from the status file below
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            errors="replace",
        )
        codes = [int(c) for c in status.read_text().split()] if status.exists() else []
        if not codes:
            raise ToolFailed(ran.returncode, ["bash", "-c", script], None, ran.stderr)
        for i, code in enumerate(codes):
            if code:
                stderr = (work / f"err{i}").read_text(errors="replace")
                raise ToolFailed(code, stages[i], None, stderr)
        os.replace(f"{tmp}.crai", f"{out}.crai")
        os.replace(tmp, out)
    finally:
        tmp.unlink(missing_ok=True)
        Path(f"{tmp}.crai").unlink(missing_ok=True)
        shutil.rmtree(work, ignore_errors=True)
    return out
