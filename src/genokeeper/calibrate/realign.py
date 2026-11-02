"""The calibration samples re-aligned from their CRAMs (HG003, both HG002 libraries): their
FASTQ were deleted after the first alignment, so the reads come back out of the old CRAM.

`cram_to_fastq` runs `samtools collate -O -u | samtools fastq -n -F 0x900 -0 /dev/null -s
/dev/null` against NCBI's unmasked FASTA, the only one an old CRAM decodes against: every
primary read, unmapped included and duplicates kept, as a gzipped pair. The pair is refused
unless its reads and bases (`fastq_counts`) equal the old CRAM's (`cram_counts`, `samtools
stats`), and `python -m genokeeper.calibrate.hg002 realign` refuses the new alignment unless
its `total_bases` equals the old CRAM's bases by the same count, never a record's stored
metrics (HG003's first record has none of its own).

`old_record` finds the CRAM a re-alignment may start from in its record, never in a typed
path; its anchor comes from the record's `cram` equalling `<work>/align/<key>/sample.cram`,
the one place the run's align step (`_aligned`, after its full platform gate) writes, which
a record naming another file cannot meet; so the read-name test the gate skips for the
re-aligned pair (`check_platform(names=False)`) was passed by these reads once. This module
imports `records` and `pipeline`, never the `hg002` CLI, which imports it.
"""

import json
import os
import secrets
import shlex
import shutil
import subprocess
from pathlib import Path

from genokeeper.calibrate import calibration_dir, hg003_dir, layout
from genokeeper.calibrate.records import (
    ALIGNED,
    Record,
    label,
    library_dir,
    need_space,
    on_this_reference,
    read_record,
    write_record,
)
from genokeeper.checksum import sha256_cached
from genokeeper.config import unmasked_fasta
from genokeeper.pipeline import GENOME_SAMPLE, IDENTITIES, RunRefused
from genokeeper.tools import ToolFailed, gated, require, run_tool

# Free disk a re-alignment needs, times the old CRAM's size: collate's level-1 temporaries
# (about 1.6x) and the gzipped pair (about 2x) exist together, and the alignment after it
# holds the pair and the new CRAM (about 3x); `need_space` keeps its margin beyond.
REALIGN_SPACE = 3.6
_CRAM = f"{GENOME_SAMPLE}.cram"
_FASTQ = ("*.fq.gz", "*.fastq.gz")


def work_name(key: str) -> str:
    """The folder under the calibration folder that `align_sample` and `run_sample` write
    `key`'s steps in: HG002's `hg002-<key>`, HG003's `hg003`, HG004's `hg004-pcrfree`."""
    if key == "hg003":
        return "hg003"  # hg003_dir's
    return "hg004-pcrfree" if key == "hg004" else f"hg002-{key}"


def cram_counts(cram: Path, fasta: Path, *, threads: int) -> tuple[int, int]:
    """The primary reads and their bases in `cram` (`samtools stats` "raw total sequences"
    and "total length"), decoded against `fasta`; `RunRefused` when `fasta` is missing."""
    if not fasta.is_file():
        raise RunRefused(
            f"missing input: the FASTA {cram} decodes against, {fasta}; run"
            " scripts/bootstrap.sh"
        )
    out = run_tool(
        gated(
            "dna:qc",
            ["samtools", "stats", "-@", str(threads), "--reference", str(fasta)]
            + [str(cram)],
        )
    )
    found: dict[str, int] = {}
    for line in out.splitlines():
        fields = line.split("\t")
        if fields[0] == "SN" and fields[1] in ("raw total sequences:", "total length:"):
            found[fields[1]] = int(fields[2])
    return found["raw total sequences:"], found["total length:"]


def fastq_counts(pair: tuple[Path, Path]) -> tuple[int, int]:
    """The reads and bases of both files of `pair`, as `seqtk size` counts them, summed."""
    reads = bases = 0
    for fastq in pair:
        n, b = run_tool(gated("dna:downsample", ["seqtk", "size", str(fastq)])).split()
        reads, bases = reads + int(n), bases + int(b)
    return reads, bases


def check_counts(expected: tuple[int, int], found: tuple[int, int], what: str) -> None:
    """`RunRefused` unless the FASTQ's reads and bases (`found`) are the CRAM's
    (`expected`)."""
    if found != expected:
        (r, b), (r0, b0) = found, expected
        raise RunRefused(
            f"{what}: the FASTQ holds {r} reads and {b} bases, the CRAM {r0} and {b0}"
        )


def cram_to_fastq(
    cram: Path, fasta: Path, out: tuple[Path, Path], *, threads: int
) -> tuple[Path, Path]:
    """Every primary read of `cram` (decoded against `fasta`), unmapped included and
    duplicates kept, as the gzipped pair `out`, mates in one order and named as in the CRAM:
    `samtools collate -O -u | samtools fastq -n -F 0x900 -0 /dev/null -s /dev/null` as one
    `bash -c` behind the core gate (`dna:downsample`). collate's temporaries go to a hidden
    `.collate-<hex>/` folder beside `out` on the lab's disk, removed after; the pair lands
    under hidden names renamed only when every stage exited 0, both or neither. Raises
    `RunRefused` without `REALIGN_SPACE` times the CRAM's size free, `ToolFailed` when a
    stage fails, `ToolMissing` without samtools; a killed run's leftovers are removed."""
    require("samtools")
    folder = out[0].parent
    folder.mkdir(parents=True, exist_ok=True)
    for pattern in (".collate-*", *(f".{o.name}.*" for o in out)):
        for old in folder.glob(pattern):
            shutil.rmtree(old) if old.is_dir() else old.unlink()
    need_space(
        folder, int(REALIGN_SPACE * cram.stat().st_size), f"re-aligning {cram.name}"
    )
    scratch = folder / f".collate-{secrets.token_hex(4)}"
    scratch.mkdir()
    tag = secrets.token_hex(4)
    # samtools fastq compresses a file whose name ends in .gz.
    tmps = [o.with_name(f".{o.name}.{tag}.gz") for o in out]
    t = str(threads)
    collate = ["samtools", "collate", "-O", "-u", "-@", t, "--reference", str(fasta)]
    collate += ["-T", f"{scratch}/", str(cram)]
    fastq = ["samtools", "fastq", "-n", "-F", "0x900", "-0", "/dev/null", "-s"]
    fastq += ["/dev/null", "-@", t, "-1", str(tmps[0]), "-2", str(tmps[1]), "-"]
    script = f"set -o pipefail; {shlex.join(collate)} | {shlex.join(fastq)}"
    command = gated("dna:downsample", ["bash", "-c", script])
    try:
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
        for tmp, final in zip(tmps, out, strict=True):
            os.replace(tmp, final)
    finally:
        for tmp in tmps:
            tmp.unlink(missing_ok=True)
        shutil.rmtree(scratch, ignore_errors=True)
    return out


def hg003_pair(lab: Path) -> tuple[Path, Path]:
    """The two FASTQ names HG003's identities (`IDENTITIES`, which its first alignment
    wrote) hold in `layout(lab)["downsampled"]`, in name order, there: the names its
    re-aligned pair takes. `RunRefused` naming what was found for any other count."""
    folder = layout(lab)["downsampled"]
    identities = hg003_dir(lab) / IDENTITIES
    known = json.loads(identities.read_text()) if identities.is_file() else {}
    names = sorted(
        Path(p).name
        for p in known
        if Path(p).parent == folder.resolve() and any(Path(p).match(g) for g in _FASTQ)
    )
    if len(names) != 2:
        raise RunRefused(
            f"HG003's identities {identities} name {len(names)} FASTQ in {folder}"
            f" ({', '.join(names) or 'none'}), not the pair its first alignment read"
        )
    return folder / names[0], folder / names[1]


def old_record(lab: Path, key: str) -> Record:
    """The record a re-alignment of `key` (`pcrfree`, `pcr` or `hg003`) starts from: its
    `aligned.json` (3a's legacy shape included) whose `cram` is
    `<work>/align/<record key>/sample.cram`, `<work>` the folder `align_sample` writes in
    (`work_name`). A record on this reference is returned as it is. HG003 with no record
    (2b never wrote one) starts from its single `<hg003>/align/<key>/sample.cram`, whose
    record (`metrics` its bases from `cram_counts`, on the unmasked FASTA) is written
    first. `RunRefused` for anything else, naming what is missing or what was found."""
    work = calibration_dir(lab) / work_name(key)
    record = read_record(lab, key)
    if record is not None:
        if on_this_reference(lab, record):
            return record
        expected = work / "align" / record.key / _CRAM
        if Path(record.cram) != expected:
            raise RunRefused(
                f"{label(key)}'s record {library_dir(lab, key) / ALIGNED} names the CRAM"
                f" {record.cram}, not {expected}, where `align` writes it: nothing to"
                " re-align from"
            )
        if not expected.is_file():
            raise RunRefused(
                f"{label(key)}'s CRAM {expected} is gone: nothing to re-align"
            )
        return record
    if key != "hg003":
        raise RunRefused(
            f"{label(key)} has no record at {library_dir(lab, key) / ALIGNED}: nothing to"
            f" re-align from; run `python -m genokeeper.calibrate.hg002 fetch {key}`, then"
            f" `align {key}`"
        )
    # A killed align step leaves its staging folder hidden; only a renamed one is complete.
    crams = sorted(
        p
        for p in (work / "align").glob(f"*/{_CRAM}")
        if not p.parent.name.startswith(".")
    )
    if len(crams) != 1:
        found = (
            f"{len(crams)} CRAMs ({', '.join(map(str, crams))})" if crams else "no CRAM"
        )
        raise RunRefused(
            f"HG003 has no record and {found} under {work / 'align'}: one is needed to"
            " re-align from"
        )
    (cram,) = crams
    fasta = unmasked_fasta(lab)
    _, bases = cram_counts(cram, fasta, threads=os.cpu_count() or 4)
    record = Record(
        cram=str(cram),
        key=cram.parent.name,
        metrics={"total_bases": bases},
        gvcf=None,
        reference_sha256=sha256_cached(fasta, calibration_dir(lab) / IDENTITIES),
        accepted=True,
        snp_f1=None,
        replaced=None,
    )
    write_record(lab, key, record)
    return record
