"""HG002's two public DNBSEQ-T7 libraries (ENA, PRJNA1465120: `pcrfree` and `pcr`) and HG004's
PCR-free one (`hg004`, CNGB's CNP0008319), downloaded resumably through their md5, downsampled
by seqtk to a fraction of the read pairs, then aligned like the user's genome.

`download` fetches a file whole in fixed-size byte ranges over many connections (ENA serves one
connection slowly) into a hidden part file with a journal of finished chunks, so a killed run
resumes where it stopped; a chunk holding the storage gateway's in-band error is fetched
again, and the md5 of the whole part gates the rename. `downsample_fastq` then
runs `gzip -dc | seqtk sample | bgzip` over it, and the raw file is deleted. `align` keeps the
CRAM's place in `aligned.json`, so the calibration finds it after the downsampled FASTQ is
deleted. The records and their folders are `genokeeper.calibrate.records`', re-exported here.
`realign` re-aligns a library aligned on another reference from its own CRAM
(`genokeeper.calibrate.realign`), and `clean` deletes what a current record replaced.

    python -m genokeeper.calibrate.hg002 fetch [--connections N] {hg004,pcr,pcrfree}
    python -m genokeeper.calibrate.hg002 align {hg004,pcr,pcrfree}
    python -m genokeeper.calibrate.hg002 realign [--threads N] {hg003,pcr,pcrfree}
    python -m genokeeper.calibrate.hg002 clean {hg003,hg004,pcr,pcrfree}
    python -m genokeeper.calibrate.hg002 accept hg004   # HG004's SNP F1 against its truth
    python -m genokeeper.calibrate.hg002 layout    # key, sample, aligned.json, the pair
"""

import argparse
import glob
import hashlib
import mmap
import os
import re
import secrets
import shlex
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path
from typing import IO

from genokeeper.calibrate import hg004_truth, layout, realign, rtg_path
from genokeeper.calibrate.measure import count_cells
from genokeeper.calibrate.records import (
    ALIGNED,
    HG004_MIN_SNP_F1,
    MARGIN_BYTES,
    RECORDED,
    Record,
    Replaced,
    hg002_dir,
    label,
    library_dir,
    need_space,
    on_this_reference,
    read_aligned,
    read_called,
    read_record,
    realign_hint,
    reference_sha256,
    write_record,
)
from genokeeper.calibrate.sample import (
    SEED,
    align_sample,
    downsample,
    run_sample,
    sample_fraction,
)
from genokeeper.config import Paths, load_paths, unmasked_fasta
from genokeeper.genotype.reliability import StrataIndex
from genokeeper.masking import load_mask
from genokeeper.pipeline import RunRefused, calling_inputs, recorded_bases
from genokeeper.regions import read_bed, subtract
from genokeeper.sex import Sex
from genokeeper.tools import ToolFailed, ToolMissing, gated, require, run_tool

__all__ = [
    "ALIGNED",
    "HG004_MIN_SNP_F1",
    "MARGIN_BYTES",
    "RECORDED",
    "Record",
    "hg002_dir",
    "library_dir",
    "need_space",
    "read_aligned",
    "read_called",
    "read_record",
    "write_record",
]

CHUNK_BYTES = 64 * 2**20
BGZIP_THREADS = 4
_CURL = ["curl", "-fsS", "--retry", "5", "--retry-all-errors"]
# A connection that stalls for two minutes is cut and retried.
_STALL = ["--connect-timeout", "30", "--speed-limit", "10000", "--speed-time", "120"]

_ENA = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
# The gateway's in-band error, observed 2026-10-27: it overwrites part of a correct 206 body.
_GATEWAY_ERROR = b"<Error><Code>"
_TRIES = 5
# A chunk's tries when the body comes back short (a redirect saved as the body).
_SHORT_TRIES = 20


@dataclass(frozen=True)
class Library:
    """A public library of a GIAB sample: its key, run (ENA's or CNGB's), files and their
    md5s, its base count (None where the archive publishes none: `fetch` counts the pair),
    the sample and its sex, and whether DeepVariant calls it (`calls`) or it is aligned
    only."""

    name: str
    run: str
    urls: tuple[str, str]
    md5s: tuple[str, str]
    bases: int | None
    sample: str = "HG002"
    sex: Sex = "male"
    calls: bool = False


def _library(
    name: str,
    run: str,
    md5s: tuple[str, str],
    bases: int,
    *,
    sample: str = "HG002",
    sex: Sex = "male",
    calls: bool = False,
) -> Library:
    folder = f"{_ENA}/{run[:6]}/0{run[-2:]}/{run}"
    return Library(
        name,
        run,
        (f"{folder}/{run}_1.fastq.gz", f"{folder}/{run}_2.fastq.gz"),
        md5s,
        bases,
        sample,
        sex,
        calls,
    )


# md5s and base counts from ENA's filereport for PRJNA1465120, read 2026-10-24. 3a's and
# 3b's calibrations iterate HG002's two libraries; DeepVariant calls the PCR-free one only.
LIBRARIES: dict[str, Library] = {
    "pcrfree": _library(
        "pcrfree",
        "SRR38517577",
        ("c018da93350a7b0313b9edb27ea56aa5", "0073c089fad9d9d9f4af2975389af638"),
        294758313900,
        calls=True,
    ),
    "pcr": _library(
        "pcr",
        "SRR38517578",
        ("f5cffb64f16efeca3b06f602b1d0d171", "d4a2a7484c27a4829ee961eb3a0a22eb"),
        151918693321,
    ),
}
# HG004's PCR-free DNBSEQ-T7 run in CNGB's CNP0008319 (sample CNS1414629, NA24143, female),
# from the same lane as 2b's HG003; the md5s are CNGBdb's metadata table's (read 2026-11-01),
# which gives no base count. ENA's PRJNA1465120 is no source: the run under BioSample HG004
# (SRR38517571, library HG003-PCRFREE) aligned as male on 2026-11-01, and its HG004-PCRFREE
# run holds no reads. `accept hg004` checks the calls against HG004's truth.
_CNGB_HG004 = (
    "https://ftp.cngb.org/pub/CNSA/data7/CNP0008319/CNS1414629/CNX1267210/CNR1441269"
)
HG004 = Library(
    "hg004",
    "CNR1441269",
    (
        f"{_CNGB_HG004}/E200029822_L01_UDB-452_1.fq.gz",
        f"{_CNGB_HG004}/E200029822_L01_UDB-452_2.fq.gz",
    ),
    ("179cb9ee9f6baa77ba7e23ea472365f5", "1a37dc8d44c2ec77583eb7508a1b9d88"),
    None,
    sample="HG004",
    sex="female",
    calls=True,
)
SAMPLE_LIBRARIES: dict[str, Library] = {**LIBRARIES, "hg004": HG004}


def library_sex(key: str) -> Sex:
    """The sex of the sample behind `key` (HG003, the trio's father, is male)."""
    return "male" if key == "hg003" else SAMPLE_LIBRARIES[key].sex


def downsampled_pair(lab: Path, library: str) -> tuple[Path, Path]:
    """Where `fetch` puts the library's downsampled FASTQ pair."""
    folder = library_dir(lab, library) / "downsampled"
    urls = SAMPLE_LIBRARIES[library].urls
    first, second = (folder / url.rsplit("/", 1)[1] for url in urls)
    return first, second


def _size(url: str) -> int:
    """The file's size, from its headers (the last Content-Length after redirects)."""
    head = run_tool(["curl", "-fsSIL", *_STALL, "--retry", "5", url])
    sizes = re.findall(r"^content-length:\s*(\d+)", head, re.IGNORECASE | re.MULTILINE)
    if not sizes:
        raise ValueError(f"{url} gave no Content-Length")
    return int(sizes[-1])


def _journaled(journal: Path, part: Path, size: int) -> set[int] | None:
    """The chunk indices the journal records as written into `part`, or None when the
    journal or part is missing or was made for another size or chunk size."""
    if not journal.is_file() or not part.is_file() or part.stat().st_size != size:
        return None
    lines = journal.read_text().split("\n")
    if lines[0] != f"{size} {CHUNK_BYTES}":
        return None
    # The last element is empty, or a line a kill cut short.
    return {int(line) for line in lines[1:-1] if line.isdigit()}


def _gateway_errors(path: Path) -> list[int]:
    """The offsets of `_GATEWAY_ERROR` in the (non-empty) file at `path`."""
    with path.open("rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as m:
        found: list[int] = []
        at = m.find(_GATEWAY_ERROR)
        while at != -1:
            found.append(at)
            at = m.find(_GATEWAY_ERROR, at + 1)
    return found


def _pwrite_all(fd: int, data: bytes, offset: int) -> None:
    """Write all of `data` at `offset`: `os.pwrite` may write less than it is given."""
    view = memoryview(data)
    while view:
        written = os.pwrite(fd, view, offset)
        view = view[written:]
        offset += written


def _check_range(headers: str, first: int, last: int, size: int) -> str | None:
    """Why the last response in `headers` (curl's `-D` dump) is not the bytes
    `first`-`last` of a `size`-byte file, or None when it is."""
    blocks = [b for b in re.split(r"\r?\n\r?\n", headers) if b.startswith("HTTP/")]
    if not blocks:
        return "no response headers"
    head = blocks[-1]
    status = head.split(None, 2)[1]
    if status != "206":
        return f"HTTP {status}, not 206"
    got = re.search(
        r"^content-range:\s*bytes (\d+)-(\d+)/(\d+)", head, re.IGNORECASE | re.MULTILINE
    )
    if got is None or tuple(map(int, got.groups())) != (first, last, size):
        return f"Content-Range {got.group(0) if got else 'missing'}"
    return None


class _Download:
    """Missing chunks of `url` fetched with `curl`, each into its own file, then written into
    the part at its offset and journaled; the transfers in flight are killed when stopped."""

    def __init__(
        self, url: str, size: int, part_fd: int, journal: IO[str], window: Path
    ):
        self.url, self.size, self.fd, self.journal, self.window = (
            url,
            size,
            part_fd,
            journal,
            window,
        )
        self.live: set[subprocess.Popen[str]] = set()
        self.lock = threading.Lock()
        self.stopped = False

    def get(self, i: int) -> None:
        first = i * CHUNK_BYTES
        last = min(first + CHUNK_BYTES, self.size) - 1
        dest = self.window / str(i)
        headers = self.window / f"{i}.headers"
        command = [*_CURL, *_STALL, "-r", f"{first}-{last}", "-D", str(headers)]
        command += ["-o", str(dest), self.url]
        gateway = 0
        for _ in range(_SHORT_TRIES):
            self.curl(command)
            # A short body is fetched again: CNGB's HTTP answers about two in five
            # requests with a redirect to its HTTPS, observed 2026-11-01.
            if dest.stat().st_size != last - first + 1:
                continue
            # file:// has no response headers to check.
            if self.url.startswith(("http://", "https://")):
                bad = _check_range(
                    headers.read_text(errors="replace"), first, last, self.size
                )
                if bad:
                    raise ValueError(f"{self.url} bytes {first}-{last}: {bad}")
            if not _gateway_errors(dest):
                break
            gateway += 1
            if gateway == _TRIES:
                raise ValueError(
                    f"{self.url} bytes {first}-{last} (chunk {i}): the gateway's"
                    f" in-band error in {_TRIES} tries"
                )
        else:
            raise ValueError(
                f"{self.url} bytes {first}-{last}: the server sent another size in"
                f" {_SHORT_TRIES} tries"
            )
        with dest.open("rb") as f:
            offset = first
            while block := f.read(8 * 2**20):
                _pwrite_all(self.fd, block, offset)
                offset += len(block)
        dest.unlink()
        headers.unlink(missing_ok=True)
        os.fsync(self.fd)
        with self.lock:  # the part holds the chunk before the journal says so
            self.journal.write(f"{i}\n")
            self.journal.flush()
            os.fsync(self.journal.fileno())

    def curl(self, command: list[str]) -> None:
        with self.lock:
            if self.stopped:
                raise ToolFailed(-15, command, None, "stopped")
            proc = subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                errors="replace",
            )
            self.live.add(proc)
        _, err = proc.communicate()
        with self.lock:
            self.live.discard(proc)
        if proc.returncode != 0:
            raise ToolFailed(proc.returncode, command, None, err)

    def stop(self) -> None:
        with self.lock:
            self.stopped = True
            for proc in self.live:
                proc.terminate()


class _Progress:
    """One stderr line about every 1% of the chunks: chunks and GB done, MB/s since the last
    line, and the time left at that rate."""

    def __init__(self, name: str, size: int, chunks: int, done: int, done_bytes: int):
        self.name, self.size, self.chunks = name, size, chunks
        self.done, self.bytes = done, done_bytes
        self.step = max(1, chunks // 100)
        self.mark = (time.monotonic(), done_bytes)
        self.line(f"{done} of {chunks} chunks already in the journal", 0.0)

    def line(self, note: str, rate: float) -> None:
        left = self.size - self.bytes
        eta = f", {left / rate / 3600:.1f} h left" if rate > 0 else ""
        print(
            f"hg002 {self.name}: {note}; {self.bytes / 1e9:.1f} of {self.size / 1e9:.1f} GB"
            f", {rate / 1e6:.1f} MB/s{eta}",
            file=sys.stderr,
            flush=True,
        )

    def chunk(self, nbytes: int) -> None:
        self.done += 1
        self.bytes += nbytes
        if self.done % self.step and self.done != self.chunks:
            return
        now = time.monotonic()
        then, was = self.mark
        self.mark = (now, self.bytes)
        rate = (self.bytes - was) / (now - then) if now > then else 0.0
        self.line(f"{self.done}/{self.chunks} chunks", rate)


def download(url: str, md5: str, out: Path, *, connections: int = 24) -> Path:
    """Fetch `url` whole to `out`, `connections` byte ranges at a time, resumably: ranges
    land in a hidden part file beside `out`, and a journal records each finished chunk
    (fsynced), so a rerun after a kill or a failed range fetches only the missing chunks (a
    part or journal made for another size starts over). `out` appears only once the whole
    part has the md5 `md5`. A chunk holding the gateway's in-band error is fetched again, up
    to 5 times, and journaled chunks holding it are taken out of the journal and fetched
    again. On a mismatch the part is kept as `<out>.mismatch-<its md5>` (the next call
    removes it) and the journal removed, and the error gives both md5s and the chunks this
    run fetched and the journal had. An `out` already there was verified and is kept.
    Raises `ValueError` (md5, size, the gateway's error), `ToolFailed` (a transfer),
    `RunRefused` (too little free disk), `ToolMissing` without curl."""
    if out.exists():
        return out
    require("curl")
    out.parent.mkdir(parents=True, exist_ok=True)
    part, journal_path = (out.with_name(f"{out.name}.{x}") for x in ("part", "journal"))
    window = out.with_name(f"{out.name}.chunks")
    # Files kept by an earlier md5 mismatch go, and their disk counts below.
    for kept in out.parent.glob(f"{glob.escape(out.name)}.mismatch-*"):
        kept.unlink()
    size = _size(url)
    chunks = -(-size // CHUNK_BYTES)
    done = _journaled(journal_path, part, size)
    if done is None:
        part.unlink(missing_ok=True)
        journal_path.unlink(missing_ok=True)
        done = set()
    # Earlier code journaled chunks without looking for the gateway's error.
    marked = {at // CHUNK_BYTES for at in _gateway_errors(part)} if done else set()
    if marked & done:
        done -= marked
        tmp = journal_path.with_name(f"{journal_path.name}.tmp")
        tmp.write_text(
            "".join(f"{x}\n" for x in [f"{size} {CHUNK_BYTES}", *sorted(done)])
        )
        os.replace(tmp, journal_path)
    missing = [i for i in range(chunks) if i not in done]

    def length(i: int) -> int:
        return min(CHUNK_BYTES, size - i * CHUNK_BYTES)

    left = sum(length(i) for i in missing)
    need_space(out.parent, left, f"{out.name} ({size // 10**9} GB)")
    shutil.rmtree(window, ignore_errors=True)  # what a kill left
    window.mkdir()
    fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        if not done:
            os.ftruncate(fd, size)
            journal_path.write_text(f"{size} {CHUNK_BYTES}\n")
        with journal_path.open("a") as journal:
            job = _Download(url, size, fd, journal, window)
            progress = _Progress(out.name, size, chunks, len(done), size - left)
            with ThreadPoolExecutor(max_workers=connections) as pool:
                futures = {pool.submit(job.get, i): i for i in missing}
                try:
                    for future in as_completed(futures):
                        future.result()
                        progress.chunk(length(futures[future]))
                except BaseException:
                    job.stop()
                    pool.shutdown(wait=True, cancel_futures=True)
                    raise
        digest = hashlib.md5()
        with part.open("rb") as f:
            while block := f.read(8 * 2**20):
                digest.update(block)
    finally:
        os.close(fd)
        shutil.rmtree(window, ignore_errors=True)
    if digest.hexdigest() != md5:
        # Kept for inspection; the next `download` removes it.
        mismatch = out.with_name(f"{out.name}.mismatch-{digest.hexdigest()}")
        os.replace(part, mismatch)
        journal_path.unlink(missing_ok=True)
        raise ValueError(
            f"md5 mismatch for {url}: {digest.hexdigest()}, not {md5}"
            f" ({len(missing)} chunks downloaded in this run, {len(done)} resumed from"
            f" the journal); the file is kept as {mismatch}"
        )
    os.replace(part, out)
    journal_path.unlink(missing_ok=True)
    return out


def downsample_fastq(
    raw: Path, out: Path, *, fraction: float, seed: int = SEED
) -> Path:
    """`raw` (gzip FASTQ) through `gzip -dc | seqtk sample -s seed - fraction | bgzip` into a
    hidden temporary beside `out`, renamed to `out` only when all three exited 0 (seqtk alone
    exits 0 on a gzip cut short). Raises `ToolFailed`, `ToolMissing` and `RunRefused` (too
    little free disk for the output); hidden leftovers of an earlier run are removed."""
    require("gzip", "seqtk", "bgzip")
    out.parent.mkdir(parents=True, exist_ok=True)
    for old in out.parent.glob(f".{out.name}.*"):
        shutil.rmtree(old) if old.is_dir() else old.unlink()
    need_space(out.parent, int(raw.stat().st_size * fraction), f"{out.name}")
    tmp = out.with_name(f".{out.name}.{secrets.token_hex(4)}")
    sample = shlex.join(["seqtk", "sample", "-s", str(seed), "-", repr(fraction)])
    script = (
        f"set -o pipefail; gzip -dc < {shlex.quote(str(raw))} | {sample}"
        f" | bgzip -c -@ {BGZIP_THREADS} > {shlex.quote(str(tmp))}"
    )
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
        os.replace(tmp, out)
    finally:
        tmp.unlink(missing_ok=True)
    return out


def fetch_one(
    url: str, md5: str, raw: Path, out: Path, *, fraction: float, connections: int = 24
) -> Path:
    """`download` the file to the hidden `raw`, `downsample_fastq` it into `out`, delete
    `raw`. An interrupted run resumes the download, or redoes only the downsample when `raw`
    is complete."""
    download(url, md5, raw, connections=connections)
    downsample_fastq(raw, out, fraction=fraction)
    raw.unlink()
    return out


def fetch(
    lab: Path, library: Library, *, target_bases: int, connections: int = 24
) -> tuple[Path, Path]:
    """The library's downsampled pair (`downsampled_pair`), each file keeping the fraction
    that brings the pair to `target_bases`, one file's raw FASTQ on disk at a time. A file
    already there is not fetched again, and nothing is fetched for a library with a record
    (its pair was deleted after alignment; one on another reference is re-aligned from its
    CRAM, never read from ENA again). A library without a base count has both raw files
    downloaded into `raw/` first, its pair's bases counted by `sample.downsample`, and both
    raw files deleted once the pair has landed."""
    folder = library_dir(lab, library.name)
    outs = downsampled_pair(lab, library.name)
    if read_record(lab, library.name) is not None:
        return outs
    if library.bases is None:
        if not all(out.exists() for out in outs):
            raws = [folder / "raw" / out.name for out in outs]
            for url, md5, raw in zip(library.urls, library.md5s, raws, strict=True):
                download(url, md5, raw, connections=connections)
            downsample(raws[0], raws[1], outs[0].parent, target_bases=target_bases)
            for raw in raws:
                raw.unlink()
        return outs
    fraction = sample_fraction(library.bases, target_bases)
    for url, md5, out in zip(library.urls, library.md5s, outs, strict=True):
        if not out.exists():
            raw = folder / f".raw-{out.name}"
            fetch_one(url, md5, raw, out, fraction=fraction, connections=connections)
    return outs


def _fetch_command(args: argparse.Namespace) -> int:
    paths = load_paths()
    record = read_record(paths.lab, args.library)
    if record is not None:
        path = library_dir(paths.lab, args.library) / ALIGNED
        if on_this_reference(paths.lab, record):
            print(
                f"{label(args.library)}: aligned on this reference ({path}); nothing fetched"
            )
        else:
            print(
                f"{label(args.library)}: aligned on another reference ({path}); nothing"
                f" fetched: {realign_hint(paths.lab, args.library)}"
            )
        return 0
    target = recorded_bases(paths.reports)
    if target is None:
        raise RunRefused("no recorded bases for your run: run genokeeper run first")
    fetch(
        paths.lab,
        SAMPLE_LIBRARIES[args.library],
        target_bases=target,
        connections=args.connections,
    )
    return 0


def _calls(key: str) -> bool:
    """Whether DeepVariant calls `key`'s library (HG003 is called; HG002 pcr is aligned
    only)."""
    return key == "hg003" or SAMPLE_LIBRARIES[key].calls


def _align(
    paths: Paths,
    key: str,
    pair: tuple[Path, Path],
    *,
    threads: int | None,
    from_cram: bool,
) -> Record:
    """`key`'s pair through `run_sample` when its library is called, else `align_sample`,
    under its work folder (`realign.work_name`) as its sample's sex; the record of it, on
    this reference, accepted unless it is HG004's (`accept hg004` scores that one first)."""
    name, sex, gvcf = realign.work_name(key), library_sex(key), None
    if _calls(key):
        run = run_sample(
            paths,
            name=name,
            fastq1=pair[0],
            fastq2=pair[1],
            sex=sex,
            threads=threads,
            from_cram=from_cram,
        )
        cram, cram_key, metrics, gvcf = run.cram, run.key, run.metrics, str(run.gvcf)
    else:
        cram, cram_key, metrics = align_sample(
            paths,
            name=name,
            fastq1=pair[0],
            fastq2=pair[1],
            sex=sex,
            threads=threads,
            from_cram=from_cram,
        )
    return Record(
        cram=str(cram),
        key=cram_key,
        metrics=metrics,
        gvcf=gvcf,
        reference_sha256=reference_sha256(paths.lab),
        accepted=key != "hg004",
        snp_f1=None,
        replaced=None,
    )


def _align_command(args: argparse.Namespace) -> int:
    paths = load_paths()
    key = args.library
    old = read_record(paths.lab, key)
    if old is not None and on_this_reference(paths.lab, old):
        path = library_dir(paths.lab, key) / ALIGNED
        print(f"{label(key)}: aligned on this reference ({path}); nothing aligned")
        return 0
    pair = downsampled_pair(paths.lab, key)
    if not all(p.is_file() for p in pair):
        fix = (
            f"run `python -m genokeeper.calibrate.hg002 fetch {key}`"
            if old is None
            else realign_hint(paths.lab, key)  # `fetch` skips any record
        )
        raise RunRefused(
            f"no downsampled pair for {label(key)} in {pair[0].parent}; {fix}"
        )
    record = _align(paths, key, pair, threads=None, from_cram=False)
    write_record(paths.lab, key, record)
    return 0


def _realign_command(args: argparse.Namespace) -> int:
    """`key` re-aligned from the CRAM its record names (`realign.old_record`): its reads
    written back as a pair, counted against the old CRAM, aligned (and called) on this
    reference, the new CRAM's bases counted against the old CRAM's again, and the record
    rewritten naming what it replaced. On a refusal the old record and CRAM stay as they
    are; a pair already written is reused once its counts match."""
    paths = load_paths()
    lab, key = paths.lab, args.library
    old = realign.old_record(lab, key)
    if on_this_reference(lab, old):
        print(
            f"{label(key)}: aligned on this reference"
            f" ({library_dir(lab, key) / ALIGNED}); nothing to re-align"
        )
        return 0
    calling_inputs(lab)  # the masked FASTA and its index refuse here, not hours later
    pair = realign.hg003_pair(lab) if key == "hg003" else downsampled_pair(lab, key)
    threads = args.threads or os.cpu_count() or 4
    cram, fasta = Path(old.cram), unmasked_fasta(lab)
    expected = realign.cram_counts(cram, fasta, threads=threads)
    if not all(p.is_file() for p in pair):
        realign.cram_to_fastq(cram, fasta, pair, threads=threads)
    try:
        realign.check_counts(expected, realign.fastq_counts(pair), label(key))
    except RunRefused:
        for p in pair:
            p.unlink(missing_ok=True)
        raise
    record = _align(paths, key, pair, threads=args.threads, from_cram=True)
    if record.metrics["total_bases"] != expected[1]:
        raise RunRefused(
            f"{label(key)}: the new CRAM {record.cram} holds"
            f" {record.metrics['total_bases']:.0f} bases, the old CRAM {expected[1]}; the"
            " old record is kept"
        )
    replaced = Replaced(cram=str(cram.parent), pair=[str(p) for p in pair])
    write_record(lab, key, replace(record, replaced=replaced))
    print(
        f"{label(key)}: re-aligned on this reference into {record.cram}; `clean {key}`"
        f" deletes {cram.parent} and the pair"
    )
    return 0


def _footprint(path: Path) -> int:
    """The bytes of the file or folder at `path`."""
    if path.is_dir():
        return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
    return path.stat().st_size


def _clean_command(args: argparse.Namespace) -> int:
    """What a current record replaced deleted: the old CRAM's step folder and the pair it
    was written into (HG002, HG003), or HG004's downsampled pair (it replaced nothing);
    the GB freed printed. Refuses a record on another reference, and a record whose
    `replaced` holds its own CRAM."""
    lab, key = load_paths().lab, args.library
    record = read_record(lab, key)
    if record is None or not on_this_reference(lab, record):
        fix = "align" if key == "hg004" else "realign"
        raise RunRefused(
            f"{label(key)} is not aligned on this reference; run `python -m"
            f" genokeeper.calibrate.hg002 {fix} {key}` first: clean deletes only what a"
            " current record replaced"
        )
    if key == "hg004":
        targets = list(downsampled_pair(lab, key))
    elif record.replaced is not None:
        targets = [Path(record.replaced["cram"]), *map(Path, record.replaced["pair"])]
    else:
        targets = []
    own = Path(record.cram)
    for target in targets:
        if target == own or target in own.parents:
            raise RunRefused(
                f"{label(key)}'s record names {target}, which holds its own CRAM {own},"
                " as replaced; nothing deleted"
            )
    there = [t for t in targets if t.exists()]
    if not there:
        print(f"{label(key)}: nothing left to clean")
        return 0
    freed = sum(_footprint(t) for t in there)
    for target in there:
        shutil.rmtree(target) if target.is_dir() else target.unlink()
    print(
        f"{label(key)}: {freed / 1e9:.1f} GB freed ({', '.join(str(t) for t in there)})"
    )
    return 0


# Where `hg004_snp_f1` writes the benchmark it scores over: HG004's BED minus the mask.
HG004_SCORED_BED = "hg004_benchmark_unmasked.bed"


def hg004_snp_f1(lab: Path, gvcf: Path) -> float:
    """The genotype-exact SNP F1 of HG004's gVCF against GIAB's HG004 v4.2.1 truth
    (`count_cells`, on `layout(lab)["sdf"]`), over its benchmark BED minus the masked
    intervals (`load_mask()`), written under HG004's record folder. The strata do not enter
    the F1, so none are read. `RunRefused` naming the bootstrap on a missing input."""
    vcf, bench = hg004_truth(lab)
    sdf = layout(lab)["sdf"]
    for what, path, flag in (
        ("GIAB's HG004 v4.2.1 truth", vcf, "--with-hg004"),
        ("its index", Path(f"{vcf}.tbi"), "--with-hg004"),
        ("its benchmark BED", bench, "--with-hg004"),
        ("rtg's SDF of the GRCh38 FASTA", sdf, "--with-masked-reference"),
    ):
        if not path.exists():
            raise RunRefused(
                f"missing input: {what} at {path}; run scripts/bootstrap.sh {flag}"
            )
    scored = library_dir(lab, "hg004") / HG004_SCORED_BED
    tmp = scored.with_name(f".{scored.name}.{secrets.token_hex(4)}")
    tmp.write_text(
        "".join(
            f"{c}\t{s}\t{e}\n" for c, s, e in subtract(read_bed(bench), load_mask())
        )
    )
    os.replace(tmp, scored)
    _, f1 = count_cells(
        gvcf, vcf, scored, sdf=sdf, strata=StrataIndex({}), rtg=rtg_path(lab)
    )
    return f1["snp"]


def _accept_command(args: argparse.Namespace) -> int:
    """HG004's record rewritten with its SNP F1 and `accepted` (at least
    `HG004_MIN_SNP_F1`); exit 2 when not accepted. Run again, it scores again."""
    lab = load_paths().lab
    record = read_record(lab, args.library)
    if record is None or record.gvcf is None or not Path(record.gvcf).is_file():
        raise RunRefused(
            f"{label(args.library)} is not called; run"
            f" `python -m genokeeper.calibrate.hg002 align {args.library}`"
        )
    if not on_this_reference(lab, record):
        raise RunRefused(
            f"{label(args.library)} is aligned on another reference;"
            f" {realign_hint(lab, args.library)}"
        )
    f1 = hg004_snp_f1(lab, Path(record.gvcf))
    accepted = f1 >= HG004_MIN_SNP_F1
    write_record(lab, args.library, replace(record, snp_f1=f1, accepted=accepted))
    if not accepted:
        print(
            f"hg002 accept: HG004's calls do not match its truth (SNP F1 {f1:.2f});"
            " check that the run is HG004's",
            file=sys.stderr,
        )
        return 2
    print(
        f"HG004 accepted: SNP F1 {f1:.4f} against its truth, {HG004_MIN_SNP_F1} needed"
    )
    return 0


def _layout_command(_args: argparse.Namespace) -> int:
    """`key<TAB>sample<TAB>aligned.json<TAB>pair 1<TAB>pair 2`, one `RECORDED` key per
    line; `-` for a pair with no fixed name (HG003's)."""
    lab = load_paths().lab
    for key in RECORDED:
        if key in SAMPLE_LIBRARIES:
            sample = SAMPLE_LIBRARIES[key].sample
            first, second = (str(p) for p in downsampled_pair(lab, key))
        else:
            sample, first, second = label(key), "-", "-"
        print(f"{key}\t{sample}\t{library_dir(lab, key) / ALIGNED}\t{first}\t{second}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m genokeeper.calibrate.hg002")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, handler, keys in (
        ("fetch", _fetch_command, SAMPLE_LIBRARIES),
        ("align", _align_command, SAMPLE_LIBRARIES),
        ("realign", _realign_command, [*LIBRARIES, "hg003"]),
        ("clean", _clean_command, RECORDED),
    ):
        p = sub.add_parser(name)
        p.add_argument("library", choices=sorted(keys))
        if name == "fetch":
            p.add_argument("--connections", type=int, default=24)
        if name == "realign":
            p.add_argument("--threads", type=int, help="default: this machine's CPUs")
        p.set_defaults(handler=handler)
    accept = sub.add_parser("accept")
    accept.add_argument("library", choices=["hg004"])
    accept.set_defaults(handler=_accept_command)
    sub.add_parser("layout").set_defaults(handler=_layout_command)
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (RunRefused, ValueError, ToolFailed, ToolMissing) as e:
        print(f"hg002 {args.command}: {' '.join(str(e).split())}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
