"""File checksums for run inputs and artifact caches."""

import hashlib
import json
import os
from pathlib import Path

_CHUNK = 1 << 20


def sha256_file(path: Path) -> str:
    """The hex sha256 of the file at `path`, read in 1 MiB chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    """The hex md5 of the file at `path`, read in 1 MiB chunks: an identity the calibration
    compares (the checkpoint's and reference's md5, as their publishers list them), not a
    security check."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as f:
        while chunk := f.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def md5_cached(path: Path, cache: Path) -> str:
    """`md5_file(path)`, kept in `cache` beside the file's sha256 as `sha256_cached` keeps
    that (and on the same anchor)."""
    return _cached_digest(path, cache, "md5")


def sha256_cached(path: Path, cache: Path) -> str:
    """`sha256_file(path)`, kept in the JSON file `cache` under the resolved path with the
    file's size, mtime and inode, and read from there while all three are unchanged.

    A delivered FASTQ is about 100 GB: hashing it whole takes minutes, every run. A file
    replaced by another (a new inode), rewritten (another mtime) or truncated (another size)
    is hashed again; so is one the cache does not know. The cache holds only file metadata
    and digests, written whole through a temp file and a rename.

    Its anchor comes from the filesystem's size, mtime and inode of the file, which a
    delivery copied over, re-downloaded or rewritten cannot leave all unchanged, because a
    new file takes a new inode and a write sets the mtime; only a deliberate reset of the
    mtime on an in-place rewrite of equal size defeats it, which no step of this lab does."""
    return _cached_digest(path, cache, "sha256")


def _cached_digest(path: Path, cache: Path, algorithm: str) -> str:
    """The `algorithm` ("sha256" or "md5") digest of `path` through `cache`: an entry whose
    size, mtime and inode still match keeps its digests, and gains one it lacks; any other
    entry is replaced."""
    stat = path.stat()
    key = str(path.resolve())
    known: dict[str, dict[str, object]] = (
        json.loads(cache.read_text()) if cache.is_file() else {}
    )
    entry: dict[str, object] = {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }
    found = known.get(key)
    if found is not None and {k: found.get(k) for k in entry} == entry:
        if algorithm in found:
            return str(found[algorithm])
        entry = found
    digest = sha256_file(path) if algorithm == "sha256" else md5_file(path)
    known[key] = {**entry, algorithm: digest}
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_name(f".{cache.name}.{os.getpid()}")
    tmp.write_text(json.dumps(known, indent=1, sort_keys=True))
    tmp.replace(cache)
    return digest
