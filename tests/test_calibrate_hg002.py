"""HG002's two public libraries: streamed whole through their md5 while downsampled on the
fly (parallel byte ranges, fed in order), landing nothing when a check fails."""

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from genokeeper import tools
from genokeeper.calibrate import (
    build,
    calibration_dir,
    hg002,
    hg003_dir,
    hg004_truth,
    layout,
)
from genokeeper.calibrate.hg002 import (
    HG004,
    LIBRARIES,
    RECORDED,
    SAMPLE_LIBRARIES,
    Library,
    download,
    downsampled_pair,
    fetch,
    fetch_one,
    hg002_dir,
    library_dir,
)
from genokeeper.calibrate.records import (
    ALIGNED,
    Record,
    read_aligned,
    read_called,
    read_record,
    write_record,
)
from genokeeper.calibrate.sample import SampleRun
from genokeeper.checksum import sha256_file
from genokeeper.config import Paths, grch38_fasta, unmasked_fasta
from genokeeper.pipeline import RunRefused
from genokeeper.regions import read_bed
from genokeeper.tools import ToolFailed

SRA = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
CNGB = "https://ftp.cngb.org/pub/CNSA/data7/CNP0008319/CNS1414629/CNX1267210/CNR1441269"


@pytest.fixture(autouse=True)
def _no_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")


def _fastq(path: Path, mate: int, reads: int = 200) -> bytes:
    with gzip.open(path, "wt") as f:
        for i in range(reads):
            f.write(f"@read{i}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    return path.read_bytes()


def _names(fastq: Path) -> list[str]:
    with gzip.open(fastq, "rt") as f:
        return [line.split("/")[0] for line in f.read().splitlines()[0::4]]


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def _source(tmp_path: Path, name: str = "r_1.fastq.gz", reads: int = 200) -> Path:
    fq = tmp_path / "src" / name
    fq.parent.mkdir(exist_ok=True)
    _fastq(fq, 1, reads)
    return fq


# A curl that logs each ranged request (its first byte), fails the one starting at
# $FAIL_FROM (as a dropped connection does), and otherwise runs the real curl.
SHIM = """#!/bin/bash
prev=; first=
for a; do
  [ "$prev" = -r ] && first=${a%%-*}
  prev=$a
done
if [ -n "$first" ]; then
  echo "$first" >> "$CURL_LOG"
  [ "$first" = "${FAIL_FROM:-x}" ] && exit 22
fi
exec REAL_CURL "$@"
"""


def _shim(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    real = shutil.which("curl")
    assert real is not None
    shim = tmp_path / "shim" / "curl"
    shim.parent.mkdir()
    shim.write_text(SHIM.replace("REAL_CURL", real))
    shim.chmod(0o755)
    log = tmp_path / "curl.log"
    monkeypatch.setenv("PATH", f"{shim.parent}:{os.environ['PATH']}")
    monkeypatch.setenv("CURL_LOG", str(log))
    return log


def test_download_resumes_and_fetches_only_missing_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    log = _shim(tmp_path, monkeypatch)
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 400)
    fq = _source(tmp_path, reads=2000)
    size = fq.stat().st_size
    chunks = -(-size // 400)
    assert chunks > 6
    out = tmp_path / "out" / ".raw"
    monkeypatch.setenv("FAIL_FROM", str(3 * 400))
    with pytest.raises(ToolFailed):
        download(fq.as_uri(), _md5(fq), out, connections=1)
    assert not out.exists()
    journal = out.with_name(".raw.journal").read_text().splitlines()
    assert journal == [f"{size} 400", "0", "1", "2"]  # one connection: in order
    assert out.with_name(".raw.part").stat().st_size == size
    log.write_text("")
    monkeypatch.delenv("FAIL_FROM")
    download(fq.as_uri(), _md5(fq), out, connections=4)
    assert sorted(int(x) for x in log.read_text().split()) == [
        i * 400 for i in range(3, chunks)
    ]
    assert out.read_bytes() == fq.read_bytes()
    assert list(out.parent.iterdir()) == [out]  # part, journal and chunks are gone
    log.write_text("")
    download(
        fq.as_uri(), "0" * 32, out
    )  # present: verified when it landed, not refetched
    assert log.read_text() == ""


# Serves $SERVE_FILE as an http file: `-I` answers its Content-Length, a ranged request writes
# those bytes and, as $MODE says, the headers of a real 206, of a 200, or of another range.
# The range from $INJECT_FROM (logged to $INJECT_LOG) gets $INJECT_XML over the middle of its
# body, on its first request or, with INJECT=always, on every one. The range from
# $REDIRECT_FROM (logged to $REDIRECT_LOG) is answered with a 301 and its short body twice
# or, with REDIRECT=always, every time.
HTTP_CURL = """#!/bin/bash
range=; out=; hdr=; head=0; prev=
for a; do
  case "$prev" in -r) range=$a ;; -o) out=$a ;; -D) hdr=$a ;; esac
  [ "$a" = -fsSIL ] && head=1
  prev=$a
done
size=$(wc -c < "$SERVE_FILE" | tr -d ' ')
if [ "$head" = 1 ]; then printf 'HTTP/1.1 200 OK\\r\\nContent-Length: %s\\r\\n\\r\\n' "$size"; exit 0; fi
a=${range%-*}; b=${range#*-}
if [ "$a" = "${REDIRECT_FROM:-x}" ]; then
  echo "$a" >> "$REDIRECT_LOG"
  if [ "$REDIRECT" = always ] || [ "$(wc -l < "$REDIRECT_LOG")" -le 2 ]; then
    printf moved > "$out"
    printf 'HTTP/1.1 301 Moved Permanently\\r\\nContent-Length: 5\\r\\n\\r\\n' > "$hdr"
    exit 0
  fi
fi
tail -c +$((a + 1)) "$SERVE_FILE" | head -c $((b - a + 1)) > "$out"
if [ "$a" = "${INJECT_FROM:-x}" ]; then
  echo "$a" >> "$INJECT_LOG"
  if [ "$INJECT" = always ] || [ "$(wc -l < "$INJECT_LOG")" -eq 1 ]; then
    printf %s "$INJECT_XML" | dd of="$out" bs=1 seek=$(((b - a + 1 - ${#INJECT_XML}) / 2)) conv=notrunc 2>/dev/null
  fi
fi
case "$MODE" in
  ok) printf 'HTTP/1.1 206 Partial Content\\r\\nContent-Range: bytes %s-%s/%s\\r\\n\\r\\n' $a $b $size ;;
  full) printf 'HTTP/1.1 200 OK\\r\\nContent-Length: %s\\r\\n\\r\\n' $((b - a + 1)) ;;
  shifted) printf 'HTTP/1.1 206 Partial Content\\r\\nContent-Range: bytes %s-%s/%s\\r\\n\\r\\n' $((a + 1)) $((b + 1)) $size ;;
esac > "$hdr"
"""


def _http(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fq: Path, mode: str) -> None:
    fake = tmp_path / "fakebin" / "curl"
    fake.parent.mkdir()
    fake.write_text(HTTP_CURL)
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake.parent}:{os.environ['PATH']}")
    monkeypatch.setenv("SERVE_FILE", str(fq))
    monkeypatch.setenv("MODE", mode)


@pytest.mark.parametrize("mode", ["ok", "full", "shifted"])
def test_download_requires_the_range_asked_for(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    """A chunk of the right size that the server answered with a 200, or for another range,
    is refused and not journaled."""
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 400)
    fq = _source(tmp_path, reads=500)
    _http(tmp_path, monkeypatch, fq, mode)
    out = tmp_path / "out" / ".raw"
    url = "https://example.test/r_1.fastq.gz"
    if mode == "ok":
        assert (
            download(url, _md5(fq), out, connections=3).read_bytes() == fq.read_bytes()
        )
        return
    with pytest.raises(ValueError, match=r"bytes \d+-\d+: (HTTP 200|Content-Range)"):
        download(url, _md5(fq), out, connections=1)
    assert out.with_name(".raw.journal").read_text().splitlines()[1:] == []
    assert not out.exists()


@pytest.mark.parametrize("redirect", ["twice", "always"])
def test_download_refetches_a_redirected_chunk(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, redirect: str
) -> None:
    """A range answered with a redirect's short body (CNGB's HTTP) is fetched again; after
    20 short answers the chunk is refused, never written or journaled."""
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 1000)
    fq = _source(tmp_path, reads=5000)
    _http(tmp_path, monkeypatch, fq, "ok")
    log = tmp_path / "redirect.log"
    for key, value in (("REDIRECT", redirect), ("REDIRECT_FROM", "2000")):
        monkeypatch.setenv(key, value)
    monkeypatch.setenv("REDIRECT_LOG", str(log))
    out = tmp_path / "out" / ".raw"
    url = "http://example.test/r_1.fastq.gz"
    if redirect == "twice":
        assert (
            download(url, _md5(fq), out, connections=1).read_bytes() == fq.read_bytes()
        )
        assert log.read_text().split() == ["2000"] * 3
        return
    with pytest.raises(
        ValueError, match=r"bytes 2000-2999: .*another size in 20 tries"
    ):
        download(url, _md5(fq), out, connections=1)
    assert log.read_text().split() == ["2000"] * 20
    assert out.with_name(".raw.journal").read_text().split("\n")[1:-1] == ["0", "1"]


# What ENA's storage gateway writes over ~202 bytes of a 206 body when it loses its backend.
GATEWAY_ERROR = (
    "<Error><Code>ConnectionClosedException</Code><Message>Premature end of"
    " Content-Length delimited message body (expected: 97,332,938,379; received:"
    " 14,225,573)</Message><ErrorMessage/><RequestId/></Error>"
)


@pytest.mark.parametrize("inject", ["once", "always"])
def test_download_refetches_a_chunk_holding_the_gateway_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, inject: str
) -> None:
    """A right-sized 206 for the range asked for whose body holds the gateway's error is
    fetched again; after 5 such answers the chunk is refused, never written or journaled."""
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 1000)
    fq = _source(tmp_path, reads=5000)
    _http(tmp_path, monkeypatch, fq, "ok")
    log = tmp_path / "inject.log"
    for key, value in (("INJECT", inject), ("INJECT_FROM", "2000")):
        monkeypatch.setenv(key, value)
    monkeypatch.setenv("INJECT_LOG", str(log))
    monkeypatch.setenv("INJECT_XML", GATEWAY_ERROR)
    out = tmp_path / "out" / ".raw"
    url = "https://example.test/r_1.fastq.gz"
    if inject == "once":
        assert (
            download(url, _md5(fq), out, connections=1).read_bytes() == fq.read_bytes()
        )
        assert log.read_text().split() == ["2000", "2000"]
        return
    with pytest.raises(ValueError, match=r"bytes 2000-2999 \(chunk 2\).*5 tries"):
        download(url, _md5(fq), out, connections=1)
    assert log.read_text().split() == ["2000"] * 5
    assert out.with_name(".raw.journal").read_text().split("\n")[1:-1] == ["0", "1"]
    assert out.with_name(".raw.part").read_bytes()[2000:3000] == bytes(1000)


def test_download_refetches_journaled_chunks_holding_the_gateway_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A part whose journaled chunk 2 holds the error (written by earlier code): chunk 2 is
    taken out of the journal and is the only one fetched again."""
    log = _shim(tmp_path, monkeypatch)
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 1000)
    fq = _source(tmp_path, reads=5000)
    size = fq.stat().st_size
    chunks = -(-size // 1000)
    data = bytearray(fq.read_bytes())
    data[2400 : 2400 + len(GATEWAY_ERROR)] = GATEWAY_ERROR.encode()
    out = tmp_path / "out" / ".raw"
    out.parent.mkdir()
    out.with_name(".raw.part").write_bytes(bytes(data))
    journal = out.with_name(".raw.journal")
    journal.write_text(f"{size} 1000\n" + "".join(f"{i}\n" for i in range(chunks)))
    monkeypatch.setenv("FAIL_FROM", "2000")
    with pytest.raises(ToolFailed):
        download(fq.as_uri(), _md5(fq), out, connections=1)
    assert journal.read_text().split("\n")[1:-1] == [
        str(i) for i in range(chunks) if i != 2
    ]
    monkeypatch.delenv("FAIL_FROM")
    log.write_text("")
    assert download(fq.as_uri(), _md5(fq), out).read_bytes() == fq.read_bytes()
    assert log.read_text().split() == ["2000"]


def test_download_survives_a_short_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 400)
    fq = _source(tmp_path, reads=500)
    real = os.pwrite
    calls: list[int] = []

    def half_first(fd: int, data: bytes, offset: int) -> int:
        calls.append(len(data))
        return (
            real(fd, data[: len(data) // 2], offset)
            if len(calls) == 1
            else real(fd, data, offset)
        )

    monkeypatch.setattr(os, "pwrite", half_first)
    out = download(fq.as_uri(), _md5(fq), tmp_path / "out" / ".raw", connections=1)
    assert out.read_bytes() == fq.read_bytes()
    assert len(calls) > 1 and calls[1] == calls[0] - calls[0] // 2


def test_download_starts_over_for_another_size(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 400)
    fq = _source(tmp_path, reads=500)
    out = tmp_path / "out" / ".raw"
    out.parent.mkdir()
    out.with_name(".raw.part").write_bytes(b"x" * 10)
    out.with_name(".raw.journal").write_text("10 400\n0\n")
    (out.parent / ".raw.chunks").mkdir()
    download(fq.as_uri(), _md5(fq), out, connections=3)
    assert out.read_bytes() == fq.read_bytes()


def test_download_md5_mismatch_keeps_the_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(hg002, "CHUNK_BYTES", 400)
    fq = _source(tmp_path, reads=500)
    out = tmp_path / "out" / ".raw"
    with pytest.raises(
        ValueError, match=rf"md5 mismatch.*{_md5(fq)}, not {'0' * 32}.*\d+ chunks"
    ):
        download(fq.as_uri(), "0" * 32, out, connections=2)
    # The journal goes, the file stays for inspection, under a hidden name.
    kept = out.with_name(f".raw.mismatch-{_md5(fq)}")
    assert [p.name for p in out.parent.iterdir()] == [kept.name]
    assert kept.read_bytes() == fq.read_bytes()
    err = capsys.readouterr().err
    assert "0 of " in err and "chunks already in the journal" in err and "MB/s" in err
    # A later download removes the leftover first (its disk counts), and lands the file.
    download(fq.as_uri(), _md5(fq), out, connections=2)
    assert [p.name for p in out.parent.iterdir()] == [".raw"]


def test_a_test_never_reads_the_hosts_free_disk(tmp_path: Path) -> None:
    """The download and re-alignment tests reach `need_space`; each sees conftest's stub
    (`ample_free_disk`), never the free space of the machine running it, which on a small
    CI runner is under the 20 GB margin."""
    assert shutil.disk_usage(tmp_path).free == 10**12


def test_download_refuses_without_free_space(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fq = _source(tmp_path)
    monkeypatch.setattr(hg002.shutil, "disk_usage", lambda _p: SimpleNamespace(free=5))
    out = tmp_path / "out" / ".raw"
    with pytest.raises(RunRefused, match=r"needs \d+ GB free"):
        download(fq.as_uri(), _md5(fq), out)
    assert list(out.parent.iterdir()) == []
    raw = tmp_path / "raw.fastq.gz"
    raw.write_bytes(fq.read_bytes())
    with pytest.raises(RunRefused, match="free"):
        hg002.downsample_fastq(raw, tmp_path / "ds" / "o.gz", fraction=0.5)
    assert raw.exists() and list((tmp_path / "ds").iterdir()) == []


@pytest.mark.tools
def test_fetch_one_keeps_pairs_and_deletes_raw(tmp_path: Path) -> None:
    outs = []
    for mate in (1, 2):
        fq = tmp_path / "src" / f"r_{mate}.fastq.gz"
        fq.parent.mkdir(exist_ok=True)
        _fastq(fq, mate)
        folder = tmp_path / "lib"
        out = folder / "downsampled" / fq.name
        stale = out.parent / f".{out.name}.deadbeef"  # an earlier design's temporary
        (stale.parent).mkdir(parents=True, exist_ok=True)
        stale.write_bytes(b"x")
        (out.parent / f".{out.name}.deadbeef.chunks").mkdir()
        outs.append(
            fetch_one(
                fq.as_uri(), _md5(fq), folder / f".raw-{fq.name}", out, fraction=0.5
            )
        )
    kept = _names(outs[0])
    assert 60 <= len(kept) <= 140
    assert kept == _names(outs[1])
    assert sorted(p.name for p in (tmp_path / "lib").iterdir()) == ["downsampled"]
    assert sorted(p.name for p in outs[0].parent.iterdir()) == [
        "r_1.fastq.gz",
        "r_2.fastq.gz",
    ]


@pytest.mark.tools
def test_cut_stream_lands_nothing_and_keeps_the_raw(tmp_path: Path) -> None:
    fq = tmp_path / "src" / "r_1.fastq.gz"
    fq.parent.mkdir()
    data = _fastq(fq, 1)
    fq.write_bytes(data[: len(data) // 2])  # a gzip cut short, whose md5 is then right
    raw = tmp_path / "lib" / ".raw-r_1.fastq.gz"
    out = tmp_path / "lib" / "downsampled" / fq.name
    with pytest.raises(ToolFailed):
        fetch_one(fq.as_uri(), _md5(fq), raw, out, fraction=0.5)
    assert not out.exists()
    assert list(out.parent.iterdir()) == []  # no hidden temporary
    assert (
        raw.read_bytes() == fq.read_bytes()
    )  # complete: a rerun redoes the downsample


@pytest.mark.tools
def test_chunks_and_connections_give_one_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fq = _source(tmp_path, reads=2000)
    results = []
    for name, chunk, connections in (("one", 2**30, 1), ("many", 997, 8)):
        monkeypatch.setattr(hg002, "CHUNK_BYTES", chunk)
        out = tmp_path / name / "downsampled" / fq.name
        fetch_one(
            fq.as_uri(),
            _md5(fq),
            tmp_path / name / ".raw",
            out,
            fraction=0.5,
            connections=connections,
        )
        results.append(out.read_bytes())
    assert results[0] == results[1]


def _library(tmp_path: Path) -> Library:
    src = tmp_path / "src"
    src.mkdir(exist_ok=True)
    urls = ((src / "lib_1.fastq.gz").as_uri(), (src / "lib_2.fastq.gz").as_uri())
    return Library("pcr", "SRRX", urls, ("a" * 32, "b" * 32), bases=1000)


def test_fetch_skips_present(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    lab = tmp_path / "lab"
    lib = _library(tmp_path)
    monkeypatch.setitem(SAMPLE_LIBRARIES, "pcr", lib)
    first, second = downsampled_pair(lab, "pcr")
    first.parent.mkdir(parents=True)
    first.write_bytes(b"there")
    streamed: list[tuple[str, str, Path, Path, float]] = []

    def fake(
        url: str, md5: str, raw: Path, out: Path, *, fraction: float, **_: object
    ) -> Path:
        streamed.append((url, md5, raw, out, fraction))
        out.write_bytes(b"new")
        return out

    monkeypatch.setattr(hg002, "fetch_one", fake)
    got = fetch(lab, lib, target_bases=250)
    assert got == (first, second)
    raw = hg002_dir(lab) / "pcr" / ".raw-lib_2.fastq.gz"
    assert streamed == [(lib.urls[1], lib.md5s[1], raw, second, 0.25)]
    assert first.read_bytes() == b"there"


def test_fetch_counts_a_library_without_bases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No published base count: both raw files land in raw/, the pair is downsampled
    from their counted bases, then both raw files go; a rerun fetches nothing."""
    lab = tmp_path / "lab"
    lib = replace(_library(tmp_path), bases=None)
    monkeypatch.setitem(SAMPLE_LIBRARIES, "pcr", lib)
    outs = downsampled_pair(lab, "pcr")
    raw = hg002_dir(lab) / "pcr" / "raw"
    got: list[tuple[str, str, Path]] = []

    def fake_download(url: str, md5: str, out: Path, **_: object) -> Path:
        got.append((url, md5, out))
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b"raw")
        return out

    def fake_downsample(
        first: Path, second: Path, out_dir: Path, *, target_bases: int
    ) -> tuple[Path, Path]:
        assert (first, second) == (raw / outs[0].name, raw / outs[1].name)
        assert (out_dir, target_bases) == (outs[0].parent, 250)
        out_dir.mkdir(parents=True, exist_ok=True)
        for out in outs:
            out.write_bytes(b"kept")
        return outs

    monkeypatch.setattr(hg002, "download", fake_download)
    monkeypatch.setattr(hg002, "downsample", fake_downsample)
    assert fetch(lab, lib, target_bases=250) == outs
    assert got == [
        (lib.urls[0], lib.md5s[0], raw / outs[0].name),
        (lib.urls[1], lib.md5s[1], raw / outs[1].name),
    ]
    assert list(raw.iterdir()) == []
    assert fetch(lab, lib, target_bases=250) == outs
    assert len(got) == 2


def test_fetch_skips_aligned_library(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = tmp_path / "lab"
    lib = _library(tmp_path)
    monkeypatch.setitem(SAMPLE_LIBRARIES, "pcr", lib)
    (hg002_dir(lab) / "pcr").mkdir(parents=True)
    # 3a's record: aligned on the unmasked reference; `realign` works from its CRAM.
    (hg002_dir(lab) / "pcr" / "aligned.json").write_text(
        json.dumps({"cram": "pcr.cram", "key": "k", "metrics": {}})
    )

    def refuse(*_a: object, **_k: object) -> Path:
        raise AssertionError("downloaded an aligned library again")

    monkeypatch.setattr(hg002, "fetch_one", refuse)
    fetch(lab, lib, target_bases=250)
    assert not downsampled_pair(lab, "pcr")[0].exists()


def test_align_without_the_pair_refuses_naming_fetch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    assert hg002.main(["align", "pcr"]) == 2
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1
    assert "hg002 fetch pcr" in err[0]


def test_libraries_pinned() -> None:
    assert LIBRARIES == {
        "pcrfree": Library(
            "pcrfree",
            "SRR38517577",
            (
                f"{SRA}/SRR385/077/SRR38517577/SRR38517577_1.fastq.gz",
                f"{SRA}/SRR385/077/SRR38517577/SRR38517577_2.fastq.gz",
            ),
            ("c018da93350a7b0313b9edb27ea56aa5", "0073c089fad9d9d9f4af2975389af638"),
            294758313900,
            calls=True,
        ),
        "pcr": Library(
            "pcr",
            "SRR38517578",
            (
                f"{SRA}/SRR385/078/SRR38517578/SRR38517578_1.fastq.gz",
                f"{SRA}/SRR385/078/SRR38517578/SRR38517578_2.fastq.gz",
            ),
            ("f5cffb64f16efeca3b06f602b1d0d171", "d4a2a7484c27a4829ee961eb3a0a22eb"),
            151918693321,
        ),
    }


# --- slice 3c: the masked reference and the library table ---------------------------------


def _fasta(lab: Path) -> Path:
    fasta = grch38_fasta(lab)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text(">chr1\nACGTNNNNAC\n")
    return fasta


def _current(lab: Path, key: str, **fields: object) -> Record:
    """A record of `key` on the lab's masked FASTA, its CRAM on disk, written."""
    folder = library_dir(lab, key)
    folder.mkdir(parents=True, exist_ok=True)
    cram = folder / f"{key}.cram"
    cram.write_bytes(b"")
    record = Record(
        cram=str(cram),
        key="k" * 8,
        metrics={"mean_depth": 30.0, "total_bases": 1.0e11},
        gvcf=None,
        reference_sha256=sha256_file(_fasta(lab)),
        accepted=True,
        snp_f1=None,
        replaced=None,
    )
    record = Record(**{**record.__dict__, **fields})
    write_record(lab, key, record)
    return record


def test_grch38_fasta_is_the_masked_set(tmp_path: Path) -> None:
    assert grch38_fasta(tmp_path) == (
        tmp_path / "reference" / "GRCh38_no_alt_maskedGRC.fa"
    )
    assert unmasked_fasta(tmp_path) == tmp_path / "reference" / "GRCh38_no_alt.fa"
    assert grch38_fasta(tmp_path) != unmasked_fasta(tmp_path)


def test_read_aligned_refuses_another_reference(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    record = _current(lab, "pcr", reference_sha256="0" * 64)
    assert read_record(lab, "pcr") == record
    with pytest.raises(RunRefused) as e:
        read_aligned(lab, "pcr")
    said = str(e.value)
    assert "HG002 pcr is aligned on another reference" in said
    assert "python -m genokeeper.calibrate.hg002 realign pcr" in said


def test_read_aligned_returns_a_current_record(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    record = _current(lab, "pcrfree")
    assert read_aligned(lab, "pcrfree") == (
        Path(record.cram),
        record.key,
        record.metrics,
    )
    with pytest.raises(RunRefused, match="hg002 align pcrfree"):
        read_called(lab, "pcrfree")  # aligned, never called
    gvcf = library_dir(lab, "pcrfree") / "calls.g.vcf.gz"
    gvcf.write_bytes(b"")
    _current(lab, "pcrfree", gvcf=str(gvcf))
    assert read_called(lab, "pcrfree") == gvcf


def test_legacy_record_reads_and_refuses(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    folder = hg002_dir(lab) / "pcrfree"
    folder.mkdir(parents=True)
    cram = folder / "hg002-pcrfree.cram"
    cram.write_bytes(b"")
    metrics = {"mean_depth": 30.0, "total_bases": 1.0e11}
    # As 3a's `align` wrote it, and as the lab holds it today.
    (folder / ALIGNED).write_text(
        json.dumps({"cram": str(cram), "key": "abc", "metrics": metrics})
    )
    assert read_record(lab, "pcrfree") == Record(
        cram=str(cram),
        key="abc",
        metrics=metrics,
        gvcf=None,
        reference_sha256=None,
        accepted=True,
        snp_f1=None,
        replaced=None,
    )
    with pytest.raises(RunRefused, match="hg002 realign pcrfree"):
        read_aligned(lab, "pcrfree")


def test_records_imports_no_cli() -> None:
    probe = (
        "import sys, genokeeper.calibrate.records;"
        " print('genokeeper.calibrate.hg002' in sys.modules)"
    )
    done = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    )
    assert done.stdout.strip() == "False"


def test_read_aligned_refuses_unaccepted_hg004(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    _current(lab, "hg004", accepted=False, snp_f1=0.71)
    with pytest.raises(RunRefused) as e:
        read_aligned(lab, "hg004")
    said = str(e.value)
    assert "HG004 is not accepted" in said
    assert "0.71" in said and "0.99" in said
    _current(lab, "hg004", accepted=False, snp_f1=None)
    with pytest.raises(RunRefused, match="hg002 accept hg004"):
        read_aligned(lab, "hg004")


def test_hg004_pinned() -> None:
    assert SAMPLE_LIBRARIES["hg004"] == HG004
    assert HG004 == Library(
        "hg004",
        "CNR1441269",
        (
            f"{CNGB}/E200029822_L01_UDB-452_1.fq.gz",
            f"{CNGB}/E200029822_L01_UDB-452_2.fq.gz",
        ),
        ("179cb9ee9f6baa77ba7e23ea472365f5", "1a37dc8d44c2ec77583eb7508a1b9d88"),
        None,
        "HG004",
        "female",
        True,
    )
    assert SAMPLE_LIBRARIES == {**LIBRARIES, "hg004": HG004}
    assert RECORDED == ("pcrfree", "pcr", "hg003", "hg004")
    assert [hg002.library_sex(k) for k in RECORDED] == [
        "male",
        "male",
        "male",
        "female",
    ]


def test_layout_lists_every_recorded_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    assert hg002.main(["layout"]) == 0
    rows = [line.split("\t") for line in capsys.readouterr().out.splitlines()]
    assert [(r[0], r[1]) for r in rows] == [
        ("pcrfree", "HG002"),
        ("pcr", "HG002"),
        ("hg003", "HG003"),
        ("hg004", "HG004"),
    ]
    assert rows[2] == ["hg003", "HG003", str(hg003_dir(lab) / ALIGNED), "-", "-"]
    hg004 = calibration_dir(lab) / "hg004" / "pcrfree"
    assert rows[3] == [
        "hg004",
        "HG004",
        str(hg004 / ALIGNED),
        str(hg004 / "downsampled" / "E200029822_L01_UDB-452_1.fq.gz"),
        str(hg004 / "downsampled" / "E200029822_L01_UDB-452_2.fq.gz"),
    ]
    assert rows[0][2] == str(hg002_dir(lab) / "pcrfree" / ALIGNED)


def test_fetch_skips_a_recorded_library(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    _current(lab, "hg004")

    def refuse(*_a: object, **_k: object) -> Path:
        raise AssertionError("downloaded a recorded library again")

    monkeypatch.setattr(hg002, "fetch_one", refuse)
    monkeypatch.setattr(hg002, "download", refuse)
    monkeypatch.setattr(hg002, "recorded_bases", refuse)  # never reads the reports
    assert hg002.main(["fetch", "hg004"]) == 0
    captured = capsys.readouterr()
    assert len((captured.out + captured.err).strip().splitlines()) == 1
    first, second = downsampled_pair(lab, "hg004")
    assert first.parent.parent == library_dir(lab, "hg004")
    assert second.parent == first.parent
    assert not first.exists() and not second.exists()


def test_fetch_before_a_run_says_not_yet(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Before `genokeeper run` recorded the bases, `fetch` exits `NOT_YET` (the bootstrap
    goes on to its later steps) and names the order that works, not a bare rerun."""
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path / "lab"))
    monkeypatch.setattr(hg002, "recorded_bases", lambda _reports: None)
    assert hg002.main(["fetch", "pcrfree"]) == hg002.NOT_YET
    said = capsys.readouterr().err
    assert "run genokeeper run, then rerun scripts/bootstrap.sh" in said


class _Stop(Exception):
    """Raised by a fake after the per-library loop: the calls before it are the test."""


def test_calibrations_take_each_library_sex(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = tmp_path / "lab"
    paths = Paths(repo=tmp_path, raw=tmp_path / "raw", reports=tmp_path / "r", lab=lab)
    _fasta(lab)
    monkeypatch.setattr(
        build, "LIBRARIES", {"pcrfree": LIBRARIES["pcrfree"], "hg004": HG004}
    )
    asked: list[str] = []

    def aligned(_lab: Path, library: str) -> tuple[Path, str, dict[str, float]]:
        asked.append(library)
        return tmp_path / f"{library}.cram", library, {"mean_depth": 30.0}

    monkeypatch.setattr(build, "read_aligned", aligned)
    sexes: dict[str, object] = {}

    def stage3(work: Path, **kwargs: object) -> tuple[object, dict[str, object]]:
        (key,) = [k for k in RECORDED if library_dir(lab, k) == work]
        sexes[key] = kwargs["sex"]
        found = SimpleNamespace(failed={}, flank=None, repeats={}, rfc1=None)
        return found, {"repeats": SimpleNamespace(key=key)}

    def stop(*_a: object, **_k: object) -> object:
        raise _Stop

    vcf, bench = build.hg002_truth(lab)
    vcf.parent.mkdir(parents=True)
    for p in (vcf, Path(f"{vcf}.tbi"), bench):
        p.write_bytes(b"")
    monkeypatch.setattr(build, "stage3_versions", lambda *_a: {})
    monkeypatch.setattr(build, "stage3_steps", stage3)
    monkeypatch.setattr(build, "repeat_steps", stage3)
    monkeypatch.setattr(build, "truth_repeats", stop)
    with pytest.raises(_Stop):
        build.build_stage3(paths, threads=1)
    assert sexes == {
        "pcrfree": "male",
        "pcr": "male",
        "hg003": "male",
        "hg004": "female",
    }

    males: list[object] = []

    def deletions(work: Path, **kwargs: object) -> tuple[object, None]:
        assert work == library_dir(lab, asked[-1])
        males.append(kwargs["male"])
        return SimpleNamespace(called=[], failed={}), None

    truth, sv_bench = build.truth_files(lab)
    for p in (truth, sv_bench, *build.sv_files(lab).required()):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"")
    monkeypatch.setattr(build, "deletion_versions", lambda *_a: {})
    monkeypatch.setattr(build, "read_deletions", lambda *_a, **_k: [])
    monkeypatch.setattr(build, "deletion_steps", deletions)
    monkeypatch.setattr(
        build, "load_screen", lambda *_a: SimpleNamespace(population=None)
    )
    monkeypatch.setattr(build, "judge_library", lambda *_a, **_k: ([], []))
    monkeypatch.setattr(build, "expected_true_carried", stop)
    with pytest.raises(_Stop):
        build.build_deletions(paths, threads=1)
    assert males == [True, False]


def test_accept_hg004_at_f1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """`accept hg004` scores HG004's gVCF against its truth over its BED minus the mask:
    below 0.99 the record says not accepted and the command exits 2 naming ENA's label;
    at 0.996 it is accepted and `read_aligned` returns it."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    gvcf = library_dir(lab, "hg004") / "hg004.g.vcf.gz"
    gvcf.parent.mkdir(parents=True)
    gvcf.write_bytes(b"")
    _current(lab, "hg004", gvcf=str(gvcf), accepted=False, snp_f1=None)
    vcf, bench = hg004_truth(lab)
    bench.parent.mkdir(parents=True)
    for p in (vcf, Path(f"{vcf}.tbi")):
        p.write_bytes(b"")
    bench.write_text("chr1\t0\t1000\nchr2\t0\t50\n")
    layout(lab)["sdf"].mkdir(parents=True)
    monkeypatch.setattr(hg002, "load_mask", lambda: [("chr1", 100, 200)])
    scored: list[list[tuple[str, int, int]]] = []
    snp = [0.71]

    def cells(
        called: Path, truth: Path, bed: Path, **_: object
    ) -> tuple[dict[str, object], dict[str, float]]:
        assert (called, truth) == (gvcf, vcf)
        scored.append(read_bed(bed))
        return {}, {"snp": snp[0], "indel": 0.9}

    monkeypatch.setattr(hg002, "count_cells", cells)
    assert hg002.main(["accept", "hg004"]) == 2
    assert "check that the run is HG004's" in capsys.readouterr().err
    record = read_record(lab, "hg004")
    assert record is not None
    assert (record.accepted, record.snp_f1) == (False, 0.71)
    assert scored == [[("chr1", 0, 100), ("chr1", 200, 1000), ("chr2", 0, 50)]]
    with pytest.raises(RunRefused, match="0.71"):
        read_aligned(lab, "hg004")
    snp[0] = 0.996
    assert hg002.main(["accept", "hg004"]) == 0
    record = read_record(lab, "hg004")
    assert record is not None
    assert (record.accepted, record.snp_f1) == (True, 0.996)
    assert read_aligned(lab, "hg004")[0] == Path(record.cram)


def test_hg004_align_records_unaccepted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`align hg004` calls HG004 under `hg004-pcrfree` as a female and records it not
    accepted, with no F1: no calibration reads it before `accept hg004` scores it."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    fasta = _fasta(lab)
    pair = downsampled_pair(lab, "hg004")
    pair[0].parent.mkdir(parents=True)
    for p in pair:
        p.write_bytes(b"reads")
    asked: list[dict[str, object]] = []
    cram, gvcf = tmp_path / "sample.cram", tmp_path / "sample.g.vcf.gz"
    for p in (cram, gvcf):
        p.write_bytes(b"")

    def run(_paths: object, **kwargs: object) -> SampleRun:
        asked.append(kwargs)
        return SampleRun(cram=cram, key="k", metrics={"total_bases": 9.0}, gvcf=gvcf)

    def never(*_a: object, **_k: object) -> object:
        raise AssertionError("HG004 is called, not only aligned")

    monkeypatch.setattr(hg002, "run_sample", run)
    monkeypatch.setattr(hg002, "align_sample", never)
    assert hg002.main(["align", "hg004"]) == 0
    (kwargs,) = asked
    assert (kwargs["name"], kwargs["sex"]) == ("hg004-pcrfree", "female")
    assert (kwargs["fastq1"], kwargs["fastq2"]) == pair
    assert read_record(lab, "hg004") == Record(
        cram=str(cram),
        key="k",
        metrics={"total_bases": 9.0},
        gvcf=str(gvcf),
        reference_sha256=sha256_file(fasta),
        accepted=False,
        snp_f1=None,
        replaced=None,
    )
    with pytest.raises(RunRefused, match="hg002 accept hg004"):
        read_aligned(lab, "hg004")


def test_align_reruns_nothing_on_a_current_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """`align` on a record already on this reference prints one line and exits 0, the
    record untouched: an accepted HG004 stays accepted."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    _current(lab, "hg004", gvcf=str(tmp_path / "x.g.vcf.gz"), snp_f1=0.996)
    before = (library_dir(lab, "hg004") / ALIGNED).read_bytes()

    def never(*_a: object, **_k: object) -> object:
        raise AssertionError("aligned again")

    monkeypatch.setattr(hg002, "run_sample", never)
    assert hg002.main(["align", "hg004"]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == 1
    assert (library_dir(lab, "hg004") / ALIGNED).read_bytes() == before


def _named(text: str) -> list[tuple[str, str]]:
    """The `hg002` commands `text` names, each proved to exist: its subcommand and key
    pass the parser (`--help` exits 0; a key outside its choices exits 2)."""
    found = re.findall(r"`python -m genokeeper\.calibrate\.hg002 (\w+) (\w+)`", text)
    assert found, text
    for command, key in found:
        with pytest.raises(SystemExit) as e:
            hg002.main([command, key, "--help"])
        assert e.value.code == 0, (command, key)
    return found


def test_hints_name_commands_that_exist(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Every hint names a command that exists and does something in that state."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    _fasta(lab)

    def refused(call: object, *args: object) -> str:
        with pytest.raises(RunRefused) as e:
            call(*args)  # type: ignore[operator]
        return str(e.value)

    # HG003 without a record, or current and not called: `realign hg003`, never `align`.
    assert _named(refused(read_aligned, lab, "hg003")) == [("realign", "hg003")]
    _current(lab, "hg003")
    said = refused(read_called, lab, "hg003")
    assert "aside" in said and _named(said) == [("realign", "hg003")]
    # HG004 on another reference: its record moved aside, then `align hg004`.
    gvcf = tmp_path / "x.g.vcf.gz"
    gvcf.write_bytes(b"")
    _current(lab, "hg004", gvcf=str(gvcf), reference_sha256="0" * 64)
    aside = f"move {library_dir(lab, 'hg004') / ALIGNED} aside, then run"
    for said in (
        refused(read_aligned, lab, "hg004"),
        refused(hg002._accept_command, SimpleNamespace(library="hg004")),
    ):
        assert aside in said and _named(said) == [("align", "hg004")]
    capsys.readouterr()
    assert hg002.main(["fetch", "hg004"]) == 0
    said = capsys.readouterr().out
    assert aside in said and _named(said) == [("align", "hg004")]
    # No pair and a record on another reference: `fetch` would fetch nothing.
    _current(lab, "pcr", reference_sha256="0" * 64)
    said = refused(hg002._align_command, SimpleNamespace(library="pcr"))
    assert _named(said) == [("realign", "pcr")]


def test_clean_removes_what_each_record_replaced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    # HG003, re-aligned: its old CRAM's step folder and the pair made from it go.
    old = hg003_dir(lab) / "align" / "old"
    old.mkdir(parents=True)
    (old / "sample.cram").write_bytes(b"x" * 1000)
    pair = [layout(lab)["downsampled"] / f"r_{i}.fq.gz" for i in (1, 2)]
    pair[0].parent.mkdir(parents=True)
    for p in pair:
        p.write_bytes(b"y" * 500)
    record = _current(
        lab, "hg003", replaced={"cram": str(old), "pair": [str(p) for p in pair]}
    )
    assert hg002.main(["clean", "hg003"]) == 0
    assert "GB freed" in capsys.readouterr().out
    assert not old.exists() and not any(p.exists() for p in pair)
    assert Path(record.cram).is_file()
    assert hg002.main(["clean", "hg003"]) == 0  # nothing left
    assert "nothing" in capsys.readouterr().out
    # HG004 replaced nothing: its downsampled pair goes, no CRAM is touched.
    hg004 = _current(lab, "hg004")
    for p in downsampled_pair(lab, "hg004"):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"z")
    assert hg002.main(["clean", "hg004"]) == 0
    assert not any(p.exists() for p in downsampled_pair(lab, "hg004"))
    assert Path(hg004.cram).is_file()
    # A record on another reference: refused, nothing deleted.
    kept = calibration_dir(lab) / "hg002-pcr" / "align" / "old"
    kept.mkdir(parents=True)
    _current(
        lab,
        "pcr",
        reference_sha256="0" * 64,
        replaced={"cram": str(kept), "pair": []},
    )
    assert hg002.main(["clean", "pcr"]) == 2
    assert "realign pcr" in capsys.readouterr().err
    assert kept.is_dir()
    # A record that names its own CRAM's folder as replaced: refused, nothing deleted.
    own = _current(lab, "pcrfree")
    write_record(
        lab,
        "pcrfree",
        Record(
            **{
                **own.__dict__,
                "replaced": {"cram": str(Path(own.cram).parent), "pair": []},
            }
        ),
    )
    assert hg002.main(["clean", "pcrfree"]) == 2
    assert Path(own.cram).is_file()
