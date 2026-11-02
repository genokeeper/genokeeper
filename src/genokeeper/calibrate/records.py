"""The calibration samples' alignment records and their folders: one `aligned.json` per
library key (`RECORDED`), which `genokeeper.calibrate.hg002` writes and every calibration reads
through `read_aligned` and `read_called`.

A record names the CRAM, its align step's key and QC metrics, the gVCF of a called library,
and the sha256 of the FASTA the CRAM was written against; a record on another FASTA (3a's,
written before the masked reference, carries none) or a library not accepted against its
truth (HG004) is refused, so no calibration reads it. This module imports no CLI:
`genokeeper.calibrate.hg002` re-exports its names, and the re-alignment imports it without
importing `hg002`.
"""

import json
import os
import secrets
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TypedDict

from genokeeper.calibrate import calibration_dir, hg003_dir
from genokeeper.checksum import sha256_cached
from genokeeper.config import grch38_fasta
from genokeeper.pipeline import IDENTITIES, RunRefused

ALIGNED = "aligned.json"
MARGIN_BYTES = 20 * 10**9  # free disk kept beyond what a file needs
# Every key with a record, in the order `hg002 layout` lists them.
RECORDED = ("pcrfree", "pcr", "hg003", "hg004")
# HG004's run must reach this SNP F1 against HG004's own truth before any calibration reads
# it: another sample's calls (ENA's mislabelled trio runs) would score far below.
HG004_MIN_SNP_F1 = 0.99
_HG002 = ("pcrfree", "pcr")
_CLI = "python -m genokeeper.calibrate.hg002"


def hg002_dir(lab: Path) -> Path:
    """HG002's downsampled FASTQ and alignment records, under the calibration folder."""
    return calibration_dir(lab) / "hg002"


def library_dir(lab: Path, key: str) -> Path:
    """The folder holding `key`'s record, downsampled pair and steps: HG002's libraries
    under `hg002_dir`, HG004's PCR-free library under `hg004/pcrfree`, HG003's own folder."""
    if key in _HG002:
        return hg002_dir(lab) / key
    if key == "hg004":
        return calibration_dir(lab) / "hg004" / "pcrfree"
    if key == "hg003":
        return hg003_dir(lab)
    raise ValueError(f"no calibration library {key!r}; one of {', '.join(RECORDED)}")


def label(key: str) -> str:
    """The library as the messages name it: `HG002 pcr`, `HG003`, `HG004`."""
    return f"HG002 {key}" if key in _HG002 else key.upper()


def need_space(folder: Path, need: int, what: str) -> None:
    """`RunRefused` unless `folder`'s disk has `need` bytes and `MARGIN_BYTES` to spare."""
    free = shutil.disk_usage(folder).free
    if free < need + MARGIN_BYTES:
        raise RunRefused(
            f"{what} needs {(need + MARGIN_BYTES) // 10**9} GB free in {folder},"
            f" {free // 10**9} GB is"
        )


class Replaced(TypedDict):
    """What a re-alignment replaced: the old CRAM and the pair made from it."""

    cram: str
    pair: list[str]


@dataclass(frozen=True)
class Record:
    """One library's `aligned.json`. `reference_sha256` is None only in 3a's legacy shape,
    which never names the FASTA it was written against."""

    cram: str
    key: str
    metrics: dict[str, float]
    gvcf: str | None
    reference_sha256: str | None
    accepted: bool
    snp_f1: float | None
    replaced: Replaced | None


def write_record(lab: Path, key: str, record: Record) -> Path:
    """`record` as `key`'s `aligned.json`, written whole through a hidden temporary and a
    rename."""
    path = library_dir(lab, key) / ALIGNED
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{ALIGNED}.{secrets.token_hex(4)}")
    tmp.write_text(json.dumps(asdict(record)))
    os.replace(tmp, path)
    return path


def read_record(lab: Path, key: str) -> Record | None:
    """`key`'s record, or None when it has none. 3a's shape, `{"cram", "key", "metrics"}`,
    reads with no gVCF, no reference digest, accepted, no F1 and nothing replaced."""
    path = library_dir(lab, key) / ALIGNED
    if not path.is_file():
        return None
    data = json.loads(path.read_text())
    return Record(
        cram=data["cram"],
        key=data["key"],
        metrics=data["metrics"],
        gvcf=data.get("gvcf"),
        reference_sha256=data.get("reference_sha256"),
        accepted=data.get("accepted", True),
        snp_f1=data.get("snp_f1"),
        replaced=data.get("replaced"),
    )


def reference_sha256(lab: Path) -> str:
    """The sha256 of the lab's GRCh38 FASTA (`grch38_fasta`), the digest a current record
    carries; `RunRefused` naming the bootstrap when the FASTA is not there."""
    fasta = grch38_fasta(lab)
    if not fasta.is_file():
        raise RunRefused(
            f"missing input: the masked GRCh38 FASTA at {fasta}; run"
            " scripts/bootstrap.sh --with-masked-reference"
        )
    return sha256_cached(fasta, calibration_dir(lab) / IDENTITIES)


def on_this_reference(lab: Path, record: Record) -> bool:
    """Whether `record`'s CRAM was written against the lab's GRCh38 FASTA."""
    return record.reference_sha256 is not None and record.reference_sha256 == (
        reference_sha256(lab)
    )


def _first(key: str) -> str:
    """The command that makes `key`'s first record: HG003's is re-aligned from its CRAM."""
    return f"`{_CLI} realign hg003`" if key == "hg003" else f"`{_CLI} align {key}`"


def aside(lab: Path, key: str) -> str:
    """Remake `key`'s record from scratch: `align` and `realign` keep a current one."""
    return f"move {library_dir(lab, key) / ALIGNED} aside, then run {_first(key)}"


def realign_hint(lab: Path, key: str) -> str:
    """The command that puts `key`'s record on this reference: `realign`, which HG004 has
    not (its record moved aside, `align hg004` aligns its pair, `fetch hg004` first)."""
    return aside(lab, key) if key == "hg004" else f"run `{_CLI} realign {key}`"


def _current(lab: Path, key: str) -> Record:
    """`key`'s record when a calibration may read it: aligned, on this reference and
    accepted; `RunRefused` naming the command that makes it so."""
    record = read_record(lab, key)
    if record is None:
        raise RunRefused(f"{label(key)} is not aligned; run {_first(key)}")
    if not on_this_reference(lab, record):
        raise RunRefused(
            f"{label(key)} is aligned on another reference; {realign_hint(lab, key)}"
        )
    if not record.accepted:
        if record.snp_f1 is None:
            raise RunRefused(
                f"{label(key)} is not accepted: not yet scored against its truth; run"
                f" `{_CLI} accept {key}`"
            )
        raise RunRefused(
            f"{label(key)} is not accepted: SNP F1 {record.snp_f1:.2f} against its"
            f" truth, below {HG004_MIN_SNP_F1}"
        )
    if not Path(record.cram).is_file():
        raise RunRefused(
            f"{label(key)}'s CRAM {record.cram} is gone; {aside(lab, key)}"
        )
    return record


def read_aligned(lab: Path, key: str) -> tuple[Path, str, dict[str, float]]:
    """`key`'s CRAM, its align step's key and QC metrics, as `align` recorded them;
    `RunRefused` when it is not aligned (or its CRAM is gone), aligned on another
    reference, or not accepted, naming the command that fixes it."""
    record = _current(lab, key)
    return Path(record.cram), record.key, record.metrics


def read_called(lab: Path, key: str) -> Path:
    """The gVCF of `key`'s record, refused as `read_aligned` refuses, and when the record
    has no gVCF (or it is gone), naming `align`."""
    record = _current(lab, key)
    if record.gvcf is None or not Path(record.gvcf).is_file():
        raise RunRefused(f"{label(key)} is not called; {aside(lab, key)}")
    return Path(record.gvcf)
