"""SMN1/SMN2 copy number from two callers (slice 3a): SMNCopyNumberCaller and Parascopy.

Each runs offline in its image (`containers.IMAGES`, built from docker/) behind the core gate
as `dna:smn`, on the CRAM, its index and the FASTA (with its `.fai`) mounted read-only and one
output directory mounted read-write; Parascopy also reads its data from `$lab/parascopy/`
(`parascopy_data`). The readers give the results as the grading reads them.

Parascopy's two steps are separate (slice 3d): its depth (`parascopy_depth`) is computed once
per CRAM and read by every copy-number step (`parascopy_cn`), 3a's SMN1 model and 3d's paralog
models alike.
"""

import gzip
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genokeeper.containers import run_container

SMNCALLER_VERSION = "1.1.2"
PARASCOPY_VERSION = "1.19.0"
LABEL = "dna:smn"
# Parascopy's names for the homology table and the SMN1 model, inside its unpacked data
# (`scripts/bootstrap.sh`, step_parascopy).
_TABLE = Path("homology_table") / "GRCh38.bed.gz"
_MODELS = Path("models_GRCh38_1KGP") / "EUR"
_MODEL = "SMN1.gz"
# The sample name Parascopy gives the CRAM, as its read group does.
SAMPLE = "sample"
# SMN1's and SMN2's exon 7 c.840 sites (GRCh38, 1-based; SMNCopyNumberCaller's
# SMN_SNP_38.txt): the Parascopy row holding a site is the gene's.
SMN1_SITE = ("chr5", 70951946)
SMN2_SITE = ("chr5", 70076526)


@dataclass(frozen=True)
class SmnCaller:
    """SMNCopyNumberCaller's copy numbers (None: not reported), its count of the g.27134T>G
    silent-carrier marker and its own confidence."""

    smn1: int | None
    smn2: int | None
    g27134tg: int | None
    confidence: str | None


@dataclass(frozen=True)
class ParascopyCall:
    """Parascopy's copy numbers at SMN1 and SMN2 (None: not reported), its quality and
    whether the call passed its filters."""

    smn1: int | None
    smn2: int | None
    qual: float | None
    passed: bool


def parascopy_data(lab: Path) -> tuple[Path, Path]:
    """The homology table and the EUR models folder in the lab's unpacked Parascopy data;
    `FileNotFoundError` names the bootstrap when either is absent."""
    root = lab / "parascopy"
    table, models = root / _TABLE, root / _MODELS
    for needed in (table, models / _MODEL):
        if not needed.is_file():
            raise FileNotFoundError(f"{needed}: run scripts/bootstrap.sh")
    return table, models


def parascopy_md5s(lab: Path) -> str:
    """The Parascopy data's identity: the homology table's own md5 list (`GRCh38.md5` beside
    it) and each tarball's md5 the bootstrap checked before unpacking it (its `.unpacked.*`
    marker); `FileNotFoundError` names the bootstrap when one is absent."""
    root = lab / "parascopy"
    listed = root / _TABLE.with_name("GRCh38.md5")
    markers = sorted(root.glob(".unpacked.*"))
    if not listed.is_file() or not markers:
        raise FileNotFoundError(
            f"{root}: no Parascopy data md5s; run scripts/bootstrap.sh"
        )
    parts = [" ".join(listed.read_text().split())]
    parts += [
        f"{m.name.removeprefix('.unpacked.')}:{m.read_text().strip()}" for m in markers
    ]
    return "; ".join(parts)


def _inputs(cram: Path, fasta: Path, out_dir: Path) -> list[tuple[Path, str, bool]]:
    return [
        (cram.parent, "/cram", True),
        (fasta.parent, "/ref", True),
        (out_dir, "/out", False),
    ]


def run_smn_caller(cram: Path, fasta: Path, out_dir: Path, *, threads: int) -> Path:
    """SMNCopyNumberCaller on `cram` (GRCh38), its JSON written as `out_dir/smn.json`."""
    (out_dir / "manifest.txt").write_text(f"/cram/{cram.name}\n")
    run_container(
        "smncaller",
        [
            "smn_caller.py", "--manifest", "/out/manifest.txt", "--genome", "38",
            "--reference", f"/ref/{fasta.name}", "--prefix", "smn",
            "--outDir", "/out", "--threads", str(threads),
        ],
        mounts=_inputs(cram, fasta, out_dir),
        label=LABEL,
    )  # fmt: skip
    return out_dir / "smn.json"


def _count(value: Any) -> int | None:
    return None if value is None or value == "None" else int(value)


def read_smn_caller(path: Path) -> SmnCaller:
    """The one sample's copy numbers in SMNCopyNumberCaller's JSON at `path`; the text
    `None` (a count it did not give) reads as None. Its confidence is the `Info` it gives
    the call."""
    (sample,) = json.loads(path.read_text()).values()
    info = sample.get("Info")
    return SmnCaller(
        smn1=_count(sample.get("SMN1")),
        smn2=_count(sample.get("SMN2")),
        g27134tg=_count(sample.get("g27134TG_CN")),
        confidence=None if info in (None, "None") else str(info),
    )


def _reads(cram: Path, fasta: Path, threads: int) -> list[str]:
    return [
        "-i", f"/cram/{cram.name}::{SAMPLE}", "-f", f"/ref/{fasta.name}",
        "-@", str(threads),
    ]  # fmt: skip


def parascopy_depth(cram: Path, fasta: Path, out_dir: Path, *, threads: int) -> Path:
    """Parascopy's background depth of `cram` (GRCh38, `-g hg38`) in `out_dir/depth`, behind
    `dna:smn`; returns its `depth.csv`."""
    run_container(
        "parascopy",
        [
            "parascopy", "depth", *_reads(cram, fasta, threads), "-g", "hg38",
            "-o", "/out/depth",
        ],
        mounts=_inputs(cram, fasta, out_dir),
        label=LABEL,
    )  # fmt: skip
    return out_dir / "depth" / "depth.csv"


def parascopy_cn(
    cram: Path,
    fasta: Path,
    depth: Path,
    out_dir: Path,
    *,
    lab: Path,
    threads: int,
    models: Sequence[str],
    label: str,
) -> Path:
    """Parascopy's copy number on `cram` with each of the EUR `models` (file stems, `SMN1`,
    `STRC`, ...) in one `cn-using` behind `label`, on `depth` (`parascopy_depth`'s file,
    another step's: its folder mounted read-only as `/depth`), in `out_dir/cn`; returns the
    paralog-specific table `out_dir/cn/res.paralog.bed.gz`. A model absent from the lab's
    data is refused, naming the bootstrap."""
    table, folder = parascopy_data(lab)
    root = lab / "parascopy"
    paths = [folder / f"{m}.gz" for m in models]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"{path}: run scripts/bootstrap.sh")
    run_container(
        "parascopy",
        [
            "parascopy", "cn-using",
            *(str(Path("/parascopy") / p.relative_to(root)) for p in paths),
            *_reads(cram, fasta, threads),
            "-t", str(Path("/parascopy") / table.relative_to(root)),
            "-d", f"/depth/{depth.name}", "-o", "/out/cn",
        ],
        mounts=[
            *_inputs(cram, fasta, out_dir),
            (root, "/parascopy", True),
            (depth.parent, "/depth", True),
        ],
        label=label,
    )  # fmt: skip
    return out_dir / "cn" / "res.paralog.bed.gz"


def _row_at(rows: list[list[str]], site: tuple[str, int]) -> list[str] | None:
    """The first row of Parascopy's table holding `site` (1-based; a row's start is 0-based)."""
    chrom, pos = site
    return next(
        (r for r in rows if r[0] == chrom and int(r[1]) < pos <= int(r[2])), None
    )


def known(value: str) -> str | None:
    """A cell of Parascopy's table, None for its `?`."""
    return None if value == "?" else value


def paralog_rows(path: Path) -> list[list[str]]:
    """The rows of Parascopy's paralog table at `path` (gzip when it ends in `.gz`), its
    `#` header lines skipped. Columns: chrom start end sample cn_filter cn cn_qual
    agcn_filter agcn agcn_qual main_region (`?` for an unknown value)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]


def read_parascopy(path: Path) -> ParascopyCall:
    """Parascopy's copy numbers at SMN1 and SMN2 in its paralog table at `path`
    (`paralog_rows`): each gene's row is the one holding its exon 7 site. `?` and a missing
    row read as None; `qual` and `passed` are the SMN1 row's."""
    rows = paralog_rows(path)
    smn1, smn2 = _row_at(rows, SMN1_SITE), _row_at(rows, SMN2_SITE)
    # Columns: chrom start end sample cn_filter cn cn_qual agcn_filter agcn agcn_qual ...
    cn1 = None if smn1 is None else known(smn1[5])
    cn2 = None if smn2 is None else known(smn2[5])
    qual = None if smn1 is None else known(smn1[6])
    return ParascopyCall(
        smn1=None if cn1 is None else int(cn1),
        smn2=None if cn2 is None else int(cn2),
        qual=None if qual is None else float(qual),
        passed=smn1 is not None and smn1[4] == "PASS",
    )
