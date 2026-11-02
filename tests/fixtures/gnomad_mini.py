"""A miniature gnomAD v4.1 joint sites VCF, a GRCh38 reference it agrees with, and a fake
bucket that serves both the VCF and its object metadata to `build_chromosome` through the
`opener` seam, recording every request. Everything is written at test time.
"""

import base64
import hashlib
import io
import json
import random
import urllib.error
from dataclasses import dataclass
from email.message import Message
from pathlib import Path
from urllib.request import Request

import pysam

from genokeeper.evidence.frequency import FAF_FIELD, gnomad_url, metadata_url

ANC_FIELD = "fafmax_faf95_max_gen_anc_joint"
CHROMS = ("chr21", "chr22")
LENGTH = 300
# A C at 100, then six As (101-106) and a G: a deletion of one A anywhere in the run
# left-aligns to 100 CA>C.
RUN_START = 101


def sequence(chrom: str) -> str:
    """The synthetic contig: random bases from a seed of the contig's name, with the A run."""
    bases = random.Random(chrom).choices("ACGT", k=LENGTH)
    bases[RUN_START - 2] = "C"
    bases[RUN_START - 1 : RUN_START + 5] = "A" * 6
    bases[RUN_START + 5] = "G"
    return "".join(bases)


def base(chrom: str, pos: int) -> str:
    return sequence(chrom)[pos - 1]


def other(chrom: str, pos: int) -> str:
    """An ALT base that differs from the reference at `pos`."""
    return "T" if base(chrom, pos) != "T" else "G"


@dataclass(frozen=True)
class Site:
    """One joint VCF record; `faf=None` leaves the field missing."""

    chrom: str
    pos: int
    ref: str
    alt: str
    faf: float | None
    filter: str = "PASS"
    gen_anc: str = "nfe"


def snv(chrom: str, pos: int, faf: float | None, **kw: str) -> Site:
    return Site(chrom, pos, base(chrom, pos), other(chrom, pos), faf, **kw)


def write_reference(path: Path, chroms: tuple[str, ...] = CHROMS) -> Path:
    """The contigs as an indexed FASTA at `path`."""
    path.write_text("".join(f">{c}\n{sequence(c)}\n" for c in chroms))
    pysam.faidx(str(path))
    return path


def vcf_bytes(tmp: Path, chrom: str, sites: list[Site]) -> bytes:
    """The bgzipped joint sites VCF of `sites` (sorted here), as the bucket would serve it."""
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom},length={LENGTH}>",
        *(
            f'##FILTER=<ID={f},Description="{f}">'
            for f in ("EXOMES_FILTERED", "GENOMES_FILTERED", "BOTH_FILTERED")
        ),
        f'##INFO=<ID={FAF_FIELD},Number=A,Type=Float,Description="faf95 max">',
        f'##INFO=<ID={ANC_FIELD},Number=A,Type=String,Description="its group">',
        '##INFO=<ID=vep,Number=.,Type=String,Description="Allele|Consequence|SYMBOL">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    rows = []
    for s in sorted(sites, key=lambda s: s.pos):
        info = [
            f"{ANC_FIELD}={s.gen_anc}",
            f"vep={s.alt}|missense_variant|GENE1",
        ]
        if s.faf is not None:
            info.insert(0, f"{FAF_FIELD}={s.faf!r}")
        rows.append(
            f"{s.chrom}\t{s.pos}\t.\t{s.ref}\t{s.alt}\t.\t{s.filter}\t{';'.join(info)}"
        )
    plain = tmp / f"{chrom}.vcf"
    plain.write_text("\n".join([*header, *rows]) + "\n")
    packed = tmp / f"{chrom}.vcf.bgz"
    pysam.tabix_compress(str(plain), str(packed), force=True)
    return packed.read_bytes()


def md5_b64(data: bytes) -> str:
    """The md5 as the bucket's JSON API gives it: base64 of the digest."""
    return base64.b64encode(hashlib.md5(data).digest()).decode()


class _Served(io.BytesIO):
    """A response serving `data[start:stop]`, then raising `drop` (None: it just ends);
    `delivered` counts the bytes read from it."""

    def __init__(
        self,
        data: bytes,
        start: int,
        stop: int,
        drop: type[Exception] | None,
        status: int,
        headers: Message,
    ) -> None:
        super().__init__(data[start : max(start, stop)])
        self.drop, self.status, self.headers, self.delivered = drop, status, headers, 0

    def read(self, size: int | None = -1) -> bytes:
        chunk = super().read(size)
        self.delivered += len(chunk)
        if not chunk and self.drop is not None:
            raise self.drop("the connection broke")
        return chunk


type Fault = tuple[int, type[Exception] | None]


class Bucket:
    """A fake `urlopen` for the gnomAD bucket: serves `files[chrom]` at `gnomad_url(chrom)`
    and its size and md5 at `metadata_url(chrom)`; 404 for any other URL. Every request is
    kept in `requests`, and each file response beside it in `log`. A `Range: bytes=<n>-`
    gets 206 from byte n, or 416 past the end; `answer` other than "right" answers it
    wrongly: "ignored" (200, the whole file), "elsewhere" (206 from byte n - 1) or
    "other-total" (206 naming another size). A chromosome in `cut` drops every response
    (ConnectionResetError) at that byte; `faults[chrom]` lists one (byte, exception or None
    for a clean early end) per response, used up in order."""

    def __init__(self, files: dict[str, bytes]) -> None:
        self.files = files
        self.requests: list[Request] = []
        self.log: list[tuple[Request, _Served]] = []
        self.cut: dict[str, int] = {}
        self.faults: dict[str, list[Fault]] = {}
        self.answer = "right"

    def __call__(self, request: Request, timeout: float | None = None) -> io.BytesIO:
        self.requests.append(request)
        for chrom, data in self.files.items():
            if request.full_url == gnomad_url(chrom):
                return self._file(request, chrom, data)
            if request.full_url == metadata_url(chrom):
                meta = {"size": str(len(data)), "md5Hash": md5_b64(data)}
                return io.BytesIO(json.dumps(meta).encode())
        raise urllib.error.HTTPError(
            request.full_url, 404, "Not Found", Message(), None
        )

    def _file(self, request: Request, chrom: str, data: bytes) -> _Served:
        headers, status, start = Message(), 200, 0
        asked = request.get_header("Range")
        if asked is not None and self.answer != "ignored":
            start = int(asked.removeprefix("bytes=").removesuffix("-"))
            if start >= len(data):
                raise urllib.error.HTTPError(
                    request.full_url, 416, "Range Not Satisfiable", Message(), None
                )
            start -= self.answer == "elsewhere"
            total = len(data) + (self.answer == "other-total")
            headers["Content-Range"] = f"bytes {start}-{len(data) - 1}/{total}"
            status = 206
        stop, drop = len(data), None
        if chrom in self.cut:
            stop, drop = self.cut[chrom], ConnectionResetError
        elif self.faults.get(chrom):
            stop, drop = self.faults[chrom].pop(0)
        served = _Served(data, start, stop, drop, status, headers)
        self.log.append((request, served))
        return served
