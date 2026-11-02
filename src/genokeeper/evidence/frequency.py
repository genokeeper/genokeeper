"""Population frequency for the BA1/BS1 rules: gnomAD's filtering allele frequency table and
the ClinGen SVI BA1 exception list.

The table (`FafTable`) is the slim local copy of gnomAD v4.1's joint sites VCFs the bootstrap
builds under `<lab>/gnomad/4.1/` (`scripts/bootstrap.sh --with-gnomad`; track design,
A1): one `faf-<chrom>.parquet` per gnomAD chromosome, columns `vid VARCHAR`, `faf95_max
DOUBLE` (`fafmax_faf95_max_joint`) and `gen_anc VARCHAR` (`fafmax_faf95_max_gen_anc_joint`),
holding the variants at or above the manifest's `floor`, and a `manifest.json`:
`{"release": "4.1", "field": "fafmax_faf95_max_joint", "floor": <float>, "chromosomes":
{<chrom>: {"url", "size", "md5", "rows"}}}`. A vid absent from the table is below the floor
or absent from gnomAD: it gets no frequency. `FafTable.open` refuses a table whose manifest
does not list every gnomAD chromosome or whose file for one is missing; its anchor comes
from GNOMAD_CHROMS, the 24 chromosome names written in this module, which an interrupted or
partial build cannot move because the manifest must list each one and each file must exist.

`build_chromosome` (`genokeeper build-faf`) builds one chromosome's file: it streams the
bucket's whole `gnomad.joint.v4.1.sites.<chrom>.vcf.bgz` once, hashing every chunk while
`bcftools norm -f <the lab's GRCh38 FASTA> -m -both` and `bcftools query` keep the rows at or
above the floor (`faf_floor`: the lowest BS1 in config/grading.toml) whose FILTER is not
BOTH_FILTERED, and renames the file into place, rewriting the manifest, only when the stream's
size and md5 equal the bucket's object metadata. With the two-tier floor (`FAF_TWO_TIER`,
`TwoTier`), the rows below 0.01 are kept only inside the padded spans of the `[bs1_by_gene]`
genes (`resources/bs1_gene_spans.tsv`); the manifest records the tiers, and a table built
under another floor or other tiers is rebuilt whole. Every request is a GET of a URL
`gnomad_url` or `metadata_url` builds; the first for a URL carries no `Range` header, and
after a read error, a timeout, an early end or a 5xx or 429 answer the same process reopens
the same URL with exactly `Range: bytes=<N>-`, N being the bytes it has already received and
fed on (design decision, 2026-10-19). Its anchor comes from GNOMAD_CHROMS, the chromosome
names written in this module, and from the moment the network broke, which anyone reading the
bucket's access logs cannot move, because no input of the user's (genome, run or store)
reaches a URL, a header or N, and both functions refuse a name outside the list. A process
that dies restarts its chromosome from the first byte: the md5 state and the bcftools pipe
live only in it.

The exception list (`resources/ba1_exceptions.tsv`) names the alleles ClinGen exempts from
BA1 although they are common (Ghosh et al. 2018; HFE C282Y, BTD D444H among them); the grading
exempts them from BA1 and BS1 in every gene (track design, A2). Its anchor comes from the
ClinVar variation IDs the published table gives, resolved to ClinVar's own GRCh38 records,
which no input of a run can move because the list is a file in the package, read as it is.
"""

import argparse
import base64
import contextlib
import fcntl
import hashlib
import http.client
import json
import re
import shutil
import ssl
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, cast

import duckdb

from genokeeper.checksum import sha256_file
from genokeeper.config import Grading, load_grading, load_paths, unmasked_fasta
from genokeeper.tools import ToolFailed, ToolMissing, gated
from genokeeper.variants.key import make_vid

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

GNOMAD_RELEASE = "4.1"
FAF_FIELD = "fafmax_faf95_max_joint"
ANC_FIELD = "fafmax_faf95_max_gen_anc_joint"
GNOMAD_CHROMS = (*(f"chr{n}" for n in range(1, 23)), "chrX", "chrY")
GNOMAD_BUCKET = "gcp-public-data--gnomad"
_STORAGE = "https://storage.googleapis.com"
MANIFEST = "manifest.json"
_CHROM_KEYS = frozenset({"url", "size", "md5", "rows"})
# The manifest keys a chromosome's entry is valid under: one built under others is stale.
_POLICY_KEYS = ("release", "field", "floor", "two_tier")
# Track design A1: one tier, unless chr22's file projects the whole table above 20 GB; then
# True, and every chromosome is rebuilt under the two tiers.
FAF_TWO_TIER = False
TWO_TIER_GENOME_FLOOR = 0.01
BS1_GENE_SPANS_FILE = "bs1_gene_spans.tsv"
# Around each gene body: ClinVar files promoter and near-gene variants under the gene too.
SPAN_PAD = 10_000
_CHUNK = 1 << 20
# Seconds without a byte before a read fails: a healthy stream from the bucket never pauses a
# minute (the clock runs only while waiting on the socket, not while bcftools holds us up), and
# a stalled connection, ESTABLISHED at 0 B/s, is resumed after one.
_TIMEOUT = 60
# Resumes of one chromosome's stream, `_WAIT` x the attempt seconds apart (at most a minute):
# about a quarter of an hour of outage in all, before the chromosome fails.
_RESUMES = 20
_WAIT = 5
# What a broken connection raises; an early end is the other drop.
_DROPS = (
    TimeoutError,
    ConnectionError,
    ssl.SSLError,
    http.client.IncompleteRead,
    urllib.error.URLError,
)
_TRANSIENT = frozenset({429, *range(500, 600)})
BA1_EXCEPTIONS_FILE = "ba1_exceptions.tsv"
_COLUMNS = ["vid", "gene", "hgvs", "source"]


def parse_ba1_exceptions(text: str) -> frozenset[str]:
    """The vids of an exception list: `#` lines, a `vid gene hgvs source` header, one row
    per allele. Refuses a row of another width, a vid `make_vid` would not write, or a
    duplicate."""
    vids: set[str] = set()
    header = False
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if not header:
            if fields != _COLUMNS:
                raise ValueError(f"line {number}: expected the header {_COLUMNS}")
            header = True
            continue
        if len(fields) != len(_COLUMNS):
            raise ValueError(f"line {number}: expected {len(_COLUMNS)} columns")
        vid = fields[0]
        build, chrom, pos, ref, alt = vid.split(":")
        if (
            build != "GRCh38"
            or not pos.isdigit()
            or make_vid(chrom, int(pos), ref, alt) != vid
        ):
            raise ValueError(f"line {number}: {vid!r} is not a GRCh38 vid")
        if vid in vids:
            raise ValueError(f"line {number}: duplicate {vid}")
        vids.add(vid)
    if not vids:
        raise ValueError("no exception rows")
    return frozenset(vids)


def load_ba1_exceptions() -> frozenset[str]:
    """The GRCh38 vids of the packaged ClinGen SVI BA1 exception list."""
    text = files("genokeeper").joinpath("resources", BA1_EXCEPTIONS_FILE).read_text()
    return parse_ba1_exceptions(text)


def faf_folder(lab: Path) -> Path:
    """Where the lab keeps the FAF table: `<lab>/gnomad/4.1/`."""
    return lab / "gnomad" / GNOMAD_RELEASE


def faf_file(folder: Path, chrom: str) -> Path:
    """One chromosome's table file: `<folder>/faf-<chrom>.parquet`."""
    return folder / f"faf-{chrom}.parquet"


class FafMissing(RuntimeError):
    """The FAF table is absent or incomplete: a missing or unreadable manifest, one of
    another release, field or shape, a gnomAD chromosome it does not list, or a listed
    chromosome's missing file."""


@dataclass(frozen=True)
class FafTable:
    """An opened FAF table: its folder, its manifest and the manifest's sha256."""

    folder: Path
    manifest: dict[str, Any]
    sha256: str

    @classmethod
    def open(cls, folder: Path) -> FafTable:
        """Open the table in `folder`; raise `FafMissing` unless its manifest is this
        release's, lists every chromosome of GNOMAD_CHROMS and each one's file exists."""
        path = folder / MANIFEST
        try:
            loaded: object = json.loads(path.read_text())
        except (OSError, ValueError) as e:
            raise FafMissing(f"no readable manifest {path}: {e}") from e
        if not isinstance(loaded, dict):
            raise FafMissing(f"{path}: the manifest is not an object")
        data = cast(dict[str, Any], loaded)
        for key, value in (("release", GNOMAD_RELEASE), ("field", FAF_FIELD)):
            if data.get(key) != value:
                raise FafMissing(f"{path}: {key} {data.get(key)!r}, expected {value!r}")
        floor = data.get("floor")
        if isinstance(floor, bool) or not isinstance(floor, int | float):
            raise FafMissing(f"{path}: floor {floor!r} is not a number")
        listed = data.get("chromosomes")
        if not isinstance(listed, dict):
            raise FafMissing(f"{path}: chromosomes is not an object")
        chromosomes = cast(dict[str, Any], listed)
        missing = [
            c
            for c in GNOMAD_CHROMS
            if not isinstance(chromosomes.get(c), dict)
            or not _CHROM_KEYS <= set(chromosomes[c])
        ]
        if missing:
            raise FafMissing(f"{path} does not list {', '.join(missing)}")
        absent = [
            faf_file(folder, c).name
            for c in GNOMAD_CHROMS
            if not faf_file(folder, c).is_file()
        ]
        if absent:
            raise FafMissing(f"{folder} lacks {', '.join(absent)}")
        return cls(folder=folder, manifest=data, sha256=sha256_file(path))

    @property
    def manifest_path(self) -> Path:
        return self.folder / MANIFEST

    def lookup(self, vids: Collection[str]) -> dict[str, float]:
        """`faf95_max` of each of `vids` the table holds, in one DuckDB join of the vids
        against the 24 chromosome files; a vid it lacks is absent from the result."""
        if not vids:
            return {}
        paths = [str(faf_file(self.folder, c)) for c in GNOMAD_CHROMS]
        with duckdb.connect() as con:
            found = con.execute(
                "SELECT f.vid, max(f.faf95_max) FROM read_parquet(?) f"
                " JOIN (SELECT DISTINCT unnest(?::VARCHAR[]) AS vid) w ON w.vid = f.vid"
                " GROUP BY f.vid",
                [paths, list(vids)],
            ).fetchall()
        return {vid: float(faf) for vid, faf in found}

    def uncovered(self, grading: Grading) -> str | None:
        """Why the table misses a BS1 of `grading`, or None: a floor above the lowest
        (`faf_floor`), or two tiers other than `two_tier(grading)`'s."""
        floor, tiers = self.manifest["floor"], self.manifest.get("two_tier")
        if floor > faf_floor(grading):
            return f"floor {floor:g} is above the lowest BS1 {faf_floor(grading):g}"
        if tiers is not None and tiers != two_tier(grading).record():
            return "its two tiers are not the [bs1_by_gene] genes' spans"
        return None

    def describe(self) -> str:
        """One line naming the source, the field, the floor, the size and the manifest's
        sha256: no genotype, no position."""
        chromosomes = cast(dict[str, dict[str, Any]], self.manifest["chromosomes"])
        variants = sum(int(c["rows"]) for c in chromosomes.values())
        return (
            f"gnomAD v{self.manifest['release']} {self.manifest['field']}, floor"
            f" {self.manifest['floor']:g}, {variants} variants in {len(GNOMAD_CHROMS)}"
            f" chromosome files (manifest sha256 {self.sha256[:12]})"
        )


# --- building the table -----------------------------------------------------------------

# `urllib.request.urlopen` or a stand-in: (request, timeout=seconds) -> a readable response.
type Opener = Callable[..., Any]


class StreamMismatch(RuntimeError):
    """A stream whose size or md5 differs from the bucket's object metadata, or a resume the
    bucket did not answer from the bytes already received."""


def _object(chrom: str) -> str:
    if chrom not in GNOMAD_CHROMS:
        raise ValueError(
            f"{chrom!r} is not a gnomAD chromosome ({', '.join(GNOMAD_CHROMS)})"
        )
    return (
        f"release/{GNOMAD_RELEASE}/vcf/joint/"
        f"gnomad.joint.v{GNOMAD_RELEASE}.sites.{chrom}.vcf.bgz"
    )


def gnomad_url(chrom: str) -> str:
    """The whole joint sites VCF of one gnomAD chromosome in the public bucket."""
    return f"{_STORAGE}/{GNOMAD_BUCKET}/{_object(chrom)}"


def metadata_url(chrom: str) -> str:
    """The bucket's JSON API record of that file: its size and base64 md5 alone."""
    name = urllib.parse.quote(_object(chrom), safe="")
    return f"{_STORAGE}/storage/v1/b/{GNOMAD_BUCKET}/o/{name}?fields=size,md5Hash"


def _get(url: str, start: int | None = None) -> urllib.request.Request:
    """A GET of `url` with no body and no header of ours, or, resuming a stream after
    `start` bytes, exactly `Range: bytes=<start>-`."""
    headers = {} if start is None else {"Range": f"bytes={start}-"}
    return urllib.request.Request(url, headers=headers, method="GET")


def object_metadata(
    chrom: str, *, opener: Opener = urllib.request.urlopen
) -> tuple[int, str]:
    """The size and the base64 md5 the bucket lists for `gnomad_url(chrom)`."""
    with opener(_get(metadata_url(chrom)), timeout=_TIMEOUT) as response:
        data = json.load(response)
    try:
        return int(data["size"]), str(data["md5Hash"])
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(
            f"{metadata_url(chrom)}: no size and md5Hash in {data!r}"
        ) from e


def faf_floor(grading: Grading) -> float:
    """The table's floor: the lowest BS1 in the grading config, the fallback's and every
    gene's. A variant below it can trigger neither BA1 nor any BS1."""
    return min(grading.bs1, *grading.bs1_by_gene.values())


@dataclass(frozen=True)
class TwoTier:
    """Track design A1's two-tier floor: rows at or above `genome_floor` everywhere, and
    down to the table's floor only inside `spans` (gene -> chrom, start, end; 1-based,
    inclusive)."""

    genome_floor: float
    spans: Mapping[str, tuple[str, int, int]]

    def record(self) -> dict[str, Any]:
        """The manifest's `two_tier` value."""
        return {
            "genome_floor": self.genome_floor,
            "spans": {gene: list(span) for gene, span in sorted(self.spans.items())},
        }


def load_gene_spans() -> dict[str, tuple[str, int, int]]:
    """The packaged GRCh38 gene bodies (`resources/bs1_gene_spans.tsv`, NCBI Gene)."""
    text = files("genokeeper").joinpath("resources", BS1_GENE_SPANS_FILE).read_text()
    spans: dict[str, tuple[str, int, int]] = {}
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    for line in lines[1:]:
        gene, chrom, start, end, _ = line.split("\t")
        spans[gene] = (chrom, int(start), int(end))
    return spans


def two_tier(grading: Grading) -> TwoTier:
    """The two tiers for `grading`: 0.01 genome-wide, and the padded span of every
    `[bs1_by_gene]` gene; refuse a gene the span file lacks."""
    known = load_gene_spans()
    missing = sorted(set(grading.bs1_by_gene) - set(known))
    if missing:
        raise ValueError(f"{BS1_GENE_SPANS_FILE} has no span for {', '.join(missing)}")
    return TwoTier(
        genome_floor=TWO_TIER_GENOME_FLOOR,
        spans={
            gene: (
                known[gene][0],
                max(1, known[gene][1] - SPAN_PAD),
                known[gene][2] + SPAN_PAD,
            )
            for gene in sorted(grading.bs1_by_gene)
        },
    )


def _policy(floor: float, tiers: TwoTier | None) -> dict[str, Any]:
    policy: dict[str, Any] = {
        "release": GNOMAD_RELEASE,
        "field": FAF_FIELD,
        "floor": floor,
    }
    if tiers is not None:
        policy["two_tier"] = tiers.record()
    return policy


def _listed(folder: Path, policy: dict[str, Any]) -> dict[str, Any]:
    """The manifest's chromosome entries when it was built under `policy`, else none."""
    try:
        loaded: object = json.loads((folder / MANIFEST).read_text())
    except OSError, ValueError:
        return {}
    if not isinstance(loaded, dict):
        return {}
    data = cast(dict[str, Any], loaded)
    if any(data.get(k) != policy.get(k) for k in _POLICY_KEYS):
        return {}
    listed = data.get("chromosomes")
    return cast(dict[str, Any], listed) if isinstance(listed, dict) else {}


def stale_chromosomes(
    folder: Path,
    *,
    floor: float,
    two_tier: TwoTier | None = None,
    metadata: Mapping[str, tuple[int, str]] | None = None,
) -> list[str]:
    """The chromosomes to (re)build, in GNOMAD_CHROMS order: every one when the manifest is
    absent or was built under another floor or other tiers; else each it lacks or lists
    without its file, and, given the bucket's `metadata` (size and md5 per chromosome), each
    whose recorded size or md5 differs."""
    listed = _listed(folder, _policy(floor, two_tier))
    stale: list[str] = []
    for chrom in GNOMAD_CHROMS:
        entry = listed.get(chrom)
        if (
            not isinstance(entry, dict)
            or not _CHROM_KEYS <= set(cast(dict[str, Any], entry))
            or not faf_file(folder, chrom).is_file()
            or (
                metadata is not None
                and chrom in metadata
                and (entry["size"], entry["md5"]) != tuple(metadata[chrom])
            )
        ):
            stale.append(chrom)
    return stale


def build_chromosome(
    chrom: str,
    folder: Path,
    *,
    floor: float,
    fasta: Path,
    opener: Opener = urllib.request.urlopen,
    expected: tuple[int, str] | None = None,
    two_tier: TwoTier | None = None,
) -> int:
    """Stream `gnomad_url(chrom)` once into `folder/faf-<chrom>.parquet` and its manifest
    entry, and return the rows kept. `expected` is the bucket's (size, base64 md5), fetched
    when None; `fasta` is the GRCh38 FASTA the alleles are normalised against; the bcftools
    pipeline runs behind the core budget gate (`tools.gated`). A stream that is cut, or
    whose size or md5 differs (`StreamMismatch`), leaves no file and no entry; a chromosome
    already there stays as it was."""
    url = gnomad_url(chrom)
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        raise ToolMissing("bcftools")
    size, md5 = (
        expected if expected is not None else object_metadata(chrom, opener=opener)
    )
    folder.mkdir(parents=True, exist_ok=True)
    out = faf_file(folder, chrom)
    tmp = out.with_name(f"{out.name}.tmp")
    rows_tsv = out.with_name(f"{out.name}.tsv.tmp")
    try:
        got_size, got_md5 = _stream(
            url,
            rows_tsv,
            opener=opener,
            commands=_commands(bcftools, fasta, floor),
            total=size,
        )
        if got_size != size:
            raise StreamMismatch(
                f"{url}: {got_size} bytes, the bucket lists size {size}"
            )
        if got_md5 != md5:
            raise StreamMismatch(f"{url}: md5 {got_md5}, the bucket lists {md5}")
        rows = _write_parquet(rows_tsv, tmp, two_tier)
        tmp.replace(out)
    finally:
        rows_tsv.unlink(missing_ok=True)
        tmp.unlink(missing_ok=True)
    entry = {"url": url, "size": size, "md5": md5, "rows": rows}
    _record(folder, chrom, entry, _policy(floor, two_tier))
    return rows


def _commands(bcftools: str, fasta: Path, floor: float) -> tuple[list[str], list[str]]:
    """`bcftools norm` (split, left-align against `fasta`) and `bcftools query` (the kept
    rows as chrom, pos, ref, alt, faf95, its group).

    `-c x` leaves out a record whose REF disagrees with `fasta` (a gnomAD REF on an N of
    NCBI's no-alt FASTA), which would otherwise abort the whole chromosome's stream. It
    cannot silently empty the table: its anchor comes from the md5 the bootstrap checks the
    FASTA against (NCBI's md5checksums.txt), which a wrong or truncated reference cannot
    move, because the bootstrap refuses a FASTA whose md5 differs.

    The gate wraps `norm` only: the pipeline is one unit of work and takes one slot. Gating
    `query` too makes it wait for the slot `norm` holds while `norm` blocks on the full
    pipe, so the two deadlock until the gate's timeout. Its label is `dna:faf`, not
    `dna:annotate`: the stream is network-bound and its weight is one core."""
    keep = f'INFO/{FAF_FIELD}>={floor!r} && FILTER!="BOTH_FILTERED"'
    columns = f"%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/{FAF_FIELD}\\t%INFO/{ANC_FIELD}\\n"
    norm = ["norm", "-f", str(fasta), "-c", "x", "-m", "-both", "-Ou", "-"]
    return (
        gated("dna:faf", [bcftools, *norm]),
        [bcftools, "query", "-i", keep, "-f", columns, "-"],
    )


def _stream(
    url: str,
    sink: Path,
    *,
    opener: Opener,
    commands: tuple[list[str], list[str]],
    total: int,
) -> tuple[int, str]:
    """GET `url` whole into `norm | query > sink`, chunk by chunk; return the bytes read
    and their base64 md5. A connection that breaks, stalls, ends before `total` bytes or is
    answered 5xx or 429 is reopened from the bytes received, up to `_RESUMES` times. Raises
    what the stream raises, killing both tools, `StreamMismatch` for a wrong answer to a
    resume, or `ToolFailed` when either tool exits non-zero."""
    norm_cmd, query_cmd = commands
    digest = hashlib.md5(usedforsecurity=False)
    size = 0
    with (
        sink.open("wb") as out,
        tempfile.TemporaryFile() as norm_err,
        tempfile.TemporaryFile() as query_err,
    ):
        norm = subprocess.Popen(
            norm_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=norm_err
        )
        stdin, stdout = cast(IO[bytes], norm.stdin), cast(IO[bytes], norm.stdout)
        query = subprocess.Popen(query_cmd, stdin=stdout, stdout=out, stderr=query_err)
        stdout.close()  # query holds the read end
        resumes = 0
        try:
            while True:
                try:
                    request = _get(url, size if resumes else None)
                    with opener(request, timeout=_TIMEOUT) as response:
                        if resumes:
                            _check_resume(url, response, size, total)
                        while chunk := response.read(_CHUNK):
                            digest.update(chunk)
                            size += len(chunk)
                            stdin.write(chunk)
                    if size >= total:
                        break
                    raise ConnectionError(f"the stream ended at byte {size} of {total}")
                except BrokenPipeError:
                    raise
                except _DROPS as e:
                    # 5xx and 429 are transient, met as a drop; another HTTP error is not.
                    if (
                        isinstance(e, urllib.error.HTTPError)
                        and e.code not in _TRANSIENT
                    ):
                        if resumes:
                            _check_resume(url, e, size, total)
                        raise
                    if resumes == _RESUMES:
                        raise
                    resumes += 1
                    time.sleep(min(_WAIT * resumes, 60))
        except BrokenPipeError:
            pass  # norm stopped early: its exit status and stderr say why, below
        except BaseException:
            for tool in (norm, query):
                tool.kill()
                tool.wait()
            raise
        finally:
            with contextlib.suppress(BrokenPipeError):
                stdin.close()
        for cmd, tool, err in (
            (norm_cmd, norm, norm_err),
            (query_cmd, query, query_err),
        ):
            if code := tool.wait():
                err.seek(0)
                raise ToolFailed(code, cmd, None, err.read().decode(errors="replace"))
    return size, base64.b64encode(digest.digest()).decode()


def _check_resume(url: str, response: Any, start: int, total: int) -> None:
    """Refuse an answer to `Range: bytes=<start>-` other than a 206 from byte `start` of an
    object of `total` bytes: it would feed bytes twice, skip some, or splice another file."""
    found = response.headers.get("Content-Range") or ""
    match = re.fullmatch(r"bytes (\d+)-\d+/(\d+)", found)
    if (
        response.status != 206
        or match is None
        or match.groups() != (f"{start}", f"{total}")
    ):
        raise StreamMismatch(
            f"{url}: resuming at byte {start} of size {total}, the bucket answered"
            f" {response.status} {found!r}"
        )


def _sql_path(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _write_parquet(rows_tsv: Path, out: Path, tiers: TwoTier | None) -> int:
    """The query's rows as the table's columns, written to `out`; the rows written."""
    keep = "TRUE"
    with duckdb.connect() as con:
        if tiers is not None:
            con.execute("CREATE TABLE span (chrom VARCHAR, lo BIGINT, hi BIGINT)")
            con.executemany(
                "INSERT INTO span VALUES (?, ?, ?)",
                [list(s) for s in tiers.spans.values()],
            )
            keep = (
                f"r.faf >= {tiers.genome_floor!r} OR EXISTS (SELECT 1 FROM span s"
                " WHERE s.chrom = r.chrom AND r.pos BETWEEN s.lo AND s.hi)"
            )
        # The vid as make_vid writes it: gnomAD's contigs already carry the `chr` prefix.
        con.execute(
            "COPY (SELECT 'GRCh38:' || r.chrom || ':' || r.pos || ':' || upper(r.ref)"
            " || ':' || upper(r.alt) AS vid, r.faf AS faf95_max, r.gen_anc"
            f" FROM read_csv({_sql_path(rows_tsv)}, delim = '\t', header = false,"
            " quote = '', escape = '', nullstr = '.', auto_detect = false, columns ="
            " {'chrom': 'VARCHAR', 'pos': 'BIGINT', 'ref': 'VARCHAR', 'alt': 'VARCHAR',"
            " 'faf': 'DOUBLE', 'gen_anc': 'VARCHAR'}) r"
            f" WHERE {keep}) TO {_sql_path(out)} (FORMAT parquet)"
        )
        found = con.execute(
            f"SELECT count(*) FROM read_parquet({_sql_path(out)})"
        ).fetchone()
    return int(found[0]) if found is not None else 0


def _record(
    folder: Path, chrom: str, entry: dict[str, Any], policy: dict[str, Any]
) -> None:
    """Write `chrom`'s entry into the manifest through a temporary file and a rename. The
    table's identity is the manifest's sha256 alone, so every file written or replaced
    rewrites it; a manifest of another floor or other tiers starts over with `chrom`. An
    exclusive lock held from the read to the rename keeps concurrent builds' entries."""
    with (folder / ".manifest.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        listed = {**_listed(folder, policy), chrom: entry}
        data = {
            **policy,
            "chromosomes": {c: listed[c] for c in GNOMAD_CHROMS if c in listed},
        }
        path = folder / MANIFEST
        tmp = path.with_name(f".{MANIFEST}.tmp")
        tmp.write_text(json.dumps(data, indent=1) + "\n")
        tmp.replace(path)


def _build_faf(args: argparse.Namespace) -> int:
    paths = load_paths()
    grading = load_grading()
    floor = faf_floor(grading)
    tiers = two_tier(grading) if FAF_TWO_TIER else None
    folder = faf_folder(paths.lab)
    if args.check:
        for chrom in stale_chromosomes(folder, floor=floor, two_tier=tiers):
            print(chrom)
        return 0
    fasta = unmasked_fasta(paths.lab)
    if not fasta.is_file() or not fasta.with_name(f"{fasta.name}.fai").is_file():
        print(
            f"genokeeper build-faf: missing input: no GRCh38 FASTA with its .fai at {fasta};"
            " run scripts/bootstrap.sh",
            file=sys.stderr,
        )
        return 2
    if shutil.which("bcftools") is None:
        print(f"genokeeper build-faf: {ToolMissing('bcftools')}", file=sys.stderr)
        return 2
    doing = "the bucket's object metadata"
    try:
        metadata: dict[str, tuple[int, str]] = {}
        if args.chrom:
            chroms = list(dict.fromkeys(args.chrom))
        else:
            metadata = {c: object_metadata(c) for c in GNOMAD_CHROMS}
            chroms = stale_chromosomes(
                folder, floor=floor, two_tier=tiers, metadata=metadata
            )
        for chrom in chroms:
            doing = gnomad_url(chrom)
            print(f"genokeeper build-faf: streaming {doing}", file=sys.stderr)
            rows = build_chromosome(
                chrom,
                folder,
                floor=floor,
                fasta=fasta,
                expected=metadata.get(chrom),
                two_tier=tiers,
            )
            print(f"{chrom}\t{rows} variants")
    except (OSError, ValueError, StreamMismatch, ToolFailed) as e:
        print(f"genokeeper build-faf: {doing}: {e}", file=sys.stderr)
        return 1
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `build-faf` sub-command."""
    p = subparsers.add_parser(
        "build-faf",
        help="stream gnomAD v4.1 joint (0.88 TB) into the local FAF table (a few GB)",
    )
    p.add_argument(
        "--chrom",
        action="append",
        choices=GNOMAD_CHROMS,
        metavar="CHROM",
        help="build this chromosome (repeatable); default: each the manifest lacks or"
        " whose bucket md5 differs",
    )
    p.add_argument(
        "--check",
        action="store_true",
        help="list the chromosomes the local table lacks, one per line; no network",
    )
    p.set_defaults(handler=_build_faf)
