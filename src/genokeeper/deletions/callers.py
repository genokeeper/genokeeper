"""DELLY 2.7.0 and Dysgu 1.9.0 on the CRAM, and their deletions read as `Deletion`s (slice 3b).

Both run offline in their images (`containers.IMAGES`, built from docker/) behind the core
gate, on the CRAM, its index and the FASTA (with its `.fai`) mounted read-only and one output
directory mounted read-write. DELLY parallelises over samples only, so one sample runs on one
core: `run_delly` runs one container per chromosome it is given (`exons.PRIMARY`, or without
chrY for a female sample: DELLY 2.7.0 exits 1 on a shard with too few reads to estimate its
library), each behind its own `dna:delly` gate (`delly_shards`, which slice 3d's `delly cnv`
shards through too), and concatenates their BCFs natively.
`run_dysgu` shards the same way, one chromosome after another, to bound its memory. bcftools
runs natively (`run_tool`), as 3a's wrappers run samtools.

Coordinates are `Deletion`'s: `pos` is the base before the deletion, `end` its last deleted
base. A sequence-resolved record's end comes from its REF (`POS + len(REF) - 1`), whatever its
`INFO/END` says: Dysgu 1.9.0 sets `END` (and `SVLEN`) one base past its own REF, on 355 of 355
of HG002's chr20 deletions (plan Task 9, measured 2026-10-28); a symbolic record's is its
`INFO/END`.

Of another type (slice 3d): an insertion is read at its point (`end = POS + 1`) with its
inserted length, `len(ALT) - len(REF)` for a sequence-resolved record and `|INFO/SVLEN|` for a
symbolic one; a duplication's end is its `INFO/END`, sequence-resolved or not (Dysgu writes
both). A Dysgu insertion whose length is inferred, not assembled (`INFO/LPREC=0`), and an
INSurVeyor one with `INFO/INCOMPLETE_ASSEMBLY` have only a lower bound (`lower_bound`); such
an INSurVeyor record carries no `SVLEN`, and its assembled part (`INFO/SVINSSEQ`) is the
bound.
"""

import gzip
from collections.abc import Collection, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genokeeper.containers import run_container
from genokeeper.deletions.screen import Deletion
from genokeeper.tools import run_tool

DELLY_LABEL = "dna:delly"
DYSGU_LABEL = "dna:dysgu"
# The GRCh38 exclude list inside DELLY's image (docker/delly.Dockerfile, pinned by sha256).
DELLY_TEMPLATE = "/opt/delly/human.hg38.excl.tsv"
# Per type, a record's length as `read_deletions` reads it. A deletion: `INFO/END - POS` for
# a symbolic record, `len(REF) - 1` for a sequence-resolved one (not Dysgu's `SVLEN`, one
# more). An insertion: `|INFO/SVLEN|` symbolic, `len(ALT) - len(REF)` sequence-resolved. A
# duplication: `INFO/END - POS`.
_LENGTH = {
    "DEL": '((ALT~"^<" && INFO/END-POS>={n}) || (ALT!~"^<" && strlen(REF)-1>={n}))',
    "INS": '((ALT~"^<" && abs(INFO/SVLEN)>={n})'
    ' || (ALT!~"^<" && strlen(ALT)-strlen(REF)>={n}))',
    "DUP": "INFO/END-POS>={n}",
}
# INSurVeyor writes no `SVLEN` on an incomplete assembly (11 of HG002's 92 chr20 calls, Task
# 6): its assembled part, `SVINSSEQ`, is the length's lower bound (`_event`).
_INSURVEYOR_PART = '(ALT~"^<" && INFO/SVLEN="." && strlen(INFO/SVINSSEQ)>={n})'
# Per caller and type, what `deletions_only` keeps besides the length. DELLY: FILTER PASS
# and the sample's FT PASS (its LowQual marks PE and SR support under 3 reads or MAPQ 20,
# and GQ under 15; modvcf.h at v2.7.0), IMPRECISE calls kept. Dysgu: FILTER PASS (PROB >=
# 0.45, its paired-end default; no DNBSEQ-specific value is published). INSurVeyor: FILTER
# PASS (its out.pass.vcf.gz holds nothing else).
_KEEP = {
    ("delly", "DEL"): 'INFO/SVTYPE="DEL" && FILTER="PASS" && FMT/FT="PASS" && GT="alt"',
    ("dysgu", "DEL"): 'INFO/SVTYPE="DEL" && FILTER="PASS" && GT="alt"',
    ("dysgu", "INS"): 'INFO/SVTYPE="INS" && FILTER="PASS" && GT="alt"',
    ("dysgu", "DUP"): 'INFO/SVTYPE="DUP" && FILTER="PASS" && GT="alt"',
    ("insurveyor", "INS"): 'INFO/SVTYPE="INS" && FILTER="PASS" && GT="alt"',
}
_CALLERS = frozenset(caller for caller, _ in _KEEP)


def _inputs(cram: Path, fasta: Path, out_dir: Path) -> list[tuple[Path, str, bool]]:
    return [
        (cram.parent, "/cram", True),
        (fasta.parent, "/ref", True),
        (out_dir, "/out", False),
    ]


def _lengths(fasta: Path) -> dict[str, str]:
    """The FASTA's contig names, in its `.fai`'s order, and their lengths."""
    lines = Path(f"{fasta}.fai").read_text().splitlines()
    return dict(line.split("\t")[:2] for line in lines if line)


def _contigs(fasta: Path) -> list[str]:
    """The FASTA's contig names, in its `.fai`'s order."""
    return list(_lengths(fasta))


def _exclude(
    template: str,
    contigs: Sequence[str],
    own: str,
    path: Path,
    extra: Sequence[str] = (),
) -> None:
    """`template`'s lines, the `extra` ones, and a whole-contig line (a bare name, as
    DELLY's `-x` reads it) for every contig but `own`."""
    lines = [template.rstrip("\n"), *extra, *(c for c in contigs if c != own)]
    path.write_text("\n".join(lines) + "\n")


def delly_shards(
    cram: Path,
    fasta: Path,
    out_dir: Path,
    *,
    threads: int,
    shards: Collection[str],
    folder: str,
    command: Sequence[str],
    label: str,
    extra: Sequence[str] = (),
) -> list[Path]:
    """DELLY's `command` (its subcommand and options before `-x`) on `cram`, one container
    per chromosome of the FASTA in `shards` (in its `.fai`'s order), at most `threads` at a
    time, each behind its own `label` gate in its own `out_dir/<folder>/<chrom>/` with its
    own exclude list (DELLY's template, the `extra` lines, every other contig); their BCFs
    in that order. A shard that fails fails the step: the error is raised once the running
    shards finish."""
    contigs = _contigs(fasta)
    folders = [out_dir / folder / c for c in contigs if c in shards]
    run_container(
        "delly",
        ["cp", DELLY_TEMPLATE, "/out/exclude.template.tsv"],
        mounts=[(out_dir, "/out", False)],
    )
    template = (out_dir / "exclude.template.tsv").read_text()

    def shard(path: Path) -> Path:
        path.mkdir(parents=True, exist_ok=True)
        _exclude(template, contigs, path.name, path / "exclude.tsv", extra)
        run_container(
            "delly",
            [
                *command, "-x", "/out/exclude.tsv", "-o", "/out/delly.bcf",
                f"/cram/{cram.name}",
            ],
            mounts=_inputs(cram, fasta, path),
            label=label,
        )  # fmt: skip
        return path / "delly.bcf"

    pool = ThreadPoolExecutor(max_workers=threads)
    try:
        return list(pool.map(shard, folders))
    finally:
        pool.shutdown(cancel_futures=True)


def run_delly(
    cram: Path, fasta: Path, out_dir: Path, *, threads: int, shards: Collection[str]
) -> Path:
    """DELLY's deletions on `cram` (GRCh38), sharded per chromosome (`delly_shards`, in
    `out_dir/delly-shards/`); their BCFs concatenated in the `.fai`'s order into
    `out_dir/delly.vcf.gz`, bgzipped and indexed. A shard that fails fails the step, and no
    file is written."""
    bcfs = delly_shards(
        cram,
        fasta,
        out_dir,
        threads=threads,
        shards=shards,
        folder="delly-shards",
        command=["delly", "sr", "-t", "DEL", "-g", f"/ref/{fasta.name}"],
        label=DELLY_LABEL,
    )
    out = out_dir / "delly.vcf.gz"
    run_tool(
        ["bcftools", "concat", "--no-version", "-Oz", "-o", str(out), *map(str, bcfs)]
    )
    run_tool(["bcftools", "index", "-f", "-t", str(out)])
    return out


def run_dysgu(
    cram: Path, fasta: Path, out_dir: Path, *, threads: int, shards: Collection[str]
) -> Path:
    """Dysgu's calls on `cram` (GRCh38; `--mode pe`, its default, and no `dysgu filter` on one
    sample), one container per chromosome of the FASTA in `shards` (in its `.fai`'s order),
    one after another on `threads`, each limited to its chromosome (`--search`) in its
    own `out_dir/dysgu-shards/<chrom>/` and kept to that chromosome's records (a
    translocation's other end may name another, and `bcftools concat` refuses a chromosome
    in two blocks; no deletion is lost); concatenated in that order, sorted, bgzipped and
    indexed natively as `out_dir/dysgu.vcf.gz`, every record unfiltered, kept in the step's
    folder: slice 3d's `dysgu_ins_dup` step reads its insertions and duplications there
    (`pipeline.dysgu_ins_dup_step`). One whole-genome run outgrew Docker's
    memory on a 30x library (21.8 GiB on HG002's); a chromosome took 1.0 to 1.9 GiB
    on HG002's (plan Task 15, measured 2026-10-29). A shard that fails fails the step."""
    bcfs = [
        _dysgu_shard(cram, fasta, out_dir / "dysgu-shards" / c, length, threads)
        for c, length in _lengths(fasta).items()
        if c in shards
    ]
    plain = out_dir / "dysgu.vcf"
    run_tool(
        [
            "bcftools",
            "concat",
            "--no-version",
            "-Ov",
            "-o",
            str(plain),
            *map(str, bcfs),
        ]
    )
    return sort_index(plain, out_dir / "dysgu.vcf.gz")


def _dysgu_shard(
    cram: Path, fasta: Path, shard: Path, length: str, threads: int
) -> Path:
    """Dysgu's calls on the one chromosome `shard.name` (`length` bases long), its work
    directory `shard/dysgu-work`, those on that chromosome as `shard/dysgu.bcf`."""
    shard.mkdir(parents=True, exist_ok=True)
    (shard / "search.bed").write_text(f"{shard.name}\t0\t{length}\n")
    run_container(
        "dysgu",
        [
            "dysgu", "run", "-p", str(threads), "-x", "--search", "/out/search.bed",
            "-o", "/out/dysgu.vcf", f"/ref/{fasta.name}", "/out/dysgu-work",
            f"/cram/{cram.name}",
        ],
        mounts=_inputs(cram, fasta, shard),
        label=DYSGU_LABEL,
    )  # fmt: skip
    own = shard / "dysgu.bcf"
    run_tool(
        [
            "bcftools", "view", "--no-version", "-t", shard.name, "-Ob", "-o", str(own),
            str(shard / "dysgu.vcf"),
        ]
    )  # fmt: skip
    return own


def sort_index(plain: Path, out: Path) -> Path:
    """`plain` sorted (its temporary files under `out`'s folder, not the system's), bgzipped
    and indexed natively as `out`; `plain` removed."""
    temp = f"{out.parent}/.{out.name}.sort."
    run_tool(["bcftools", "sort", "-T", temp, "-Oz", "-o", str(out), str(plain)])
    run_tool(["bcftools", "index", "-f", "-t", str(out)])
    plain.unlink()
    return out


def filter_expression(caller: str, min_length: int, svtype: str = "DEL") -> str:
    """The bcftools expression `deletions_only` keeps `caller`'s records of `svtype` by
    (`_KEEP` and the length)."""
    if (caller, svtype) not in _KEEP:
        raise ValueError(f"no {svtype} filter for caller {caller!r}")
    length = _LENGTH[svtype]
    if caller == "insurveyor":
        length = f"({length} || {_INSURVEYOR_PART})"
    return f"{_KEEP[caller, svtype]} && {length.format(n=min_length)}"


def caller_min_length(caller: str, min_length: int) -> int:
    """The length `caller`'s deletions are kept from under the screen's `min_length`:
    Dysgu's at it, DELLY's at half of it, so a DELLY call a little short can still confirm
    a Dysgu one (`truvari.match_callers`' `--sizefilt`); any other caller's at it."""
    if caller not in _CALLERS:
        raise ValueError(f"no deletion filter for caller {caller!r}")
    return min_length // 2 if caller == "delly" else min_length


def caller_filter(caller: str, min_length: int) -> str:
    """`filter_expression` at `caller_min_length`: what the run's caller step keys on and
    the deletion calibration records among its versions."""
    return filter_expression(caller, caller_min_length(caller, min_length))


def deletions_only(
    vcf: Path, out: Path, *, caller: str, min_length: int, svtype: str = "DEL"
) -> Path:
    """`caller`'s records of `svtype` in `vcf` of at least `min_length` bp (as
    `read_deletions` measures them) with a non-reference GT and its PASS filters (`_KEEP`),
    as the sorted, bgzipped and indexed VCF `out`."""
    kept = out.with_name(f"{out.name}.kept.bcf")
    expression = filter_expression(caller, min_length, svtype)
    run_tool(
        [
            "bcftools", "view", "--no-version", "-i", expression,
            "-Ob", "-o", str(kept), str(vcf),
        ]
    )  # fmt: skip
    return sort_index(kept, out)


def _info(text: str) -> dict[str, str]:
    """A record's INFO as name -> value; a flag's value is empty."""
    pairs = (item.partition("=") for item in text.split(";"))
    return {name: value for name, _, value in pairs}


def depth_ratio(values: Mapping[str, str]) -> float | None:
    """DELLY's read count inside the deletion over its two flanks' together (`RC`, `RCL`,
    `RCR` at v2.7.0); None when a field is missing or a flank count is 0."""
    try:
        inside, left, right = (int(values[k]) for k in ("RC", "RCL", "RCR"))
    except KeyError, ValueError:
        return None
    return None if left == 0 or right == 0 else inside / (left + right)


def _of_type(svtype: str, info: Mapping[str, str], ref: str, alt: str) -> bool:
    """Whether a record is of `svtype`: its `SVTYPE`; a sequence-resolved deletion also by a
    REF longer than its ALT, an insertion without `SVTYPE` by an ALT longer than its REF."""
    declared = info.get("SVTYPE")
    resolved = not alt.startswith("<")
    if svtype == "DEL":
        return declared == "DEL" or (resolved and len(alt) < len(ref))
    if svtype == "INS":
        return declared == "INS" or (
            declared is None and resolved and len(alt) > len(ref)
        )
    return declared == svtype


def _event(
    svtype: str, start: int, ref: str, alt: str, info: Mapping[str, str]
) -> tuple[int, int | None, bool]:
    """A record's end, `svlen` and `lower_bound` as the module docstring reads them."""
    symbolic = alt.startswith("<")
    if svtype == "DEL":
        return (int(info["END"]) if symbolic else start + len(ref) - 1), None, False
    if svtype == "DUP":
        return int(info["END"]), None, False
    bound = info.get("LPREC") == "0" or "INCOMPLETE_ASSEMBLY" in info
    if not symbolic:
        return start + 1, len(alt) - len(ref), bound
    if info.get("SVLEN", ".") not in (".", ""):
        return start + 1, abs(int(info["SVLEN"])), bound
    part = info.get("SVINSSEQ", ".")
    return start + 1, (None if part in (".", "") else len(part)), True


def vcf_records(
    vcf: Path,
) -> Iterator[tuple[list[str], dict[str, str], dict[str, str]]]:
    """Each record of `vcf` (gzip when it ends in `.gz`): its first eight columns (ALT cut
    to its first allele), its INFO and its first sample's FORMAT values."""
    opener = gzip.open if vcf.suffix == ".gz" else open
    with opener(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys, sample = fields[8:10]
            values = dict(zip(keys.split(":"), sample.split(":"), strict=False))
            cols = [*fields[:4], fields[4].split(",")[0], *fields[5:8]]
            yield cols, _info(fields[7]), values


def read_records(
    vcf: Path, *, caller: str, svtype: str = "DEL", pass_only: bool = False
) -> list[tuple[Deletion, dict[str, str]]]:
    """Every record of `svtype` in `vcf` (`vcf_records`; `_of_type`), with its INFO,
    unfiltered unless `pass_only` (FILTER PASS only): POS, its end, REF's first base, the
    first sample's GT, `callers = {caller: "<ID> QUAL <QUAL>"}`, for DELLY its
    `depth_ratio`, and an insertion's `svlen` and `lower_bound`."""
    found: list[tuple[Deletion, dict[str, str]]] = []
    for cols, info, values in vcf_records(vcf):
        chrom, pos, vid, ref, alt, qual, filt, _ = cols
        if not _of_type(svtype, info, ref, alt) or (pass_only and filt != "PASS"):
            continue
        start = int(pos)
        end, svlen, lower_bound = _event(svtype, start, ref, alt, info)
        deletion = Deletion(
            chrom=chrom,
            pos=start,
            end=end,
            ref=ref[0],
            gt=values["GT"],
            callers={caller: f"{vid} QUAL {qual}"},
            depth_ratio=depth_ratio(values) if caller == "delly" else None,
            svtype=svtype,
            svlen=svlen,
            lower_bound=lower_bound,
        )
        found.append((deletion, info))
    return found


def read_deletions(vcf: Path, *, caller: str, svtype: str = "DEL") -> list[Deletion]:
    """Every record of `svtype` in `vcf`, as `read_records` reads it."""
    return [d for d, _ in read_records(vcf, caller=caller, svtype=svtype)]
