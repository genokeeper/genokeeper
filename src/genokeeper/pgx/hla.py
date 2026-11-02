"""HLA-A and HLA-B drug alleles where OptiType and T1K agree (track design, slice 1b,
"Outside calls").

The reads: the MHC region (`MHC_REGION`, GRCh38) and the unplaced unmapped reads, as pairs,
from the CRAM by samtools natively (`extract_hla_reads`). Two callers read them, each
offline in its pinned BioContainers image behind the core gate as `outside.LABEL`:
OptiType 1.5.0 (GitHub tag v1.5.0, commit 407872e0; RazerS3 and the GLPK solver in its
image) and T1K 1.0.10 (`--preset hla-wgs`) with its index built from IPD-IMGT/HLA
`IMGT_HLA_VERSION` (`build_t1k_index`, run by the bootstrap). Alleles are compared at two
fields (`two_fields`); T1K's alleles of quality 0 or less are dropped.

For each PGx allele (`PGX_HLA`) a caller's call is positive (one of its alleles is it),
negative (none can be), or unknown (an allele dropped, or several equivalent alleles one
of which is it). The two must agree, positive or negative, at every one of the gene's PGx
alleles; otherwise the gene is unresolved, the alleles in dispute named. One outside line
per gene: PharmCAT finds each allele's status from a diplotype, not from a phenotype naming
several (its run found no recommendation for `*15:02 negative; *57:01 positive; ...`), so a
resolved gene is given OptiType's diplotype, whose status at every PGx allele is the agreed
one; the two diplotypes may still differ at an allele no recommendation reads.
"""

import csv
import re
import shutil
import sys
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Literal

from genokeeper.containers import PLATFORM, run_container
from genokeeper.pgx.outside import LABEL, CallerFailed, OutsideCall
from genokeeper.tools import gated, run_tool

MHC_REGION = ("chr6", 28510120, 33480577)
PGX_HLA = {"HLA-A": ("*31:01",), "HLA-B": ("*15:02", "*57:01", "*58:01")}
OPTITYPE, T1K = "OptiType", "T1K"
# The IPD-IMGT/HLA release T1K's index is built from (bootstrap step_t1k).
IMGT_HLA_VERSION = "3.66.0"
T1K_INDEX = "hla_dna_seq.fa"
OPTITYPE_PREFIX = "optitype"
T1K_PREFIX = "t1k"
FASTQ = ("hla_1.fq", "hla_2.fq")
# The genes OptiType types, and the only ones read from T1K.
CLASS_I = ("HLA-A", "HLA-B", "HLA-C")
# An allele name's first two fields, and an expression suffix ending the second (`N`).
_TWO_FIELDS = re.compile(r"\*(\d+:\d+)([A-Z]?)(?=$|:)")

# A caller's call of one gene: two alleles at two fields, an allele T1K dropped None, and
# equivalent alleles joined by `,`.
type HlaCall = tuple[str | None, str | None]


def t1k_index(lab: Path) -> Path:
    """Where the bootstrap files T1K's HLA index (`step_t1k`)."""
    return lab / "t1k" / f"IPD-IMGT-HLA-{IMGT_HLA_VERSION}" / T1K_INDEX


def two_fields(name: str) -> str:
    """`HLA-B*57:01:01:02` (or `B*57:01:01`) as `*57:01`; `*01:11N` keeps its `N`."""
    match = _TWO_FIELDS.search(name)
    if match is None:
        raise ValueError(f"not an HLA allele name: {name!r}")
    return f"*{match[1]}{match[2]}"


def extract_hla_reads(cram: Path, fasta: Path, out_dir: Path) -> tuple[Path, Path]:
    """The pairs of `cram` (indexed, written against `fasta`) with a read in the MHC
    region or both unmapped, without secondary, supplementary or duplicate records, as
    `out_dir/hla_1.fq` and `hla_2.fq`. One samtools pipeline behind one gate; the files
    are renamed into place only once it has finished."""
    out_dir.mkdir(parents=True, exist_ok=True)
    chrom, start, end = MHC_REGION
    finals = [out_dir / name for name in FASTQ]
    parts = [path.with_name(path.name + ".part") for path in finals]
    script = (
        "set -euo pipefail; "
        'samtools view -u -h -M -F 0xD00 -T "$1" "$2" "$3" "*" '
        '| samtools collate -u -O - "$4/collate" '
        '| samtools fastq -n -0 /dev/null -s /dev/null -1 "$5" -2 "$6" -'
    )
    work = out_dir / ".collate"
    work.mkdir(exist_ok=True)
    try:
        run_tool(
            gated(
                LABEL,
                [
                    "bash", "-c", script, "extract-hla",
                    str(fasta), str(cram), f"{chrom}:{start}-{end}", str(work),
                    str(parts[0]), str(parts[1]),
                ],
            )
        )  # fmt: skip
        for part, final in zip(parts, finals, strict=True):
            part.replace(final)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        for part in parts:
            part.unlink(missing_ok=True)
    return finals[0], finals[1]


def read_optitype(path: Path) -> dict[str, HlaCall]:
    """OptiType's result (`<prefix>_result.tsv`): each gene's two alleles at two fields."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise CallerFailed(f"{path}: no OptiType result")
    row = rows[0]
    calls: dict[str, HlaCall] = {}
    for letter in ("A", "B", "C"):
        alleles = [row.get(f"{letter}{i}") or "" for i in (1, 2)]
        if all(alleles):
            calls[f"HLA-{letter}"] = (two_fields(alleles[0]), two_fields(alleles[1]))
    return calls


def _t1k_allele(names: str, quality: str) -> str | None:
    if names == "." or float(quality) <= 0:
        return None
    return ",".join(dict.fromkeys(two_fields(n) for n in names.split(",")))


def read_t1k(
    path: Path, *, genes: Collection[str] | None = CLASS_I
) -> dict[str, HlaCall]:
    """T1K's genotypes (`<prefix>_genotype.tsv`): each gene's alleles at two fields, its
    equivalent alleles that share two fields made one, one of quality 0 or less dropped,
    one distinct allele a homozygote; a gene with none kept is left out. By default only the
    class I genes OptiType calls are read, for PharmCAT; `genes=None` reads every `HLA-`
    gene T1K names at two fields (slice 3f's class II), skipping one named at one field
    (`MICB*034`, as its run showed)."""
    calls: dict[str, HlaCall] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 8 or not fields[1].isdigit() or fields[1] == "0":
            continue
        gene = fields[0]
        if not (gene.startswith("HLA-") if genes is None else gene in genes):
            continue
        try:
            first = _t1k_allele(fields[2], fields[4])
            second = first if fields[1] == "1" else _t1k_allele(fields[5], fields[7])
        except ValueError:
            if genes is not None:
                raise
            continue
        if first is not None or second is not None:
            calls[gene] = (first, second)
    return calls


def run_optitype(fq1: Path, fq2: Path, out_dir: Path) -> dict[str, HlaCall]:
    """OptiType on the paired reads (both in one directory), offline with GLPK; its result
    and coverage plot stay in `out_dir`. Raises `CallerFailed` without a result."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "optitype",
        [
            "optitype", "run", "-i", f"/hla-in/{fq1.name}", "-i", f"/hla-in/{fq2.name}",
            "--dna", "-o", "/optitype-out", "-p", OPTITYPE_PREFIX, "--solver", "glpk",
        ],
        mounts=[(fq1.parent, "/hla-in", True), (out_dir, "/optitype-out", False)],
        platform=PLATFORM,
        label=LABEL,
    )  # fmt: skip
    result = out_dir / f"{OPTITYPE_PREFIX}_result.tsv"
    if not result.is_file():
        raise CallerFailed(f"OptiType wrote no {result.name} in {out_dir}")
    return read_optitype(result)


def run_t1k(fq1: Path, fq2: Path, out_dir: Path, *, index: Path) -> dict[str, HlaCall]:
    """T1K on the paired reads (both in one directory) with the HLA `index`
    (`t1k_index`), as WGS, offline; its files stay in `out_dir`. Raises `CallerFailed`
    without a genotype table."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "t1k",
        [
            "run-t1k", "-1", f"/hla-in/{fq1.name}", "-2", f"/hla-in/{fq2.name}",
            "--preset", "hla-wgs", "-f", f"/t1k-index/{index.name}",
            "-o", T1K_PREFIX, "--od", "/t1k-out",
        ],
        mounts=[
            (fq1.parent, "/hla-in", True),
            (index.parent, "/t1k-index", True),
            (out_dir, "/t1k-out", False),
        ],
        platform=PLATFORM,
        label=LABEL,
    )  # fmt: skip
    result = out_dir / f"{T1K_PREFIX}_genotype.tsv"
    if not result.is_file():
        raise CallerFailed(f"T1K wrote no {result.name} in {out_dir}")
    return read_t1k(result)


def _status(call: HlaCall | None, allele: str) -> str:
    if call is None:
        return "unknown"
    if any(a is not None and a.split(",") == [allele] for a in call):
        return "positive"
    if any(a is None or allele in a.split(",") for a in call):
        return "unknown"
    return "negative"


def status(
    call: HlaCall | None, allele: str
) -> Literal["positive", "negative", "unknown"]:
    """A caller's call at `allele`: a two-field allele (`*57:01`) as `_status` reads it, or a
    one-field group (`*27`), whose members are the names starting `*27:`. Positive when one
    allele's every equivalent name is a member, unknown when only some are or an allele was
    dropped (None), negative otherwise."""
    group = ":" not in allele

    def member(name: str) -> bool:
        return name.startswith(f"{allele}:") if group else name == allele

    if call is None:
        return "unknown"
    names = [None if a is None else a.split(",") for a in call]
    if any(n is not None and all(member(x) for x in n) for n in names):
        return "positive"
    if any(n is None or any(member(x) for x in n) for n in names):
        return "unknown"
    return "negative"


def shown(call: HlaCall | None) -> str:
    return "no call" if call is None else "/".join(a or "?" for a in call)


def hla_outside(
    opt: Mapping[str, HlaCall], t1k: Mapping[str, HlaCall]
) -> list[OutsideCall]:
    """One outside call per gene of `PGX_HLA` from OptiType's and T1K's calls (see the
    module docstring)."""
    found: list[OutsideCall] = []
    for gene, alleles in PGX_HLA.items():
        o, t = opt.get(gene), t1k.get(gene)
        callers = {OPTITYPE: shown(o), T1K: shown(t)}
        statuses: list[str] = []
        disputed: list[str] = []
        for allele in alleles:
            so, st = _status(o, allele), _status(t, allele)
            agreed = so == st and so != "unknown"
            statuses.append(f"{allele} {so if agreed else 'unresolved'}")
            if not agreed:
                disputed.append(f"{allele} ({OPTITYPE} {so}, {T1K} {st})")
        silent = [name for name, call in ((OPTITYPE, o), (T1K, t)) if call is None]
        reason = None
        if silent:
            reason = f"{' and '.join(silent)} made no call"
        elif disputed:
            reason = f"{OPTITYPE} and {T1K} differ on {', '.join(disputed)}"
        resolved = reason is None and o is not None
        found.append(
            OutsideCall(
                gene=gene,
                diplotype="/".join(a or "" for a in o) if resolved and o else None,
                phenotype="; ".join(statuses),
                activity_score=None,
                callers=callers,
                agreed=resolved,
                unresolved_reason=reason,
            )
        )
    return found


def build_t1k_index(dat: Path, out_dir: Path) -> Path:
    """T1K's HLA index built from IPD-IMGT/HLA's `hla.dat` into `out_dir`, offline in T1K's
    pinned image; the index file's path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "t1k",
        [
            "t1k-build.pl",
            "-d",
            f"/t1k-dat/{dat.name}",
            "-o",
            "/t1k-index",
            "--prefix",
            "hla",
        ],
        mounts=[(dat.parent, "/t1k-dat", True), (out_dir, "/t1k-index", False)],
        platform=PLATFORM,
        label=LABEL,
    )
    index = out_dir / T1K_INDEX
    if not index.is_file():
        raise CallerFailed(f"t1k-build.pl wrote no {T1K_INDEX} in {out_dir}")
    return index


def main(argv: Sequence[str] | None = None) -> int:
    """`python -m genokeeper.pgx.hla build-index DAT OUT_DIR` (the bootstrap's step_t1k)."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3 or args[0] != "build-index":
        print(
            "usage: python -m genokeeper.pgx.hla build-index DAT OUT_DIR",
            file=sys.stderr,
        )
        return 2
    print(build_t1k_index(Path(args[1]), Path(args[2])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
