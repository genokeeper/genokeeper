"""Runs of homozygosity (slice 3f, "ROH").

`bcftools roh -G30 --AF-file` (in the `truvari` image's bcftools 1.25, offline, behind the
core gate as `dna:roh`; `run_roh`) reads genotypes at the AF file's sites: the
superpopulation pgsc_calc names, one file per superpopulation built from the HGDP+1kGP panel
(`typing_files.build_roh_af`). It skips a row whose ALT is `.` or `<*>`, so every row it
reads must carry the AF file's REF and ALT with a reference call written `0/0` (ruling R7):

- `sites_vcf`: our genotypes, from the gVCF through the call layer (`gvcf.calls_at`) at each
  AF allele and written by pgsc_calc's writer (`pgsc.write_target_vcf`), each graded by the
  caller's `level` (a low or ungraded call, or an unknown one, is `./.`);
- `truth_sites_vcf`: a truth's genotypes at the AF alleles inside its benchmark, an absent
  site `0/0`, one where the truth calls another allele `./.` (as our calls write it).

Both take the same sites (`af_alleles`: autosomes, positions with one ALT in the AF file), so
the run and the calibration feed `run_roh` one kind of input. `inside` cuts our sites to a
benchmark, and `spike` makes an artificial isodisomy the calibration must find. `read_runs` reads the runs (`RG` lines) and
`summarise` gives the total length and count of long runs, the longest run and the runs long
enough to suggest isodisomy. Heterodisomy leaves no run, and nothing here speaks of the
parents.
"""

import gzip
import itertools
from collections import defaultdict
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.config import RohThresholds
from genokeeper.containers import run_container
from genokeeper.genotype.gvcf import Allele, Called, calls_at
from genokeeper.genotype.reliability import StrataIndex
from genokeeper.prs.pgsc import write_target_vcf
from genokeeper.store.db import Level

# The core gate's label for `bcftools roh`.
LABEL = "dna:roh"
AUTOSOMES = frozenset(f"chr{n}" for n in range(1, 23))
# Alleles per `calls_at` call: a chromosome's AF sites run to millions.
BATCH = 200_000
_HEADER = (
    "##fileformat=VCFv4.2\n"
    # `genomelab`, the name before the rename: this text reaches a step key.
    "##source=genomelab (truth genotypes at the AF file's alleles)\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


@dataclass(frozen=True)
class Run:
    """One run of homozygosity: its chromosome, first and last base, length (bp) and
    bcftools' quality."""

    chrom: str
    start: int
    end: int
    length: int
    quality: float


@dataclass(frozen=True)
class RohSummary:
    """The total length and the count of runs of `min_run` or more, the longest run, and
    each run of `isodisomy_min` or more."""

    total_long: int
    n_long: int
    longest: Run | None
    isodisomy: tuple[Run, ...]


def af_alleles(af: Path) -> Iterator[Allele]:
    """The AF file's alleles (`chr<N>`, position, `REF,ALT`, frequency; bgzipped, sorted),
    on chr1-22, leaving out every position the file lists with more than one ALT."""
    with gzip.open(af, "rt", encoding="utf-8") as handle:
        rows = (line.rstrip("\n").split("\t") for line in handle if line.strip())
        for (chrom, pos), group in itertools.groupby(rows, key=lambda r: (r[0], r[1])):
            site = list(group)
            if chrom not in AUTOSOMES or len(site) != 1:
                continue
            ref, alt = site[0][2].split(",")
            yield chrom, int(pos), ref.upper(), alt.upper()


def sites_vcf(
    gvcf: Path,
    af: Path,
    fasta: Path,
    out: Path,
    *,
    haploid: Callable[[str, int], bool],
    scratch: Path,
    level: Callable[[Allele, Called], Level | None],
) -> Path:
    """Our genotypes at the AF file's alleles (`af_alleles`) as `out` (gzipped VCF): each
    allele's call from DeepVariant's indexed `gvcf` on `fasta` (`calls_at`, with `haploid`
    and `scratch` as it takes them; read one chromosome at a time, in batches of `BATCH`),
    graded by `level`, written by `pgsc.write_target_vcf`.

    `out` holds the user's genotypes: the caller puts it under the run's work directory."""

    def graded() -> Iterator[tuple[Allele, Called | None, Level | None]]:
        for _, chromosome in itertools.groupby(af_alleles(af), key=lambda a: a[0]):
            for batch in itertools.batched(chromosome, BATCH):
                found = calls_at(
                    gvcf, batch, fasta=fasta, haploid=haploid, scratch=scratch
                )
                for allele in batch:
                    called = found[allele]
                    yield allele, called, level(allele, called)

    return write_target_vcf(graded(), out)


type _TruthRecord = tuple[str, tuple[str, ...], tuple[int | None, ...]]


def truth_gt(records: Sequence[_TruthRecord], ref: str, alt: str) -> str:
    """The truth's genotype for `ref>alt` from its records at the position (REF, ALTs,
    GT): `0/0` when none calls an allele there, the allele's copies when one record calls
    only it and the reference, else `./.` (another allele called, or a no-call)."""
    called = [r for r in records if any(a != 0 for a in r[2])]
    if not called:
        return "0/0"
    if len(called) > 1:
        return "./."
    rec_ref, alts, gt = called[0]
    if rec_ref != ref or alt not in alts or None in gt:
        return "./."
    index = alts.index(alt) + 1
    if any(a not in (0, index) for a in gt):
        return "./."
    copies = 2 if len(gt) == 1 else gt.count(index)
    return ("0/0", "0/1", "1/1")[min(copies, 2)]


def truth_sites_vcf(truth: Path, bench: Path, af: Path, out: Path) -> Path:
    """A truth's genotypes (`truth`: bgzipped, indexed VCF) at the AF file's alleles
    (`af_alleles`) inside the benchmark BED `bench`, with the AF file's REF and ALT, as `out`
    (gzipped VCF): an AF allele the truth does not call is `0/0`, one beside another called
    allele `./.`; sites outside the benchmark are left out."""
    intervals = StrataIndex.from_beds({"benchmark": bench})
    out.parent.mkdir(parents=True, exist_ok=True)
    with (
        pysam.VariantFile(str(truth)) as vcf,
        gzip.open(out, "wt", encoding="utf-8", compresslevel=1) as handle,
    ):
        handle.write(_HEADER)
        contigs = set(vcf.header.contigs)
        for chrom, alleles in itertools.groupby(af_alleles(af), key=lambda a: a[0]):
            calls: dict[int, list[_TruthRecord]] = defaultdict(list)
            for rec in vcf.fetch(chrom) if chrom in contigs else ():
                gt = tuple(rec.samples[0].get("GT") or (None,))
                alts = tuple(a.upper() for a in rec.alts or ())
                calls[rec.pos].append(((rec.ref or "").upper(), alts, gt))
            for _, pos, ref, alt in alleles:
                if intervals.flags_at(chrom, pos):
                    gt = truth_gt(calls.get(pos, ()), ref, alt)
                    handle.write(
                        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
                    )
    return out


_HETEROZYGOUS = frozenset({"0/1", "1/0", "0|1", "1|0"})


def spike(vcf: Path, region: tuple[str, int, int], out: Path) -> Path:
    """`vcf` (gzipped) as `out` with every heterozygous genotype inside `region` (chrom,
    first and last base, 1-based) rewritten homozygous ALT `1/1`: an artificial isodisomy
    the calibration must find. Every other row is copied unchanged."""
    chrom, start, end = region
    with (
        gzip.open(vcf, "rt", encoding="utf-8") as src,
        gzip.open(out, "wt", encoding="utf-8", compresslevel=1) as dst,
    ):
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if (
                not line.startswith("#")
                and fields[0] == chrom
                and start <= int(fields[1]) <= end
                and fields[9] in _HETEROZYGOUS
            ):
                fields[9] = "1/1"
                line = "\t".join(fields) + "\n"
            dst.write(line)
    return out


def inside(vcf: Path, bench: Path, out: Path) -> Path:
    """`vcf` (gzipped) as `out` with the header and only the rows inside the BED `bench`, as
    `truth_sites_vcf` leaves its rows out: `bcftools roh` refuses `-t` and `-T` beside
    `--AF-file`, and `-r` and `-R` need an indexed file."""
    intervals = StrataIndex.from_beds({"benchmark": bench})
    with (
        gzip.open(vcf, "rt", encoding="utf-8") as src,
        gzip.open(out, "wt", encoding="utf-8", compresslevel=1) as dst,
    ):
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            chrom, pos, _ = line.split("\t", 2)
            if intervals.flags_at(chrom, int(pos)):
                dst.write(line)
    return out


def run_roh(vcf: Path, af: Path, out: Path) -> Path:
    """`bcftools roh -G30 --AF-file af` on `vcf` (`sites_vcf`'s or `truth_sites_vcf`'s:
    autosomal biallelic SNVs) into `out` (`-O r`, its `RG` lines), in the `truvari` image,
    offline, behind the gate as `LABEL`."""
    out.parent.mkdir(parents=True, exist_ok=True)
    run_container(
        "truvari",
        [
            "bcftools", "roh", "-G30", "--AF-file", f"/af/{af.name}",
            "-O", "r", "-o", f"/out/{out.name}", f"/in/{vcf.name}",
        ],
        mounts=[
            (vcf.parent, "/in", True),
            (af.parent, "/af", True),
            (out.parent, "/out", False),
        ],
        label=LABEL,
    )  # fmt: skip
    return out


def read_runs(path: Path) -> list[Run]:
    """The runs of `bcftools roh -O r`'s output: its `RG` lines (sample, chromosome, start,
    end, length, markers, quality)."""
    runs: list[Run] = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] != "RG":
            continue
        _, _, chrom, start, end, length, _, quality = fields[:8]
        runs.append(Run(chrom, int(start), int(end), int(length), float(quality)))
    return runs


def summarise(runs: Sequence[Run], *, grading: RohThresholds) -> RohSummary:
    """The runs' summary: total length and count of runs of `min_run` or more, the longest
    run of all (None without runs), and the runs of `isodisomy_min` or more."""
    long = [r for r in runs if r.length >= grading.min_run]
    return RohSummary(
        total_long=sum(r.length for r in long),
        n_long=len(long),
        longest=max(runs, key=lambda r: r.length, default=None),
        isodisomy=tuple(r for r in runs if r.length >= grading.isodisomy_min),
    )
