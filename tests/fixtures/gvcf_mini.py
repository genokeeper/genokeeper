"""A miniature DeepVariant 1.10 gVCF on a small synthetic reference, built in a test's tmp_path.

Every record shape here was read off a real gVCF: DeepVariant 1.10.0
(`google/deepvariant@sha256:962e5a83…`, the run `genokeeper.call.deepvariant.run_deepvariant`
makes) with the T7 checkpoint, on public synthetic reads (a random 64 kb reference with
planted variants, aligned by `genokeeper.align.align`, called as a male with the PARs moved
into the mini contigs; Task 13 step 1 of the findings track plan, run on
2026-10-19). What it showed:

- a reference block: `REF <*>`, QUAL 0, FILTER `.`, `END=`, FORMAT `GT:GQ:MIN_DP:PL`, GT `0/0`,
  GQ at most 50 (one value per block, as `gvcf_gq_binsize=1` keeps them), no DP;
- a block over no reads is still `0/0`, with GQ 1, MIN_DP 0 and PL `0,0,0`;
- a block at a site the reads favour non-reference without a candidate: GT `./.`, GQ 0;
- every variant record ends its ALTs with `<*>`, AD and VAF carry a 0 for it, FORMAT
  `GT:GQ:DP:AD:VAF:PL`; FILTER is PASS, RefCall or NoCall;
- a RefCall's GT is `0/0` (also on a haploid contig);
- a NoCall's GT is `./.` with a GQ (a haploid site whose reads look heterozygous);
- all alleles at one position in one record: `A C,G,<*>` at `1/2`, or `0/1` when a third
  base's reads did not make a genotype; an SNV and a deletion as `TGG AGG,T,<*>`; two
  homopolymer deletions as `CAA C,CA,<*>` (the `CA` allele is `CA>C` once normalised);
- a deletion's REF span has no block: the next block starts after it; an SNV inside a called
  deletion (the other haplotype) is a record of its own;
- a male's non-PAR chrX: a hom call is written `1/1` and every block `0/0` (PL `0,9989,…`).

DeepVariant wrote no `*` allele (an SNV on the haplotype opposite a deletion came out as its
own record); `star_record` is hand-written from the VCF 4.2 specification (section 5.5) for
the call layer's defensive rule.
"""

import random
from pathlib import Path

import pysam

CONTIGS = {"chr1": 4000, "chrX": 4000}
SAMPLE = "mini"
# A homopolymer: C, eight A, G at chr1:1000-1009.
HOMOPOLYMER = ("chr1", 1000, "CAAAAAAAAG")
# A haploid region: a male's chrX beyond this position (the mini PAR is 1-1000).
PAR_END = 1000

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
##FILTER=<ID=NoCall,Description="Site has depth=0 resulting in no call.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
##FORMAT=<ID=MED_DP,Number=1,Type=Integer,Description="Median DP observed within the GVCF block rounded to the nearest integer.">
##DeepVariant_version=1.10.0
"""


def _sequences(seed: int = 13) -> dict[str, str]:
    rng = random.Random(seed)
    seqs = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in CONTIGS.items()}
    contig, pos, run = HOMOPOLYMER
    seqs[contig][pos - 1 : pos - 1 + len(run)] = list(run)
    return {c: "".join(s) for c, s in seqs.items()}


SEQUENCES = _sequences()


def ref(contig: str, pos: int, length: int = 1) -> str:
    """The reference bases at 1-based `pos`."""
    return SEQUENCES[contig][pos - 1 : pos - 1 + length]


def other_base(base: str, k: int = 0) -> str:
    """The `k`th base (0-2) that is not `base`."""
    return [b for b in "ACGT" if b != base][k]


def write_fasta(path: Path) -> Path:
    """The mini reference, with its .fai."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for name, seq in SEQUENCES.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def block(
    contig: str,
    start: int,
    end: int,
    *,
    gq: int = 50,
    min_dp: int = 17,
    gt: str = "0/0",
) -> str:
    """A reference block over `start..end` (1-based, inclusive), as DeepVariant writes it."""
    pl = "0,0,0" if min_dp == 0 else f"0,{gq + 1},{10 * gq + 9}"
    return (
        f"{contig}\t{start}\t.\t{ref(contig, start)}\t<*>\t0\t.\tEND={end}"
        f"\tGT:GQ:MIN_DP:PL\t{gt}:{gq}:{min_dp}:{pl}"
    )


def variant(
    contig: str,
    pos: int,
    ref_len: int,
    alts: list[str],
    *,
    gt: str,
    gq: int,
    ad: list[int],
    filt: str = "PASS",
    qual: float = 40.0,
) -> str:
    """A variant record: REF is the reference's `ref_len` bases at `pos`, ALTs end with `<*>`,
    AD lists REF and `alts` (the `<*>` 0 is added), VAF follows from AD."""
    depth = sum(ad)
    vaf = ",".join(f"{n / depth:.6g}" if depth else "0" for n in ad[1:])
    n_alleles = len(alts) + 2
    pl = ",".join(["0"] + ["990"] * (n_alleles * (n_alleles + 1) // 2 - 1))
    return (
        f"{contig}\t{pos}\t.\t{ref(contig, pos, ref_len)}\t{','.join([*alts, '<*>'])}"
        f"\t{qual}\t{filt}\t.\tGT:GQ:DP:AD:VAF:PL"
        f"\t{gt}:{gq}:{depth}:{','.join(map(str, [*ad, 0]))}:{vaf},0:{pl}"
    )


def star_record(contig: str, pos: int, alt: str, *, gt: str, gq: int = 30) -> str:
    """A record with a `*` allele (an upstream deletion on one haplotype). Hand-written: the
    VCF 4.2 shape, not observed from DeepVariant 1.10."""
    return variant(contig, pos, 1, [alt, "*"], gt=gt, gq=gq, ad=[0, 10, 10])


def write_gvcf(path: Path, records: list[str]) -> Path:
    """`records` (sorted here) as a bgzipped, tabix-indexed gVCF; returns its path."""
    order = list(CONTIGS)
    rows = sorted(
        records, key=lambda r: (order.index(r.split("\t")[0]), int(r.split("\t")[1]))
    )
    contigs = "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CONTIGS.items())
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + SAMPLE + "\n"
    plain = path.with_suffix("")
    plain.parent.mkdir(parents=True, exist_ok=True)
    plain.write_text(HEADER + contigs + columns + "".join(f"{r}\n" for r in rows))
    out = pysam.tabix_index(str(plain), preset="vcf", force=True)
    return Path(out)
