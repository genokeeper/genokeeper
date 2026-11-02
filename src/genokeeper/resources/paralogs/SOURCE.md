# Paralog loci: HG002's copy-number truth

`hg002_truth.tsv` holds the copy numbers slice 3d's paralog calibration checks Parascopy and
KILDA against (closing design, "Paralog loci through Parascopy and Gauchian" and "Lp(a)"; the
slice 3d plan's Decisions, "HG002's truth"). HG002, HG003 and HG004 are public GIAB samples.

## STRC

- HG002: GIAB v5.0q (`HG002_GRCh38_v5.0q_stvar.vcf.gz` and its benchmark BED, the files 3b
  pinned). STRC's MANE v1.5 Select span, chr15:43,599,563–43,618,800, lies wholly inside the
  structural-variant and small-variant benchmark regions, and the 16 records over it are SNVs and
  indels of at most 10 bp: no structural variant, two copies.
- HG003 and HG004: two copies each, from research note R2 (STRC: deletion-negative in Sentieon's
  statement on Illumina and Ultima data).

## RHD, CYP2A6, AMY1, C4 and KIV-2

- Assembly: the HG002 Q100 v1.1 diploid assembly, `hg002v1.1.fasta.gz`, from
  `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/`, 1,769,801,343
  bytes, sha256 `28953e5ebaacd445c9e32bbabd8db103a2c10a08af13121159a1527addf84803` (computed on
  the fetch of 2026-10-30; the server publishes no checksum). Contigs `chr*_MATERNAL` are
  haplotype 1, `chr*_PATERNAL` haplotype 2.
- Method: each query, cut from GRCh38 (no-alt analysis set), aligned with minimap2 2.28-r1209
  `-x asm5 -c -N 500 -p 0.1 --secondary=yes` to the whole assembly, and to GRCh38 as a control;
  identity is matching bases over the alignment block, coverage the aligned share of the query.
  - RHD and CYP2A6: their MANE v1.5 Select spans, hits at 98% identity or more over 90% or more
    of the span (each differs from its paralog, RHCE and CYP2A7, by more than 2%: those copies
    align at 96% or less, or not at all). One copy on each haplotype.
  - AMY1: AMY1A's MANE span at 99% or more, which counts every AMY1 copy and no AMY2 copy (the
    GRCh38 control finds its own three, AMY1A, AMY1B and AMY1C). Three on each haplotype.
  - C4: C4A's MANE span at 98% or more, counting C4A and C4B together (the control finds two).
    Two on each haplotype.
  - KIV-2: KILDA's `KIV2_hg38.bed` (tag v1.2.3) is the whole 35 kb array, chr6:160,611,535–
    160,646,545, so the query is one unit, chr6:160,617,661–160,623,206 (5,546 bp, from one
    unit's start to the next). Hits at 95% identity or more over 80% or more of the unit: one
    unit type aligns over 83–88% only (a divergent 640 bp head), and each array ends in a
    1,470 bp partial unit (27%), which is not counted. The control finds six, the count KILDA's
    reference k-mer file assumes (`KIV2_hg38_kmers_6copies_specific.tsv`). 14 maternal, 24
    paternal.
- UGT2B17 is not in the table: the slice does not cover it (Parascopy's `UGT2B15` model covers
  21% of its MANE span).

# genes.GRCh38.tsv

The MANE v1.5 Select spans of slice 3d's paralog genes (`deletions.exons.PARALOG_SPAN_GENES`),
which the depth screen's owned regions and the paralog readers read (`load_spans`,
`locus_span`).

- Source: the MANE GTF `resources/deletions/SOURCE.md` names (MANE v1.5, sha256
  `71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66`).
- Rebuild: `uv run python -m genokeeper.deletions.exons spans <the .gtf.gz>`, which refuses a
  file with another sha256.
- The cut: each gene's `transcript` row tagged `MANE_Select` on chr1-22, chrX or chrY, start
  0-based half-open. MANE v1.5 also places a C4A Select transcript on `chr6_GL000254v2_alt`;
  it is not read, so chr6's C4A span stands. `GBA1_ex9-11` is the union of GBA1's MANE
  Select exons 9 to 11 (ENST00000368373.8), the part GBAP1 copies.
- `OWNED_CYP2D` (in `deletions/exons.py`, not this file) is CYP2D6's span from 5 kb before it
  to 25 kb after it, holding CYP2D7 and CYP2D8P: NCBI Gene 1564 and 1568, GRCh38
  chr22:42,139,576-42,144,483 and chr22:42,149,883-42,155,016 (1-based), read from
  E-utilities esummary on 2026-10-30.
