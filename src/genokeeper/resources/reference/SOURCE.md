# masked_intervals.GRCh38.tsv

The GRC false duplications GIAB masks in its GRCh38 no-alt analysis set, as the intervals the
masked FASTA sets to N. Every alignment, call and calibration runs on that FASTA
(`genokeeper.config.grch38_fasta`); the benchmark subtractions and HG004's acceptance subtract
these intervals.

- URL: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
- Release: GIAB maskedGRC_exclusions_v2 (2021-07), 885,955,190 bytes
- sha256 of the .gz: 74482c74f2a47e8908404028d5305167636c9942fba4129ea0d70c5f04cab5ad
  (GIAB publishes none: computed on the copy fetched 2026-10-30, whose size equals the
  server's Content-Length)
- Compared with: NCBI's GRCh38 no-alt analysis set
  (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`, md5-verified by the bootstrap's
  `grch38_fasta` step), the lab's `genokeeper.config.unmasked_fasta`
- Rebuild: unpack the .gz above, `samtools faidx` it, then `uv run python -m genokeeper.masking
  derive <unmasked .fa> <masked .fa>`.

Why derived and not GIAB's BED: GIAB's `GRC_exclusions_T2Tv2.bed` (2023-02) post-dates the v2
FASTA (2021-07), so it need not describe what the FASTA masks. `genokeeper.masking
masked_intervals` compares the two FASTAs contig by contig (same contigs and lengths, in the
same order), upper-cased (soft masking is case, not sequence): a position the masked FASTA
sets to N where the unmasked one has a base is masked, an N in both is not, and any other
difference refuses. The runs are merged, 0-based half-open (BED), one per row after the
`chrom start end` header.

Derived 2026-10-30: 12 intervals, 1,237,765 bp, on chr9, chr13 (2), chr16, chr21 (5), chrX
and all of `chr11_KI270721v1_random` and `chr22_KI270734v1_random` but each one's first
base, which the masked FASTA keeps.

The bootstrap's `masked_reference` step derives the mask again from the lab's two FASTAs
(`python -m genokeeper.masking verify`) and refuses unless it equals this file, so the mask
every check reads is the one in git. `tests/test_masking.py`'s
`test_no_region_resource_touches_the_mask` checks that no `*.GRCh38.json` or `*.GRCh38.tsv`
region file under `resources/`, and no stage-3 window (`genokeeper.masking.WINDOWS`), lies in
a masked block, where it would read zero depth.

# false_duplications_correct_copy.GRCh38.tsv

The copies GRCh38 keeps of the same GRC false duplications: the correct copy of each masked
block. GIAB's v4.2.1 truths were made on the unmasked reference, where a read of such a haplotype
could go to the false copy, so they can lack a real variant on the correct copy; on the masked
reference every read lands there and the caller reports what the truth lacks. So HG003's
`small_variants` and HG004's SNP F1 leave these intervals out as well as the mask. A run reads
HG003's benchmark without them too, so HG002's cells rate a site there inside v5.0q's benchmark
and any other is uncalibrated; v5.0q was held valid on them by a measurement (BR-030).

- URL: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz
- Release: GIAB genome stratifications v3.6, retrieved 2026-11-05, 201 bytes
- sha256 of the .gz: 44d07f373f6841a99fb85bd712e535e3e1364a033fd37cfd42e9e17988277717
  (computed on that copy; GIAB's md5 list was not fetched)
- sha256 of the decompressed BED: e77d1572ecb059e1f6fa3684ce60e4b9235740749618969ea7a27370e90240eb
- Rebuild: decompress the .gz above under the three `#` lines and the `chrom start end` header
  this file has; its rows are the BED's, unchanged. The file's sha256 is the `correct_copies`
  version of the `small_variants` and `hg002_cells` tables, so an edit stales them.

Checked 2026-11-05 against `GRCh38_no_alt.fa`: 12 intervals, 1,237,625 bp, none overlapping
the mask. Each of the mask's 12 blocks has one, the interval holding 86-100% of the block's
31-mers on one strand (the next best holds at most 1.3%); chrX's block has its copy on chr8 and
chr13's first on chr16, so `chrX` has none.
