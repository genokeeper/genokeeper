# clonal_spans.GRCh38.tsv

The spans of the nine clonal-haematopoiesis genes slice 3c's allele-fraction guard reads
(`genokeeper.genotype.clonal`, `reliability.guard_allele_fraction`), cut from NCBI's MANE
(Matched Annotation from NCBI and EMBL-EBI) GTF, the file slice 3b's exons come from.

- URL: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz
- Release: MANE v1.5 (GRCh38)
- sha256 of the GTF: 71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66
- Rebuild: download the file above, then `uv run python -m genokeeper.genotype.clonal cut
  <the .gtf.gz>`. The command refuses a file with another sha256, and a listed gene without
  a MANE Select transcript.

The cut: for each gene of `config/grading.toml`'s `[af_guard] clonal_genes` (TP53, CHEK2,
DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1, SRSF2), the span of the `exon` rows of its
transcript tagged `MANE_Select` on a primary chromosome, from the first exon's start to the
last exon's end (introns and UTRs inside). Columns: gene, transcript (with its version),
chromosome (`chr` names), start (0-based, half-open, as BED), end. 9 rows.
