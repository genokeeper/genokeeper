# exons.GRCh38.tsv

The exons the slice 3b deletion screen reads, cut from NCBI's MANE (Matched Annotation from
NCBI and EMBL-EBI) GTF.

- URL: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz
- Release: MANE v1.5 (GRCh38)
- sha256 of the GTF: 71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66
- Rebuild: download the file above, then `uv run python -m genokeeper.deletions.exons cut
  <the .gtf.gz>`. The command refuses a file with another sha256.

The cut: every `exon` row of a transcript tagged `MANE_Select` or `MANE_Plus_Clinical` whose
`gene_name` is one of the 192 genes of the ACMG SF v3.3 list (84), the ACMG carrier Tier 3
list (113), 7 of them on both, and the screen genes of `genelists/screen_genes.tsv` (slice
3d: DPYD, whose MANE Select NM_000110.4 has 23 exons, and CHEK2, whose exons 9 and 10 the
del5395 founder deletion removes), on a primary chromosome (chr1-22, chrX, chrY): MANE
v1.5 places SLC37A4's only transcript on `chr11_KN196481v1_fix`, where no caller shard runs,
so SLC37A4 has no row and the screen does not cover it. Columns: gene, transcript (as the
GTF gives it, with its version), the transcript's own exon number, chromosome (`chr` names),
start (0-based, half-open, as BED), end, and the MANE tag. A gene with a Select and a Plus Clinical
transcript has the exons of both, each with its own numbering.
