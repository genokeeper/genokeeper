# cds.GRCh38.tsv

The coding blocks slice 3c's per-gene coverage reads (`genokeeper.coverage`), cut from NCBI's
MANE (Matched Annotation from NCBI and EMBL-EBI) GTF, the file slice 3b's exons come from.

- URL: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz
- Release: MANE v1.5 (GRCh38)
- sha256 of the GTF: 71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66
- Rebuild: download the file above, then `uv run python -m genokeeper.coverage cut <the
  .gtf.gz>`. The command refuses a file with another sha256.

The cut: every `CDS` row of a transcript tagged `MANE_Select` or `MANE_Plus_Clinical` whose
`gene_name` is one of the coverage genes (`coverage.coverage_genes`, 204): the ACMG SF v3.3
and ACMG carrier Tier 3 genes (190), the eleven `[repeat_loci]` genes of
`config/grading.toml`, SMN1 and SMN2, and CBS, KCNE1, CRYAA and U2AF1 (the reference-limited
genes the masked GRCh38 makes readable; three are on neither list), on a primary chromosome
(chr1-22, chrX, chrY). As in 3b's exons, SLC37A4 has no row: MANE v1.5 places its only
transcript on `chr11_KN196481v1_fix`. Columns are 3b's (`deletions.exons.COLUMNS`): gene,
transcript (with its version), the transcript's own exon number of the block, chromosome
(`chr` names), start (0-based, half-open, as BED), end, and the MANE tag. 4,631 rows.
