#!/usr/bin/env bash
# Install the tools and databases `genokeeper run` needs under $GENOKEEPER_HOME (else
# $GENOMELAB_HOME, its name before the rename; default ~/genome-lab).
#
#   scripts/bootstrap.sh                 # install and download whatever is missing
#   scripts/bootstrap.sh --with-gnomad   # also build the gnomAD FAF table
#   scripts/bootstrap.sh --with-calibration  # also fetch HG003 and its truth (stage 2b)
#   scripts/bootstrap.sh --with-hg002    # also stream and downsample HG002 (stage 3's calibration)
#   scripts/bootstrap.sh --with-vep      # also fetch VEP's cache and slice 1e's score files
#   scripts/bootstrap.sh --with-sv       # also fetch gnomAD-SV (slice 3b), unless it is in prefetch/3b/
#   scripts/bootstrap.sh --with-masked-reference  # also fetch the masked GRCh38 (slice 3c,
#                                        # 886 MB), unless it is in prefetch/3c/
#   scripts/bootstrap.sh --with-hg004    # also fetch HG004's truth and stream its library (3c)
#   scripts/bootstrap.sh --with-typing   # also fetch slice 3f's HHV-6 and hmmix files (3.5 GB),
#                                        # unless they are in prefetch/3f/
#   scripts/bootstrap.sh --with-publications  # also fetch LitVar2's variants file and ClinVar's
#                                        # citations (about 2.2 GB, publications only)
#   scripts/bootstrap.sh --check         # list what is missing; reads only, exits 0
#   scripts/bootstrap.sh --refresh STEP  # fetch STEP's newest release, whatever the manifest
#                                        # says, and nothing else (`genokeeper refresh` calls it)
#
# The flags combine in any order; --check with --with-gnomad, --with-calibration,
# --with-hg002, --with-vep, --with-sv, --with-masked-reference, --with-hg004, --with-typing
# or --with-publications still only lists. --refresh STEP takes no
# other flag: it runs that one step, ignoring its manifest rows (a step that is done is
# fetched again), files the new release under its dated folder, prints
# `refreshed: STEP <dated folder>` on stdout and exits 0; any failure exits non-zero and
# leaves the previous folders and manifest rows as they were. STEP is one of: clinvar
# clingen clingen_dosage clingen_regions cpic pharmgkb gwas medlineplus litvar
# clinvar_citations gencc. Every one fetches a whole file from a URL written in this script.
#
# Steps, each skipped when its outputs exist with the sha256 recorded in
# $GENOKEEPER_HOME/bootstrap-manifest.tsv (step<TAB>path<TAB>sha256):
#   tools       bcftools, samtools, bgzip and tabix (brew installs only the missing ones);
#   clinvar     the ClinVar GRCh38 weekly VCF with its .tbi, md5-verified and read through
#               by the run's own ClinVar reader (genokeeper.panel.build) before it is filed,
#               in clinvar/<##fileDate as YYYYMMDD>/clinvar_GRCh38.vcf.gz;
#   clingen     ClinGen's own gene-disease validity export (CSV), checked by genokeeper's
#               loader, in clingen/<its FILE CREATED date as YYYYMMDD>/gene-validity.csv,
#               filed only when a line other than FILE CREATED differs from the newest copy;
#   grch38_fasta  NCBI's GRCh38 no-alt analysis set, md5-verified against its directory's
#               md5checksums.txt, unpacked with a .fai to genokeeper.config.unmasked_fasta
#               (the one definition of the path, where `genokeeper build-faf` reads it and
#               CRAMs written before the masked reference decode);
#   tools38     bwa 0.7.19 and mosdepth 0.3.14 (brew), rtg-tools 3.13 (its release zip,
#               sha256-pinned, with Homebrew's openjdk) in tools/rtg-tools-3.13/; each
#               version checked (another bwa or mosdepth version is a warning); no
#               manifest row;
#   giab38      the GIAB v3.6 GRCh38 BEDs genokeeper.regions lists ([strata] and the XY,
#               gap and mappability sources), md5-verified and tabix-indexed as published,
#               and the region BEDs it derives (chrY X-degenerate, autosomes), in
#               giab/v3.6-GRCh38/;
#   dv_model    the DeepVariant Complete Genomics T7 checkpoint's four files, md5-pinned, in
#               deepvariant/complete-t7/ (files already there with the pinned md5 are kept);
#   cgd         the Clinical Genomic Database's CGD.txt.gz, fetched whole, in
#               cgd/<its Last-Modified date as YYYYMMDD>/;
#   masked_reference  slice 3c: GIAB's GRCh38 no-alt analysis set with the GRC false
#               duplications masked (maskedGRC_exclusions_v2, 886 MB, sha256-pinned), moved
#               from prefetch/3c/ when the copy there has the pin, else (only with
#               --with-masked-reference) downloaded through a resumed .part; unpacked to
#               genokeeper.config.grch38_fasta with a .fai, its mask derived against the
#               unmasked FASTA (the grch38_fasta step comes first) into
#               genokeeper.masking.masked_bed and refused unless it equals the committed
#               mask, bwa-indexed beside it (`python -m genokeeper.masking index`), and rtg's
#               SDF of it built (`rtg format`, stamped with its sha256) where
#               `python -m genokeeper.calibrate layout` names it; recorded only once all
#               are there. After tools38 and pharmcat.
#   pharmcat    PharmCAT 3.4.0's release assets pharmcat_positions_3.4.0.vcf.bgz and its .csi,
#               sha256-pinned to the digests GitHub publishes, in pharmcat/3.4.0/;
#   cpic        CPIC's gene-drug pairs and MT-RNR1's allele names, JSON rows from CPIC's
#               API, fetched whole and never resumed, checked by genokeeper's reader, in
#               cpic/<fetch date as YYYYMMDD, UTC>/, filed only when a file differs from
#               the newest copy;
#   pharmgkb    ClinPGx's annotations download (summaryAnnotations.zip), fetched
#               whole and never resumed, checked by genokeeper's reader, in
#               pharmgkb/<its CREATED date as YYYYMMDD>/;
#   gwas        the GWAS Catalog's associations with their ontology annotations (one zipped
#               TSV), fetched whole and never resumed, checked by genokeeper's reader, in
#               gwas/<its member's build date as YYYYMMDD>/;
#   medlineplus MedlinePlus Genetics' summaries (ghr-summaries.xml), fetched whole and never
#               resumed, checked by genokeeper's loader, in medlineplus/<fetch date as
#               YYYYMMDD, UTC>/;
#   images      every image in genokeeper.containers.IMAGES, pulled by its digest (Docker
#               refuses other content), or built from its Dockerfile under docker/ where no
#               image of the pinned version exists for the platform it runs on
#               (genokeeper.containers.BUILT; each IMAGES comment says why); no manifest row
#               (Docker's store);
#   t1k         IPD-IMGT/HLA 3.66.0's hla.dat.zip, fetched whole and sha256-pinned, and T1K's
#               HLA index built from it in T1K's image (genokeeper.pgx.hla.t1k_index), in
#               t1k/IPD-IMGT-HLA-3.66.0/;
#   pgs         slice 1d's polygenic scores, in pgs/ (genokeeper.prs.pgsc names each path):
#               Nextflow 25.10.0's -dist launcher, its nf-schema 2.0.0 plugin unpacked under
#               pgs/nxf-home/ (Nextflow runs offline) and pgsc_calc at the commit its tag
#               v2.3.0 points to, each sha256-pinned, with a link to Homebrew's openjdk@21;
#               the chosen scores' harmonised GRCh38 scoring files (genokeeper.prs.scores),
#               fetched whole and checked against the .md5 the PGS Catalog publishes
#               beside each, fetched with it; the HGDP+1kGP
#               reference panel (16 GB, md5-pinned, no manifest row: like a FASTQ, it is
#               under its final name only once its md5 matched), moved from
#               prefetch/1d/pgsc_calc/ when the copy there has the pinned md5, else
#               downloaded; the panel's GRCh38 variant file (every variant, for the GWAS
#               associations' alleles) and its ancestry positions, filtered by plink2 in
#               pgsc_calc's image (python -m genokeeper.prs.pgsc ancestry-positions). After
#               images: the filter runs in a pulled image.
#   clingen_dosage  ClinGen's dosage sensitivity curations (slice 1e's haploinsufficient
#               genes), fetched whole and never resumed, checked by genokeeper's reader, in
#               clingen-dosage/<the date it names as YYYYMMDD>/, filed only when a line
#               other than that date differs from the newest copy;
#   clingen_regions  ClinGen's region curations and recurrent CNV regions (slice 3d's large
#               depth events; approved by the user, A5), both fetched whole and never
#               resumed, checked by genokeeper's readers, in clingen-regions/<the region
#               list's date as YYYYMMDD>/, filed only when either file changed;
#   gencc       GenCC's gene-disease validity submissions in their new format (slice 3f),
#               fetched whole and never resumed (the export is built per request, 20 a day),
#               checked by genokeeper's reader, in gencc/<its newest submitted_run_date as
#               YYYYMMDD>/gencc-submissions.tsv, filed only when it differs from the newest
#               copy; the copy in prefetch/3f/ is moved in when it has the pinned size and
#               sha256 (not under --refresh, which fetches the newest);
#   litvar, clinvar_citations  only with --with-publications or --refresh (they are no
#               grading input): LitVar2's whole variants file and ClinVar's var_citations.txt,
#               fetched whole and never resumed (LitVar2 is rebuilt in place), checked by
#               genokeeper.evidence.publications, in litvar/<its Last-Modified date as
#               YYYYMMDD, UTC>/ and clinvar-citations/<the same>/;
#   vep         slice 1e (genokeeper.annotate.vep names each path), only with --with-vep
#               (without it the step says what is missing, and --check lists it): Ensembl's
#               indexed VEP 116 GRCh38 cache (27.6 GB), REVEL v1.3, AlphaMissense,
#               SpliceVault's GRCh38 data with its .tbi, gnomAD v4.1's base-level pext and
#               UCSC's NMD Escape track for MANE, each sha256-pinned, moved from
#               prefetch/1e/ when the copy there has the pin, else downloaded through a
#               resumed .part; then, in vep/116/, the cache unpacked, REVEL sorted on its
#               GRCh38 column and tabix-indexed, AlphaMissense indexed, pext as a
#               tabix-indexed BED with each gene's largest value, and the NMD track as BED
#               by bigBedToBed in its pinned image. After images.
#   sv          slice 3b (genokeeper.deletions.files names each path): gnomAD v4.1's SV sites VCF
#               (1.74 GB) with its .tbi in sv/, each sha256-pinned and moved from prefetch/3b/
#               when the copy there has the pin, else downloaded through a resumed .part, and
#               then only with --with-sv (without it the step says what is missing, and --check
#               lists it). After images.
#   typing      slice 3f (genokeeper.typing_files names each path): NCBI's HHV-6A (NC_001664.4)
#               and HHV-6B (NC_000898.1) FASTA files, concatenated into hhv6/hhv6.fa and
#               bwa-indexed; hmmix 0.9.2's four hg38 files (Zenodo 11212339) and 22 archaic
#               BCFs with their .csi (Zenodo 13368126) in hmmix/0.9.2/ (3.5 GB in all); each
#               sha256-pinned and moved from prefetch/3f/ when the copy there has the pin,
#               else downloaded through a resumed .part, and then only with --with-typing
#               (without it the step says what is missing, and --check lists it). Then,
#               when their outputs are missing, the supervised-admixture reference in
#               ancestry/ and the ROH allele-frequency files in roh/, built from the HGDP+1kGP
#               panel in pgsc_calc's plink2 image and fastmixture's
#               (python -m genokeeper.typing_files build-panel, build-af). After images, pgs.
#   panelapp    Genomics England PanelApp's Intellectual disability panel (285, 2.8 MB),
#               one JSON file, version- and sha256-pinned, in panelapp/285/<version>/. Its
#               URL serves the panel's current version: a newer one fails the pin and
#               is not filed, and neither that nor a failed fetch stops the bootstrap, as
#               the run then records the layer as not run;
#   gnomad_faf  no manifest row: the table keeps its own manifest. Only with --with-gnomad,
#               `genokeeper build-faf` streams gnomAD v4.1 joint (0.88 TB, keeps a few GB)
#               into gnomad/4.1/ for the chromosomes the table lacks or whose bucket md5
#               changed; without the flag the step says the table is missing, and --check
#               lists the missing chromosomes.
#   hg003       stage 2b's calibration inputs, only with --with-calibration (without it the
#               step says what is missing, and --check lists it), into the places
#               `python -m genokeeper.calibrate layout` names: GIAB's HG003 v4.2.1 GRCh38
#               benchmark (VCF, .tbi, BED; md5-pinned) and HG003's
#               DNBSEQ-T7 FASTQ (147 GB, md5-pinned, no manifest row: a file under its final
#               name was verified, as it is renamed there only after its md5 matched). A
#               transfer cut part-way leaves <name>.part, which the rerun resumes with
#               `curl -C -`; a FASTQ is not fetched again once its downsampled pair exists.
#
#   parascopy   stage 3's Parascopy data, on every bootstrap (52 MB): Zenodo record 15019940's
#               GRCh38_v1.7.tar.gz (homology table) and models_GRCh38_1KGP_v1.7.tar.gz,
#               md5-pinned, kept in parascopy/ and unpacked beside the tarballs; a cut
#               transfer leaves <name>.part, which the rerun resumes;
#   hg002       stage 3's calibration inputs, only with --with-hg002 (without it the step says
#               what is missing, with the download size, and --check lists it): GIAB's HG002
#               v5.0q GRCh38 small-variant truth (VCF, .tbi, benchmark BED; md5-pinned) in
#               giab/HG002-v5.0q-GRCh38/, then each public HG002 library (ENA, PRJNA1465120)
#               through `python -m genokeeper.calibrate.hg002 fetch`: each FASTQ downloaded
#               whole over 24 connections into a hidden file (resumed from a journal of
#               finished chunks after a kill), md5-checked, downsampled, then deleted (300 GB
#               read, about 170 GB kept; one raw file, about 150 GB, on disk at a time), a
#               library skipped when its downsampled pair exists or it was aligned, and every
#               library left for a rerun before a `genokeeper run` recorded the bases (the
#               bootstrap goes on). No manifest row for the libraries: a file under its final
#               name was verified.
#   hg002_sv    slice 3b's truth: GIAB's HG002 v5.0q GRCh38 structural-variant VCF, .tbi and
#               benchmark BED (50 MB, sha256-pinned, moved from prefetch/3b/ or downloaded)
#               into the same giab/HG002-v5.0q-GRCh38/, recorded under its own step
#               hg002_sv_truth, only with --with-hg002 (without it the step says so).
#   hg004       slice 3c's HG004 inputs, only with --with-hg004 (without it the step says
#               what is missing, with the sizes, and --check lists it): GIAB's HG004 v4.2.1
#               GRCh38 benchmark (VCF, .tbi, BED; 163 MB, sha256-pinned, moved from
#               prefetch/3c/ or downloaded) in giab/HG004-v4.2.1-GRCh38/, recorded as
#               hg004_truth, then HG004's PCR-free library through
#               `python -m genokeeper.calibrate.hg002 fetch hg004` (144 GB from CNGB) unless its
#               downsampled pair or its record is there (left for a rerun, as HG002's, before
#               a `genokeeper run` recorded the bases).
#
# brew runs with HOMEBREW_NO_ANALYTICS=1. Sourcing the script defines its functions without
# running any step.

set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
lab=${GENOKEEPER_HOME:-${GENOMELAB_HOME:-$HOME/genome-lab}}
manifest=$lab/bootstrap-manifest.tsv
check=0
with_gnomad=0
with_calibration=0
with_hg002=0
with_vep=0
with_sv=0
with_masked_reference=0
with_hg004=0
with_typing=0
with_publications=0
refresh=
refresh_steps="clinvar clingen clingen_dosage clingen_regions cpic pharmgkb gwas medlineplus litvar clinvar_citations gencc"
clinvar_base=https://ftp.ncbi.nlm.nih.gov/pub/clinvar
clingen_url=https://search.clinicalgenome.org/kb/gene-validity/download
clingen_file=gene-validity.csv
grch38_base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
grch38_name=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# Versions and checksums are pinned here: a file that differs is refused, a bwa or mosdepth
# version that differs is a warning.
bwa_version=0.7.19
mosdepth_version=0.3.14
rtg_version=3.13
rtg_url=https://github.com/RealTimeGenomics/rtg-tools/releases/download/$rtg_version/rtg-tools-$rtg_version-nojre.zip
# sha256 of the release zip, computed on its download of 2026-10-19 (GitHub publishes none).
rtg_sha256=bc5c6badb07d7e20d1c5c557bd6d571a022bbd9f58fa1e3840bcff9431a18f96
giab38_base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
giab38_md5s=GRCh38-genome-stratifications-md5s.txt
# The DeepVariant Complete Genomics T7 checkpoint: file:md5, the md5s the bucket listing
# (storage/v1/b/deepvariant/o?prefix=complete-case-study-testdata/complete-t7/) gave as
# md5Hash on 2026-10-19. A file whose md5 differs is fetched again; a fetched one that still
# differs is refused, so a re-uploaded checkpoint never changes the caller silently.
dv_model_base=https://storage.googleapis.com/deepvariant/complete-case-study-testdata/complete-t7
dv_model_files="example_info.json:15ae88dfbe01fe3fcd624592d5db3e67 model.example_info.json:0af4100a148992efdfe705d48f010623 weights-51-0.995354.ckpt.data-00000-of-00001:c178f108e2a63121563dd3899ff0b3b1 weights-51-0.995354.ckpt.index:5cb901452d3974ae343fd01358941408"
cgd_url=https://research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz
cgd_file=CGD.txt.gz
# PharmCAT's positions (genokeeper.pgx.pharmcat.positions_path), name:sha256. The sha256s are
# the asset digests GitHub's release API gives for tag v3.4.0 (published 2026-07-14, read
# 2026-10-20), so a re-uploaded asset is refused.
pharmcat_version=3.4.0
pharmcat_base=https://github.com/PharmGKB/PharmCAT/releases/download/v$pharmcat_version
pharmcat_files="pharmcat_positions_$pharmcat_version.vcf.bgz:01a28fe31675c9e8a8b18bcab4afb82a62abdfea6606856aeca0e0378e7d6b58 pharmcat_positions_$pharmcat_version.vcf.bgz.csi:eb536edd898613977e7429469e4a8fceb31afc85a839f98423272f90a9491ac4"
# CPIC's current data from its API (PostgREST), name:query: the gene-drug pairs with their
# CPIC levels, and MT-RNR1's allele names (genokeeper.evidence.pgx.CPIC_FILES). CPIC updates
# them in place and publishes no checksum, so each fetch is filed under its date with its
# sha256 in the manifest. A query names a gene or none, never anything of the user's.
cpic_api=https://api.cpicpgx.org/v1
cpic_files="cpic_pairs.json:pair_view?select=genesymbol,drugname,cpiclevel MT-RNR1_alleles.json:allele?genesymbol=eq.MT-RNR1&select=name"
# ClinPGx's (formerly PharmGKB's) annotations, CC BY-NC-SA 4.0 (the live file; the frozen
# clinicalAnnotations.zip was CC BY-SA 4.0, and genokeeper still reads it). Rebuilt in place
# under one URL with no published checksum, so each fetch is filed under the CREATED date it
# names with its sha256 in the manifest.
pharmgkb_url=https://api.clinpgx.org/v1/download/file/data/summaryAnnotations.zip
pharmgkb_file=summaryAnnotations.zip
# The GWAS Catalog's associations with their ontology annotations (genokeeper.evidence.gwas).
# The latest release is rebuilt in place under one URL, and its release directory publishes
# no checksum (read 2026-10-21), so each fetch is filed under the build date its member
# carries, with its sha256 in the manifest.
gwas_file=gwas-catalog-associations_ontology-annotated-full.zip
gwas_url=https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/$gwas_file
# MedlinePlus Genetics' gene and condition summaries (genokeeper.evidence.medlineplus), public
# domain. Rebuilt daily under one URL with no published checksum, so each fetch is filed under
# its date with its sha256 in the manifest.
medlineplus_url=https://medlineplus.gov/download/ghr-summaries.xml
medlineplus_file=ghr-summaries.xml
# The publication sources (genokeeper.evidence.publications): LitVar2's whole variants file
# (about 2 GB, public domain, rebuilt in place) and ClinVar's citations table. Each is one
# URL, the same for every user: no request names a variant. The per-rsID LitVar endpoint is
# never called.
litvar_url=https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz
litvar_file=litvar2_variants.json.gz
citations_url=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
citations_file=var_citations.txt
# IPD-IMGT/HLA's hla.dat, which T1K's HLA index (genokeeper.pgx.hla.t1k_index) is built from:
# release 3.66.0 (tag v3.66.0-alpha, 2026-10-14). The sha256 was computed on its download of
# 2026-10-20, whose md5 matched the release's md5checksum.txt
# (d854cf6b3b63a12e3de178aa411c7d35), so a re-uploaded release is refused.
imgt_hla_version=3.66.0
imgt_hla_url=https://raw.githubusercontent.com/ANHIG/IMGTHLA/v$imgt_hla_version-alpha/hla.dat.zip
imgt_hla_sha256=509dd60b30666e70b217ec499f10d27a575427c5015da601e0ac5269182840fd
# Slice 1d (genokeeper.prs.pgsc): Nextflow, the sha256 GitHub's release API gives for the
# -dist asset of v25.10.0 (published 2025-10-22); nf-schema, whose zip has no published
# sha256: computed on its download of 2026-10-21, whose sha512 matched the release's
# nf-schema-2.0.0-meta.json (ecca8c50...); pgsc_calc, GitHub's archive of the commit the tag
# v2.3.0 points to, its sha256 computed on two downloads of 2026-10-21 that agreed. All
# read 2026-10-21. Nextflow 25.10 runs on Java 17 to 25: Homebrew's openjdk is 27.
nextflow_version=25.10.0
nextflow_url=https://github.com/nextflow-io/nextflow/releases/download/v$nextflow_version/nextflow-$nextflow_version-dist
nextflow_sha256=294376ec555695ee0b92e21477600f97c113f2d1ed3fb1f480daf2ee439c4626
nf_schema_version=2.0.0
nf_schema_url=https://github.com/nextflow-io/nf-schema/releases/download/$nf_schema_version/nf-schema-$nf_schema_version.zip
nf_schema_sha256=8c68d1a8bcfa670758c8f591bb8f885a4b181d0b4aa90e1316c8798647705c98
pgsc_calc_version=2.3.0
pgsc_calc_commit=72ee54f2119d2fe65ac671327342dfc412bd0a3b
pgsc_calc_url=https://github.com/PGScatalog/pgsc_calc/archive/$pgsc_calc_commit.tar.gz
pgsc_calc_sha256=88c354e24ff383d5955386e6f8dee8ee2772ef4277388c64aa8bcbcaaf848ce3
nextflow_java=openjdk@21
# StellarPGx (MIT), GitHub's archive of master on 2026-11-09, the sha256 of those bytes.
stellarpgx_commit=1952a4494fd6e9f809be87baa0a56c379e1d757b
stellarpgx_url=https://github.com/SBIMB/StellarPGx/archive/$stellarpgx_commit.tar.gz
stellarpgx_sha256=dc7df7667861f94977f4ec29221e76e1d7f60b3684de1e05a0b6abb3f214b97b
# The chosen scores (genokeeper.prs.scores), whose harmonised GRCh38 scoring files the PGS
# Catalog re-harmonises in place, each with its md5 beside it (<file>.md5).
pgs_scores_base=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores
pgs_ids="PGS000018 PGS005168 PGS002765 PGS003867 PGS003978 PGS005331
PGS005240 PGS003852 PGS004886 PGS012590 PGS000903 PGS002724
PGS001790 PGS003972 PGS005217 PGS001796 PGS001792 PGS005113
PGS001789 PGS002768 PGS002378 PGS004156 PGS004606 PGS004884
PGS019980 PGS001794 PGS012555 PGS004873 PGS000754 PGS002726
PGS004254 PGS004253 PGS000820 PGS001787 PGS002755 PGS001788
PGS000705 PGS004614 PGS019967 PGS019973 PGS019983 PGS004885
PGS000134 PGS002802 PGS000882"
# The research scores (genokeeper.prs.scores.load_research_scores): scores outside the chosen
# set, a separate step so that adding one never re-fetches a chosen file.
pgs_research_ids="PGS002746 PGS003753 PGS002790 PGS003724 PGS004427 PGS003573
PGS002231 PGS002786 PGS000133 PGS002342 PGS004451 PGS002209
PGS000205 PGS000908 PGS004161 PGS003492 PGS004799 PGS002143
PGS002254 PGS002255 PGS002168 PGS002193 PGS002133 PGS003485
PGS002201 PGS002813 PGS002033"
# The HGDP+1kGP reference panel (15,953,388,718 B, 2023-10-05), md5 as the PGS Catalog's
# resources/md5s.txt publishes it (read 2026-10-21).
pgs_panel_name='pgsc_HGDP+1kGP_v1.tar.zst'
pgs_panel_url=https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/$pgs_panel_name
pgs_panel_md5=2bd3586e9a09700f9bcbe62ba4dcee49
# Its GRCh38 variant file, every variant, kept beside the ancestry positions for the GWAS
# associations' alleles (genokeeper.prs.pgsc.PANEL_PVAR_NAME).
pgs_panel_pvar=GRCh38_HGDP+1kGP_panel.pvar.zst
# HG003 (NA24149, GIAB's Ashkenazi father) on the original DNBSEQ-T7 (flowcell E200029822),
# PCR-free PE150: CNGBdb project CNP0008319, run CNR1441268 ("WGS of NA24149: PCR-Free").
# Each FASTQ as name:md5:size. The folder is the run's download link in CNGBdb's own API
# (db.cngb.org/cnsa/ajax/run_tableview/?q=CNR1441268); the md5s are CNGBdb's, the run's
# file_md5 and file2_md5 in the project's metadata table
# (ftp.cngb.org/pub/CNSA/data7/public_info/CNP0008319/metadata_CNP0008319_experiment.tsv);
# the sizes the server's Content-Length. All read 2026-10-20.
hg003_fastq_base=https://ftp.cngb.org/pub/CNSA/data7/CNP0008319/CNS1414628/CNX1267209/CNR1441268
hg003_fastq_files="E200029822_L01_UDB-451_1.fq.gz:5a1593a2ca3f72d5d89d7e0888e84933:70700990194 E200029822_L01_UDB-451_2.fq.gz:0e011c6e9ea4f9024d91b107be7ee968:75959390481"
# GIAB's HG003 v4.2.1 GRCh38 benchmark (chr1-22), name:md5. GIAB publishes no checksum for
# these files (unchanged since December 2020): the md5s were computed on their download of
# 2026-10-20, so a file served differently later is refused.
hg003_truth_base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38
hg003_truth_files="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz:ea98d8b9de762e050e34666462bf3452 HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi:769599d9863a3290ea1da97297dbbc91 HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed:26f7a0a353c653e4e2979654586110b1"
# Parascopy's data (Zenodo record 15019940, 52 MB), name:md5, as Zenodo's file list gives them.
parascopy_base=https://zenodo.org/records/15019940/files
parascopy_files="GRCh38_v1.7.tar.gz:a95bf674f43317d3a4c1b8ddbb140945 models_GRCh38_1KGP_v1.7.tar.gz:244110d8fa883cf11334527ec2383498"
# GIAB's HG002 v5.0q GRCh38 small-variant truth (T2T Q100 based), name:md5. The md5s were
# transcribed from GIAB's directory by a research run and are confirmed on the first real
# download, which refuses a mismatch.
hg002_truth_base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q
hg002_truth_files="HG002_GRCh38_v5.0q_smvar.vcf.gz:c71acc71069bf7cd7f51eb8fb0c1a1ab HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi:58cab2a06e29b74a5bdc190ead079083 HG002_GRCh38_v5.0q_smvar.benchmark.bed:3366858af85875cbb3c579412ffa4c56"
# ClinGen's dosage sensitivity curations (genokeeper.evidence.clingen.DOSAGE_URL).
clingen_dosage_file=ClinGen_gene_curation_list_GRCh38.tsv
clingen_dosage_url=https://ftp.clinicalgenome.org/$clingen_dosage_file
# ClinGen's region curations and recurrent CNV regions (genokeeper.evidence.clingen
# REGIONS_FILE and RECURRENT_FILE), both CC0.
clingen_regions_file=ClinGen_region_curation_list_GRCh38.tsv
clingen_regions_url=https://ftp.clinicalgenome.org/$clingen_regions_file
clingen_recurrent_file=ClinGen_recurrent_CNV_V2.1-hg38.bed
clingen_recurrent_url=https://ftp.clinicalgenome.org/$clingen_recurrent_file
# GenCC's submissions export, new format (genokeeper.check.GENCC_URL), and the copy fetched
# into prefetch/3f/ on 2026-10-30 (approved by the user): its size and sha256.
gencc_url=https://thegencc.org/download/action/submissions-export-tsv?format=new
gencc_file=gencc-submissions.tsv
gencc_size=25164448
gencc_sha256=6b099f2497837c8bb0b62930bcfff13c952869dd52635511da181521ce1cde17
# Genomics England PanelApp's Intellectual disability panel (genokeeper.sources.panelapp):
# version 11.93, 2,867 genes (approved by the user on 2026-11-08). The API lists an entity's
# `evidence` in no fixed order, so one version differs byte for byte between two fetches: the
# sha256 is of the file sorted (`step_panelapp`), as two fetches of 2026-11-09 agree.
panelapp_version=11.93
panelapp_url='https://panelapp.genomicsengland.co.uk/api/v1/panels/285/?format=json'
panelapp_sha256=99bc5d75640874c78b2a932caedfb4a689db8b6ba5fe97923763104dd3436c8e
# Slice 1e's downloads, one per line: name, prefetch/1e/ subfolder, size, sha256, URL. Each
# sha256 was computed on the download of 2026-10-21 to 2026-10-23 that matched what its
# publisher gives: the size Ensembl's FTP lists and its CHECKSUMS (BSD sum 56036 26996736)
# for the VEP cache; GCS's md5 for AlphaMissense (n9Fnc18Wobh9pus+TCX8tQ==) and pext
# (hezEjdhZvz8MZNclYcmWiw==); Ensembl's CHECKSUMS for SpliceVault (57027 863191 and
# 55851 731); the size the server gave for REVEL and UCSC's NMD track, whose publishers give
# no checksum. So a file re-published with other bytes is refused.
vep_release=116
vep_cache_tarball=homo_sapiens_vep_116_GRCh38.tar.gz
revel_zip=revel-v1.3_all_chromosomes.zip
alphamissense_file=AlphaMissense_hg38.tsv.gz
splicevault_file=SpliceVault_data_GRCh38.tsv.gz
pext_tsv=gnomad.pext.gtex_v10.base_level.tsv.gz
nmd_bb=nmdEscMane.bb
vep_downloads="$vep_cache_tarball vep 27644657162 014c7dd9bb5ad06665866d62eb80f31ca761197bcb9b59280300676f996e600d https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/$vep_cache_tarball
$revel_zip revel 667102707 815801d832e4ab76922dddfd901d7f24a755df57078a7138aa14f615da86a8a2 https://rothsj06.dmz.hpc.mssm.edu/$revel_zip
$alphamissense_file alphamissense 642961469 0516cfd71c0767ac8f9c469252d429000e94e02c008b6e3a46d4b4646fcd3475 https://storage.googleapis.com/dm_alphamissense/$alphamissense_file
$splicevault_file splicevault 883907365 04d487f484dc2dac6d0b98f7a417cf8f7a5c98550afd602cde547d07af3a1de6 https://ftp.ensembl.org/pub/current_variation/SpliceVault/$splicevault_file
$splicevault_file.tbi splicevault 747762 571852b7a0bcadeddc7d0a7bd429ad7b5ec1e40c90de0943b1cc8d3c22738cb0 https://ftp.ensembl.org/pub/current_variation/SpliceVault/$splicevault_file.tbi
$pext_tsv pext 198452726 7e56d4613c9e59a459bb50d53cc43698114a5598f95c5d64d465889beaab0fd8 https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/pext/$pext_tsv
$nmd_bb nmd 3153487 1089835208c6cbac25e7cb072b16f7e4d48d99e1bfaf6aa3e8bcfef7df6da784 https://hgdownload.soe.ucsc.edu/gbdb/hg38/nmd/$nmd_bb"

# Slice 3b's downloads, one per line: name, size, sha256, URL. Each sha256 and size were
# computed on the bytes in prefetch/3b/ on 2026-10-27 and matched against the publisher's
# digest that day (GCS's md5 for gnomAD-SV, GIAB's md5 listing for the v5.0q SV files), and
# each URL's Content-Length equals its size. So a file re-published with other bytes is refused.
sv_base=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv
sv_downloads="gnomad.v4.1.sv.sites.vcf.gz 1743544245 26a70c309a853c77c043c8d710016f97a326d4a55cc89a08bc45384e00de2b04 $sv_base/gnomad.v4.1.sv.sites.vcf.gz
gnomad.v4.1.sv.sites.vcf.gz.tbi 511659 d41b761f20107b26279a5e02ebb178a31a710f7b47f2ea7c37753d21cd893022 $sv_base/gnomad.v4.1.sv.sites.vcf.gz.tbi"
hg002_sv_downloads="HG002_GRCh38_v5.0q_stvar.vcf.gz 48040012 d66b2d2496ff5418763813d3195599007dbff950514d10e5bc27b6c8b76e34b8 $hg002_truth_base/HG002_GRCh38_v5.0q_stvar.vcf.gz
HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi 1554436 57898d17e44ecd5c286aeb863f27e613244396a3ad4b04ca79c8add647b8a81d $hg002_truth_base/HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi
HG002_GRCh38_v5.0q_stvar.benchmark.bed 230609 2f75ce942e1dd9e1a443e4e04fac640104aae5e23fba24d7d3f2c1f1f9985b00 $hg002_truth_base/HG002_GRCh38_v5.0q_stvar.benchmark.bed"
# Slice 3c's downloads. The sizes are those the plan read from each server's Content-Length;
# GIAB publishes no checksum for these files, so each sha256 was computed on the bytes
# fetched into prefetch/3c/ on 2026-10-30, whose sizes equal those, and a file re-published
# with other bytes is refused. The masked FASTA's pin is also genokeeper.masking.MASKED_SHA256.
masked_name=GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
masked_size=885955190
masked_sha256=74482c74f2a47e8908404028d5305167636c9942fba4129ea0d70c5f04cab5ad
masked_url=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/$masked_name
hg004_truth_base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38
hg004_truth_downloads="HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz 148948118 13c113d93eb2884ebb97573cbe558554d3fb0a8748a6d2263cf4616312ac14ff $hg004_truth_base/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi 1653057 05925a6b43cd6894b65a47a3b9612fee634d6a7164fdd02b5f4890f40e36cb9a $hg004_truth_base/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed 13081043 88d1c926fdc8abd9c39b3e9fe2af3fc43dd67d284690a3dc6127b369d96bccad $hg004_truth_base/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# Slice 3f's typing downloads, one per line: path under $GENOKEEPER_HOME, size, sha256, URL.
# The HHV-6 FASTA files are NCBI efetch's, fetched on 2026-10-30 (NCBI publishes no digest, so
# their pins vouch for that first approved fetch only; the bootstrap concatenates them into
# hhv6/hhv6.fa). hmmix's four hg38 files (Zenodo 11212339, CC BY 4.0) and its 22 archaic BCFs
# with their .csi (Zenodo 13368126, CC BY 4.0) were each checked against Zenodo's md5 and
# hashed on 2026-10-30. Each file waits in prefetch/3f/ under the same layout (hhv6 files at
# the top, then hmmix/ and archaic/), or is downloaded with --with-typing. A file served with
# other bytes is refused.
ncbi_efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'
zenodo_hmmix=https://zenodo.org/api/records/11212339/files
zenodo_archaic=https://zenodo.org/api/records/13368126/files
hhv6_parts="NC_001664.4 NC_000898.1"
typing_downloads="hhv6/NC_001664.4.fa 161748 89d048cf2c599507cc649f76a48beb56fb80549737758cc2f41dd7bded408e3b $ncbi_efetch&id=NC_001664.4&rettype=fasta&retmode=text
hhv6/NC_000898.1.fa 164482 98f40eb9f85220c595c3a023d9c38b1f0060414d0a1fc2986f1759c16d981eef $ncbi_efetch&id=NC_000898.1&rettype=fasta&retmode=text
hmmix/0.9.2/hg38_Outgroup_1000g_HGDP.txt 1587317664 ea2056d3d8c06b7560de35c07b19f3d5655c27c5fc191198f83394d5a57a0e31 $zenodo_hmmix/hg38_Outgroup_1000g_HGDP.txt/content
hmmix/0.9.2/hg38_ancestral.tar.gz 870051319 56a62f96c8061983417c503b8223ae1ad833f6bbfaa70d0dee7437f3b86d288d $zenodo_hmmix/hg38_ancestral.tar.gz/content
hmmix/0.9.2/hg38_strick_callability_mask.bed 98186619 99db6c80c00d0892d86f727d7eaedd197c3754b55f05c1da94f63e4755506834 $zenodo_hmmix/hg38_strick_callability_mask.bed/content
hmmix/0.9.2/hg38_mutationrate.bed 87107 0e371655e59ab8fff0be31b20d56ff7be043aada213c73c41c800617ead9851a $zenodo_hmmix/hg38_mutationrate.bed/content
hmmix/0.9.2/archaic/individuals_highcov.1.bcf 73785350 87ad9e641c6215a734bec58bc55106a477651fc5683ffaf9f6a6c73b143ee296 $zenodo_archaic/individuals_highcov.1.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.1.bcf.csi 134834 162e06e8cae08d589eea22be9ed7050ce14a069e0efbd8d06d360641b80e3e2d $zenodo_archaic/individuals_highcov.1.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.2.bcf 83256012 ae5a91556756cc467de2e112a5446ad35656f3b331bcbe25751cd14eb80986c6 $zenodo_archaic/individuals_highcov.2.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.2.bcf.csi 144911 e6a30f40951732ac16f1dc2baca0c87348e6e8e043695c0b54274a83f6d96fa5 $zenodo_archaic/individuals_highcov.2.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.3.bcf 70265195 ce2941a19fbaab3f4a9858a0b22b7481b0a4ff1ace9069357a4f89d5d5303eea $zenodo_archaic/individuals_highcov.3.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.3.bcf.csi 121531 2a483eb4814665e11f1b77e1c4a09dd532b333e4428c6fa75ed3d53eaf763572 $zenodo_archaic/individuals_highcov.3.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.4.bcf 68631619 9c290aace5a2a7582561f45ab3e6712cefaeb8e555eaae2df0c3699f2431844a $zenodo_archaic/individuals_highcov.4.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.4.bcf.csi 116833 a956504e87107d501c02aa16632563260e2682b10f374b25a648c61219893e52 $zenodo_archaic/individuals_highcov.4.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.5.bcf 63418060 b8dd8ac768d9e1f8758b0b40439b72e93b1105e1d2cd0741c4f59c71e7d495e9 $zenodo_archaic/individuals_highcov.5.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.5.bcf.csi 109112 53b877637ae7b955941dc9e14b6082883b9092c9dfb9f225cf387baeb7813f14 $zenodo_archaic/individuals_highcov.5.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.6.bcf 58273678 9780301cc90ef65036e5969ad1976915aa06922fd433032d491309b39ecd163a $zenodo_archaic/individuals_highcov.6.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.6.bcf.csi 102928 2347539da625cb70111082920e7f55c3c3f38e97552e0f09057cb6a76f06ad76 $zenodo_archaic/individuals_highcov.6.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.7.bcf 51348217 1fad31548086deab8755bff7971bf6ecce6f4eb577b7a2e150a07b25f6d721f0 $zenodo_archaic/individuals_highcov.7.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.7.bcf.csi 93028 d18c1207fe55776b4bb8f244ca8cab35b522ca6288fa1e4241ec2b74a3110660 $zenodo_archaic/individuals_highcov.7.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.8.bcf 52869550 dd0f2c6d8cf540d8a2dae77c82be2438bd5641f0bb7525190e989e4d3550f8ae $zenodo_archaic/individuals_highcov.8.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.8.bcf.csi 87267 e3380d927a8a72f9085d39b57b5638a80ba3c5e94022a19517735a3dd665a5b9 $zenodo_archaic/individuals_highcov.8.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.9.bcf 40138540 070767de84051747cbbfc5e120d417be9e6ef14ea5ad3737af384f9a8ae7058b $zenodo_archaic/individuals_highcov.9.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.9.bcf.csi 67558 a9801399620d95b1ff4a4e32ae5386cc239356d29dcd51cccc8f6a0280a8ddd4 $zenodo_archaic/individuals_highcov.9.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.10.bcf 47157094 0def4f52ea413e6ebec9e4b0ac857e925da9ef2e148ae923262fc69c240a4fe8 $zenodo_archaic/individuals_highcov.10.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.10.bcf.csi 78940 13a496db9f34656bfc22f70e90c31fa9940fe7d49d5526e0602e3986af67a821 $zenodo_archaic/individuals_highcov.10.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.11.bcf 45086009 346a6785b97ffd28c088fa26153368ed5465b8776c89e6afebc4197c7487abf8 $zenodo_archaic/individuals_highcov.11.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.11.bcf.csi 79230 fa7d5d7e114042175a6232b7b5d88feadb7f8649b14de5d88b9a488d1c1e23bf $zenodo_archaic/individuals_highcov.11.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.12.bcf 44896875 f7f3252ec6d5c981ba8e60c32020a9832c827e52b30df7845a37cfc2cbf928d3 $zenodo_archaic/individuals_highcov.12.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.12.bcf.csi 79606 821c47acb7f6aa4be838be6531dc23bc2067b690903ad9c28ca84818e8634259 $zenodo_archaic/individuals_highcov.12.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.13.bcf 36013962 b313f760083ffa3d6f17754fad4146b6eb14e2a251745659e3b54bc47a876e8f $zenodo_archaic/individuals_highcov.13.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.13.bcf.csi 58579 5f150f7d13485cd43a975fb1b929b89c438dd949eab839f74a39d5d6f7d61eae $zenodo_archaic/individuals_highcov.13.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.14.bcf 29950306 82f41d834957f98e443c2fbdbaee1370fdaa652308c9f9e481e79976b188da07 $zenodo_archaic/individuals_highcov.14.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.14.bcf.csi 53877 d163099e524392b826530c87ae09d32ffdf02499bf62b57cb9bc1b4d485f224d $zenodo_archaic/individuals_highcov.14.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.15.bcf 27520613 f6d13a5457ce8f7e1b48c5e315813dbd759756f7925552a0c1cc0d7c5766be05 $zenodo_archaic/individuals_highcov.15.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.15.bcf.csi 47585 9fc8a90027a758441990ee8625f6a3d29addb74f6da6eb987562f893b5e01048 $zenodo_archaic/individuals_highcov.15.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.16.bcf 28771938 df380b1bd142ce2cb33d7f8e4f6a1338d493a96402ca6a5872b40c7157ead871 $zenodo_archaic/individuals_highcov.16.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.16.bcf.csi 45653 00342c1840da8e2deef1f044bd5972759da497f7d46116d3c39863af5f5a4075 $zenodo_archaic/individuals_highcov.16.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.17.bcf 24771514 c7f32fbfc5427f8403674e19a8db331597b8e5969de8b195a63a14aba57add31 $zenodo_archaic/individuals_highcov.17.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.17.bcf.csi 46087 b14b21542294652ab20dabe044de5eab04c659cac542e86663884ec7d4c4da92 $zenodo_archaic/individuals_highcov.17.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.18.bcf 28649110 4bbaf1b9318908c3dad368816bad2fdec7a8d934bfa6a963f12a66b3e0f856de $zenodo_archaic/individuals_highcov.18.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.18.bcf.csi 45396 ef329005855b1dba968dfc789f9ed39a0a57460b39af275383854eafd010f3e2 $zenodo_archaic/individuals_highcov.18.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.19.bcf 16042692 24ea6440390f7bc61fa8bea1c1e99eab630e2f67b3ac2fc6087c1ae12b18e52c $zenodo_archaic/individuals_highcov.19.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.19.bcf.csi 32121 bf05bddc9d43cdd3a4f4aecc1de3529a26df8e8d1124260dc516f559dce54340 $zenodo_archaic/individuals_highcov.19.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.20.bcf 22594270 95fc1aad9f9eaa52b5e414609df196aca8d017aa8bb32f650702a975b9a8db39 $zenodo_archaic/individuals_highcov.20.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.20.bcf.csi 36252 501b67787872559929980a6a525240307aa79992ecfc988a2fa4dc608c6d8752 $zenodo_archaic/individuals_highcov.20.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.21.bcf 13211088 6163756f4021d4814f423ea7e30eaee4d0d63fc0e57635b3706831467a6e5c0c $zenodo_archaic/individuals_highcov.21.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.21.bcf.csi 20650 7e5a1366e6a013c22c95a461194dff3792d8e9dca1d82d95ae1e8a636318a9d6 $zenodo_archaic/individuals_highcov.21.bcf.csi/content
hmmix/0.9.2/archaic/individuals_highcov.22.bcf 12082465 7ce1662bdb3f6885f83ae825d40072efacd749f8a230ddd36e30c225f3f3d506 $zenodo_archaic/individuals_highcov.22.bcf/content
hmmix/0.9.2/archaic/individuals_highcov.22.bcf.csi 20064 08d874eabeb8b98f6c35348be8cf9df3a0ff850b723dc55c1b36a9ebbd72868b $zenodo_archaic/individuals_highcov.22.bcf.csi/content"

die() { echo "bootstrap: $*" >&2; exit 1; }
# Progress goes to stderr: stdout carries --check's list and function results.
say() { echo "bootstrap: $*" >&2; }
missing() { echo "missing: $*"; }

sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
md5of() {
  if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}

# recorded STEP PATH: the manifest has a row for STEP and PATH, and the file has its sha256.
recorded() {
  local h
  [ -f "$manifest" ] && [ -f "$2" ] || return 1
  h=$(awk -F'\t' -v s="$1" -v p="$2" '$1 == s && $2 == p { print $3 }' "$manifest")
  [ -n "$h" ] && [ "$(sha "$2")" = "$h" ]
}

# step_done STEP [PATH...]: the manifest has rows for STEP, each row's file has its recorded
# sha256, and every PATH given is among those rows.
step_done() {
  local step=$1 rows p
  shift
  [ -z "$refresh" ] || return 1 # --refresh fetches whatever the manifest says
  [ -f "$manifest" ] || return 1
  rows=$(awk -F'\t' -v s="$step" '$1 == s { print $2 }' "$manifest")
  [ -n "$rows" ] || return 1
  while IFS= read -r p; do
    recorded "$step" "$p" || return 1
  done <<< "$rows"
  # A here-string, not `printf | grep -q`: grep exits at its first match, and a printf still
  # writing then dies of SIGPIPE, which pipefail turns into "not done" (Linux's bash writes
  # each line apart; the 45 scoring files hit it).
  for p in "$@"; do
    grep -qxF -- "$p" <<< "$rows" || return 1
  done
}

# file_day FILE: the UTC day of FILE's modification time, which `curl -R` set to the server's
# Last-Modified: the name of the folder it is filed under, in the check's own scheme.
file_day() { date -u -r "$1" +%Y%m%d; }

# refreshed STEP FOLDER: under --refresh, the one stdout line `genokeeper refresh` reads.
refreshed() {
  if [ -n "$refresh" ]; then echo "refreshed: $1 $2"; fi
}

# record STEP PATH...: replace STEP's manifest rows with these files and their sha256.
record() {
  local step=$1 tmp p
  shift
  mkdir -p "$lab"
  tmp=$(mktemp "$lab/.manifest.XXXXXX")
  if [ -f "$manifest" ]; then
    awk -F'\t' -v s="$step" '$1 != s' "$manifest" > "$tmp"
  fi
  for p in "$@"; do
    printf '%s\t%s\t%s\n' "$step" "$p" "$(sha "$p")" >> "$tmp"
  done
  mv "$tmp" "$manifest"
}

# download CURL_ARG...: curl with the bootstrap's retries and bounds, the one place it runs:
# a connection not made in a minute, or a transfer under 1 KB/s for 5 minutes, fails and is
# retried, so a stalled server cannot block a bootstrap or a refresh.
download() {
  curl -fL --retry 5 --retry-delay 10 --connect-timeout 60 --speed-limit 1024 --speed-time 300 "$@"
}

# fetch URL DEST: download through DEST.part, resuming a partial file.
fetch() {
  say "downloading $1"
  download -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}

# genokeeper_run COMMAND...: COMMAND (python, genokeeper) in the project's environment.
genokeeper_run() {
  if [ "$check" = 1 ]; then
    # `uv run --no-sync` would still create an empty .venv.
    [ -x "$repo/.venv/bin/python" ] || return 1
    uv run --project "$repo" --no-sync --quiet "$@"
  else
    uv run --project "$repo" --frozen --quiet "$@"
  fi
}
genokeeper_python() { genokeeper_run python "$@"; }

step_tools() {
  local tool formula formulas=""
  for tool in bcftools samtools bgzip tabix; do
    command -v "$tool" >/dev/null && continue
    case "$tool" in
      bcftools | samtools) formula=$tool ;;
      *) formula=htslib ;;
    esac
    if [ "$check" = 1 ]; then
      missing "tool $tool (brew install $formula)"
    else
      case " $formulas " in *" $formula "*) ;; *) formulas="$formulas $formula" ;; esac
    fi
  done
  if [ -n "$formulas" ]; then
    # shellcheck disable=SC2086 # one word per formula
    HOMEBREW_NO_ANALYTICS=1 brew install $formulas
  fi
  if ! command -v uv >/dev/null; then
    [ "$check" = 1 ] || die "uv is not on PATH"
    missing "tool uv"
  fi
}

# The GRCh38 weekly VCF of the current release: its .md5 names the dated file the weekly link
# points to, and its ##fileDate files it under clinvar/<YYYYMMDD>/.
step_clinvar() {
  local work base md5 path name date dest
  step_done clinvar && return
  if [ "$check" = 1 ]; then
    missing "ClinVar GRCh38 weekly VCF and .tbi ($clinvar_base/vcf_GRCh38/weekly/)"
    return
  fi
  work=$lab/clinvar/.download
  # A refresh never resumes onto a part a failed run left: the release may have moved on.
  if [ -n "$refresh" ]; then rm -rf "$work"; fi
  mkdir -p "$work"
  base=$clinvar_base/vcf_GRCh38/weekly
  rm -f "$work/GRCh38.md5.part" # the weekly .md5 is another file each week
  fetch "$base/clinvar.vcf.gz.md5" "$work/GRCh38.md5"
  read -r md5 path < "$work/GRCh38.md5" || die "empty ClinVar GRCh38 .md5"
  name=$(basename "$path")
  case "$name" in
    clinvar_[0-9]*.vcf.gz) ;;
    *) die "unexpected name in the ClinVar GRCh38 .md5: $name" ;;
  esac
  # Under the release's own name: a part cut in an earlier week is never resumed onto it.
  fetch "$base/$name" "$work/$name"
  fetch "$base/$name.tbi" "$work/$name.tbi"
  [ "$(md5of "$work/$name")" = "$md5" ] || die "md5 mismatch for $base/$name"
  tabix -l "$work/$name" > /dev/null || die "unreadable index $base/$name.tbi"
  date=$(bcftools view -h "$work/$name" | sed -n 's/^##fileDate=//p' | sed 's/-//g')
  [[ "$date" =~ ^[0-9]{8}$ ]] || die "no ##fileDate in $base/$name"
  # The run reads this file with genokeeper.panel.build's reader: a release that reader
  # rejects (a record without its ID, no classified record) must never become the newest.
  genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.panel.build import check_clinvar38
check_clinvar38(Path(sys.argv[1]))' "$work/$name" 2> /dev/null || {
    rm -rf "$work"
    die "the ClinVar release $name is not one genokeeper's reader takes"
  }
  dest=$lab/clinvar/$date
  mkdir -p "$dest"
  mv "$work/$name" "$dest/clinvar_GRCh38.vcf.gz"
  mv "$work/$name.tbi" "$dest/clinvar_GRCh38.vcf.gz.tbi"
  rm -rf "$work"
  record clinvar "$dest/clinvar_GRCh38.vcf.gz" "$dest/clinvar_GRCh38.vcf.gz.tbi"
  say "ClinVar $date ($name) -> $dest"
  refreshed clinvar "$dest"
}

# unchanged_copy FOLDER FILE NEW: the newest dated folder under FOLDER, when its FILE holds
# the same curations as NEW. ClinGen writes the day of the build into both of its files, so
# a download equal but for that line is not filed again; GenCC's export, which has no such
# line, only when every line is equal.
unchanged_copy() {
  local d newest=
  for d in "$1"/[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]; do
    [ ! -f "$d/$2" ] || newest=$d
  done
  [ -n "$newest" ] || return 1
  genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.clingen import same_curations
sys.exit(0 if same_curations(Path(sys.argv[1]), Path(sys.argv[2])) else 1)' "$newest/$2" "$3" || return 1
  echo "$newest"
}

# The export is generated on request, dated by its own "FILE CREATED" line (the day of the
# request), and filed only when its curations changed; a partial download is never resumed,
# since the next request may serve another day's file.
step_clingen() {
  local work date dest
  step_done clingen && return
  if [ "$check" = 1 ]; then
    missing "ClinGen gene-disease validity export ($clingen_url)"
    return
  fi
  work=$lab/clingen/.download
  mkdir -p "$work"
  rm -f "$work/$clingen_file.part"
  fetch "$clingen_url" "$work/$clingen_file"
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.clingen import gene_validity_date, load_gene_validity
path = Path(sys.argv[1])
load_gene_validity(path)
print(gene_validity_date(path))' "$work/$clingen_file") || {
    rm -rf "$work"
    die "the ClinGen gene validity export from $clingen_url does not parse"
  }
  [[ "$date" =~ ^[0-9]{8}$ ]] || die "no FILE CREATED date in the ClinGen gene validity export"
  if dest=$(unchanged_copy "$lab/clingen" "$clingen_file" "$work/$clingen_file"); then
    rm -rf "$work"
    record clingen "$dest/$clingen_file"
    say "ClinGen gene validity unchanged since $dest"
    refreshed clingen "$dest"
    return
  fi
  dest=$lab/clingen/$date
  mkdir -p "$dest"
  mv "$work/$clingen_file" "$dest/$clingen_file"
  rm -rf "$work"
  record clingen "$dest/$clingen_file"
  say "ClinGen gene validity $date -> $dest/$clingen_file"
  refreshed clingen "$dest"
}

step_grch38_fasta() {
  local dest work md5
  dest=$(unmasked_fasta_path) || {
    [ "$check" = 1 ] && { missing "project environment for the GRCh38 FASTA path (uv sync)"; return; }
    die "cannot read the GRCh38 FASTA path from genokeeper.config"
  }
  step_done grch38_fasta "$dest" "$dest.fai" && return
  if [ "$check" = 1 ]; then
    missing "GRCh38 no-alt FASTA $dest with .fai (from $grch38_base/$grch38_name)"
    return
  fi
  work=$(dirname "$dest")/.grch38
  mkdir -p "$work"
  fetch "$grch38_base/md5checksums.txt" "$work/md5checksums.txt"
  md5=$(awk -v n="./$grch38_name" '$2 == n { print $1 }' "$work/md5checksums.txt")
  [ -n "$md5" ] || die "$grch38_name is not in $grch38_base/md5checksums.txt"
  [ -f "$work/$grch38_name" ] || fetch "$grch38_base/$grch38_name" "$work/$grch38_name"
  if [ "$(md5of "$work/$grch38_name")" != "$md5" ]; then
    rm -f "$work/$grch38_name"
    die "md5 mismatch for $grch38_base/$grch38_name"
  fi
  gzip -dc "$work/$grch38_name" > "$dest.tmp"
  mv "$dest.tmp" "$dest"
  samtools faidx "$dest"
  rm -rf "$work"
  record grch38_fasta "$dest" "$dest.fai"
  say "GRCh38 no-alt FASTA -> $dest"
}

# tool_version TOOL: the version TOOL reports, or nothing. bwa prints it to stderr and exits
# 1 without a command.
tool_version() {
  case "$1" in
    bwa) { bwa 2>&1 || true; } | sed -n 's/^Version: \([0-9][0-9.]*\).*/\1/p' ;;
    mosdepth) { mosdepth --version 2>&1 || true; } | awk '$1 == "mosdepth" { print $2 }' ;;
  esac
}
pinned_version() {
  case "$1" in
    bwa) echo "$bwa_version" ;;
    mosdepth) echo "$mosdepth_version" ;;
  esac
}

rtg_dir() { echo "$lab/tools/rtg-tools-$rtg_version"; }
rtg_ok() {
  local dir
  dir=$(rtg_dir)
  [ -x "$dir/rtg" ] || return 1
  # awk reads the whole output: `grep -q` could stop early and fail the pipe with SIGPIPE.
  "$dir/rtg" version 2>/dev/null |
    awk -v v="$rtg_version" '/RTG Tools/ { for (i = 1; i <= NF; i++) if ($i == v) f = 1 }
                             END { exit !f }'
}

# warn_version TOOL: a warning when TOOL is not the version tested here. Not a refusal: no
# adversary hides behind a version (the user after `brew upgrade` is the one it would stop),
# and every run records the version and keys its align and QC caches on it.
warn_version() {
  local v
  v=$(tool_version "$1")
  [ "$v" = "$(pinned_version "$1")" ] ||
    say "warning: $1 ${v:-of an unknown version} is not the tested $(pinned_version "$1")"
}

# bwa and mosdepth from Homebrew, rtg-tools from its release zip; another bwa or mosdepth
# version is a warning (warn_version), never upgraded or replaced.
step_tools38() {
  local tool formulas=""
  for tool in bwa mosdepth; do
    if command -v "$tool" >/dev/null; then
      warn_version "$tool"
    elif [ "$check" = 1 ]; then
      missing "tool $tool $(pinned_version "$tool") (brew install $tool)"
    else
      formulas="$formulas $tool"
    fi
  done
  if [ -n "$formulas" ]; then
    # shellcheck disable=SC2086 # one word per formula
    HOMEBREW_NO_ANALYTICS=1 brew install $formulas
    for tool in $formulas; do warn_version "$tool"; done
  fi
  step_rtg
}

# No manifest row: `rtg version` is the check. rtg.cfg is written before the first run, so
# rtg never prompts, and it turns off crash reports and usage logging (both send data to Real
# Time Genomics; rtg-tools reads the calibration's calls in stage 2b).
step_rtg() {
  local dir work java
  rtg_ok && return
  dir=$(rtg_dir)
  if [ "$check" = 1 ]; then
    missing "rtg-tools $rtg_version in $dir (Homebrew's openjdk and $rtg_url)"
    return
  fi
  java=$(HOMEBREW_NO_ANALYTICS=1 brew --prefix openjdk 2>/dev/null)/bin/java
  if [ ! -x "$java" ]; then
    HOMEBREW_NO_ANALYTICS=1 brew install openjdk
    java=$(HOMEBREW_NO_ANALYTICS=1 brew --prefix openjdk)/bin/java
  fi
  work=$lab/tools/.download
  rm -rf "$work"
  mkdir -p "$work"
  fetch "$rtg_url" "$work/rtg.zip"
  [ "$(sha "$work/rtg.zip")" = "$rtg_sha256" ] || {
    rm -rf "$work"
    die "sha256 mismatch for $rtg_url"
  }
  unzip -q "$work/rtg.zip" -d "$work"
  printf 'RTG_JAVA="%s"\nRTG_TALKBACK=false\nRTG_USAGE=\n' "$java" \
    > "$work/rtg-tools-$rtg_version/rtg.cfg"
  rm -rf "$dir"
  mv "$work/rtg-tools-$rtg_version" "$dir"
  rm -rf "$work"
  rtg_ok || die "rtg-tools in $dir does not report version $rtg_version"
  say "rtg-tools $rtg_version -> $dir"
}

# The GIAB v3.6 GRCh38 BEDs as published (already bgzipped, chr names), tabix-indexed, and
# the region BEDs genokeeper.regions derives from them.
step_giab38() {
  local dir names derived name rel md5 outputs listed=0
  outputs=()
  names=$(genokeeper_python -m genokeeper.regions downloads) &&
    derived=$(genokeeper_python -m genokeeper.regions derived) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.regions (uv sync)"; return; }
    die "cannot list the GIAB GRCh38 files from genokeeper.regions"
  }
  dir=$(genokeeper_python -c '
from genokeeper.config import load_paths
from genokeeper.regions import giab38_dir
print(giab38_dir(load_paths().lab))') || die "cannot read the GIAB GRCh38 path"
  for name in $names; do
    outputs+=("$dir/$name" "$dir/$name.tbi")
  done
  for name in $derived; do
    outputs+=("$dir/$name")
  done
  step_done giab38 "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    for name in $names; do
      if ! recorded giab38 "$dir/$name" || ! recorded giab38 "$dir/$name.tbi"; then
        missing "GIAB v3.6 GRCh38 stratification $name"
        listed=1
      fi
    done
    for name in $derived; do
      recorded giab38 "$dir/$name" && continue
      missing "region BED $dir/$name (derived)"
      listed=1
    done
    [ "$listed" = 1 ] || missing "GIAB v3.6 GRCh38 stratifications: stale manifest rows"
    return
  fi
  mkdir -p "$dir/.download"
  fetch "$giab38_base/$giab38_md5s" "$dir/.download/$giab38_md5s"
  for name in $names; do
    # A file already recorded is kept: a lab from before a name was added (slice 3c's
    # three chrX BEDs) fetches the new names only. --refresh fetches every one again.
    if [ -z "$refresh" ] && recorded giab38 "$dir/$name" && recorded giab38 "$dir/$name.tbi"; then
      continue
    fi
    read -r md5 rel < <(awk -v n="$name" '{ b = $2; sub(/.*\//, "", b) } b == n { print $1, $2 }' \
      "$dir/.download/$giab38_md5s") || die "$name is not in $giab38_base/$giab38_md5s"
    fetch "$giab38_base/$rel" "$dir/.download/$name"
    [ "$(md5of "$dir/.download/$name")" = "$md5" ] || die "md5 mismatch for $giab38_base/$rel"
    mv "$dir/.download/$name" "$dir/$name"
    tabix -f -p bed "$dir/$name" || die "$giab38_base/$rel is not a bgzipped BED"
  done
  rm -rf "$dir/.download"
  genokeeper_python -m genokeeper.regions derive "$dir" > /dev/null ||
    die "deriving the region BEDs in $dir failed"
  record giab38 "${outputs[@]}"
  say "GIAB v3.6 GRCh38 stratifications and derived regions -> $dir"
}

unmasked_fasta_path() {
  genokeeper_python -c '
from genokeeper.config import load_paths, unmasked_fasta
print(unmasked_fasta(load_paths().lab))'
}

step_dv_model() {
  local dir entry name md5 outputs
  outputs=()
  dir=$lab/deepvariant/complete-t7
  for entry in $dv_model_files; do outputs+=("$dir/${entry%%:*}"); done
  step_done dv_model "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    for entry in $dv_model_files; do
      name=${entry%%:*}
      recorded dv_model "$dir/$name" ||
        missing "DeepVariant T7 checkpoint file $name ($dv_model_base/$name)"
    done
    return
  fi
  mkdir -p "$dir/.download"
  for entry in $dv_model_files; do
    name=${entry%%:*}
    md5=${entry#*:}
    # On disk with the pinned md5 (fetched before): recognised, not fetched again.
    [ -f "$dir/$name" ] && [ "$(md5of "$dir/$name")" = "$md5" ] && continue
    fetch "$dv_model_base/$name" "$dir/.download/$name"
    if [ "$(md5of "$dir/.download/$name")" != "$md5" ]; then
      rm -rf "$dir/.download"
      die "md5 mismatch for $dv_model_base/$name"
    fi
    mv "$dir/.download/$name" "$dir/$name"
  done
  rm -rf "$dir/.download"
  record dv_model "${outputs[@]}"
  say "DeepVariant T7 checkpoint -> $dir"
}

# The Clinical Genomic Database's own download, fetched whole and never resumed (the next
# request may serve a newer file), filed under the server's Last-Modified date.
step_cgd() {
  local work date dest first
  step_done cgd && return
  if [ "$check" = 1 ]; then
    missing "CGD gene file ($cgd_url)"
    return
  fi
  work=$lab/cgd/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $cgd_url"
  download -R -o "$work/$cgd_file" "$cgd_url"
  first=$(gzip -dc "$work/$cgd_file" 2>/dev/null | sed -n 1p) || first=
  case "$first" in
    "#GENE"$'\t'"HGNC ID"$'\t'*) ;;
    *) rm -rf "$work"; die "$cgd_url is not the CGD table (first line: ${first:0:60})" ;;
  esac
  date=$(file_day "$work/$cgd_file")
  dest=$lab/cgd/$date
  mkdir -p "$dest"
  mv "$work/$cgd_file" "$dest/$cgd_file"
  rm -rf "$work"
  record cgd "$dest/$cgd_file"
  say "CGD $date -> $dest/$cgd_file"
}

# PharmCAT's positions file and its index, kept only with the sha256 pinned above: its anchor
# comes from the release digests written in this script, which a re-uploaded or tampered asset
# cannot move because the comparison runs on the bytes on this disk before they are filed.
step_pharmcat() {
  local dir entry name pin outputs
  outputs=()
  dir=$lab/pharmcat/$pharmcat_version
  for entry in $pharmcat_files; do outputs+=("$dir/${entry%%:*}"); done
  step_done pharmcat "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    for entry in $pharmcat_files; do
      name=${entry%%:*}
      recorded pharmcat "$dir/$name" ||
        missing "PharmCAT $pharmcat_version positions file $name ($pharmcat_base/$name)"
    done
    return
  fi
  rm -rf "$dir/.download"
  mkdir -p "$dir/.download"
  for entry in $pharmcat_files; do
    name=${entry%%:*}
    pin=${entry#*:}
    fetch "$pharmcat_base/$name" "$dir/.download/$name"
    if [ "$(sha "$dir/.download/$name")" != "$pin" ]; then
      rm -rf "$dir/.download"
      die "sha256 mismatch for $pharmcat_base/$name"
    fi
  done
  for entry in $pharmcat_files; do
    mv "$dir/.download/${entry%%:*}" "$dir/${entry%%:*}"
  done
  rm -rf "$dir/.download"
  record pharmcat "${outputs[@]}"
  say "PharmCAT $pharmcat_version positions -> $dir"
}

# CPIC's data, fetched whole and never resumed (the next request may serve newer rows),
# filed under the fetch date when it differs from the newest copy. Only rows genokeeper's reader takes are recorded: its anchor
# comes from that reader run on the bytes on this disk, which a maintenance page, an error
# answer or an empty one cannot pass because none is a non-empty list of the rows it reads.
# The step is done once the newest dated folder holds the files cpic_files names.
step_cpic() {
  local work entry name dest names entries newest=
  read -ra entries <<< "$cpic_files"
  for dest in "$lab"/cpic/2*/; do newest=${dest%/}; done
  names=()
  for entry in "${entries[@]}"; do names+=("$newest/${entry%%:*}"); done
  step_done cpic "${names[@]}" && return
  if [ "$check" = 1 ]; then
    missing "CPIC gene-drug pairs and MT-RNR1 alleles ($cpic_api/)"
    return
  fi
  work=$lab/cpic/.download
  rm -rf "$work"
  mkdir -p "$work"
  names=()
  for entry in "${entries[@]}"; do
    name=${entry%%:*}
    say "downloading $cpic_api/${entry#*:}"
    download -o "$work/$name" "$cpic_api/${entry#*:}"
    names+=("$name")
  done
  genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.pgx import CPIC_FILES, load_cpic_alleles, load_cpic_levels
folder = Path(sys.argv[1])
load_cpic_levels(folder / CPIC_FILES["pairs"])
load_cpic_alleles(folder / CPIC_FILES["MT-RNR1"])' "$work" 2> /dev/null || {
    rm -rf "$work"
    die "$cpic_api/ answered with no CPIC rows genokeeper reads"
  }
  # CPIC dates nothing: rows equal to the newest copy's are no new release.
  dest=$newest
  for name in "${names[@]}"; do cmp -s "$work/$name" "$newest/$name" || dest=; done
  if [ -n "$dest" ]; then
    rm -rf "$work"
    record cpic "${names[@]/#/$dest/}"
    say "CPIC data unchanged since $dest"
    refreshed cpic "$dest"
    return
  fi
  dest=$lab/cpic/$(date -u +%Y%m%d)
  mkdir -p "$dest"
  for name in "${names[@]}"; do mv "$work/$name" "$dest/$name"; done
  rm -rf "$work"
  record cpic "${names[@]/#/$dest/}"
  say "CPIC data -> $dest"
  refreshed cpic "$dest"
}

# ClinPGx's annotations, fetched whole and never resumed (the next request may serve
# a newer build), filed under the CREATED date the download itself names. Only the download
# is recorded: its anchor comes from genokeeper's reader run on the bytes on this disk, which a
# maintenance page or another file cannot pass because the reader refuses a file without the
# two annotation tables and the CREATED member.
step_pharmgkb() {
  local work date dest
  step_done pharmgkb && return
  if [ "$check" = 1 ]; then
    missing "ClinPGx clinical annotations ($pharmgkb_url)"
    return
  fi
  work=$lab/pharmgkb/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $pharmgkb_url"
  download -o "$work/$pharmgkb_file" "$pharmgkb_url"
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.pharmgkb import read_clinical_annotations
print(read_clinical_annotations(Path(sys.argv[1])).created)' "$work/$pharmgkb_file" 2> /dev/null) || {
    rm -rf "$work"
    die "$pharmgkb_url is not the ClinPGx clinical annotations download"
  }
  [[ "$date" =~ ^[0-9]{8}$ ]] || { rm -rf "$work"; die "no CREATED date in $pharmgkb_url"; }
  dest=$lab/pharmgkb/$date
  mkdir -p "$dest"
  mv "$work/$pharmgkb_file" "$dest/$pharmgkb_file"
  rm -rf "$work"
  record pharmgkb "$dest/$pharmgkb_file"
  say "ClinPGx clinical annotations $date -> $dest/$pharmgkb_file"
  refreshed pharmgkb "$dest"
}

# MedlinePlus Genetics' summaries, fetched whole and never resumed (MedlinePlus rebuilds the
# file daily, so a resumed part could splice two builds), filed under the fetch date. Only a
# file genokeeper's loader reads is recorded: its anchor comes from that loader run on the
# bytes on this disk, which a maintenance page or another file cannot pass because the
# loader refuses a file whose root is not `summaries` or that holds no gene summary.
step_medlineplus() {
  local work dest
  step_done medlineplus && return
  if [ "$check" = 1 ]; then
    missing "MedlinePlus Genetics summaries ($medlineplus_url)"
    return
  fi
  work=$lab/medlineplus/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $medlineplus_url"
  download -o "$work/$medlineplus_file" "$medlineplus_url"
  genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.medlineplus import load_medlineplus
load_medlineplus(Path(sys.argv[1]))' "$work/$medlineplus_file" 2> /dev/null || {
    rm -rf "$work"
    die "$medlineplus_url is not the MedlinePlus Genetics summaries file"
  }
  dest=$lab/medlineplus/$(date -u +%Y%m%d)
  mkdir -p "$dest"
  mv "$work/$medlineplus_file" "$dest/$medlineplus_file"
  rm -rf "$work"
  record medlineplus "$dest/$medlineplus_file"
  say "MedlinePlus Genetics summaries -> $dest/$medlineplus_file"
  refreshed medlineplus "$dest"
}

# image_pinned REF: the local image REF names lists REF among its digests (a pulled image)
# or its tags (an image built here, tagged with its inputs' sha256).
image_pinned() {
  docker image inspect \
    --format '{{range .RepoDigests}}{{println .}}{{end}}{{range .RepoTags}}{{println .}}{{end}}' \
    "$1" 2>/dev/null | grep -xF -- "$1" > /dev/null
}

# The GWAS Catalog's associations, fetched whole and never resumed (the next request may
# serve a newer release). Its anchor is genokeeper's reader, which refuses a file that is not
# the download before anything is filed, so a cut transfer or an error page cannot pass.
step_gwas() {
  local work date dest
  step_done gwas && return
  if [ "$check" = 1 ]; then
    missing "GWAS Catalog associations ($gwas_url)"
    return
  fi
  work=$lab/gwas/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $gwas_url"
  download -o "$work/$gwas_file" "$gwas_url"
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.gwas import load_associations
print(load_associations(Path(sys.argv[1])).release)' "$work/$gwas_file" 2> /dev/null) || {
    rm -rf "$work"
    die "$gwas_url is not the GWAS Catalog associations download"
  }
  [[ "$date" =~ ^[0-9]{8}$ ]] || { rm -rf "$work"; die "no build date in $gwas_url"; }
  dest=$lab/gwas/$date
  mkdir -p "$dest"
  mv "$work/$gwas_file" "$dest/$gwas_file"
  rm -rf "$work"
  record gwas "$dest/$gwas_file"
  say "GWAS Catalog associations $date -> $dest/$gwas_file"
  refreshed gwas "$dest"
}

# Every image in genokeeper.containers.IMAGES, pulled by its digest: Docker resolves
# `repo@sha256:<digest>` by content and refuses other bytes, so nothing re-checks the pull.
# An image genokeeper.containers.BUILT names (none of its version for its platform) is built
# from its Dockerfile under docker/ instead, once docker/ holds the inputs its name records,
# and the build checks the Dockerfile's own pins. No manifest row: Docker's store is the record.
step_images() {
  local images platform how name ref names
  # name, reference, platform ("-": pinned by a multi-platform index, pulled for the host)
  # and whether it is pulled or built.
  images=$(genokeeper_python -c '
from genokeeper.containers import BUILT, IMAGES, PLATFORMS
for name, ref in IMAGES.items():
    print(name, ref, PLATFORMS.get(name, "-"), "build" if name in BUILT else "pull")') || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.containers (uv sync)"; return; }
    die "cannot read the image table from genokeeper.containers"
  }
  if ! docker info > /dev/null 2>&1; then
    names=$(awk '{ printf "%s%s", sep, $1; sep = ", " }' <<< "$images")
    if [ "$check" = 1 ]; then
      missing "Docker (start Docker Desktop; the images $names cannot be checked)"
      return
    fi
    die "Docker is not running: start Docker Desktop and rerun (images: $names)"
  fi
  while read -r name ref platform how; do
    image_pinned "$ref" && continue
    if [ "$check" = 1 ]; then
      missing "container image $name ($ref)"
      continue
    fi
    if [ "$how" = build ]; then
      genokeeper_python -m genokeeper.containers build "$name" >&2 ||
        die "container image $name: the build from docker/ failed"
    elif [ "$platform" = - ]; then
      docker pull "$ref" >&2
    else
      docker pull --platform "$platform" "$ref" >&2
    fi
  done <<< "$images"
}

# t1k_build DAT OUT_DIR: T1K's HLA index from hla.dat, offline in T1K's pinned image.
t1k_build() { genokeeper_python -m genokeeper.pgx.hla build-index "$1" "$2" >&2; }

# T1K's HLA index, built from IPD-IMGT/HLA's hla.dat.zip, fetched whole (never resumed) and
# kept only with the sha256 pinned above: its anchor comes from that pin, which a re-uploaded
# release cannot move because the comparison runs on the bytes on this disk before the index
# is built from them. After step_images: the build runs in T1K's image.
step_t1k() {
  local dir work index
  dir=$lab/t1k/IPD-IMGT-HLA-$imgt_hla_version
  index=$dir/hla_dna_seq.fa
  step_done t1k "$dir/hla.dat.zip" "$index" && return
  if [ "$check" = 1 ]; then
    missing "T1K HLA index from IPD-IMGT/HLA $imgt_hla_version ($imgt_hla_url)"
    return
  fi
  work=$dir/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $imgt_hla_url"
  download -o "$work/hla.dat.zip" "$imgt_hla_url"
  if [ "$(sha "$work/hla.dat.zip")" != "$imgt_hla_sha256" ]; then
    rm -rf "$work"
    die "sha256 mismatch for $imgt_hla_url"
  fi
  unzip -q -o "$work/hla.dat.zip" hla.dat -d "$work" 2> /dev/null || {
    rm -rf "$work"
    die "$imgt_hla_url holds no hla.dat"
  }
  t1k_build "$work/hla.dat" "$work/index" || {
    rm -rf "$work"
    die "building T1K's HLA index failed"
  }
  mv "$work/hla.dat.zip" "$dir/hla.dat.zip"
  mv "$work/index/"* "$dir/"
  rm -rf "$work"
  record t1k "$dir/hla.dat.zip" "$index"
  say "T1K HLA index (IPD-IMGT/HLA $imgt_hla_version) -> $index"
}

# Nextflow, its nf-schema plugin and pgsc_calc, fetched whole (never resumed) and kept only
# with the sha256 pinned above: their anchor comes from those pins, which a re-uploaded asset
# or a moved tag cannot pass because the comparison runs on the bytes on this disk before
# anything is unpacked. The plugin is unpacked under pgs/nxf-home/plugins/, where Nextflow
# finds it offline.
step_pgs_tools() {
  local dir work java entry url pin nextflow plugin pipeline
  dir=$lab/pgs
  nextflow=$dir/nextflow/nextflow-$nextflow_version-dist
  plugin=$dir/nf-schema-$nf_schema_version.zip
  pipeline=$dir/pgsc_calc-$pgsc_calc_commit.tar.gz
  step_done pgs_tools "$nextflow" "$plugin" "$pipeline" && [ -x "$dir/java/bin/java" ] && return
  if [ "$check" = 1 ]; then
    missing "Nextflow $nextflow_version, nf-schema $nf_schema_version, pgsc_calc $pgsc_calc_version ($pgsc_calc_commit) and $nextflow_java"
    return
  fi
  java=$(HOMEBREW_NO_ANALYTICS=1 brew --prefix "$nextflow_java")
  [ -x "$java/bin/java" ] || HOMEBREW_NO_ANALYTICS=1 brew install "$nextflow_java"
  work=$dir/.download
  rm -rf "$work"
  mkdir -p "$work" "$dir/nextflow" "$dir/pgsc_calc"
  ln -sfn "$java" "$dir/java"
  for entry in "$nextflow_url $nextflow_sha256" "$nf_schema_url $nf_schema_sha256" \
    "$pgsc_calc_url $pgsc_calc_sha256"; do
    read -r url pin <<< "$entry"
    say "downloading $url"
    download -o "$work/${url##*/}" "$url"
    [ "$(sha "$work/${url##*/}")" = "$pin" ] || { rm -rf "$work"; die "sha256 mismatch for $url"; }
  done
  rm -rf "$dir/nxf-home/plugins/nf-schema-$nf_schema_version" "$dir/pgsc_calc/$pgsc_calc_commit"
  mkdir -p "$dir/nxf-home/plugins/nf-schema-$nf_schema_version" "$dir/pgsc_calc/$pgsc_calc_commit"
  # -o: the release zip lists classes/META-INF/ twice, which would make unzip ask.
  unzip -q -o "$work/${nf_schema_url##*/}" -d "$dir/nxf-home/plugins/nf-schema-$nf_schema_version"
  tar -xzf "$work/${pgsc_calc_url##*/}" -C "$dir/pgsc_calc/$pgsc_calc_commit" --strip-components 1
  mv "$work/${nextflow_url##*/}" "$nextflow"
  chmod +x "$nextflow"
  mv "$work/${nf_schema_url##*/}" "$plugin"
  mv "$work/${pgsc_calc_url##*/}" "$pipeline"
  rm -rf "$work"
  record pgs_tools "$nextflow" "$plugin" "$pipeline"
  say "Nextflow $nextflow_version, nf-schema $nf_schema_version, pgsc_calc $pgsc_calc_version -> $dir"
}

# StellarPGx, CYP2D6's third caller, fetched whole and unpacked only with the sha256 pinned above.
step_stellarpgx() {
  local dir archive
  dir=$lab/stellarpgx/$stellarpgx_commit
  step_done stellarpgx "$dir/main.nf" && return
  if [ "$check" = 1 ]; then
    missing "StellarPGx ($stellarpgx_commit)"
    return
  fi
  mkdir -p "$lab/stellarpgx"
  archive=$lab/stellarpgx/.download.tar.gz
  say "downloading $stellarpgx_url"
  download -o "$archive" "$stellarpgx_url"
  [ "$(sha "$archive")" = "$stellarpgx_sha256" ] || { rm -f "$archive"; die "sha256 mismatch for $stellarpgx_url"; }
  rm -rf "$dir"
  mkdir -p "$dir"
  tar -xzf "$archive" -C "$dir" --strip-components 1
  rm -f "$archive"
  record stellarpgx "$dir/main.nf"
  say "StellarPGx -> $dir"
}

# The chosen scores' scoring files, fetched whole (never resumed) and kept only with the md5
# the PGS Catalog publishes beside each, fetched with it: their anchor comes from that md5,
# which a cut or garbled transfer cannot pass because the comparison runs on the bytes on
# this disk before they are filed. It is no defence against EBI's own server serving other
# bytes (the .md5 comes from it too); a re-harmonised file is filed, and the run records
# each file's sha256 (run.inputs.prs.scores), so the change is never silent.
# pgs_fetch_scores STEP WHAT IDS: those scoring files, under their own step in the manifest.
pgs_fetch_scores() {
  local step=$1 what=$2 ids=$3 dir id name url outputs=()
  dir=$lab/pgs/scores
  for id in $ids; do outputs+=("$dir/${id}_hmPOS_GRCh38.txt.gz"); done
  step_done "$step" "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    missing "${#outputs[@]} $what ($pgs_scores_base/)"
    return
  fi
  rm -rf "$dir/.download"
  mkdir -p "$dir/.download"
  for id in $ids; do
    name=${id}_hmPOS_GRCh38.txt.gz
    url=$pgs_scores_base/$id/ScoringFiles/Harmonized/$name
    say "downloading $url"
    download -o "$dir/.download/$name" "$url"
    download -o "$dir/.download/$name.md5" "$url.md5"
    [ "$(md5of "$dir/.download/$name")" = "$(awk '{ print $1; exit }' "$dir/.download/$name.md5")" ] || {
      rm -rf "$dir/.download"
      die "md5 mismatch for $url against $url.md5"
    }
  done
  for id in $ids; do
    mv "$dir/.download/${id}_hmPOS_GRCh38.txt.gz" "$dir/${id}_hmPOS_GRCh38.txt.gz"
  done
  rm -rf "$dir/.download"
  record "$step" "${outputs[@]}"
  say "${#outputs[@]} $what -> $dir"
}

step_pgs_scores() {
  pgs_fetch_scores pgs_scores "PGS Catalog scoring files" "$pgs_ids"
}

step_pgs_research_scores() {
  pgs_fetch_scores pgs_research_scores "PGS Catalog research scoring files" "$pgs_research_ids"
}

# pgs_ancestry_build PANEL OUT_DIR: the panel's ancestry positions, by plink2 offline in
# pgsc_calc's pinned image.
pgs_ancestry_build() {
  genokeeper_python -m genokeeper.prs.pgsc ancestry-positions "$1" "$2" >&2
}

# The HGDP+1kGP panel under its final name only once its md5 is the published one (its
# anchor: the md5 pinned above, which a cut or substituted download cannot pass because the
# comparison runs on the bytes on this disk before the rename); then its ancestry positions.
step_pgs_panel() {
  local dir panel prefetch pvar
  dir=$lab/pgs
  panel=$dir/$pgs_panel_name
  prefetch=$lab/prefetch/1d/pgsc_calc/$pgs_panel_name
  pvar=$dir/ancestry/GRCh38_HGDP+1kGP_ancestry.pvar.zst
  full=$dir/ancestry/$pgs_panel_pvar
  [ -f "$panel" ] && step_done pgs_ancestry "$pvar" "$full" && return
  if [ "$check" = 1 ]; then
    [ -f "$panel" ] || missing "HGDP+1kGP reference panel $panel ($pgs_panel_url, 16 GB)"
    missing "ancestry positions $pvar and the panel's variants $full (built from the panel)"
    return
  fi
  mkdir -p "$dir"
  if [ ! -f "$panel" ]; then
    if [ -f "$prefetch" ] && [ "$(md5of "$prefetch")" = "$pgs_panel_md5" ]; then
      mv "$prefetch" "$panel"
      say "HGDP+1kGP panel moved from $prefetch"
    else
      fetch_md5 "$pgs_panel_url" "$panel" "$pgs_panel_md5"
    fi
  fi
  pgs_ancestry_build "$panel" "$dir/ancestry" || die "building the panel's ancestry positions failed"
  record pgs_ancestry "$pvar" "$full"
  say "HGDP+1kGP panel and its ancestry positions -> $dir"
}

step_pgs() {
  step_pgs_tools
  step_stellarpgx
  step_pgs_scores
  step_pgs_research_scores
  step_pgs_panel
}

# No manifest row: `genokeeper build-faf --check` reads the table's own manifest.
step_gnomad_faf() {
  local stale chrom
  stale=$(genokeeper_run genokeeper build-faf --check) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper build-faf (uv sync)"; return; }
    die "genokeeper build-faf --check failed"
  }
  if [ "$check" = 1 ]; then
    for chrom in $stale; do missing "gnomAD FAF table $chrom (run with --with-gnomad)"; done
    return
  fi
  if [ "$with_gnomad" = 1 ]; then
    # It also asks the bucket for each file's md5: a changed chromosome is rebuilt.
    genokeeper_run genokeeper build-faf
    say "gnomAD FAF table built"
  elif [ -n "$stale" ]; then
    say "gnomAD FAF table missing: run with --with-gnomad (streams 0.88 TB, keeps a few GB)"
  fi
}

# fetch_md5 URL DEST MD5: download through DEST.part (resuming one a cut transfer left) and
# rename it to DEST only when its md5 is MD5; a mismatch removes it and refuses. DEST under
# its final name is therefore a verified file: its anchor comes from the md5 written in this
# script, which a server serving other bytes, or a transfer resumed onto a changed file,
# cannot move because the comparison runs on the bytes on this disk.
fetch_md5() {
  say "downloading $1"
  download -C - -o "$2.part" "$1"
  if [ "$(md5of "$2.part")" != "$3" ]; then
    rm -f "$2.part"
    die "md5 mismatch for $1; removed, rerun to fetch it again"
  fi
  mv "$2.part" "$2"
}

# fetch_truth STEP BASE DIR FILES: each name:md5 of FILES from BASE into DIR through
# fetch_md5 (one on disk with its pinned md5, fetched before, is kept), recorded as STEP.
fetch_truth() {
  local entry name md5 paths=()
  mkdir -p "$3/.download"
  for entry in $4; do
    name=${entry%%:*}
    md5=${entry#*:}
    paths+=("$3/$name")
    [ -f "$3/$name" ] && [ "$(md5of "$3/$name")" = "$md5" ] && continue
    fetch_md5 "$2/$name" "$3/.download/$name" "$md5"
    mv "$3/.download/$name" "$3/$name"
  done
  rm -rf "$3/.download"
  record "$1" "${paths[@]}"
}

# The FASTA's sha256 as the manifest records it, whichever step wrote the row, or nothing.
fasta_recorded() {
  [ -f "$manifest" ] || return 0
  awk -F'\t' -v p="$1" '$2 == p { print $3; exit }' "$manifest"
}

# Stage 2b's calibration inputs (see the header); only with --with-calibration.
step_hg003() {
  local layout fastq_dir ds_dir truth_dir entry name md5 size
  local wanted="" fastq_bytes=0 truth_ok=1 outputs=()
  layout=$(genokeeper_python -m genokeeper.calibrate layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.calibrate (uv sync)"; return; }
    die "cannot read the calibration paths from genokeeper.calibrate"
  }
  fastq_dir=$(awk -F'\t' '$1 == "fastq" { print $2 }' <<< "$layout")
  ds_dir=$(awk -F'\t' '$1 == "downsampled" { print $2 }' <<< "$layout")
  truth_dir=$(awk -F'\t' '$1 == "truth" { print $2 }' <<< "$layout")
  for entry in $hg003_truth_files; do outputs+=("$truth_dir/${entry%%:*}"); done
  if ! step_done hg003_truth "${outputs[@]}"; then
    truth_ok=0
    for entry in $hg003_truth_files; do
      name=${entry%%:*}
      recorded hg003_truth "$truth_dir/$name" ||
        wanted+="GIAB HG003 v4.2.1 truth file $name ($hg003_truth_base/$name)"$'\n'
    done
  fi
  # Once downsampled, the full FASTQ (`genokeeper calibrate` deletes it then) is not wanted again;
  # nor once HG003 has a record (its aligned.json beside the FASTQ folder) or a CRAM under
  # align/, which `python -m genokeeper.calibrate.hg002 realign hg003` re-aligns from.
  local crams=("${fastq_dir%/*}"/align/*/sample.cram)
  if [ ! -f "${fastq_dir%/*}/aligned.json" ] && [ ! -f "${crams[0]}" ]; then
    for entry in $hg003_fastq_files; do
      [ -f "$ds_dir/${entry%%:*}" ] || { ds_dir=; break; }
    done
  fi
  if [ -z "$ds_dir" ]; then
    for entry in $hg003_fastq_files; do
      IFS=: read -r name md5 size <<< "$entry"
      [ -f "$fastq_dir/$name" ] && [ "$(wc -c < "$fastq_dir/$name" | tr -d ' ')" = "$size" ] &&
        continue
      fastq_bytes=$((fastq_bytes + size))
      wanted+="HG003 FASTQ $fastq_dir/$name ($hg003_fastq_base/$name, $size B)"$'\n'
    done
  fi
  [ -n "$wanted" ] || return 0
  if [ "$check" = 1 ]; then
    while IFS= read -r entry; do [ -z "$entry" ] || missing "$entry"; done <<< "$wanted"
    return
  fi
  if [ "$with_calibration" != 1 ]; then
    say "calibration inputs missing: run with --with-calibration (downloads $((fastq_bytes / 1000000000)) GB of HG003 FASTQ, deleted once downsampled, and GIAB's HG003 benchmark)"
    return
  fi
  if [ "$truth_ok" = 0 ]; then
    fetch_truth hg003_truth "$hg003_truth_base" "$truth_dir" "$hg003_truth_files"
    say "GIAB HG003 v4.2.1 GRCh38 benchmark -> $truth_dir"
  fi
  if [ -z "$ds_dir" ]; then
    mkdir -p "$fastq_dir"
    for entry in $hg003_fastq_files; do
      IFS=: read -r name md5 size <<< "$entry"
      [ -f "$fastq_dir/$name" ] && [ "$(wc -c < "$fastq_dir/$name" | tr -d ' ')" = "$size" ] &&
        continue
      fetch_md5 "$hg003_fastq_base/$name" "$fastq_dir/$name" "$md5"
    done
    say "HG003 FASTQ -> $fastq_dir"
  fi
}

# Parascopy's data (see the header): fetched through fetch_md5, unpacked into a scratch folder
# and moved beside the tarballs only once it unpacked whole. A marker holding the tarball's
# md5 records the unpacking, so a tarball replaced later is unpacked again.
step_parascopy() {
  local dir=$lab/parascopy entry name md5 tmp member
  local wanted="" outputs=()
  for entry in $parascopy_files; do
    outputs+=("$dir/${entry%%:*}" "$dir/.unpacked.${entry%%:*}")
  done
  step_done parascopy "${outputs[@]}" && return
  for entry in $parascopy_files; do
    name=${entry%%:*}
    recorded parascopy "$dir/$name" && recorded parascopy "$dir/.unpacked.$name" && continue
    wanted+="Parascopy data $name ($parascopy_base/$name?download=1)"$'\n'
  done
  if [ "$check" = 1 ]; then
    while IFS= read -r entry; do [ -z "$entry" ] || missing "$entry"; done <<< "$wanted"
    return
  fi
  mkdir -p "$dir"
  for entry in $parascopy_files; do
    name=${entry%%:*}
    md5=${entry#*:}
    [ -f "$dir/$name" ] && [ "$(md5of "$dir/$name")" = "$md5" ] ||
      fetch_md5 "$parascopy_base/$name?download=1" "$dir/$name" "$md5"
    tmp=$dir/.unpack
    rm -rf "$tmp"
    mkdir "$tmp"
    tar -xzf "$dir/$name" -C "$tmp" || { rm -rf "$tmp"; die "cannot unpack $dir/$name"; }
    while IFS= read -r member; do
      rm -rf "${dir:?}/$(basename "$member")"
      mv "$member" "$dir/"
    done < <(find "$tmp" -mindepth 1 -maxdepth 1)
    rm -rf "$tmp"
    printf '%s\n' "$md5" > "$dir/.unpacked.$name"
  done
  record parascopy "${outputs[@]}"
  say "Parascopy data -> $dir"
}

# `genokeeper.calibrate.hg002 fetch KEY` for the reference library NAME; status 1 when the
# fetch exits 3 (no `genokeeper run` has recorded the bases it downsamples to yet), which
# leaves the library for a rerun after that run and lets the bootstrap go on to the steps
# the run needs; any other failure stops the bootstrap.
fetch_library() {
  local code=0
  genokeeper_python -m genokeeper.calibrate.hg002 fetch "$2" >&2 || code=$?
  case $code in
    0) return 0 ;;
    3) say "$1 waits for your first genokeeper run, which records the bases it is downsampled to; rerun with this flag after it"; return 1 ;;
    *) die "fetching $1 failed; rerun to continue with the files not yet downsampled" ;;
  esac
}

# Stage 3's HG002 inputs (see the header); only with --with-hg002. The libraries come from
# `genokeeper.calibrate.hg002 fetch`, which skips a downsampled file already there; of its
# layout only HG002's rows are this step's (HG003 and HG004 have their own steps).
step_hg002() {
  local layout truth_dir entry name lib sample aligned first second
  local wanted="" truth_ok=1 libs=() outputs=()
  layout=$(genokeeper_python -m genokeeper.calibrate.hg002 layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.calibrate.hg002 (uv sync)"; return; }
    die "cannot read the HG002 paths from genokeeper.calibrate.hg002"
  }
  truth_dir=$lab/giab/HG002-v5.0q-GRCh38
  for entry in $hg002_truth_files; do outputs+=("$truth_dir/${entry%%:*}"); done
  if ! step_done hg002_truth "${outputs[@]}"; then
    truth_ok=0
    for entry in $hg002_truth_files; do
      name=${entry%%:*}
      recorded hg002_truth "$truth_dir/$name" ||
        wanted+="GIAB HG002 v5.0q truth file $name ($hg002_truth_base/$name)"$'\n'
    done
  fi
  while IFS=$'\t' read -r lib sample aligned first second; do
    [ "$sample" = HG002 ] || continue
    [ -f "$aligned" ] || { [ -f "$first" ] && [ -f "$second" ]; } && continue
    libs+=("$lib")
    wanted+="HG002 library $lib: downsampled pair $first and $second (downloaded from ENA)"$'\n'
  done <<< "$layout"
  [ -n "$wanted" ] || return 0
  if [ "$check" = 1 ]; then
    while IFS= read -r entry; do [ -z "$entry" ] || missing "$entry"; done <<< "$wanted"
    return
  fi
  if [ "$with_hg002" != 1 ]; then
    say "HG002 inputs missing: run with --with-hg002 (reads 300 GB from ENA, keeps about 170 GB of downsampled FASTQ, and fetches GIAB's HG002 v5.0q truth)"
    return
  fi
  if [ "$truth_ok" = 0 ]; then
    fetch_truth hg002_truth "$hg002_truth_base" "$truth_dir" "$hg002_truth_files"
    say "GIAB HG002 v5.0q GRCh38 truth -> $truth_dir"
  fi
  for lib in ${libs[@]+"${libs[@]}"}; do
    say "downloading and downsampling HG002 $lib (one raw FASTQ at a time; a rerun resumes)"
    fetch_library "HG002 $lib" "$lib" || return 0
  done
}

# Slice 3b's HG002 truth: the v5.0q structural-variant files beside 3a's small-variant ones,
# under their own step so that a lab holding 3a's truth is not told to fetch it again.
step_hg002_sv() {
  local truth_dir entry name size pin url outputs=()
  truth_dir=$lab/giab/HG002-v5.0q-GRCh38
  while read -r name size pin url; do outputs+=("$truth_dir/$name"); done <<< "$hg002_sv_downloads"
  step_done hg002_sv_truth "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    while read -r name size pin url; do
      recorded hg002_sv_truth "$truth_dir/$name" ||
        missing "GIAB HG002 v5.0q SV truth file $name ($url)"
    done <<< "$hg002_sv_downloads"
    return
  fi
  if [ "$with_hg002" != 1 ]; then
    say "HG002 SV truth missing: run with --with-hg002 (downloads $(downloads_size "$hg002_sv_downloads" 1000000 0) MB into $truth_dir)"
    return
  fi
  obtain_downloads "$hg002_sv_downloads" "$truth_dir" "$lab/prefetch/3b"
  record hg002_sv_truth "${outputs[@]}"
  say "GIAB HG002 v5.0q GRCh38 SV truth -> $truth_dir"
}

# ClinGen's dosage sensitivity curations, fetched whole and never resumed (ClinGen rebuilds
# the file in place daily), filed under the date it names when its curations changed. Only a file genokeeper's reader
# takes is recorded: its anchor comes from that reader run on the bytes on this disk, which
# a maintenance page or another file cannot pass because the reader refuses a file without
# its date line, its header with the gene and haploinsufficiency columns, and a gene row.
step_clingen_dosage() {
  local work date dest
  step_done clingen_dosage && return
  if [ "$check" = 1 ]; then
    missing "ClinGen dosage sensitivity curations ($clingen_dosage_url)"
    return
  fi
  work=$lab/clingen-dosage/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $clingen_dosage_url"
  download -o "$work/$clingen_dosage_file" "$clingen_dosage_url"
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.clingen import load_dosage
print(load_dosage(Path(sys.argv[1])).date)' "$work/$clingen_dosage_file" 2> /dev/null) || {
    rm -rf "$work"
    die "$clingen_dosage_url is not ClinGen's gene curation list"
  }
  if dest=$(unchanged_copy "$lab/clingen-dosage" "$clingen_dosage_file" "$work/$clingen_dosage_file"); then
    rm -rf "$work"
    record clingen_dosage "$dest/$clingen_dosage_file"
    say "ClinGen dosage sensitivity unchanged since $dest"
    refreshed clingen_dosage "$dest"
    return
  fi
  dest=$lab/clingen-dosage/$date
  mkdir -p "$dest"
  mv "$work/$clingen_dosage_file" "$dest/$clingen_dosage_file"
  rm -rf "$work"
  record clingen_dosage "$dest/$clingen_dosage_file"
  say "ClinGen dosage sensitivity $date -> $dest/$clingen_dosage_file"
  refreshed clingen_dosage "$dest"
}

# ClinGen's region curations and recurrent CNV regions, fetched whole and never resumed and
# filed together under the region list's date when either changed, as clingen_dosage files
# its list. Only a pair genokeeper's readers take is recorded: its anchor comes from those
# readers run on the bytes on this disk, which a maintenance page or another file cannot
# pass because they refuse a list without its date line, its header with the region name
# and location columns and a row, and a BED without a nine-column recurrent (255,111,0) row.
step_clingen_regions() {
  local work date dest
  step_done clingen_regions && return
  if [ "$check" = 1 ]; then
    missing "ClinGen region curations and recurrent CNV regions ($clingen_regions_url, $clingen_recurrent_url)"
    return
  fi
  work=$lab/clingen-regions/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $clingen_regions_url and $clingen_recurrent_url"
  download -o "$work/$clingen_regions_file" "$clingen_regions_url"
  download -o "$work/$clingen_recurrent_file" "$clingen_recurrent_url"
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.clingen import load_recurrent, load_regions, regions_date
load_regions(Path(sys.argv[1]))
load_recurrent(Path(sys.argv[2]))
print(regions_date(Path(sys.argv[1])))' "$work/$clingen_regions_file" "$work/$clingen_recurrent_file" 2> /dev/null) || {
    rm -rf "$work"
    die "$clingen_regions_url and $clingen_recurrent_url is not ClinGen's region curation list and recurrent CNV BED"
  }
  if dest=$(unchanged_copy "$lab/clingen-regions" "$clingen_regions_file" "$work/$clingen_regions_file") &&
    cmp -s "$dest/$clingen_recurrent_file" "$work/$clingen_recurrent_file"; then
    rm -rf "$work"
    record clingen_regions "$dest/$clingen_regions_file" "$dest/$clingen_recurrent_file"
    say "ClinGen regions unchanged since $dest"
    refreshed clingen_regions "$dest"
    return
  fi
  dest=$lab/clingen-regions/$date
  # The BED changed under an unchanged region list: today's folder, never the old one.
  [ ! -e "$dest" ] || dest=$lab/clingen-regions/$(date -u +%Y%m%d)
  mkdir -p "$dest"
  mv "$work/$clingen_regions_file" "$dest/$clingen_regions_file"
  mv "$work/$clingen_recurrent_file" "$dest/$clingen_recurrent_file"
  rm -rf "$work"
  record clingen_regions "$dest/$clingen_regions_file" "$dest/$clingen_recurrent_file"
  say "ClinGen regions $date -> $dest"
  refreshed clingen_regions "$dest"
}

# GenCC's submissions (slice 3f), fetched whole with the server's Last-Modified on the file
# (curl -R; the daily check compares it when GenCC gives no ETag) and never resumed, or
# moved from prefetch/3f/ when that copy has its pin; filed under its newest
# submitted_run_date only when it differs from the newest copy. Only a file genokeeper's
# reader takes is recorded: its anchor comes from that reader run on the bytes on this
# disk, which an error page or a cut transfer cannot pass because the reader refuses a file
# without its columns, a submission row or a known classification title.
step_gencc() {
  local work date dest prefetch=$lab/prefetch/3f/$gencc_file
  step_done gencc && return
  if [ "$check" = 1 ]; then
    missing "GenCC gene-disease submissions ($gencc_url)"
    return
  fi
  work=$lab/gencc/.download
  rm -rf "$work"
  mkdir -p "$work"
  if [ -z "$refresh" ] && prefetched "$prefetch" "$gencc_size" "$gencc_sha256"; then
    mv "$prefetch" "$work/$gencc_file"
    say "GenCC submissions moved from $prefetch"
  else
    say "downloading $gencc_url"
    download -R -o "$work/$gencc_file" "$gencc_url"
  fi
  date=$(genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence.gencc import gencc_date, load_gencc
path = Path(sys.argv[1])
load_gencc(path)
print(gencc_date(path))' "$work/$gencc_file" 2> /dev/null) || {
    rm -rf "$work"
    die "$gencc_url is not GenCC's submissions export"
  }
  if dest=$(unchanged_copy "$lab/gencc" "$gencc_file" "$work/$gencc_file"); then
    rm -rf "$work"
    record gencc "$dest/$gencc_file"
    say "GenCC submissions unchanged since $dest"
    refreshed gencc "$dest"
    return
  fi
  dest=$lab/gencc/$date
  mkdir -p "$dest"
  mv "$work/$gencc_file" "$dest/$gencc_file"
  rm -rf "$work"
  record gencc "$dest/$gencc_file"
  say "GenCC submissions $date -> $dest/$gencc_file"
  refreshed gencc "$dest"
}

# PanelApp's panel, filed only with the sha256 pinned above: its anchor is that pin, which
# another version of the panel or a cut download cannot pass, as the comparison runs on the
# bytes on this disk before the rename. Its URL serves the panel's current version, so a miss
# there is expected once the panel moves on: it warns and leaves the rest of the bootstrap
# running, and the run reports the research layer as not run.
step_panelapp() {
  local dir=$lab/panelapp/285/$panelapp_version part
  local dest=$dir/panel.json
  step_done panelapp "$dest" && return
  if [ "$check" = 1 ]; then
    missing "PanelApp's Intellectual disability panel $panelapp_version ($panelapp_url)"
    return
  fi
  mkdir -p "$dir"
  part=$dest.part
  if ! download -o "$part" "$panelapp_url"; then
    rm -f "$part"
    say "could not fetch $panelapp_url: the PanelApp research layer stays not run"
    return
  fi
  if ! genokeeper_python -c '
import json, sys
with open(sys.argv[1], encoding="utf-8") as handle:
    panel = json.load(handle)
for kind in ("genes", "strs", "regions"):
    for entity in panel.get(kind, []):
        entity["evidence"].sort()
with open(sys.argv[1], "w", encoding="utf-8") as handle:
    json.dump(panel, handle, sort_keys=True, ensure_ascii=False)
' "$part"; then
    rm -f "$part"
    say "could not read $panelapp_url as a panel: the PanelApp research layer stays not run"
    return
  fi
  if [ "$(sha "$part")" != "$panelapp_sha256" ]; then
    rm -f "$part"
    say "sha256 mismatch for $panelapp_url, not panel $panelapp_version: the PanelApp research layer stays not run"
    return
  fi
  mv "$part" "$dest"
  record panelapp "$dest"
  say "PanelApp panel $panelapp_version -> $dest"
}

# A publication source, fetched whole and never resumed (-C - would splice two builds: LitVar2
# is rebuilt in place), filed under the server's Last-Modified date in UTC, which curl -R
# puts on the file. Only a file genokeeper.evidence.publications reads is recorded: its
# anchor is that reader run on the bytes on this disk, which an error page or a cut
# transfer cannot pass.
# publication_step STEP FOLDER FILE URL CHECKER LABEL
publication_step() {
  local step=$1 folder=$2 file=$3 url=$4 checker=$5 label=$6 work date dest
  step_done "$step" && return
  if [ "$check" = 1 ]; then
    missing "$label ($url)"
    return
  fi
  work=$lab/$folder/.download
  rm -rf "$work"
  mkdir -p "$work"
  say "downloading $url"
  download -R -D "$work/headers" -o "$work/$file" "$url"
  grep -qi '^last-modified:' "$work/headers" || {
    rm -rf "$work"
    die "$url answered with no Last-Modified date"
  }
  genokeeper_python -c '
import sys
from pathlib import Path
from genokeeper.evidence import publications
getattr(publications, sys.argv[2])(Path(sys.argv[1]))' "$work/$file" "$checker" 2> /dev/null || {
    rm -rf "$work"
    die "$url is not the $label file genokeeper reads"
  }
  date=$(file_day "$work/$file")
  [[ "$date" =~ ^[0-9]{8}$ ]] || { rm -rf "$work"; die "no Last-Modified date for $url"; }
  dest=$lab/$folder/$date
  mkdir -p "$dest"
  mv "$work/$file" "$dest/$file"
  rm -rf "$work"
  record "$step" "$dest/$file"
  say "$label $date -> $dest/$file"
  refreshed "$step" "$dest"
}

step_litvar() {
  publication_step litvar litvar "$litvar_file" "$litvar_url" check_litvar "LitVar2 variants"
}

step_clinvar_citations() {
  publication_step clinvar_citations clinvar-citations "$citations_file" "$citations_url" \
    check_citations "ClinVar citations"
}

# prefetched PATH SIZE SHA256: PATH is a file of that size with that sha256.
prefetched() {
  [ -f "$1" ] && [ "$(wc -c < "$1" | tr -d ' ')" = "$2" ] && [ "$(sha "$1")" = "$3" ]
}

# obtain PREFETCH NAME SIZE SHA256 URL DEST: NAME at DEST with the pinned sha256, moved from
# PREFETCH when the copy there has it, else downloaded through DEST.part (a cut transfer
# resumed) and renamed only once its sha256 matches. DEST under its final name is therefore
# a verified file: its anchor comes from the sha256 written in this script, which a server
# serving other bytes, a prefetched copy cut short, or a transfer resumed onto a changed
# file cannot move, because the comparison runs on the bytes on this disk.
obtain() {
  local prefetch=$1 name=$2 size=$3 pin=$4 url=$5 dest=$6
  if prefetched "$prefetch" "$size" "$pin"; then
    mv "$prefetch" "$dest"
    say "$name moved from $prefetch"
    return
  fi
  say "downloading $url"
  download -C - -o "$dest.part" "$url"
  if [ "$(sha "$dest.part")" != "$pin" ]; then
    rm -f "$dest.part"
    die "sha256 mismatch for $url; removed, rerun to fetch it again"
  fi
  mv "$dest.part" "$dest"
}

# vep_obtain NAME SUBFOLDER SIZE SHA256 URL DEST: obtain from prefetch/1e/SUBFOLDER/.
vep_obtain() { obtain "$lab/prefetch/1e/$2/$1" "$1" "$3" "$4" "$5" "$6"; }

# downloads_ready DOWNLOADS DEST_DIR PREFETCH_DIR: every row (name, size, sha256, URL) is
# already in DEST_DIR with its pin, or in PREFETCH_DIR with it, so obtaining it downloads
# nothing.
downloads_ready() {
  local name size pin url
  while read -r name size pin url; do
    prefetched "$2/$name" "$size" "$pin" ||
      prefetched "$3/$name" "$size" "$pin" || return 1
  done <<< "$1"
}

# obtain_downloads DOWNLOADS DEST_DIR PREFETCH_DIR: each row not yet in DEST_DIR with its
# pin, obtained from PREFETCH_DIR or downloaded.
obtain_downloads() {
  local name size pin url
  mkdir -p "$2"
  while read -r name size pin url; do
    prefetched "$2/$name" "$size" "$pin" ||
      obtain "$3/$name" "$name" "$size" "$pin" "$url" "$2/$name"
  done <<< "$1"
}

# downloads_size DOWNLOADS DIVISOR DECIMALS: the rows' total bytes over DIVISOR.
downloads_size() {
  awk -v d="$2" -v p="$3" '{ s += $2 } END { printf "%.*f", p, s / d }' <<< "$1"
}

# Slice 1e's VEP cache and score files (see the header); only with --with-vep. The
# downloads stay in vep/116/.download/ until every output is built, so a rerun after a
# failure fetches nothing again.
step_vep() {
  local layout dir info revel am sv pext pext_max nmd work entry tab rel
  local name sub size pin url bytes=0 outputs=()
  layout=$(genokeeper_python -m genokeeper.annotate.vep layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.annotate.vep (uv sync)"; return; }
    die "cannot read slice 1e's paths from genokeeper.annotate.vep"
  }
  # VepFiles.required(), in its order: the cache's info.txt, REVEL, AlphaMissense and
  # SpliceVault each with its .tbi, pext with its .tbi, pext's gene maxima, the NMD track.
  while IFS= read -r entry; do outputs+=("$entry"); done <<< "$layout"
  info=${outputs[0]} revel=${outputs[1]} am=${outputs[3]} sv=${outputs[5]}
  pext=${outputs[7]} pext_max=${outputs[9]} nmd=${outputs[10]} dir=${revel%/*}
  step_done vep "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    for entry in "${outputs[@]}"; do
      recorded vep "$entry" || missing "VEP $vep_release file $entry (run with --with-vep)"
    done
    return
  fi
  if [ "$with_vep" != 1 ]; then
    while read -r name sub size pin url; do bytes=$((bytes + size)); done <<< "$vep_downloads"
    say "VEP $vep_release cache and score files missing: run with --with-vep (downloads $((bytes / 1000000000)) GB into $dir)"
    return
  fi
  work=$dir/.download
  mkdir -p "$work"
  while read -r name sub size pin url; do
    [ -f "$work/$name" ] && [ "$(sha "$work/$name")" = "$pin" ] && continue
    vep_obtain "$name" "$sub" "$size" "$pin" "$url" "$work/$name"
  done <<< "$vep_downloads"
  rm -rf "$dir/.cache"
  mkdir -p "$dir/.cache"
  tar -xzf "$work/$vep_cache_tarball" -C "$dir/.cache"
  rel=${info#"$dir/cache/"}
  [ -f "$dir/.cache/$rel" ] || die "$vep_cache_tarball holds no $rel"
  # REVEL as its VEP plugin reads it (REVEL.pm's GRCh38 recipe): tab-separated, the header
  # marked with #, the rows without a GRCh38 position left out and sorted on it.
  tab=$(printf '\t')
  unzip -p "$work/$revel_zip" | tr ',' '\t' |
    awk -v h="$work/revel.header" 'NR == 1 { print "#" $0 > h; next } $3 != "."' |
    LC_ALL=C sort -T "$work" -t "$tab" -k1,1 -k3,3n > "$work/revel.sorted"
  cat "$work/revel.header" "$work/revel.sorted" | bgzip -c > "$work/revel.tsv.gz"
  tabix -f -s 1 -b 3 -e 3 "$work/revel.tsv.gz"
  tabix -f -s 1 -b 2 -e 2 -S 1 "$work/$alphamissense_file"
  # pext as BED (chrom, start, end, gene_id, exp_prop_mean; NaN rows left out) and each
  # gene's largest exp_prop_mean.
  gzip -dc "$work/$pext_tsv" |
    awk -F'\t' -v m="$work/pext.max" 'NR > 1 && $4 != "NaN" {
        split($1, l, ":"); print l[1] "\t" l[2] - 1 "\t" l[2] "\t" $2 "\t" $4
        if (!($2 in top) || $4 + 0 > top[$2]) top[$2] = $4 + 0
      } END { for (g in top) print g "\t" top[g] > m }' |
    bgzip -c > "$work/pext.bed.gz"
  tabix -f -p bed "$work/pext.bed.gz"
  LC_ALL=C sort -o "$work/pext.max" "$work/pext.max"
  genokeeper_python -m genokeeper.annotate.vep nmd-bed "$work/$nmd_bb" "$work/nmd.bed" >&2 ||
    die "converting $nmd_bb to BED failed"
  rm -rf "$dir/cache"
  mv "$dir/.cache" "$dir/cache"
  mv "$work/revel.tsv.gz" "$revel"
  mv "$work/revel.tsv.gz.tbi" "$revel.tbi"
  mv "$work/$alphamissense_file" "$am"
  mv "$work/$alphamissense_file.tbi" "$am.tbi"
  mv "$work/$splicevault_file" "$sv"
  mv "$work/$splicevault_file.tbi" "$sv.tbi"
  mv "$work/pext.bed.gz" "$pext"
  mv "$work/pext.bed.gz.tbi" "$pext.tbi"
  mv "$work/pext.max" "$pext_max"
  mv "$work/nmd.bed" "$nmd"
  rm -rf "$work"
  record vep "${outputs[@]}"
  say "VEP $vep_release cache and slice 1e's score files -> $dir"
}

# Slice 3b's gnomAD-SV (see the header). Files in prefetch/3b/ with their pins need no flag,
# moving them downloads nothing; anything else waits for --with-sv.
step_sv() {
  local layout entry dir
  local outputs=()
  layout=$(genokeeper_python -m genokeeper.deletions.files layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.deletions.files (uv sync)"; return; }
    die "cannot read slice 3b's paths from genokeeper.deletions.files"
  }
  while IFS= read -r entry; do outputs+=("$entry"); done <<< "$layout"
  dir=${outputs[0]%/*}
  step_done sv "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    for entry in "${outputs[@]}"; do
      recorded sv "$entry" || missing "gnomAD-SV file $entry (run with --with-sv)"
    done
    return
  fi
  if [ "$with_sv" != 1 ] && ! downloads_ready "$sv_downloads" "$dir" "$lab/prefetch/3b"; then
    say "gnomAD-SV missing: run with --with-sv (downloads $(downloads_size "$sv_downloads" 1000000000 2) GB into $dir)"
    return
  fi
  obtain_downloads "$sv_downloads" "$dir" "$lab/prefetch/3b"
  for entry in "${outputs[@]}"; do
    [ -f "$entry" ] || die "the sv downloads hold no $entry"
  done
  record sv "${outputs[@]}"
  say "gnomAD-SV -> $dir"
}

# Slice 3c's masked GRCh38 (see the header). In prefetch/3c/ with its pin it needs no flag.
# The .gz stays in reference/.masked/ until every output is built, so a rerun after a
# failure fetches nothing again; nothing is recorded before the SDF exists.
step_masked_reference() {
  local paths fasta unmasked bed sdf work gz fasta_sha ext outputs=()
  paths=$(genokeeper_python -c '
from genokeeper.calibrate import layout
from genokeeper.config import load_paths, unmasked_fasta
from genokeeper.masking import masked_bed
lab = load_paths().lab
print(layout(lab)["fasta"], unmasked_fasta(lab), masked_bed(lab), layout(lab)["sdf"], sep="\n")') || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.masking (uv sync)"; return; }
    die "cannot read the masked reference's paths from genokeeper"
  }
  { read -r fasta; read -r unmasked; read -r bed; read -r sdf; } <<< "$paths"
  outputs=("$fasta" "$fasta.fai" "$bed")
  for ext in amb ann bwt pac sa; do outputs+=("$fasta.$ext"); done
  outputs+=("$sdf/done" "$sdf/fasta.sha256")
  step_done masked_reference "${outputs[@]}" && return
  if [ "$check" = 1 ]; then
    missing "masked GRCh38 $fasta with its .fai, mask, bwa index and rtg SDF (run with --with-masked-reference; $masked_url)"
    return
  fi
  work=$(dirname "$fasta")/.masked
  gz=$work/$masked_name
  if [ "$with_masked_reference" != 1 ] && ! prefetched "$gz" "$masked_size" "$masked_sha256" &&
    ! prefetched "$lab/prefetch/3c/$masked_name" "$masked_size" "$masked_sha256"; then
    say "masked GRCh38 missing: run with --with-masked-reference (downloads $(downloads_size "$masked_name $masked_size" 1000000 0) MB into $(dirname "$fasta"))"
    return
  fi
  [ -n "$(fasta_recorded "$unmasked")" ] && [ -f "$unmasked.fai" ] ||
    die "no recorded GRCh38 FASTA at $unmasked: the grch38_fasta step comes first"
  rtg_ok || die "rtg-tools $rtg_version is missing from $(rtg_dir): the tools38 step comes first"
  mkdir -p "$work"
  prefetched "$gz" "$masked_size" "$masked_sha256" ||
    obtain "$lab/prefetch/3c/$masked_name" "$masked_name" "$masked_size" "$masked_sha256" \
      "$masked_url" "$gz"
  # Verified in reference/.masked/, renamed into place only once verified: a refused FASTA
  # never sits at the path every consumer reads.
  gzip -dc "$gz" > "$work/masked.fa"
  samtools faidx "$work/masked.fa"
  # Its anchor comes from the mask committed in git (genokeeper.masking.MASK_RESOURCE),
  # which a FASTA served with other bases cannot move, because `verify` compares the bases
  # on this disk with the unmasked FASTA's and refuses any difference but a base set to N,
  # then refuses unless the masked runs equal the committed ones.
  genokeeper_python -m genokeeper.masking verify "$unmasked" "$work/masked.fa" "$bed" >&2 ||
    die "the mask of $fasta differs from the committed one; nothing recorded"
  mv "$work/masked.fa.fai" "$fasta.fai"
  mv "$work/masked.fa" "$fasta"
  genokeeper_python -m genokeeper.masking index "$fasta" >&2 || die "bwa index of $fasta failed"
  fasta_sha=$(sha "$fasta")
  rm -rf "$sdf.tmp"
  "$(rtg_dir)/rtg" format -o "$sdf.tmp" "$fasta" >&2
  printf '%s\n' "$fasta_sha" > "$sdf.tmp/fasta.sha256"
  rm -rf "$sdf"
  mv "$sdf.tmp" "$sdf"
  rm -rf "$work"
  record masked_reference "${outputs[@]}"
  say "masked GRCh38 with its mask, bwa index and rtg SDF -> $fasta"
}

# Slice 3c's HG004 inputs (see the header); only with --with-hg004. Its library comes from
# `genokeeper.calibrate.hg002 fetch hg004`, which never fetches a recorded library again.
step_hg004() {
  local layout truth_dir row lib sample aligned first second name size pin url
  local wanted="" truth_ok=1 pair_ok=1 outputs=()
  layout=$(genokeeper_python -m genokeeper.calibrate layout) &&
    row=$(genokeeper_python -m genokeeper.calibrate.hg002 layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.calibrate (uv sync)"; return; }
    die "cannot read HG004's paths from genokeeper.calibrate"
  }
  truth_dir=$(awk -F'\t' '$1 == "hg004_truth" { print $2 }' <<< "$layout")
  IFS=$'\t' read -r lib sample aligned first second < <(awk -F'\t' '$1 == "hg004"' <<< "$row")
  while read -r name size pin url; do outputs+=("$truth_dir/$name"); done <<< "$hg004_truth_downloads"
  if ! step_done hg004_truth "${outputs[@]}"; then
    truth_ok=0
    while read -r name size pin url; do
      recorded hg004_truth "$truth_dir/$name" ||
        wanted+="GIAB HG004 v4.2.1 truth file $name ($url)"$'\n'
    done <<< "$hg004_truth_downloads"
  fi
  if ! [ -f "$aligned" ] && ! { [ -f "$first" ] && [ -f "$second" ]; }; then
    pair_ok=0
    wanted+="HG004 library $lib: downsampled pair $first and $second (downloaded from ENA)"$'\n'
  fi
  [ -n "$wanted" ] || return 0
  if [ "$check" = 1 ]; then
    while IFS= read -r row; do [ -z "$row" ] || missing "$row"; done <<< "$wanted"
    return
  fi
  if [ "$with_hg004" != 1 ]; then
    say "HG004 inputs missing: run with --with-hg004 (downloads $(( $(downloads_size "$hg004_truth_downloads" 1 0) / 1000000 )) MB of GIAB's HG004 truth into $truth_dir, downloads 144 GB of HG004's FASTQ from CNGB and keeps its downsampled pair)"
    return
  fi
  if [ "$truth_ok" = 0 ]; then
    obtain_downloads "$hg004_truth_downloads" "$truth_dir" "$lab/prefetch/3c"
    record hg004_truth "${outputs[@]}"
    say "GIAB HG004 v4.2.1 GRCh38 truth -> $truth_dir"
  fi
  if [ "$pair_ok" = 0 ]; then
    say "downloading and downsampling HG004 (one raw FASTQ at a time; a rerun resumes)"
    fetch_library HG004 hg004 || return 0
  fi
}

# Slice 3f's typing inputs (see the header). A file in prefetch/3f/ with its pin needs no
# flag, moving it downloads nothing; anything else waits for --with-typing.

# typing_prefetch PATH: where prefetch/3f/ keeps the file the row PATH (under the lab) names.
typing_prefetch() {
  case $1 in
    hmmix/*/archaic/*) echo "$lab/prefetch/3f/archaic/${1##*/}" ;;
    hmmix/*) echo "$lab/prefetch/3f/hmmix/${1##*/}" ;;
    *) echo "$lab/prefetch/3f/${1##*/}" ;;
  esac
}

# typing_ready: every row of $typing_downloads is already in place with its pin, or in
# prefetch/3f/ with it, so obtaining it downloads nothing.
typing_ready() {
  local path size pin url
  while read -r path size pin url; do
    prefetched "$lab/$path" "$size" "$pin" ||
      prefetched "$(typing_prefetch "$path")" "$size" "$pin" || return 1
  done <<< "$typing_downloads"
}

# typing_obtain: each row not yet in place with its pin, obtained from prefetch/3f/ or
# downloaded.
typing_obtain() {
  local path size pin url
  while read -r path size pin url; do
    mkdir -p "$lab/${path%/*}"
    prefetched "$lab/$path" "$size" "$pin" ||
      obtain "$(typing_prefetch "$path")" "${path##*/}" "$size" "$pin" "$url" "$lab/$path"
  done <<< "$typing_downloads"
}

# The tools of step_typing, apart so a test can stand in for each.
typing_bwa_index() { bwa index "$1" >&2; }
typing_build_panel() { genokeeper_python -m genokeeper.typing_files build-panel "$lab" >&2; }
typing_build_af() { genokeeper_python -m genokeeper.typing_files build-af "$lab" >&2; }

# Slice 3f's typing inputs: HHV-6A and B concatenated into hhv6/hhv6.fa and bwa-indexed,
# hmmix's files, then the supervised-admixture reference and the per-superpopulation AF files
# built from the HGDP+1kGP panel when their outputs are missing. After pgs (the panel) and
# tools38 (bwa).
step_typing() {
  local layout entry fasta part parts=() build_panel=0 build_af=0
  local outputs=()
  layout=$(genokeeper_python -m genokeeper.typing_files layout) || {
    [ "$check" = 1 ] && { missing "project environment for genokeeper.typing_files (uv sync)"; return; }
    die "cannot read slice 3f's paths from genokeeper.typing_files"
  }
  while IFS= read -r entry; do outputs+=("$entry"); done <<< "$layout"
  fasta=${outputs[0]}
  for part in $hhv6_parts; do parts+=("${fasta%/*}/$part.fa"); done
  step_done typing "${outputs[@]}" "${parts[@]}" && return
  if [ "$check" = 1 ]; then
    for entry in "${outputs[@]}"; do
      recorded typing "$entry" || missing "typing file $entry (run with --with-typing)"
    done
    return
  fi
  if [ "$with_typing" != 1 ] && ! typing_ready; then
    say "typing inputs missing: run with --with-typing (downloads $(downloads_size "$typing_downloads" 1000000000 1) GB into $lab)"
    return
  fi
  typing_obtain
  # The two genomes in one FASTA, in the order of $hhv6_parts.
  [ -f "$fasta" ] || {
    cat "${parts[@]}" > "$fasta.part"
    mv "$fasta.part" "$fasta"
  }
  for entry in "${outputs[@]:1:5}"; do
    [ -f "$entry" ] || { typing_bwa_index "$fasta" || die "bwa index of $fasta failed"; break; }
  done
  for entry in "${outputs[@]}"; do
    [ -f "$entry" ] || case $entry in
      "$lab"/ancestry/*) build_panel=1 ;;
      "$lab"/roh/*) build_af=1 ;;
    esac
  done
  if [ "$build_panel" = 1 ] || [ "$build_af" = 1 ]; then
    [ -f "$lab/pgs/$pgs_panel_name" ] || die "the HGDP+1kGP panel $lab/pgs/$pgs_panel_name is missing: run the pgs step"
  fi
  if [ "$build_panel" = 1 ]; then
    typing_build_panel || die "building the supervised-admixture reference failed"
  fi
  if [ "$build_af" = 1 ]; then
    typing_build_af || die "building the ROH allele frequencies failed"
  fi
  for entry in "${outputs[@]}"; do
    [ -f "$entry" ] || die "the typing step made no $entry"
  done
  record typing "${outputs[@]}" "${parts[@]}"
  say "typing inputs -> $lab"
}

usage() {
  echo "usage: $0 [--check] [--with-gnomad] [--with-calibration] [--with-hg002] [--with-vep] [--with-sv] [--with-masked-reference] [--with-hg004] [--with-typing] [--with-publications]" >&2
  echo "       $0 --refresh STEP   (STEP: $refresh_steps)" >&2
}

parse_args() {
  local arg
  while [ "$#" -gt 0 ]; do
    arg=$1
    shift
    case "$arg" in
      --check) check=1 ;;
      --with-gnomad) with_gnomad=1 ;;
      --with-calibration) with_calibration=1 ;;
      --with-hg002) with_hg002=1 ;;
      --with-vep) with_vep=1 ;;
      --with-sv) with_sv=1 ;;
      --with-masked-reference) with_masked_reference=1 ;;
      --with-hg004) with_hg004=1 ;;
      --with-typing) with_typing=1 ;;
      --with-publications) with_publications=1 ;;
      --refresh)
        [ "$#" -gt 0 ] || { usage; exit 64; }
        case "$1" in "" | *[!a-z_]*) usage; exit 64 ;; esac # one step name, one word
        case " $refresh_steps " in *" $1 "*) refresh=$1 ;; *) usage; exit 64 ;; esac
        shift
        ;;
      *) usage; exit 64 ;;
    esac
  done
  # A refresh fetches; --check only lists.
  if [ -n "$refresh" ] && [ "$check" = 1 ]; then usage; exit 64; fi
}

main() {
  parse_args "$@"
  [ "$check" = 1 ] || mkdir -p "$lab"
  if [ -n "$refresh" ]; then
    "step_$refresh"
    return
  fi
  step_tools
  step_tools38
  step_clinvar
  step_clingen
  step_giab38
  step_grch38_fasta
  step_dv_model
  step_cgd
  step_pharmcat
  step_masked_reference
  step_cpic
  step_pharmgkb
  step_gwas
  step_medlineplus
  step_clingen_dosage
  step_clingen_regions
  step_gencc
  step_panelapp
  if [ "$with_publications" = 1 ]; then
    step_litvar
    step_clinvar_citations
  fi
  step_images
  step_t1k
  step_pgs
  step_typing
  step_vep
  step_sv
  step_hg003
  step_parascopy
  step_hg002
  step_hg002_sv
  step_hg004
  step_gnomad_faf
  [ "$check" = 1 ] || say "done; manifest $manifest"
}

# `return` works only in a sourced file, so this never runs a step when sourced, whatever
# $0 the sourcing shell has.
if ! (return 0 2> /dev/null); then
  main "$@"
fi
