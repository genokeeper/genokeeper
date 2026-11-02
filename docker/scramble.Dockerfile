# Scramble 1.0.2 (github.com/GeneDx/scramble tag 1.0.2, commit
# 14806b4569d20e85835ace7f0e32419085d9bae3, CC BY-NC-SA 4.0, personal non-commercial use,
# read 2026-10-30; run, never redistributed), for the mobile-element insertion calls. Bioconda
# has no linux-aarch64 build of it (its recipe is skiplisted), so this image is built here, for
# linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - htslib, R with optparse, stringr, Rsamtools, Biostrings and rBLAST, BLAST+, make, a C compiler
#   and liblzma's headers by scramble.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-30). Biostrings is held below 2.77.1, where
#   pairwiseAlignment() became defunct in favour of the pwalign package (SCRAMble.R's do.meis.R
#   and do.dels.R call it unqualified): bioconductor-biostrings 2.74.0 on r-base 4.4.3, with
#   bioconductor-rsamtools 2.22.0 and bioconductor-rblast 1.2.0 (Bioconductor 3.20), the oldest
#   set bioconda builds for linux-aarch64. Biostrings 2.74.0 still forwards the call, but only
#   to an installed pwalign, so bioconductor-pwalign 1.2.0 is in the lock too; the check below
#   runs the alignment the way do.meis.R does;
# - Scramble itself by the git commit tag 1.0.2 points to: cluster_identifier is compiled from
#   its C sources against that htslib, and its cluster_analysis/ folder (SCRAMble.R and the
#   MEI consensus FASTA) is kept at /opt/scramble/cluster_analysis/. Upstream's
#   `cluster_identifier -v` prints 1.0.1 at this tag, and a call without arguments exits 1, so
#   the check below reads its usage line.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER scramble.lock /tmp/
RUN micromamba create --yes --name scramble --file /tmp/scramble.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=scramble PATH=/opt/conda/envs/scramble/bin:$PATH
ADD --keep-git-dir=false --chown=$MAMBA_USER:$MAMBA_USER https://github.com/GeneDx/scramble.git#14806b4569d20e85835ace7f0e32419085d9bae3 /tmp/scramble-src
RUN make -C /tmp/scramble-src/cluster_identifier/src \
        CC=aarch64-conda-linux-gnu-gcc \
        INCLUDES=-I/opt/conda/envs/scramble/include \
        CFLAGS="-Wall -O2 -L/opt/conda/envs/scramble/lib -Wl,-rpath,/opt/conda/envs/scramble/lib" \
    && install -m 0755 /tmp/scramble-src/cluster_identifier/src/build/cluster_identifier /opt/conda/envs/scramble/bin/
USER root
RUN mkdir -m 0755 /opt/scramble \
    && cp -r /tmp/scramble-src/cluster_analysis /opt/scramble/cluster_analysis \
    && chmod -R a+rX /opt/scramble \
    && rm -rf /tmp/scramble-src
USER $MAMBA_USER
# Offline, as every run.
RUN --network=none (cluster_identifier 2>&1 || true) | grep 'Usage: cluster_identifier' \
    && test -s /opt/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa \
    && test -s /opt/scramble/cluster_analysis/bin/SCRAMble.R \
    && Rscript -e 'library(optparse); library(rBLAST); library(Rsamtools); library(stringr); library(Biostrings); ref <- readDNAStringSet("/opt/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa"); a <- pairwiseAlignment(pattern = tolower(as.character(subseq(ref[[1]], 100, 160))), subject = ref[[1]], type = "local"); stopifnot(nchar(a) == 61, score(a) > 0, pid(a) > 99)'
