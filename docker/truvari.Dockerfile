# Truvari 5.5.0 (github.com/ACEnglish/truvari, MIT at tag v5.5.0, read 2026-10-27), for the
# benchmark of structural-variant calls against the Genome in a Bottle truth set. BioContainers'
# image is not used: this one is built here, for linux/arm64
# (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by truvari.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-27): bioconda's noarch truvari 5.5.0 pyhdfd78af_1
#   on python 3.12.16, with bcftools 1.25 and htslib 1.24 (its tabix) for the VCFs Truvari
#   reads and writes.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER truvari.lock /tmp/
RUN micromamba create --yes --name truvari --file /tmp/truvari.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=truvari PATH=/opt/conda/envs/truvari/bin:$PATH
# Offline, as every run.
RUN --network=none truvari version 2>&1 | grep '5\.5\.0' && tabix --version | head -1 && bcftools --version | head -1
