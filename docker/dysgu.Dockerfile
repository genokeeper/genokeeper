# Dysgu 1.9.0 (github.com/kcleal/dysgu, MIT at tag v1.9.0, read 2026-10-27), for the
# structural-variant calls (deletions). BioContainers' image is not used: this one is built
# here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by dysgu.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-27): bioconda's dysgu 1.9.0 py312hed307a4_0, on
#   python 3.12.16.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER dysgu.lock /tmp/
RUN micromamba create --yes --name dysgu --file /tmp/dysgu.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=dysgu PATH=/opt/conda/envs/dysgu/bin:$PATH
# Offline, as every run.
RUN --network=none dysgu --version 2>&1 | grep '1\.9\.0'
