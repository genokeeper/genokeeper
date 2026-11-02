# INSurVeyor 1.1.3 (github.com/kensung-lab/INSurVeyor tag 1.1.3, commit
# 914e65ae0a305ca6768076a0235f21239cebddd7, GPL-3.0, read 2026-10-30; run, never redistributed),
# for the insertion calls. BioContainers' image is not used: this one is built
# here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by insurveyor.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-30): bioconda's insurveyor 1.1.3 h1e84f2d_2 and
#   samtools 1.21 on htslib 1.21, with the python 3.13 the package resolves to.
# Upstream's insurveyor.py still has VERSION = "1.1.2" at that tag, so `--version` prints
# v1.1.2; the check below reads the package's version from conda's own record.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER insurveyor.lock /tmp/
RUN micromamba create --yes --name insurveyor --file /tmp/insurveyor.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=insurveyor PATH=/opt/conda/envs/insurveyor/bin:$PATH
# Offline, as every run.
RUN --network=none test -f /opt/conda/envs/insurveyor/conda-meta/insurveyor-1.1.3-h1e84f2d_2.json \
    && insurveyor.py --version 2>&1 | grep 'INSurVeyor v1\.1\.' && samtools --version > /dev/null
