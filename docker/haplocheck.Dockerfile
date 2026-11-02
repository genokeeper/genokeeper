# haplocheck 1.3.3 (github.com/genepi/haplocheck, MIT at tag v1.3.3, read 2026-10-30; still
# its newest release, published 2021-02-06), for contamination in the mitochondrial calls.
# Upstream publishes no image and BioContainers' is linux/amd64 only, so it is built here, for
# linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by haplocheck.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-30): bioconda's noarch haplocheck 1.3.3
#   h2a3209d_2 (its jar) on conda-forge's linux-aarch64 openjdk 25.0.2.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER haplocheck.lock /tmp/
RUN micromamba create --yes --name haplocheck --file /tmp/haplocheck.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=haplocheck PATH=/opt/conda/envs/haplocheck/bin:$PATH
# Offline, as every run.
RUN --network=none haplocheck --version 2>&1 | grep '1\.3\.3'
