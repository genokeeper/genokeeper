# Parascopy 1.19.0 (github.com/tprodanov/parascopy, MIT at tag v1.19.0, read 2026-10-24), for
# the SMN1/SMN2 paralog-specific copy numbers. BioContainers' image is amd64 only, so it is
# built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by parascopy.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-26): bioconda's parascopy 1.19.0 py312haa49b00_0,
#   which pins htslib 1.21 (SMNCopyNumberCaller's environment is another image).
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER parascopy.lock /tmp/
RUN micromamba create --yes --name parascopy --file /tmp/parascopy.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=parascopy PATH=/opt/conda/envs/parascopy/bin:$PATH
# Offline, as every run.
RUN --network=none parascopy --version | grep -q '^Parascopy v1\.19\.0$'
