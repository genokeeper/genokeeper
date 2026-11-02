# REViewer 0.2.7 (github.com/Illumina/REViewer, GPL-3 at tag v0.2.7, read 2026-10-24; its master
# branch has been PolyForm Strict since 2026-04, which is not used here), for the read-pile
# images of repeat expansions. BioContainers' image is amd64 only, so it is built here, for
# linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by reviewer.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-26): bioconda's reviewer 0.2.7 h739925f_0, which
#   pins htslib 1.19.1 (ExpansionHunter's pin differs, so it is another image).
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER reviewer.lock /tmp/
RUN micromamba create --yes --name reviewer --file /tmp/reviewer.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=reviewer PATH=/opt/conda/envs/reviewer/bin:$PATH
# Offline, as every run.
RUN --network=none REViewer --help > /dev/null
