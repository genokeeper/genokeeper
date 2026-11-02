# Gauchian 1.0.2 (github.com/Illumina/Gauchian tag v1.0.2, GPLv3 in LICENSE.txt and setup.py,
# read 2026-10-30; run, never redistributed), for the GBA gene's copy number and GBA/GBAP1
# conversion calls. Gauchian publishes no image, so it is built here, for linux/arm64
# (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - python 3.12 with pysam, numpy, scipy and statsmodels by gauchian.lock, an explicit lockfile
#   (`micromamba create --dry-run` for linux-aarch64 inside that image, 2026-10-30);
# - Gauchian itself by the git commit v1.0.2 points to (371a6d79599389e590038c4fecf13ec853abb8ab),
#   whose package_data ships its gauchian/data/ files (the GRCh38 *_38* ones among them).
#   Its setup.py lists pytest-runner under setup_requires, which would fetch an unpinned egg
#   from PyPI at build time; that line is deleted, as the package needs it for tests only.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER gauchian.lock /tmp/
RUN micromamba create --yes --name gauchian --file /tmp/gauchian.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=gauchian PATH=/opt/conda/envs/gauchian/bin:$PATH
ADD --keep-git-dir=false --chown=$MAMBA_USER:$MAMBA_USER https://github.com/Illumina/Gauchian.git#371a6d79599389e590038c4fecf13ec853abb8ab /tmp/gauchian-src
RUN sed -i '/setup_requires/d' /tmp/gauchian-src/setup.py \
    && pip install --no-cache-dir --no-deps --no-build-isolation /tmp/gauchian-src \
    && rm -rf /tmp/gauchian-src
# Offline, as every run.
RUN --network=none gauchian --help > /dev/null \
    && test -s /opt/conda/envs/gauchian/lib/python3.12/site-packages/gauchian/data/GBA_region_38.bed \
    && test -s /opt/conda/envs/gauchian/lib/python3.12/site-packages/gauchian/data/GBA_SNP_38.txt
