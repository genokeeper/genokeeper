# Yleaf 4.1.6 (github.com/genid/Yleaf tag v4.1.6, released 2026-10-19; LICENSE.txt is GPL-3,
# setup.py says MIT, both read at that tag 2026-10-26), for the Y haplogroup. Yleaf publishes no
# image, so it is built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned
# inputs, each checked before it is used:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - Yleaf's dependencies and samtools 1.24 by yleaf.lock, an explicit lockfile
#   (`micromamba create --dry-run` for linux-aarch64 inside that image, 2026-10-26);
# - Yleaf itself by the git commit v4.1.6 points to (git ls-remote, 2026-10-26); its YFull tree
#   and hg38 position files come with it. A run passes the chrY reference with --ref-fasta.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER yleaf.lock /tmp/
RUN micromamba create --yes --name yleaf --file /tmp/yleaf.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=yleaf PATH=/opt/conda/envs/yleaf/bin:$PATH
ENV YLEAF_NO_UPDATE_CHECK=1
ADD --keep-git-dir=false --chown=$MAMBA_USER:$MAMBA_USER https://github.com/genid/Yleaf.git#c167c4e8951623a53c370e7c300a0cfc7f411b62 /tmp/yleaf-src
RUN pip install --no-cache-dir --no-deps --no-build-isolation /tmp/yleaf-src \
    && rm -rf /tmp/yleaf-src
# Offline, as every run.
RUN --network=none Yleaf --help > /dev/null
