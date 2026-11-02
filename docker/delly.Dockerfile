# DELLY 2.7.0 (github.com/dellytools/delly, BSD-3-Clause at tag v2.7.0, read 2026-10-27), for
# the structural-variant calls (deletions). BioContainers' image is not used: this one is built
# here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by delly.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-27): bioconda's delly 2.7.0 ha658904_0 and
#   bcftools 1.25 (DELLY writes BCF), on htslib 1.24;
# - the GRCh38 exclude list, human.hg38.excl.tsv, from the repository's excludeTemplates/ at
#   the commit tag v2.7.0 points to (3d505adc17afccb4443cf97dd01f78afc1d5a267, read
#   2026-10-27), by the sha256 of the whole file; ADD refuses other bytes.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER delly.lock /tmp/
RUN micromamba create --yes --name delly --file /tmp/delly.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=delly PATH=/opt/conda/envs/delly/bin:$PATH
# ADD makes a missing parent directory with its --chmod (0644: closed to the unprivileged user),
# so the directory is made first.
USER root
RUN mkdir -m 0755 /opt/delly
USER $MAMBA_USER
ADD --checksum=sha256:caef8593f82f513694ae41b429680eab8abe39c057a119dae7760715aef006ee --chmod=0644 \
    https://raw.githubusercontent.com/dellytools/delly/3d505adc17afccb4443cf97dd01f78afc1d5a267/excludeTemplates/human.hg38.excl.tsv \
    /opt/delly/human.hg38.excl.tsv
# Offline, as every run.
RUN --network=none delly --version 2>&1 | grep '2\.7\.0' && test -s /opt/delly/human.hg38.excl.tsv
