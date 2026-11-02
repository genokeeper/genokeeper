# KILDA 1.2.3 (github.com/HCL-HUBL/KILDA tag v1.2.3, AGPL-3.0 at that tag, read 2026-10-30; run,
# never redistributed or served), for the LPA KIV-2 repeat's copy number from k-mer counts.
# KILDA publishes no image, so it is built here, for linux/arm64
# (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - jellyfish 2.3.1, samtools and python 3.12 with pandas, numpy and matplotlib by kilda.lock,
#   an explicit lockfile (`micromamba create --dry-run` for linux-aarch64 inside that image,
#   2026-10-30);
# - KILDA itself by the git commit v1.2.3 points to (b0ddf73f43636ddd7acef218bfda55b36e45fa1e):
#   its bin/kilda.py and data/ folder (the curated k-mer lists and BED files) are kept at
#   /opt/kilda/.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER kilda.lock /tmp/
RUN micromamba create --yes --name kilda --file /tmp/kilda.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=kilda PATH=/opt/conda/envs/kilda/bin:$PATH
ADD --keep-git-dir=false --chown=$MAMBA_USER:$MAMBA_USER https://github.com/HCL-HUBL/KILDA.git#b0ddf73f43636ddd7acef218bfda55b36e45fa1e /tmp/kilda-src
USER root
RUN mkdir -m 0755 /opt/kilda \
    && cp -r /tmp/kilda-src/bin /tmp/kilda-src/data /opt/kilda/ \
    && chmod -R a+rX /opt/kilda \
    && rm -rf /tmp/kilda-src
USER $MAMBA_USER
# Offline, as every run.
RUN --network=none python /opt/kilda/bin/kilda.py -V \
    && jellyfish --version | grep '2\.3\.1' \
    && test -s /opt/kilda/data/kmers_GRCh38_LPA_k31/KIV2_hg38_kmers_6copies_specific.tsv
