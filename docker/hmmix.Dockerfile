# hmmix 0.9.2 (github.com/LauritsSkov/Introgression-detection, which has no tags or releases;
# read 2026-10-30: the repository's LICENSE is MIT and the wheel's classifier says MIT, but the
# 0.9.2 wheel's own LICENSE.txt is the GPL-3.0 text, so it is treated as GPL-3.0: run, never
# redistributed), for the Neanderthal share. Neither upstream nor BioContainers publishes an
# image, so it is built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned
# inputs, each checked before it is used:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - Python 3.10.22 (hmmix requires Python >3.5,<3.11), its pip, and bioconda's bcftools 1.25
#   on htslib 1.25 and vcftools 0.1.17, which hmmix runs through its shell (its create_ingroup
#   pipes into vcftools and writes an empty obs file when vcftools is missing), by hmmix.lock,
#   an explicit lockfile (`micromamba create --dry-run` for linux-aarch64 inside that image,
#   2026-10-30);
# - hmmix itself and its dependencies (numba's aarch64 wheel among them) by the sha256s in
#   hmmix.requirements.txt.
# genomelab.containers names the image by the sha256 of this file and the two lists.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER hmmix.lock hmmix.requirements.txt /tmp/
RUN micromamba create --yes --name hmmix --file /tmp/hmmix.lock \
    && micromamba clean --all --yes \
    && /opt/conda/envs/hmmix/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/hmmix.requirements.txt
ENV ENV_NAME=hmmix PATH=/opt/conda/envs/hmmix/bin:$PATH
# Offline, as every run: numba loads, the subcommands the slice runs are offered, and both
# tools hmmix shells out to run.
RUN --network=none hmmix -h > /tmp/help.txt \
    && grep -q create_ingroup /tmp/help.txt && grep -q decode /tmp/help.txt \
    && python -c "import numba, hmm_functions" \
    && bcftools --version | head -1 && vcftools --version \
    && rm /tmp/help.txt
