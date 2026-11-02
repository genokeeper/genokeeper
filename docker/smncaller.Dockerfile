# SMNCopyNumberCaller 1.1.2 (github.com/Illumina/SMNCopyNumberCaller, Apache-2.0 at tag v1.1.2,
# read 2026-10-24), for the SMN1/SMN2 copy numbers. BioContainers' image is amd64 only, so it is
# built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by smncaller.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-26): bioconda's smncopynumbercaller
#   1.1.2 py312h7e72e81_1 (noarch), with python pinned to 3.12 because its recipe does not pin one.
# The package ships the GRCh38 data files; the last step checks they are there.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER smncaller.lock /tmp/
RUN micromamba create --yes --name smncaller --file /tmp/smncaller.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=smncaller PATH=/opt/conda/envs/smncaller/bin:$PATH
# Offline, as every run: the program starts and its GRCh38 data files are in the package.
RUN --network=none smn_caller.py --help > /dev/null \
    && for f in SMN_region_38.bed SMN_SNP_38.txt SMN_target_variant_38.txt SMN_gmm.txt; do \
        find /opt/conda/envs/smncaller -name "$f" | grep -q . || exit 1; done
