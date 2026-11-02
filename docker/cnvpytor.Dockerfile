# CNVpytor 1.3.2 (github.com/abyzovlab/CNVpytor tag v1.3.2, MIT, read 2026-10-30), for the
# read-depth copy-number calls. BioContainers' image is not used: this one is built here, for
# linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - every package by cnvpytor.lock, an explicit lockfile (`micromamba create --dry-run` for
#   linux-aarch64 inside that image, 2026-10-30): bioconda's cnvpytor 1.3.2 pyhdfd78af_0 on
#   python 3.12;
# - the GRCh38 GC-content and mask files, gc_hg38.pytor and mask_hg38.pytor, from the
#   repository's cnvpytor/data/ at the commit tag v1.3.2 points to
#   (3adc56fc6433ab8d76db1f151720cc4a4e487f20, read 2026-10-30), each by the sha256 of the
#   whole file (its git blob sha1 matched the repository's, by git hash-object); ADD refuses
#   other bytes. CNVpytor 1.3.2 reads them from its own package's data/ directory and, when
#   one is missing, downloads it from GitHub, which a run offline cannot. These two are the
#   only data the image carries; nothing else is downloaded.
# - empty placeholders for the other built-in genomes' resource files (hg19, chm13v2.0,
#   chm13v1.1, kn99; names read from cnvpytor/genome.py inside the image). CNVpytor's
#   __main__.py exits 0, doing nothing, unless every built-in genome's gc and mask file exists
#   (Genome.check_resources, existence only), `--version` included. The placeholders only get
#   past that check: a run that resolves to another genome reads an empty file, and the
#   pipeline's runs name GRCh38.
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --chown=$MAMBA_USER:$MAMBA_USER cnvpytor.lock /tmp/
RUN micromamba create --yes --name cnvpytor --file /tmp/cnvpytor.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=cnvpytor PATH=/opt/conda/envs/cnvpytor/bin:$PATH
ADD --checksum=sha256:96dc8423c95792380b8750e6c500fb6dd432b7fecf6c0a219998e3d3aa19395b --chmod=0644 \
    https://raw.githubusercontent.com/abyzovlab/CNVpytor/3adc56fc6433ab8d76db1f151720cc4a4e487f20/cnvpytor/data/gc_hg38.pytor \
    /opt/conda/envs/cnvpytor/lib/python3.12/site-packages/cnvpytor/data/gc_hg38.pytor
ADD --checksum=sha256:389ba0dc3ec0c6a9fb7c9275fcdcc85cd2a939f5313a4ac6d3eac01d51c7026c --chmod=0644 \
    https://raw.githubusercontent.com/abyzovlab/CNVpytor/3adc56fc6433ab8d76db1f151720cc4a4e487f20/cnvpytor/data/mask_hg38.pytor \
    /opt/conda/envs/cnvpytor/lib/python3.12/site-packages/cnvpytor/data/mask_hg38.pytor
RUN cd /opt/conda/envs/cnvpytor/lib/python3.12/site-packages/cnvpytor/data \
    && touch gc_hg19.pytor mask_hg19.pytor gc_chm13v2.0.pytor gc_chm13v1.1.pytor gc_kn99.pytor
# Offline, as every run.
RUN --network=none cnvpytor --version 2>&1 | grep '1\.3\.2' \
    && test -s /opt/conda/envs/cnvpytor/lib/python3.12/site-packages/cnvpytor/data/gc_hg38.pytor \
    && test -s /opt/conda/envs/cnvpytor/lib/python3.12/site-packages/cnvpytor/data/mask_hg38.pytor
