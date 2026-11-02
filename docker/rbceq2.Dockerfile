# RBCeq2 2.4.4 (github.com/limcintyre/RBCeq2 tag v2.4.4; MIT, the wheel's licenses/LICENSE and
# License-Expression, read 2026-10-30; "NOT FOR CLINICAL USE" per its README), for the blood
# group and platelet antigen alleles. Neither upstream nor BioContainers publishes an image, so
# it is built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs, each
# checked before it is used:
# - the base image by its index digest (python:3.12-slim-bookworm, Docker Hub, read
#   2026-10-20; RBCeq2 requires Python 3.12 or newer);
# - RBCeq2 itself, its ISBT allele table (rbceq2/resources/db.tsv, inside the wheel) and every
#   dependency, by the sha256s in rbceq2.requirements.txt.
# genomelab.containers names the image by the sha256 of this file and the list.
FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
COPY rbceq2.requirements.txt /tmp/
RUN python -m venv /opt/rbceq2 \
    && /opt/rbceq2/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/rbceq2.requirements.txt
ENV PATH=/opt/rbceq2/bin:$PATH
# The allele table a run reads, for the step that records its digest.
LABEL rbceq2.db=/opt/rbceq2/lib/python3.12/site-packages/rbceq2/resources/db.tsv
# Offline, as every run: it prints `rbceq2 2.4.4 (DB: <the table's version>)` on stdout.
RUN --network=none rbceq2 --version 2>/dev/null | grep '^rbceq2 2\.4\.4 (DB: ' \
    && test -s /opt/rbceq2/lib/python3.12/site-packages/rbceq2/resources/db.tsv
