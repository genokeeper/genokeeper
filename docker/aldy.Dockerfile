# Aldy 4.8.3 (github.com/0xTCG/aldy; Aldy/IURTC licence: non-commercial research only, no
# clinical or diagnostic use), for CYP2D6 and CYP2B6. BioContainers has no Aldy image, so
# it is built here from pinned inputs, each checked before it is used:
# - the base images by their index digests (Docker Hub, read 2026-10-20): python:3.12-bookworm
#   compiles Aldy's Cython extensions and mappy, python:3.12-slim-bookworm runs them;
# - Aldy's sdist from PyPI (the only file of 4.8.3, uploaded 2026-01-08) by its sha256;
# - every dependency by the sha256s in aldy.build-requirements.txt and aldy.requirements.txt.
# genomelab.containers names the image by the sha256 of this file and the two lists.
FROM python:3.12-bookworm@sha256:5560e9ab8709f459489e5b8aa696eda8a07ef821e14bb122be62d91234bfa98b AS build
COPY aldy.build-requirements.txt aldy.requirements.txt /tmp/
ADD --checksum=sha256:2612e849b35bc3c11fd0ee5752c8c732ae7d39fe095975d5249ea1666e867acd \
    https://files.pythonhosted.org/packages/source/a/aldy/aldy-4.8.3.tar.gz /tmp/
RUN python -m venv /opt/aldy \
    && /opt/aldy/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/aldy.build-requirements.txt \
    && /opt/aldy/bin/pip install --no-cache-dir --require-hashes --no-deps \
        --no-build-isolation -r /tmp/aldy.requirements.txt \
    && /opt/aldy/bin/pip install --no-cache-dir --no-deps --no-build-isolation \
        /tmp/aldy-4.8.3.tar.gz

FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
COPY --from=build /opt/aldy /opt/aldy
ENV PATH=/opt/aldy/bin:$PATH
RUN aldy --help > /dev/null
