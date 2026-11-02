# fastmixture 1.3.0 (github.com/Rosemeis/fastmixture tag v1.3.0; GPL-3.0: run, never
# redistributed: the image is built and run on this machine only; LICENSE at the tag and the
# sdist's pyproject.toml, GPL-3.0-only, read 2026-10-30), for supervised admixture proportions
# against the HGDP+1kGP panel. Neither upstream nor BioContainers publishes a linux/arm64 image,
# so it is built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs,
# each checked before it is used:
# - the base images by their index digests (Docker Hub, read 2026-10-20): python:3.12-bookworm
#   compiles fastmixture's Cython extensions with OpenMP, python:3.12-slim-bookworm runs them
#   (with the build image's libgomp, which the slim image lacks);
# - fastmixture's sdist from PyPI (the only file of 1.3.0, uploaded 2026-01-30) by its sha256;
# - its build tools and numpy by the sha256s in fastmixture.requirements.txt.
# Its setup.py compiles with -Ofast -march=native: the extensions target the CPU of the
# machine that builds the image, which is the machine that runs it.
# genomelab.containers names the image by the sha256 of this file and the list.
FROM python:3.12-bookworm@sha256:5560e9ab8709f459489e5b8aa696eda8a07ef821e14bb122be62d91234bfa98b AS build
COPY fastmixture.requirements.txt /tmp/
ADD --checksum=sha256:64c78ba77d0842257058ac9225f702afc68c3ff0d4cfba9426a92928fa248f46 \
    https://files.pythonhosted.org/packages/source/f/fastmixture/fastmixture-1.3.0.tar.gz /tmp/
RUN python -m venv /opt/fastmixture \
    && /opt/fastmixture/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/fastmixture.requirements.txt \
    && /opt/fastmixture/bin/pip install --no-cache-dir --no-deps --no-build-isolation \
        /tmp/fastmixture-1.3.0.tar.gz \
    && cp -L "$(gcc -print-file-name=libgomp.so.1)" /opt/fastmixture/libgomp.so.1

FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
COPY --from=build /opt/fastmixture /opt/fastmixture
RUN mv /opt/fastmixture/libgomp.so.1 /usr/local/lib/ && ldconfig
ENV PATH=/opt/fastmixture/bin:$PATH
# Offline, as every run: the extensions load and both modes the slice uses are offered.
RUN --network=none fastmixture --help > /tmp/help.txt \
    && grep -q -- '--supervised' /tmp/help.txt && grep -q -- '--projection' /tmp/help.txt \
    && python -c "from fastmixture import em, svd, shared" && rm /tmp/help.txt
