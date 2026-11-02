# BCyrius 2.2.0 (gitlab.com/andreassh/bcyrius; Apache-2.0; bundles Beagle 5.5, GPL-3, and
# its CYP2D6 phasing panel), for CYP2D6. Neither BioContainers nor upstream publishes an
# image, so it is built here from pinned inputs, each checked before it is used:
# - the base image by its index digest (python:3.12-slim-bookworm, Docker Hub, read
#   2026-10-20);
# - the source by its git commit (tag v2.2.0), which names its content;
# - every dependency, with the jdk4py Java runtime Beagle runs on (its `java` extra), by
#   the sha256s in bcyrius.requirements.txt.
# genomelab.containers names the image by the sha256 of this file and the list.
FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
COPY bcyrius.requirements.txt /tmp/
ADD --keep-git-dir=false \
    https://gitlab.com/andreassh/bcyrius.git#6ee6edfb73cf9c00f62ad1d6e37d10cd80a34624 \
    /tmp/bcyrius
RUN python -m venv /opt/bcyrius \
    && /opt/bcyrius/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/bcyrius.requirements.txt \
    && /opt/bcyrius/bin/pip install --no-cache-dir --no-deps --no-build-isolation \
        /tmp/bcyrius \
    && rm -rf /tmp/bcyrius
ENV PATH=/opt/bcyrius/bin:$PATH
RUN bcyrius --help > /dev/null
