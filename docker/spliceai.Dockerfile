# SpliceAI 1.3.1 (github.com/Illumina/SpliceAI; its code GPL-3.0, its models and scores
# CC BY-NC 4.0: non-commercial use only), for slice 1e's splice predictions in the listed
# genes. Illumina publishes no image and BioContainers' has no TensorFlow for linux/arm64,
# so it is built here, for linux/arm64 (`genomelab.containers.PLATFORMS`), from pinned inputs,
# each checked before it is used:
# - the base image by its index digest (python:3.11-slim-bookworm, Docker Hub, read
#   2026-10-23; TensorFlow 2.15 runs on Python 3.9 to 3.11);
# - SpliceAI itself, its five models and GENCODE annotations inside its PyPI files, and every
#   dependency, by the sha256s in spliceai.requirements.txt.
# genomelab.containers names the image by the sha256 of this file and the list.
FROM python:3.11-slim-bookworm@sha256:f06d01d4ba3f8ad5af45a3ec87df19c1d317d8e8961edc0fabda61811af6f7cf
COPY spliceai.requirements.txt /tmp/
RUN python -m venv /opt/spliceai \
    && /opt/spliceai/bin/pip install --no-cache-dir --require-hashes --no-deps \
        -r /tmp/spliceai.requirements.txt
ENV PATH=/opt/spliceai/bin:$PATH
# Loads TensorFlow and SpliceAI's models, as a run does.
RUN python -c "from spliceai.utils import Annotator; from keras.models import load_model; from pkg_resources import resource_filename; load_model(resource_filename('spliceai', 'models/spliceai1.h5'))"
