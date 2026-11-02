# ExpansionHunter, built from Illumina's source (github.com/Illumina/ExpansionHunter) at commit
# 9e93322f4ae72d450a277c43fc6dd54e5b610a95 (2022-02-02, "Motif analysis updates"), for the
# repeat-expansion loci and the RFC1 motif analysis (ehunter/locus/RFC1MotifAnalysis.cpp), read
# via the GitHub API 2026-10-30. At that commit LICENSE.txt is Apache-2.0; the licence became
# PolyForm Strict at 0c9864d229 (2022-12-08), which is not used here. The binary still reports
# `ExpansionHunter v5.0.0`: the commit did not bump the version. Neither alternative serves:
# bioconda's expansionhunter 5.0.0 (and the v5.0.0 tag) lack the analysis, and the bw2 fork's
# branch moves daily. BioContainers' image is amd64 only, so this is built for linux/arm64
# (`genomelab.containers.PLATFORMS`) from pinned inputs, each checked before it is used:
# - the base image by its index digest (mambaorg/micromamba 2.8.1, Docker Hub, read 2026-10-26;
#   2.9.0's tag was re-pushed on 2026-10-24, so 2.8.1 is the newest tag over two weeks old);
# - the toolchain and the libraries it links by expansionhunter.lock, an explicit lockfile
#   (`micromamba create --dry-run` for linux-aarch64 inside that image, conda-forge, 2026-10-30).
#   It pins cmake 3.31: CMake 4 refuses the cmake_minimum_required below 3.5 that spdlog
#   v1.6.1 and googletest 1.10.0 declare;
# - the source by the commit above, and the five tarballs the commit's CMakeLists.txt would
#   download by ExternalProject_Add, each by the sha256 of its bytes (computed on the download
#   2026-10-30; ADD refuses other bytes), so that the build itself runs offline. Boost 1.73.0 is
#   fetched from archives.boost.io, because the boostorg.jfrog.io URL the commit names now
#   redirects to a JFrog reactivation page; its sha256 is the digest Boost publishes for
#   boost_1_73_0.tar.bz2, on https://www.boost.org/releases/1.73.0/ and in
#   https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2.json (read 2026-10-30).
#   htslib 1.22.1 is from its GitHub release; spdlog
#   v1.6.1, googletest release-1.10.0 and abseil 20210324.2 are GitHub archives.
# - htslib is 1.22.1, not the one the commit's CMakeLists.txt names (1.10.2): every CRAM in
#   this lab is version 3.1, which htslib reads only from 1.12 on (the older one makes
#   ExpansionHunter abort with `[E::cram_next_slice] Failure to decode slice`). 1.22.1 is the
#   htslib that bioconda's ExpansionHunter 5.0.0 was linked against for 3a, and bioconda's
#   recipe (recipes/expansionhunter/ in bioconda-recipes, build.sh, meta.yaml and its patches,
#   read 2026-10-30) is the model for pointing the build at a newer htslib: ExpansionHunter's
#   own CMake builds it, not a system one. Its sha256 is the one bioconda's htslib 1.22.1
#   recipe pins (recipes/htslib/meta.yaml at ccef06ec45, merged 2025-07-14, read 2026-10-30),
#   which equals the sha256 of the bytes downloaded from the GitHub release 2026-10-30. That
#   anchor comes from bioconda's recipe repository, which a tampered release asset cannot
#   move because the recipe was merged in 2025 and is pinned here by its commit.
# Two stages: the build stage compiles and checks ExpansionHunter; the final stage, on the
# same base image, carries only the binary and the libraries it links from the conda
# environment (the ones `ldd` resolves under /opt/conda; the rest are the base image's own).
# genomelab.containers names the image by the sha256 of this file and the lockfile.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4 AS build
COPY --chown=$MAMBA_USER:$MAMBA_USER expansionhunter.lock /tmp/
RUN micromamba create --yes --name expansionhunter --file /tmp/expansionhunter.lock \
    && micromamba clean --all --yes
ENV ENV_NAME=expansionhunter PATH=/opt/conda/envs/expansionhunter/bin:$PATH
# The binary links the environment's libcurl, OpenSSL and libstdc++ dynamically.
ENV LD_LIBRARY_PATH=/opt/conda/envs/expansionhunter/lib
ADD --keep-git-dir=false --chown=$MAMBA_USER:$MAMBA_USER https://github.com/Illumina/ExpansionHunter.git#9e93322f4ae72d450a277c43fc6dd54e5b610a95 /opt/eh-src
ADD --checksum=sha256:3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7 --chown=$MAMBA_USER:$MAMBA_USER \
    https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 \
    /opt/eh-deps/htslib-1.22.1.tar.bz2
ADD --checksum=sha256:4eb3b8d442b426dc35346235c8733b5ae35ba431690e38c6a8263dce9fcbb402 --chown=$MAMBA_USER:$MAMBA_USER \
    https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2 \
    /opt/eh-deps/boost_1_73_0.tar.bz2
ADD --checksum=sha256:378a040d91f787aec96d269b0c39189f58a6b852e4cbf9150ccfacbe85ebbbfc --chown=$MAMBA_USER:$MAMBA_USER \
    https://github.com/gabime/spdlog/archive/refs/tags/v1.6.1.tar.gz \
    /opt/eh-deps/spdlog-v1.6.1.tar.gz
ADD --checksum=sha256:9dc9157a9a1551ec7a7e43daea9a694a0bb5fb8bec81235d8a1e6ef64c716dcb --chown=$MAMBA_USER:$MAMBA_USER \
    https://github.com/google/googletest/archive/refs/tags/release-1.10.0.tar.gz \
    /opt/eh-deps/googletest-release-1.10.0.tar.gz
ADD --checksum=sha256:59b862f50e710277f8ede96f083a5bb8d7c9595376146838b9580be90374ee1f --chown=$MAMBA_USER:$MAMBA_USER \
    https://github.com/abseil/abseil-cpp/archive/refs/tags/20210324.2.tar.gz \
    /opt/eh-deps/abseil-20210324.2.tar.gz
# Each ExternalProject_Add URL in CMakeLists.txt is pointed at the local file; the guard fails
# the build if any remote URL is left. `make` also runs the commit's unit tests (the ehunter
# project's TEST_COMMAND ctest), and the binary goes on PATH. The micromamba image's RUN shell
# activates the environment (which sets CC and CXX to its compilers) only for this build.
# The conda toolchain names its tools aarch64-conda-linux-gnu-*, while Boost.Build's bootstrap
# and htslib's Makefile look for gcc and ar by their plain names, so those are linked on PATH.
WORKDIR /opt/eh-src
ENV MAMBA_DOCKERFILE_ACTIVATE=1
RUN --network=none sed -i \
        -e 's|https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2|/opt/eh-deps/htslib-1.22.1.tar.bz2|' \
        -e 's|https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2|/opt/eh-deps/boost_1_73_0.tar.bz2|' \
        -e 's|https://github.com/gabime/spdlog/archive/refs/tags/v1.6.1.tar.gz|/opt/eh-deps/spdlog-v1.6.1.tar.gz|' \
        -e 's|https://github.com/google/googletest/archive/refs/tags/release-1.10.0.tar.gz|/opt/eh-deps/googletest-release-1.10.0.tar.gz|' \
        -e 's|https://github.com/abseil/abseil-cpp/archive/refs/tags/20210324.2.tar.gz|/opt/eh-deps/abseil-20210324.2.tar.gz|' \
        CMakeLists.txt \
    && ! grep -n 'https\?://' CMakeLists.txt \
    && mkdir /tmp/eh-tools \
    && for t in gcc g++ cc c++ ar ranlib; do ln -s "$(command -v "aarch64-conda-linux-gnu-$t")" "/tmp/eh-tools/$t"; done \
    && export PATH="/tmp/eh-tools:$PATH" \
    && mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make -j"$(nproc)" \
    && cp install/bin/ExpansionHunter /opt/conda/envs/expansionhunter/bin/
# Check, offline: the htslib just built writes a CRAM 3.1 file from a SAM file in its own
# tarball (test/ce#5.sam on test/ce.fa) and reads it back; the build fails if any of a record's
# 11 mandatory SAM fields (name to quality) differs (the reader adds MD and NM tags of its own).
# The first six bytes must be `CRAM` 3 1, so that a writer which fell back to 3.0 cannot pass.
# Its anchor comes from htslib's own bundled test data, inside the sha256-pinned tarball above,
# which a later edit of this file cannot move because ADD refuses other bytes. No CRAM or BAM
# is committed: the file is made here, in /tmp.
RUN --network=none export PATH="/tmp/eh-tools:$PATH" \
    && cd /opt/eh-src/build/htslib-prefix/src/htslib \
    && make test/test_view \
    && test/test_view -S -t test/ce.fa -C -o VERSION=3.1 'test/ce#5.sam' > /tmp/ce5.cram \
    && test "$(head -c 6 /tmp/ce5.cram | od -An -tu1 | tr -s ' ')" = " 67 82 65 77 3 1" \
    && test/test_view -D -t test/ce.fa /tmp/ce5.cram | grep -v '^@' | cut -f1-11 > /tmp/ce5.out \
    && grep -v '^@' 'test/ce#5.sam' | cut -f1-11 | cmp - /tmp/ce5.out \
    && test -s /tmp/ce5.out && echo "CRAM 3.1 round trip: $(wc -l < /tmp/ce5.out) records" \
    && rm /tmp/ce5.cram /tmp/ce5.out
ENV MAMBA_DOCKERFILE_ACTIVATE=0
# Offline, as every run. The binary is the one with the RFC1 analysis, links everything, and
# carries the version string of the htslib above (a standalone `1.22.1`, as hts_version()
# returns it), so it was not linked against another one. (`grep -c` reads to the end: under
# pipefail a `grep -q` that exits early would fail the pipe by SIGPIPE.)
RUN --network=none ExpansionHunter --version 2>&1 | grep -q 'v5\.0\.0$' \
    && grep -aq RFC1MotifAnalysis "$(command -v ExpansionHunter)" \
    && test "$(tr '\0' '\n' < "$(command -v ExpansionHunter)" | grep -cxF 1.22.1)" -ge 1 \
    && ! ldd "$(command -v ExpansionHunter)" | grep 'not found'
# The binary and the conda libraries it links, for the final stage (under /tmp: the build
# stage runs as the micromamba image's user, who cannot write /opt).
RUN --network=none mkdir -p /tmp/eh/bin /tmp/eh/lib \
    && cp "$(command -v ExpansionHunter)" /tmp/eh/bin/ \
    && ldd /tmp/eh/bin/ExpansionHunter | awk '$3 ~ /^\/opt\/conda\// {print $3}' > /tmp/eh-libs \
    && test -s /tmp/eh-libs \
    && xargs cp -L -t /tmp/eh/lib/ < /tmp/eh-libs

FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4
COPY --from=build /tmp/eh/ /opt/eh/
ENV PATH=/opt/eh/bin:$PATH LD_LIBRARY_PATH=/opt/eh/lib
# Offline, as every run: the binary finds every library it links, and runs. (`--help`'s
# exit status is not relied on: its text must name `--reads`.)
RUN --network=none ! ldd /opt/eh/bin/ExpansionHunter | grep 'not found' \
    && { ExpansionHunter --help || true; } > /tmp/eh-help.txt 2>&1 \
    && grep -q -- '--reads' /tmp/eh-help.txt \
    && rm /tmp/eh-help.txt
