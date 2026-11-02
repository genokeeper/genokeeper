# Haplogrep3 3.3.2 (github.com/genepi/haplogrep3, MIT licence) with PhyloTree Build 17's rCRS
# tree, phylotree-rcrs@17.3 (genomelab.haplogroup.TREE), for the mitochondrial haplogroup.
# Haplogrep3 publishes no image and fetches its trees from the network on first use, which
# the run's --network none forbids: so the tree is baked in here, and haplogrep3.yaml names
# no repository. Pinned inputs, each checked before it is used:
# - the base images by their index digests (Docker Hub, read 2026-10-21): busybox 1.37.0
#   unpacks, eclipse-temurin:21-jre-noble runs Haplogrep3's jar;
# - Haplogrep3's release tarball and the tree's release zip by the sha256 digests GitHub's
#   release API gives for those assets (v3.3.2, published 2026-06-09; 17.3, 2026-06-08).
# genomelab.containers names the image by the sha256 of this file and haplogrep3.yaml.
FROM busybox:1.37.0@sha256:bdf57e528e45e4433820e045b29b4597825a1c9e38353532d90a01445013f82e AS unpack
ADD --checksum=sha256:957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959 \
    https://github.com/genepi/haplogrep3/releases/download/v3.3.2/haplogrep3-3.3.2-linux.tar.gz /tmp/
ADD --checksum=sha256:5278133257f519fb94572eb9e6266ea24b34546cab32f26587700291e0103736 \
    https://github.com/genepi/phylotree-rcrs-17/releases/download/17.3/phylotree-rcrs-17.3.zip /tmp/
COPY haplogrep3.yaml /opt/haplogrep3/
RUN tar -xzf /tmp/haplogrep3-3.3.2-linux.tar.gz -C /opt/haplogrep3 haplogrep3.jar LICENSE \
    && mkdir -p /opt/haplogrep3/trees/phylotree-rcrs/17.3 \
    && unzip -q /tmp/phylotree-rcrs-17.3.zip -d /opt/haplogrep3/trees/phylotree-rcrs/17.3

FROM eclipse-temurin:21-jre-noble@sha256:35501d3f93173b420ac9001658eb68157e813316b865087e345c187fca422c4c
COPY --from=unpack /opt/haplogrep3 /opt/haplogrep3
# Offline, as every run: the image lists its one tree without a repository.
RUN --network=none java -jar /opt/haplogrep3/haplogrep3.jar trees | grep -q 'phylotree-rcrs@17.3'
