# Changelog

All notable changes to Genokeeper are listed here. Versions follow [Semantic
Versioning](https://semver.org); until 1.0, a minor version may change commands and outputs.

## [0.1.0] — unreleased

The first public version.

- From paired DNBSEQ-T7 FASTQ to GRCh38 calls on your Mac: bwa, a QC gate, DeepVariant, chrM
  with mutserve, read-backed phasing.
- Two confidences for every finding: a genotype error rate measured on GIAB's HG002, HG003 and
  HG004, and the strength of the evidence behind the claim.
- What it reports:
  - ClinVar findings and carrier status;
  - variants ClinVar has not classified, graded as predictions;
  - pharmacogenomics (PharmCAT, CYP2D6, CYP2B6, HLA);
  - polygenic scores;
  - mitochondrial and Y haplogroups;
  - SMN1/SMN2 copy number and repeat expansions;
  - deletions, copy number, insertions and paralog loci;
  - blood groups, ancestry, traits.
- A static dashboard with reliability and evidence filters.
- An agent skill that answers questions about your genome.
- Model-written explanations, opt-in.
- A monthly refresh that grades the findings again against new releases and reports what changed.
- No telemetry; a run opens no network connection (README, "Privacy and network").
- English by default; `GENOKEEPER_LANG=ru` for Russian.
