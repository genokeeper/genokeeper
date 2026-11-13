# Changelog

All notable changes to Genokeeper are listed here. Versions follow [Semantic
Versioning](https://semver.org); until 1.0, a minor version may change commands and outputs.

## [0.2.0] — 2026-11-13

- A public demo on GIAB's HG002: three pages (an overview, the dashboard, the calibration
  tables) and `genokeeper demo` (`uvx genokeeper demo`, macOS or Linux), which downloads the
  demo store, checks it against its manifest and writes a folder Claude Code can answer
  questions about with the `genome-ask` skill. The maintainer-only `python -m genokeeper.demo
  build|package` makes it.
- The dashboard's two sliders become a 3x4 matrix of reliability and evidence thresholds; the
  view follows the URL hash. The first screen opens on what stands out, the Changes tab opens
  on the changes that exist when none survive your thresholds and names which source versions
  differ between the two runs, and a generated explanation shows a short gist above its
  paragraph. The page uses the system font and no web font.
- `genokeeper explain` takes `--min-reliability` and `--min-validity` to choose which findings
  it explains, `--model`, and `--tone engaging` (a gist and a paragraph for a curious reader);
  the default tone's prompts and cache are unchanged. It prints the calls made and what they
  cost. Explanations and the MedlinePlus Genetics summaries they are grounded in ship with the
  demo.
- Drug-response grading changed. Each drug result now quotes one curator's level: for a gene
  PharmCAT calls, the guideline's strength (CPIC's gene-drug level where it reports none); for
  any other gene, ClinPGx's level of the same variant (1B strong, 2A and 2B moderate); a
  ClinVar drug-response record with neither is limited. ClinVar's review stars no longer make
  one strong, and a variant ClinPGx annotates no longer shows as a second, stronger result. A
  drug result may therefore change grade at your next `genokeeper refresh`, because the
  refresh grades your calls again by the new rule.
- The store is found outside a checkout (the installed tool works in its working folder), and
  the `genome-ask` skill ships in the wheel.

## [0.1.0] — 2026-11-11

The first public version.

- From paired DNBSEQ-T7 FASTQ to GRCh38 calls on your own machine (tested on Apple-silicon
  Macs): bwa, a QC gate, DeepVariant, chrM with mutserve, read-backed phasing.
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
