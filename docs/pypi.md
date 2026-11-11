# Genokeeper

A local, privacy-first analyser for your own whole-genome sequencing data. It takes the raw
paired reads of a 30x WGS run, aligns them to GRCh38 and calls them on your own machine, and
grades every finding with an error rate measured on public reference samples (GIAB's HG002,
HG003 and HG004). `genokeeper refresh` grades your calls again on the newest releases of ClinVar,
ClinGen, GenCC, CPIC, ClinPGx and the GWAS Catalog and shows what changed.

> **Not medical advice, not a medical device.** This is a research tool a person runs on their
> own genome; nothing it reports is a diagnosis or a basis for a medical decision.

**Install from the repository, not from this package alone.** `pip install genokeeper` gives
you the command-line tool, but building the lab it runs on (the pinned containers, reference
data and calibration samples) needs `scripts/bootstrap.sh` and the `docker/` folder of the
source checkout. Follow the
[Quick start](https://github.com/genokeeper/genokeeper#quick-start) instead:

```bash
git clone https://github.com/genokeeper/genokeeper.git
cd genokeeper
uv sync --locked
```

Requirements, supported input (DNBSEQ-T7 FASTQ in this version), privacy and every network
connection the tool makes are in the
[README](https://github.com/genokeeper/genokeeper#readme). Licence: AGPL-3.0-or-later.
