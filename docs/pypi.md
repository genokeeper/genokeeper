# Genokeeper

A local, privacy-first analyser for your own whole-genome sequencing data. It takes the raw
paired reads of a 30x WGS run, aligns them to GRCh38 and calls them on your own machine, and
grades every finding with an error rate measured on public reference samples (GIAB's HG002,
HG003 and HG004). `genokeeper refresh` grades your calls again on the newest releases of ClinVar,
ClinGen, GenCC, CPIC, ClinPGx and the GWAS Catalog and shows what changed.

> **Research and educational use only.** Not intended for diagnosis, prevention, prediction,
> prognosis or treatment; not a medical device, not CE-marked. Confirm any result in a clinical
> laboratory.

**Status:** alpha (0.2.x), one maintainer, tested only on Apple-silicon Macs; until 1.0 a minor
version may change commands and outputs. **Use:** the code is AGPL-3.0-or-later, but six
components the bootstrap fetches are non-commercial only, so the pipeline as a whole is for
non-commercial use ([Licence](https://github.com/genokeeper/genokeeper#licence)).

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

To look at what it reports before you sequence anything, `uvx genokeeper demo` (macOS or
Linux) downloads a public demo on GIAB's HG002 into `~/genokeeper-demo`, a folder Claude Code
can answer questions about.

Requirements, supported input (DNBSEQ-T7 FASTQ in this version), privacy and every network
connection the tool makes are in the
[README](https://github.com/genokeeper/genokeeper#readme). Licence: AGPL-3.0-or-later.
