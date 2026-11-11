# Contributing to Genokeeper

Thank you for helping. Genokeeper reads people's genomes, so one rule comes before all others.

## Never share genome data

No contribution, issue, pull request, commit message or log may contain a person's genome data:
genotypes, rsIDs with their alleles, findings, haplogroups, counts from a person's run, the output
of `genokeeper query`, `show`, `genotype`, `changes` or `dashboard`, or file names that carry a
kit or sample ID. Test data comes from a public sample (GIAB's HG002, HG003 and HG004, 1000
Genomes) or is invented, and says which. The genome-data guard (`scripts/check_no_genome_data.sh`)
refuses the common mistakes; it cannot read what a sentence says.

## Set up

```bash
git clone https://github.com/genokeeper/genokeeper.git && cd genokeeper
uv sync
git config core.hooksPath .githooks
```

The hooks run the genome-data guard before every commit and push. Developing needs no genome: the
test suite runs on fixtures it builds at test time. It does run some real tools, the ones CI
installs (`.github/workflows/ci.yml`):

- bcftools 1.24, samtools 1.24, htslib 1.24 (bgzip, tabix), seqtk 1.5, bwa 0.7.19 and mosdepth
  0.3.14 on `PATH`: `brew install bcftools samtools htslib seqtk bwa mosdepth` on a Mac;
- rtg-tools 3.13 at `<lab>/tools/rtg-tools-3.13/rtg`, where `<lab>` is `$GENOKEEPER_HOME` or
  `~/genome-lab`, with an `rtg.cfg` that turns off its crash reports and usage logging.
  `bash scripts/bootstrap.sh` installs it there (with everything else); CI unpacks the release
  zip by hand. Without it, the five tests in `tests/test_calibrate_measure.py` that run
  `rtg vcfeval` error;
- WhatsHap's pinned image in Docker, for the phasing tests (CI pulls it at the digest
  `genokeeper.containers` pins).

Running the pipeline itself needs the whole bootstrap
([docs/reference.md](docs/reference.md#building-the-store)).

To keep private IDs (your kit or sample ID) out of every commit, list them, one extended regular
expression per line, in `.git/info/genome-private-ids`. That file is never committed.

## Before a pull request

Run what CI runs:

```bash
uv run basedpyright
uv run ruff check
uv run ruff format --check
uv run pytest -q -n auto
bash scripts/check_no_genome_data.sh --tracked
```

- Keep a pull request to one change, with a test that fails without it.
- Write code, comments, commit messages and pull requests in English.
- A new download is pinned by checksum: a sha256 written in the script, or the md5 its
  publisher lists when that is the only digest it gives (NCBI's, ClinVar's, Zenodo's and the
  cloud buckets' listings), never an md5 computed here. HTTPS carries every download. A new tool
  runs in a container started with
  `--network none` (`run_container` in `src/genokeeper/containers.py`). Anything that sends data off
  the machine needs a row in README, "Privacy and network", and a reason it cannot run locally.
- A new data source or tool states its licence in `docs/reference.md`; terms a user must
  accept (non-commercial ones among them) also go in the README's "Terms you accept" and
  "Licence"; and if the repository carries any of its data, in NOTICE.

## Reviewing

Maintainers review a pull request from outside in a checkout with no genome data and no lab: the
git hooks, `uv run genokeeper` and an agent's project settings all run the checked-out code.

## Licence

Genokeeper is licensed under the GNU AGPL, version 3 or any later version. By contributing, you
agree that your contribution is licensed the same way.

## Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
