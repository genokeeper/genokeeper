# Genokeeper

[![CI](https://github.com/genokeeper/genokeeper/actions/workflows/ci.yml/badge.svg)](https://github.com/genokeeper/genokeeper/actions/workflows/ci.yml)
[![Licence: AGPL v3](https://img.shields.io/badge/licence-AGPL--3.0--or--later-blue.svg)](LICENSE)
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/16523/badge)](https://www.bestpractices.dev/projects/16523)

> **Not medical advice, not a medical device.** This is a research tool a person runs on their
> own genome; nothing it reports is a diagnosis or a basis for a medical decision.

A local analyser for your own whole-genome sequencing data. It takes the raw paired reads of a
30x WGS run, aligns them to GRCh38 and calls them on your Mac. It grades every finding with an
error rate measured on public reference samples (GIAB's HG002, HG003 and HG004). When you run
`genokeeper refresh`, it fetches the newest releases of ClinVar, ClinGen, GenCC, CPIC, ClinPGx
and the GWAS Catalog and grades your calls on them again. An optional daily check tells you
when a newer release is out. It never runs the refresh for you.

**Nothing leaves your machine unless you ask.** Genokeeper has no telemetry, no analytics and
no crash reporting. The pipeline never uploads your reads, genotypes or findings: every tool
that reads them runs in a container with no network, and a run opens no connection at all.
Two optional features send finding data to Anthropic, and only when you use them: the
`genome-ask` skill in your own Claude Code session, and `genokeeper explain`.
[Privacy and network](#privacy-and-network) lists every connection the tool makes and how to
check the claim yourself.

**Supported input (this version):** one pair of FASTQ files (`*.fq.gz` or `*.fastq.gz`) from an
original DNBSEQ-T7 run of about 30x, as tellmeGen delivers it, for example. Reads from another
sequencer, Illumina included, are refused with the reason, because the variant caller's model
and every calibration sample are T7's. BAM, CRAM and VCF input are not supported yet.

**Requirements:**
- a Mac with Apple silicon and 48 GB of memory or more: Docker needs 32 GB of it, and macOS
  needs the rest;
- about 1.5 TB of downloads, and 1 TB of free disk at the peak; about 0.85 TB stays after
  calibration, and about 0.5 TB once you free what may go ([Freeing disk](#freeing-disk) adds
  it up);
- several days of mostly unattended downloading and computing before the first graded result;
- optional: [Claude Code](https://claude.com/claude-code), for the `genome-ask` skill and
  `genokeeper explain`. Both use your own Claude Code login, so they count against your Claude
  plan's usage, or are billed per token if Claude Code runs on an API key.

[Quick start](#quick-start) takes a fresh Mac from nothing to the dashboard, step by step. The
maintainers' detail (each tool, each calibration measurement, how a run and a refresh behave)
is in [docs/reference.md](docs/reference.md).

## See it

Everything in this section is made from an invented genome, not anyone's results.

**Ask your agent.** The `genome-ask` skill answers questions about your genome in a Claude
Code session from `genokeeper query`. Every claim comes with how reliable the genotype is and how
strong the evidence behind it is. The skill never chooses a drug or a dose for you. Here is a
real, unedited session on the invented genome:

![A Claude Code session: the genome-ask skill answers which painkillers the invented genome's CYP2D6 and CYP2C9 results affect, with each result's reliability and evidence](docs/demo/genome-ask.gif)

Other questions it answers the same way:
- "Am I a carrier of anything my children could inherit?" It lists the carrier findings, how
  sure each genotype is, and the genes this run could not assess.
- "What do my results say about heart disease?" It gives the ClinVar findings and the
  polygenic scores, each score as a percentile with its population and its predictive power.
- "Do I have the factor V Leiden variant?" It reports the genotype at the site, from
  `genokeeper genotype`.
- "What changed since the last update?" It reports what the last refresh graded again, from
  `genokeeper changes`.

**The dashboard.** `genokeeper dashboard` writes one static page of every finding. It has a
slider for genotype reliability and one for strength of evidence, topic chips and a search box.
Below the list, the page counts what the sliders hide and why.

A whole genome yields thousands of candidate findings, and most rest on an uncertain genotype or
weak evidence. Each finding carries two measures:
- how likely its genotype is wrong, measured on reference samples;
- how strong the evidence behind its claim is.

With both sliders at zero, the invented genome shows 2,405 findings. At the default
thresholds, 18 remain.

An open row shows its genotype, its evidence and sources, and, once `genokeeper explain` has run,
a plain-language explanation. The page labels that explanation with the model that wrote it and
the sources it was given.

![The demo dashboard: both sliders rise from zero, the count falls from 2,405 findings to 18, then two rows open with explanations written by a language model](docs/demo/dashboard.gif)

To try the filters yourself, download the demo page
[docs/demo/dashboard.html](docs/demo/dashboard.html) and open it in a browser. It loads nothing
from the network.

![The demo dashboard at its default thresholds: findings with genotype-reliability and evidence bars, the topic chips and the two sliders](docs/demo/dashboard.png)

**Language.** Text is in English by default. `GENOKEEPER_LANG=ru`, or `--lang ru` on
`dashboard` and `explain`, switches the dashboard, the explanations and the update notices to
Russian. The skill answers in the language you ask in.

## Quick start

From a fresh Apple-silicon Mac to the dashboard and the agent skill, in the order that works.
The order matters. The reference samples are cut down to the number of bases your own reads
hold, so your FASTQ goes through a first `genokeeper run` before they are fetched. Times are
rough, for a recent Mac that meets the requirements; downloads take as long as your connection
makes them. Every download step can be rerun: it resumes a cut transfer and skips what is done.
Read [Terms you accept](#terms-you-accept) before step 5.

**1. Prerequisites.**
- Xcode's command line tools (`xcode-select --install`) and [Homebrew](https://brew.sh).
- uv and seqtk: `brew install uv seqtk`. The bootstrap installs neither. Nothing checks for
  seqtk until a reference sample is cut down (steps 7 and 8), and those steps stop without it.
- Rosetta (`softwareupdate --install-rosetta`): eight of the pinned images are published for
  Intel (linux/amd64) only.
- [Docker Desktop](https://www.docker.com/products/docker-desktop/), running for every step
  below. In Settings, Resources, give it 32 GB of memory or more. In Settings, General, keep
  "Use Rosetta for x86_64/amd64 emulation on Apple Silicon" on.
- Optional: [Claude Code](https://claude.com/claude-code), for step 10's skill and `explain`.

**2. Clone and install.** uv fetches the Python the project pins.

```bash
git clone https://github.com/genokeeper/genokeeper.git && cd genokeeper
uv sync
git config core.hooksPath .githooks
```

The hooks run the genome-data guard before every commit and push, so your data cannot be
committed by mistake.

**3. Choose where the lab goes.** The lab is the folder of tools, references and reference
samples; the bootstrap fills `~/genome-lab/`. To put it on another disk, set
`$GENOKEEPER_HOME` in your shell profile now. Every command reads it. Set it before step 5: the
lab records absolute paths, so moving it later is not supported.

```bash
echo 'export GENOKEEPER_HOME=/Volumes/Data/genome-lab' >> ~/.zshrc && source ~/.zshrc
```

Your run's own outputs (the store, the alignment, the dashboard) go to `reports/` inside the
checkout, so clone onto a disk with room for your FASTQ pair and about 100 GB more.

**4. Put your reads in `data/raw/` and check their platform.** Exactly two FASTQ files,
`*.fq.gz` or `*.fastq.gz`, from a DNBSEQ-T7 run. The run pairs them in name order, so read 1's
file must sort first (`…_1.fq.gz`, `…_2.fq.gz`). A delivery split into lanes is concatenated
into one pair first, the lanes in the same order on both sides, and the lane files kept outside
`data/raw/` (gzip files concatenate into a valid gzip file):

```bash
mkdir -p data/raw
cat L01_1.fq.gz L02_1.fq.gz > data/raw/genome_1.fq.gz
cat L01_2.fq.gz L02_2.fq.gz > data/raw/genome_2.fq.gz
```

Then check, in seconds and before any download, that the reads come from the platform this
version supports. With nothing else installed yet, the first run checks read 1's first 100,000
reads and stops:

```bash
uv run genokeeper run
```

On a DNBSEQ-T7 it goes on to say the gnomAD table of step 5 is missing, which is expected now.
Reads from another platform are refused with what they are instead (see
[Troubleshooting](#troubleshooting)), so you find out before any download.

If the delivery has a 23andMe-format genotype file (`*.txt`, not `*.vc.txt`), put it beside
the FASTQ: the run checks the sample's identity against it, and records "not checked" without
it.

**5. Bootstrap the lab.** It installs the command line tools through Homebrew, pulls or
builds the pinned container images, and downloads the references and knowledge bases. It
streams gnomAD's 0.88 TB once to build a table of a few GB. Hours at 100 MB/s, a day or more on
a slower line.

```bash
bash scripts/bootstrap.sh --check   # what is missing; changes nothing
bash scripts/bootstrap.sh --with-gnomad --with-vep --with-masked-reference --with-sv --with-typing
```

`--check` lists every missing file of every step, whatever `--with-…` flags you give it, and
exits 0. A line that ends "(run with --with-vep)" says which flag fetches the file; it reads
the same when you gave that flag. The lines for the reference samples (HG002, HG003, HG004 and
their truth sets) name no flag: they belong to step 7. It checks the container images only
while Docker is running.

Leave out `--with-hg002` and `--with-hg004` for now. They cut the samples down to your run's
bases, so before step 6 they only wait for that run (the bootstrap says so and goes on).

**6. First run: alignment and QC.** Step 5 must be complete: `run` refuses without the gnomAD
table that `--with-gnomad` builds, and names the flag. The run checks that the reads come from
a T7, aligns them with bwa to GRCh38 (about 5 to 6 hours), refuses a sample that fails the QC
gate, and records your run's bases. With no calibration table yet it then stops with exit code
2. That is expected: "this run has recorded its … bases, which `genokeeper calibrate`
downsamples HG003 to; the calibration table at … is not calibrated". The alignment is cached,
and step 9 does not repeat it.

```bash
uv run genokeeper run
```

**7. Fetch the reference samples.** HG003's FASTQ (147 GB, kept until step 8 cuts it down),
HG002's two libraries (300 GB read, one file of about 150 GB on disk at a time, about 170 GB of
cut-down reads kept), HG004's library (144 GB read whole, then cut down) and GIAB's truth sets
for the three:

```bash
bash scripts/bootstrap.sh --with-calibration --with-hg002 --with-hg004
```

Then align them, each like a run (hours apiece):

```bash
uv run python -m genokeeper.calibrate.hg002 align pcrfree  # HG002 PCR-free: aligned and called
uv run python -m genokeeper.calibrate.hg002 align pcr      # HG002 PCR: aligned only
uv run python -m genokeeper.calibrate.hg002 align hg004    # refused unless HG004 shows as female
uv run python -m genokeeper.calibrate.hg002 accept hg004   # accepted at SNP F1 0.99 or more
```

**8. Calibrate.** It builds every calibration table in order. The first, `small_variants`,
cuts HG003 down to your bases, deletes the 147 GB download, aligns and calls HG003, and scores
it against GIAB's truth. The others measure the structural, repeat, copy-number, sex-dosage and
typing calls on the aligned samples. A table whose inputs are missing is skipped with one line
naming the command that supplies them (exit code 2); a rerun builds only what is still stale.
Expect a day or more.

```bash
uv run genokeeper calibrate
```

A shorter path to a first result: in step 7 fetch only `--with-calibration` and skip the
`align` commands, then run `uv run genokeeper calibrate --only small_variants`. The run then
grades the small variants, and every structural, repeat, copy-number, sex-dosage and typing
statement reads "not calibrated" until the rest is done. It saves about 450 GB of downloads.

**9. Second run.** It reuses the alignment, then calls, phases, annotates and grades
everything into `reports/genomelab.duckdb`, printing the identity verdict and an accounting of
every source. Most of a day. SpliceAI is the slowest part, and its time grows with the number
of your variants in the genes it scores: the code measured one SpliceAI process at 1,224 s on
1,757 of HG003's sites, and four at 597 s (`src/genokeeper/annotate/vep.py`). A run starts one
SpliceAI process per three of its threads.

```bash
uv run genokeeper run
```

**10. Read the results.**

```bash
uv run genokeeper dashboard   # writes reports/dashboard.html; open it in a browser
uv run genokeeper query       # the findings at the default thresholds
```

For the agent skill, start Claude Code in the checkout (`claude`) and ask in plain words: "am
I a carrier of anything?" The skill sends the `genokeeper query` answers it reads to Anthropic
through your Claude Code session. The optional `genokeeper explain` sends findings to
Anthropic as well. No other step sends your data off the machine.

**11. The daily update check (optional).** It notifies you when a source has a newer release.
It never fetches or grades anything: you run `uv run genokeeper refresh` when you choose.

```bash
brew install terminal-notifier       # open its app once for the permission prompt
uv run genokeeper schedule install
```

Approve the agent once under System Settings → General → Login Items & Extensions → "Allow in
the Background".

## Commands

Each command prints its options with `--help`. The ones a new user meets:

| Command | Does |
|---|---|
| `run` | aligns, calls and grades `data/raw/`'s FASTQ pair into the store. `--fastq1` and `--fastq2` name other files (both or neither); `--clinvar38 PATH` grades on that ClinVar GRCh38 VCF instead of the newest in the lab; `--genotypes PATH` names the 23andMe-format file; `--sex male\|female` overrides the sex read from chrX depth; `--threads N` |
| `query`, `show`, `genotype` | read the store: the findings at the filters, one statement in full, the answer at a gene's or marker's sites; `--json` prints what the skill reads |
| `dashboard` | writes `reports/dashboard.html` |
| `explain` | asks Claude Code for a plain explanation of each finding not explained yet; `--dry-run` prints every prompt and the call count and sends nothing |
| `refresh` | fetches the newest knowledge-base releases and grades your calls on them; `--explain` also explains the new findings |
| `changes` | what changed between two snapshots (one per run). By default it compares the latest run with the previous refresh's run; `--base RUN_ID` compares with another snapshot, `--run RUN_ID` picks the newer one. It lists only what entered or left the trusted findings; `--all` lists every change |
| `check-updates` | the update check by hand; `--notify` sends the notification |
| `schedule install`, `schedule uninstall` | the daily update check |
| `calibrate` | rebuilds every stale calibration table; `--only NAME` builds one (`small_variants`, `hg002_cells`, `stage3`, `deletions`, `cnv`, `insertions`, `paralogs`, `foundation`, `slice3f`), `--report` prints what each measured |
| `build-faf` | builds the gnomAD table (`bootstrap.sh --with-gnomad` runs it). With no option it streams each chromosome the table lacks or whose source changed; `--chrom chr1` (repeatable) builds that chromosome; `--check` lists the chromosomes the table lacks and uses no network |

`run` and `query` print the run id (`run <id>`).

## Terms you accept

The repository's own code is AGPL (see [Licence](#licence)). The tools and data the bootstrap
fetches keep their own terms, and you accept them by running the bootstrap. No step asks you to
click through them, so read these first:

- **Aldy** calls CYP2D6 and CYP2B6. Its licence, from the Indiana University Research and
  Technology Corporation (IURTC), permits non-commercial research only. It excludes clinical
  trials and diagnostic use without IURTC's written consent, and forbids processing data for a
  third party or distributing the software. Running the lab on someone else's genome for them
  is processing for a third party. The bootstrap builds Aldy's image from `docker/`; after that,
  this prints the licence:

  ```bash
  docker run --rm --network none "$(uv run python -c 'from genokeeper.containers import IMAGES; print(IMAGES["aldy"])')" aldy license
  ```
- **SpliceAI's models and scores** are CC BY-NC 4.0: non-commercial use only. The bootstrap
  builds SpliceAI's image with them.
- **The PGS Catalog's scoring files** come under each score's own terms: cite the score's
  publication, follow its authors' restrictions and EBI's terms of use. Each score's terms are
  in `src/genokeeper/resources/pgs/scores.toml` (`licence`).
- **StellarPGx** (a third CYP2D6 caller) is MIT-licensed.
- **hmmix** (the Neanderthal share) ships the GPL-3.0 text in its wheel, while PyPI and its
  repository say MIT. It runs here for you and is never redistributed.
- **Docker Desktop** has its own terms, which it shows the first time it starts; you accept them
  there.

Other non-commercial terms apply to REVEL's scores, AlphaMissense's predictions, Scramble and
ClinPGx's annotations; [Licence](#licence) lists them, and
[docs/reference.md](docs/reference.md#building-the-store) gives every tool's terms.

## Keeping it up to date

**The knowledge bases.** Run `uv run genokeeper refresh` when you want new releases, for example
when the daily check notifies you. It fetches each source's newest file, runs `genokeeper run`
only when a grading input, a calibration table or Genokeeper's own code changed, and shows the result in the dashboard's Changes tab and in
`genokeeper changes`. It never re-aligns or re-calls your reads. The pinned tools and data
(gnomAD, the PGS Catalog, VEP, PharmCAT and others) are reported when newer, never applied: they
change only with a new version of Genokeeper.

**Genokeeper itself.**

```bash
git pull
uv sync
bash scripts/bootstrap.sh --with-gnomad --with-vep --with-masked-reference --with-sv --with-typing
uv run genokeeper calibrate
uv run genokeeper run
```

The bootstrap, run again with the flags you used, fetches or builds only what the new version
pins differently. A new version can change the pipeline (a tool, an image, a flag). Each
calibration table records the pipeline it was measured on, so a change makes it stale. `run`
then refuses while the small-variant table is stale, naming `genokeeper calibrate`; a stale
table for another kind of call leaves those statements unassessed. `genokeeper calibrate`
rebuilds only the stale tables and says "current" for the rest. If the change touched your own
steps, the run redoes only those; a changed aligner re-aligns.

## Freeing disk

Where the space goes. The numbers are approximate; each names its source.

| What | Where | About | Source |
|---|---|---|---|
| tools, references, knowledge bases (step 5) | the lab | 100 GB | the earlier README's estimate |
| the 37 container images | Docker's virtual disk | 50 GB | their sizes as `docker image inspect` reports them on a built lab |
| your FASTQ pair | `data/raw/` | its own size, F below | |
| your run: alignment, cached steps, store | `reports/` | 100 GB, 45 GB of it the alignment | the earlier README's estimate |
| HG003's download (steps 7 to 8) | the lab | 147 GB | the bootstrap |
| HG004's download (step 7, while it is cut down) | the lab | 144 GB | the bootstrap |
| HG002's cut-down reads | the lab | 170 GB | the bootstrap |
| HG003's and HG004's cut-down reads | the lab | F each | they hold your run's bases |
| four aligned reference libraries | the lab | 45 GB each | the earlier README's estimate |

HG002's 170 GB for two cut-down pairs suggests F is about 85 GB for a 30x run.

- **Peak:** in step 8, just before `calibrate` deletes HG003's download. On disk then: the lab,
  the images, your FASTQ and run, HG003's download and its cut-down reads, HG002's and HG004's
  cut-down reads, and three aligned libraries. That is about 700 GB plus 3F, so about 950 GB.
  Keep 1 TB free.
- **After calibration:** about 600 GB plus 3F, so about 850 GB.
- **After freeing what is listed below:** about 430 GB plus F, so about 0.5 TB.

May be deleted once `genokeeper calibrate` says every table is current:
- HG004's cut-down reads: `uv run python -m genokeeper.calibrate.hg002 clean hg004`.
- HG002's and HG003's cut-down reads, by hand: the lab's
  `calibration/hg002/pcrfree/downsampled/`, `calibration/hg002/pcr/downsampled/` and
  `calibration/hg003/downsampled/`. No command reads them again: the
  bootstrap fetches nothing for a library that has an alignment record, and a re-alignment
  writes the reads back from the library's alignment.
- After a re-alignment, `clean <library>` deletes what it replaced.
- Docker's build cache: `docker builder prune`. Keep the images: a run never pulls one.

Keep:
- `reports/work/`, above all your alignment. `refresh` and every later `run` reuse it; without
  it a run aligns again (about 5 to 6 hours).
- `reports/genomelab.duckdb`, the store: it holds every run's snapshot that `changes` compares.
- your FASTQ pair in `data/raw/`: every run checks it.
- the aligned reference libraries and their `aligned.json` records under the lab's
  `calibration/`: `calibrate` rebuilds a stale table from them.
- the rest of the lab.

## Uninstall

```bash
uv run genokeeper schedule uninstall          # removes the daily check, if you installed it
uv run python -c 'from genokeeper.containers import IMAGES; print("\n".join(IMAGES.values()))' \
  | xargs docker image rm                     # removes the pinned images
rm -rf "${GENOKEEPER_HOME:-$HOME/genome-lab}" # removes the lab
```

Then delete the checkout. It holds your FASTQ in `data/raw/` and your results in `reports/`, so
move out what you want to keep first. Remove the `GENOKEEPER_HOME` line from your shell profile
if you added one. The Homebrew formulas the bootstrap installed (bcftools, samtools, htslib,
bwa, mosdepth, openjdk, openjdk@21) and the ones you installed (uv, seqtk, terminal-notifier)
can go with `brew uninstall` if nothing else needs them.

## Troubleshooting

**Docker memory.** A run warns when Docker gives its containers less than 32 GiB, and the
polygenic scores' ancestry step can then fail out of memory. Raise it in Docker Desktop,
Settings, Resources. DeepVariant also sizes its work to Docker's memory, so less memory makes
it slower. Docker with 32 GB needs a Mac with 48 GB or more.

**Rosetta.** Eight images (DeepVariant, VerifyBamID2, mutserve, WhatsHap, OptiType, T1K,
bigBedToBed and StellarPGx) are Intel-only. They need Rosetta installed
(`softwareupdate --install-rosetta`) and Docker Desktop's "Use Rosetta for x86_64/amd64
emulation on Apple Silicon" on.

**Reads from another platform.** The run reads the flowcell id in the read names and the base
qualities. A DNBSEQ-G400 run (flowcell `V…`), a DNBSEQ-T7+ run (flowcell `ML…`, or qualities
binned to four levels), Illumina and other platforms are refused, naming what was found. This
version has no path for them: the variant caller's model and every calibration sample are the
original T7's.

**A QC refusal.** The run refuses before calling when a sample fails the QC gate: mean depth,
the share of bases at 20x, mapped reads, duplicates, insert size, contamination, Ti/Tv or the
het/hom ratio. The message names the metric. The bounds and their sources are in
`config/grading.toml` under `[qc]`. A sample that fails is a question for your sequencing
provider. Changing a bound is your own decision, and it changes what the next refresh re-runs.

**Disk running out.** Fetching, cutting down or re-aligning a reference sample checks the free
space first. It refuses with "… needs N GB free in …, M GB is", keeping 20 GB to spare. Other
steps fail with their tool's own error when the disk fills. Free space as
[Freeing disk](#freeing-disk) says, then run the same command again: downloads resume where they
stopped, and finished steps are not redone.

**`bootstrap --check` still lists the reference samples after step 5.** That is expected. It
lists every input of every step whatever flags you give it, and the reference samples (HG002,
HG003, HG004 and their truth sets) come in step 7. Their lines go once step 7 has fetched them.
HG002's and HG004's cut-down reads stay listed until a first run has recorded your bases.

## Privacy and network

Genokeeper has no telemetry, no analytics, no crash reporting and no account. The list below is
every connection it makes, followed by the evidence that the rest of it makes none.

**Offline by construction.** Once the bootstrap has built the lab, these commands work with the
network off: `genokeeper run`, `query`, `show`, `genotype`, `changes` and `dashboard`.
- Every tool that reads your reads or calls runs in a pinned container started with
  `--network none --pull never`. There is one place where a container starts:
  `run_container` in `src/genokeeper/containers.py`, tested in `tests/test_containers.py`.
- The pgsc_calc tasks that Nextflow starts get the same flag, with `NXF_OFFLINE` set, and each
  task's launch line is checked for the flag after the run (`src/genokeeper/prs/pgsc.py`).
- VEP runs with `--offline` on its local cache.
- `test_a_run_opens_no_network_connection` (in `tests/test_pipeline_invariants.py`) runs the
  whole pipeline on a mini genome with every name lookup and socket connection refused. The
  run completes and attempts none.
- The dashboard is one static page whose Content-Security-Policy is `default-src 'none'`, so it
  can load nothing. A PubMed link opens only when you click it.

**What goes out, and when:**

| Command | Goes to | Sends |
|---|---|---|
| `bash scripts/bootstrap.sh` | the download hosts of the tools and references ([docs/reference.md](docs/reference.md#building-the-store)) | requests for public files, each checked by its checksum; Homebrew runs with `HOMEBREW_NO_ANALYTICS=1`, and rtg-tools runs with crash reports and usage logging off |
| `genokeeper refresh` | the same public release hosts | requests for public files |
| `genokeeper check-updates`, and the daily check if you install it (`genokeeper schedule install`) | the release endpoints of ClinVar, ClinGen, GenCC, ClinPGx, the GWAS Catalog, MedlinePlus, LitVar, CGD, the PGS Catalog, gnomAD and CSpec, and GitHub's release API | version requests, with no sample, run or user identifier |
| `genokeeper explain` or `genokeeper refresh --explain` (opt-in) | Anthropic, through your own Claude Code login | for each finding and watch item (a borderline repeat, an SMN1 silent-carrier marker) not yet explained: its gene, variant, your genotype and the claim (`--dry-run` prints every prompt first) |
| the `genome-ask` skill, if you ask questions through Claude Code (opt-in) | Anthropic, through your own Claude Code session | the `genokeeper query`, `show`, `genotype` and `changes` answers the agent reads in order to answer you |

Nothing in the table runs on its own except the daily check, which you install yourself. The
last two rows are the only ones that send finding data, and only when you ask a question or run
`explain`.

**Without a language model.** Nothing needs one.
- The dashboard shows every finding with its sources, and leaves the explanation empty until
  you run `genokeeper explain`.
- `query`, `show` and `genotype` print the store's data as it is.

**Your data and git.** `data/raw/` and `reports/` never enter git: `.gitignore` excludes them,
and the genome-data guard (`scripts/check_no_genome_data.sh`, run by the hooks of step 2 and in
CI) refuses genomic files and data files outside the data folders. List your kit or sample ID,
one extended regular expression per line, in `.git/info/genome-private-ids` (never committed),
and the guard refuses it too. Claude Code sessions in the checkout cannot read `data/raw/`,
`reports/` or the lab's `calibration/` with their file tools; a lab moved by `$GENOKEEPER_HOME`
is outside that rule. [docs/reference.md](docs/reference.md#data) gives the details and limits.

**Outside Genokeeper.** Docker Desktop sends its own usage statistics unless you turn them off in
its settings.

**Check it yourself.** After the bootstrap, turn the network off and run `uv run genokeeper run`
and `uv run genokeeper dashboard`. Or watch the processes with an outbound firewall (LuLu, Little
Snitch).

## Licence

Genokeeper's own code is licensed under the GNU Affero General Public License, version 3 or any
later version ([LICENSE](LICENSE)). You may use, change and share it. If you share a changed
version, or offer it to others over a network (a hosted service included), you must make its
source available under the same licence. Running it for yourself requires nothing.
[NOTICE](NOTICE) lists the third-party data the repository carries and their terms. Among them
are Illumina's repeat catalog, STRchive's thresholds, ClinGen's and GenCC's data, MITOMAP, RBCeq2,
PharmCAT's example report, and the HLA table derived from PharmGKB.

The tools and data the bootstrap fetches keep their own licences
([docs/reference.md](docs/reference.md#building-the-store) sets them out), and each user
downloads them from their sources. Six are non-commercial only:
- Aldy;
- REVEL's scores;
- the SpliceAI models (CC BY-NC 4.0);
- AlphaMissense's predictions (CC BY-NC-SA 4.0);
- Scramble (CC BY-NC-SA 4.0);
- ClinPGx's current download (CC BY-NC-SA 4.0).

IPD-IMGT/HLA's `hla.dat` is CC BY-ND 4.0. T1K builds its index from it in the lab, and that
adaptation is never distributed. T1K is MIT-licensed and OptiType BSD-3-Clause.

So the pipeline as a whole is for non-commercial use, whatever the AGPL allows for the code. The
container images built from `docker/` must not be published, because Aldy's licence forbids
distributing it.

Not medical advice and not a medical device.
