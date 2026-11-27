# Disk, uninstall and troubleshooting

How much disk a lab takes and what may be deleted, how to remove Genokeeper, and what to do when a
step refuses or fails. The step numbers are those of the README's
[Quick start](../README.md#quick-start).

## Freeing disk

Where the space goes. The numbers are approximate; each names its source.

| What | Where | About | Source |
|---|---|---|---|
| tools, references, knowledge bases (step 5) | the lab | 100 GB | a rough estimate, not measured by the code |
| the 37 container images | Docker's virtual disk | 50 GB | their sizes as `docker image inspect` reports them on a built lab |
| your FASTQ pair | `data/raw/` | its own size, F below | |
| your run: alignment, cached steps, store | `reports/` | 100 GB, 45 GB of it the alignment | a rough estimate, not measured by the code |
| HG003's download (steps 7 to 8) | the lab | 147 GB | the bootstrap |
| HG004's download (step 7, while it is cut down) | the lab | 144 GB | the bootstrap |
| HG002's cut-down reads | the lab | 170 GB | the bootstrap |
| HG003's and HG004's cut-down reads | the lab | F each | they hold your run's bases |
| four aligned reference libraries | the lab | 45 GB each | a rough estimate, not measured by the code |

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
  `calibration/hg003/downsampled/`. While a library keeps its alignment record nothing reads
  them: the bootstrap fetches nothing for such a library, and a re-alignment writes the reads
  back from its alignment.
- After a re-alignment, `clean <library>` deletes the old alignment and the reads the new one
  was made from.
- Docker's build cache: `docker builder prune`. Keep the images: a run never pulls one.

Keep:
- `reports/work/`, above all your alignment. `refresh` and every later `run` reuse it; without
  it a run aligns again (about 5 to 6 hours, a rough estimate, not measured by the code).
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

**Docker memory.** A run warns when Docker gives its containers less than 36 GiB: the
polygenic scores' matching steps (pgsc_calc's FORMAT_SCOREFILES and MATCH_COMBINE) peaked at
34.6 GB on a 30x genome and can fail out of memory below it. Raise it in Docker Desktop,
Settings, Resources. DeepVariant sizes its work to Docker's memory, so less memory only makes
it slower. Docker with 36 GB needs a Mac with 48 GB or more; running the scores in smaller
batches to lower this is tracked in [#50](https://github.com/genokeeper/genokeeper/issues/50).

**Rosetta.** Eight images (DeepVariant, VerifyBamID2, mutserve, WhatsHap, OptiType, T1K,
bigBedToBed and StellarPGx) are Intel-only. They need Rosetta installed
(`softwareupdate --install-rosetta`) and Docker Desktop's "Use Rosetta for x86_64/amd64
emulation on Apple Silicon" on.

**Reads from another platform.** The run reads the flowcell id in the read names and the base
qualities. A DNBSEQ-G400 run (flowcell `V…`), a DNBSEQ-T7+ run (flowcell `ML…`, or qualities
binned to four levels), Illumina and other platforms are refused, naming what was found. This
version has no path for them: the variant caller's model and every calibration sample are the
original T7's.

**A QC refusal.** The QC gate judges a sample twice. After alignment, before calling, it reads
mean depth, the share of bases at 20x, mapped reads, duplicates, insert size and contamination;
a refusal there comes in step 6. After DeepVariant has called, it reads Ti/Tv and the het/hom
ratio as well; a refusal on those comes only in step 9, because the first run stops before
calling. The message names the metric. The bounds and their sources are in
`config/grading.toml` under `[qc]`. A sample that fails is a question for your sequencing
provider. Changing a bound is your own decision, and it changes what the next refresh re-runs.

**Disk running out.** Three steps check the free space first and refuse with "… needs N GB free
in …, M GB is", keeping 20 GB to spare: downloading HG002's and HG004's libraries, cutting
HG002's down, and re-aligning a reference library from its alignment. The others do not check:
HG003's 147 GB download, cutting HG003 and HG004 down, and every other step fail with their
tool's own error when the disk fills. Free space as
[Freeing disk](#freeing-disk) says, then run the same command again: downloads resume where they
stopped, and finished steps are not redone.

**`bootstrap --check` still lists the reference samples after step 5.** That is expected. It
lists the inputs of every step whatever `--with-…` flags you give it (only LitVar2's file and
ClinVar's citations wait for `--with-publications`), and the reference samples (HG002,
HG003, HG004 and their truth sets) come in step 7. Their lines go once step 7 has fetched them.
HG002's and HG004's cut-down reads stay listed until a first run has recorded your bases.
