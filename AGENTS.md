# genokeeper

A local, self-updating analyser for a person's own WGS 30x data (DNBSEQ-T7 FASTQ, re-aligned
to GRCh38 here): open-source tools run on the user's Mac, a findings store keeps one
snapshot per run, and a monthly refresh reports what changed. Personal genome data never
enters git and never leaves the machine unasked (`scripts/check_no_genome_data.sh`).

Keep this file within its budgets — at most 300 lines and 3,000 words, and the status section
below at most 50 lines — because it is loaded into every session; detail belongs in the
documents it links to.

## Current status

- Repository, guards and stayfixed adoption in place.
- The findings track is delivered: `genokeeper run` goes from the user's FASTQ to GRCh38
  calls (bwa, a QC gate, DeepVariant, mutserve, phasing), error rates calibrated on HG003,
  ClinVar claims, pharmacogenomics, polygenic scores, the haplogroups, GWAS research
  associations and, since slice 1e, VEP-predicted findings. `genokeeper dashboard` and
  `genokeeper explain` show them. Build steps: [README.md](README.md#building-the-store).
- Stage 4 (refresh): `genokeeper refresh` fetches the newest releases and grades the calls on
  them again, and the Changes tab and `genokeeper changes` show what changed
  ([README.md](README.md#refreshing)); `genokeeper check-updates`, run daily by
  `genokeeper schedule install`'s LaunchAgent, writes `~/genome-lab/refresh/check.json` and
  notifies once per newer set.
- Stage 3, slice 3a delivered (2026-10-28): SMN1/SMN2 copy number (SMNCopyNumberCaller and
  Parascopy), eleven repeat loci (ExpansionHunter on Illumina's catalog regions, STRchive's
  thresholds, REViewer images) and the Y haplogroup (Yleaf), each gated by a stage-3 table
  measured on HG002's two libraries (`genokeeper calibrate --only stage3`).
- Stage 3, slice 3b delivered (2026-10-29): deletions of 50 bp or more in the ACMG SF and
  carrier genes' exons from DELLY and Dysgu (each run one chromosome at a time), matched by
  Truvari, each class's chance of being false measured on HG002's two libraries
  (`genokeeper calibrate --only deletions`); every shown deletion needs confirmation by MLPA, a
  microarray, qPCR or long reads.
- Stage 3, slice 3c delivered (2026-11-01): GRCh38 with the GRC false duplications masked,
  `realign` for the calibration libraries, HG004, HG002 cells for chrX, chrY and v5.0q-only regions
  (their own table), and four run-time checks (sample identity, per-gene coverage, sex dosage, an
  allele-fraction guard). The `foundation` table of 2026-11-05 on HG002 PCR-free, HG002 PCR, HG003
  and HG004 sets the sex-dosage and allele-fraction bounds; identity and coverage use design
  thresholds.
- Slice 3g (2026-10-30): 34 more polygenic scores by 1d's rule read strictly, 24 conditions
  named not covered, a weak-power wording; six single-gene traits (2026-11-02: FUT2, ABCC11,
  MC1R, HERC2, AMPD1, CCR5Δ32).
- Stage 3, slice 3d calibrated (2026-11-06): depth copy number (DELLY cnv, CNVpytor),
  insertions (INSurVeyor, Dysgu, Scramble), paralog loci (Parascopy, Gauchian for GBA1, KILDA for
  Lp(a)) and the Tier 4 carrier list. AMY1 is withheld and a mobile element only Scramble called
  is not measured, each unassessed with its value under `claim.withheld`; no measured class
  reads "likely" (README.md, slice 3d's limits).
- Stage 3, slice 3e delivered (2026-11-01; `stage3` table of 2026-11-06 on HG002, HG003, HG004):
  one ExpansionHunter build (Illumina commit 9e93322f4a, htslib 1.22.1) over 24 loci, RFC1 and
  CFTR's poly-T/TG, an exclusion-only rule at the new loci, FMR1 and AFF2 by the X dosage, 22
  loci named not covered. The table withholds 12 items: not graded, each value kept for the
  agent ([which and why](README.md#building-the-store)).
- Stage 3, slice 3f delivered, its `slice3f` table built (2026-11-06, the GIAB trio): HLA, MITOMAP's
  chrM variants, HHV-6, runs of homozygosity, blood groups, ancestry, traits. Withheld, so not graded,
  the value kept: ABO, C4A, C4B, CTL2, FUT3, MAM, PIGG, XG (measurable on no library) and chrM
  heteroplasmy. UGT1A1\*28 is trusted; PharmCAT's read-backed phase sets are not (BR-036).
- Roadmap: [docs/roadmap.md](docs/roadmap.md).
- Questions about the user's genome: the `genome-ask` skill
  ([.claude/skills/genome-ask/SKILL.md](.claude/skills/genome-ask/SKILL.md)) answers from
  `uv run genokeeper query --json`.

## Personal data

The user's genotypes, carried variants, findings, allele sizes, copy numbers and haplogroups,
and any output that holds them (`genokeeper query`, `show`, `genotype`, `changes`, `explain`,
the dashboard), reach Claude only as the answer that needs them. Never write them into:

- auto-memory notes: the overlay's hooks commit and push them to GitHub at the end of the turn;
- git: files, commit messages, branch names, PR bodies, issues and review comments;
- the bug ledger (`docs/bugs/`): evidence from the user's data is counts and shapes only;
- claude.ai Artifacts or the Docs connector, the side-session notes folder, a session brief, a
  cloud session (a fork, `target: "cloud"`, `move_to_cloud`) or a scheduled or cloud routine;
- Telegram or any message, a web search, a fetch or any other service,

unless the user asks for that channel by name. These are rules the agent reads, not checks:
only the genome-data guard refuses anything, and only in git. The `.claude/settings.json` deny
rules refuse a read that names a file under `reports/`, `data/raw/` or
`~/genome-lab/calibration/`, never a command that reads without naming it (`grep -r`,
`find … -exec`, a Python subprocess). So search with the Grep tool, `rg` or `git grep`, which
skip ignored folders, and never run `grep -r` or `find` over the repository root.
<!-- stayfixed:harness:begin -->
## stayfixed

This repository is initialised with stayfixed. `stayfixed.toml` names its document paths and
its note store. `docs/bug-reports.md` is a generated index over `docs/bugs/` — never edit it by hand;
`stayfixed bugs index` renders it. `docs/roadmap.md` carries the design-and-plan trail between its
markers, which `stayfixed docs trail` rewrites. Specs, plans and research notes are not in this
public repository: they live in the maintainer's private notes, and new ones
are written there. Each gate runs advisory until it enforces:
`stayfixed adopt promote` enforces the gates that pass, and `stayfixed assess` says which would.

The `python` profile's rules are in `docs/stayfixed/rules/python.md`. Before the first command:

- Run commands through the package manager whose lockfile is committed (`uv run`, `poetry run`, `pdm run`, `pipenv run`); do not switch tools or add a second lockfile.
- After editing dependencies, re-lock and commit the lockfile in the same change; CI refuses a stale lock (`uv sync --locked`, `poetry check --lock`, `pdm install --check`, `pipenv install --deploy`).
- Run the test command CI runs: a `-m` on the command line replaces the one in `addopts`, and a misspelt marker in `-m "not x"` deselects nothing and passes.
<!-- stayfixed:harness:end -->
