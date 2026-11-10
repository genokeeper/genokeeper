# Genokeeper reference

The maintainers' detail behind the [README](../README.md): how the work was staged, what the
bootstrap fetches and builds, what each calibration table measured and on which samples, how a
run, the query, the explanations and the refresh behave, and each tool's terms. The README is
the user guide; read it first. A "stage" or a "slice" below is a step of the build history;
[roadmap.md](roadmap.md) has the plan behind them, and `BR-…` ids are entries in the bug ledger
([bug-reports.md](bug-reports.md)).

## Stages

1. **Findings with two confidences**, in four slices:
   - **1a** findings core: the DuckDB store, ClinVar claims, the query CLI and the agent
     skill (first on the delivered GRCh37 calls, now on our own GRCh38 ones);
   - **1b** pharmacogenomics: PharmCAT, CYP2D6, CYP2B6, the HLA drug alleles and
     PharmGKB's clinical annotations;
   - **1c** the dashboard and the generated explanations;
   - **1d** polygenic scores (pgsc_calc) with APOE, the mitochondrial haplogroup
     (Haplogrep3) and single-SNP research associations (the GWAS Catalog);
   - **1e** variants ClinVar has not classified: Ensembl VEP with REVEL, AlphaMissense,
     SpliceVault and gnomAD's pext, SpliceAI in the listed genes, graded as predictions.
2. **Re-alignment to GRCh38**: bwa, a QC gate, DeepVariant into a gVCF, chrM by mutserve,
   read-backed phasing (2a); genotype error rates calibrated on HG003 (2b).
3. **Structural variants, repeat expansions and Y haplogroups**, in slices: **3a** SMN1 copy
   number, eleven repeat-expansion loci and the Y haplogroup, calibrated on HG002;
   **3b** deletions of 50 bp or more in the ACMG SF and carrier genes (DELLY and Dysgu,
   matched by Truvari), calibrated on HG002; **3c** the foundation under them: the masked
   GRCh38, a sample identity check, per-gene coverage, sex-chromosome dosage and an
   allele-fraction guard; **3d** copy number from read depth, insertions
   of 50 bp or more, the paralog loci, GBA1 and Lp(a); **3e** thirteen more repeat loci,
   RFC1 and CFTR's poly-T/TG, with 22 loci stated as not covered;
   inversions and translocations are not covered.
4. **The monthly knowledge-base refresh** and its change feed.

## Building the store

In the order that works on a new lab (the README's [Quick start](../README.md#quick-start)
gives sizes and times):

```bash
bash scripts/bootstrap.sh --check          # list what is missing; changes nothing
bash scripts/bootstrap.sh --with-gnomad --with-vep --with-masked-reference --with-sv --with-typing
uv run genokeeper run                      # aligns, records your bases, stops: not calibrated
bash scripts/bootstrap.sh --with-calibration --with-hg002 --with-hg004  # downsampled to your bases
uv run python -m genokeeper.calibrate.hg002 align pcrfree  # also pcr, and hg004 (refused unless female)
uv run python -m genokeeper.calibrate.hg002 accept hg004   # accepted at SNP F1 0.99 or more
uv run genokeeper calibrate [--only NAME]  # every stale calibration table, in order
uv run genokeeper run                      # FASTQ -> GRCh38 calls -> reports/genomelab.duckdb
```

A lab whose libraries were aligned before slice 3c moves them onto the masked reference with
`realign` instead of `align`, then `clean` (slice 3c, below).

The bootstrap fills `~/genome-lab/`: the ClinVar GRCh38 weekly VCF, ClinGen's gene-disease
validity export, the CGD gene table, NCBI's GRCh38 no-alt FASTA (0.87 GB download, about
3.2 GB unpacked, md5-checked; CRAMs written before slice 3c decode against it, and the FAF
builder normalises with it), the GIAB v3.6 GRCh38
stratifications and the region BEDs derived from them, the DeepVariant T7 checkpoint, bwa,
mosdepth and rtg-tools at pinned versions, PharmCAT's positions, CPIC's gene-drug pairs and
MT-RNR1 alleles (from its JSON API), ClinPGx's
clinical annotations, the GWAS Catalog's associations (with their ontology annotations, filed
under the release's build date), MedlinePlus Genetics' gene and condition summaries (the
explanations' grounding), IPD-IMGT/HLA's `hla.dat` with T1K's HLA index built from it, the
container images (pulled by digest; Aldy's,
BCyrius's and Haplogrep3's, which no registry publishes at their pinned versions, built from
`docker/` with every input pinned, Haplogrep3's with its PhyloTree tree baked in), and
slice 1d's polygenic scores under `pgs/`: Nextflow 25.10.0 with its nf-schema plugin (it runs
offline, on Homebrew's openjdk@21), pgsc_calc at the commit its tag v2.3.0 points to, the
chosen scores' harmonised GRCh38 scoring files and the research scores' (a separate step,
so that fetching one never re-fetches a chosen file; all md5-checked against the PGS
Catalog's; slice 3g's 34 new files are 677,577,718 bytes, fetched once they were approved), the
HGDP+1kGP reference panel (16 GB, md5-checked, moved from `prefetch/1d/pgsc_calc/` when a
copy with that md5 is there) and the panel's ancestry positions, filtered from it by plink2
in pgsc_calc's own image. The
flags combine in any order. `--with-gnomad` also builds the gnomAD v4.1 FAF table
(`uv run genokeeper build-faf`): it streams the 24 joint sites VCFs once (0.88 TB, chr22 alone 15.6 GB), each checked against the bucket's md5, and
keeps a few GB of Parquet under `gnomad/4.1/`. A rerun streams only the chromosomes the table
lacks or whose bucket md5 changed; without the flag the bootstrap only says the table is
missing, and `--check` lists the missing chromosomes. `--with-calibration` fetches stage
2b's inputs: GIAB's HG003 v4.2.1 GRCh38 benchmark, HG003's DNBSEQ-T7 FASTQ from CNGBdb
(147 GB, md5-checked, resumed after a cut transfer, not fetched again once downsampled) and
an rtg SDF of the GRCh38 FASTA, under `calibration/` and beside the reference.
`--with-vep` fetches slice 1e's files into `vep/116/`, each sha256-pinned and moved from
`prefetch/1e/` when a copy with the pin is there: Ensembl's indexed VEP 116 GRCh38 cache
(27.6 GB, unpacked), REVEL v1.3 (sorted on GRCh38 and indexed for VEP's plugin),
AlphaMissense, SpliceVault's GRCh38 data, gnomAD v4.1's base-level pext (as an indexed BED
with each gene's largest value) and UCSC's NMD Escape track for MANE (as BED, by
bigBedToBed in its pinned image); without the flag it says they are missing. ClinGen's
dosage sensitivity list comes with every bootstrap, and SpliceAI's image is built from
`docker/` with its TensorFlow, models and dependencies pinned.
`--with-publications` also fetches the two publication sources (about 2.2 GB, no grading
input): LitVar2's whole variants file into `litvar/<Last-Modified date>/` and ClinVar's
`var_citations.txt` into `clinvar-citations/<date>/`. `--refresh STEP` (STEP: `clinvar clingen
clingen_dosage cpic pharmgkb gwas medlineplus litvar clinvar_citations`) fetches that one
source's newest release even when the manifest says it is done, files it under its dated folder,
prints `refreshed: STEP <folder>` and exits non-zero on any failure; `genokeeper refresh` calls
it once per source. Every one of these requests is for a whole file at a fixed URL.

Stage 3a's tools run from images the bootstrap builds from `docker/` (by hand: `uv run python -m
genokeeper.containers build <name>`), each with its inputs pinned and read with `--network none`: `expansionhunter`
(ExpansionHunter 5.0.0), `reviewer` (REViewer 0.2.7, its own image because its htslib pin
differs from ExpansionHunter's), `smncaller` (SMNCopyNumberCaller 1.1.2), `parascopy`
(Parascopy 1.19.0) and `yleaf` (Yleaf 4.1.6). Every bootstrap also fetches Parascopy's
precomputed data (52 MB, Zenodo record 15019940, md5-checked). `--with-hg002` fetches GIAB's
HG002 v5.0q GRCh38 truth and the two public HG002 libraries (ENA, PRJNA1465120): each FASTQ is
downloaded whole and resumably over parallel range requests (one raw file, about 150 GB, on
disk at a time), md5-checked, downsampled to the bases of your own run (which a first
`genokeeper run` records; before one, the libraries wait and the bootstrap goes on) and deleted, so about 170 GB of
downsampled FASTQ stay under `calibration/`. `uv run genokeeper calibrate` (or `--only
stage3`) runs both libraries through the run's own align, QC gate and the three tools and
writes the stage-3 calibration table under `calibration/`; a run refuses a table made with
other image digests, catalog or reference. A watch item (a borderline or intermediate repeat,
an SMN1 silent-carrier marker) or a non-finding of SMN, a repeat locus or a paralog locus is
neither in the feed nor in its counts: the dashboard lists them under the feed, each kind
apart, with its call and its stored reason where it has one (CFTR's poly-T/TG aside, and a
normal flagged as an exclusion in its own block).

Slice 3b's callers run from images built the same way: `delly` (DELLY 2.7.0, one container
per primary chromosome, as it parallelises over samples only), `dysgu` (Dysgu 1.9.0) and
`truvari` (Truvari 5.5.0), each with `--network none`. The bootstrap's `sv` step moves
gnomAD v4.1's SV sites VCF (1.74 GB) and its index into `sv/` from `prefetch/3b/` when the
copy there has the pinned sha256, or downloads them with `--with-sv` (without it the step
says what is missing); `--with-hg002` also fetches GIAB's HG002 v5.0q structural-variant
truth (50 MB, sha256-pinned). `uv run genokeeper calibrate` (or `--only deletions`) runs
both HG002 libraries through the run's deletion steps and writes the deletion table under
`calibration/deletions/`: each deletion class's chance of being false over the listed
genes' exons. A run without the `sv` files, or without a current deletion table, leaves the
deletion screen unassessed with the reason and goes on; the run's tools name DELLY, Dysgu
and Truvari only when the screen ran.

Slice 3d screens the same exons, with DPYD and CHEK2 added, for what the deletion screen
leaves. Losses and gains of read depth of 1 kb or more come from DELLY's `delly cnv` (the
3b image, one whole-genome run) and CNVpytor 1.3.2 (`cnvpytor`). Insertions of 50 bp or
more come from Dysgu and INSurVeyor 1.1.3 (`insurveyor`), matched by Truvari, with mobile
elements labelled by Scramble 1.0.2 (`scramble`) on the target exons only. Each goes through
3b's gene gate and wording: at validity limited, needing confirmation. A depth event of
500 kb or more is shown as "screening result, needs a chromosomal microarray" only where it
overlaps a region or gene ClinGen scores 3 for haploinsufficiency or triplosensitivity, and
ends in a segmental duplication only where it matches a ClinGen recurrent CNV region. Such an
event is shown, not calibrated: no truth set holds one. The bootstrap's `clingen_regions` step fetches ClinGen's region curations and
recurrent CNV BED (CC0), and rerunning `scripts/bootstrap.sh` refreshes them; `genokeeper
refresh` does not. The GJB6-D13S1830 deletion is a named DFNB1 carrier locus. STRC
(Parascopy, with CATSPER2's unique depth) is a carrier locus. RHD, AMY1, C4 and CYP2A6
(Parascopy, CYP2A6 checked against Aldy) are information, and RHD is a copy number, never a
blood type. GBA1 comes from Gauchian 1.0.2 (`gauchian`), its copies compared with Parascopy's
over GBA1's GBAP1-homologous exons 9-11. Lp(a)'s KIV-2 count comes from KILDA 1.2.3 (`kilda`),
as information beside the blood test. ClinGen's recessive genes at Definitive or Strong with
no dominant curation (Tier 4, 1,309 genes) join Tier 3 in the carrier gate. `uv run genokeeper
calibrate` builds the `cnv`, `insertions` and `paralogs` tables on HG002's two libraries
and the HG002/HG003/HG004 trio, and measures the `deletions` table again, because its exon set
and class features change. Each depth or insertion class's chance of being false is measured
as 3b's, gains also by Mendelian consistency on the trio. A paralog locus is withheld where any
of the four libraries gave no call, where HG002's two libraries disagree, where a library
contradicts its sample's truth (HG002's assembly-derived copy numbers, HG003's and HG004's for
STRC), or where HG002 shows more copies than its parents together; GBA1 also where HG002 is not
wild type, and Lp(a) where KILDA's count differs from HG002's by more than
`[lpa] max_kiv2_error` (a share of that count). A screen without a current table is
unassessed with its reason.

The 3d tables as measured on 2026-11-05 and 06, on the masked reference. A withheld item means
"not graded": the statement stays unassessed with the table's reason, and `claim.withheld` keeps
the run's own value at low reliability for `genokeeper query --json` and the `genome-ask` skill,
never as a result.

- **Withheld:** AMY1 ("not covered: no call in calibration": HG004 gave no usable call, while
  both of HG002's libraries gave the truth's 6). STRC, RHD, C4, CYP2A6, GBA1 and Lp(a) are graded.
- **Not measured:** a mobile element only Scramble called. Neither HG002 library had a rare one
  inside the benchmark around the target exons (0 of the 8 truth insertions there found), so the
  `mei` band has no rate and such an event reads "not calibrated: no mei calls measured", its
  value kept the same way. An insertion Dysgu or INSurVeyor called keeps its own class, and
  Scramble only labels it.
- **Large events:** none of the three genomes gave a large call, so no large band is withheld.
  A large event in a run has no rate: it is shown, not calibrated.
- **CNVpytor stays:** 356 of the 405 calls both depth callers made were false, against 1,583 of
  DELLY's 1,755.

| Band | Truth events found (HG002 PCR-free, PCR) | Chance of being false, across the classes measured |
|---|---|---|
| Depth loss, 1-5 kb | 20 and 31 of 42 spike-ins | 0.971 to 0.999 |
| Depth loss, 5-500 kb | 17 and 16 of 18 spike-ins | 0.832 to 0.965 |
| Depth gain, 1-5 kb | 2 and 7 of 40 spike-ins | 0.963 to 0.983 |
| Depth gain, 5-500 kb | 11 and 16 of 20 spike-ins | 0.915 to 0.992 |
| Insertion, 50-999 bp | 1,643 and 1,494 of 12,243 | 0.998 to 1.000 |
| Insertion, 1 kb or more | 96 and 93 of 2,113 | 0.998 to 0.999 |
| Deletion, 50-999 bp | 1,328 and 1,195 of 7,260 | 0.851 to 0.996 |
| Deletion, 1 kb or more | 86 and 86 of 197 | 0.136 to 0.997 |

Every class's point estimate stays above `[deletions] likely_below` (0.1), so every event with
a rate reads "possible", none "likely". No call is weak evidence of no event: the found column is the screen's reach. The
insertion rates are if anything high (BR-037). Eight classes have an interval wider than 0.5
and are kept as measured: the depth gains of 1-5 kb called by both callers with an end masked,
and by DELLY alone with a junction, masked or not; the depth gains of 5-500 kb with a
junction, by both callers with an end masked and by DELLY alone with none; and the deletions of
1 kb or more with an end masked by both callers with depth support, and by Dysgu alone, masked
or not.

Slice 3e's repeat loci run from the one `expansionhunter` image, now built from source:
`uv run python -m genokeeper.containers build expansionhunter` (or the bootstrap), then `uv
run genokeeper calibrate` (or `--only stage3`) re-measures the stage-3 table with the 24
loci, RFC1's motif analysis, CFTR's poly-T/TG and the trio checks on HG002, HG003 and HG004;
a run without a current table leaves every graded repeat statement unassessed, naming
`genokeeper calibrate`, and still writes the 22 not-covered statements.

The table of 2026-11-06 keeps 15 of the 24 loci (HTT, AR, DMPK, C9orf72, ATXN1, ATXN3, ATXN7,
ATN1, FXN, TCF4, PABPN1, CNBP, ATXN10, JPH3, FGF14) and CFTR's poly-T/TG
(`expansionhunter:CFTR-polyT-TG`), and withholds the 12 items below. The table has no HG002 truth
for FXN, TCF4, CNBP or FGF14, so those four rest on the libraries' agreement and the trio. A
withheld item is not graded, which is not dropped: it stays unassessed with the table's reason,
and the value the tool measured stays visible to the skill at the lowest reliability (see
[Asking about the genome](#asking-about-the-genome)). Each reason below starts "not covered:" in
the store. PCR-free and PCR are HG002's two libraries.

| Withheld | The table's reason | What was measured |
|---|---|---|
| ATXN2 | false call in calibration | Truth 22/30. PCR sized it exactly; PCR-free called 22/27 (interval 26–64), 3 off and not Mendelian. HG002's 30 is above the normal maximum of 28, so even an exact call is not a normal. |
| CACNA1A | false call in calibration | Truth 11/11. PCR-free's interval reached 19, past the normal maximum of 18; PCR called 10/11. |
| FMR1 | false call in calibration | Truth 31. PCR-free was exact; PCR called 24 with interval 19–60 and no spanning read. |
| AFF2 | false call in calibration | Truth 19. PCR-free was exact, but no allele of HG004's is within 2 repeats of it; PCR called 16 with interval 16–79 and no spanning read. |
| FMR1's two-X branch | the X rule did not hold on the trio in calibration | HG004 must be a two-X normal carrying HG002's allele. It read 30/31 on 1 spanning read, unassessed. |
| AFF2's two-X branch | the X rule did not hold on the trio in calibration | HG004 read 26/26 on no spanning read, unassessed, and neither allele is near HG002's 19. |
| NOP56 | HG002's two libraries disagree in calibration | Both 4/7. PCR-free has 2 spanning reads behind an allele, under the floor of 3, so it is unassessed; PCR is a normal on 10. |
| TBP | HG002's two libraries disagree in calibration | Both 37/36. PCR-free is a normal on 5 spanning reads; PCR has 2, unassessed. |
| PPP2R2B | the rule never gave a normal on a healthy calibration sample | All four samples read 10/10. Two alleles of one length are unassessed under the exclusion rule. |
| ZFHX3 | the rule never gave a normal on a healthy calibration sample | All four samples read 20/20, unassessed the same way. |
| THAP11 | the rule never gave a normal on a healthy calibration sample | All four samples read 29/29, unassessed the same way. |
| RFC1 | the RFC1 motif call was not a clear normal on a healthy calibration sample | All four samples read normal, "no expanded allele", but HG002 carries a benign expanded AAAAG allele (GIAB's v5.0q truth) that both libraries had to see. They missed it, so a normal here has not been shown to see an expanded allele (BR-038). |

While FMR1 and AFF2 are withheld as genes, a two-X sample's statement carries the gene's reason
rather than its branch's.

`uv run genokeeper calibrate [--only NAME] [--report] [--threads N]` is the one calibration
command: it rebuilds every stale table (one whose recorded versions differ from this
pipeline's) in dependency order, each swapped in whole, and says "current" for the rest; a
table whose inputs are missing is skipped with one line naming them and the command that
supplies them (exit 2). `--only NAME` builds that table alone (`small_variants`,
`hg002_cells`, `stage3`, `deletions`, `cnv`, `insertions`, `paralogs`, `foundation`,
`slice3f`) and
refuses while one it needs is stale. The `small_variants` table downsamples HG003's FASTQ by bases to your run's (the total
a `genokeeper run` recorded in its QC; it refuses before one has) into
`calibration/hg003/downsampled/`, runs it through the run's
own align, QC gate, sex (refusing unless male) and call steps, scores its gVCF against GIAB's truth with rtg, and writes the
calibration table (`calibration/small_variants/`): genotype error per GQ bin, row type, SNV or indel
and easy or difficult region, with the versions of the pipeline it measured (bwa, samtools,
the DeepVariant image, the checkpoint's and reference's md5, the DeepVariant flags).
`--report` prints the table, HG003's SNP and INDEL F1, and its held-out checks (odd against
even chromosomes, chr1–19 against chr20–22), each beside the limit that stops stage 2b.

Slice 3c aligns everything to GIAB's GRCh38 with the GRC false duplications masked (the
lab's `reference/GRCh38_no_alt_maskedGRC.fa`): `--with-masked-reference` moves it from
`prefetch/3c/` at its sha256 pin or downloads it, derives the masked intervals by comparing
it with NCBI's FASTA (refused unless they equal the committed copy), bwa-indexes it and
builds its rtg SDF. GIAB's v4.2.1 truths were made on the unmasked reference and can lack a
real variant on the copy of a false duplication that stays, so HG003's table and HG004's F1
leave out those correct copies as well as the mask. At run time HG002's cells rate a site on
such a copy where v5.0q's benchmark holds it, and any other site there is uncalibrated. The
libraries aligned before it are re-aligned from their own CRAMs:
`python -m genokeeper.calibrate.hg002 realign pcrfree` (and `pcr`, `hg003`) turns the old
CRAM back into a FASTQ pair, refused unless its reads and bases equal the CRAM's, and
aligns it again; `clean` then deletes what the new record replaced. On a new lab, `align
pcrfree` and `align pcr` align HG002's downsampled pairs, and `genokeeper calibrate` aligns
HG003 itself. `--with-hg004` fetches
HG004's v4.2.1 truth and downloads its PCR-free library (CNGB `CNR1441269`, from the lane
of 2b's HG003; ENA's trio runs are mislabelled), downsampled to your run's bases; `align hg004` aligns and calls it, refused unless its
sex is female, and records it not yet accepted; `accept hg004` then scores it against
HG004's truth and accepts it only when its SNP F1 reaches 0.99.
The `foundation` table then measures, on HG002's two libraries, HG003 and HG004, the X and
Y depth ratios and the chrX allele balance a normal karyotype gives, and HG002's true
heterozygotes' allele fractions; `--report` proposes the `[dosage]` and `[af_guard]`
bounds, which go into `config/grading.toml` by hand, in a commit citing the report.

The `foundation` table as measured on 2026-11-05, on the masked reference, over HG002 PCR-free,
HG002 PCR, HG003 and HG004. The three male libraries' X depth ratios run from 0.4908 to 0.5090
and their Y ratios from 0.5121 to 0.5306; HG004, female, reads 0.9847 and 0.0276, and its SNP F1
against its v4.2.1 truth is 0.9962. The chrX heterozygous share is 0.0371 on HG002 PCR-free
(69,855 sites), 0.0648 on HG003 (71,552) and 0.5912 on HG004 (100,474); HG002 PCR is aligned and
not called, so it gives the depth ratios alone. Those measurements set `[dosage]`'s `x_het_min`
at 0.33 and `y_absent_max` at 0.06. `[af_guard]`'s `min_het_vaf` of 0.25 comes from HG002
PCR-free alone: 2,138,950 true heterozygous SNVs, 1.399% of them below 0.30.

`genokeeper run` reads the paired FASTQ in `data/raw/` (`--fastq1`/`--fastq2` name them),
checks the reads come from the platform the T7 checkpoint was trained on, aligns them with
bwa, refuses a sample that fails the QC gate (`config/grading.toml` [qc]) before calling,
calls it with DeepVariant into a gVCF and chrM with mutserve, and grades every claim of the
ClinVar GRCh38 release on our own call at its allele; candidate compound heterozygotes are
phased from the reads. Then the pharmacogenes: PharmCAT 3.4.0 on our calls at each of its
positions (a call graded low or unknown goes in as a no-call), given the outside callers'
calls (CYP2D6 where BCyrius and Aldy agree on the activity score, or where StellarPGx's score equals one of theirs when they do not, CYP2B6 from Aldy, HLA-A and
HLA-B where OptiType and T1K agree on the drug alleles, MT-RNR1 from mutserve), each gene
stored as a gene call under drugs; and PharmGKB's level 1B, 2A and 2B clinical annotations
in genes PharmCAT's report gives no drug (F2 and F5 among them), located at their rsID's ClinVar allele and graded on our call
there. Then the polygenic scores, 45 scores (11 from slice 1d, 34 from slice 3g) beside 24
further conditions named not covered with the reason: pgsc_calc v2.3.0 (Nextflow, offline;
every process in its pinned image with `--network none`, checked task by task after the run,
which refuses results it cannot verify) on a VCF of our calls at every chosen score's and the
panel's ancestry positions (an unknown or low call goes in as a no-call), each score a finding with
its percentile against the population you are closest to (one below 0.75 overlap is
unassessed, what pgsc_calc computed kept under `claim.withheld`), and APOE's diplotype from our
calls at rs429358 and rs7412. Docker needs at least 32 GB of memory for it (Docker Desktop,
Settings, Resources). 27 more scores outside the chosen set (`[[research_score]]`
in `src/genokeeper/resources/pgs/scores.toml`, each with `why`) run in a second pgsc_calc run of
their own, so the chosen scores' cache stays as it was, and are stored as research: validity
limited, topic research, never actionable, returned only to `--gene`, `--marker` or `--text`
and never on the dashboard's feed. A research score says little about one person, and one
below the 0.75 overlap is still shown with its overlap: nothing is hidden, its reliability is
stated. Their target writes each allele on a row of its own: a person who carries another ALT
of a position reads as a no-call on this one's row. Then the mitochondrial
haplogroup: Haplogrep3 3.3.2 (offline, in its image) on mutserve's homoplasmic chrM calls,
on PhyloTree Build 17's rCRS tree (`phylotree-rcrs@17.3`); the chrM positions too few reads
cover are named beside it, and with more than 1% of chrM uncovered it is unassessed. Then
the research associations: the GWAS Catalog's single-SNP associations at p ≤ 5e-8 (one per
SNP, risk allele and trait set, the strongest study's), each located on the HGDP+1kGP
panel's variants (their alleles give the SNP's REF and ALTs; a risk allele that is none of
them or could be read on either strand is skipped, and so is a trait that speaks of a
dose) and graded on our call there: carrying the risk allele, also
when it is the reference, is a finding, at validity limited, with the study's odds ratio
or beta and never an absolute risk. Then the variants ClinVar has not classified: VEP 116
(offline, in Ensembl's image) on every PASS allele of ours that no ClinVar record
classifies germline, whatever its shape, SpliceAI on those in the listed genes; a predicted
loss of function in a gene ClinGen rates haploinsufficient, or in a carrier Tier 3 gene,
that passes the LoF filter's six gates (homopolymer and low mappability, MANE, NMD on the real stop with ClinVar's downstream
truncations as rescue, no start-loss, pext, SpliceVault), a missense at or above its
calibrated REVEL or AlphaMissense score (REVEL read from another transcript of the same protein
change where REVEL v1.3 does not name the transcript and its rows agree), or a splice change at SpliceAI 0.2, is a finding
at validity limited, "predicted, not classified", with a note to discuss it with a
geneticist (`config/grading.toml` [predicted]); outside the listed genes, or failing a
gate, it is research. A predicted finding never outranks a ClinVar-classified one in the
feed. Rare, probably damaging variants in the genes of Genomics England PanelApp's
Intellectual disability panel (panel 285, version 11.93, 2,867 genes, one pinned file from
the bootstrap) are stated from the same VEP records when gnomAD's filtering allele frequency is under
1%, and under the gene's own limit where it has one, as research: validity limited,
topic research, never actionable, returned only to `--gene`, `--marker` or `--text` and
never on the dashboard's feed, each `predicted, not classified`, naming the panel, the
gene's confidence level and its inheritance, and repeating no predicted finding. SpliceAI
scores only the listed genes, so a splice change in another panel gene is not scored by
it; a panel file that is missing or damaged leaves one unassessed statement and the
rest of the run as it is. It is a computer's prediction on a broad panel, which most
people carry something on, and never a diagnosis, of intellectual disability, ADHD or
autism alike. The store keeps no dosing text; PharmCAT's HTML report, whose path the run prints,
stays under `reports/work/pharmcat/` for you. `--sex male|female` overrides the sex
inferred from the X copy count (chrX non-PAR depth over the autosomes: one X male, two female);
`--threads` defaults to this machine's CPUs. Each slow step is cached under `reports/work/`
and reused while its inputs are unchanged. Each call's `p_wrong` is the error rate measured
for its cell, with its interval, on GIAB reference samples: HG003 inside its benchmark regions,
and elsewhere HG002's cells, which may have been in DeepVariant's training, inside v5.0q's
benchmark and at the locus's ploidy, once they are current. They hold a male's chrX and chrY
outside the PARs and the correct copies of the false duplications. Any other call, a non-PAR
chrX call with two X chromosomes among them, is marked "outside calibrated regions". The run refuses
to start without a complete FAF table (naming `--with-gnomad`), with a calibration table
measured under other pipeline versions (naming `genokeeper calibrate`), or without any other
input the bootstrap supplies (naming `scripts/bootstrap.sh`); with no calibration table it
aligns and measures QC, recording the bases the calibration needs, then refuses.

Slice 3c's checks run on every run. After the call, the run compares its genotypes at the
GWAS Catalog SNPs with the delivered 23andMe-format file in `data/raw/` (the one `*.txt`,
not `*.vc.txt`, with its column line; `--genotypes PATH` names another, and two candidates
refuse): below 0.98 concordance over at least 1,000 SNPs it refuses before grading; with
fewer SNPs, or no file, the identity is recorded "not checked" with the reason. This catches
a FASTQ mix-up, not a sample swap at the lab, which derived that file from the same reads.
Each heterozygous call below `[af_guard] min_het_vaf` is warned "allele fraction low", and
in the nine clonal-haematopoiesis genes its reliability is capped at low ("may be clonal,
not inherited; confirm in another tissue"). Each gene of the ACMG SF and carrier lists, the
stage-3 genes and CBS, KCNE1, CRYAA and U2AF1 gets a coverage row: the share of its MANE CDS
covered at MAPQ 20 by 10 reads per expected copy, beside the shares in difficult strata and
outside the calibrated benchmark; a gene below `[coverage] min_callable` is "not reliably
assessed" in `query --gene` and `genotype --gene`, unless a stage-3 source states it. The
sex-chromosome dosage (X and Y depth ratios, the chrX allele balance) is one statement,
`dosage:sex-chromosomes`: every non-normal pattern a "screening result, needs a karyotype",
mosaic loss of Y research information, and unassessed while the `foundation` calibration table
is missing ("not calibrated") or stale; FMR1 and AFF2, which read that dosage, are then
unassessed too.
The run prints the identity verdict and how many genes fall below the callable share, and
records the reference in its tools: `genokeeper changes` names "reference changed" (the
Changes tab, "reference genome", in Russian «референсный геном») as the cause of every change
between two runs on different references, so a call the mask moved is not read as new knowledge. It prints its accounting, one line per source (ClinVar,
PharmGKB, PharmCAT, pgsc_calc, Haplogrep3, the GWAS Catalog, VEP, stage 3a's SMN, ExpansionHunter and Yleaf, slice 3b's deletions, slice 3d's `cnv`, `insertions`, `paralogs` and `lpa`, and slice 3f's typing sources): each input row is a finding, a watch item, a non-finding, an unassessed
claim (a finding were the genotype or gene call known) or dropped with a reason. A store another schema version wrote is refused untouched (there are no
migrations): move it aside and run again. Then read the store with `genokeeper query` /
`show` (below) or ask Claude Code.

Slice 3f types what the variant track does not. `--with-typing` fetches its HHV-6 genomes and
hmmix's files (about 3.5 GB, each sha256-pinned) and builds the panel files admixture and
runs of homozygosity read; `genokeeper calibrate --only slice3f` measures its table on the
GIAB trio, and without a current table each 3f statement is unassessed "not calibrated".
OptiType and T1K type HLA: a class I allele both read is a finding (HLA-B\*27, the PharmGKB
2A/2B drug alleles, C\*06:02 and B\*51 as research associations), and the coeliac screen
needs T1K's DQ types and six DQ tag SNPs to agree; its rule-out is a non-finding flagged as
an exclusion, which the dashboard lists under the feed in its "Excluded by screening" block
with the test that would confirm it. KIR is not covered. mutserve's chrM heteroplasmies are
graded against a floor (`[mito] heteroplasmy_floor`, 10%, which the calibration can raise
but never lower), MITOMAP's Confirmed variants join ClinVar's, and haplocheck's
contamination check withholds the heteroplasmies (so does haplocheck failing, its image
missing on a lab bootstrapped before 3f among the causes). The unmapped reads are aligned to HHV-6A
and HHV-6B for inherited integrated HHV-6; bcftools roh reports runs of homozygosity;
RBCeq2 2.4.4 types every blood group system but RhD and RHCE, MNS and the platelet antigens
(its GRCh38 table holds GRCh37 positions for them), and a system with a no-call at any of
its positions is unassessed; pgsc_calc's similarity, fastmixture's supervised admixture
and hmmix's Neanderthal share describe ancestry; and eleven traits (lactase persistence,
ADH1B, ACTN3, TAS2R38, CYP1A2, and slice 3g's six single-gene traits: FUT2, ABCC11, MC1R,
HERC2, AMPD1 and CCR5Δ32) are read from our genotypes. The ancestry statements, the
blood group phenotypes, the HHV-6 band and the traits are information, not calibrated.
PharmCAT gains read-backed phase sets and UGT1A1\*28 when the table trusts them, and the G6PD
haplotypes it does not call are named; PharmGKB's level 3 and the HLA disease associations
join the GWAS Catalog as research associations, returned only to a query that names them; GenCC joins ClinGen for
gene-disease validity; and a polygenic score is shown as its decile. Nothing in 3f is
actionable.

**What the `slice3f` table measured** (`genokeeper calibrate --report` prints it). It was built
on 2026-11-06 from the GIAB trio's three called libraries (HG002 PCR-free, HG003, HG004).
Withheld means not graded, not dropped: the statement stays unassessed, and the claim keeps the
tool's value under `withheld`, at the lowest reliability, with the table's reason (a chrM
heteroplasmy may keep its level in the claim's `level` or its warnings instead). Two things are
withheld:

- **Eight blood group systems: ABO, C4A, C4B, CTL2, FUT3, MAM, PIGG and XG.** Each is measurable
  on no library: on each of the three, a position the system reads is a no-call in our calls or
  in the truth's, so RBCeq2's genotype on one could not be set against its genotype on the
  other. That is a missing measurement, not a wrong genotype. The other 37 of the 45 packaged
  systems were compared on at least one library, and no genotype differed.
- **chrM heteroplasmy.** The highest call below `[mito] min_homoplasmy` (0.9) in the unmixed
  libraries was 0.58. One at 0.20 or above leaves the table with no floor, so a MITOMAP or
  ClinVar statement on a chrM call below homoplasmy stays unassessed, with the reason "false
  heteroplasmy at 0.58 in calibration". The gate counts every call below homoplasmy in the
  unmixed libraries as false, and no truth set for chrM heteroplasmy was compared, so a real
  heteroplasmy in one of these cell lines would withhold it too: "false" is the gate's
  assumption, not a finding. Homoplasmic chrM calls and the haplogroup are not touched.

PharmCAT's read-backed phase sets are not trusted: one pair of sites was compared with HG002's
phased truth and agreed, and the table needs 20 (BR-036). PharmCAT therefore runs without them,
as before 3f; the phasing of compound-heterozygote candidates and of APOE is not gated by this
table. The other measures:

- UGT1A1\*28 is trusted: its TA repeat matched the truth on all three libraries, each inside
  its benchmark, and HG003 carries TA7.
- HLA: no gene with a published type is withheld. OptiType (class I) and T1K gave HG002's type
  at all six that have one (HLA-A, -B, -C, -DRB1, -DQA1, -DQB1), and each parent carries the
  allele HG002 got from it; T1K's other genes have no truth here. No coeliac DQ tag SNP
  differed from a truth; all six are graded on HG003 and HG004, four on HG002.
- haplocheck, which is no gate here because the table has no floor, flagged all four mixtures
  of HG003's chrM reads into HG002's (2% to 20%) and none of the three unmixed libraries.
- HHV-6: background ratios of 0, 0 and 4.0×10⁻⁵, a non-finding ceiling of 1.2×10⁻⁴ against the
  band's 0.5. The trio offers no positive, so this measures the background only.
- Runs of homozygosity: no run of 5 Mb or more is false or missed on any library, and the 12 Mb
  isodisomy spiked into chr2 was found on all three.
- The Neanderthal share: 1.06% to 1.15%, inside ±30% of the 1.18% published for NA20810.
- Ancestry: HG002's fastmixture proportions are within 0.45 percentage points of its parents'
  mean in every group. That is a consistency check with no gate, so ancestry stays information.

**Aldy's terms.** CYP2D6 (with BCyrius) and CYP2B6 are called by Aldy 4.8.3, whose licence
from the Indiana University Research and Technology Corporation (IURTC) permits use solely
for non-commercial research, excludes clinical trials and diagnostic use involving human
subjects without IURTC's written consent, and forbids providing data processing to a third
party or distributing the software (its `LICENSE.rst`, shown by `aldy license`). Each user
accepts that licence for their own non-commercial research on their own genome; running the
lab on someone else's genome for them is processing for a third party. Its pharmacogenomic
results are not a diagnosis.
ClinPGx's (PharmGKB's) annotations download (`summaryAnnotations.zip`, the one ClinPGx
serves now) is licensed CC BY-NC-SA 4.0: personal use, as here, is within it, and it asks for
attribution to ClinPGx. The frozen `clinicalAnnotations.zip` an earlier bootstrap fetched
(CC BY-SA 4.0) is still read. The downloads stay in the lab and the store, never
redistributed. One table derived from that frozen release is in the repository,
`src/genokeeper/resources/hla/alleles.tsv`, under CC BY-SA 4.0 (NOTICE).

**Polygenic score terms.** pgsc_calc is Apache-2.0. Each scoring file comes under the PGS
Catalog's terms its score records (`src/genokeeper/resources/pgs/scores.toml`, `licence`: cite
the score's publication and follow its authors' restrictions under EBI's terms of use; the
chosen ones add none). The HGDP+1kGP panel is the PGS Catalog's build of gnomAD v3.1's HGDP
and 1000 Genomes callset, whose primary data gnomAD releases under CC0 1.0. All stay in the
lab, never redistributed.

**Haplogroup terms.** Haplogrep3 is MIT-licensed, and so is genepi's packaging of PhyloTree
Build 17 (van Oven, 2016) it classifies on.

**Prediction terms.** SpliceAI's code is GPL-3.0 and its models and scores CC BY-NC 4.0;
AlphaMissense's predictions are CC BY-NC-SA 4.0 (the copyright line in its file): both are
non-commercial only, and personal use on your own genome is within that.
REVEL's scores are free for non-commercial use; VEP is Apache-2.0. The score files stay in
the lab, never redistributed, and a predicted finding is not a diagnosis.

**Stage 3a terms.** ExpansionHunter 5.0.0 (Apache-2.0 at tag v5.0.0), SMNCopyNumberCaller
1.1.2 (Apache-2.0) and Parascopy 1.19.0 (MIT) are permissive. REViewer 0.2.7 is GPL-3 at its
tag, and Yleaf 4.1.6's LICENSE.txt says GPL-3 while its setup.py says MIT: both are run here
for you, never redistributed. The repeat catalog is cut from Illumina's ExpansionHunter v5.0.0 catalog, its thresholds are STRchive v2.28.0's, and
Parascopy's data is Zenodo record 15019940's. Illumina has since moved ExpansionHunter,
REViewer, SMNCopyNumberCaller, Cyrius and SpliceAI to the PolyForm Strict License 1.0.0 on
their main branches, so a newer version of any of them is read for its licence before use.
A stage-3 result is a screening result that needs confirmation, never a diagnosis; FMR1 is assessed only by slice 3e's rule (below).

**Stage 3e terms.** ExpansionHunter runs at commit `9e93322f4a` (Illumina, 2022-02-02),
whose LICENSE.txt is Apache-2.0; the licence became PolyForm Strict later, so no newer commit
is used without reading it first. It is built here from source, offline and by checksum:
htslib 1.22.1 (not the 1.10.2 the commit names: the lab's CRAMs are version 3.1, which
htslib reads only from 1.12 on), Boost 1.73.0 from archives.boost.io (the URL the commit
names now redirects), spdlog 1.6.1, abseil 20210324.2 and googletest 1.10.0. The image stays
in the lab, never redistributed. The catalog is Illumina's v5.0.0 GRCh38 one plus four
entries cut to whole repeat units from STRchive v2.28.0's ExpansionHunter asset (PABPN1,
ZFHX3, THAP11, FGF14) and CFTR's poly-T/TG, read off GRCh38 and ClinVar; the thresholds are
STRchive v2.28.0's. Out by design: TRGT (long reads only, and a PacBio-only licence),
ScatTR (validated on simulations alone) and genome-wide repeat screening (thousands of
alleles nobody can interpret without a cohort on the same platform). Every non-normal result
is a screening result that needs confirmation, never a diagnosis; a normal at FMR1, AFF2 or
an exclusion-only locus says nothing about anyone but the sample.

**Stage 3b terms.** DELLY 2.7.0 (BSD-3-Clause at tag v2.7.0), Dysgu 1.9.0 (MIT at tag
v1.9.0) and Truvari 5.5.0 (MIT at tag v5.5.0) are permissive. The exons are cut from NCBI's
MANE v1.5 (`src/genokeeper/resources/deletions/SOURCE.md`), the population frequencies are
gnomAD v4.1's SV sites, and the truth is GIAB's HG002 v5.0q. A shown deletion counts as a
loss of function only where 1e's predicted ones do (a carrier Tier 3 or ClinGen
haploinsufficient gene); it is at validity limited and needs confirmation by MLPA, a
microarray, qPCR or long reads, never a diagnosis. Duplications are slice 3d's; SLC37A4 is
not covered (MANE v1.5 places its only transcript on a fix contig, off the primary assembly the
callers run on), and the screen can miss a deletion.

**Stage 3d terms.** CNVpytor 1.3.2 is MIT-licensed. INSurVeyor 1.1.3 is GPL-3.0 and Gauchian
1.0.2 GPLv3 at its tag (only its `master` moved to PolyForm in 2024). KILDA 1.2.3 is
AGPL-3.0, and Scramble 1.0.2 CC BY-NC-SA 4.0: personal non-commercial use, as here, is within
it. Each is run here for you, never redistributed, and KILDA is never served. ClinGen's
gene-disease validity export (the Tier 4 list, dated 2026-10-30) and its region curations and
recurrent CNV regions are CC0. Every non-normal 3d result is a screening result that needs
confirmation, never a diagnosis. Not covered: depth events under about 1 kb (most single-exon events), target
exons in segmental duplications without a Parascopy model, GSTT1, UGT2B17 (Parascopy's
UGT2B15 model covers 21% of its MANE span), where an extra whole-gene copy lies, and
inversions and translocations (short-read inversion recall is 0.12-0.17 at 30x).

**Slice 3f terms.** haplocheck 1.3.3 and RBCeq2 2.4.4 are MIT. fastmixture 1.3.0 is
GPL-3.0, and hmmix 0.9.2's wheel ships the GPL-3.0 text while PyPI and its repository say
MIT: both are run here for you, never redistributed. hmmix's helper files and the
archaic genomes are CC BY 4.0 (Zenodo records 11212339 and 13368126). MITOMAP's Confirmed
list (`src/genokeeper/resources/mito/`, read 2026-10-30) is CC BY 3.0: cite MITOMAP as its
CitingMitomap page asks. GenCC's submissions are CC0. ClinPGx's level 3 annotations come
from the same download as the levels above, CC BY-NC-SA 4.0: personal use only, never
redistributed. An HLA allele, the coeliac screen, a MITOMAP variant and HHV-6 are screening
results that need confirmation; the ancestry, blood group, HHV-6 band and trait statements
are information, never a diagnosis.

**GWAS Catalog terms.** The NHGRI-EBI GWAS Catalog's associations are available under the
general Terms of Use for EMBL-EBI Services (its licensing page, read 2026-10-21); they stay
in the lab and the store, never redistributed.

## Asking about the genome

Open Claude Code in this repository and ask in plain words ("am I a carrier of anything?",
"what is in BRCA1?"). The `genome-ask` skill (`.claude/skills/genome-ask/SKILL.md`) reads the
findings store and answers with each claim's genotype reliability, clinical validity, effect
and sources, plus what the filters hid. A claim at a site whose genotype is unknown is never a
finding: the query counts these `unassessed` claims per gene, with why the genotype is
unknown, and lists them on request. What a calibration table withholds stays one of those
claims, never a finding, a watch item, a normal or an exclusion, and the dashboard's counts
do not change; where the tool measured a value, the claim keeps it for the skill under
`withheld`, at the lowest reliability, with the table's reason, what it would mean were it
true and the test that would confirm it. Carrier and "do I have" questions go through `genokeeper
genotype`: an answer per site from its dosage and reliability (`carrier`, `not_a_carrier`,
`probably_not_a_carrier`, `unknown`; `possible_carrier` for a screening result that needs
confirmation, every stage-3 finding among them; `needs_confirmation` for another deletion,
depth event or insertion; `information` for a copy number or Lp(a); "not assessed" when the
store holds no site), and for a
gene counts, what was not looked at, and the gene's main mechanism when SNV calls cannot
see it (SMN1, HBA1/HBA2, CYP21A2, GBA1, F8, DMD, FMR1), and what the deletion screens
cannot see in it (CBS: a large deletion over gnomAD-SV's common depth-only records there,
`gene.deletion_screen_limit`); for SMN1 and the repeat genes
the gene's stage-3 statement leads the sites with its own answer (the copy number or the
repeat count against its thresholds), and the mechanism flag is cleared when that statement
is assessed. CFTR's poly-T/TG (stage 3e) is printed as its own line, never as
cystic-fibrosis carrier status. Population frequency is assessed
for a run whose versions list `gnomad_faf`: there a pathogenic claim common in gnomAD v4.1
(BA1, or its gene's BS1) is graded `conflicting`. Drug questions read `--topic drugs`
(PharmCAT's gene results, PharmGKB's annotations, ClinVar's drug responses); the skill
never gives a dose or a drug choice and sends you to a prescriber.
The same data from the command line:

```bash
uv run genokeeper query                      # default view: reliability >= medium, validity >= moderate
uv run genokeeper query --gene BRCA1 --json  # the JSON contract the skill reads
uv run genokeeper query --min-reliability low --min-validity conflicting  # everything
uv run genokeeper query --list-unassessed    # also list the sites whose genotype is unknown
uv run genokeeper show clinvar:<variation_id>  # one statement, whatever its kind
uv run genokeeper genotype --gene CFTR       # per-site answers and the gene's counts
uv run genokeeper genotype --marker rs80357906 --json
```

The store is `reports/genomelab.duckdb` by default (`--db` moves it; `--run` picks a run
other than the latest).

## Dashboard and explanations

```bash
uv run genokeeper explain --dry-run   # each prompt that would be sent, and the call count
uv run genokeeper explain             # ask Claude for the findings and watch items not explained yet
uv run genokeeper dashboard           # write reports/dashboard.html
```

`genokeeper dashboard` writes one static page of every finding (topic chips, the two sliders,
what the filter hides and why, expandable rows), built from the query alone. It loads
nothing from the network: its Content-Security-Policy allows only its own inlined script and
style. Its search box narrows the current view to findings whose gene, condition,
classification, id or drug names (what `genokeeper query --text` searches), or whose cached
explanation, contain the typed text, ignoring case and reading `_` as a space.

`genokeeper explain` writes, per finding, a short "what it means" and a plain headline, in
English unless `--lang ru` or `GENOKEEPER_LANG=ru` asks for Russian. The headline leads the
finding's row on the page above the source's own words. A stage-3
watch item (a borderline or intermediate repeat, an SMN1 silent-carrier marker) is asked for
too, as a watch item and never as a finding, and the page's stage-3 list shows its "what it
means" (not the headline) under its call, or nothing until one is written. It runs
headless Claude Code (`claude -p`, model `claude-opus-5-5`, no tools, no MCP servers, no
settings). **`explain` and `refresh --explain` are the commands that send finding data (findings and stage-3 watch items) off
the machine**; a `refresh` without the flag sends nothing, and its new findings show without
an explanation until one of the two runs. Each prompt carries to Anthropic one finding's facts — its id, topic and gene; its marker
(rsID, GRCh38 position); its genotype (the caller's `gt`, the allele and its dosage); its
validity (level, the source's own level, the reason); whether it is actionable; its effect; a
gene call's diplotype, phenotype, activity score, callers, level and warnings; the genotype
call's reliability warnings (not `p_wrong`, GQ or flags); on chrX, chrY and chrM a note on
what the dosage means; and for a stage-3 watch item its kind (`watch`), with a note in the
instructions that it is not a finding — and its sources, each as its id and text, with no version (the claim
itself, from which any absolute file path is left out, and the MedlinePlus Genetics summaries
of its gene and conditions). Nothing leaves with `explain --dry-run`, which prints those
prompts instead. An answer is kept only when it is the JSON object asked
for, with a headline, and cites at least one of its prompt's sources; a rejected or failed
one is counted and not kept. It asks four findings at a time (`--jobs` sets how many; 1 asks
one after another), and each finding asked prints one stderr line, in the order the answers
arrive and tagged with the finding's id, with its outcome, and a rejected answer's text or a
failed call's reason (an expired login: run `claude` and log in again; a usage or rate limit
shows as `API error 429 …`, after Claude Code's own retries). The command exits 1 when a call
failed; a rerun asks only the findings without an answer.
Accepted answers are cached in `reports/explanations/` (`--cache` moves it) under a sha256 of
the model, the prompt version and the prompt, so a rerun asks only for findings whose claim,
sources' text, model or prompt changed (a new release of a source alone does not). The page shows an explanation only while its key
still matches, and otherwise "No explanation generated yet" (in Russian, "Объяснение ещё не
сгенерировано"). The command refuses without
the bootstrap's MedlinePlus file.

## Refreshing

A refresh grades through `genokeeper run`, so each new run carries the research polygenic
scores and `genokeeper changes` counts them under research, listing one only when named. A
failure anywhere in the research run (a missing scoring file names `scripts/bootstrap.sh`)
leaves every research score unassessed, saying what failed: the run and the refresh go on.

```bash
uv run genokeeper refresh             # newest releases, a run when one is due, the Changes tab
uv run genokeeper refresh --explain   # the same, and explain the new findings and watch items (sends them)
uv run genokeeper changes             # what changed since the previous refresh
```

**History: once after stage 4 merged, a lab built before it ran `uv run genokeeper run` by
hand.** VEP's step key and the explanations' prompt version both changed, so that one run
annotated again with VEP and SpliceAI (about 4.5 h of SpliceAI in one process, about half that
in its shards: the README's figure, which no measurement in the code records) and `genokeeper explain` then asked every explanation
again. The refreshes after it reuse VEP's cache. `SPLICEAI_CORES` in
`src/genokeeper/annotate/vep.py` records the measurement the shards are sized by: on 1,757 of
HG003's sites (2026-10-29), one SpliceAI process took 1,224 s and four took 597 s.

`genokeeper refresh` holds `reports/work/.refresh.lock` (a second refresh is refused) and runs
under umask 077. When the analyser (a hash of the package's code, its resources and
`config/grading.toml`) differs from the latest snapshot's, it first runs `genokeeper run` on the
releases already on disk. It then runs `scripts/bootstrap.sh --refresh STEP` once for each
refreshed source (ClinVar, ClinGen validity and dosage, CPIC, ClinPGx, the GWAS Catalog,
MedlinePlus, LitVar2 and ClinVar's citations) and records each one new, unchanged or failed.
It runs `genokeeper run` again only when a grading input's newest file has a sha256 other than
the latest snapshot recorded (ClinVar, ClinGen validity and dosage, CGD, CPIC, ClinPGx, the
GWAS Catalog), when the lab built a table since (the gnomAD FAF table, or a calibration table
`genokeeper calibrate` wrote), or the analyser still differs. A run re-grades from the cached calls: it never
re-aligns, re-calls or re-runs VEP and SpliceAI for a new ClinVar release, and it never takes
a pinned component.

It then writes the PMIDs LitVar2 and ClinVar's citations list for every trusted finding and
watch-list VUS to `reports/refresh/<started>/publications.json` (one list per refresh), with `--explain` asks for the
explanations that are not cached yet, and writes the dashboard with its Changes tab. A source that failed is
named in the tab's banner with the release its rows still date from, never shown as "no
change". Retention then deletes each refreshed source's dated folder that is not its newest
and that neither of the last two refreshes' runs read (MedlinePlus, LitVar2 and the citations
keep their newest only).

Each refresh appends one line to `~/genome-lab/refresh/log.jsonl`: the start and end, each
source's release before and after with its status, the runs, and the outcome. The line holds
an error's class only; its full text, which can name variants, goes to
`reports/refresh/<started>.error.txt`. The exit code is 0 when complete, 3 when partial (a
source failed), and 1 when failed (a run refused, or another refresh under way). A failed
refresh is never a comparison base: the next refresh compares with the last complete or
partial one, and runs even when the release it needs is already on disk.

### The daily update check

```bash
brew install terminal-notifier       # once; open its app once for the permission prompt
uv run genokeeper schedule install    # writes and loads the LaunchAgent, sends a test notification
uv run genokeeper schedule uninstall  # unloads it and removes the plist
uv run genokeeper check-updates       # the check by hand; --notify announces what is news
```

`genokeeper check-updates` asks each source for its newest version at a URL written in
`src/genokeeper/check.py` (a checksum file, a header, a listing or a small JSON document; no
request depends on your data) and writes `~/genome-lab/refresh/check.json`: per source its
newest version, the one here, and newer, same or failed. A source it could not check, or one
whose answer is not the shape of that source's version (a captive portal's page, a garbled
tag), is "could not check", never "up to date". A pinned component (PharmCAT, IPD-IMGT/HLA,
Ensembl VEP, Haplogrep3, STRchive, ExpansionHunter (its repeat catalog), Yleaf, PGS Catalog, gnomAD, CGD) and a ClinGen CSpec specification of a
`[bs1_by_gene]` gene (`src/genokeeper/resources/cspec_versions.tsv`) are reported when newer,
never applied; edit a threshold by hand with its source, and the specification's version
beside it.

`genokeeper schedule install` writes `~/Library/LaunchAgents/local.genokeeper.check-updates.plist`
(10:00 every day and at login, absolute paths, the main checkout as the project), loads it
with `launchctl bootout`, `enable` and `bootstrap`, and sends one test notification; it
removes the agent again and fails when launchd refuses a step or the notifier's exit status
is not 0. Approve the agent once under System Settings → General → Login Items & Extensions
→ "Allow in the Background"; this shows it allowed:

```bash
sudo sfltool dumpbtm | grep genokeeper
```

The project was called genomelab before 2026-11. An agent installed then
(`local.genomelab.check-updates`) runs a command that no longer exists: `genokeeper schedule
install` removes it and installs the new one, so run it once after updating. Each
`GENOKEEPER_*` variable is still read under its old name, `GENOMELAB_*`, when it is unset.

A notification names source names, versions and days, never a gene or a variant, and fires
once per new set of newer versions: when the last refresh is 30 days old or more and a
refreshed source is newer, or when a pinned component or a specification has a version not
announced before. The dashboard and `genokeeper changes` say "the update check has not run for
N days" past 35.

## Data

Raw files live in `data/raw/` and **never enter git**: `.gitignore` excludes them, and
`scripts/check_no_genome_data.sh` refuses genomic paths, data files (JSON, HTML, SVG, text,
archives, anything but source, docs and config) outside the data folders (`tests/fixtures/`,
`src/genokeeper/resources/`, `docker/`, `docs/demo/`), files over 5 MB, the private IDs you list in
`.git/info/genome-private-ids` (your kit or sample ID, one extended regular expression per
line; the file is never committed), and a VCF that names a sample other than a public one (`HG…`, `NA…`: the pipeline names every
sample `sample`, so a fixture made from HG002's run renames its column to `HG002`). It runs as
the pre-commit hook, as the pre-push hook over the outgoing commits, and in CI. After cloning,
enable the hooks:

```bash
git config core.hooksPath .githooks
```

Claude Code sessions in this repository cannot read the genome data with their file tools:
`.claude/settings.json` denies reads of `data/raw/` and `reports/` (in this checkout and,
for worktree sessions, wherever on disk they sit) and of `~/genome-lab/calibration/`. Its
anchor is the harness's permission check, which an implementer subagent cannot change
without the user's approval prompt. It guards against accidental reads, not against a
determined shell: a script that opens the files itself is not stopped, and neither is a
search that reads without naming the file (`grep -r`, `find … -exec`); the `genome-ask`
skill's privacy rule sends searches through the Grep tool, `rg` or `git grep` instead. Settings files are
strict JSON, so this paragraph carries the rule's reason; `$GENOKEEPER_HOME` moves the lab
out of the rule's reach.

A typical delivery (tellmeGen's) holds:

| File | Content |
|---|---|
| `*_1.fq.gz`, `*_2.fq.gz` | raw paired reads: the only files `genokeeper run` reads |
| `*.md.cram` + `.crai` | the delivery's reads aligned to GRCh37/hg19 (not used) |
| `*.vcf.gz` | the delivery's DeepVariant 1.10 calls on GRCh37 (not used) |
| `*.vc.txt`, `*.txt` | 23andMe-format genotypes |

## Layout

```
data/raw/     raw delivery (never in git)
reports/      run outputs (never in git: they contain genotypes)
docs/demo/    the README's demo, made from an invented genome
scripts/      guards, setup and run scripts
```

Tools and annotation databases are installed outside the repository, under `~/genome-lab/`.

