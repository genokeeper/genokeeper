"""The run: stage 2a's GRCh38 pipeline (`run_genome`, `genokeeper run`).

`run_genome` aligns the delivered FASTQ to the lab's GRCh38 FASTA, gates the sample on QC
before calling, calls it with DeepVariant into a gVCF, and grades every claim of a ClinVar
GRCh38 release (`panel.build.build_claims`) on our own call at its allele
(`genotype.gvcf.calls_at`; on chrM, mutserve's, `genokeeper.mito`).

It works in two phases. The first does every slow step outside the store, each one's artifact
cached under `reports/work/<step>/<sha256 of the step's inputs>/` and reused while those
inputs are unchanged. The second writes the run into the store in one transaction: the run
row; then each source (`genokeeper.sources`) collects its variants, calls, statements and
accounting, the compound-heterozygote promotion runs over all of them as a pure step with
the candidate pairs' read-backed phases, and `store_sources` writes them; then a check that
every statement's kind agrees with the call stored beside it (`check_calls`). The per-site
rows go through a `StoreWriter`, in bulk. An interrupted run leaves no rows, so a query of
the latest run never reads half of one. Every run writes its own complete rows under a new
`run_id`; `run.inputs` records each step's input hash and whether its artifact was reused.
The snapshot is exported after the commit.

Gene-disease validity comes from ClinGen's own export (`<lab>/clingen/<YYYYMMDD>/`, the
newest), and inheritance for any gene from the CGD table (`<lab>/cgd/<YYYYMMDD>/`), both
required. Population frequency comes from the lab's gnomAD FAF table (`<lab>/gnomad/4.1/`,
`genokeeper.evidence.frequency`), which every run requires and opens before any step;
`run.inputs` records its manifest's sha256. The source looks its vids up once, when its
statements are graded in the store phase, which is never cached: no step's key holds the
table, and a rebuilt table changes the next run's grades without rebuilding a step.

Genotype error comes from the calibration table (`<lab>/calibration/small_variants/`,
written by `genokeeper calibrate`, `genokeeper.calibrate.tables`): each call's `p_wrong` is the
rate measured on HG003 for its cell, inside HG003's benchmark regions. Every run requires the
table and refuses, before any step, one measured under other pipeline versions than its own
(`calibration_versions`, `open_calibration`); like the FAF table, it is in no step's key.

Pharmacogenomics (slice 1b) comes after the calls, in every run, from the bootstrap's
slice-1b files (`pgx_inputs`, refused before any step when one is missing): PharmCAT on our
calls at each of its positions (`pharmcat_step`) given the outside callers' calls
(`outside_steps`, `outside_calls`), stored by `PgxSource`; and PharmGKB's clinical
annotations located at their ClinVar alleles and graded on our calls there, stored by
`PharmGkbSource`. CPIC's gene-drug levels also make a carried ClinVar drug_response
actionable. No dosing text enters the store; PharmCAT's HTML report stays in its step's
directory under `reports/work/pharmcat/` for the user. An outside caller that fails (exits
non-zero, or writes no result its reader takes) leaves the genes it serves unresolved, the
failure as the reason, and the run goes on; nothing of its step is cached.

Polygenic scores (slice 1d) come next, from the bootstrap's slice-1d files (`prs_inputs`,
refused before any step when one is missing): our calls at every chosen score's and the
HGDP+1kGP panel's ancestry positions (`prs_genotypes_step`), graded and written as
pgsc_calc's target VCF, and pgsc_calc on it (`pgsc_step`, `genokeeper.prs.pgsc`), whose
results a task it cannot show ran offline in a pinned image refuses; stored by `PrsSource`
with APOE's diplotype from the graded ClinVar claims at its two SNPs (a double
heterozygote phased with the candidate pairs). The mitochondrial haplogroup comes from
Haplogrep3 on mutserve's homoplasmic chrM calls (`genokeeper.haplogroup`), stored by
`MitoSource`; the research associations from the GWAS Catalog's download, located on the
panel's variants (`gwas_alleles_step`) and our calls there (`genotypes_step`),
stored by `GwasSource`.

Variants ClinVar has not classified (slice 1e) come last, from the bootstrap's slice-1e
files and ClinGen's dosage list (`predicted_inputs`, refused before any step when one is
missing): every PASS allele of our VCF that no record of the ClinVar release classifies
germline, whatever its shape (`clinvar_vids`: filtered at grading, after the VEP cache),
through VEP and SpliceAI (`vep_step`, `genokeeper.annotate.vep`), our calls at the
candidates (`genotypes_step`), stored by `PredictedSource` as predictions
(`genokeeper.evidence.predicted`); two carrier ones in one Tier 3 gene join the candidate
pairs.

Slice 3c checks the run itself. Right after the call, our genotypes at the located GWAS
SNPs are compared with the delivered 23andMe-format file (`genokeeper.sample_identity`), and
a mismatch refuses before the claims are graded (a FASTQ mix-up check: the provider derived
that file from the same reads). Every graded call passes the allele-fraction guard
(`reliability.guard_allele_fraction`, capped in the clonal genes' spans); each coverage
gene's callability goes into the store's `coverage` table, and the sex-chromosome dosage
is one statement (`DosageSource`), both from `foundation_steps`.

Slice 3f types what the reads say beyond the variants: haplocheck's contamination check
and GenCC's validity in every run, MITOMAP's confirmed chrM variants and eleven traits (3f's
five, 3g's six) on the run's calls, and, with a current slice 3f calibration table
(`genokeeper calibrate`), HLA, HHV-6, runs of homozygosity, blood groups and ancestry
(`typing_steps`), each gated statement unassessed with the reason without one.

`genokeeper run` runs it on this machine's files under umask 077, and exits 2, with a
one-line message, on a missing tool or input (naming the bootstrap), a platform the T7
checkpoint was not trained on, a failed QC gate, an ambiguous sex, Docker not running or an
image missing, PharmCAT, pgsc_calc or Haplogrep3 failing, or the containment check refusing
pgsc_calc's results. It prints the run id, the accounting table (unassessed
claims beside findings, watch items and non-findings), the QC metrics, the population frequency table it graded with and the path
of PharmCAT's HTML report: counts only, never a genotype.
"""

import argparse
import csv
import fcntl
import functools
import gzip
import hashlib
import importlib.resources
import itertools
import json
import math
import os
import re
import secrets
import shutil
import subprocess
import sys
import tempfile
from collections.abc import (
    Callable,
    Collection,
    Iterable,
    Iterator,
    Mapping,
    Sequence,
)
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, replace
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Concatenate, cast

import duckdb
import pysam

from genokeeper import __version__
from genokeeper import coverage as gene_coverage
from genokeeper.align import PlatformRefused, align, check_platform
from genokeeper.ancestry.admixture import (
    pruned_alleles,
    read_groups,
    run_projection,
    to_bed,
)
from genokeeper.ancestry.archaic import (
    ArchaicShare,
    archaic_share,
    read_segments,
    run_hmmix,
)
from genokeeper.annotate.vep import CACHE_VERSION as VEP_CACHE_VERSION
from genokeeper.annotate.vep import (
    VEP_INPUT,
    VEP_VERSION,
    VepFiles,
    VepRecord,
    annotate,
    locus,
    vep_files,
    vep_input,
)
from genokeeper.bloodgroups import (
    POSITIONS_FILE,
    SystemCall,
    rbceq2_positions,
    read_rbceq2,
    run_rbceq2,
)
from genokeeper.bloodgroups import uncovered as uncovered_systems
from genokeeper.calibrate import BENCH_BED, layout, table_dir
from genokeeper.calibrate.classes import ClassTable, open_class_table
from genokeeper.calibrate.deletions import BANDS as DELETION_BANDS
from genokeeper.calibrate.deletions import (
    DELETIONS_UNION_STEP_VERSION,
    DELLY_STEP_VERSION,
    DYSGU_STEP_VERSION,
    deletion_versions,
    open_deletions,
    read_population,
)
from genokeeper.calibrate.hg002_cells import (
    HAPLOID,
    ONLY,
    REGIONS,
    cells_dir,
    open_cells,
)
from genokeeper.calibrate.slice3d import (
    ALDY_CYP2A6_STEP_VERSION,
    CNV_STEP_VERSION,
    CNVPYTOR_STEP_VERSION,
    DYSGU_INS_DUP_STEP_VERSION,
    GAUCHIAN_STEP_VERSION,
    INSERTIONS_UNION_STEP_VERSION,
    INSURVEYOR_STEP_VERSION,
    KILDA_STEP_VERSION,
    PARALOGS,
    PARASCOPY_DEPTH_STEP_VERSION,
    PARASCOPY_LOCI_STEP_VERSION,
    SCRAMBLE_STEP_VERSION,
    ParalogTable,
    cnv_versions,
    insertion_versions,
    open_paralog_table,
    paralog_versions,
)
from genokeeper.calibrate.stage3 import (
    SMN_KEY,
    YLEAF_KEY,
    hg002_truth,
    open_stage3,
    stage3_versions,
    two_x_key,
)
from genokeeper.calibrate.table import MANIFEST, CalibrationTable
from genokeeper.calibrate.tables import NOT_CALIBRATED, Stale, open_table
from genokeeper.call.deepvariant import (
    CALLING_FLAGS,
    CHECKPOINT_FILES,
    checkpoint_path,
    run_deepvariant,
)
from genokeeper.checksum import md5_cached, sha256_cached, sha256_file
from genokeeper.cnv.callers import (
    expected_ploidy,
    read_cnvpytor,
    read_delly_cnv,
    read_junctions,
    run_cnvpytor,
    run_delly_cnv,
)
from genokeeper.cnv.screen import (
    DEPTH_CALLERS,
    MATCH_OVERLAP,
    NamedCnv,
    gene_spans,
    load_named,
    merge_depth,
    named_deletions,
    owned_regions,
    screen_cnv,
    segdup_exons,
)
from genokeeper.config import (
    GRADING_FILE,
    MASKED_REFERENCE,
    DeletionThresholds,
    Grading,
    Paths,
    grch38_fasta,
    load_grading,
    load_paths,
)
from genokeeper.containers import IMAGES, DockerUnavailable, image_missing
from genokeeper.coverage import CoverageRow
from genokeeper.deletions.callers import (
    caller_filter,
    caller_min_length,
    deletions_only,
    filter_expression,
    run_delly,
    run_dysgu,
)
from genokeeper.deletions.exons import (
    GBA1_HOMOLOGY,
    PRIMARY,
    load_exons,
    load_spans,
    locus_span,
    targets,
)
from genokeeper.deletions.files import sv_files
from genokeeper.deletions.screen import (
    CallClass,
    Deletion,
    Screen,
    Shown,
    max_distance,
    population_spans,
    screen_all,
)
from genokeeper.deletions.truvari import BETWEEN_CALLERS, match_callers, read_union
from genokeeper.evidence.cgd import load_cgd, newest_cgd
from genokeeper.evidence.clingen import (
    DOSAGE_FILE,
    RECURRENT_FILE,
    REGIONS_FILE,
    DosageList,
    gene_level,
    gene_validity_date,
    load_curations,
    load_dosage,
    load_recurrent,
    load_regions,
    regions_date,
)
from genokeeper.evidence.clingen import FILE_NAME as GENE_VALIDITY_FILE
from genokeeper.evidence.clinvar import clnsig_class
from genokeeper.evidence.frequency import (
    FafMissing,
    FafTable,
    faf_folder,
    load_ba1_exceptions,
)
from genokeeper.evidence.gencc import gencc_date, load_gencc, newest_gencc
from genokeeper.evidence.gencc import merge as merge_gencc
from genokeeper.evidence.gwas import DOWNLOAD as GWAS_DOWNLOAD
from genokeeper.evidence.gwas import (
    Association,
    load_associations,
    panel_alleles,
)
from genokeeper.evidence.gwas import locate as locate_gwas
from genokeeper.evidence.pgx import CPIC_FILES, load_cpic_alleles, load_cpic_levels
from genokeeper.evidence.pharmgkb import DOWNLOADS as PHARMGKB_DOWNLOADS
from genokeeper.evidence.pharmgkb import (
    PHARMGKB_DROP_REASONS,
    ClinicalAnnotations,
    clinvar_alleles,
    left_by_pharmcat,
    locate,
    read_clinical_annotations,
)
from genokeeper.evidence.predicted import homopolymer_runs, load_escapes, predicted_rule
from genokeeper.genelists import LISTED, GeneList, load_gene_lists, screen_overrides
from genokeeper.genotype.clonal import ClonalIndex
from genokeeper.genotype.gvcf import Allele, Called, calls_at
from genokeeper.genotype.reliability import (
    ExtraCells,
    Reliability,
    StrataIndex,
    depth_medians,
    grade_called,
    guard_allele_fraction,
    on_haploid_contig,
)
from genokeeper.haplocheck import OUT_NAME as HAPLOCHECK_OUT
from genokeeper.haplocheck import RAW_SUFFIX as HAPLOCHECK_RAW
from genokeeper.haplocheck import Haplocheck, demoted, read_haplocheck, run_haplocheck
from genokeeper.haplogroup import (
    HAPLOGREP_VERSION,
    TREE,
    VCF_NAME,
    Haplogroup,
    read_haplogroup,
    run_haplogrep,
    too_uncovered,
    uncovered,
    write_chrm_vcf,
)
from genokeeper.haplogroup import OUTPUT as HAPLOGROUPS
from genokeeper.hhv6 import (
    Hhv6Call,
    align_viral,
    extract_unmapped,
    read_coverage,
    viral_coverage,
)
from genokeeper.hhv6 import classify as classify_hhv6
from genokeeper.identity import analyser_identity
from genokeeper.insertions import (
    INSERTION_CALLERS,
    INSERTION_WORDING,
    INSURVEYOR_VERSION,
    MEI_SCOPE,
    SCRAMBLE_MARGIN,
    SCRAMBLE_VERSION,
    add_scramble,
    insertion_union,
    read_scramble,
    run_insurveyor,
    run_scramble,
    scramble_regions,
)
from genokeeper.masking import masked_bed
from genokeeper.mito import (
    CALLS_FILE,
    VCF_FILE,
    MitoCall,
    call_chrm,
    chrm_coverage,
    grade_mito,
    mito_called,
    read_calls,
)
from genokeeper.panel.build import Claim, ClaimSet, build_claims, clinvar_vids
from genokeeper.paralogs import (
    GAUCHIAN_VERSION,
    KILDA_VERSION,
    LOCI,
    PARALOGS_LABEL,
    GauchianCall,
    PsCall,
    aldy_copies,
    depth_ratio,
    read_gauchian,
    read_kilda,
    read_locus,
    run_aldy_cyp2a6,
    run_gauchian,
    run_kilda,
    unique_part,
)
from genokeeper.pgx.cyp2b6 import cyp2b6_outside
from genokeeper.pgx.cyp2d6 import (
    ALDY,
    ALDY_VERSION,
    BCYRIUS,
    BCYRIUS_ID,
    BCYRIUS_VERSION,
    agree,
    read_aldy,
    read_bcyrius,
    run_aldy,
    run_bcyrius,
)
from genokeeper.pgx.hla import (
    FASTQ,
    IMGT_HLA_VERSION,
    MHC_REGION,
    OPTITYPE,
    OPTITYPE_PREFIX,
    PGX_HLA,
    T1K,
    T1K_PREFIX,
    HlaCall,
    extract_hla_reads,
    hla_outside,
    read_optitype,
    read_t1k,
    run_optitype,
    run_t1k,
    t1k_index,
)
from genokeeper.pgx.hla_typing import TAG_ALLELES, TAG_RSIDS
from genokeeper.pgx.mtrnr1 import mtrnr1_outside
from genokeeper.pgx.outside import (
    CallerCall,
    CallerFailed,
    OutsideCall,
    outside_lines,
)
from genokeeper.pgx.pharmcat import (
    PHARMCAT_VERSION,
    REPORT_JSON,
    UGT1A1_28,
    PharmcatFailed,
    gene_regions,
    parse_report,
    pharmcat_positions,
    positions_path,
    run_pharmcat,
    write_pharmcat_vcf,
)
from genokeeper.phasing import Pair, Phase, phase_candidates, phase_pairs, phase_sites
from genokeeper.prs.pgsc import (
    APOE_E2,
    APOE_E4,
    NEXTFLOW_VERSION,
    NF_SCHEMA_VERSION,
    PGSC_CALC_COMMIT,
    PGSC_CALC_VERSION,
    PGSC_LABELS,
    ContainmentRefused,
    PgscFailed,
    PopSimilarity,
    ancestry_pvar,
    java_home,
    nextflow_path,
    nxf_home,
    panel_path,
    panel_pvar,
    parse_results,
    parse_similarity,
    pgsc_calc_dir,
    run_pgsc_calc,
    scoring_file,
    target_positions,
    write_target_vcf,
)
from genokeeper.prs.scores import Score, load_not_covered, load_scores
from genokeeper.qc import judge, mosdepth, run_qc, vcf_stats
from genokeeper.query import GENOTYPE_UNKNOWN
from genokeeper.regions import (
    QC_AUTOSOMES,
    SEX_REGIONS,
    Interval,
    giab38_dir,
    merge,
    read_bed,
    read_named,
    subtract,
)
from genokeeper.repeats.expansionhunter import (
    EH_VERSION,
    REVIEWER_VERSION,
    catalog_sha256,
    flank_ratio,
    review,
    run_expansionhunter,
)
from genokeeper.repeats.expansionhunter import read_calls as read_repeat_calls
from genokeeper.repeats.expansionhunter import region as repeat_region
from genokeeper.repeats.loci import RepeatCall, Rfc1Call
from genokeeper.repeats.rfc1 import RFC1_LOCUS, read_rfc1
from genokeeper.roh import Run, read_runs, run_roh, sites_vcf, summarise
from genokeeper.sample_identity import (
    IDENTITY_TOOL,
    Verdict,
    compare,
    find_delivered,
    read_delivered,
)
from genokeeper.sample_identity import judge as judge_identity
from genokeeper.sex import Dosage, Sex, SexRefused, copies, measure_ratios, x_het_share
from genokeeper.sex import infer_sex as infer_sex_from_cram
from genokeeper.smn import _MODEL as SMN_MODEL
from genokeeper.smn import LABEL as SMN_LABEL
from genokeeper.smn import (
    PARASCOPY_VERSION,
    SMNCALLER_VERSION,
    ParascopyCall,
    SmnCaller,
    parascopy_cn,
    parascopy_data,
    parascopy_depth,
    parascopy_md5s,
    read_parascopy,
    read_smn_caller,
    run_smn_caller,
)
from genokeeper.sources import (
    Collected,
    Row,
    RunContext,
    Source,
    promote,
    write_collected,
)
from genokeeper.sources.ancestry import AncestrySource
from genokeeper.sources.bloodgroups import BloodGroupSource
from genokeeper.sources.cnv import CNVPYTOR_VERSION, CnvSource
from genokeeper.sources.deletions import (
    DELLY_VERSION,
    DYSGU_VERSION,
    TRUVARI_VERSION,
    DeletionSource,
    Rate,
)
from genokeeper.sources.dosage import DosageSource
from genokeeper.sources.gwas import GWAS_TOOL, GwasSource
from genokeeper.sources.gwas import graded as graded_gwas
from genokeeper.sources.hhv6 import Hhv6Source
from genokeeper.sources.hla import HlaSource
from genokeeper.sources.mito import MitoSource
from genokeeper.sources.mitomap import MitomapSource
from genokeeper.sources.panel import ClaimsSource
from genokeeper.sources.paralogs import (
    GBA1_CONTRADICTED,
    LpaSource,
    ParalogSource,
    gba1_contradicted,
)
from genokeeper.sources.pgx import PgxSource
from genokeeper.sources.pharmgkb import PHARMGKB_RESEARCH_TOOL, PharmGkbSource
from genokeeper.sources.predicted import PredictedSource, clinvar_truncations
from genokeeper.sources.prs import PrsSource
from genokeeper.sources.repeats import CFTR_GENE, REPEATS_TOOL, RepeatSource, XDosage
from genokeeper.sources.roh import RohSource
from genokeeper.sources.smn import SMN_TOOL, SmnSource
from genokeeper.sources.traits import TRAIT_ALLELES, TraitSource
from genokeeper.sources.yleaf import YLEAF_TOOL, YleafSource
from genokeeper.store.db import (
    Kind,
    Level,
    RowType,
    StoreTooNew,
    StoreTooOld,
    check_schema_version,
    connect,
    export_snapshots,
    start_run,
    transaction,
)
from genokeeper.store.writer import StoreWriter, sweep_rows
from genokeeper.tools import ToolFailed, ToolMissing, run_tool
from genokeeper.typing_files import typing_files
from genokeeper.yhaplogroup import FOLDER as YLEAF_FOLDER
from genokeeper.yhaplogroup import PREDICTION as YLEAF_PREDICTION
from genokeeper.yhaplogroup import YLEAF_VERSION, YHaplogroup, read_yleaf, run_yleaf

if TYPE_CHECKING:
    from genokeeper.calibrate.slice3f import Slice3fTable
    from genokeeper.cli import Subparsers

# The store keeps its name from before the rename to genokeeper: the user's runs are in it.
STORE_NAME = "genomelab.duckdb"
# Under reports/work: held (flock) by the one run that may write the store and the cache.
LOCK_NAME = ".run.lock"
_CLINVAR_DATE = re.compile(r"\d{8}")

# Each cached step's version, part of its artifact's key: bump when the step's output
# changes (another tool option, a new field, a fix in how it is built), so an artifact the
# old code built is rebuilt rather than reused. A step's key also holds the versions of the
# tools it runs and the keys of the steps it reads, so a bump rebuilds those after it too.
# The genokeeper version is not in a key: a release that changes no step's output rebuilds
# nothing.
ALIGN_STEP_VERSION = "1"
QC_STEP_VERSION = "2"  # 2: total_bases
SEX_STEP_VERSION = "1"
CALL_STEP_VERSION = "1"
CALL_QC_STEP_VERSION = "1"
GENOME_DEPTH_STEP_VERSION = "1"
CLAIMS_STEP_VERSION = "3"  # 2: structural_variant; 3: consequence (MC)
GENOTYPES_STEP_VERSION = "1"
MITO_STEP_VERSION = "2"  # 2: mutserve's VCF too, for haplocheck (slice 3f)
HAPLOCHECK_STEP_VERSION = "1"
HAPLOGROUP_STEP_VERSION = "1"
PHASES_STEP_VERSION = "2"  # 2: the run's het SNVs between a pair's sites
PGX_ALLELES_STEP_VERSION = "2"  # 2: the level-3 annotations' rsIDs (slice 3f)
PGX_GENOTYPES_STEP_VERSION = "2"  # 2: the level-3 annotations' alleles (slice 3f)
PGX_PHASES_STEP_VERSION = "1"
BCYRIUS_STEP_VERSION = "1"
ALDY_STEP_VERSION = "1"
HLA_READS_STEP_VERSION = "1"
OPTITYPE_STEP_VERSION = "1"
T1K_STEP_VERSION = "1"
PHARMCAT_STEP_VERSION = "2"  # 2: phase sets and trusted alleles (slice 3f)
PRS_GENOTYPES_STEP_VERSION = "1"
PGSC_STEP_VERSION = "3"  # 2: --min_overlap 0; 3: no --keep_multiallelic
GWAS_ALLELES_STEP_VERSION = "1"
GWAS_GENOTYPES_STEP_VERSION = "1"
CLINVAR_VIDS_STEP_VERSION = "3"  # 3: every germline-classified record, any shape
VEP_STEP_VERSION = "5"  # 5: records per variant and gene (BR-012, BR-013)
PREDICTED_GENOTYPES_STEP_VERSION = "1"
TYPING_GENOTYPES_STEP_VERSION = "1"
REPEATS_STEP_VERSION = "1"
SMN_STEP_VERSION = "2"  # 2: Parascopy's depth from the shared `parascopy_depth` step
YLEAF_STEP_VERSION = "1"
COVERAGE_STEP_VERSION = "1"
DOSAGE_STEP_VERSION = "1"


class RunRefused(RuntimeError):
    """The run cannot start or go on: a missing tool or input, another run under way, a
    store another process holds or another schema version wrote, a platform the T7
    checkpoint was not trained on, a failed QC gate, an ambiguous sex, a FASTQ that does
    not match the delivered genotypes, or Docker not running. `genokeeper run` exits 2 with
    the message."""


def _one_line(error: BaseException) -> str:
    return " ".join(str(error).split())


@dataclass(frozen=True)
class _Step:
    """A cached artifact: the hash of its step's inputs, its path, and whether it was
    reused rather than built by this run."""

    key: str
    path: Path
    reused: bool

    def record(self) -> dict[str, object]:
        return {"key": self.key, "path": str(self.path), "reused": self.reused}


def _key(*parts: str) -> str:
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def _cached(
    work: Path,
    step: str,
    version: str,
    parts: Sequence[str],
    name: str,
    build: Callable[[Path], None],
) -> _Step:
    """`<work>/<step>/<key>/<name>`, where `build(directory)` writes `name` and anything
    beside it; built in a fresh directory renamed into place, so a directory that exists
    is complete and is reused as it is, unless `name` no longer resolves there (a symlink
    whose target moved): that directory is removed and built again.

    The key hashes the step's `version` (its `*_STEP_VERSION`) and the step's inputs, not
    the code nor the genokeeper version: after changing how a step builds its artifact, bump
    its step version.
    """
    key = _key(step, version, *parts)
    final = work / step / key
    if final.is_dir():
        if (final / name).exists():
            return _Step(key=key, path=final / name, reused=True)
        shutil.rmtree(final)
    staging = work / step / f".building-{secrets.token_hex(4)}"
    staging.mkdir(parents=True)
    try:
        build(staging)
        if not (staging / name).exists():
            raise RuntimeError(f"step {step} wrote no {name}")
        try:
            staging.rename(final)
        except OSError:
            # Another run renamed its complete directory into place first: use that one.
            if not (final / name).exists():
                raise
            shutil.rmtree(staging, ignore_errors=True)
            return _Step(key=key, path=final / name, reused=True)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return _Step(key=key, path=final / name, reused=False)


# What a step stages under the work folder, or a tool beside the CRAM it reads; a killed
# run leaves them, and no run reuses one.
_STAGING = (
    "*/.building-*",
    ".pharmcat-input-*",
    ".pgsc-input-*",
    "align/*/.qc-*",
    "align/*/.sex-*",
    "align/*/.phase-*",
)


@contextmanager
def run_lock(work: Path) -> Iterator[None]:
    """Hold `<work>/.run.lock` for the block, having removed a killed run's staging; refuse
    (`RunRefused`) while another run holds it."""
    work.mkdir(parents=True, exist_ok=True)
    lock = work / LOCK_NAME
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RunRefused(
                f"another run is under way (it holds {lock}); wait for it to finish"
            ) from None
        # Its anchor comes from this flock, which a live run cannot lose because the kernel
        # drops it only when that run's descriptor closes; every staging writer holds it.
        for folder in sorted(p for glob in _STAGING for p in work.glob(glob)):
            shutil.rmtree(folder, ignore_errors=True)
        yield
    finally:
        os.close(fd)


def _single_writer[**P, R](
    run: Callable[Concatenate[Paths, P], R],
) -> Callable[Concatenate[Paths, P], R]:
    """`run(paths, ...)` under `paths.reports/work/.run.lock`: one run at a time writes
    the store and the step cache."""

    @functools.wraps(run)
    def locked(paths: Paths, *args: P.args, **kwargs: P.kwargs) -> R:
        work = paths.reports / "work"
        with run_lock(work):
            # A killed run's row files hold the user's rows; no other run is writing now.
            for failure in sweep_rows(work):
                print(
                    f"genokeeper run: a killed run's rows were not removed: {failure}",
                    file=sys.stderr,
                )
            return run(paths, *args, **kwargs)

    return locked


# --- checks -------------------------------------------------------------------------------


def _check_files(files: Mapping[str, Path]) -> None:
    for what, path in files.items():
        if not path.is_file():
            raise RunRefused(f"missing input {what}: {path}")


def _in_use(store: Path, error: duckdb.Error) -> RunRefused:
    return RunRefused(
        f"the findings store {store} is in use by another process; close it (a"
        f" genokeeper query, a DuckDB shell) and run again: {_one_line(error)}"
    )


def _check_store(store: Path) -> None:
    """Refuse (`RunRefused`), before any step runs, a store another schema version wrote
    or another process holds for writing; read-only, so the store is left as it is. A
    process holding it read-only (a `genokeeper query`) passes here, as DuckDB lets a second
    reader in, and is refused only when the run opens the store to write. A missing store
    is fine: the run creates it."""
    if not store.is_file():
        return
    try:
        con = duckdb.connect(str(store), read_only=True)
    except duckdb.IOException as e:
        raise _in_use(store, e) from e
    try:
        check_schema_version(con, store)
    except (StoreTooNew, StoreTooOld) as e:
        raise RunRefused(str(e)) from e
    finally:
        con.close()


def _gene_lists_key(lists: Mapping[str, GeneList]) -> str:
    return json.dumps(
        {name: [gl.version, gl.genes] for name, gl in lists.items()}, sort_keys=True
    )


# --- the store ----------------------------------------------------------------------------


class CallMismatch(RuntimeError):
    """A run's statement kind contradicts the call stored beside it: a bug, never data."""


def check_calls(con: duckdb.DuckDBPyConnection, run_id: str) -> None:
    """Raise `CallMismatch` when a variant finding of the run sits on an unknown, missing or
    zero-copy call, or an unassessed claim on a known one; counts only, never a genotype.

    A GWAS association (`sources.gwas`) is checked instead against its risk allele's copies,
    recounted from the call stored at its vid with the function that graded it
    (`sources.gwas.graded`, ruling R26): its kind and its claim's `risk_allele_copies` must
    both be what that call gives, and the call must be of a countable shape. So it refuses
    an association unassessed on a known call (Task 34's run), one whose call is not
    countable (a reference risk allele beside another allele, or a risk allele that is not
    exactly the vid's REF or ALT), and a kind or count its call contradicts. Its anchor
    comes from the vid, whose REF `make_vid` took from an allele checked against the GRCh38
    FASTA, and from the call row stored at that vid: the claim gives only its risk allele,
    which must equal the vid's REF or ALT, and its count and kind, which are compared with
    those `graded` recounts from the call's GT and dosage, never taken from the claim.

    Every statement is graded on the call stored beside it, so this never fires on a
    correct run."""
    row = con.execute(
        "WITH graded AS (SELECT s.kind, c.vid AS called, c.dosage"
        " FROM statement s LEFT JOIN call c ON c.run_id = s.run_id AND c.vid = s.vid"
        " WHERE s.run_id = ? AND s.vid IS NOT NULL AND s.source IS DISTINCT FROM ?)"
        " SELECT count(*) FILTER (WHERE kind = 'finding'"
        " AND (called IS NULL OR dosage IN ('unknown', '0'))),"
        " count(*) FILTER (WHERE kind = 'unassessed'"
        " AND called IS NOT NULL AND dosage <> 'unknown')"
        " FROM graded",
        [run_id, GWAS_TOOL],
    ).fetchone()
    findings, unassessed = row if row is not None else (0, 0)
    associations = con.execute(
        "SELECT s.vid, s.claim->>'risk_allele',"
        " CAST(s.claim->>'risk_allele_copies' AS INTEGER), s.kind, c.gt, c.dosage"
        " FROM statement s LEFT JOIN call c ON c.run_id = s.run_id AND c.vid = s.vid"
        " WHERE s.run_id = ? AND s.source = ?",
        [run_id, GWAS_TOOL],
    ).fetchall()
    gwas = sum(
        graded_gwas(vid, risk, gt, dosage) != (claimed, kind)
        for vid, risk, claimed, kind, gt, dosage in associations
    )
    if findings or unassessed or gwas:
        raise CallMismatch(
            f"run {run_id}: {findings} variant findings on an unknown, missing or zero-copy"
            f" call and {unassessed} unassessed claims on a known call, and {gwas} GWAS"
            " associations whose risk-allele copies the call contradicts or cannot count;"
            " the run is rolled back"
        )


def store_sources(
    writer: StoreWriter,
    con: duckdb.DuckDBPyConnection,
    ctx: RunContext,
    sources: Sequence[Source],
    phases: Mapping[Pair, Phase] | None = None,
    gauchian: GauchianCall | None = None,
) -> None:
    """The run's statements into the store: each source collects its rows, the statements
    Gauchian's GBA1 call contradicts become unassessed (`gba1_unassessed`), the
    compound-heterozygote promotion runs across all of them with the candidate pairs'
    read-backed `phases` (none: every pair unphased), and they are written with one
    accounting row per tool. A new source is one more item in `sources`."""
    collected = gba1_unassessed([source.collect(ctx) for source in sources], gauchian)
    promoted = promote.promote_compound_hets(
        collected, ctx.lists, ctx.inheritance, ctx.haploid, phases
    )
    write_collected(writer, con, ctx.run_id, promoted)


def _tool_version(tool: str) -> str:
    first = run_tool([tool, "--version"]).splitlines()[0]
    return first.split()[-1]


def _header_value(vcf: Path, key: str) -> str | None:
    with pysam.VariantFile(str(vcf)) as f:
        for rec in f.header.records:
            if rec.key == key and rec.value is not None:
                return str(rec.value)
    return None


def _clinvar_release(clinvar: Path) -> str:
    """The ClinVar release, `##fileDate` as YYYYMMDD (the statements' source version)."""
    date = (_header_value(clinvar, "fileDate") or "").replace("-", "")
    if not _CLINVAR_DATE.fullmatch(date):
        raise RunRefused(f"{clinvar} has no ##fileDate (YYYY-MM-DD or YYYYMMDD)")
    return date


def _file(path: Path, sha: str) -> dict[str, str]:
    return {"path": str(path), "sha256": sha}


# --- the run --------------------------------------------------------------------------------

GENOME_TOOLS = ("bcftools", "samtools", "bwa", "mosdepth")
# The read group's sample name in the CRAM and in DeepVariant's VCFs: never the user's.
GENOME_SAMPLE = "sample"
# Under reports/work: each file's sha256 with the size, mtime and inode it was taken at.
IDENTITIES = "identities.json"
CHRM = "chrM"


def _bwa_version() -> str:
    """bwa's version from its usage text (`Version: 0.7.19-r1273`), which it prints on
    stderr while exiting 1."""
    done = subprocess.run(
        ["bwa"], check=False, capture_output=True, text=True, stdin=subprocess.DEVNULL
    )
    for line in done.stderr.splitlines():
        if line.startswith("Version:"):
            return line.split()[1]
    raise RunRefused("bwa printed no version; reinstall it (brew install bwa)")


def _check_genome_tools() -> None:
    for tool in GENOME_TOOLS:
        if shutil.which(tool) is None:
            raise RunRefused(str(ToolMissing(tool)))


def _bootstrap_missing(what: str) -> RunRefused:
    return RunRefused(f"missing input: {what}; run scripts/bootstrap.sh")


def genome_tool_versions() -> dict[str, str]:
    """The versions of the genome tools a run uses: bcftools, samtools, bwa and mosdepth.
    Refuses (`RunRefused`) when one is not on PATH."""
    _check_genome_tools()
    versions = {tool: _tool_version(tool) for tool in ("bcftools", "samtools")}
    versions["bwa"] = _bwa_version()
    versions["mosdepth"] = _tool_version("mosdepth")
    return versions


def calling_inputs(lab: Path) -> tuple[Path, Path]:
    """The lab's GRCh38 FASTA, its .fai and bwa index beside it, and the T7 checkpoint
    prefix; refuses (`RunRefused`) what is missing, naming the bootstrap."""
    fasta = grch38_fasta(lab)
    for what, path in (
        ("the GRCh38 FASTA", fasta),
        ("its .fai", Path(f"{fasta}.fai")),
        ("its bwa index", Path(f"{fasta}.bwt")),
    ):
        if not path.is_file():
            raise _bootstrap_missing(f"no {what} at {path}")
    checkpoint = checkpoint_path(lab)
    absent = [n for n in CHECKPOINT_FILES if not (checkpoint.parent / n).is_file()]
    if absent:
        raise _bootstrap_missing(
            f"T7 checkpoint files {', '.join(absent)} in {checkpoint.parent}"
        )
    return fasta, checkpoint


def calibration_versions(
    lab: Path, *, tools: Mapping[str, str], identities: Path
) -> dict[str, str]:
    """The six versions a calibration table is measured under (`VERSION_KEYS`), as this
    machine's pipeline has them: bwa's and samtools' from `tools` (`genome_tool_versions`),
    the DeepVariant image digest written in the code, each T7 checkpoint file's md5, the
    DeepVariant arguments that decide the calls (`CALLING_FLAGS`; the shard count follows
    the machine and is not one), and the GRCh38 FASTA's md5. The digests are cached in
    `identities` (`md5_cached`). `genokeeper calibrate` writes these into the table it
    measures; `genokeeper run` computes its own and compares them (`open_calibration`)."""
    fasta, checkpoint = calling_inputs(lab)
    return {
        "bwa": tools["bwa"],
        "samtools": tools["samtools"],
        "deepvariant_image": IMAGES["deepvariant"],
        "checkpoint_md5": ",".join(
            f"{name}:{md5_cached(checkpoint.parent / name, identities)}"
            for name in CHECKPOINT_FILES
        ),
        "make_examples_flags": " ".join(CALLING_FLAGS),
        "reference_md5": md5_cached(fasta, identities),
    }


@dataclass(frozen=True)
class Calibration:
    """The calibration table a run grades with, HG003's benchmark regions it was measured
    in, HG002's cells where HG003's truth does not reach (`hg002_cells`, empty when they
    are not current), and what the run records of them."""

    table: CalibrationTable
    benchmark: StrataIndex
    record: dict[str, object]
    extra: tuple[ExtraCells, ...] = ()


def open_calibration(lab: Path, versions: Mapping[str, str]) -> Calibration | None:
    """The lab's calibration table (`tables.folder(lab, "small_variants")`) and HG003's
    benchmark BED, or None with no table there (the run refuses after its QC step, which
    records the bases the calibration needs); refuse (`RunRefused`) a table measured under
    other `versions` than the run's own (`calibration_versions`) with `tables.open_table`'s
    reason, naming each that differs, a table that is unreadable, and a missing benchmark
    BED. HG002's cells are opened beside it (`open_cells`) and never refuse:
    `record["hg002_cells"]` holds each part's path and the sha256 of its manifest and of
    its `regions.bed`, or the reason they are not used.

    The version check's anchor comes from the versions the run computes itself (the tools
    on PATH, the digest and flags in the code, the md5 of the checkpoint and FASTA files on
    disk), which a table measured under another pipeline cannot move, because the run never
    reads its own versions from the table. No one gains by fooling it: it stops an honest
    mistake, a bwa, image, checkpoint or reference changed after calibrating, which would
    grade calls with rates measured on another pipeline (track design, "2b. Calibration":
    "A run with any of these different refuses the table")."""
    found = open_table(lab, "small_variants", versions=versions)
    folder = table_dir(lab)
    if isinstance(found, Stale):
        if found.reason == NOT_CALIBRATED:
            return None
        raise RunRefused(f"the calibration table at {folder} is {found.reason}")
    bench = layout(lab)["truth"] / BENCH_BED
    if not bench.is_file():
        raise RunRefused(
            f"missing input: no HG003 benchmark BED at {bench};"
            " run scripts/bootstrap.sh --with-calibration"
        )
    try:
        table = CalibrationTable.read(folder)
    except FileNotFoundError as e:
        raise RunRefused(
            f"missing input: no calibration table at {folder} ({_one_line(e)});"
            " run `genokeeper calibrate`"
        ) from e
    except (ValueError, KeyError) as e:
        raise RunRefused(
            f"unreadable calibration table at {folder} ({_one_line(e)});"
            " run `genokeeper calibrate` again"
        ) from e
    extra, reason = open_cells(lab, versions)
    cells: dict[str, object] = (
        {"reason": reason}
        if reason is not None
        else {
            part: {
                "path": str(cells_dir(lab) / part),
                "manifest_sha256": sha256_file(cells_dir(lab) / part / MANIFEST),
                "regions_sha256": sha256_file(cells_dir(lab) / part / REGIONS),
            }
            for part in (HAPLOID, ONLY)
        }
    )
    return Calibration(
        table=table,
        benchmark=StrataIndex.from_beds({"benchmark": bench}),
        record={
            "path": str(folder),
            "manifest_sha256": sha256_file(folder / MANIFEST),
            "benchmark": _file(bench, sha256_file(bench)),
            "hg002_cells": cells,
        },
        extra=extra,
    )


def align_step(
    work: Path,
    *,
    fastq1: Path,
    fastq2: Path,
    fasta: Path,
    shas: Mapping[str, str],
    versions: Mapping[str, str],
    threads: int,
) -> _Step:
    """The run's alignment (`align`) of the pair, cached under `<work>/align/<key>/`, the
    key holding the sha256 of both FASTQ and the FASTA (`shas` "fastq1", "fastq2",
    "fasta"), bwa's and samtools' `versions` and the read group's sample name. The run
    and the calibration sample (`genokeeper.calibrate.sample`) both align through here."""

    def build(out: Path) -> None:
        align(
            fastq1,
            fastq2,
            fasta,
            out / f"{GENOME_SAMPLE}.cram",
            sample=GENOME_SAMPLE,
            threads=threads,
        )

    return _cached(
        work,
        "align",
        ALIGN_STEP_VERSION,
        [shas["fastq1"], shas["fastq2"], shas["fasta"], versions["bwa"]]
        + [versions["samtools"], GENOME_SAMPLE],
        f"{GENOME_SAMPLE}.cram",
        build,
    )


def call_step(
    work: Path,
    *,
    aligned: _Step,
    fasta: Path,
    fasta_sha: str,
    checkpoint: Path,
    identities: Path,
    sex: Sex,
    grading: Grading,
    threads: int,
) -> _Step:
    """The run's DeepVariant call (`run_deepvariant`) of the `aligned` CRAM into a VCF and a
    gVCF beside it, cached under `<work>/call/<key>/`, the key holding the alignment's key,
    the image digest, the FASTA's sha256, the checkpoint's identity, the sex and, for a
    male, [par]. The run and the calibration sample both call through here, so HG003 is
    called with the arguments a run gives a sample of its sex. Raises `DockerUnavailable`
    when Docker is not running."""
    male = sex == "male"

    def build(out: Path) -> None:
        run_deepvariant(
            aligned.path,
            fasta,
            out,
            checkpoint=checkpoint,
            male=male,
            grading=grading,
            threads=threads,
        )

    return _cached(
        work,
        "call",
        CALL_STEP_VERSION,
        [aligned.key, IMAGES["deepvariant"], fasta_sha]
        + [_checkpoint_id(checkpoint, identities), str(sex)]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
        f"{GENOME_SAMPLE}.g.vcf.gz",
        build,
    )


def qc_step(
    work: Path,
    *,
    aligned: _Step,
    fasta: Path,
    lab: Path,
    versions: Mapping[str, str],
    grading: Grading,
    threads: int,
) -> _Step:
    """The QC metrics of the `aligned` CRAM before calling (`run_qc`) as `qc.json`, cached
    under `<work>/qc/<key>/`, the key holding the alignment's key, mosdepth's and samtools'
    `versions`, VerifyBamID2's image, the QC autosomes BED and the depth MAPQ. The run and
    the calibration sample both measure through here; `qc_gate` judges the metrics.
    Refuses (`RunRefused`) when Docker is not running."""

    def build(out: Path) -> None:
        measured = run_qc(
            aligned.path, fasta, grading=grading, lab=lab, threads=threads
        )
        (out / "qc.json").write_text(json.dumps(measured))

    try:
        return _cached(
            work,
            "qc",
            QC_STEP_VERSION,
            [aligned.key, versions["mosdepth"], versions["samtools"]]
            + [IMAGES["verifybamid2"], sha256_file(giab38_dir(lab) / QC_AUTOSOMES)]
            + [str(grading.qc.depth_mapq)],
            "qc.json",
            build,
        )
    except DockerUnavailable as e:
        raise RunRefused(_one_line(e)) from e


def sex_step(
    work: Path,
    *,
    aligned: _Step,
    fasta: Path,
    lab: Path,
    versions: Mapping[str, str],
    grading: Grading,
    threads: int,
) -> _Step:
    """The sex the `aligned` CRAM shows (`infer_sex`) as `sex.json`, cached under
    `<work>/sex/<key>/`, the key holding the alignment's key, mosdepth's version, the sex
    regions BED, [sex] and `[dosage] tolerance` (`sex_from_ratios` reads it). The run
    (given no sex) and the calibration sample both infer through here. Refuses
    (`RunRefused`) a sex the depth does not decide."""

    def build(out: Path) -> None:
        found = infer_sex_from_cram(
            aligned.path, grading=grading, fasta=fasta, lab=lab, threads=threads
        )
        (out / "sex.json").write_text(json.dumps({"sex": found}))

    regions = sha256_file(giab38_dir(lab) / SEX_REGIONS)
    try:
        return _cached(
            work,
            "sex",
            SEX_STEP_VERSION,
            [aligned.key, versions["mosdepth"], regions]
            + [json.dumps(asdict(grading.sex), sort_keys=True)]
            + [str(grading.dosage.tolerance)],
            "sex.json",
            build,
        )
    except SexRefused as e:
        raise RunRefused(_one_line(e)) from e


def recorded_bases(reports: Path) -> int | None:
    """The total bases the newest QC of the user's run recorded (`reports/work/qc/`), or
    None before a run measured them: the target `genokeeper calibrate` downsamples HG003
    to. An aggregate count, read from the run's own cache."""
    qc = sorted((reports / "work" / "qc").glob("*/qc.json"), key=os.path.getmtime)
    bases = [json.loads(p.read_text()).get("total_bases") for p in qc]
    return next((int(b) for b in reversed(bases) if b is not None), None)


def call_qc_step(work: Path, *, called: _Step, bcftools: str) -> _Step:
    """Ti/Tv and het/hom of the `called` VCF (`vcf_stats`) as `call_qc.json`, cached under
    `<work>/call_qc/<key>/`; the run and the calibration sample both measure through
    here."""
    vcf = called.path.with_name(f"{GENOME_SAMPLE}.vcf.gz")

    def build(out: Path) -> None:
        (out / "call_qc.json").write_text(json.dumps(vcf_stats(vcf)))

    return _cached(
        work,
        "call_qc",
        CALL_QC_STEP_VERSION,
        [called.key, bcftools],
        "call_qc.json",
        build,
    )


def qc_gate(metrics: Mapping[str, float], grading: Grading) -> None:
    """Refuse (`RunRefused`) a sample whose metrics fall outside [qc], naming each."""
    if failed := judge(dict(metrics), grading.qc):
        raise RunRefused("QC gate: " + "; ".join(failed))


def _called_json(called: Called) -> dict[str, Any]:
    return asdict(called)


def _called(data: dict[str, Any]) -> Called:
    ad = data["ad"]
    return Called(
        gt=data["gt"],
        dosage=data["dosage"],
        row_type=RowType(data["row_type"]),
        gq=data["gq"],
        dp=data["dp"],
        ad=None if ad is None else tuple(ad),
        vaf=data["vaf"],
        warnings=tuple(data["warnings"]),
    )


def _as_allele(value: list[Any]) -> Allele:
    chrom, pos, ref, alt = value
    return str(chrom), int(pos), str(ref), str(alt)


def _claim_set(path: Path) -> ClaimSet:
    data = json.loads(path.read_text())
    return ClaimSet(
        claims=[Claim(**c) for c in data["claims"]],
        input_rows=data["input_rows"],
        dropped=data["dropped"],
    )


def genome_calls(
    gvcf: Path,
    claims: Sequence[Claim],
    *,
    fasta: Path,
    haploid: Callable[[str, int], bool],
    scratch: Path,
) -> list[Called | None]:
    """The call at each claim's allele from our gVCF, in `claims` order, read one
    chromosome at a time (`calls_at` keeps every record near its claims in memory), its
    temporary files under `scratch`. A claim on chrM gets None: DeepVariant's chrM calls
    are not used (`genokeeper.mito`)."""
    by_chrom: dict[str, list[Allele]] = {}
    for c in claims:
        if c.chrom != CHRM:
            by_chrom.setdefault(c.chrom, []).append((c.chrom, c.pos, c.ref, c.alt))
    found: dict[Allele, Called] = {}
    for alleles in by_chrom.values():
        found |= calls_at(gvcf, alleles, fasta=fasta, haploid=haploid, scratch=scratch)
    return [
        None if c.chrom == CHRM else found[(c.chrom, c.pos, c.ref, c.alt)]
        for c in claims
    ]


def _checkpoint_id(checkpoint: Path, identities: Path) -> str:
    return _key(
        *(
            f"{name}:{sha256_cached(checkpoint.parent / name, identities)}"
            for name in CHECKPOINT_FILES
        )
    )


# --- pharmacogenomics (slice 1b) ------------------------------------------------------------


def _newest_dated(folder: Path, names: Sequence[str]) -> Path | None:
    """The newest `<folder>/<YYYYMMDD>/` holding every one of `names`."""
    dated = sorted(
        (d for d in folder.glob("*") if _CLINVAR_DATE.fullmatch(d.name)), reverse=True
    )
    return next((d for d in dated if all((d / n).is_file() for n in names)), None)


def newest_pharmgkb(lab: Path) -> Path | None:
    """The ClinPGx annotations download of the newest `<lab>/pharmgkb/<YYYYMMDD>/` holding
    either download; in one folder the live one."""
    folders = [
        (folder, name)
        for name in PHARMGKB_DOWNLOADS
        if (folder := _newest_dated(lab / "pharmgkb", [name])) is not None
    ]
    if not folders:
        return None
    folder, name = max(folders, key=lambda found: found[0].name)
    return folder / name


@dataclass(frozen=True)
class PgxInputs:
    """The bootstrap's slice-1b files a run reads, each read before any step: PharmCAT's
    positions file and its alleles, the newest CPIC folder with CPIC's gene-drug levels and
    MT-RNR1's allele names, the newest ClinPGx clinical annotations download and what it
    holds, and T1K's HLA index."""

    positions: Path
    alleles: list[Allele]
    cpic: Path
    cpic_levels: dict[tuple[str, str], str]
    mtrnr1_alleles: list[str]
    pharmgkb: Path
    annotations: ClinicalAnnotations
    t1k_index: Path

    def record(self) -> dict[str, object]:
        """What `run.inputs` records of them: each file's path and sha256."""
        return {
            "positions": _file(self.positions, sha256_file(self.positions)),
            "cpic": {
                name: _file(self.cpic / name, sha256_file(self.cpic / name))
                for name in CPIC_FILES.values()
            },
            "pharmgkb": {
                **_file(self.pharmgkb, sha256_file(self.pharmgkb)),
                "created": self.annotations.created,
            },
            "t1k_index": _file(self.t1k_index, sha256_file(self.t1k_index)),
        }


def pgx_inputs(lab: Path) -> PgxInputs:
    """The lab's slice-1b inputs (`PgxInputs`); refuse (`RunRefused`) one that is missing,
    naming the bootstrap, or unreadable."""
    positions = positions_path(lab)
    if not positions.is_file():
        raise _bootstrap_missing(
            f"no PharmCAT {PHARMCAT_VERSION} positions file {positions}"
        )
    cpic = _newest_dated(lab / "cpic", list(CPIC_FILES.values()))
    if cpic is None:
        raise _bootstrap_missing(
            f"no CPIC data in {lab}/cpic/<YYYYMMDD>/ ({', '.join(CPIC_FILES.values())})"
        )
    download = newest_pharmgkb(lab)
    if download is None:
        raise _bootstrap_missing(
            f"no ClinPGx clinical annotations {lab}/pharmgkb/<YYYYMMDD>/"
            f"{' or '.join(PHARMGKB_DOWNLOADS)}"
        )
    index = t1k_index(lab)
    if not index.is_file():
        raise _bootstrap_missing(f"no T1K HLA index {index}")
    try:
        found = PgxInputs(
            positions=positions,
            alleles=pharmcat_positions(positions),
            cpic=cpic,
            cpic_levels=load_cpic_levels(cpic / CPIC_FILES["pairs"]),
            mtrnr1_alleles=load_cpic_alleles(cpic / CPIC_FILES["MT-RNR1"]),
            pharmgkb=download,
            annotations=read_clinical_annotations(download),
            t1k_index=index,
        )
    except (ValueError, OSError, KeyError) as e:
        raise RunRefused(f"unreadable slice-1b input: {_one_line(e)}") from e
    return found


def _pgx_step(
    work: Path,
    step: str,
    version: str,
    parts: Sequence[str],
    name: str,
    build: Callable[[Path], None],
) -> _Step:
    """`_cached` for PharmCAT's step: Docker not running or an image missing, or PharmCAT
    failing or writing no report, refuses the run (`RunRefused`); the finished steps stay
    cached."""
    try:
        return _cached(work, step, version, parts, name, build)
    except (DockerUnavailable, ToolFailed, PharmcatFailed) as e:
        raise RunRefused(_one_line(e)) from e


def _failure(error: BaseException) -> str:
    """A caller's failure in one line: its exit status and the last line of its stderr,
    or the error's message."""
    if not isinstance(error, ToolFailed):
        return _one_line(error)
    said = [line.strip() for line in (error.stderr or "").splitlines() if line.strip()]
    return f"exit {error.returncode}" + "".join(f": {line}" for line in said[-1:])


# The MHC reads' step, as a failure's name: OptiType and T1K read them.
MHC_READS = "MHC reads"


def outside_steps(
    work: Path,
    *,
    aligned: _Step,
    fasta: Path,
    fasta_sha: str,
    samtools: str,
    index: Path,
    hla_only: bool = False,
) -> tuple[dict[str, _Step], dict[str, str]]:
    """The outside callers' cached steps on the `aligned` CRAM: BCyrius (`bcyrius`), Aldy
    on CYP2D6 and CYP2B6 (`aldy`), the MHC-region and unmapped reads (`hla_reads`, keyed by
    samtools' version and the region), and OptiType and T1K on those (`optitype`; `t1k`
    with the HLA `index`'s sha256); with `hla_only` the last three alone (slice 3f's
    `typing_steps`). Each caller's key holds its image.

    A step that fails (exit non-zero, no result, or one its reader cannot take) is left
    out, and its failure returned (caller -> why), so the genes it serves are unresolved and
    the run goes on; `_cached` keeps nothing of it, so a rerun runs it again. Docker not
    running or an image missing refuses the run (`RunRefused`)."""
    cram = aligned.path
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def step(
        caller: str,
        name: str,
        version: str,
        parts: list[str],
        result: str,
        build: Callable[[Path], None],
    ) -> None:
        try:
            steps[name] = _cached(work, name, version, parts, result, build)
        except DockerUnavailable as e:
            raise RunRefused(_one_line(e)) from e
        except (ToolFailed, CallerFailed, ValueError, KeyError, IndexError) as e:
            failed[caller] = _failure(e)

    def build_bcyrius(out: Path) -> None:
        run_bcyrius(cram, fasta, out)

    def build_aldy(out: Path) -> None:
        for gene in ("CYP2D6", "CYP2B6"):
            run_aldy(cram, fasta, out, gene)

    def build_hla_reads(out: Path) -> None:
        extract_hla_reads(cram, fasta, out)

    if not hla_only:
        step(
            BCYRIUS,
            "bcyrius",
            BCYRIUS_STEP_VERSION,
            [aligned.key, IMAGES["bcyrius"], fasta_sha],
            f"{BCYRIUS_ID}.tsv",
            build_bcyrius,
        )
        step(
            ALDY,
            "aldy",
            ALDY_STEP_VERSION,
            [aligned.key, IMAGES["aldy"], fasta_sha],
            "CYP2B6.aldy",
            build_aldy,
        )
    step(
        MHC_READS,
        "hla_reads",
        HLA_READS_STEP_VERSION,
        [aligned.key, samtools, fasta_sha, json.dumps(MHC_REGION)],
        FASTQ[1],
        build_hla_reads,
    )
    reads = steps.get("hla_reads")
    if reads is None:
        return steps, failed
    fq1, fq2 = reads.path.with_name(FASTQ[0]), reads.path

    def build_optitype(out: Path) -> None:
        run_optitype(fq1, fq2, out)

    def build_t1k(out: Path) -> None:
        run_t1k(fq1, fq2, out, index=index)

    step(
        OPTITYPE,
        "optitype",
        OPTITYPE_STEP_VERSION,
        [reads.key, IMAGES["optitype"]],
        f"{OPTITYPE_PREFIX}_result.tsv",
        build_optitype,
    )
    step(
        T1K,
        "t1k",
        T1K_STEP_VERSION,
        [reads.key, IMAGES["t1k"], sha256_file(index)],
        f"{T1K_PREFIX}_genotype.tsv",
        build_t1k,
    )
    return steps, failed


def _failed(gene: str, failed: Mapping[str, str], *callers: str) -> OutsideCall | None:
    """`gene` unresolved when one of the `callers` it needs failed, None when none did."""
    why = {c: failed[c] for c in callers if c in failed}
    if not why:
        return None
    reason = "; ".join(f"{c} failed ({w})" for c, w in why.items())
    return OutsideCall(
        gene, None, None, None, dict.fromkeys(why, "failed"), False, reason
    )


def outside_calls(
    steps: Mapping[str, _Step],
    failed: Mapping[str, str],
    pgx: PgxInputs,
    *,
    mito_calls: Sequence[MitoCall],
    coverage: Mapping[int, int],
    grading: Grading,
) -> list[OutsideCall]:
    """The outside calls PharmCAT is given (track design, slice 1b, "Outside calls"), one
    per gene, from the callers' cached results (`outside_steps`) and mutserve's: CYP2D6
    agreed by BCyrius and Aldy, CYP2B6 from Aldy, HLA-A and HLA-B where OptiType and T1K
    agree, MT-RNR1 from mutserve's levels. A gene whose caller `failed` is unresolved with
    the failure as its reason."""

    def aldy(gene: str) -> CallerCall:
        return read_aldy(steps["aldy"].path.with_name(f"{gene}.aldy"))

    hla = [
        g for gene in PGX_HLA if (g := _failed(gene, failed, MHC_READS, OPTITYPE, T1K))
    ]
    return [
        _failed("CYP2D6", failed, BCYRIUS, ALDY)
        or agree(read_bcyrius(steps["bcyrius"].path), aldy("CYP2D6"), "CYP2D6"),
        _failed("CYP2B6", failed, ALDY) or cyp2b6_outside(aldy("CYP2B6")),
        *(
            hla
            or hla_outside(
                read_optitype(steps["optitype"].path), read_t1k(steps["t1k"].path)
            )
        ),
        mtrnr1_outside(
            mito_calls,
            pgx.mtrnr1_alleles,
            min_homoplasmy=grading.mito.min_homoplasmy,
            coverage=coverage,
            min_dp=grading.min_dp,
        ),
    ]


def pharmcat_step(
    work: Path,
    *,
    calls: Mapping[Allele, Called],
    levels: Mapping[Allele, Level],
    pgx: PgxInputs,
    outside: Sequence[OutsideCall],
    male: bool,
    grading: Grading,
    phases: Mapping[Allele, tuple[int, int]] | None = None,
    trusted: frozenset[Allele] = frozenset(),
) -> _Step:
    """PharmCAT's run (`run_pharmcat`) on its input VCF, written from our `calls` at every
    position of `pgx.alleles` graded by `levels` (`write_pharmcat_vcf`, with the phase sets
    `phases` and the `trusted` alleles the slice 3f table earns), with the `outside` calls
    (`-po`), cached under `<work>/pharmcat/<key>/` with the HTML report beside the JSON one,
    for the user. The key holds the image, PharmCAT's version and the sha256 of both
    files: an input that changed, a grade the calibration table moved among them, rebuilds
    it."""
    with tempfile.TemporaryDirectory(dir=work, prefix=".pharmcat-input-") as tmp:
        vcf = write_pharmcat_vcf(
            calls,
            pgx.alleles,
            Path(tmp) / "pharmcat.vcf",
            male=male,
            levels=levels,
            grading=grading,
            phases=phases,
            trusted=trusted,
        )
        vcf_text = vcf.read_text(encoding="utf-8")
    lines = outside_lines(outside)

    def build(out: Path) -> None:
        given = out / "input"
        given.mkdir()
        (given / "pharmcat.vcf").write_text(vcf_text, encoding="utf-8")
        po = None
        if lines:
            po = given / "outside.tsv"
            po.write_text(lines, encoding="utf-8")
        run_pharmcat(given / "pharmcat.vcf", po, out)

    return _pgx_step(
        work,
        "pharmcat",
        PHARMCAT_STEP_VERSION,
        [IMAGES["pharmcat"], PHARMCAT_VERSION, _key(vcf_text), _key(lines)],
        REPORT_JSON,
        build,
    )


# --- polygenic scores (slice 1d) -------------------------------------------------------------

# The prs_genotypes step's artifact: one JSON line per target position, the allele and our
# call there (millions of them: never one JSON document).
PRS_GENOTYPES = "genotypes.jsonl.gz"
type Grade = Callable[[Allele, dict[str, Any] | None], tuple[Called, Reliability]]


@dataclass(frozen=True)
class PrsInputs:
    """The bootstrap's slice-1d files a run reads, each checked before any step: the chosen
    scores (`prs.scores`) with their scoring files, the HGDP+1kGP panel and its ancestry
    positions (Nextflow, pgsc_calc, the plugin and Java are checked, not read)."""

    scores: tuple[Score, ...]
    files: dict[str, Path]
    panel: Path
    ancestry: Path

    @functools.cached_property
    def shas(self) -> dict[str, str]:
        """Each scoring file's sha256 by its PGS ID, and the ancestry positions' as
        `ancestry`."""
        return {
            **{pgs_id: sha256_file(path) for pgs_id, path in self.files.items()},
            "ancestry": sha256_file(self.ancestry),
        }

    def record(self, panel_sha: str) -> dict[str, object]:
        """What `run.inputs` records of them: each file's path and sha256."""
        return {
            "scores": {i: _file(p, self.shas[i]) for i, p in self.files.items()},
            "panel": _file(self.panel, panel_sha),
            "ancestry_positions": _file(self.ancestry, self.shas["ancestry"]),
        }


def prs_inputs(lab: Path) -> PrsInputs:
    """The lab's slice-1d inputs (`PrsInputs`); refuse (`RunRefused`) one that is missing,
    naming the bootstrap."""
    scores = load_scores()
    files = {s.pgs_id: scoring_file(lab, s.pgs_id) for s in scores}
    main_nf = pgsc_calc_dir(lab) / "main.nf"
    java = java_home(lab) / "bin" / "java"
    for what, path in (
        *((f"PGS Catalog scoring file {p}", p) for p in files.values()),
        (f"HGDP+1kGP reference panel {panel_path(lab)}", panel_path(lab)),
        (f"ancestry positions {ancestry_pvar(lab)}", ancestry_pvar(lab)),
        (f"Nextflow {NEXTFLOW_VERSION} {nextflow_path(lab)}", nextflow_path(lab)),
        (f"pgsc_calc {PGSC_CALC_VERSION} {main_nf}", main_nf),
        (f"Java for Nextflow {java}", java),
    ):
        if not path.is_file():
            raise _bootstrap_missing(f"no {what}")
    plugin = nxf_home(lab) / "plugins" / f"nf-schema-{NF_SCHEMA_VERSION}"
    if not plugin.is_dir():
        raise _bootstrap_missing(f"no Nextflow plugin nf-schema {plugin}")
    return PrsInputs(scores, files, panel_path(lab), ancestry_pvar(lab))


def prs_genotypes_step(
    work: Path,
    *,
    gvcf: Path,
    call: _Step,
    prs: PrsInputs,
    fasta: Path,
    fasta_sha: str,
    haploid: Callable[[str, int], bool],
    bcftools: str,
) -> _Step:
    """Our calls at every target position (`target_positions`: each score's and each
    ancestry position), read one chromosome at a time (`calls_at` keeps every record near
    its alleles in memory), cached by the call step's key and the scoring files', ancestry
    positions' and FASTA's sha256."""

    def build(out: Path) -> None:
        positions = target_positions(
            [prs.files[s.pgs_id] for s in prs.scores], prs.ancestry, fasta=fasta
        )
        with gzip.open(out / PRS_GENOTYPES, "wt", encoding="utf-8") as handle:
            for _, alleles in itertools.groupby(positions, key=lambda a: a[0]):
                chromosome = list(alleles)
                # Its staging directory under reports/work: the deny rules cover it.
                found = calls_at(
                    gvcf, chromosome, fasta=fasta, haploid=haploid, scratch=out
                )
                for allele in chromosome:
                    line = json.dumps([*allele, _called_json(found[allele])])
                    handle.write(f"{line}\n")

    return _cached(
        work,
        "prs_genotypes",
        PRS_GENOTYPES_STEP_VERSION,
        [call.key, *prs.shas.values(), fasta_sha, bcftools],
        PRS_GENOTYPES,
        build,
    )


def _prs_sites(genotypes: Path, grade: Grade) -> Iterator[tuple[Allele, Called, Level]]:
    """Each target position's allele, our call there and its level, streamed."""
    with gzip.open(genotypes, "rt", encoding="utf-8") as handle:
        for line in handle:
            chrom, pos, ref, alt, data = json.loads(line)
            allele = (chrom, pos, ref, alt)
            called, reliability = grade(allele, data)
            yield allele, called, reliability.level


def pgsc_step(
    work: Path,
    *,
    genotypes: Path,
    grade: Grade,
    prs: PrsInputs,
    panel_sha: str,
    lab: Path,
) -> _Step:
    """pgsc_calc's run (`run_pgsc_calc`) on its target VCF, written from our calls at every
    target position (`genotypes`) graded by `grade` (`write_target_vcf`: an unknown or
    low call goes in as a no-call), cached under `<work>/pgsc/<key>/` with its results.
    The key holds the target VCF's, the scoring files' and the panel's sha256, pgsc_calc's
    commit, Nextflow's version and the images: a grade the calibration table moved
    rebuilds it. Docker not running, pgsc_calc failing, or a task
    the containment check cannot verify refuses the run (`RunRefused`), and nothing of the
    step is kept: a failed Nextflow's stdout and `.nextflow.log` go to
    `<work>/pgsc-failed/`, which the next failure replaces and the refusal names."""
    with tempfile.TemporaryDirectory(dir=work, prefix=".pgsc-input-") as tmp:
        target = write_target_vcf(
            _prs_sites(genotypes, grade), Path(tmp) / "target.vcf.gz"
        )

        def build(out: Path) -> None:
            given = out / "input"
            given.mkdir()
            shutil.copy(target, given / target.name)
            scores = [prs.files[s.pgs_id] for s in prs.scores]
            try:
                run_pgsc_calc(given / target.name, scores, prs.panel, out, lab=lab)
            except ToolFailed as e:
                kept = work / "pgsc-failed"
                shutil.rmtree(kept, ignore_errors=True)
                kept.mkdir()
                (kept / "nextflow.stdout").write_text(e.stdout or "")
                if (out / ".nextflow.log").is_file():
                    shutil.copy(out / ".nextflow.log", kept)
                stdout = kept / "nextflow.stdout"
                raise RunRefused(f"{_one_line(e)} Nextflow's output: {stdout}") from e

        parts = [
            sha256_file(target),
            *prs.shas.values(),
            panel_sha,
            PGSC_CALC_COMMIT,
            NEXTFLOW_VERSION,
            *(IMAGES[name] for name in PGSC_LABELS.values()),
        ]
        try:
            return _cached(work, "pgsc", PGSC_STEP_VERSION, parts, "results", build)
        except (DockerUnavailable, ToolFailed, ContainmentRefused) as e:
            raise RunRefused(_one_line(e)) from e


# --- research associations (slice 1d) ----------------------------------------------------

GWAS_ALLELES = "alleles.json.gz"
GENOTYPES = "genotypes.jsonl.gz"


def gwas_inputs(lab: Path) -> tuple[Path, Path]:
    """The newest GWAS Catalog associations download (`<lab>/gwas/<YYYYMMDD>/`) and the
    panel's variants its SNPs are located on (`panel_pvar`); refuse (`RunRefused`) one that
    is missing, naming the bootstrap."""
    folder = _newest_dated(lab / "gwas", [GWAS_DOWNLOAD])
    if folder is None:
        raise _bootstrap_missing(
            f"no GWAS Catalog associations {lab}/gwas/<YYYYMMDD>/{GWAS_DOWNLOAD}"
        )
    if not panel_pvar(lab).is_file():
        raise _bootstrap_missing(f"no HGDP+1kGP panel variants {panel_pvar(lab)}")
    return folder / GWAS_DOWNLOAD, panel_pvar(lab)


def gwas_alleles_step(
    work: Path, *, catalog: Path, pvar: Path, fasta: Path, shas: Sequence[str]
) -> _Step:
    """The catalog's single-SNP associations (`evidence.gwas.load_associations`), each
    located at its panel allele (`locate`), with the download's rows, where the rest went
    and its release; cached by the catalog's, the panel variants' and the FASTA's sha256
    (`shas`). A download its reader refuses refuses the run."""

    def build(out: Path) -> None:
        found = load_associations(catalog)
        panel = panel_alleles(pvar, {(a.chrom, a.pos) for a in found.associations})
        located, dropped = locate_gwas(found.associations, panel, fasta=fasta)
        document = {
            "release": found.release,
            "input_rows": found.input_rows,
            "dropped": {r: n + dropped[r] for r, n in found.dropped.items()},
            "located": [[asdict(a), list(allele)] for a, allele in located],
        }
        with gzip.open(out / GWAS_ALLELES, "wt", encoding="utf-8") as handle:
            json.dump(document, handle)

    try:
        return _cached(
            work, "gwas_alleles", GWAS_ALLELES_STEP_VERSION, shas, GWAS_ALLELES, build
        )
    except ValueError as e:
        raise RunRefused(_one_line(e)) from e


def _gwas_association(data: dict[str, Any]) -> Association:
    ci = data["ci"]
    return Association(
        **data | {"efo_ids": tuple(data["efo_ids"]), "ci": ci and (ci[0], ci[1])}
    )


def genotypes_step(
    work: Path,
    step: str,
    version: str,
    *,
    alleles: list[Allele],
    keyed_by: _Step,
    gvcf: Path,
    fasta: Path,
    haploid: Callable[[str, int], bool],
    parts: Sequence[str],
) -> _Step:
    """Our calls at every allele of `alleles` (sorted), read one chromosome at a time, as
    `step` (the research associations' located alleles, slice 1d; VEP's candidates, slice
    1e), cached by the step that chose them (`keyed_by`) and the call `parts`."""

    def build(out: Path) -> None:
        with gzip.open(out / GENOTYPES, "wt", encoding="utf-8") as handle:
            for _, group in itertools.groupby(alleles, key=lambda a: a[0]):
                chromosome = list(group)
                # Its staging directory under reports/work: the deny rules cover it.
                found = calls_at(
                    gvcf, chromosome, fasta=fasta, haploid=haploid, scratch=out
                )
                for allele in chromosome:
                    line = json.dumps([*allele, _called_json(found[allele])])
                    handle.write(f"{line}\n")

    return _cached(work, step, version, [keyed_by.key, *parts], GENOTYPES, build)


# --- variants ClinVar has not classified (slice 1e) ---------------------------------------

CLINVAR_VIDS = "vids.txt.gz"
ANNOTATED = "annotated.json.gz"


def predicted_inputs(lab: Path) -> tuple[VepFiles, Path]:
    """The bootstrap's slice-1e files (`annotate.vep.vep_files`) and the newest ClinGen
    dosage sensitivity list (`<lab>/clingen-dosage/<YYYYMMDD>/`); refuse (`RunRefused`)
    one that is missing, naming the bootstrap step that fetches it."""
    files = vep_files(lab)
    lacking = [str(p) for p in files.required() if not p.is_file()]
    if lacking:
        raise RunRefused(
            f"missing input: no VEP {VEP_CACHE_VERSION} file {lacking[0]} (and"
            f" {len(lacking) - 1} more); run scripts/bootstrap.sh --with-vep"
        )
    folder = _newest_dated(lab / "clingen-dosage", [DOSAGE_FILE])
    if folder is None:
        raise _bootstrap_missing(
            f"no ClinGen dosage sensitivity list {lab}/clingen-dosage/<YYYYMMDD>/"
            f"{DOSAGE_FILE}"
        )
    return files, folder / DOSAGE_FILE


def vep_step(
    work: Path,
    *,
    vcf: Path,
    fasta: Path,
    files: VepFiles,
    listed: Collection[str],
    threads: int,
    parts: Sequence[str],
) -> _Step:
    """VEP's records (`annotate.vep.annotate`, SpliceAI on the `listed` genes) for every
    PASS allele of our VCF (`vep_input`), with the vids VEP was given; cached by `parts`,
    which hold no ClinVar release: the run's grading filters the classified alleles.
    Docker not running, an image missing or VEP or SpliceAI failing refuses the run."""

    def build(out: Path) -> None:
        given = vep_input(vcf, out / VEP_INPUT, fasta=fasta)
        records = annotate(
            out,
            genes_for_spliceai=listed,
            files=files,
            fasta=fasta,
            threads=threads,
        )
        document = {"vids": given, "records": [asdict(r) for r in records]}
        with gzip.open(out / ANNOTATED, "wt", encoding="utf-8") as handle:
            json.dump(document, handle)

    try:
        return _cached(work, "vep", VEP_STEP_VERSION, parts, ANNOTATED, build)
    except (DockerUnavailable, ToolFailed) as e:
        raise RunRefused(f"VEP: {_one_line(e)}") from e


# --- stage 3 (slice 3a) ---------------------------------------------------------------------

EH_VCF = "eh.vcf"
FLANK_JSON = "flank.json"
# The `failed` key of RFC1's motif call, which fails apart from its step.
RFC1_FAILED = "rfc1"
# The locus whose GC-rich flanks tell a PCR-amplified library (spec, "Library type").
FLANK_GENE = "C9orf72"


@dataclass(frozen=True)
class RepeatResults:
    """What the `repeats` step gave on one CRAM: ExpansionHunter's calls by VARID, RFC1's
    motif call, REViewer's images (gene -> SVG path) and C9orf72's flank ratio (nan when
    the step failed). A failure is in `failed`: the step's under `expansionhunter` (its
    fields None), the motif call's under `rfc1` (the step kept)."""

    repeats: dict[str, RepeatCall] | None
    rfc1: Rfc1Call | None
    images: dict[str, str]
    flank: float
    failed: dict[str, str]


@dataclass(frozen=True)
class Stage3Results:
    """What the three stage-3 steps gave on one CRAM: ExpansionHunter's calls by VARID,
    RFC1's motif call, REViewer's images (gene -> SVG path), both SMN callers' results,
    Yleaf's haplogroup and its markers file, and C9orf72's flank ratio (nan when
    ExpansionHunter's step failed). A step that failed is in `failed` (tool -> one line)
    and its fields are None."""

    repeats: dict[str, RepeatCall] | None
    rfc1: Rfc1Call | None
    images: dict[str, str]
    caller: SmnCaller | None
    parascopy: ParascopyCall | None
    yleaf: YHaplogroup | None
    yleaf_out: Path | None
    flank: float
    failed: dict[str, str]


def repeat_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    sex: Sex,
    grading: Grading,
    mean_depth: float,
    threads: int,
) -> tuple[RepeatResults, dict[str, _Step]]:
    """The cached `repeats` step on `cram` (its align step's key `cram_key`):
    ExpansionHunter, REViewer on every locus with an allele above its normal range, and
    C9orf72's flank ratio over `mean_depth`; its key holds `cram_key`, both images, the
    catalog's and the FASTA's sha256, the sex, the mean depth and [repeat_loci]. Then RFC1's
    motif call from its JSON (`read_rfc1`; None for a no-call), whose failure (the analysis
    missing) is `failed["rfc1"]` and keeps the step. The run and the stage-3 calibration (all four trio samples) run through here.

    A tool that fails, or writes what its reader cannot take, leaves the step out and is in
    `failed["expansionhunter"]`; nothing of it is cached. Docker not running or an image
    missing refuses (`RunRefused`)."""
    loci = grading.repeats.loci
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def build_repeats(out: Path) -> None:
        vcf = run_expansionhunter(cram, fasta, out, sex=sex, threads=threads)
        calls = read_repeat_calls(vcf, out / "eh.json")
        above = [
            locus.locus_id
            for locus in loci.values()
            if locus.variant_id in calls
            and any(
                a.repeats > locus.normal_max for a in calls[locus.variant_id].alleles
            )
        ]
        if above:
            review(out, vcf=vcf, fasta=fasta, locus_ids=above)
        ratio = flank_ratio(
            cram, fasta, region=repeat_region(loci[FLANK_GENE]), mean_depth=mean_depth
        )
        (out / FLANK_JSON).write_text(json.dumps({"flank": ratio}))

    repeats: dict[str, RepeatCall] | None = None
    rfc1: Rfc1Call | None = None
    images: dict[str, str] = {}
    flank = math.nan
    found = _cached_step(
        work,
        steps,
        failed,
        "repeats",
        REPEATS_STEP_VERSION,
        [cram_key, catalog_sha256(), fasta_sha, str(sex)]
        + [IMAGES["expansionhunter"], IMAGES["reviewer"], str(mean_depth)]
        + [json.dumps({g: asdict(locus) for g, locus in loci.items()}, sort_keys=True)],
        EH_VCF,
        build_repeats,
        tool=REPEATS_TOOL,
    )
    if found is not None:
        folder = found.path.parent
        repeats = read_repeat_calls(found.path, folder / "eh.json")
        images = {
            gene: str(folder / f"rev.{locus.locus_id}.svg")
            for gene, locus in loci.items()
            if (folder / f"rev.{locus.locus_id}.svg").is_file()
        }
        flank = float(json.loads((folder / FLANK_JSON).read_text())["flank"])
        try:
            rfc1 = read_rfc1(folder / "eh.json", repeats)
        except ValueError as e:
            failed[RFC1_FAILED] = _failure(e)
    results = RepeatResults(
        repeats=repeats, rfc1=rfc1, images=images, flank=flank, failed=failed
    )
    return results, steps


def stage3_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    sex: Sex,
    mean_depth: float,
    lab: Path,
    grading: Grading,
    threads: int,
) -> tuple[Stage3Results, dict[str, _Step]]:
    """The three stage-3 steps on `cram` (its align step's key `cram_key`), each cached:
    `repeats` (`repeat_steps`), `smn` (SMNCopyNumberCaller and Parascopy's SMN1 model on the
    shared `parascopy_depth` step's depth, which slice 3d's paralog loci read too) and, for a
    male, `yleaf`. Each key holds `cram_key`, its images, the catalog's sha256, the FASTA's
    sha256 and the sex; `smn`'s also the Parascopy data md5s and the depth step's key. The
    run and the stage-3 calibration of HG002's two libraries both run through here.

    A tool that fails, or writes what its reader cannot take, leaves its step out and is in
    the results' `failed`; nothing of it is cached. Docker not running or an image missing
    refuses (`RunRefused`)."""
    found_repeats, steps = repeat_steps(
        work,
        cram=cram,
        cram_key=cram_key,
        fasta=fasta,
        fasta_sha=fasta_sha,
        sex=sex,
        grading=grading,
        mean_depth=mean_depth,
        threads=threads,
    )
    common = [cram_key, catalog_sha256(), fasta_sha, str(sex)]
    failed = dict(found_repeats.failed)

    step = partial(_cached_step, work, steps, failed)

    def build_smn(out: Path, depth: Path) -> None:
        read_smn_caller(run_smn_caller(cram, fasta, out, threads=threads))
        found = parascopy_cn(
            cram,
            fasta,
            depth,
            out,
            lab=lab,
            threads=threads,
            models=[Path(SMN_MODEL).stem],
            label=SMN_LABEL,
        )
        read_parascopy(found)

    def build_yleaf(out: Path) -> None:
        read_yleaf(run_yleaf(cram, fasta, out, threads=threads))

    caller: SmnCaller | None = None
    parascopy: ParascopyCall | None = None

    def make_depth() -> _Step:
        parascopy_data(lab)  # refused before the depth step runs, as 3a's run was
        return parascopy_depth_step(
            work,
            cram=cram,
            cram_key=cram_key,
            fasta=fasta,
            fasta_sha=fasta_sha,
            threads=threads,
        )

    depth = _attempt(steps, failed, PARASCOPY_DEPTH, make_depth, tool=SMN_TOOL)
    found = (
        None
        if depth is None
        else step(
            "smn",
            SMN_STEP_VERSION,
            [*common, IMAGES["smncaller"], IMAGES["parascopy"], parascopy_md5s(lab)]
            + [depth.key],
            "smn.json",
            lambda out: build_smn(out, depth.path),
            tool=SMN_TOOL,
        )
    )
    if found is not None:
        caller = read_smn_caller(found.path)
        parascopy = read_parascopy(found.path.parent / "cn" / "res.paralog.bed.gz")
    yleaf: YHaplogroup | None = None
    yleaf_out: Path | None = None
    if sex == "male":
        found = step(
            "yleaf",
            YLEAF_STEP_VERSION,
            [*common, IMAGES["yleaf"]],
            f"{YLEAF_FOLDER}/{YLEAF_PREDICTION}",
            build_yleaf,
            tool=YLEAF_TOOL,
        )
        if found is not None:
            yleaf = read_yleaf(found.path)
            out_file = found.path.parent / cram.stem / f"{cram.stem}.out"
            yleaf_out = out_file if out_file.is_file() else None
    results = Stage3Results(
        repeats=found_repeats.repeats,
        rfc1=found_repeats.rfc1,
        images=found_repeats.images,
        caller=caller,
        parascopy=parascopy,
        yleaf=yleaf,
        yleaf_out=yleaf_out,
        flank=found_repeats.flank,
        failed=failed,
    )
    return results, steps


# --- stage 3b: deletions ------------------------------------------------------------------

DELETIONS_UNION = "both/tp-base.vcf.gz"
SV_MISSING = "gnomAD-SV missing: run scripts/bootstrap.sh --with-sv"


@dataclass(frozen=True)
class DeletionResults:
    """The two callers' union (`read_union`), or None when a step failed; each failed step
    (`delly`, `dysgu`, `truvari`) -> one line."""

    called: list[Deletion] | None
    failed: dict[str, str]


# In a female sample's caller keys: its shards leave chrY out (`deletion_steps`).
NO_CHRY = "shards without chrY"


def deletion_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    male: bool,
) -> tuple[DeletionResults, dict[str, _Step]]:
    """The three deletion steps on `cram` (its align step's key `cram_key`), each cached:
    `delly` and `dysgu` at the same time (DELLY's chromosome shards one core each, at most
    `threads` at a time, and Dysgu's chromosome shards one after another on `threads`,
    clamped to this machine's CPUs; the core gate shares the cores), each output through
    `deletions_only` (each from `caller_min_length`: Dysgu at `min_length`, DELLY at half of
    it, so a DELLY call a little short can still confirm a Dysgu one), then
    `deletions_union`, Truvari between them. Each caller's key holds `cram_key`, its image,
    the FASTA's sha256 and its filter (`caller_filter`, as the deletion table's versions
    hold it); the union's both callers' keys, Truvari's image, its options and
    `min_length`. The run and the deletion calibration both run through here.

    The callers shard over the primary chromosomes (`PRIMARY`), chrY only for a `male`
    sample: on a female one's chrY DELLY finds too few reads to estimate its library and
    fails the step. A female sample's caller keys also hold `NO_CHRY`, as the sex may be
    given rather than read from the CRAM; a male one's keys are as before.

    A step that fails is in the results' `failed` and nothing after it runs; nothing of it
    is cached. Docker not running or an image missing refuses (`RunRefused`)."""
    t = grading.deletions
    dysgu_threads = max(1, min(threads, os.cpu_count() or 1))
    shards = [c for c in PRIMARY if male or c != "chrY"]
    runs = {
        "delly": (
            DELLY_STEP_VERSION,
            lambda out: run_delly(cram, fasta, out, threads=threads, shards=shards),
        ),
        "dysgu": (
            DYSGU_STEP_VERSION,
            lambda out: run_dysgu(
                cram, fasta, out, threads=dysgu_threads, shards=shards
            ),
        ),
    }
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def caller(tool: str) -> None:
        version, run = runs[tool]
        name = f"{tool}.deletions.vcf.gz"
        filters = caller_filter(tool, t.min_length)
        floor = caller_min_length(tool, t.min_length)

        def build(out: Path) -> None:
            deletions_only(run(out), out / name, caller=tool, min_length=floor)

        parts = [
            cram_key,
            IMAGES[tool],
            fasta_sha,
            filters,
            *([] if male else [NO_CHRY]),
        ]
        _cached_step(work, steps, failed, tool, version, parts, name, build)

    _side_by_side(*(partial(caller, tool) for tool in runs))
    if failed:
        return DeletionResults(None, failed), steps

    def build_union(out: Path) -> None:
        both = match_callers(
            steps["dysgu"].path, steps["delly"].path, out, min_length=t.min_length
        )
        read_union(both)

    union = _cached_step(
        work,
        steps,
        failed,
        "deletions_union",
        DELETIONS_UNION_STEP_VERSION,
        [steps["delly"].key, steps["dysgu"].key, IMAGES["truvari"]]
        + [" ".join(BETWEEN_CALLERS), str(t.min_length)],
        DELETIONS_UNION,
        build_union,
        tool="truvari",
    )
    if union is None:
        return DeletionResults(None, failed), steps
    called = read_union(union.path.parent)
    return DeletionResults(called, failed), steps


def _by_chrom(intervals: Iterable[Interval]) -> dict[str, list[Interval]]:
    found: dict[str, list[Interval]] = {}
    for interval in intervals:
        found.setdefault(interval[0], []).append(interval)
    return found


def load_screen(
    lab: Path,
    grading: Grading,
    called: Sequence[Deletion],
    work: Path,
    *,
    svtype: str = "DEL",
    callers: tuple[str, str] = ("dysgu", "delly"),
    thresholds: DeletionThresholds | None = None,
) -> Screen:
    """What the screen reads for `called` (of `svtype`, the union of `callers`): the listed
    genes' exons widened by `exon_margin`, the mask (GIAB's segdup and low-mappability BEDs
    of [strata], merged), the repeats (its low_complexity BED) and gnomAD-SV's PASS records
    of `svtype` with `AN` at least `min_gnomad_an` overlapping a called event, an insertion's
    within `max_distance` (`read_population`, its temporary file under `work`), under
    `thresholds` (default `[deletions]`, 3b's screen)."""
    t, giab = thresholds or grading.deletions, giab38_dir(lab)
    mask = merge(
        [
            *read_bed(giab / grading.strata["segdup"]),
            *read_bed(giab / grading.strata["low_mappability"]),
        ]
    )
    spans = [
        span for d in called for span in population_spans([d], pad=max_distance(d, t))
    ]
    population = read_population(
        sv_files(lab).gnomad, spans, t.min_gnomad_an, work, svtype=svtype
    )
    return Screen(
        targets=targets(load_exons(), t.exon_margin),
        mask=_by_chrom(mask),
        repeats=_by_chrom(read_bed(giab / grading.strata["low_complexity"])),
        population=population,
        thresholds=t,
        callers=callers,
    )


# --- slice 3c: per-gene coverage and sex-chromosome dosage -----------------------------------

COVERAGE_JSON = "coverage.json"
DOSAGE_JSON = "dosage.json"
_NOT_DELIVERED = "no delivered genotype file in data/raw"


@dataclass(frozen=True)
class FoundationResults:
    """Slice 3c's two measurements on the run's CRAM: each coverage gene's callability
    (`genokeeper.coverage`) and the sex-chromosome dosage (`genokeeper.sex`). A step that
    failed is in `failed` (step -> one line) and its field is None."""

    coverage: list[CoverageRow] | None
    dosage: Dosage | None
    failed: dict[str, str]


def _failed_step(failed: Mapping[str, str], step: str) -> dict[str, str]:
    """What the run records of a slice-3c step: its failure, or nothing."""
    return {"failed": failed[step]} if step in failed else {}


def _contigs(fasta: Path) -> dict[str, int]:
    """Contig -> length, from the FASTA's `.fai`."""
    lengths: dict[str, int] = {}
    for line in Path(f"{fasta}.fai").read_text().splitlines():
        name, length = line.split("\t")[:2]
        lengths[name] = int(length)
    return lengths


def foundation_steps(
    work: Path,
    *,
    aligned: _Step,
    called: _Step,
    fasta: Path,
    lab: Path,
    grading: Grading,
    benchmarks: Sequence[Path],
    owners: Mapping[str, str],
    male: bool,
    threads: int,
) -> tuple[FoundationResults, dict[str, _Step]]:
    """Slice 3c's two steps on the `aligned` CRAM, each cached. `coverage`: the coverage
    genes' MANE CDS +/- `[coverage] cds_margin` (`coverage.cds_bed`), one `qc.mosdepth`
    call per expected copy count (`expected_copies`) at `min_depth_per_copy` x copies, and
    each gene's row (`coverage_rows`): the shares inside the low-mappability and segdup
    strata and outside the `benchmarks` (read minus the lab's mask) beside its callable
    share, and its `owners` source. Its key holds the align key, the sex (and, for a male,
    [par]), the CDS resource's, those BEDs' and the mask's sha256, [coverage], the owners and
    mosdepth's version. `dosage`: the X and Y depth ratios (`sex.measure_ratios`) and the
    chrX heterozygous share from the `called` step's VCF (`sex.x_het_share`); its key holds
    the align and call keys, the sex regions BED's sha256, `[sex] mapq`, `[dosage] het_low,
    het_high, min_x_sites`, `[reliability] min_dp` and mosdepth's version.

    A step whose tool fails, or whose input file is missing, is in the results' `failed`;
    nothing of it is cached and the run goes on."""
    t = grading.coverage
    giab = giab38_dir(lab)
    cram = aligned.path
    version = _tool_version("mosdepth")
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}
    difficult = [giab / grading.strata[s] for s in ("low_mappability", "segdup")]

    def build_coverage(out: Path) -> None:
        genes = gene_coverage.cds_bed(
            gene_coverage.load_cds(),
            margin=t.cds_margin,
            contigs=_contigs(fasta),
            out=out / "cds.bed",
        )
        copies = gene_coverage.expected_copies(genes, male=male, par=grading.par)
        depth: dict[str, tuple[float, float]] = {}
        for n in sorted(set(copies.values())):
            group = sorted(
                (i, gene)
                for gene, ivs in genes.items()
                if copies[gene] == n
                for i in ivs
            )
            bed = out / f"cds.{n}.bed"
            bed.write_text("".join(f"{c}\t{a}\t{b}\t{g}\n" for (c, a, b), g in group))
            depth |= mosdepth(
                cram,
                fasta,
                bed,
                out / f"depth.{n}",
                mapq=t.mapq,
                threads=threads,
                thresholds=t.min_depth_per_copy * n,
            )
        benchmark = subtract(
            [i for b in benchmarks for i in read_bed(b)], read_bed(masked_bed(lab))
        )
        rows = gene_coverage.coverage_rows(
            depth,
            genes,
            difficult=[i for b in difficult for i in read_bed(b)],
            benchmark=benchmark,
            owners=owners,
            min_callable=t.min_callable,
        )
        (out / COVERAGE_JSON).write_text(json.dumps([asdict(r) for r in rows]))

    rows: list[CoverageRow] | None = None
    try:
        cds = importlib.resources.files("genokeeper").joinpath(
            "resources", "coverage", gene_coverage.RESOURCE
        )
        steps["coverage"] = _cached(
            work,
            "coverage",
            COVERAGE_STEP_VERSION,
            [aligned.key, "male" if male else "female"]
            + [json.dumps(grading.par if male else {}, sort_keys=True)]
            + [hashlib.sha256(cds.read_bytes()).hexdigest()]
            + [sha256_file(p) for p in [*difficult, *benchmarks, masked_bed(lab)]]
            + [
                json.dumps(asdict(t), sort_keys=True),
                json.dumps(owners, sort_keys=True),
            ]
            + [version],
            COVERAGE_JSON,
            build_coverage,
        )
    except (ToolFailed, FileNotFoundError) as e:
        failed["coverage"] = _failure(e)
    else:
        found = json.loads(steps["coverage"].path.read_text())
        rows = [CoverageRow(**row) for row in found]

    regions = giab / SEX_REGIONS
    vcf = called.path.with_name(f"{GENOME_SAMPLE}.vcf.gz")

    def build_dosage(out: Path) -> None:
        x, y = measure_ratios(
            cram, grading=grading, fasta=fasta, lab=lab, threads=threads
        )
        share, sites = x_het_share(
            vcf, xnonpar=read_named(regions, "xnonpar"), grading=grading
        )
        (out / DOSAGE_JSON).write_text(json.dumps(asdict(Dosage(x, y, share, sites))))

    dosage: Dosage | None = None
    d = grading.dosage
    try:
        steps["dosage"] = _cached(
            work,
            "dosage",
            DOSAGE_STEP_VERSION,
            [aligned.key, called.key, sha256_file(regions), str(grading.sex.mapq)]
            + [json.dumps([d.het_low, d.het_high, d.min_x_sites])]
            + [str(grading.min_dp), version],
            DOSAGE_JSON,
            build_dosage,
        )
    except (ToolFailed, FileNotFoundError, SexRefused) as e:
        failed["dosage"] = _failure(e)
    else:
        dosage = Dosage(**json.loads(steps["dosage"].path.read_text()))
    return FoundationResults(coverage=rows, dosage=dosage, failed=failed), steps


# --- slice 3d: copy number, insertions, paralog loci ------------------------------------

# The shared depth step's file (`parascopy_depth_step`), and the files the `dysgu_ins_dup`
# step keeps beside each other from 3b's Dysgu run.
PARASCOPY_DEPTH = "parascopy_depth"
DEPTH_CSV = "depth/depth.csv"
DYSGU_RAW = "dysgu.vcf.gz"
DYSGU_INS = "dysgu.ins.vcf.gz"
DYSGU_DUP = "dysgu.dup.vcf.gz"
DYSGU_RAW_GONE = "Dysgu's raw calls are gone: remove its step folder and rerun"
UNIQUE_JSON = "unique.json"
CNVPYTOR_DROPPED = (
    "CNVpytor dropped by the calibration: depth calls are DELLY's alone, so no large event"
    " is shown"
)
REGIONS_MISSING = (
    "large events not screened: ClinGen region lists missing; run scripts/bootstrap.sh"
)
SLICE3D_PARTS = ("cnv", "insertions", PARALOGS)


def _attempt(
    steps: dict[str, _Step],
    failed: dict[str, str],
    name: str,
    make: Callable[[], _Step],
    tool: str | None = None,
) -> _Step | None:
    """The step `make()` builds or reuses, recorded in `steps` under `name`; a failure of
    the deletion steps' class (a tool exiting non-zero or missing, a reader refusing its output, a
    missing or unwritable file, `_cached`'s "wrote no" `RuntimeError`) is recorded in
    `failed` (`tool`, default the step -> one line) instead and nothing of it is cached, so
    its part is unassessed and the run goes on. Docker not running or an image missing
    refuses the run (`RunRefused`), as every container step does. Every stage-3 part runs
    its steps through here."""
    try:
        steps[name] = make()
    except DockerUnavailable as e:
        raise RunRefused(_one_line(e)) from e
    except RunRefused:
        raise
    except (ToolFailed, ToolMissing, ValueError, KeyError, OSError, RuntimeError) as e:
        failed[tool or name] = _failure(e)
        return None
    return steps[name]


def _cached_step(
    work: Path,
    steps: dict[str, _Step],
    failed: dict[str, str],
    name: str,
    version: str,
    parts: Sequence[str],
    result: str,
    build: Callable[[Path], None],
    tool: str | None = None,
) -> _Step | None:
    """`_attempt` of the cached step `name` (`_cached`)."""
    return _attempt(
        steps,
        failed,
        name,
        lambda: _cached(work, name, version, parts, result, build),
        tool,
    )


def _side_by_side(*runs: Callable[[], object]) -> None:
    """`runs` at the same time (the core gate shares the cores); Docker not running
    refuses once all finished."""
    with ThreadPoolExecutor(max_workers=max(len(runs), 1)) as pool:
        futures = [pool.submit(run) for run in runs]
    for future in futures:
        future.result()


def parascopy_depth_step(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    threads: int,
) -> _Step:
    """Parascopy's background depth of `cram` (`smn.parascopy_depth`), once per CRAM: the
    `depth.csv` 3a's SMN step and 3d's `parascopy_loci` both read. Keyed by `cram_key`,
    Parascopy's image and the FASTA's sha256; raises what the tool raises."""
    return _cached(
        work,
        PARASCOPY_DEPTH,
        PARASCOPY_DEPTH_STEP_VERSION,
        [cram_key, IMAGES["parascopy"], fasta_sha],
        DEPTH_CSV,
        lambda out: _ignore(parascopy_depth(cram, fasta, out, threads=threads)),
    )


def _ignore(_: object) -> None:
    """A build step's result left on disk, where `_cached` looks for it."""


def dysgu_ins_dup_step(work: Path, *, dysgu: _Step, grading: Grading) -> _Step:
    """Dysgu's insertions and duplications out of 3b's `dysgu` step: its raw `dysgu.vcf.gz`
    (kept in that step's folder, `deletions.callers.run_dysgu`) kept to its PASS `INS`
    records of `[insertions] min_length` or more (`deletions_only`, `dysgu.ins.vcf.gz`) and
    its PASS `DUP` records of that length (`dysgu.dup.vcf.gz`, the depth gains' junctions)
    beside it. Keyed by 3b's step key and both filters; a raw file that is gone fails the
    step (`FileNotFoundError`, `DYSGU_RAW_GONE`)."""
    raw = dysgu.path.parent / DYSGU_RAW
    floor = grading.insertions.min_length

    def build(out: Path) -> None:
        if not raw.is_file():
            raise FileNotFoundError(f"{raw}: {DYSGU_RAW_GONE}")
        deletions_only(
            raw, out / DYSGU_DUP, caller="dysgu", min_length=floor, svtype="DUP"
        )
        deletions_only(
            raw, out / DYSGU_INS, caller="dysgu", min_length=floor, svtype="INS"
        )

    return _cached(
        work,
        "dysgu_ins_dup",
        DYSGU_INS_DUP_STEP_VERSION,
        [dysgu.key]
        + [filter_expression("dysgu", floor, svtype) for svtype in ("INS", "DUP")],
        DYSGU_INS,
        build,
    )


@dataclass(frozen=True)
class CnvResults:
    """The depth events (`merge_depth`), or None when a step failed; Dysgu's DUP records
    (the gains' junctions); each failed step -> one line."""

    events: list[Deletion] | None
    junctions: list[Deletion]
    failed: dict[str, str]


def cnv_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    male: bool,
    masked: Path,
    cnvpytor: bool,
    junctions: _Step,
) -> tuple[CnvResults, dict[str, _Step]]:
    """The depth callers' steps on `cram`, side by side as `deletion_steps` runs its
    callers: `dellycnv` (`delly cnv` from `[cnv] min_length`, sharded as 3b's DELLY, with
    3c's `masked` intervals excluded and, for a `male`, chrX and chrY haploid; a female's
    shards leave chrY out) and, when `cnvpytor` (the `cnv` table kept it), `cnvpytor` at
    `[cnv] bin_size`. Each key holds `cram_key`, its image, the FASTA's sha256 and its
    options (DELLY's also the masked BED's sha256 and the shards). The events are DELLY's
    read against the sample's ploidy and CNVpytor's merged (`merge_depth` at
    `MATCH_OVERLAP`); the junctions Dysgu's DUP records in `junctions`' `dysgu.dup.vcf.gz`.
    A step that fails is in `failed` and no events are given."""
    t = grading.cnv
    shards = [c for c in PRIMARY if male or c != "chrY"]
    haploid = ("chrX", "chrY") if male else ()
    ploidy = expected_ploidy(grading, male=male)
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def build_delly(out: Path) -> None:
        vcf = run_delly_cnv(
            cram,
            fasta,
            out,
            threads=threads,
            shards=shards,
            haploid=haploid,
            masked=masked,
            min_length=t.min_length,
        )
        read_delly_cnv(vcf, ploidy=ploidy)

    def build_cnvpytor(out: Path) -> None:
        calls = run_cnvpytor(cram, fasta, out, threads=threads, bin_size=t.bin_size)
        read_cnvpytor(calls, min_length=t.min_length)

    common = [cram_key, fasta_sha]
    runs = [
        lambda: _cached_step(
            work,
            steps,
            failed,
            "dellycnv",
            CNV_STEP_VERSION,
            [*common, IMAGES["delly"], str(t.min_length), sha256_file(masked)]
            + [",".join(shards), ",".join(haploid)],
            "dellycnv.vcf.gz",
            build_delly,
        )
    ]
    if cnvpytor:
        runs.append(
            lambda: _cached_step(
                work,
                steps,
                failed,
                "cnvpytor",
                CNVPYTOR_STEP_VERSION,
                [*common, IMAGES["cnvpytor"], str(t.bin_size)],
                "calls.tsv",
                build_cnvpytor,
            )
        )
    _side_by_side(*runs)
    if failed:
        return CnvResults(None, [], failed), steps
    try:
        delly = read_delly_cnv(steps["dellycnv"].path, ploidy=ploidy)
        pytor = (
            read_cnvpytor(steps["cnvpytor"].path, min_length=t.min_length)
            if cnvpytor
            else None
        )
        found = read_junctions(junctions.path.parent / DYSGU_DUP)
    except (ValueError, KeyError, OSError) as e:
        return CnvResults(None, [], {"cnv": _failure(e)}), steps
    events = merge_depth(delly, pytor, min_overlap=MATCH_OVERLAP)
    return CnvResults(events, found, {}), steps


@dataclass(frozen=True)
class InsertionResults:
    """The insertions called (the Dysgu+INSurVeyor union with Scramble's mobile elements,
    `add_scramble`), or None when a step failed; each failed step -> one line."""

    called: list[Deletion] | None
    failed: dict[str, str]


def insertion_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    dysgu_ins: _Step,
    regions: Path,
) -> tuple[InsertionResults, dict[str, _Step]]:
    """The insertion steps on `cram`: `insurveyor` and `scramble` (on `regions`, the screen
    genes' exons +/- `SCRAMBLE_MARGIN`, `scramble_regions`) side by side, each keyed by
    `cram_key`, its image and the FASTA's sha256 (Scramble's also the regions' sha256);
    then `insertions_union`, Truvari between Dysgu's kept insertions (`dysgu_ins`,
    `dysgu_ins_dup_step`) and INSurVeyor's at `[insertions] min_length`, keyed by both
    steps' keys, Truvari's image and options. The union labelled by Scramble's mobile
    elements within `[insertions] max_distance` (`add_scramble`). A step that fails is in
    `failed` and nothing is called."""
    t = grading.insertions
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def build_insurveyor(out: Path) -> None:
        run_insurveyor(cram, fasta, out, threads=threads)

    def build_scramble(out: Path) -> None:
        read_scramble(run_scramble(cram, fasta, out, regions=regions, threads=threads))

    common = [cram_key, fasta_sha]
    _side_by_side(
        lambda: _cached_step(
            work,
            steps,
            failed,
            "insurveyor",
            INSURVEYOR_STEP_VERSION,
            [*common, IMAGES["insurveyor"]],
            "insurveyor.vcf.gz",
            build_insurveyor,
        ),
        lambda: _cached_step(
            work,
            steps,
            failed,
            "scramble",
            SCRAMBLE_STEP_VERSION,
            [*common, IMAGES["scramble"], sha256_file(regions)],
            "meis.txt",
            build_scramble,
        ),
    )
    if failed:
        return InsertionResults(None, failed), steps
    insurveyor = steps["insurveyor"]

    def build_union(out: Path) -> None:
        insertion_union(dysgu_ins.path, insurveyor.path, out, min_length=t.min_length)

    union = _cached_step(
        work,
        steps,
        failed,
        "insertions_union",
        INSERTIONS_UNION_STEP_VERSION,
        [dysgu_ins.key, insurveyor.key, IMAGES["truvari"]]
        + [" ".join(BETWEEN_CALLERS), str(t.min_length)],
        DELETIONS_UNION,
        build_union,
    )
    if union is None:
        return InsertionResults(None, failed), steps
    try:
        called = read_union(union.path.parent, INSERTION_CALLERS, svtype="INS")
        meis = read_scramble(steps["scramble"].path)
    except (ValueError, KeyError, OSError) as e:
        return InsertionResults(None, {"insertions": _failure(e)}), steps
    return InsertionResults(
        add_scramble(called, meis, max_distance=t.max_distance), {}
    ), steps


@dataclass(frozen=True)
class ParalogResults:
    """What the paralog steps gave: Parascopy's call per locus (`read_locus`), CATSPER2's
    depth over the run's mean (STRC's unique part; None: not measured), Gauchian's GBA1 call,
    KILDA's KIV-2 count, Aldy's CYP2A6 copies; each failed step -> one line."""

    calls: dict[str, PsCall]
    unique_ratio: float | None
    gauchian: GauchianCall | None
    kiv2: float | None
    aldy_cyp2a6: int | None
    failed: dict[str, str]


PARALOG_STEPS = ("parascopy_loci", "gauchian", "kilda", "aldy_cyp2a6")


def paralog_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    lab: Path,
    threads: int,
    mean_depth: float,
    depth: _Step,
    run: Collection[str] = PARALOG_STEPS,
) -> tuple[ParalogResults, dict[str, _Step]]:
    """The paralog steps in `run` (all by default) on `cram`, each keyed by `cram_key`, its
    image and the FASTA's sha256: `parascopy_loci` (every model of `LOCI` in one
    `parascopy_cn` behind `dna:paralogs` on `depth`'s `depth.csv`, and CATSPER2's unique
    depth over `mean_depth`; also keyed by the depth step's key, the Parascopy data md5s,
    the loci and their spans and the mean depth), `gauchian`, `kilda` (on LPA's span) and
    `aldy_cyp2a6`. A step that fails, or is not in `run`, leaves its fields empty; a
    failure is in `failed`."""
    spans = load_spans()
    loci = {name: locus_span(spans, locus.genes) for name, locus in LOCI.items()}
    strc = LOCI["STRC"]
    unique = strc.unique or "CATSPER2"
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def build_loci(out: Path) -> None:
        bed = parascopy_cn(
            cram,
            fasta,
            depth.path,
            out,
            lab=lab,
            threads=threads,
            models=[locus.model for locus in LOCI.values()],
            label=PARALOGS_LABEL,
        )
        for name, at in LOCI.items():
            read_locus(bed, at, loci[name])
        _, models = parascopy_data(lab)
        parts = unique_part(spans[unique], models / f"{strc.model}.gz")
        ratio = depth_ratio(cram, fasta, parts, mean_depth)
        (out / UNIQUE_JSON).write_text(json.dumps({"ratio": ratio}))

    def build_gauchian(out: Path) -> None:
        tsv = run_gauchian(cram, fasta, out, threads=threads)
        read_gauchian(tsv, tsv.with_suffix(".json"))

    def build_kilda(out: Path) -> None:
        read_kilda(run_kilda(cram, fasta, out, threads=threads, span=spans["LPA"]))

    def build_aldy(out: Path) -> None:
        run_aldy_cyp2a6(cram, fasta, out)

    common = [cram_key, fasta_sha]

    def step(name: str, *args: Any) -> _Step | None:
        return _cached_step(work, steps, failed, name, *args) if name in run else None

    found = step(
        "parascopy_loci",
        PARASCOPY_LOCI_STEP_VERSION,
        [*common, IMAGES["parascopy"], depth.key, parascopy_md5s(lab), str(mean_depth)]
        + [json.dumps({n: [asdict(LOCI[n]), loci[n]] for n in LOCI}, sort_keys=True)],
        "cn/res.paralog.bed.gz",
        build_loci,
    )
    calls: dict[str, PsCall] = {}
    ratio: float | None = None
    if found is not None:
        calls = {n: read_locus(found.path, LOCI[n], loci[n]) for n in LOCI}
        ratio = json.loads((found.path.parent.parent / UNIQUE_JSON).read_text())[
            "ratio"
        ]
    found = step(
        "gauchian",
        GAUCHIAN_STEP_VERSION,
        [*common, IMAGES["gauchian"]],
        "gauchian.tsv",
        build_gauchian,
    )
    gauchian = (
        None
        if found is None
        else read_gauchian(found.path, found.path.with_suffix(".json"))
    )
    found = step(
        "kilda",
        KILDA_STEP_VERSION,
        [*common, IMAGES["kilda"], json.dumps(spans["LPA"])],
        "kilda-output/kilda_kiv2.CNs",
        build_kilda,
    )
    kiv2 = None if found is None else read_kilda(found.path)
    found = step(
        "aldy_cyp2a6",
        ALDY_CYP2A6_STEP_VERSION,
        [*common, IMAGES["aldy"]],
        "CYP2A6.aldy",
        build_aldy,
    )
    diplotype = None if found is None else read_aldy(found.path).diplotype
    copies = None if diplotype is None else aldy_copies(diplotype)
    results = ParalogResults(calls, ratio, gauchian, kiv2, copies, failed)
    return results, steps


@dataclass(frozen=True)
class Slice3dTables:
    """The three slice-3d calibration tables the run grades with, each None with why not
    (`reasons`, part -> reason)."""

    cnv: ClassTable | None
    insertions: ClassTable | None
    paralogs: ParalogTable | None
    reasons: dict[str, str]
    versions: dict[str, dict[str, str]]


def open_slice3d_tables(lab: Path, fasta_sha: str, grading: Grading) -> Slice3dTables:
    """The `cnv`, `insertions` and `paralogs` tables when each was measured under the
    versions the run computes itself (`cnv_versions`, `insertion_versions`,
    `paralog_versions`), through 3c's `tables.open_table`; otherwise None and why ("not
    calibrated: run `genokeeper calibrate`", a stale key, or an input missing). Opened before
    alignment, so a part without a current table never starts a step. 3b's deletion table
    is opened apart from these (C6)."""
    reasons: dict[str, str] = {}
    versions: dict[str, dict[str, str]] = {}
    found: dict[str, ClassTable | ParalogTable | None] = {}
    for name, compute in (
        ("cnv", cnv_versions),
        ("insertions", insertion_versions),
        (PARALOGS, paralog_versions),
    ):
        try:
            versions[name] = compute(lab, fasta_sha, grading)
        except FileNotFoundError as e:
            found[name], reasons[name] = None, f"{name} not run: {_one_line(e)}"
            continue
        if name == PARALOGS:
            table, reason = open_paralog_table(lab, versions[name])
        else:
            table, reason = open_class_table(lab, name, versions[name])
        found[name] = table
        if reason is not None:
            reasons[name] = reason
    return Slice3dTables(
        cnv=cast(ClassTable | None, found["cnv"]),
        insertions=cast(ClassTable | None, found["insertions"]),
        paralogs=cast(ParalogTable | None, found[PARALOGS]),
        reasons=reasons,
        versions=versions,
    )


@dataclass(frozen=True)
class ScreenedDeletions:
    """3b's screen as slice 3d reads it: the callers' union, the shown deletions, the screen
    they were shown under (None: 3b did not run) and each 3b class's chance of being
    false."""

    called: Sequence[Deletion]
    shown: Sequence[Shown]
    screen: Screen | None
    p_wrong: Callable[[CallClass], Rate | None]


@dataclass(frozen=True)
class Slice3dResults:
    """The four slice-3d sources, the 3b deletions read depth also supports ((chrom, pos,
    end) -> the depth callers, 3b's `DeletionSource.depth_support`), Gauchian's GBA1 call
    (for `gba1_unassessed`; None where the calibration withholds GBA1), the parts that ran and what the run records of the ClinGen
    region lists."""

    cnv: CnvSource
    insertions: DeletionSource
    paralogs: ParalogSource
    lpa: LpaSource
    joined: dict[tuple[str, int, int], str]
    gauchian: GauchianCall | None
    ran: frozenset[str]
    regions: dict[str, object]


def _failed_line(failed: Mapping[str, str]) -> str:
    return "; ".join(f"{step} failed: {why}" for step, why in failed.items())


def _clingen_regions(
    lab: Path,
) -> tuple[
    list[tuple[str, int, int, str]], list[tuple[str, int, int, str]], dict[str, object]
]:
    """ClinGen's region curations and recurrent CNV regions from the newest
    `<lab>/clingen-regions/<YYYYMMDD>/` holding both (bootstrap's `clingen_regions`), and
    what the run records of them; without them, no lists and the reason."""
    folder = _newest_dated(lab / "clingen-regions", [REGIONS_FILE, RECURRENT_FILE])
    if folder is None:
        return [], [], {"reason": REGIONS_MISSING}
    curations, recurrent = folder / REGIONS_FILE, folder / RECURRENT_FILE
    try:
        record: dict[str, object] = {
            "regions": {
                **_file(curations, sha256_file(curations)),
                "date": regions_date(curations),
            },
            "recurrent": _file(recurrent, sha256_file(recurrent)),
        }
        return load_regions(curations), load_recurrent(recurrent), record
    except (ValueError, UnicodeDecodeError, KeyError) as e:
        reason = f"large events not screened: ClinGen region lists unreadable ({e})"
        return [], [], {"reason": " ".join(reason.split())}


def slice3d_steps(
    work: Path,
    *,
    tables: Slice3dTables,
    lab: Path,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    male: bool,
    mean_depth: float,
    dysgu: _Step | None,
    deletions: ScreenedDeletions,
    dosage: DosageList,
) -> tuple[Slice3dResults, dict[str, _Step]]:
    """Slice 3d on `cram`, after 3b's screen, each part only with its table current (else
    unassessed with the table's reason, its steps never started):

    - copy number: `dysgu_ins_dup_step` (on 3b's `dysgu` step; without it the part is
      unassessed naming why 3b's Dysgu did not run) and `cnv_steps`, CNVpytor only when the
      table kept it; the events screened (`screen_cnv`) against 3b's shown deletions
      (`JOINED`), 3c's masked intervals, the regions a locus-specific caller owns, ClinGen's
      region lists (missing: the large bands withheld, `REGIONS_MISSING`) and the named
      loci, which 3b's rare deletions also reach (`named_deletions`);
    - insertions: `insertion_steps` on Dysgu's kept insertions, screened by 3b's
      `screen_all` with `[insertions]`;
    - paralogs and Lp(a): `parascopy_depth_step` (shared with 3a's SMN step) and
      `paralog_steps`; a failed step leaves the loci it serves unassessed naming it.

    Returns the four sources, the `joined` map 3b's source takes as `depth_support`, and
    the steps."""
    steps: dict[str, _Step] = {}
    reasons = dict(tables.reasons)
    dysgu_reason = "3b's Dysgu calls missing: the deletion screen did not run"
    ins_dup: _Step | None = None
    if dysgu is not None and (tables.cnv is not None or tables.insertions is not None):
        raw = dysgu
        failed: dict[str, str] = {}
        ins_dup = _attempt(
            steps,
            failed,
            "dysgu_ins_dup",
            lambda: dysgu_ins_dup_step(work, dysgu=raw, grading=grading),
        )
        if ins_dup is None:
            dysgu_reason = _failed_line(failed)
    ran: set[str] = set()

    # Copy number from read depth.
    cnv_table = tables.cnv
    regions, recurrent, regions_record = _clingen_regions(lab)
    cnvpytor = cnv_table is not None and (
        cnv_table.extra.get("cnvpytor", {}).get("decision") != "dropped"
    )
    screen = None
    called = 0
    segdups: dict[str, tuple[int, ...]] = {}
    named: list[tuple[NamedCnv, Shown]] = []
    withheld: dict[str, str] = (
        {}
        if cnv_table is None
        else dict(cnv_table.extra.get("large", {}).get("withheld", {}))
    )
    if cnv_table is not None and "cnv" not in reasons:
        if ins_dup is None:
            reasons["cnv"] = dysgu_reason
        else:
            masked = masked_bed(lab)
            found, cnv_ran = cnv_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                grading=grading,
                threads=threads,
                male=male,
                masked=masked,
                cnvpytor=cnvpytor,
                junctions=ins_dup,
            )
            steps |= cnv_ran
            if found.events is None:
                reasons["cnv"] = _failed_line(found.failed)
            else:
                try:
                    s = load_screen(
                        lab, grading, found.events, work, callers=DEPTH_CALLERS
                    )
                    # A gain is common by gnomAD-SV's duplications, a loss by its
                    # deletions (`is_common` matches the event's own type).
                    gains = read_population(
                        sv_files(lab).gnomad,
                        [(d.chrom, d.pos, d.end) for d in found.events],
                        grading.deletions.min_gnomad_an,
                        work,
                        svtype="DUP",
                    )
                    s = replace(s, population=[*s.population, *gains])
                    spans = load_spans()
                    _, model_dir = parascopy_data(lab)
                    owned = owned_regions(
                        spans,
                        [model_dir / f"{locus.model}.gz" for locus in LOCI.values()],
                    )
                    masked_intervals = read_bed(masked)
                except (ValueError, KeyError, OSError) as e:
                    reasons["cnv"] = f"cnv failed: {_failure(e)}"
                else:
                    loci = load_named()
                    screen = screen_cnv(
                        found.events,
                        [x.deletion for x in deletions.shown],
                        found.junctions,
                        s,
                        grading.cnv,
                        masked=masked_intervals,
                        owned=owned,
                        spans=gene_spans(s.targets),
                        dosage=dosage,
                        regions=regions,
                        recurrent=recurrent,
                        named=loci,
                    )
                    if deletions.screen is not None:
                        named = named_deletions(
                            deletions.called, loci, deletions.screen
                        )
                    called = len(found.events)
                    segdups = segdup_exons(
                        s.targets, [m for v in s.mask.values() for m in v], owned
                    )
                    ran.add("cnv")
    if "reason" in regions_record:
        withheld |= {
            band: str(regions_record["reason"]) for band in ("loss_large", "gain_large")
        }

    def depth_p_wrong(c: CallClass) -> Rate | None:
        # A named 3b deletion keeps its 3b class; a depth event's band is the cnv table's.
        if c.band in DELETION_BANDS:
            return deletions.p_wrong(c)
        return None if cnv_table is None else cnv_table.p_wrong(c)

    cnv_source = CnvSource(
        shown=[] if screen is None else screen.shown,
        large=[] if screen is None else screen.large,
        named=[] if screen is None else [*named, *screen.named],
        called=called,
        dropped={} if screen is None else screen.dropped,
        not_run=None if "cnv" in ran else reasons.get("cnv", NOT_CALIBRATED),
        p_wrong=depth_p_wrong,
        haploinsufficient=dosage.haploinsufficient,
        triplosensitive=dosage.triplosensitive,
        placement_overrides=screen_overrides(),
        segdup_exons=segdups,
        large_withheld=withheld,
        notes=() if cnvpytor or cnv_table is None else (CNVPYTOR_DROPPED,),
    )

    # Insertions.
    ins_table = tables.insertions
    shown: list[Shown] = []
    ins_dropped: dict[str, int] = {}
    ins_called = 0
    if ins_table is not None and "insertions" not in reasons:
        if ins_dup is None:
            reasons["insertions"] = dysgu_reason
        else:
            exons = targets(load_exons(), SCRAMBLE_MARGIN)
            bed = scramble_regions(exons, work / "scramble.regions.bed")
            found_ins, ins_ran = insertion_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                grading=grading,
                threads=threads,
                dysgu_ins=ins_dup,
                regions=bed,
            )
            steps |= ins_ran
            if found_ins.called is None:
                reasons["insertions"] = _failed_line(found_ins.failed)
            else:
                ins_called = len(found_ins.called)
                s = load_screen(
                    lab,
                    grading,
                    found_ins.called,
                    work,
                    svtype="INS",
                    callers=INSERTION_CALLERS,
                    thresholds=grading.insertions,
                )
                shown, ins_dropped = screen_all(found_ins.called, s)
                ran.add("insertions")
    insertion_source = DeletionSource(
        shown=shown,
        called=ins_called,
        dropped=ins_dropped,
        not_run=None
        if "insertions" in ran
        else reasons.get("insertions", NOT_CALIBRATED),
        p_wrong=(lambda _: None) if ins_table is None else ins_table.p_wrong,
        haploinsufficient=dosage.haploinsufficient,
        recall={} if ins_table is None else ins_table.recall(),
        wording=INSERTION_WORDING,
        placement_overrides=screen_overrides(),
        screen_claim={"not_covered": [MEI_SCOPE]},
    )

    # Paralog loci and Lp(a).
    paralog_table = tables.paralogs
    paralogs = ParalogResults({}, None, None, None, None, {})
    if paralog_table is not None and PARALOGS not in reasons:
        failed = {}
        depth = _attempt(
            steps,
            failed,
            PARASCOPY_DEPTH,
            lambda: parascopy_depth_step(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                threads=threads,
            ),
        )
        if depth is None:
            reasons[PARALOGS] = _failed_line(failed)
        else:
            paralogs, found_paralogs = paralog_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                lab=lab,
                threads=threads,
                mean_depth=mean_depth,
                depth=depth,
            )
            steps |= found_paralogs
            ran.add(PARALOGS)
    table_withheld = {} if paralog_table is None else paralog_table.withheld
    # A failed step leaves the loci it serves unassessed, naming it.
    serves = {
        "parascopy_loci": tuple(LOCI),
        "gauchian": ("GBA1",),
        "aldy_cyp2a6": ("CYP2A6",),
    }
    locus_withheld = {n: why for n, why in table_withheld.items() if n in LOCI}
    for step, why in paralogs.failed.items():
        for name in serves.get(step, ()):
            locus_withheld.setdefault(name, f"{step} failed: {why}")
    not_run = None if PARALOGS in ran else reasons.get(PARALOGS, NOT_CALIBRATED)
    paralog_source = ParalogSource(
        calls=paralogs.calls,
        unique_ratio=paralogs.unique_ratio,
        gauchian=paralogs.gauchian,
        aldy_cyp2a6=paralogs.aldy_cyp2a6,
        withheld=locus_withheld,
        not_run=not_run,
    )
    kilda = paralogs.failed.get("kilda")
    lpa_source = LpaSource(
        kiv2=paralogs.kiv2,
        withheld=table_withheld.get("LPA"),
        not_run=not_run if kilda is None else f"kilda failed: {kilda}",
    )
    results = Slice3dResults(
        cnv=cnv_source,
        insertions=insertion_source,
        paralogs=paralog_source,
        lpa=lpa_source,
        joined={} if screen is None else screen.joined,
        # A GBA1 the calibration withheld contradicts nothing (`gba1_unassessed`).
        gauchian=None if "GBA1" in table_withheld else paralogs.gauchian,
        ran=frozenset(ran),
        regions=regions_record,
    )
    return results, steps


def gba1_unassessed(
    collected: Sequence[Collected], gauchian: GauchianCall | None
) -> list[Collected]:
    """`collected` with every statement Gauchian contradicts in GBA1's exons 9-11 made
    unassessed (`GBA1_CONTRADICTED`; Decisions, "GBA1 outside the handled cases"): a
    ClinVar or predicted statement there that `gba1_contradicted` flags, on a call that
    carries its allele (Gauchian wild type agrees with a reference call). Its genotype is
    then unknown, as a mismapped pseudogene read leaves it: every row at its vid has its
    call's dosage unknown with a `GENOTYPE_UNKNOWN` warning naming why, and every finding or
    watch item there, and every GWAS association (its copies with it), is unassessed, so the
    store keeps an unassessed statement only beside an unknown call and no finding on one
    (`check_calls`); a non-finding (a benign classification) stays one. Gauchian not run:
    `collected` as it is."""
    if gauchian is None:
        return list(collected)
    exons = load_spans()[GBA1_HOMOLOGY[0]]
    calls: dict[object, Mapping[str, object]] = {}
    for c in collected:
        for row in c.rows:
            if row.call is not None:
                calls.setdefault(row.call["vid"], row.call)
    contradicted = {
        row.statement["vid"]
        for c in collected
        for row in c.rows
        if gba1_contradicted(row.statement, gauchian, exons)
        and calls.get(row.statement["vid"], {}).get("dosage")
        not in (None, "0", "unknown")
    }
    if not contradicted:
        return list(collected)
    warning = f"{GENOTYPE_UNKNOWN}: {GBA1_CONTRADICTED}"

    def judged(row: Row) -> Row:
        if row.statement.get("vid") not in contradicted:
            return row
        gwas = row.statement.get("source") == GWAS_TOOL
        call = row.call
        if call is not None:
            warnings = [str(w) for w in cast(list[object], call.get("warnings") or [])]
            call = {
                **call,
                "dosage": "unknown",
                "level": str(Level.LOW),
                "warnings": [warning, *warnings],
            }
        if row.statement.get("kind") == Kind.NON_FINDING and not gwas:
            return replace(row, call=call)
        statement = {
            **row.statement,
            "kind": str(Kind.UNASSESSED),
            "reason": GBA1_CONTRADICTED,
            "actionable": False,
        }
        if gwas:
            claim = json.loads(str(statement["claim"]))
            statement["claim"] = json.dumps({**claim, "risk_allele_copies": None})
        return replace(row, call=call, statement=statement)

    return [replace(c, rows=[judged(r) for r in c.rows]) for c in collected]


# --- slice 3f: the typing steps -----------------------------------------------------------

HHV6_STEP_VERSION = "1"
ROH_SITES_STEP_VERSION = "1"
ROH_STEP_VERSION = "1"
BLOODGROUPS_STEP_VERSION = "1"
ADMIXTURE_STEP_VERSION = "1"
ARCHAIC_STEP_VERSION = "1"
# The typing steps' names in `TypingResults.failed` ("<tool> failed: <line>").
HHV6_STEP, ROH_STEP, RBCEQ2_STEP = "HHV-6", "ROH", "RBCeq2"
ADMIXTURE_STEP, ARCHAIC_STEP = "fastmixture", "hmmix"
HHV6_COVERAGE = "coverage.tsv"
ROH_SITES, ROH_OUT = "sites.vcf.gz", "roh.txt"
RBCEQ2_SITES = "rbceq2_sites.vcf"
BLOODGROUPS_JSON, ADMIXTURE_JSON, ARCHAIC_JSON = (
    "bloodgroups.json",
    "admixture.json",
    "archaic.json",
)


def _discard(_result: object) -> None:
    """A step's build that keeps a tool's files and not what it returns."""


@dataclass(frozen=True)
class TypingResults:
    """What slice 3f's typing steps gave on one library: OptiType's and T1K's calls (T1K's
    every `HLA-` gene, `read_t1k(genes=None)`), the HHV-6 call, the runs of homozygosity,
    RBCeq2's calls with the systems no-called at a written position (R10) and the alleles
    left out of its input (R14), the admixture proportions and the archaic share. A step
    that failed is in `failed` (step -> one line) and its field is None; ROH without a
    population is not run and not failed."""

    optitype: dict[str, HlaCall] | None = None
    t1k: dict[str, HlaCall] | None = None
    hhv6: Hhv6Call | None = None
    roh: list[Run] | None = None
    bloodgroups: dict[str, SystemCall] | None = None
    bloodgroup_uncovered: dict[str, tuple[int, ...]] = field(default_factory=dict)
    bloodgroup_unchecked: dict[str, tuple[str, ...]] = field(default_factory=dict)
    admixture: dict[str, float] | None = None
    archaic: ArchaicShare | None = None
    failed: dict[str, str] = field(default_factory=dict)


def typing_steps(
    work: Path,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    vcf: Path,
    gvcf: Path,
    calls_key: str,
    mean_depth: float,
    population: str | None,
    male: bool,
    lab: Path,
    grading: Grading,
    threads: int,
    table: Slice3fTable | None,
    level: Callable[[Allele, Called], Level | None],
    tools: Mapping[str, str],
) -> tuple[TypingResults, dict[str, _Step]]:
    """Slice 3f's typing steps on one library (its CRAM `cram`, align key `cram_key`; its
    PASS calls `vcf` and indexed `gvcf`), each cached under `work`. The run and the slice 3f
    calibration both run through here.

    - `hla_reads`, `optitype` and `t1k`: `outside_steps`' own (`hla_only`), so a run
      reuses the MHC reads and both callers' results;
    - `hhv6`: the unmapped pairs aligned to HHV-6A and HHV-6B, the unique regions' depth,
      classified against half `mean_depth` and the table's background (none without one, so
      no ratio is a non-finding);
    - `roh_sites` and `roh` (only with a `population`, the superpopulation pgsc_calc names):
      our genotypes at that population's AF alleles graded by `level` (`roh.sites_vcf`),
      then `bcftools roh`;
    - `bloodgroups`: RBCeq2 on our calls at its positions (`bloodgroups.rbceq2_positions`),
      written by PharmCAT's writer with each call graded by `level`;
    - `admixture`: our calls at the panel's pruned SNPs, projected by fastmixture;
    - `archaic`: hmmix on `vcf`.

    `calls_key` stands for the calls and their grading (the run's: its call step and its
    calibration table; the calibration's: the call step and its own grade): every step that
    reads graded calls hashes it. Each key also holds the images and the inputs' sha256
    (`tools` gives samtools' and bwa's versions). A tool that fails, or writes what its reader
    cannot take, leaves its step out and is in `failed`; nothing of it is cached. Docker not
    running refuses (`RunRefused`)."""
    files = typing_files(lab)
    identities = work / IDENTITIES
    steps: dict[str, _Step] = {}
    failed: dict[str, str] = {}

    def sha(path: Path) -> str:
        return sha256_cached(path, identities)

    def haploid(chrom: str, pos: int) -> bool:
        return on_haploid_contig(chrom, pos, grading=grading, male=male)

    def graded(alleles: Sequence[Allele], scratch: Path) -> dict[Allele, Called]:
        scratch.mkdir(parents=True, exist_ok=True)
        try:
            return calls_at(
                gvcf, alleles, fasta=fasta, haploid=haploid, scratch=scratch
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    def step[T](
        tool: str,
        name: str,
        version: str,
        parts: Callable[[], list[str]],
        result: str,
        build: Callable[[Path], None],
        read: Callable[[Path], T],
    ) -> T | None:
        try:
            steps[name] = _cached(work, name, version, parts(), result, build)
            return read(steps[name].path)
        except DockerUnavailable as e:
            raise RunRefused(_one_line(e)) from e
        except RunRefused:
            raise
        except (
            ToolFailed,
            RuntimeError,
            ValueError,
            KeyError,
            IndexError,
            OSError,
        ) as e:
            steps.pop(name, None)
            failed[tool] = _failure(e)
            return None

    # HLA: outside_steps' own three steps (a run reuses its results), read here.
    found, hla_failed = outside_steps(
        work,
        aligned=_Step(cram_key, cram, False),
        fasta=fasta,
        fasta_sha=fasta_sha,
        samtools=tools["samtools"],
        index=t1k_index(lab),
        hla_only=True,
    )
    steps |= found
    failed |= hla_failed

    def typed(
        name: str, tool: str, read: Callable[[Path], dict[str, HlaCall]]
    ) -> dict[str, HlaCall] | None:
        try:
            return read(steps[name].path) if name in steps else None
        except (ValueError, KeyError, IndexError, OSError) as e:
            steps.pop(name)
            failed[tool] = _failure(e)
            return None

    optitype = typed("optitype", OPTITYPE, read_optitype)
    t1k = typed("t1k", T1K, lambda path: read_t1k(path, genes=None))

    def build_hhv6(out: Path) -> None:
        fq1, fq2 = extract_unmapped(cram, fasta, out / "reads")
        bam = align_viral(fq1, fq2, files.hhv6_fasta, out / "hhv6.bam", threads=threads)
        shutil.rmtree(out / "reads", ignore_errors=True)
        viral_coverage(bam, out / HHV6_COVERAGE)

    background = 0.0 if table is None else (table.hhv6_background or 0.0)
    hhv6 = step(
        HHV6_STEP,
        "hhv6",
        HHV6_STEP_VERSION,
        lambda: (
            [cram_key, fasta_sha, sha(files.hhv6_fasta)]
            + [tools["samtools"], tools["bwa"]]
        ),
        HHV6_COVERAGE,
        build_hhv6,
        lambda path: classify_hhv6(
            read_coverage(path),
            mean_depth / 2,
            grading=grading.hhv6,
            background=background,
        ),
    )

    roh: list[Run] | None = None
    if population is not None:
        af = files.af_file(population)

        def build_sites(out: Path) -> None:
            sites_vcf(
                gvcf,
                af,
                fasta,
                out / ROH_SITES,
                haploid=haploid,
                scratch=out / "scratch",
                level=level,
            )
            shutil.rmtree(out / "scratch", ignore_errors=True)

        sites = step(
            ROH_STEP,
            "roh_sites",
            ROH_SITES_STEP_VERSION,
            lambda: [calls_key, fasta_sha, sha(af), population],
            ROH_SITES,
            build_sites,
            lambda path: path,
        )
        if sites is not None:
            sites_key = steps["roh_sites"].key
            roh = step(
                ROH_STEP,
                "roh",
                ROH_STEP_VERSION,
                lambda: [sites_key, sha(af), IMAGES["truvari"]],
                ROH_OUT,
                lambda out: _discard(run_roh(sites, af, out / ROH_OUT)),
                read_runs,
            )

    def build_bloodgroups(out: Path) -> None:
        positions = rbceq2_positions(fasta)
        calls = graded(positions.alleles, out / "scratch")
        levels = {a: level(a, c) or Level.LOW for a, c in calls.items()}
        written = write_pharmcat_vcf(
            calls,
            positions.alleles,
            out / RBCEQ2_SITES,
            male=male,
            levels=levels,
            grading=grading,
        )
        found = read_rbceq2(run_rbceq2(written, out))
        document = {
            "calls": {s: asdict(c) for s, c in found.items()},
            "uncovered": uncovered_systems(written, positions),
            "unchecked": dict(positions.unchecked),
        }
        (out / BLOODGROUPS_JSON).write_text(json.dumps(document))

    blood = step(
        RBCEQ2_STEP,
        "bloodgroups",
        BLOODGROUPS_STEP_VERSION,
        lambda: (
            [calls_key, IMAGES["rbceq2"], fasta_sha, sha256_file(POSITIONS_FILE)]
            + [str(male), json.dumps(grading.par if male else {}, sort_keys=True)]
        ),
        BLOODGROUPS_JSON,
        build_bloodgroups,
        lambda path: json.loads(path.read_text()),
    )

    reference_p = files.reference_p()
    groups_file = files.panel_ref / "groups.json"
    bim = files.panel_ref / "pruned.bim"

    def build_admixture(out: Path) -> None:
        alleles = pruned_alleles(bim)
        calls = graded(alleles, out / "scratch")
        target = write_target_vcf(
            ((a, calls[a], level(a, calls[a])) for a in alleles), out / "target.vcf.gz"
        )
        bed = to_bed(target, out / "bed" / "target.bed")
        k = len(read_groups(groups_file))
        found = run_projection(bed, reference_p, k, out / "projection")
        (out / ADMIXTURE_JSON).write_text(json.dumps(found))

    admixture = step(
        ADMIXTURE_STEP,
        "admixture",
        ADMIXTURE_STEP_VERSION,
        lambda: (
            [calls_key, IMAGES["pgsc_plink2"], IMAGES["fastmixture"]]
            + [sha(bim), sha(reference_p), sha(groups_file)]
        ),
        ADMIXTURE_JSON,
        build_admixture,
        lambda path: {g: float(v) for g, v in json.loads(path.read_text()).items()},
    )

    hmmix_inputs = [p for p in files.required() if p.is_relative_to(files.hmmix)]

    def build_archaic(out: Path) -> None:
        share = archaic_share(read_segments(run_hmmix(vcf, files.hmmix, out / "hmmix")))
        shutil.rmtree(out / "hmmix" / "anc", ignore_errors=True)
        (out / ARCHAIC_JSON).write_text(json.dumps(asdict(share)))

    archaic = step(
        ARCHAIC_STEP,
        "archaic",
        ARCHAIC_STEP_VERSION,
        lambda: [sha(vcf), IMAGES["hmmix"], *(sha(p) for p in hmmix_inputs)],
        ARCHAIC_JSON,
        build_archaic,
        lambda path: ArchaicShare(**json.loads(path.read_text())),
    )
    results = TypingResults(
        optitype=optitype,
        t1k=t1k,
        hhv6=hhv6,
        roh=roh,
        bloodgroups=None
        if blood is None
        else {s: SystemCall(**c) for s, c in blood["calls"].items()},
        bloodgroup_uncovered={}
        if blood is None
        else {s: tuple(p) for s, p in blood["uncovered"].items()},
        bloodgroup_unchecked={}
        if blood is None
        else {s: tuple(n) for s, n in blood["unchecked"].items()},
        admixture=admixture,
        archaic=archaic,
        failed=failed,
    )
    return results, steps


# `run.inputs["gencc"]` without a GenCC download (the run grades on ClinGen's alone).
GENCC_MISSING = "no GenCC download: run scripts/bootstrap.sh"
# `run.inputs["pgx"]["phase"]` when the slice 3f table does not trust read-backed phase.
PHASE_UNTRUSTED = "not trusted: the slice 3f calibration did not earn read-backed phase"


def _grading_key(grading: Grading) -> str:
    """The grading as one deterministic string: a sha256 over its every value."""

    def plain(value: object) -> object:
        if isinstance(value, set | frozenset):
            return sorted(str(v) for v in cast(Iterable[object], value))
        return str(value)

    return _key(json.dumps(asdict(grading), sort_keys=True, default=plain))


@_single_writer
def run_genome(
    paths: Paths,
    *,
    fastq1: Path,
    fastq2: Path,
    clinvar38: Path,
    sex: Sex | None = None,
    grading: Grading | None = None,
    gene_validity: Path | None = None,
    faf: FafTable | None = None,
    cgd: Path | None = None,
    threads: int | None = None,
    delivered: Path | None = None,
) -> str:
    """Run stage 2a into `paths.reports/genomelab.duckdb` and return the run id: align the
    paired FASTQ to the lab's GRCh38 FASTA, gate the sample on QC, call it with DeepVariant
    into a gVCF, and grade every claim of the ClinVar GRCh38 release `clinvar38` on our own
    call at its allele: from the gVCF, or on chrM from mutserve's levels (`genokeeper.mito`).
    Candidate compound heterozygotes are phased from the reads (`genokeeper.phasing`) before
    the promotion judges them.

    The steps are cached (`_cached`), each keyed by its inputs: align (the FASTQ by their
    full sha256, kept with path, size, mtime and inode in `reports/work/identities.json`),
    qc, sex (when not given), call, call_qc, depth, claims, genotypes, mito, haplogroup
    (unless too much of chrM is uncovered), gwas_alleles and gwas_genotypes; slice 1b's
    pgx_alleles (PharmGKB's rsIDs in the ClinVar release), pgx_genotypes (our calls at them
    and at PharmCAT's positions), bcyrius, aldy, hla_reads, optitype, t1k and pharmcat
    (keyed by its input files' content); slice 1d's prs_genotypes (our calls at every
    score and ancestry position) and pgsc (keyed by its target VCF's content); slice 1e's
    clinvar_vids, vep (keyed by our calls, both images, the listed genes and the slice
    1e files' sha256; no ClinVar release, whose filter runs at grading) and
    predicted_genotypes; slice 3c's coverage and dosage (`foundation_steps`, a failure
    leaving its result out with the reason); slice 3b's delly, dysgu and deletions_union
    (`deletion_steps`,
    only with the sv step's files and a current deletion table, else the screen is
    unassessed with the reason, as it is when a step fails); slice 3d's parascopy_depth
    (shared with smn), dysgu_ins_dup, dellycnv, cnvpytor, insurveyor, scramble,
    insertions_union, parascopy_loci, gauchian, kilda and aldy_cyp2a6 (`slice3d_steps`, each
    part only with its table current, `open_slice3d_tables`, else unassessed with the
    reason); and phases (with candidate pairs, APOE's and two predicted carrier LoFs'
    among them). Then the gene results (`PgxSource`), the
    annotations (`PharmGkbSource`), the scores with APOE (`PrsSource`), the haplogroup
    (`MitoSource`), the research associations (`GwasSource`) and the predictions
    (`PredictedSource`) are stored beside ClinVar's claims, the deletions
    (`DeletionSource`, with the depth support slice 3d finds) after stage 3's, then slice
    3d's copy number, insertions, paralog loci and Lp(a); a ClinVar or predicted call in
    GBA1's exons 9-11 that Gauchian contradicts is unassessed (`gba1_unassessed`). The QC metrics are judged on every run against `grading`'s
    [qc], so a refusal keeps them cached and a rerun with another threshold measures
    nothing again.

    Slice 3f: haplocheck (`haplocheck`, after mito, whose VCF it reads) in every run, its
    flag demoting each chrM heteroplasmy; the DQ tag SNPs' and the traits' genotypes
    (`typing_genotypes`); GenCC's validity for the genes ClinGen does not curate (the
    newest `<lab>/gencc/<YYYYMMDD>/`, else ClinGen's alone, recorded why); the level-3
    annotations as research (a second `PharmGkbSource`); the chrM floor, the larger of
    `[mito] heteroplasmy_floor` and the slice 3f table's; and with a current slice 3f table
    (`calibrate.slice3f.open_slice3f`) PharmCAT's phase sets (`pgx_phases`) and UGT1A1*28 where the table
    trusts them, and `typing_steps` after the deletion screen (HLA, HHV-6, ROH, blood
    groups, ancestry). Without that table each statement it gates is unassessed with the
    reason; HLA, MITOMAP, HHV-6, ROH, the blood groups, ancestry and the traits are stored
    after the dosage statement.

    `grading=None` reads config/grading.toml, whose [par] and [strata] grade the run;
    `gene_validity` defaults to the newest ClinGen export, `cgd` to the newest CGD table,
    `faf` to the lab's FAF table, `threads` to this machine's CPUs. `delivered` is the
    delivered genotype file the sample's identity is checked against (None: "not
    checked", recorded with the reason); below `[identity] min_compared` SNPs compared it
    is not checked either, and below `min_concordance` the run refuses after the call
    step, before the claims step. Raises `RunRefused`,
    before any step, on a store another schema version wrote or another process holds, a
    missing tool or input (each naming the bootstrap; slice 1b's among them), a calibration
    table measured under other pipeline versions (naming `genokeeper calibrate`), or another
    run holding the lock;
    and then on a platform the T7 checkpoint was not trained on, a failed QC gate, no
    calibration table (after the QC step has recorded the bases), an
    ambiguous sex, Docker not running, a pharmacogenomic caller that wrote no result, or
    VEP or SpliceAI failing (the finished steps stay cached). Raises `CallMismatch`, rolling the run back, if
    a statement's kind contradicts its stored call (a bug).
    """
    grading_input: dict[str, object] = (
        {"given": True}
        if grading is not None
        else {
            "path": str(GRADING_FILE),
            "sha256": sha256_file(GRADING_FILE),
            "present": True,
        }
    )
    store = paths.reports / STORE_NAME
    _check_store(store)
    faf = faf if faf is not None else open_faf_table(paths.lab)
    ba1_exceptions = load_ba1_exceptions()
    grading = grading if grading is not None else load_grading()
    # Its anchor comes from the grading's BS1 values, which a table built under an older
    # config cannot move, because its manifest records the floor and tiers it was built on.
    if (gap := faf.uncovered(grading)) is not None:
        raise RunRefused(
            f"the gnomAD FAF table at {faf.folder} misses a BS1 of the grading ({gap});"
            " run genokeeper build-faf"
        )
    threads = threads if threads is not None else (os.cpu_count() or 4)
    versions = genome_tool_versions()
    bcftools = versions["bcftools"]
    lab = paths.lab
    _check_files({"FASTQ 1": fastq1, "FASTQ 2": fastq2, "ClinVar GRCh38": clinvar38})
    if delivered is not None:
        _check_files({"delivered genotypes": delivered})
    try:
        clonal = ClonalIndex.load(grading.af_guard.clonal_genes)
    except ValueError as e:
        raise RunRefused(_one_line(e)) from e
    fasta, checkpoint = calling_inputs(lab)
    giab = giab38_dir(lab)
    beds = {flag: giab / name for flag, name in grading.strata.items()}
    regions = [giab / QC_AUTOSOMES] + ([giab / SEX_REGIONS] if sex is None else [])
    lacking = sorted(str(p) for p in [*beds.values(), *regions] if not p.is_file())
    if lacking:
        raise _bootstrap_missing(f"GIAB GRCh38 BEDs {', '.join(lacking)}")
    validity_csv = (
        gene_validity if gene_validity is not None else newest_gene_validity(lab)
    )
    _check_files({"ClinGen gene validity": validity_csv})
    try:
        validity_date = gene_validity_date(validity_csv)
        curations = load_curations(validity_csv)
    except (ValueError, UnicodeDecodeError, csv.Error) as e:
        raise RunRefused(f"unreadable ClinGen gene validity {validity_csv}: {e}") from e
    # GenCC's submissions for the genes ClinGen does not curate (slice 3f); without a
    # download the run grades on ClinGen's alone and records why.
    gencc = newest_gencc(lab)
    validity_source: dict[str, str] = {}
    gencc_input: dict[str, object] = {"present": False, "reason": GENCC_MISSING}
    if gencc is not None:
        try:
            curations, validity_source = merge_gencc(curations, load_gencc(gencc))
            gencc_input = {
                **_file(gencc, sha256_file(gencc)),
                "date": gencc_date(gencc),
                "present": True,
            }
        except (ValueError, UnicodeDecodeError, csv.Error) as e:
            raise RunRefused(f"unreadable GenCC export {gencc}: {_one_line(e)}") from e
    validity_by_gene = {gene: gene_level(found) for gene, found in curations.items()}
    cgd = cgd if cgd is not None else newest_cgd(lab)
    if cgd is None or not cgd.is_file():
        raise _bootstrap_missing(
            f"no CGD table {cgd or f'{lab}/cgd/<YYYYMMDD>/CGD.txt.gz'}"
        )
    try:
        inheritance = load_cgd(cgd)
    except ValueError as e:
        raise RunRefused(_one_line(e)) from e
    pgx = pgx_inputs(lab)
    prs = prs_inputs(lab)
    catalog, variants = gwas_inputs(lab)
    annotation, dosage_list = predicted_inputs(lab)
    try:
        dosage = load_dosage(dosage_list)
    except (ValueError, UnicodeDecodeError) as e:
        raise RunRefused(f"unreadable ClinGen dosage list: {_one_line(e)}") from e
    release = _clinvar_release(clinvar38)
    lists = load_gene_lists()
    work = paths.reports / "work"
    identities = work / IDENTITIES
    shas = {
        "fastq1": sha256_cached(fastq1, identities),
        "fastq2": sha256_cached(fastq2, identities),
        "fasta": sha256_cached(fasta, identities),
        "clinvar38": sha256_file(clinvar38),
        "cgd": sha256_file(cgd),
        "pharmgkb": sha256_file(pgx.pharmgkb),
        "pgs_panel": sha256_cached(prs.panel, identities),
        "gwas": sha256_file(catalog),
        "panel_pvar": sha256_cached(variants, identities),
    }
    calibrated = calibration_versions(lab, tools=versions, identities=identities)
    calibration = open_calibration(lab, calibrated)
    # Without the Parascopy data stage 3 is unassessed, naming the bootstrap; the run goes on.
    try:
        stage3_calibrated = stage3_versions(lab, shas["fasta"], grading)
    except FileNotFoundError as e:
        stage3_table, stage3_reason = None, f"stage 3 not run: {_one_line(e)}"
    else:
        stage3_table, stage3_reason = open_stage3(lab, stage3_calibrated)
    # Without the sv step's files, a current deletion table or the files its versions hash,
    # the deletion screen is unassessed with the reason; the run goes on.
    deletion_table, deletion_reason = None, SV_MISSING
    if all(p.is_file() for p in sv_files(lab).required()):
        try:
            deletion_calibrated = deletion_versions(lab, shas["fasta"], grading)
        except FileNotFoundError as e:
            deletion_reason = f"deletions not run: {_one_line(e)}"
        else:
            deletion_table, deletion_reason = open_deletions(lab, deletion_calibrated)
    # Slice 3d's three tables, apart from 3b's (C6): a part without a current one is
    # unassessed with the reason and none of its steps starts.
    slice3d_tables = open_slice3d_tables(lab, shas["fasta"], grading)
    # Slice 3f's typing steps run only with a current slice 3f table; without one each
    # statement it gates is unassessed with the reason, and the run goes on.
    # Imported here: the calibration's modules import the pipeline.
    from genokeeper.calibrate.slice3f import CHRM_KEY, HHV6_KEY, ROH_KEY, open_slice3f

    slice3f_calibrated, slice3f_table, slice3f_reason = open_slice3f(lab, grading)
    withheld3f = {} if slice3f_table is None else slice3f_table.withheld
    steps: dict[str, _Step] = {}

    # Before hours of alignment: reads from another platform are refused here.
    try:
        platform = check_platform(fastq1)
    except PlatformRefused as e:
        raise RunRefused(_one_line(e)) from e

    steps["align"] = align_step(
        work,
        fastq1=fastq1,
        fastq2=fastq2,
        fasta=fasta,
        shas=shas,
        versions=versions,
        threads=threads,
    )
    cram = steps["align"].path
    steps["qc"] = qc_step(
        work,
        aligned=steps["align"],
        fasta=fasta,
        lab=lab,
        versions=versions,
        grading=grading,
        threads=threads,
    )
    metrics: dict[str, float] = json.loads(steps["qc"].path.read_text())
    qc_gate(metrics, grading)

    if calibration is None:
        raise RunRefused(
            f"this run has recorded its {metrics['total_bases']:.0f} bases, which"
            " `genokeeper calibrate` downsamples HG003 to; the calibration table at"
            f" {table_dir(lab)} is {NOT_CALIBRATED}"
        )

    if sex is None:
        steps["sex"] = sex_step(
            work,
            aligned=steps["align"],
            fasta=fasta,
            lab=lab,
            versions=versions,
            grading=grading,
            threads=threads,
        )
        sex = cast(Sex, json.loads(steps["sex"].path.read_text())["sex"])
        decided_by = "X copy count: chrX non-PAR depth over the autosomal (mosdepth)"
    else:
        decided_by = "given"
    male = sex == "male"

    def haploid(chrom: str, pos: int) -> bool:
        return on_haploid_contig(chrom, pos, grading=grading, male=male)

    try:
        steps["call"] = call_step(
            work,
            aligned=steps["align"],
            fasta=fasta,
            fasta_sha=shas["fasta"],
            checkpoint=checkpoint,
            identities=identities,
            sex=sex,
            grading=grading,
            threads=threads,
        )
    except DockerUnavailable as e:
        raise RunRefused(_one_line(e)) from e
    gvcf = steps["call"].path
    vcf = gvcf.with_name(f"{GENOME_SAMPLE}.vcf.gz")
    steps["call_qc"] = call_qc_step(work, called=steps["call"], bcftools=bcftools)
    metrics |= json.loads(steps["call_qc"].path.read_text())
    qc_gate(metrics, grading)

    def build_depth(out: Path) -> None:
        (out / "depth.json").write_text(json.dumps(depth_medians(vcf)))

    steps["depth"] = _cached(
        work,
        "depth",
        GENOME_DEPTH_STEP_VERSION,
        [steps["call"].key, bcftools],
        "depth.json",
        build_depth,
    )
    medians: dict[str, float] = json.loads(steps["depth"].path.read_text())

    # Research associations: the catalog's located at the panel's alleles, our calls there;
    # graded below, and first the sample's identity on them.
    steps["gwas_alleles"] = gwas_alleles_step(
        work,
        catalog=catalog,
        pvar=variants,
        fasta=fasta,
        shas=[shas["gwas"], shas["panel_pvar"], shas["fasta"]],
    )
    with gzip.open(steps["gwas_alleles"].path, "rt", encoding="utf-8") as handle:
        gwas = json.load(handle)
    gwas_located = [
        (_gwas_association(a), _as_allele(allele)) for a, allele in gwas["located"]
    ]
    gwas_alleles = sorted({allele for _, allele in gwas_located})
    steps["gwas_genotypes"] = genotypes_step(
        work,
        "gwas_genotypes",
        GWAS_GENOTYPES_STEP_VERSION,
        alleles=gwas_alleles,
        keyed_by=steps["gwas_alleles"],
        gvcf=gvcf,
        fasta=fasta,
        haploid=haploid,
        parts=[steps["call"].key, str(sex), bcftools]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
    )
    with gzip.open(steps["gwas_genotypes"].path, "rt", encoding="utf-8") as handle:
        gwas_data: dict[Allele, dict[str, Any]] = {
            (chrom, pos, ref, alt): data
            for chrom, pos, ref, alt, data in map(json.loads, handle)
        }

    # The sample's identity: our calls at the located SNPs against the delivered file; a
    # mismatch refuses before the claims are graded.
    identity: dict[str, object]
    if delivered is None:
        identity = {"verdict": str(Verdict.NOT_CHECKED), "reason": _NOT_DELIVERED}
        identity_text = f"not checked: {_NOT_DELIVERED}"
    else:
        ours = [
            (a.rsid, allele[2], allele[3], _called(gwas_data[allele]).dosage)
            for a, allele in gwas_located
        ]
        checked = compare(read_delivered(delivered), ours)
        verdict, text = judge_identity(checked, grading.identity)
        if verdict is Verdict.MISMATCH:
            raise RunRefused(text)
        matched = verdict is Verdict.MATCHED
        identity = {
            **_file(delivered, sha256_file(delivered)),
            "verdict": str(verdict),
            "compared": checked.compared,
            "concordant": checked.concordant,
            "skipped": checked.skipped,
            "reason": None if matched else text,
        }
        identity_text = text if matched else f"not checked: {text}"

    def build_claims_step(out: Path) -> None:
        found = build_claims(clinvar38, lists, fasta=fasta)
        (out / "claims.json").write_text(
            json.dumps(
                {
                    "input_rows": found.input_rows,
                    "dropped": found.dropped,
                    "claims": [asdict(c) for c in found.claims],
                }
            )
        )

    steps["claims"] = _cached(
        work,
        "claims",
        CLAIMS_STEP_VERSION,
        [shas["clinvar38"], shas["fasta"], _gene_lists_key(lists), bcftools],
        "claims.json",
        build_claims_step,
    )
    claims = _claim_set(steps["claims"].path)

    def build_genotypes(out: Path) -> None:
        # Its staging directory under reports/work: the deny rules cover it.
        found = genome_calls(
            gvcf, claims.claims, fasta=fasta, haploid=haploid, scratch=out
        )
        (out / "genotypes.json").write_text(
            json.dumps([None if c is None else _called_json(c) for c in found])
        )

    steps["genotypes"] = _cached(
        work,
        "genotypes",
        GENOTYPES_STEP_VERSION,
        [steps["claims"].key, steps["call"].key, str(sex), bcftools]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
        "genotypes.json",
        build_genotypes,
    )
    from_gvcf = json.loads(steps["genotypes"].path.read_text())

    def build_mito(out: Path) -> None:
        call_chrm(cram, fasta, out, threads=threads)

    try:
        steps["mito"] = _cached(
            work,
            "mito",
            MITO_STEP_VERSION,
            [steps["align"].key, IMAGES["mutserve"], shas["fasta"]],
            CALLS_FILE,
            build_mito,
        )
    except DockerUnavailable as e:
        raise RunRefused(_one_line(e)) from e
    mito_calls = read_calls(steps["mito"].path.parent)
    coverage = chrm_coverage(steps["mito"].path.parent)
    homoplasmy = grading.mito.min_homoplasmy

    # haplocheck on mutserve's VCF, every run: its flag only demotes (each heteroplasmy
    # unknown), so it acts with or without a slice 3f table. A failure, its image missing
    # among them (a lab bootstrapped before 3f), is recorded and demotes the same way
    # (rulings R20, R21); Docker not running refuses. Nothing of a failure is cached.
    def build_haplocheck(out: Path) -> None:
        run_haplocheck(steps["mito"].path.with_name(VCF_FILE), out)

    flagged: Haplocheck | None = None
    check_failed: str | None = None
    try:
        steps["haplocheck"] = _cached(
            work,
            "haplocheck",
            HAPLOCHECK_STEP_VERSION,
            [steps["mito"].key, IMAGES["haplocheck"]],
            f"{HAPLOCHECK_OUT}{HAPLOCHECK_RAW}",
            build_haplocheck,
        )
        flagged = read_haplocheck(steps["haplocheck"].path)
    except (ToolFailed, RuntimeError, ValueError, KeyError, OSError) as e:
        if isinstance(e, DockerUnavailable) and not image_missing("haplocheck"):
            raise RunRefused(_one_line(e)) from e
        steps.pop("haplocheck", None)
        check_failed = f"haplocheck failed: {_failure(e)}"
    haplocheck_input: dict[str, object] = {
        "status": "failed"
        if flagged is None
        else "contaminated"
        if flagged.contaminated
        else "not contaminated",
        "level": None if flagged is None else flagged.level,
        "major": None if flagged is None else flagged.major,
        "minor": None if flagged is None else flagged.minor,
        "reason": check_failed,
    }
    # The run's chrM floor: the config's, raised (never lowered) by the table's.
    floor = grading.mito.heteroplasmy_floor
    if slice3f_table is not None and slice3f_table.chrm_floor is not None:
        floor = max(floor, slice3f_table.chrm_floor)
    heteroplasmy_withheld = withheld3f.get(CHRM_KEY) or demoted(flagged, check_failed)

    # The haplogroup, unless too much of chrM is uncovered to state one.
    holes = uncovered(coverage, min_dp=grading.min_dp)
    haplogroup: Haplogroup | None = None
    if not too_uncovered(holes):

        def build_haplogroup(out: Path) -> None:
            vcf = write_chrm_vcf(
                mito_calls,
                coverage,
                out / VCF_NAME,
                min_homoplasmy=homoplasmy,
                min_dp=grading.min_dp,
            )
            run_haplogrep(vcf, out, tree=TREE)

        try:
            steps["haplogroup"] = _cached(
                work,
                "haplogroup",
                HAPLOGROUP_STEP_VERSION,
                [steps["mito"].key, IMAGES["haplogrep3"], TREE, str(homoplasmy)]
                + [str(grading.min_dp)],
                HAPLOGROUPS,
                build_haplogroup,
            )
            haplogroup = read_haplogroup(steps["haplogroup"].path)
        except (DockerUnavailable, ToolFailed, ValueError) as e:
            raise RunRefused(f"Haplogrep3: {_one_line(e)}") from e
    mito_source = MitoSource(haplogroup, holes, grading.min_dp)
    index = StrataIndex.from_beds(beds)

    def grade(
        allele: Allele, data: dict[str, Any] | None
    ) -> tuple[Called, Reliability]:
        """The call at `allele` from the gVCF's `data`, or on chrM (None) from
        mutserve's levels, and its grade."""
        chrom, pos, ref, alt = allele
        if data is None:  # chrM, from mutserve
            call = mito_called(
                allele,
                mito_calls,
                coverage,
                min_dp=grading.min_dp,
                min_homoplasmy=homoplasmy,
                floor=floor,
                heteroplasmy_withheld=heteroplasmy_withheld,
            )
            return call, grade_mito(
                call, chrom, pos, strata=index, min_homoplasmy=homoplasmy
            )
        call = _called(data)
        reliability = grade_called(
            call,
            chrom,
            pos,
            strata=index,
            medians=medians,
            grading=grading,
            indel=len(ref) != len(alt),
            haploid=haploid(chrom, pos),
            table=calibration.table,
            benchmark=calibration.benchmark,
            extra=calibration.extra,
        )
        return call, guard_allele_fraction(
            call,
            reliability,
            clonal=clonal.at(chrom, pos),
            haploid=haploid(chrom, pos),
            bound=grading.af_guard.min_het_vaf,
        )

    graded: list[tuple[Claim, Called, Reliability]] = []
    for claim, data in zip(claims.claims, from_gvcf, strict=True):
        call, reliability = grade((claim.chrom, claim.pos, claim.ref, claim.alt), data)
        graded.append((claim, call, reliability))
    source = ClaimsSource(claims, graded)

    def call_level(allele: Allele, called: Called) -> Level | None:
        """`grade`'s level of a call read elsewhere (the typing steps' inputs)."""
        return grade(allele, _called_json(called))[1].level

    # Slice 3f's genotypes, graded as every call: the six DQ tag SNPs and the traits' sites.
    typing_alleles = sorted({*TAG_ALLELES, *TRAIT_ALLELES})
    steps["typing_genotypes"] = genotypes_step(
        work,
        "typing_genotypes",
        TYPING_GENOTYPES_STEP_VERSION,
        alleles=typing_alleles,
        keyed_by=steps["call"],
        gvcf=gvcf,
        fasta=fasta,
        haploid=haploid,
        parts=[json.dumps(typing_alleles), str(sex), bcftools]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
    )
    with gzip.open(steps["typing_genotypes"].path, "rt", encoding="utf-8") as handle:
        typing_graded = {
            (chrom, pos, ref, alt): grade((chrom, pos, ref, alt), data)
            for chrom, pos, ref, alt, data in map(json.loads, handle)
        }
    # A tag SNP's dosage only from a graded call: an unknown or low one backs no exclusion.
    tag_dosages = {
        rsid: called.dosage
        if called.dosage != "unknown" and reliability.level not in (None, Level.LOW)
        else "unknown"
        for rsid, allele in zip(TAG_RSIDS, TAG_ALLELES, strict=True)
        for called, reliability in [typing_graded[allele]]
    }

    # Pharmacogenomics: PharmGKB's rsIDs located in the ClinVar release, our calls at
    # them and at every PharmCAT position, the outside callers, then PharmCAT.
    research = pgx.annotations.research  # the level-3 annotations, slice 3f's
    rsids = sorted({a.variant for a in [*pgx.annotations.annotations, *research]})

    def build_pgx_alleles(out: Path) -> None:
        found = clinvar_alleles(clinvar38, rsids, fasta=fasta)
        (out / "alleles.json").write_text(json.dumps(found))

    steps["pgx_alleles"] = _cached(
        work,
        "pgx_alleles",
        PGX_ALLELES_STEP_VERSION,
        [shas["clinvar38"], shas["pharmgkb"], shas["fasta"], bcftools],
        "alleles.json",
        build_pgx_alleles,
    )
    rsid_alleles = {
        rsid: [_as_allele(a) for a in found]
        for rsid, found in json.loads(steps["pgx_alleles"].path.read_text()).items()
    }
    # Every annotation's allele is called; which PharmCAT answers is known after its run.
    candidates, _ = locate([*pgx.annotations.annotations, *research], rsid_alleles)
    pgx_alleles = list(dict.fromkeys([*pgx.alleles, *(a for _, a in candidates)]))
    from_gvcf_pgx = [a for a in pgx_alleles if a[0] != CHRM]

    def build_pgx_genotypes(out: Path) -> None:
        # Its staging directory under reports/work: the deny rules cover it.
        found = calls_at(gvcf, from_gvcf_pgx, fasta=fasta, haploid=haploid, scratch=out)
        (out / "genotypes.json").write_text(
            json.dumps([[list(a), _called_json(found[a])] for a in from_gvcf_pgx])
        )

    steps["pgx_genotypes"] = _cached(
        work,
        "pgx_genotypes",
        PGX_GENOTYPES_STEP_VERSION,
        [steps["call"].key, json.dumps(from_gvcf_pgx), str(sex), bcftools]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
        "genotypes.json",
        build_pgx_genotypes,
    )
    pgx_data: dict[Allele, dict[str, Any]] = {
        _as_allele(a): data
        for a, data in json.loads(steps["pgx_genotypes"].path.read_text())
    }
    pgx_graded = {
        a: grade(a, None if a[0] == CHRM else pgx_data[a]) for a in pgx_alleles
    }
    callers, failed = outside_steps(
        work,
        aligned=steps["align"],
        fasta=fasta,
        fasta_sha=shas["fasta"],
        samtools=versions["samtools"],
        index=pgx.t1k_index,
    )
    steps |= callers
    outside = outside_calls(
        steps, failed, pgx, mito_calls=mito_calls, coverage=coverage, grading=grading
    )
    # Slice 3f: UGT1A1*28's low-graded call and read-backed phase sets, each only when the
    # slice 3f table earns it; otherwise PharmCAT runs as before 3f.
    trusted = (
        frozenset(UGT1A1_28)
        if slice3f_table is not None and slice3f_table.ugt1a1_trusted
        else frozenset[Allele]()
    )
    ugt1a1_unread = not trusted and not all(
        a in pgx_graded and pgx_graded[a][1].level not in (None, Level.LOW)
        for a in UGT1A1_28
    )
    pgx_phases: dict[Allele, tuple[int, int]] | None = None
    phase_reason = slice3f_reason
    if slice3f_table is not None and not slice3f_table.phase_trusted:
        phase_reason = PHASE_UNTRUSTED
    elif slice3f_table is not None:
        het = phase_candidates({a: pgx_graded[a][0] for a in pgx.alleles}, haploid)
        regions = [(c, start, end) for _, c, start, end in gene_regions(pgx.positions)]

        def build_pgx_phases(out: Path) -> None:
            found = phase_sites(cram, fasta, vcf, het, regions)
            (out / "phases.json").write_text(
                json.dumps([[list(a), *found[a]] for a in het if a in found])
            )

        try:
            steps["pgx_phases"] = _cached(
                work,
                "pgx_phases",
                PGX_PHASES_STEP_VERSION,
                [steps["align"].key, steps["call"].key, IMAGES["whatshap"]]
                + [shas["fasta"], json.dumps(het), json.dumps(regions)],
                "phases.json",
                build_pgx_phases,
            )
        except DockerUnavailable as e:
            raise RunRefused(_one_line(e)) from e
        except (ToolFailed, RuntimeError, ValueError, OSError) as e:
            phase_reason = f"WhatsHap failed: {_failure(e)}"
        else:
            phase_reason = None
            pgx_phases = {
                _as_allele(a): (ps, haplotype)
                for a, ps, haplotype in json.loads(steps["pgx_phases"].path.read_text())
            }
    steps["pharmcat"] = pharmcat_step(
        work,
        calls={a: pgx_graded[a][0] for a in pgx.alleles},
        levels={a: pgx_graded[a][1].level for a in pgx.alleles},
        pgx=pgx,
        outside=outside,
        male=male,
        grading=grading,
        phases=pgx_phases,
        trusted=trusted,
    )
    results = parse_report(steps["pharmcat"].path)
    left, answered = left_by_pharmcat(
        pgx.annotations.annotations, {r.gene for r in results if r.drugs}
    )
    located, unlocated = locate(left, rsid_alleles)
    unlocated["pharmcat_gene"] += answered
    pgx_source = PgxSource(
        results=results,
        outside={call.gene: call for call in outside},
        cpic_levels=pgx.cpic_levels,
        version=PHARMCAT_VERSION,
        ugt1a1_unread=ugt1a1_unread,
    )
    pharmgkb_source = PharmGkbSource(
        located=[(a, allele, *pgx_graded[allele]) for a, allele in located],
        input_rows=pgx.annotations.input_rows,
        dropped={
            reason: pgx.annotations.dropped[reason] + unlocated[reason]
            for reason in PHARMGKB_DROP_REASONS
        },
        version=pgx.annotations.created,
    )
    # The level-3 annotations as research (slice 3f), left and located as the first.
    research_left, research_answered = left_by_pharmcat(
        research, {r.gene for r in results if r.drugs}
    )
    research_located, research_dropped = locate(research_left, rsid_alleles)
    research_dropped["pharmcat_gene"] += research_answered
    research_source = PharmGkbSource(
        located=[(a, allele, *pgx_graded[allele]) for a, allele in research_located],
        input_rows=len(research),
        dropped=research_dropped,
        version=pgx.annotations.created,
        tool=PHARMGKB_RESEARCH_TOOL,
    )

    # Polygenic scores: our calls at every score and ancestry position, then pgsc_calc on
    # them, graded; APOE from the graded ClinVar claims at its two SNPs.
    steps["prs_genotypes"] = prs_genotypes_step(
        work,
        gvcf=gvcf,
        call=steps["call"],
        prs=prs,
        fasta=fasta,
        fasta_sha=shas["fasta"],
        haploid=haploid,
        bcftools=bcftools,
    )
    steps["pgsc"] = pgsc_step(
        work,
        genotypes=steps["prs_genotypes"].path,
        grade=grade,
        prs=prs,
        panel_sha=shas["pgs_panel"],
        lab=lab,
    )
    try:
        prs_results = parse_results(steps["pgsc"].path)
    except (PgscFailed, KeyError, ValueError) as e:
        raise RunRefused(f"unreadable pgsc_calc results: {_one_line(e)}") from e
    unreported = sorted(set(prs.files) - set(prs_results))
    if unreported:
        raise RunRefused(
            f"pgsc_calc reported nothing of {', '.join(unreported)}"
            f" ({steps['pgsc'].path})"
        )
    # The superpopulation pgsc_calc finds the user closest to: ROH's AF file, slice 3f's
    # similarity statement.
    similarity: PopSimilarity | None = None
    similarity_failed: str | None = None
    try:
        similarity = parse_similarity(steps["pgsc"].path)
    except (PgscFailed, KeyError, ValueError) as e:
        similarity_failed = f"pgsc_calc's population similarity: {_one_line(e)}"
    apoe = [
        g
        for g in graded
        if (g[0].chrom, g[0].pos, g[0].ref, g[0].alt) in (APOE_E4, APOE_E2)
    ]

    # Research associations, graded on our calls at their located alleles.
    gwas_graded = {allele: grade(allele, data) for allele, data in gwas_data.items()}
    gwas_source = GwasSource(
        located=[(a, allele, *gwas_graded[allele]) for a, allele in gwas_located],
        input_rows=gwas["input_rows"],
        dropped=gwas["dropped"],
        version=gwas["release"],
    )

    # Variants ClinVar has not classified: every PASS allele of our VCF that no ClinVar
    # record classifies germline, whatever its shape, through VEP
    # (SpliceAI on the listed genes); our calls at the candidates.
    def build_clinvar_vids(out: Path) -> None:
        found = sorted(clinvar_vids(clinvar38, fasta=fasta))
        with gzip.open(out / CLINVAR_VIDS, "wt", encoding="utf-8") as handle:
            handle.write("".join(f"{vid}\n" for vid in found))

    steps["clinvar_vids"] = _cached(
        work,
        "clinvar_vids",
        CLINVAR_VIDS_STEP_VERSION,
        [shas["clinvar38"], shas["fasta"], bcftools],
        CLINVAR_VIDS,
        build_clinvar_vids,
    )
    listed = sorted({gene for name in LISTED for gene in lists[name].genes})
    steps["vep"] = vep_step(
        work,
        vcf=vcf,
        fasta=fasta,
        files=annotation,
        listed=listed,
        threads=threads,
        parts=[steps["call"].key, IMAGES["vep"], IMAGES["spliceai"], json.dumps(listed)]
        + [shas["fasta"], bcftools]
        + [sha256_cached(p, identities) for p in annotation.required()],
    )
    with gzip.open(steps["vep"].path, "rt", encoding="utf-8") as handle:
        annotated = json.load(handle)
    # ClinVar's filter is here, after the cache: a release classifies other alleles, VEP's
    # records of them stay.
    with gzip.open(steps["clinvar_vids"].path, "rt", encoding="utf-8") as handle:
        classified = frozenset(handle.read().split())
    # chrM is mutserve's (`genokeeper.mito`), not DeepVariant's: none of its records is graded
    # or counted.
    candidates = [
        record
        for record in (VepRecord(**data) for data in annotated["records"])
        if record.vid not in classified
        and predicted_rule(record, grading.predicted) is not None
        and locus(record.vid)[0] != CHRM
    ]
    predicted_alleles = sorted({locus(r.vid) for r in candidates})
    steps["predicted_genotypes"] = genotypes_step(
        work,
        "predicted_genotypes",
        PREDICTED_GENOTYPES_STEP_VERSION,
        alleles=predicted_alleles,
        keyed_by=steps["vep"],
        gvcf=gvcf,
        fasta=fasta,
        haploid=haploid,
        parts=[json.dumps(predicted_alleles), steps["call"].key, str(sex), bcftools]
        + [json.dumps(grading.par if male else {}, sort_keys=True)],
    )
    with gzip.open(steps["predicted_genotypes"].path, "rt", encoding="utf-8") as handle:
        predicted_graded = {
            (chrom, pos, ref, alt): grade((chrom, pos, ref, alt), data)
            for chrom, pos, ref, alt, data in map(json.loads, handle)
        }
    runs = homopolymer_runs(fasta, [r.vid for r in candidates])
    predicted_source = PredictedSource(
        sites=[(r, *predicted_graded[locus(r.vid)], runs[r.vid]) for r in candidates],
        input_rows=sum(
            1
            for vid in annotated["vids"]
            if vid not in classified and locus(vid)[0] != CHRM
        ),
        escapes=load_escapes(annotation.nmd),
        truncations=clinvar_truncations(claims.claims),
        haploinsufficient=dosage.haploinsufficient,
    )

    # Stage 3: repeat expansions, SMN1 and the Y haplogroup on the CRAM, run only with a
    # current stage-3 calibration table; without one each is unassessed with its reason.
    stage3 = Stage3Results(None, None, {}, None, None, None, None, math.nan, {})
    if stage3_table is not None:
        stage3, found = stage3_steps(
            work,
            cram=cram,
            cram_key=steps["align"].key,
            fasta=fasta,
            fasta_sha=shas["fasta"],
            sex=sex,
            mean_depth=metrics["mean_depth"],
            lab=lab,
            grading=grading,
            threads=threads,
        )
        steps |= found

    def not_run(tool: str) -> str | None:
        if stage3_reason is not None:
            return stage3_reason
        failed = stage3.failed.get(tool)
        return None if failed is None else f"{tool} failed: {failed}"

    withheld = {} if stage3_table is None else stage3_table.withheld
    smn_source = SmnSource(
        caller=stage3.caller,
        parascopy=stage3.parascopy,
        not_run=not_run(SMN_TOOL),
        withheld=withheld.get(SMN_KEY),
        min_qual=grading.smn.min_parascopy_qual,
    )
    yleaf_source = YleafSource(
        found=stage3.yleaf,
        not_run=not_run(YLEAF_TOOL),
        withheld=withheld.get(YLEAF_KEY),
        min_qc=grading.yleaf.min_qc_score,
        male=male,
    )

    # Slice 3c: each coverage gene's callability and the sex-chromosome dosage. A gene a
    # stage-3 source states is that source's, whatever its callable share.
    owners = {gene: REPEATS_TOOL for gene in grading.repeats.loci} | {
        "SMN1": SMN_TOOL,
        "SMN2": SMN_TOOL,
    }
    benchmarks = [layout(lab)["truth"] / BENCH_BED]
    if calibration.extra:
        benchmarks.append(hg002_truth(lab)[1])
    foundation, found = foundation_steps(
        work,
        aligned=steps["align"],
        called=steps["call"],
        fasta=fasta,
        lab=lab,
        grading=grading,
        benchmarks=benchmarks,
        owners=owners,
        male=male,
        threads=threads,
    )
    steps |= found
    # The [dosage] bounds come from the foundation table: not current, not graded.
    dosage_failed = foundation.failed.get("dosage")
    dosage_table = open_table(lab, "foundation")
    dosage_source = DosageSource(
        dosage=foundation.dosage,
        not_run=(
            f"dosage failed: {dosage_failed}"
            if dosage_failed is not None
            else dosage_table.reason
            if isinstance(dosage_table, Stale)
            else None
        ),
        thresholds=grading.dosage,
    )

    # The repeat source after 3c's dosage: the X-linked loci follow the run's X dosage,
    # read only under a calibrated [dosage] (else None: "not calibrated"). Warning:
    # `[dosage] calibrated_on` is free text in config/grading.toml, never checked against a
    # calibration report or the masked reference, so an edit that sets it without
    # calibrating lets an unmeasured tolerance decide FMR1's and AFF2's copies.
    loci = grading.repeats.loci
    ratios = foundation.dosage
    x_dosage = (
        XDosage(
            copies(ratios.x_ratio, grading.dosage.tolerance),
            ratios.x_ratio,
            ratios.y_ratio > grading.dosage.y_absent_max,
        )
        if ratios is not None and grading.dosage.calibrated_on
        else None
    )
    graded_keys = {*loci, RFC1_LOCUS, CFTR_GENE}
    repeat_withheld = {g: why for g, why in withheld.items() if g in graded_keys}
    if not male:
        # A two-X sample is graded by the X-linked loci's two-X branch, which the table
        # withholds apart (`<gene>:two_x`).
        for gene in loci:
            branch = withheld.get(two_x_key(gene))
            if branch is not None:
                repeat_withheld.setdefault(gene, branch)
    # The user's library looks PCR-amplified when its flank ratio is nearer HG002's PCR one.
    pcr_like = stage3_table is not None and abs(
        stage3.flank - stage3_table.flank["pcr"]
    ) < abs(stage3.flank - stage3_table.flank["pcrfree"])
    rfc1_failed = stage3.failed.get(RFC1_FAILED)
    repeat_source = RepeatSource(
        calls=stage3.repeats,
        not_run=not_run(REPEATS_TOOL),
        loci=loci,
        withheld=repeat_withheld,
        caps={} if stage3_table is None else stage3_table.caps,
        images=stage3.images,
        pcr_like=pcr_like,
        insert=metrics["median_insert"],
        male=male,
        min_allele_reads=grading.repeats.min_allele_reads,
        min_spanning=grading.exclusion.min_spanning,
        x_dosage=x_dosage,
        rfc1=stage3.rfc1,
        rfc1_not_run=None if rfc1_failed is None else f"rfc1 failed: {rfc1_failed}",
    )

    # Stage 3b: deletions from DELLY and Dysgu, matched by Truvari, run only with a current
    # deletion table; without one, or with a caller failed, the screen is unassessed.
    called: list[Deletion] = []
    shown: list[Shown] = []
    dropped: dict[str, int] = {}
    deletion_screen: Screen | None = None
    if deletion_table is not None:
        deletions, found = deletion_steps(
            work,
            cram=cram,
            cram_key=steps["align"].key,
            fasta=fasta,
            fasta_sha=shas["fasta"],
            grading=grading,
            threads=threads,
            male=male,
        )
        steps |= found
        if deletions.called is None:
            deletion_reason = "; ".join(
                f"{tool} failed: {why}" for tool, why in deletions.failed.items()
            )
        else:
            called = deletions.called
            deletion_screen = load_screen(lab, grading, called, work)
            shown, dropped = screen_all(called, deletion_screen)
    screened = deletion_table is not None and deletion_reason is None
    deletion_p_wrong = (
        deletion_table.p_wrong if deletion_table is not None else lambda _: None
    )

    # Slice 3d: copy number from read depth, insertions, the paralog loci and Lp(a), each
    # run only with its table current; 3b's deletions are read, never waited on.
    slice3d, found3d = slice3d_steps(
        work,
        tables=slice3d_tables,
        lab=lab,
        cram=cram,
        cram_key=steps["align"].key,
        fasta=fasta,
        fasta_sha=shas["fasta"],
        grading=grading,
        threads=threads,
        male=male,
        mean_depth=metrics["mean_depth"],
        dysgu=steps.get("dysgu"),
        deletions=ScreenedDeletions(called, shown, deletion_screen, deletion_p_wrong),
        dosage=dosage,
    )
    # A step stage 3 already built this run (the shared Parascopy depth) keeps its record.
    steps |= {name: step for name, step in found3d.items() if name not in steps}
    deletion_source = DeletionSource(
        shown=shown,
        called=len(called),
        dropped=dropped,
        not_run=deletion_reason,
        p_wrong=deletion_p_wrong,
        haploinsufficient=dosage.haploinsufficient,
        recall={} if deletion_table is None else deletion_table.recall(),
        placement_overrides=screen_overrides(),
        depth_support=slice3d.joined,
    )

    # Slice 3f: HLA, HHV-6, ROH, blood groups and ancestry, typed only with a current
    # slice 3f table (after the deletion screen and pgsc_calc, whose population ROH reads);
    # the HLA steps are `outside_steps`' own, reused. A failed tool leaves its statements
    # unassessed with its line.
    typing = TypingResults()
    if slice3f_table is not None:
        typing, found = typing_steps(
            work,
            cram=cram,
            cram_key=steps["align"].key,
            fasta=fasta,
            fasta_sha=shas["fasta"],
            vcf=vcf,
            gvcf=gvcf,
            # The calls and their grading: the call step, the calibration and the grading.
            calls_key=_key(
                steps["call"].key,
                json.dumps(calibration.record, sort_keys=True),
                _grading_key(grading),
            ),
            mean_depth=metrics["mean_depth"],
            population=None if similarity is None else similarity.population,
            male=male,
            lab=lab,
            grading=grading,
            threads=threads,
            table=slice3f_table,
            level=call_level,
            tools=versions,
        )
        # outside_steps' HLA steps keep their own record (this run's, not "reused").
        steps |= {name: step for name, step in found.items() if name not in steps}

    def typing_not_run(*tools: str) -> str | None:
        if slice3f_reason is not None:
            return slice3f_reason
        failed = [
            f"{t} failed: {typing.failed[t]}" for t in tools if t in typing.failed
        ]
        return "; ".join(failed) or None

    no_population = None if similarity is not None else f"not run: {similarity_failed}"
    typing_sources: list[Source] = [
        HlaSource(
            optitype=typing.optitype,
            t1k=typing.t1k,
            tags=tag_dosages,
            not_run=typing_not_run(MHC_READS, OPTITYPE, T1K),
            withheld=withheld3f,
        ),
        MitomapSource(
            calls=mito_calls,
            coverage=coverage,
            clinvar_vids=frozenset(
                c.vid
                for c in claims.claims
                if c.chrom == CHRM and clnsig_class(c.clnsig) == "pathogenic"
            ),
            floor=floor,
            min_homoplasmy=homoplasmy,
            min_dp=grading.min_dp,
            heteroplasmy_withheld=heteroplasmy_withheld,
            strata=index,
        ),
        Hhv6Source(
            call=typing.hhv6,
            not_run=typing_not_run(HHV6_STEP),
            withheld=withheld3f.get(HHV6_KEY),
        ),
        RohSource(
            summary=None
            if typing.roh is None
            else summarise(typing.roh, grading=grading.roh),
            not_run=typing_not_run(ROH_STEP) or no_population,
            withheld=withheld3f.get(ROH_KEY),
        ),
        BloodGroupSource(
            calls=typing.bloodgroups,
            not_run=typing_not_run(RBCEQ2_STEP),
            withheld=withheld3f,
            uncovered=typing.bloodgroup_uncovered,
            unchecked=typing.bloodgroup_unchecked,
        ),
        AncestrySource(
            similarity=similarity,
            admixture=typing.admixture,
            archaic=typing.archaic,
            not_run={
                name: why
                for name, why in {
                    "similarity": similarity_failed,
                    "admixture": typing_not_run(ADMIXTURE_STEP),
                    "neanderthal": typing_not_run(ARCHAIC_STEP),
                }.items()
                if why is not None
            },
            withheld=withheld3f,
        ),
        TraitSource(graded=typing_graded),
    ]

    # Read-backed phasing of the candidate pairs, through the run's heterozygous SNVs between
    # each pair's sites: the rows the store phase collects again, collected here under a
    # placeholder run id to find them.
    pending = RunContext(
        run_id="pending",
        lists=lists,
        grading=grading,
        curations=curations,
        gene_level=validity_by_gene,
        faf=faf,
        ba1_exceptions=ba1_exceptions,
        inheritance=inheritance,
        male=male,
        haploid=haploid,
        source_version=release,
        validity_source=validity_source,
    )
    pairs = promote.candidate_pairs(
        [source.collect(pending), predicted_source.collect(pending)],
        lists,
        inheritance,
        haploid,
    )
    # An APOE double heterozygote's diplotype rests on its phase (sites 138 bp apart).
    apoe_dosages = {(c.chrom, c.pos, c.ref, c.alt): d.dosage for c, d, _ in apoe}
    if apoe_dosages.get(APOE_E4) == apoe_dosages.get(APOE_E2) == "1":
        pairs = sorted({*pairs, (APOE_E4, APOE_E2)})

    def build_phases(out: Path) -> None:
        found = phase_pairs(cram, fasta, vcf, pairs)
        (out / "phases.json").write_text(
            json.dumps([[list(a), list(b), str(found[(a, b)])] for a, b in pairs])
        )

    phases: dict[Pair, Phase] = {}
    if pairs:
        try:
            steps["phases"] = _cached(
                work,
                "phases",
                PHASES_STEP_VERSION,
                [
                    steps["align"].key,
                    steps["call"].key,
                    IMAGES["whatshap"],
                    shas["fasta"],
                ]
                + [json.dumps(pairs)],
                "phases.json",
                build_phases,
            )
        except DockerUnavailable as e:
            raise RunRefused(_one_line(e)) from e
        for a, b, phase in json.loads(steps["phases"].path.read_text()):
            phases[(_as_allele(a), _as_allele(b))] = Phase(phase)

    prs_source = PrsSource(
        scores=prs.scores,
        results=prs_results,
        apoe=apoe,
        apoe_phase=phases.get((APOE_E4, APOE_E2)),
        not_covered=load_not_covered(),
    )
    inputs: dict[str, object] = {
        "fastq1": {**_file(fastq1, shas["fastq1"]), "size": fastq1.stat().st_size},
        "fastq2": {**_file(fastq2, shas["fastq2"]), "size": fastq2.stat().st_size},
        "fasta": _file(fasta, shas["fasta"]),
        "clinvar38": _file(clinvar38, shas["clinvar38"]),
        "strata": {flag: _file(p, sha256_file(p)) for flag, p in beds.items()},
        "platform": asdict(platform),
        "qc": {"metrics": metrics},
        "sex": {"value": sex, "decided_by": decided_by},
        "gene_validity": {
            **_file(validity_csv, sha256_file(validity_csv)),
            "date": validity_date,
            "present": True,
        },
        "grading": grading_input,
        "faf": {
            **_file(faf.manifest_path, faf.sha256),
            "release": faf.manifest["release"],
            "field": faf.manifest["field"],
            "floor": faf.manifest["floor"],
        },
        "cgd": {**_file(cgd, shas["cgd"]), "present": True},
        "gencc": gencc_input,
        "pgx": {
            **pgx.record(),
            "phase": phase_reason,
            "ugt1a1_28_trusted": bool(trusted),
        },
        "prs": prs.record(shas["pgs_panel"]),
        "gwas": {
            **_file(catalog, shas["gwas"]),
            "release": gwas["release"],
            "panel_variants": _file(variants, shas["panel_pvar"]),
        },
        "vep": {
            str(p.relative_to(annotation.folder)): _file(
                p, sha256_cached(p, identities)
            )
            for p in annotation.required()
        },
        "clingen_dosage": {
            **_file(dosage_list, sha256_file(dosage_list)),
            "date": dosage.date,
        },
        "calibration": calibration.record,
        "stage3": {
            "table": None if stage3_table is None else stage3_table.versions,
            "reason": stage3_reason,
        },
        "deletions": {
            "table": None if deletion_table is None else deletion_table.versions,
            "reason": deletion_reason,
        },
        "slice3f": {"table": slice3f_calibrated, "reason": slice3f_reason},
        "haplocheck": haplocheck_input,
        "identity": identity,
        "coverage": _failed_step(foundation.failed, "coverage"),
        "dosage": _failed_step(foundation.failed, "dosage"),
        "slice3d": {
            **{
                part: {
                    "table": slice3d_tables.versions.get(part)
                    if part not in slice3d_tables.reasons
                    else None,
                    "reason": slice3d_tables.reasons.get(part),
                    "ran": part in slice3d.ran,
                }
                for part in SLICE3D_PARTS
            },
            "clingen_regions": slice3d.regions,
        },
        "phased_pairs": {
            str(phase): sum(1 for p in phases.values() if p == phase) for phase in Phase
        },
        "steps": {name: step.record() for name, step in steps.items()},
    }
    tools = {
        "genokeeper": __version__,
        "analyser": analyser_identity(),
        **versions,
        "clinvar": release,
        "clingen_gene_validity": validity_date,
        "gnomad_faf": f"{faf.manifest['release']}, floor {faf.manifest['floor']}",
        "deepvariant_image": IMAGES["deepvariant"],
        "verifybamid2_image": IMAGES["verifybamid2"],
        "mutserve_image": IMAGES["mutserve"],
        "haplogrep3": HAPLOGREP_VERSION,
        "haplogrep3_tree": TREE,
        "haplogrep3_image": IMAGES["haplogrep3"],
        "whatshap_image": IMAGES["whatshap"],
        "pharmcat": PHARMCAT_VERSION,
        "pharmcat_image": IMAGES["pharmcat"],
        "bcyrius": BCYRIUS_VERSION,
        "bcyrius_image": IMAGES["bcyrius"],
        "aldy": ALDY_VERSION,
        "aldy_image": IMAGES["aldy"],
        "optitype_image": IMAGES["optitype"],
        "t1k_image": IMAGES["t1k"],
        "pgsc_calc": f"{PGSC_CALC_VERSION} ({PGSC_CALC_COMMIT})",
        "nextflow": NEXTFLOW_VERSION,
        **{f"{name}_image": IMAGES[name] for name in PGSC_LABELS.values()},
        "ipd_imgt_hla": IMGT_HLA_VERSION,
        "cpic": pgx.cpic.name,
        "pharmgkb_clinical_annotations": pgx.annotations.created,
        "gwas_catalog": gwas["release"],
        "vep": VEP_VERSION,
        "vep_cache": VEP_CACHE_VERSION,
        "vep_image": IMAGES["vep"],
        "spliceai_image": IMAGES["spliceai"],
        "clingen_dosage": dosage.date,
        "expansionhunter": EH_VERSION,
        "reviewer": REVIEWER_VERSION,
        "smncaller": SMNCALLER_VERSION,
        "parascopy": PARASCOPY_VERSION,
        "yleaf": YLEAF_VERSION,
        **{
            f"{name}_image": IMAGES[name]
            for name in (
                "expansionhunter",
                "reviewer",
                "smncaller",
                "parascopy",
                "yleaf",
            )
        },
        **{
            f"{name}_image": IMAGES[name]
            for name in ("haplocheck", "rbceq2", "fastmixture", "hmmix")
        },
        # ROH's bcftools is the truvari image's, named when ROH ran.
        **({"truvari_image": IMAGES["truvari"]} if "roh" in steps else {}),
        # The deletion tools only when the screen ran (`query --json`'s `versions`).
        **(
            {
                "delly": DELLY_VERSION,
                "dysgu": DYSGU_VERSION,
                "truvari": TRUVARI_VERSION,
                **{f"{n}_image": IMAGES[n] for n in ("delly", "dysgu", "truvari")},
            }
            if screened
            else {}
        ),
        # Slice 3d's own tools only for the parts that ran (3b's DELLY, Dysgu and Truvari
        # keep meaning 3b's screen ran); `delly cnv` is the key the genome-ask skill reads
        # to tell "screened, none found" from "not run".
        **(
            {
                "delly cnv": DELLY_VERSION,
                **(
                    {"cnvpytor": CNVPYTOR_VERSION, "cnvpytor_image": IMAGES["cnvpytor"]}
                    if "cnvpytor" in found3d
                    else {}
                ),
            }
            if "cnv" in slice3d.ran
            else {}
        ),
        **(
            {
                "insurveyor": INSURVEYOR_VERSION,
                "scramble": SCRAMBLE_VERSION,
                **{f"{n}_image": IMAGES[n] for n in ("insurveyor", "scramble")},
            }
            if "insertions" in slice3d.ran
            else {}
        ),
        **(
            {
                "gauchian": GAUCHIAN_VERSION,
                "kilda": KILDA_VERSION,
                **{f"{n}_image": IMAGES[n] for n in ("gauchian", "kilda")},
            }
            if PARALOGS in slice3d.ran
            else {}
        ),
        # The six versions the calibration table was checked against (bwa, samtools and
        # the DeepVariant image among the keys above, with the same values).
        **calibrated,
        **{name: gl.version for name, gl in lists.items()},
        "identity": f"{IDENTITY_TOOL}: {identity_text}",
        "reference": MASKED_REFERENCE,
    }
    deepvariant = _header_value(vcf, "DeepVariant_version")
    if deepvariant is not None:
        tools["deepvariant"] = deepvariant

    try:
        con = connect(store)
    except duckdb.IOException as e:
        raise _in_use(store, e) from e
    except (StoreTooNew, StoreTooOld) as e:
        raise RunRefused(str(e)) from e
    writer = StoreWriter(con, work)
    try:
        with transaction(con), writer:
            run_id = start_run(con, inputs, tools)
            ctx = RunContext(
                run_id=run_id,
                lists=lists,
                grading=grading,
                curations=curations,
                gene_level=validity_by_gene,
                faf=faf,
                ba1_exceptions=ba1_exceptions,
                inheritance=inheritance,
                male=male,
                haploid=haploid,
                source_version=release,
                cpic_levels=pgx.cpic_levels,
                validity_source=validity_source,
            )
            store_sources(
                writer,
                con,
                ctx,
                [
                    source,
                    pharmgkb_source,
                    research_source,
                    pgx_source,
                    prs_source,
                    mito_source,
                    gwas_source,
                    predicted_source,
                    repeat_source,
                    smn_source,
                    yleaf_source,
                    deletion_source,
                    dosage_source,
                    slice3d.cnv,
                    slice3d.insertions,
                    slice3d.paralogs,
                    slice3d.lpa,
                    *typing_sources,
                ],
                phases,
                gauchian=slice3d.gauchian,
            )
            for row in foundation.coverage or ():
                writer.add("coverage", gene_coverage.coverage_fields(run_id, row))
            check_calls(con, run_id)
        # The run's snapshot, and any an earlier export left out by failing after its commit.
        export_snapshots(con, paths.reports / "snapshots")
    finally:
        con.close()
    return run_id


# --- command line -------------------------------------------------------------------------


def newest_gene_validity(lab: Path) -> Path:
    """The ClinGen gene validity export of the newest `<lab>/clingen/<YYYYMMDD>/`
    holding one."""
    if (folder := _newest_dated(lab / "clingen", [GENE_VALIDITY_FILE])) is None:
        raise RunRefused(
            f"missing input: no {lab}/clingen/<YYYYMMDD>/{GENE_VALIDITY_FILE};"
            " run scripts/bootstrap.sh"
        )
    return folder / GENE_VALIDITY_FILE


def open_faf_table(lab: Path) -> FafTable:
    """The lab's gnomAD FAF table (`<lab>/gnomad/4.1/`); refuse (`RunRefused`) one that
    is absent or incomplete, naming the bootstrap step that builds it."""
    folder = faf_folder(lab)
    try:
        return FafTable.open(folder)
    except FafMissing as e:
        raise RunRefused(
            f"missing input: no complete gnomAD FAF table at {folder} ({_one_line(e)});"
            " run scripts/bootstrap.sh --with-gnomad"
        ) from e


def accounting_table(con: duckdb.DuckDBPyConnection, run_id: str) -> str:
    """The run's accounting, one line per tool: counts and drop reasons only."""
    rows = con.execute(
        "SELECT tool, input_rows, findings, watch, non_findings, unassessed, dropped"
        " FROM accounting WHERE run_id = ? ORDER BY tool",
        [run_id],
    ).fetchall()
    lines = ["tool\tinput_rows\tfindings\twatch\tnon_findings\tunassessed\tdropped"]
    for tool, input_rows, findings, watch, non_findings, unassessed, dropped in rows:
        reasons = {k: v for k, v in json.loads(dropped).items() if v}
        detail = ", ".join(f"{k} {v}" for k, v in sorted(reasons.items()))
        total = sum(reasons.values())
        lines.append(
            f"{tool}\t{input_rows}\t{findings}\t{watch}\t{non_findings}\t{unassessed}"
            f"\t{total}" + (f" ({detail})" if detail else "")
        )
    return "\n".join(lines)


def newest_clinvar38(lab: Path) -> Path:
    """The GRCh38 VCF of the newest `<lab>/clinvar/<YYYYMMDD>/` holding one."""
    if (folder := _newest_dated(lab / "clinvar", ["clinvar_GRCh38.vcf.gz"])) is None:
        raise _bootstrap_missing(f"no {lab}/clinvar/<YYYYMMDD>/clinvar_GRCh38.vcf.gz")
    return folder / "clinvar_GRCh38.vcf.gz"


def fastq_pair(raw: Path) -> tuple[Path, Path]:
    """The delivered paired FASTQ: the two `*.fq.gz` or `*.fastq.gz` files in `raw`, in
    name order; refuse any other count."""
    found = sorted({*raw.glob("*.fq.gz"), *raw.glob("*.fastq.gz")})
    if len(found) != 2:
        raise RunRefused(
            f"missing input: expected two FASTQ files (*.fq.gz, *.fastq.gz) in {raw},"
            f" found {len(found)}; name them with --fastq1 and --fastq2"
        )
    return found[0], found[1]


def _genome_command(args: argparse.Namespace) -> int:
    # Everything a run writes (the store, snapshots, work artefacts) is the user's alone.
    previous = os.umask(0o077)
    try:
        return _genome(args)
    finally:
        os.umask(previous)


def _genome(args: argparse.Namespace) -> int:
    paths = load_paths()
    try:
        if (args.fastq1 is None) != (args.fastq2 is None):
            raise RunRefused("--fastq1 and --fastq2 go together: give both, or neither")
        pair = (
            (args.fastq1, args.fastq2)
            if args.fastq1 is not None and args.fastq2 is not None
            else fastq_pair(paths.raw)
        )
        try:
            delivered = args.genotypes or find_delivered(paths.raw)
        except ValueError as e:
            raise RunRefused(
                f"{_one_line(e)}; name the one to compare with --genotypes"
            ) from e
        faf = open_faf_table(paths.lab)
        run_id = run_genome(
            paths,
            fastq1=pair[0],
            fastq2=pair[1],
            clinvar38=args.clinvar38 or newest_clinvar38(paths.lab),
            sex=args.sex,
            faf=faf,
            threads=args.threads,
            delivered=delivered,
        )
    except RunRefused as e:
        print(f"genokeeper run: {e}", file=sys.stderr)
        return 2
    with duckdb.connect(str(paths.reports / STORE_NAME), read_only=True) as con:
        table = accounting_table(con, run_id)
        row = con.execute(
            "SELECT inputs, tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
        covered = con.execute(
            "SELECT count(*), count(*) FILTER (WHERE NOT reliable) FROM coverage"
            " WHERE run_id = ?",
            [run_id],
        ).fetchone()
    inputs = json.loads(row[0]) if row is not None else {}
    tools = json.loads(row[1]) if row is not None else {}
    metrics = inputs.get("qc", {}).get("metrics", {})
    print(f"run {run_id}")
    print(table)
    print("QC: " + ", ".join(f"{k} {v:.4g}" for k, v in sorted(metrics.items())))
    print(f"population frequency: {faf.describe()}")
    print(f"identity: {tools.get('identity')}")
    failed = inputs.get("coverage", {}).get("failed")
    genes, below = covered if covered is not None else (0, 0)
    min_callable = load_grading().coverage.min_callable
    print(
        f"coverage: not measured: {failed}"
        if failed is not None
        else f"coverage: {genes} genes, {below} below {min_callable:.0%} callable"
    )
    report = inputs.get("steps", {}).get("pharmcat")
    if report is not None:
        # The user's to read; it holds PharmCAT's recommendations, which the store does not.
        print(f"PharmCAT report: {Path(report['path']).with_suffix('.html')}")
    return 0


def register(subparsers: Subparsers) -> None:
    """Add the `run` sub-command."""
    g = subparsers.add_parser(
        "run",
        help="align, QC, call and grade every ClinVar GRCh38 claim into the store",
    )
    g.add_argument("--fastq1", type=Path, help="default: the first of data/raw/*.fq.gz")
    g.add_argument(
        "--fastq2", type=Path, help="default: the second of data/raw/*.fq.gz"
    )
    g.add_argument(
        "--clinvar38",
        type=Path,
        help="default: the newest <lab>/clinvar/<date>/clinvar_GRCh38.vcf.gz",
    )
    g.add_argument(
        "--sex",
        choices=["male", "female"],
        help="default: the X copy count (one X male, two female), from chrX depth",
    )
    g.add_argument("--threads", type=int, help="default: this machine's CPUs")
    g.add_argument(
        "--genotypes",
        type=Path,
        help="the delivered 23andMe-format genotype file the identity check compares"
        " with; default: the one data/raw/*.txt with its column line",
    )
    g.set_defaults(handler=_genome_command)
