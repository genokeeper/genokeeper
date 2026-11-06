"""The calibration tables' builders and versions (`ROWS`, which `genokeeper.calibrate.tables`
registers) and `genokeeper calibrate`, the one command over them: HG003 through the run's
pipeline, scored against GIAB's truth, into the calibration table every `genokeeper run`
grades with (track design, "2b. Calibration"), and the tables slice 3 measures on HG002.

`build_table` runs the downsampled HG003 pair through `run_sample` (aligned, QC-gated and
called as the user's sample is, each step cached) and records it as `hg003`'s
`aligned.json`, or, once that record is on this reference (HG003 re-aligned from its CRAM by
`python -m genokeeper.calibrate.hg002 realign hg003`), reads its gVCF from there with no
FASTQ; it scores the gVCF per chromosome
(`count_cells`), and writes the table with the versions of this machine's pipeline
(`pipeline.calibration_versions`, the same function the run compares with) to
`<lab>/calibration/small_variants/` (`tables.write_table`), whole or not at all. A missing
input refuses before anything runs.

`build_hg002_cells` scores HG002 PCR-free's gVCF where HG003's truth does not reach
(`genokeeper.calibrate.hg002_cells`) into `<lab>/calibration/hg002_cells/`.

`build_stage3` runs HG002's two aligned libraries through the run's stage-3
steps (`pipeline.stage3_steps`), HG003 and HG004 through its repeat step alone
(`pipeline.repeat_steps`), and writes the stage-3 table (`genokeeper.calibrate.stage3`) to `<lab>/calibration/stage3/`,
whole or not at all; a sample whose tool failed writes none.

`build_deletions` runs the same two libraries through the run's deletion
steps (`pipeline.deletion_steps`), compares their rare calls with HG002's v5.0q structural-
variant truth by Truvari, and writes the deletion table (`genokeeper.calibrate.deletions`) to
`<lab>/calibration/deletions/`, whole or not at all; a library whose step failed writes none.

`report` prints the table, the versions it was measured under, and the held-out checks the
spec asks for: cells measured on the odd chromosomes checked on the even ones, and on
chr1-19 checked on chr20-22 (`held_out`). A held-out cell is outside when its errors fall
outside the 95% posterior-predictive interval of its training cell for its own units, with
the dispersion between the training chromosomes; the part fails when more than 10% of the
checked cells are (Task 21 reads the report).

`build_cnv`, `build_insertions` and `build_paralogs` measure slice 3d's three tables
(`genokeeper.calibrate.slice3d`) on the calibration libraries' CRAMs (`slice3d.TRIO`; HG002's
two alone for insertions) through the run's 3d steps, each cached in the library's folder
under its run lock, and write each whole or not at all; a step failing on any library
refuses that table alone. `judge_library` is 3b's judging against v5.0q, which
`build_deletions` and the 3d builders share.

Each table is written with its report (`report.txt`), which `--report` prints.

    genokeeper calibrate [--only NAME] [--report] [--threads N]
"""

import argparse
import hashlib
import itertools
import json
import math
import os
import shutil
from collections import defaultdict
from collections.abc import Callable, Collection, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pysam

from genokeeper.bloodgroups import DEPENDS as BLOODGROUP_DEPENDS
from genokeeper.bloodgroups import (
    POSITIONS_FILE,
    SystemCall,
    rbceq2_positions,
    read_rbceq2,
    run_rbceq2,
    uncovered,
)
from genokeeper.bloodgroups import holes as bloodgroup_holes
from genokeeper.bloodgroups import systems as bloodgroup_systems
from genokeeper.calibrate import (
    BENCH_BED,
    TRUTH_VCF,
    hg003_dir,
    layout,
    rtg_path,
    table_dir,
)
from genokeeper.calibrate.classes import (
    OUTSIDE_BENCHMARK,
    ClassTable,
    outside_benchmark,
    write_class_table,
)
from genokeeper.calibrate.classes import measure as measure_classes
from genokeeper.calibrate.classes import summary as summarise_classes
from genokeeper.calibrate.deletions import (
    BANDS,
    DeletionTable,
    band_of,
    class_of,
    deletion_versions,
    deletions_dir,
    expected_true_carried,
    geometry,
    read_population,
    summary,
)
from genokeeper.calibrate.deletions import measure as measure_deletions
from genokeeper.calibrate.deletions import report as deletions_report
from genokeeper.calibrate.deletions import write_table as write_deletions
from genokeeper.calibrate.foundation import build_foundation, foundation_versions
from genokeeper.calibrate.hg002 import (
    LIBRARIES,
    Record,
    hg002_dir,
    library_dir,
    library_sex,
    read_aligned,
    read_called,
    read_record,
    write_record,
)
from genokeeper.calibrate.hg002_cells import (
    build_cells,
    cells_dir,
    cells_report,
    lab_regions,
    write_cells,
)
from genokeeper.calibrate.measure import Cells, count_cells
from genokeeper.calibrate.records import (
    RECORDED,
    label,
    on_this_reference,
    reference_sha256,
)
from genokeeper.calibrate.records import label as lib_label
from genokeeper.calibrate.repeat_rules import (
    SampleRepeats,
    measure_repeat_rules,
    merge_withheld,
)
from genokeeper.calibrate.sample import downsample, run_sample
from genokeeper.calibrate.slice3d import (
    CNV_BANDS,
    CNV_FEATURES,
    INSERTION_BANDS,
    INSERTION_FEATURES,
    PARALOGS,
    TRIO,
    ParalogTable,
    classify_gain_truth,
    cnv_versions,
    cnvpytor_kept,
    insertion_versions,
    large_withheld,
    load_kiv2_truth,
    load_truth,
    paralog_measure,
    paralog_versions,
    pool_large,
    spike_recall,
    trio_judged,
    write_paralog_table,
)
from genokeeper.calibrate.slice3d import report as slice3d_report
from genokeeper.calibrate.slice3f import (
    ARCHAIC_KEY,
    CHRM_KEY,
    GRID,
    HHV6_KEY,
    PUBLISHED_ARCHAIC,
    ROH_KEY,
    SPIKE,
    Slice3fTable,
    calls_key,
    consistency_ancestry,
    grid_step,
    measure_archaic,
    measure_bloodgroups,
    measure_chrm_floor,
    measure_haplocheck,
    measure_hhv6,
    measure_hla,
    measure_phase,
    measure_roh,
    measure_ugt1a1,
    mixture_level,
    slice3f_versions,
    ta_genotype,
    write_slice3f,
)
from genokeeper.calibrate.slice3f import LIBRARIES as SLICE3F_LIBRARIES
from genokeeper.calibrate.slice3f import truths as slice3f_truths
from genokeeper.calibrate.spikein import pick_spikes, write_spiked
from genokeeper.calibrate.stage3 import (
    RULES,
    Stage3Table,
    hg002_truth,
    measure,
    stage3_dir,
    stage3_versions,
    truth_repeats,
)
from genokeeper.calibrate.stage3 import report as stage3_report
from genokeeper.calibrate.stage3 import write_table as write_stage3
from genokeeper.calibrate.table import (
    GENOTYPE,
    KINDS,
    VERSION_KEYS,
    CalibrationTable,
    CellKey,
    Counts,
    beta_binomial_interval,
    dispersion,
)
from genokeeper.calibrate.tables import REPORT, TABLES, calibrate, write_table
from genokeeper.checksum import sha256_cached, sha256_file
from genokeeper.cnv.screen import (
    DEPTH_CALLERS,
    MATCH_OVERLAP,
    CnvScreen,
    LargeEvent,
    cnv_band,
    gene_spans,
    owned_regions,
    screen_cnv,
)
from genokeeper.config import (
    CnvThresholds,
    DeletionThresholds,
    Grading,
    Paths,
    grch38_fasta,
    load_grading,
    load_paths,
)
from genokeeper.containers import IMAGES, DockerUnavailable
from genokeeper.deletions.callers import read_deletions, vcf_records
from genokeeper.deletions.exons import PRIMARY, Target, load_exons, load_spans, targets
from genokeeper.deletions.files import sv_files, truth_files
from genokeeper.deletions.screen import (
    CallClass,
    Deletion,
    Screen,
    _Index,  # pyright: ignore[reportPrivateUsage]
    classify,
    is_short,
    rare,
)
from genokeeper.deletions.screen import key as class_key
from genokeeper.deletions.truvari import (
    AGAINST_TRUTH,
    against_truth,
    judged,
    write_calls,
)
from genokeeper.evidence.clingen import DOSAGE_FILE, load_dosage
from genokeeper.genotype.gvcf import Allele, Called, calls_at
from genokeeper.genotype.reliability import StrataIndex, on_haploid_contig
from genokeeper.haplocheck import Haplocheck, run_haplocheck
from genokeeper.insertions import INSERTION_CALLERS, SCRAMBLE_MARGIN, scramble_regions
from genokeeper.masking import (
    correct_copies_sha256,
    masked_bed,
    without_false_duplications,
)
from genokeeper.mito import VCF_FILE, MitoCall, call_chrm
from genokeeper.paralogs import LOCI, GauchianCall
from genokeeper.pgx.hla_typing import DQ25_TAG, TAG_ALLELES, TAG_RSIDS, load_hg002_truth
from genokeeper.pgx.pharmcat import (
    UGT1A1_28,
    gene_regions,
    pharmcat_positions,
    positions_path,
    write_pharmcat_vcf,
)
from genokeeper.phasing import het_snvs, phase_sites
from genokeeper.pipeline import (
    GENOME_SAMPLE,
    IDENTITIES,
    PARASCOPY_DEPTH,
    RunRefused,
    Stage3Results,
    TypingResults,
    _attempt,  # pyright: ignore[reportPrivateUsage]
    _cached,
    _clingen_regions,  # pyright: ignore[reportPrivateUsage]
    _newest_dated,  # pyright: ignore[reportPrivateUsage]
    _Step,  # pyright: ignore[reportPrivateUsage]
    calibration_versions,
    cnv_steps,
    deletion_steps,
    dysgu_ins_dup_step,
    fastq_pair,
    genome_tool_versions,
    insertion_steps,
    load_screen,
    paralog_steps,
    parascopy_depth_step,
    recorded_bases,
    repeat_steps,
    run_lock,
    stage3_steps,
    typing_steps,
)
from genokeeper.regions import Interval, giab38_dir, merge, read_bed, subtract
from genokeeper.repeats.expansionhunter import region
from genokeeper.roh import (
    Run,
    inside,
    read_runs,
    run_roh,
    spike,
    truth_gt,
    truth_sites_vcf,
)
from genokeeper.sex import Sex
from genokeeper.smn import parascopy_data
from genokeeper.sources.paralogs import (
    _usable,  # pyright: ignore[reportPrivateUsage]
)
from genokeeper.store.db import Level, RowType
from genokeeper.tools import ToolFailed, run_tool
from genokeeper.typing_files import typing_files

if TYPE_CHECKING:
    from genokeeper.cli import Subparsers

# GIAB's HG003 (NA24149) is the Ashkenazi trio's father, HG002 (NA24385) its son.
HG003_SEX: Sex = "male"
# A male's haploid chromosomes: one haplotype of truth there.
HAPLOID = ("chrX", "chrY")
# The held-out splits of the spec ("Checks"): (training, held out) chromosome names.
AUTOSOMES = tuple(f"chr{n}" for n in range(1, 23))
SPLITS = (
    (
        "odd chromosomes",
        tuple(c for c in AUTOSOMES if int(c[3:]) % 2),
        "even chromosomes",
        tuple(c for c in AUTOSOMES if not int(c[3:]) % 2),
    ),
    ("chr1-19", AUTOSOMES[:19], "chr20-22", AUTOSOMES[19:]),
)
# The spec's limit: more held-out cells outside their interval fails the part.
MAX_OUTSIDE = 0.10
# The spec's other limits: HG003's SNP or INDEL F1 below these fails the part.
MIN_F1 = {"snp": 0.995, "indel": 0.99}


def merged(
    by_chrom: Mapping[str, Cells], chroms: Collection[str] | None = None
) -> Cells:
    """The counts of `chroms` (None: every chromosome) summed per cell."""
    total: Cells = {}
    for chrom, cells in by_chrom.items():
        if chroms is None or chrom in chroms:
            for key, found in cells.items():
                total.setdefault(key, Counts()).merge(found)
    return total


@dataclass(frozen=True)
class HeldOut:
    """A held-out check of genotype error: the held-out cells checked against the training
    table, those outside their interval, and those the training half never measured."""

    checked: int
    outside: list[CellKey]
    unmeasured: list[CellKey]


def held_out(
    by_chrom: Mapping[str, Cells], *, train: Collection[str], test: Collection[str]
) -> HeldOut:
    """The cells measured on the `train` chromosomes checked on the `test` ones: a test
    cell is outside when its errors fall outside the central 95% of what its training
    cell predicts for its units, the beta-binomial under the training cell's Jeffreys
    posterior Beta(errors + 1/2, total - errors + 1/2). Both are taken at effective size,
    every count divided by the training chromosomes' dispersion phi (`dispersion`), and
    the interval scaled back by phi: the test half's own spread between chromosomes is
    then inside it with the training half's, so a table right in every cell, its rates
    varying between chromosomes as much as on the training ones, has about 5% of its
    cells outside. Only the training chromosomes set phi."""
    trained = merged(by_chrom, train)
    checked, outside, unmeasured = 0, [], []
    tested = sorted(
        (key, found)
        for (kind, key), found in merged(by_chrom, test).items()
        if kind == GENOTYPE and found.total > 0
    )
    for key, found in tested:
        cell = trained.get((GENOTYPE, key))
        if cell is None or cell.total == 0:
            unmeasured.append(key)
            continue
        checked += 1
        phi = dispersion(
            (cells[(GENOTYPE, key)].errors, cells[(GENOTYPE, key)].total)
            for chrom, cells in by_chrom.items()
            if chrom in train and (GENOTYPE, key) in cells
        )
        a, b = cell.errors / phi + 0.5, (cell.total - cell.errors) / phi + 0.5
        lo, hi = beta_binomial_interval(max(1, round(found.total / phi)), a, b)
        if not lo * phi <= found.errors <= hi * phi:
            outside.append(key)
    return HeldOut(checked=checked, outside=outside, unmeasured=unmeasured)


def report(table: CalibrationTable, by_chrom: Mapping[str, Cells]) -> str:
    """The table (each error kind's cells), its versions, HG003's SNP and INDEL F1, and
    both held-out checks, each limit of the part named beside what it judges."""
    lines = ["versions measured under:"]
    lines += [f"  {key}: {table.versions[key]}" for key in VERSION_KEYS]
    columns = (
        "row_type\tallele_class\tregion\tgq_bin\terrors\ttotal\trate\tlo\thi\tdepth"
    )
    for kind in KINDS:
        lines += ["", f"{kind} error per cell:", columns]
        for key in table.measured(errors=kind):
            rate = table.rate_for(key, errors=kind)
            depth = "" if rate.depth is None else f"{rate.depth:.1f}"
            lines.append(
                f"{key.row_type}\t{key.allele_class}\t{key.region}\t{key.gq_bin}"
                f"\t{rate.errors}\t{rate.total}\t{rate.rate:.3g}\t{rate.lo:.3g}"
                f"\t{rate.hi:.3g}\t{depth}"
            )
    lines.append("")
    for allele_type, limit in MIN_F1.items():
        f1 = table.f1.get(allele_type, math.nan)  # nan: not measured, so below
        below = "" if f1 >= limit else "; below the limit"
        lines.append(
            f"{allele_type.upper()} F1 {f1:.4f} (the part stops below {limit}){below}"
        )
    for train_name, train, test_name, test in SPLITS:
        found = held_out(by_chrom, train=train, test=test)
        share = len(found.outside) / found.checked if found.checked else 0.0
        lines.append(
            f"held out, {train_name} -> {test_name}: {found.checked} cells checked,"
            f" {len(found.outside)} outside their interval ({share:.1%}),"
            f" {len(found.unmeasured)} not measured in the training half"
            + (f"; above the {MAX_OUTSIDE:.0%} limit" if share > MAX_OUTSIDE else "")
        )
    return "\n".join(lines)


def _bootstrap_missing(what: str) -> RunRefused:
    return RunRefused(
        f"missing input: {what}; run scripts/bootstrap.sh --with-calibration"
    )


def build_table(
    paths: Paths,
    *,
    grading: Grading | None = None,
    threads: int | None = None,
) -> tuple[CalibrationTable, dict[str, Cells]]:
    """Measure HG003 and write the calibration table; return it with the per-chromosome
    counts it was built from. With `hg003`'s record on this reference and called (`realign
    hg003`), its gVCF is scored as it is (`read_called`): nothing is aligned and no FASTQ is
    read. A record on another reference is `realign hg003`'s alone (its count checks), so
    it refuses naming that command. Without a record (a fresh lab) the pair is the two
    FASTQ in the downsampled folder (`layout`), written there first by `downsample` from
    HG003's full FASTQ to the bases the user's run recorded (`recorded_bases`), the full
    FASTQ deleted once both files have landed; HG003 must show as male (`run_sample`), and
    its record is written after the call. `grading` defaults to config/grading.toml (its
    [strata] are the difficult regions, its [qc] the gate), `threads` to this machine's
    CPUs.

    Raises `RunRefused` on a missing input (naming the bootstrap) before anything runs, and
    as `run_sample` does; nothing is written then, or when scoring fails."""
    lab = paths.lab
    where = layout(lab)
    for what, path in (
        ("GIAB's HG003 truth", where["truth"] / TRUTH_VCF),
        ("HG003's benchmark BED", where["truth"] / BENCH_BED),
        ("rtg's SDF of the GRCh38 FASTA", where["sdf"]),
        ("rtg-tools", rtg_path(lab)),
    ):
        if not path.exists():
            raise _bootstrap_missing(f"{what} at {path}")
    _need_mask(lab)
    record = read_record(lab, "hg003")
    pair: tuple[Path, Path] | None
    if record is not None and not on_this_reference(lab, record):
        raise RunRefused(
            "HG003 is aligned on another reference; run"
            " `python -m genokeeper.calibrate.hg002 realign hg003`"
        )
    if record is not None and record.gvcf is not None:
        pair = None  # scored from its record
    elif any(where["downsampled"].glob("*q.gz")):
        pair = fastq_pair(where["downsampled"])
    elif (bases := recorded_bases(paths.reports)) is None:
        raise RunRefused(
            "missing input: no run of your FASTQ has recorded its bases, which"
            " HG003 is downsampled to; run `genokeeper run` first"
        )
    else:
        full = fastq_pair(where["fastq"])
        try:
            pair = downsample(*full, where["downsampled"], target_bases=bases)
        except ValueError as e:
            raise RunRefused(" ".join(str(e).split())) from e
        # Both downsampled files landed: nothing reads the full FASTQ again.
        for fastq in full:
            fastq.unlink()
    grading = grading if grading is not None else load_grading()
    strata = _strata(lab, grading)
    if pair is None:
        gvcf = read_called(lab, "hg003")
    else:
        run = run_sample(
            paths,
            name=hg003_dir(lab).name,
            fastq1=pair[0],
            fastq2=pair[1],
            sex=HG003_SEX,
            grading=grading,
            threads=threads,
        )
        gvcf = run.gvcf
        write_record(
            lab,
            "hg003",
            Record(
                cram=str(run.cram),
                key=run.key,
                metrics=run.metrics,
                gvcf=str(run.gvcf),
                reference_sha256=reference_sha256(lab),
                accepted=True,
                snp_f1=None,
                replaced=None,
            ),
        )
    versions = calibration_versions(
        lab, tools=genome_tool_versions(), identities=hg003_dir(lab) / IDENTITIES
    )
    # HG003 is scored on the masked reference: its benchmark loses the masked blocks and
    # their correct copies, where the truth, made on the unmasked reference, can lack real
    # variants (BR-030).
    bench = without_false_duplications(
        lab, where["truth"] / BENCH_BED, hg003_dir(lab) / "benchmark.masked.bed"
    )
    by_chrom, f1 = count_cells(
        gvcf,
        where["truth"] / TRUTH_VCF,
        bench,
        sdf=where["sdf"],
        strata=strata,
        rtg=rtg_path(lab),
        threads=threads,
    )
    table = CalibrationTable.from_counts(
        merged(by_chrom), versions, f1=f1, by_chrom=by_chrom
    )

    def write(staging: Path) -> None:
        table.write(staging)
        (staging / REPORT).write_text(report(table, by_chrom) + "\n")

    write_table(lab, "small_variants", write, versions)
    return table, by_chrom


def _need_mask(lab: Path) -> None:
    mask = masked_bed(lab)
    if not mask.is_file():
        raise RunRefused(
            f"missing input: the mask's BED at {mask}; run scripts/bootstrap.sh"
            " --with-masked-reference"
        )


def _strata(lab: Path, grading: Grading) -> StrataIndex:
    """The grading's GIAB strata, the difficult regions; `RunRefused` naming the bootstrap
    when one is absent."""
    giab = giab38_dir(lab)
    beds = {flag: giab / name for flag, name in grading.strata.items()}
    lacking = sorted(str(p) for p in beds.values() if not p.is_file())
    if lacking:
        raise RunRefused(
            f"missing input: GIAB GRCh38 BEDs {', '.join(lacking)};"
            " run scripts/bootstrap.sh"
        )
    return StrataIndex.from_beds(beds)


def build_hg002_cells(lab: Path, threads: int) -> None:
    """HG002's cells (`genokeeper.calibrate.hg002_cells`): HG002 PCR-free's gVCF scored in
    v5.0q's benchmark where HG003's truth does not reach, written whole to `cells_dir`
    with the versions of this machine's pipeline, and its report written. Without a
    called, current PCR-free record or v5.0q's truth it raises `RunRefused` "HG002 cells
    not measured: <what is missing>; run <the command that supplies it>", as on any other
    missing input, and writes nothing."""
    truth, bench = hg002_truth(lab)
    try:
        gvcf = read_called(lab, "pcrfree")
    except RunRefused as e:
        raise RunRefused(f"HG002 cells not measured: {' '.join(str(e).split())}") from e
    absent = [p for p in (truth, Path(f"{truth}.tbi"), bench) if not p.is_file()]
    if absent:
        raise RunRefused(
            f"HG002 cells not measured: no {absent[0]}; run scripts/bootstrap.sh"
            " --with-hg002"
        )
    where = layout(lab)
    for what, path in (
        ("HG003's benchmark BED", where["truth"] / BENCH_BED),
        ("rtg's SDF of the GRCh38 FASTA", where["sdf"]),
        ("rtg-tools", rtg_path(lab)),
    ):
        if not path.exists():
            raise _bootstrap_missing(f"{what} at {path}")
    _need_mask(lab)
    strata = _strata(lab, load_grading())
    versions = hg002_cells_versions(lab)
    parts, overlap = build_cells(
        lab,
        gvcf=gvcf,
        sdf=where["sdf"],
        strata=strata,
        versions={key: versions[key] for key in VERSION_KEYS},
        rtg=rtg_path(lab),
        threads=threads,
    )
    said = cells_report(parts, overlap)

    def write(staging: Path) -> None:
        write_cells(staging, parts, lab_regions(lab))
        (staging / REPORT).write_text(said + "\n")

    write_table(lab, "hg002_cells", write, versions)
    print(f"HG002 cells: {cells_dir(lab)}")


def _lower(*levels: Level | None) -> Level:
    """The lowest of `levels` (None ignored)."""
    order = list(Level)  # high, medium, low
    return max((lv for lv in levels if lv is not None), key=order.index)


def build_stage3(
    paths: Paths, *, grading: Grading | None = None, threads: int | None = None
) -> Stage3Table:
    """Measure the trio (`RECORDED`, each aligned by `python -m genokeeper.calibrate.hg002
    align`) through the run's stage-3 steps, each under `library_dir(lab, key)` and its run
    lock with its sex (`library_sex`): HG002's two libraries (`pcrfree`, `pcr`) through
    `stage3_steps`, HG003 and HG004 through `repeat_steps` alone. Then 3a's `measure` on HG002's libraries
    against GIAB's v5.0q truth and the repeat rules on all four (`measure_repeat_rules`,
    kept under `measured[RULES]`); the withheld keys merged (`merge_withheld`), each locus
    capped at the lower of both, and the stage-3 table written with this machine's versions
    (`stage3_versions`, the function the run compares with). `grading` defaults to
    config/grading.toml (its [repeat_loci] and [exclusion]), `threads` to this machine's
    CPUs.

    Raises `RunRefused` on a missing input (naming the bootstrap, or for a sample not
    aligned the align command, before any step runs), Docker not running, or a tool failing
    on any sample; nothing is written then."""
    lab = paths.lab
    vcf, bench = hg002_truth(lab)
    fasta = grch38_fasta(lab)
    for what, path in (
        ("GIAB's HG002 v5.0q truth", vcf),
        ("its index", Path(f"{vcf}.tbi")),
        ("its benchmark BED", bench),
        ("the GRCh38 FASTA", fasta),
    ):
        if not path.is_file():
            raise RunRefused(
                f"missing input: {what} at {path}; run scripts/bootstrap.sh --with-hg002"
            )
    grading = grading if grading is not None else load_grading()
    threads = threads if threads is not None else (os.cpu_count() or 4)
    fasta_sha = sha256_cached(fasta, hg002_dir(lab) / IDENTITIES)
    try:
        versions = stage3_versions(lab, fasta_sha, grading)
    except FileNotFoundError as e:
        raise RunRefused(f"missing input: {' '.join(str(e).split())}") from e
    # Every sample's record first: a sample not aligned refuses before hours of steps.
    aligned = {key: read_aligned(lab, key) for key in RECORDED}
    results: dict[str, Stage3Results] = {}
    samples: dict[str, SampleRepeats] = {}
    for key, (cram, cram_key, metrics) in aligned.items():
        work = library_dir(lab, key)
        sex = library_sex(key)
        with run_lock(work):
            if key in LIBRARIES:
                found, _ = stage3_steps(
                    work,
                    cram=cram,
                    cram_key=cram_key,
                    fasta=fasta,
                    fasta_sha=fasta_sha,
                    sex=sex,
                    mean_depth=metrics["mean_depth"],
                    lab=lab,
                    grading=grading,
                    threads=threads,
                )
                results[key] = found
            else:
                found, _ = repeat_steps(
                    work,
                    cram=cram,
                    cram_key=cram_key,
                    fasta=fasta,
                    fasta_sha=fasta_sha,
                    sex=sex,
                    grading=grading,
                    mean_depth=metrics["mean_depth"],
                    threads=threads,
                )
            if found.failed:
                said = "; ".join(
                    f"{tool} failed: {why}" for tool, why in found.failed.items()
                )
                raise RunRefused(f"{label(key)}: {said}; no stage 3 table written")
        samples[key] = SampleRepeats(
            sex=sex, calls=found.repeats or {}, rfc1=found.rfc1
        )
    loci = grading.repeats.loci
    truth = {
        gene: truth_repeats(
            vcf,
            bench,
            region=region(locus),
            motif=locus.motif,
            copies=1 if locus.chrom in HAPLOID else 2,
        )
        for gene, locus in loci.items()
    }
    withheld, caps, measured = measure(results, truth, loci)
    rules_withheld, rules_caps, rules = measure_repeat_rules(
        samples,
        loci,
        min_spanning=grading.exclusion.min_spanning,
        min_allele_reads=grading.repeats.min_allele_reads,
    )
    table = Stage3Table(
        versions=versions,
        withheld=merge_withheld(withheld, rules_withheld),
        caps={
            gene: _lower(caps.get(gene), rules_caps.get(gene))
            for gene in caps.keys() | rules_caps.keys()
        },
        flank={library: found.flank for library, found in results.items()},
        measured=measured | {RULES: rules},
    )

    def write(staging: Path) -> None:
        write_stage3(staging, table)
        (staging / REPORT).write_text(stage3_report(table) + "\n")

    write_table(lab, "stage3", write, versions)
    return table


def judge_library(
    work: Path,
    calls: Sequence[Deletion],
    truth: Sequence[Deletion],
    s: Screen,
    bench: Path,
    fasta: Path,
    *,
    classify: Callable[[Deletion], CallClass],
    svtype: str = "DEL",
) -> tuple[list[tuple[CallClass, int, bool]], list[tuple[int, bool]]]:
    """One library's rare calls (not common in `s`'s gnomAD-SV records) against the truth's
    rare events of `svtype` inside the benchmark BED `bench`, by Truvari (`against_truth`)
    in the folder `work`, emptied first: each call's class (`classify`), length and whether
    the truth confirms it, and each truth event's length and whether a call found it (the
    caller bands them). The same test against the same records on both sides, so a call is
    never counted a miss of a common truth event it was never meant to report."""
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    ours = [(d, class_key(classify(d))) for d in rare(calls, s)]
    theirs = [(d, "truth") for d in rare(truth, s)]
    folder = against_truth(
        write_calls(ours, fasta, work / "calls.vcf.gz"),
        write_calls(theirs, fasta, work / "truth.vcf.gz"),
        bench,
        work,
    )
    judged_calls, judged_truth = judged(folder, svtype=svtype)
    return [(class_of(name), n, ok) for name, n, ok in judged_calls], judged_truth


def build_deletions(
    paths: Paths, *, grading: Grading | None = None, threads: int | None = None
) -> DeletionTable:
    """Measure HG002's two libraries (each aligned by `python -m genokeeper.calibrate.hg002
    align`) through the run's deletion steps, under `<lab>/calibration/hg002/<library>/` and
    its run lock: the union's rare calls (not common in gnomAD-SV), each with its class,
    against the truth's rare deletions (the same test against the same gnomAD-SV records, so
    a call is never counted a miss of a common truth deletion it was never meant to report)
    inside the v5.0q benchmark; then each band's terms and likelihoods (`measure`), its prior
    (`expected_true_carried`), and the table with this machine's versions
    (`deletion_versions`, the function the run compares with). No stop rule (design decision
    5). `grading` defaults to config/grading.toml, `threads` to this machine's CPUs.

    Raises `RunRefused` on a missing input (naming the bootstrap or the align command), Docker
    not running, or a step failing on either library; nothing is written then."""
    lab = paths.lab
    truth, bench, fasta = _sv_inputs(lab)
    grading = grading if grading is not None else load_grading()
    t = grading.deletions
    threads = threads if threads is not None else (os.cpu_count() or 4)
    fasta_sha = sha256_cached(fasta, hg002_dir(lab) / IDENTITIES)
    versions = _or_refused(lambda: deletion_versions(lab, fasta_sha, grading))
    # Truvari drops a truth deletion under its `--sizemin`, so none of those can reach a result.
    floor = int(AGAINST_TRUTH[AGAINST_TRUTH.index("--sizemin") + 1])
    truth_all = [d for d in read_deletions(truth, caller="truth") if d.length >= floor]
    calls: dict[str, list[tuple[CallClass, int, bool]]] = {}
    found: dict[str, list[tuple[int, bool]]] = {}
    s = None
    for library in LIBRARIES:
        cram, cram_key, _ = read_aligned(lab, library)
        work = library_dir(lab, library)
        with run_lock(work):
            results, _ = deletion_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                grading=grading,
                threads=threads,
                male=library_sex(library) == "male",
            )
            if results.called is None:
                raise _refused(library, results.failed, "deletion")
            s = load_screen(lab, grading, [*results.called, *truth_all], work)
            # The screen's `SHORT` cut first: DELLY's calls from half of `min_length` only
            # confirm Dysgu's, and are never counted as calls.
            shown = [d for d in results.called if not is_short(d, t)]
            calls[library], found[library] = judge_library(
                work / "deletions-truth",
                shown,
                truth_all,
                s,
                bench,
                fasta,
                classify=partial(classify, s=s),
            )
    assert s is not None
    priors = {band: expected_true_carried(lab, grading, s, band) for band in BANDS}
    benchmark = read_bed(bench)
    bands, likelihoods = measure_deletions(
        calls, found, geometry(s.targets, benchmark), priors, t
    )
    table = DeletionTable(
        versions=versions,
        bands=bands,
        likelihoods=likelihoods,
        measured={lib: summary(calls[lib], found[lib], t) for lib in LIBRARIES},
        extra={
            OUTSIDE_BENCHMARK: outside_benchmark(
                s.targets, dict.fromkeys(BANDS, benchmark)
            )
        },
    )

    def write(staging: Path) -> None:
        write_deletions(staging, table)
        (staging / REPORT).write_text(deletions_report(table) + "\n")

    write_table(lab, "deletions", write, versions)
    return table


# --- slice 3d: the cnv, insertions and paralogs tables ---------------------------------


def _refused(library: str, failed: Mapping[str, str], table: str) -> RunRefused:
    """A step failing on `library` refuses `table`, naming the library and each tool."""
    said = "; ".join(f"{tool} failed: {why}" for tool, why in failed.items())
    return RunRefused(f"{label(library)}: {said}; no {table} table written")


def _sv_inputs(lab: Path) -> tuple[Path, Path, Path]:
    """v5.0q's SV truth, its benchmark BED and the FASTA, gnomAD-SV checked beside them;
    `RunRefused` naming the bootstrap for the first one missing (3b's and 3d's inputs)."""
    truth, bench = truth_files(lab)
    fasta = grch38_fasta(lab)
    for what, path, flag in (
        ("GIAB's HG002 v5.0q SV truth", truth, "--with-hg002"),
        ("its benchmark BED", bench, "--with-hg002"),
        ("the GRCh38 FASTA", fasta, ""),
        *(("gnomAD-SV", p, "--with-sv") for p in sv_files(lab).required()),
    ):
        if not path.is_file():
            command = f"scripts/bootstrap.sh {flag}".strip()
            raise RunRefused(f"missing input: {what} at {path}; run {command}")
    return truth, bench, fasta


def _or_refused[T](make: Callable[[], T]) -> T:
    """`make()`, a file it reads missing refused as a missing input (a table's versions,
    Parascopy's data)."""
    try:
        return make()
    except FileNotFoundError as e:
        raise RunRefused(f"missing input: {' '.join(str(e).split())}") from e


def _aligned(
    lab: Path, keys: Collection[str]
) -> dict[str, tuple[Path, str, dict[str, float]]]:
    """Each library's CRAM, its key and metrics, all read before any step runs
    (`read_aligned`: a library not aligned refuses, naming 3c's align command)."""
    return {key: read_aligned(lab, key) for key in keys}


def _dysgu_ins_dup(
    work: Path,
    library: str,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    table: str,
) -> tuple[list[Deletion], _Step]:
    """3b's deletion steps on the library (cached, as 3b's calibration left them) and
    `dysgu_ins_dup` on its Dysgu step: the callers' union and that step; a step failing
    refuses `table`, naming the library and the tool."""
    results, steps = deletion_steps(
        work,
        cram=cram,
        cram_key=cram_key,
        fasta=fasta,
        fasta_sha=fasta_sha,
        grading=grading,
        threads=threads,
        male=library_sex(library) == "male",
    )
    if results.called is None:
        raise _refused(library, results.failed, table)
    failed: dict[str, str] = {}
    step = _attempt(
        steps,
        failed,
        "dysgu_ins_dup",
        partial(dysgu_ins_dup_step, work, dysgu=steps["dysgu"], grading=grading),
    )
    if step is None:
        raise _refused(library, failed, table)
    return results.called, step


def _depth_events(
    work: Path,
    library: str,
    *,
    cram: Path,
    cram_key: str,
    fasta: Path,
    fasta_sha: str,
    grading: Grading,
    threads: int,
    masked: Path,
    junctions: _Step,
) -> tuple[dict[bool, list[Deletion]], list[Deletion]]:
    """The run's depth steps (`cnv_steps`) on `cram`: the events with CNVpytor merged (key
    True) and DELLY's alone (False, its cached step read again), and Dysgu's DUP records;
    a step failing refuses the `cnv` table, naming the library and the tool."""
    events: dict[bool, list[Deletion]] = {}
    junctions_found: list[Deletion] = []
    for cnvpytor in (True, False):
        found, _ = cnv_steps(
            work,
            cram=cram,
            cram_key=cram_key,
            fasta=fasta,
            fasta_sha=fasta_sha,
            grading=grading,
            threads=threads,
            male=library_sex(library) == "male",
            masked=masked,
            cnvpytor=cnvpytor,
            junctions=junctions,
        )
        if found.events is None:
            raise _refused(library, found.failed, "cnv")
        events[cnvpytor], junctions_found = found.events, found.junctions
    return events, junctions_found


type Calls = dict[str, list[tuple[CallClass, int, bool]]]


@dataclass(frozen=True)
class _Depth:
    """One library's depth events (CNVpytor merged: True; DELLY's alone: False), Dysgu's
    DUP junctions, 3b's calls of `min_length` or more, and on HG002 the spike-ins and the
    events found on its spiked copy."""

    events: dict[bool, list[Deletion]]
    junctions: list[Deletion]
    deletions: list[Deletion]
    spikes: list[Deletion]
    spiked: dict[bool, list[Deletion]]


def _with_sequence(vcf: Path) -> list[tuple[Deletion, str]]:
    """The truth's sequence-resolved insertions, each with its inserted sequence (the ALT
    without its anchor base), for `classify_gain_truth`."""
    alts = {(c[0], int(c[1]), c[2]): c[4] for c, _, _ in vcf_records(vcf)}
    found: list[tuple[Deletion, str]] = []
    for d in read_deletions(vcf, caller="truth", svtype="INS"):
        alt = alts[d.chrom, d.pos, d.callers["truth"].split(" ")[0]]
        if not alt.startswith("<"):
            found.append((d, alt[1:]))
    return found


def _depth_band(
    svtype: str, t: CnvThresholds
) -> Callable[[int, DeletionThresholds], str | None]:
    """A gnomAD-SV record's depth band by its length (`cnv_band`), for the priors."""

    def band(length: int, _: DeletionThresholds) -> str | None:
        if length < t.min_length:
            return None
        return cnv_band(Deletion("", 0, length, "N", "", {}, None, svtype=svtype), t)

    return band


def _contigs(fasta: Path) -> dict[str, int]:
    """The FASTA's primary chromosomes (`PRIMARY`) and their lengths, from its `.fai`."""
    rows = [line.split("\t") for line in Path(f"{fasta}.fai").read_text().splitlines()]
    return {row[0]: int(row[1]) for row in rows if row and row[0] in PRIMARY}


def _class_band(c: CallClass) -> str:
    return c.band


def _svtype(band: str) -> str:
    """A depth band's event type: a loss band's deletions, a gain band's duplications."""
    return "DEL" if band.startswith("loss") else "DUP"


def _at(d: Deletion) -> tuple[str, int, int, str]:
    return d.chrom, d.pos, d.end, d.svtype


def _class_in(
    classed: Mapping[tuple[str, int, int, str], CallClass], d: Deletion
) -> CallClass:
    """The class the screen gave the event (one per position and type)."""
    return classed[_at(d)]


def build_cnv(lab: Path, threads: int) -> None:
    """The `cnv` table (the plan's Decisions): per library of `TRIO` (`pcr` read as its CRAM
    alone), under its `library_dir` and run lock, with its own sex, 3b's deletion steps
    (cached), `dysgu_ins_dup` and the run's depth steps; on HG002's two libraries also a
    spiked copy of the CRAM (`write_spiked`, its depth steps alone, the copy deleted after).
    Every library's events go through the run's depth screen (`screen_cnv`) over whole
    chromosomes, so the calls measured are those it would class. HG002's rare calls of the
    measured bands are judged in the benchmark (`judge_library`: losses against v5.0q's
    deletions of `[cnv] min_length` or more, gains against `classify_gain_truth`), gains
    outside it by the trio (`trio_judged`); recall per band is the spike-ins found
    (`spike_recall`); priors gnomAD-SV's rare deletions and duplications over the targets;
    a loss band's geometry the benchmark's, a gain band's the whole genome's. The large
    events of the four libraries with CNVpytor are pooled per genome (`pool_large`,
    `large_withheld`); when `cnvpytor_kept` drops CNVpytor the bands are measured once more
    on DELLY's events alone. Raises `RunRefused` on a missing input (ClinGen's region lists
    and dosage list among them), a library not aligned, or a step failing on any library;
    nothing is written then."""
    truth, bench, fasta = _sv_inputs(lab)
    _need_mask(lab)
    regions, recurrent, record = _clingen_regions(lab)
    if "reason" in record:
        raise RunRefused(f"cnv table not measured: {record['reason']}")
    listed = _newest_dated(lab / "clingen-dosage", [DOSAGE_FILE])
    if listed is None:
        raise RunRefused(
            f"missing input: no ClinGen dosage sensitivity list {lab}/clingen-dosage/"
            f"<YYYYMMDD>/{DOSAGE_FILE}; run scripts/bootstrap.sh"
        )
    dosage = load_dosage(listed / DOSAGE_FILE)
    grading = load_grading()
    t, d_t = grading.cnv, grading.deletions
    fasta_sha = _fasta_sha(lab)
    versions = _or_refused(lambda: cnv_versions(lab, fasta_sha, grading))
    _, models = _or_refused(lambda: parascopy_data(lab))
    aligned = _aligned(lab, TRIO)
    masked = masked_bed(lab)
    masked_intervals = read_bed(masked)
    owned = owned_regions(
        load_spans(), [models / f"{locus.model}.gz" for locus in LOCI.values()]
    )
    giab = giab38_dir(lab)
    mask = [read_bed(giab / grading.strata[n]) for n in ("segdup", "low_mappability")]
    losses, gains = pick_spikes(
        targets(load_exons(), d_t.exon_margin),
        [*mask[0], *mask[1], *masked_intervals],
        owned,
    )
    data: dict[str, _Depth] = {}
    common = {"fasta": fasta, "fasta_sha": fasta_sha, "grading": grading}
    for key, (cram, cram_key, _) in aligned.items():
        work = library_dir(lab, key)
        with run_lock(work):
            called, junctions = _dysgu_ins_dup(
                work,
                key,
                cram=cram,
                cram_key=cram_key,
                threads=threads,
                table="cnv",
                **common,
            )
            events, dups = _depth_events(
                work,
                key,
                cram=cram,
                cram_key=cram_key,
                threads=threads,
                masked=masked,
                junctions=junctions,
                **common,
            )
            spiked: dict[bool, list[Deletion]] = {}
            if TRIO[key] == "HG002":
                out = work / "spikein" / "spiked.cram"
                spikes = json.dumps([cram_key, [asdict(d) for d in [*losses, *gains]]])
                try:
                    write_spiked(cram, fasta, losses, gains, out, threads=threads)
                    spiked, _ = _depth_events(
                        out.parent,
                        key,
                        cram=out,
                        cram_key=hashlib.sha256(spikes.encode()).hexdigest(),
                        threads=threads,
                        masked=masked,
                        junctions=junctions,
                        **common,
                    )
                finally:
                    out.unlink(missing_ok=True)
                    Path(f"{out}.crai").unlink(missing_ok=True)
        data[key] = _Depth(
            events=events,
            junctions=dups,
            deletions=[d for d in called if not is_short(d, d_t)],
            spikes=[*losses, *gains] if TRIO[key] == "HG002" else [],
            spiked=spiked,
        )
    truth_losses = [
        d for d in read_deletions(truth, caller="truth") if d.length >= t.min_length
    ]
    gain_truth = classify_gain_truth(_with_sequence(truth), fasta)
    everything = [
        *truth_losses,
        *gain_truth,
        *(
            d
            for lib in data.values()
            for found in (*lib.events.values(), *lib.spiked.values(), lib.deletions)
            for d in found
        ),
    ]
    s = load_screen(lab, grading, everything, hg002_dir(lab), callers=DEPTH_CALLERS)
    # A gain is common by gnomAD-SV's duplications, a loss by its deletions, as in the run.
    dups_common = read_population(
        sv_files(lab).gnomad,
        [(d.chrom, d.pos, d.end) for d in everything],
        d_t.min_gnomad_an,
        hg002_dir(lab),
        svtype="DUP",
    )
    s = replace(s, population=[*s.population, *dups_common])
    contigs = _contigs(fasta)
    # The screen over whole chromosomes: every event it would class, not only those over a
    # target exon, is a call measured (the run's screen past its `NO_EXON` drop).
    whole = replace(s, targets={c: [Target(c, 0, n, ())] for c, n in contigs.items()})
    bench_at = _Index((c, a, b, (a, b)) for c, a, b in read_bed(bench))

    def in_bench(d: Deletion) -> bool:
        return any(
            a <= d.pos and d.end <= b for a, b in bench_at.hits(d.chrom, d.pos, d.end)
        )

    def screened(lib: _Depth, events: Sequence[Deletion]) -> CnvScreen:
        return screen_cnv(
            events,
            rare(lib.deletions, s),
            lib.junctions,
            whole,
            t,
            masked=masked_intervals,
            owned=owned,
            spans=gene_spans(whole.targets),
            dosage=dosage,
            regions=regions,
            recurrent=recurrent,
            named=(),
        )

    def judge(
        cnvpytor: bool,
    ) -> tuple[Calls, dict[str, list[tuple[str, bool]]], dict[str, list[LargeEvent]]]:
        calls: Calls = {}
        recall: dict[str, list[tuple[str, bool]]] = {}
        large: dict[str, list[LargeEvent]] = {}
        parents = [
            lib.events[cnvpytor] for k, lib in data.items() if TRIO[k] != "HG002"
        ]
        for key, lib in data.items():
            found = screened(lib, lib.events[cnvpytor])
            large[key] = found.large
            if TRIO[key] != "HG002":
                continue
            shown = [x for x in found.shown if x.call_class.band in CNV_BANDS]
            class_of_event = partial(
                _class_in, {_at(x.deletion): x.call_class for x in shown}
            )
            events = [x.deletion for x in shown]
            gains = [d for d in events if d.svtype == "DUP"]
            out = (
                library_dir(lab, key) / "cnv-truth" / ("both" if cnvpytor else "delly")
            )
            with run_lock(library_dir(lab, key)):
                lost, _ = judge_library(
                    out / "loss",
                    [d for d in events if d.svtype == "DEL"],
                    truth_losses,
                    s,
                    bench,
                    fasta,
                    classify=class_of_event,
                )
                gained, _ = judge_library(
                    out / "gain",
                    [d for d in gains if in_bench(d)],
                    gain_truth,
                    s,
                    bench,
                    fasta,
                    classify=class_of_event,
                    svtype="DUP",
                )
            # The screen dropped the common gains already.
            trio = trio_judged(
                [d for d in gains if not in_bench(d)],
                parents,
                min_overlap=MATCH_OVERLAP,
                common=lambda _: False,
            )
            calls[key] = [
                *lost,
                *gained,
                *((class_of_event(d), d.length, ok) for d, ok in trio),
            ]
            spiked = screened(lib, lib.spiked[cnvpytor])
            recall[key] = spike_recall(
                [x.deletion for x in spiked.shown],
                lib.spikes,
                t,
                min_overlap=MATCH_OVERLAP,
            )
        return calls, recall, large

    calls, recall, large = judge(True)
    kept, decision = cnvpytor_kept(calls)
    pooled = pool_large(large, min_overlap=MATCH_OVERLAP)
    if not kept:
        calls, recall, _ = judge(False)
    benchmark, genome = read_bed(bench), [(c, 0, n) for c, n in contigs.items()]
    measured_in = {b: benchmark if b.startswith("loss") else genome for b in CNV_BANDS}
    geometries = {b: geometry(s.targets, region) for b, region in measured_in.items()}
    priors = {
        b: expected_true_carried(
            lab, grading, s, b, svtype=_svtype(b), band_of=_depth_band(_svtype(b), t)
        )
        for b in CNV_BANDS
    }
    bands, likelihoods = measure_classes(
        calls, recall, geometries, priors, CNV_FEATURES, band=_class_band
    )
    table = ClassTable(
        versions=versions,
        bands=bands,
        likelihoods=likelihoods,
        measured={
            lib: summarise_classes(
                calls[lib],
                recall[lib],
                CNV_FEATURES,
                CNV_BANDS,
                band=_class_band,
                key=class_key,
            )
            for lib in calls
        },
        features=CNV_FEATURES,
        extra={
            "cnvpytor": decision,
            "large": {"calls": pooled, "withheld": large_withheld(pooled)},
            OUTSIDE_BENCHMARK: outside_benchmark(s.targets, measured_in),
        },
    )
    path = write_class_table(
        lab, "cnv", table, report=slice3d_report(table, None, None)
    )
    print(f"cnv calibration table: {path}, CNVpytor {decision['decision']}")


def _by_position(
    calls: Sequence[Deletion],
    truth: Sequence[Deletion],
    region: Sequence[Interval],
    *,
    distance: int,
    classify: Callable[[Deletion], CallClass],
) -> tuple[list[tuple[CallClass, int, bool]], list[tuple[str | None, bool]]]:
    """Scramble's mobile elements judged by position (C20: their length is unknown), those
    in `region` (where Scramble looked, inside the benchmark): each true when a truth
    insertion lies within `distance`; each truth insertion there found when a call does."""
    at = _Index((c, a, b, None) for c, a, b in region)
    calls = [d for d in calls if at.covers(d.chrom, d.pos)]
    truth = [d for d in truth if at.covers(d.chrom, d.pos)]

    def near(d: Deletion, others: Sequence[Deletion]) -> bool:
        return any(
            o.chrom == d.chrom and abs(o.pos - d.pos) <= distance for o in others
        )

    return (
        [(classify(d), d.length, near(d, truth)) for d in calls],
        [("mei", near(d, calls)) for d in truth],
    )


def build_insertions(lab: Path, threads: int) -> None:
    """The `insertions` table: per HG002 library (`pcr` read as its CRAM alone), under its
    `library_dir` and run lock, 3b's deletion steps (cached), `dysgu_ins_dup` and the run's
    insertion steps; the union's rare calls of `[insertions] min_length` or more judged
    against v5.0q's rare insertions in the benchmark (`judge_library`, bands `short` and
    `long` by length), Scramble's lone mobile elements (band `mei`) by position inside the
    benchmark where Scramble looked (`_by_position`). Priors: gnomAD-SV's rare insertions
    over the targets, of the band's length (`mei`: any). The `mei` band's geometry is
    Scramble's regions within the benchmark, the others' the benchmark. Raises `RunRefused`
    on a missing input, a library not aligned, or a step failing on either library;
    nothing is written then."""
    truth, bench, fasta = _sv_inputs(lab)
    grading = load_grading()
    t = grading.insertions
    fasta_sha = _fasta_sha(lab)
    versions = _or_refused(lambda: insertion_versions(lab, fasta_sha, grading))
    aligned = _aligned(lab, LIBRARIES)
    floor = int(AGAINST_TRUTH[AGAINST_TRUTH.index("--sizemin") + 1])
    truth_all = [
        d
        for d in read_deletions(truth, caller="truth", svtype="INS")
        if d.length >= floor
    ]
    near = targets(load_exons(), SCRAMBLE_MARGIN)
    looked = merge((x.chrom, x.start, x.end) for ts in near.values() for x in ts)
    benchmark = read_bed(bench)
    scanned = subtract(looked, subtract(looked, benchmark))
    calls: dict[str, list[tuple[CallClass, int, bool]]] = {}
    found: dict[str, list[tuple[str | None, bool]]] = {}
    s = None
    for library, (cram, cram_key, _) in aligned.items():
        work = library_dir(lab, library)
        with run_lock(work):
            _, dysgu_ins = _dysgu_ins_dup(
                work,
                library,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                grading=grading,
                threads=threads,
                table="insertions",
            )
            results, _ = insertion_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                grading=grading,
                threads=threads,
                dysgu_ins=dysgu_ins,
                regions=scramble_regions(near, work / "scramble.regions.bed"),
            )
            if results.called is None:
                raise _refused(library, results.failed, "insertions")
            s = load_screen(
                lab,
                grading,
                [*results.called, *truth_all],
                work,
                svtype="INS",
                callers=INSERTION_CALLERS,
                thresholds=t,
            )
            of = partial(classify, s=s)
            classed = [(d, of(d).band == "mei") for d in results.called]
            meis = [d for d, mei in classed if mei]
            union = [d for d, mei in classed if not mei and d.length >= t.min_length]
            judged_calls, judged_truth = judge_library(
                work / "insertions-truth",
                union,
                truth_all,
                s,
                bench,
                fasta,
                classify=of,
                svtype="INS",
            )
        mei_calls, mei_truth = _by_position(
            rare(meis, s),
            rare([d for d in truth_all if d.length >= t.min_length], s),
            scanned,
            distance=t.max_distance,
            classify=of,
        )
        calls[library] = [*judged_calls, *mei_calls]
        found[library] = [*((band_of(n, t), ok) for n, ok in judged_truth), *mei_truth]
    assert s is not None

    def by_length(n: int, _: DeletionThresholds) -> str | None:
        return band_of(n, t)

    def any_length(n: int, _: DeletionThresholds) -> str | None:
        return "mei" if n >= t.min_length else None

    priors = {
        b: expected_true_carried(
            lab,
            grading,
            s,
            b,
            svtype="INS",
            band_of=any_length if b == "mei" else by_length,
        )
        for b in INSERTION_BANDS
    }
    regions = {b: scanned if b == "mei" else benchmark for b in INSERTION_BANDS}
    geometries = {b: geometry(s.targets, region) for b, region in regions.items()}
    bands, likelihoods = measure_classes(
        calls, found, geometries, priors, INSERTION_FEATURES, band=_class_band
    )
    table = ClassTable(
        versions=versions,
        bands=bands,
        likelihoods=likelihoods,
        measured={
            lib: summarise_classes(
                calls[lib],
                found[lib],
                INSERTION_FEATURES,
                INSERTION_BANDS,
                band=_class_band,
                key=class_key,
            )
            for lib in LIBRARIES
        },
        features=INSERTION_FEATURES,
        extra={OUTSIDE_BENCHMARK: outside_benchmark(s.targets, regions)},
    )
    path = write_class_table(
        lab, "insertions", table, report=slice3d_report(None, table, None)
    )
    print(f"insertions calibration table: {path}")


def build_paralogs(lab: Path, threads: int) -> None:
    """The `paralogs` table: per library of `TRIO` (`pcr` read as its CRAM alone), under
    its `library_dir` and run lock, the shared Parascopy depth step (cached, as 3a's SMN
    calibration may have left it) and those of the run's paralog steps a verdict reads
    (Parascopy's loci on every library, Gauchian and KILDA on HG002's two); each locus's usable copy
    number (as the run reads it), Gauchian's GBA1 call and KILDA's KIV-2 count measured
    against HG002's truth and the trio (`paralog_measure`). Raises `RunRefused` on a
    missing input, a library not aligned, or a step failing on any library; nothing is
    written then."""
    fasta = grch38_fasta(lab)
    if not fasta.is_file():
        raise RunRefused(
            f"missing input: the GRCh38 FASTA at {fasta}; run scripts/bootstrap.sh"
        )
    grading = load_grading()
    fasta_sha = _fasta_sha(lab)
    versions = _or_refused(lambda: paralog_versions(lab, fasta_sha, grading))
    _or_refused(lambda: parascopy_data(lab))
    aligned = _aligned(lab, TRIO)
    copies: dict[str, dict[str, int | None]] = {}
    gauchian: dict[str, GauchianCall | None] = {}
    kiv2: dict[str, float | None] = {}
    for key, (cram, cram_key, metrics) in aligned.items():
        work = library_dir(lab, key)
        with run_lock(work):
            steps: dict[str, _Step] = {}
            failed: dict[str, str] = {}
            depth = _attempt(
                steps,
                failed,
                PARASCOPY_DEPTH,
                partial(
                    parascopy_depth_step,
                    work,
                    cram=cram,
                    cram_key=cram_key,
                    fasta=fasta,
                    fasta_sha=fasta_sha,
                    threads=threads,
                ),
            )
            if depth is None:
                raise _refused(key, failed, PARALOGS)
            results, _ = paralog_steps(
                work,
                cram=cram,
                cram_key=cram_key,
                fasta=fasta,
                fasta_sha=fasta_sha,
                lab=lab,
                threads=threads,
                mean_depth=metrics["mean_depth"],
                depth=depth,
                run=("parascopy_loci", "gauchian", "kilda")
                if TRIO[key] == "HG002"
                else ("parascopy_loci",),
            )
        if results.failed:
            raise _refused(key, results.failed, PARALOGS)
        copies[key] = {}
        for name in LOCI:
            usable = _usable(results.calls.get(name), name, grading.paralogs)
            copies[key][name] = usable if isinstance(usable, int) else None
        gauchian[key], kiv2[key] = results.gauchian, results.kiv2
    withheld, measured = paralog_measure(
        copies, load_truth(), gauchian, kiv2, load_kiv2_truth().total, grading.lpa
    )
    table = ParalogTable(versions=versions, withheld=withheld, measured=measured)
    path = write_paralog_table(lab, table, report=slice3d_report(None, None, table))
    print(f"paralogs calibration table: {path}, {len(withheld)} withheld")


# --- slice 3f ------------------------------------------------------------------------------

# ROH on the trio reads the EUR AF file: Ashkenazi Jews are not a panel population.
TRIO_POPULATION = "EUR"
MIXTURE_SEED = 7
SLICE3F_STEP_VERSION = "2"  # 2: ROH cut to the benchmark before bcftools roh
# The calibration's grade in the typing steps' calls key (`deepvariant_level`).
LEVEL_NAME = "deepvariant-gq"
_UNKNOWN = Called("./.", "unknown", RowType.NOCALL, None, None, None, None, ())
_DOSAGE = {"0/0": "0", "0/1": "1", "1/1": "2"}


def deepvariant_level(grading: Grading) -> Callable[[Allele, Called], Level | None]:
    """The calibration's grade of a trio library's call: low for an unknown call or one
    whose DeepVariant GQ gives an error above `p_low`, medium otherwise. This is the run's
    own grade where no calibration cell rates a call (`reliability.grade_called`'s GQ
    fallback): the trio has no small-variant table of its own to grade by, and only low
    versus not low decides what the typing steps write."""

    def level(_allele: Allele, called: Called) -> Level | None:
        if called.dosage == "unknown" or called.gq is None:
            return Level.LOW
        return Level.LOW if 10 ** (-called.gq / 10) > grading.p_low else Level.MEDIUM

    return level


def _truth_records(
    vcf: pysam.VariantFile, chrom: str, pos: int
) -> list[tuple[str, tuple[str, ...], tuple[int | None, ...]]]:
    if chrom not in vcf.header.contigs:
        return []
    return [
        (
            (rec.ref or "").upper(),
            tuple(a.upper() for a in rec.alts or ()),
            tuple(rec.samples[0].get("GT") or (None,)),
        )
        for rec in vcf.fetch(chrom, pos - 1, pos)
        if rec.pos == pos
    ]


def truth_calls(
    vcf: Path,
    bench: Path,
    alleles: Sequence[Allele],
    *,
    haploid: Callable[[str, int], bool],
) -> dict[Allele, Called]:
    """Each allele's genotype in a truth VCF (bgzipped, indexed) as a call: inside the
    benchmark BED `bench`, its copies (`roh.truth_gt`: none where the truth calls nothing
    at the position, unknown beside another called allele), one copy at most on a
    `haploid` locus; outside the benchmark, unknown."""
    intervals = StrataIndex.from_beds({"benchmark": bench})
    found: dict[Allele, Called] = {}
    with pysam.VariantFile(str(vcf)) as handle:
        for allele in alleles:
            chrom, pos, ref, alt = allele
            if not intervals.flags_at(chrom, pos):
                found[allele] = _UNKNOWN
                continue
            gt = truth_gt(_truth_records(handle, chrom, pos), ref, alt)
            dosage = _DOSAGE.get(gt, "unknown")
            if haploid(chrom, pos):
                dosage = {"0": "0", "2": "1"}.get(dosage, "unknown")
                gt = {"0": "0", "1": "1"}.get(dosage, ".")
            found[allele] = (
                _UNKNOWN
                if dosage == "unknown"
                else Called(gt, dosage, RowType.PASS, 99, None, None, None, ())
            )
    return found


def truth_ta(vcf: Path, pos: int = UGT1A1_28[0][1]) -> str | None:
    """The truth's UGT1A1 TA-repeat genotype at chr2:`pos` (`slice3f.ta_genotype`):
    TA6/TA6 where it calls nothing, each called allele's repeats from its length against
    REF; None for a no-call or several records."""
    with pysam.VariantFile(str(vcf)) as handle:
        records = [r for r in _truth_records(handle, "chr2", pos) if any(r[2])]
    if not records:
        return ta_genotype((6, 6))
    if len(records) > 1 or None in records[0][2]:
        return None
    ref, alts, gt = records[0]
    return ta_genotype(
        [
            6 if i == 0 else 6 + (len(alts[i - 1]) - len(ref)) // 2
            for i in gt
            if i is not None
        ]
    )


def ours_ta(calls: Mapping[Allele, Called]) -> str | None:
    """Our TA-repeat genotype from the calls at PharmCAT's three UGT1A1*28 ALTs (TA5,
    TA7, TA8); None when one is unknown or they give more than two alleles."""
    copies: list[int] = []
    for allele in UGT1A1_28:
        called = calls.get(allele, _UNKNOWN)
        if called.dosage == "unknown":
            return None
        _, _, ref, alt = allele
        copies += [6 + (len(alt) - len(ref)) // 2] * int(called.dosage)
    if len(copies) > 2:
        return None
    return ta_genotype([*copies, *[6] * (2 - len(copies))])


def truth_haplotypes(
    vcf: Path, bench: Path, sites: Sequence[Allele]
) -> dict[Allele, int]:
    """Each site the phased truth holds heterozygous (`0|1` or `1|0`) inside its benchmark:
    the haplotype its ALT is on."""
    intervals = StrataIndex.from_beds({"benchmark": bench})
    found: dict[Allele, int] = {}
    with pysam.VariantFile(str(vcf)) as handle:
        for site in sites:
            chrom, pos, ref, alt = site
            if not intervals.flags_at(chrom, pos) or chrom not in handle.header.contigs:
                continue
            for rec in handle.fetch(chrom, pos - 1, pos):
                sample = rec.samples[0]
                gt = tuple(sample.get("GT") or ())
                alts = tuple(a.upper() for a in rec.alts or ())
                if (
                    rec.pos == pos
                    and sample.phased
                    and (rec.ref or "").upper() == ref
                    and alt in alts
                    and sorted(gt) == [0, alts.index(alt) + 1]
                ):
                    found[site] = gt.index(alts.index(alt) + 1)
    return found


def compare_phase(
    ours: Mapping[Allele, tuple[int, int]], truth: Mapping[Allele, int]
) -> tuple[int, int]:
    """Agreeing and disagreeing cis/trans judgements against the phased truth: each two
    sites next to each other in one of our phase sets, both phased in the truth."""
    agree = disagree = 0
    sets: dict[int, list[Allele]] = defaultdict(list)
    for site, (ps, _) in ours.items():
        if site in truth:
            sets[ps].append(site)
    for members in sets.values():
        members.sort(key=lambda a: (a[0], a[1]))
        for a, b in itertools.pairwise(members):
            same = (ours[a][1] == ours[b][1]) == (truth[a] == truth[b])
            agree, disagree = agree + same, disagree + (not same)
    return agree, disagree


def _samtools(*args: str) -> str:
    return run_tool(["samtools", *args])


def _chrm_reads(cram: Path, fasta: Path, out: Path) -> int:
    """`cram`'s chrM reads as the indexed BAM `out`; their count."""
    _samtools("view", "-b", "-T", str(fasta), "-o", str(out), str(cram), "chrM")
    _samtools("index", str(out))
    return int(_samtools("view", "-c", str(out)).strip())


def _mito_json(calls: Sequence[MitoCall], found: Haplocheck) -> str:
    return json.dumps(
        {"calls": [asdict(c) for c in calls], "haplocheck": asdict(found)}
    )


def _read_mito(path: Path) -> tuple[list[MitoCall], Haplocheck]:
    data = json.loads(path.read_text())
    return [MitoCall(**c) for c in data["calls"]], Haplocheck(**data["haplocheck"])


def _systems(
    calls: Mapping[str, SystemCall] | None, holes: Mapping[str, tuple[int, ...]]
) -> dict[str, SystemCall | None]:
    """Each packaged system's RBCeq2 call, None where a position it reads was a no-call (R10,
    and BR-029's `bloodgroups.DEPENDS`) or RBCeq2 gave nothing."""
    return {
        s: None
        if calls is None or bloodgroup_holes(holes, s)
        else calls.get(s, SystemCall(None, None))
        for s in bloodgroup_systems()
    }


@dataclass(frozen=True)
class _Library:
    """What one trio library gave the measures."""

    typing: TypingResults
    tags: dict[str, bool | None]
    rs2187668: dict[str, str | None]
    ours_ta: str | None
    truth_ta: str | None
    ta_inside: bool
    roh: dict[str, list[Run]]
    truth_blood: dict[str, SystemCall | None]
    mito: tuple[list[MitoCall], Haplocheck]


type _Record = tuple[Path, str, dict[str, float], Path]


# The errors a calibration step raises when its tool fails or Docker cannot start it.
_STEP_ERRORS = (DockerUnavailable, ToolFailed, ValueError, RuntimeError, OSError)


def _slice3f_refused(what: str, error: BaseException) -> RunRefused:
    """`error` as the refusal of the slice 3f table, naming `what` failed (a library's key
    or a step); a `RunRefused` as it is."""
    if isinstance(error, RunRefused):
        return error
    said = " ".join(str(error).split())
    name = lib_label(what) if what in SLICE3F_LIBRARIES else what
    return RunRefused(f"{name}: {said}; no slice 3f table written")


def _measure_library(
    lab: Path,
    key: str,
    record: _Record,
    truth: tuple[Path, Path],
    *,
    grading: Grading,
    threads: int,
    tools: Mapping[str, str],
    fasta_sha: str,
) -> _Library:
    """One trio library through the typing steps and the calibration's own steps, under
    its folder's run lock; `RunRefused` when a typing tool failed."""
    cram, cram_key, metrics, gvcf = record
    truth_vcf, bench = truth
    fasta = grch38_fasta(lab)
    male = library_sex(key) == "male"
    level = deepvariant_level(grading)
    af = typing_files(lab).af_file(TRIO_POPULATION)
    work = library_dir(lab, key)

    def haploid(chrom: str, pos: int) -> bool:
        return on_haploid_contig(chrom, pos, grading=grading, male=male)

    def build_roh(out: Path) -> None:
        bed = out / "benchmark.bed"
        if key == "pcrfree":
            shutil.copy(bench, bed)
        else:  # v4.2.1, which can lack real variants on the correct copies (BR-030)
            without_false_duplications(lab, bench, bed)
        sites = inside(steps["roh_sites"].path, bed, out / "ours.vcf.gz")
        given = {
            "ours": sites,
            "truth": truth_sites_vcf(truth_vcf, bed, af, out / "truth.vcf.gz"),
            "spiked": spike(sites, SPIKE, out / "spiked.vcf.gz"),
        }
        for name, vcf in given.items():
            run_roh(vcf, af, out / f"{name}.txt")

    def build_truth_blood(out: Path) -> None:
        where = rbceq2_positions(fasta)
        calls = truth_calls(truth_vcf, bench, where.alleles, haploid=haploid)
        written = write_pharmcat_vcf(
            calls,
            where.alleles,
            out / "rbceq2_sites.vcf",
            male=male,
            levels=dict.fromkeys(calls, Level.HIGH),
            grading=grading,
        )
        found = _systems(
            read_rbceq2(run_rbceq2(written, out)), uncovered(written, where)
        )
        document = {s: None if c is None else asdict(c) for s, c in found.items()}
        (out / "truth_bloodgroups.json").write_text(json.dumps(document))

    def build_mito(out: Path) -> None:
        calls = call_chrm(cram, fasta, out, threads=threads)
        found = run_haplocheck(out / VCF_FILE, out / "haplocheck")
        (out / "mito.json").write_text(_mito_json(calls, found))

    with run_lock(work):
        typing, steps = typing_steps(
            work,
            cram=cram,
            cram_key=cram_key,
            fasta=fasta,
            fasta_sha=fasta_sha,
            vcf=gvcf.with_name(f"{GENOME_SAMPLE}.vcf.gz"),
            gvcf=gvcf,
            calls_key=f"{calls_key(gvcf)}:{LEVEL_NAME}:{grading.p_low}",
            mean_depth=metrics["mean_depth"],
            population=TRIO_POPULATION,
            male=male,
            lab=lab,
            grading=grading,
            threads=threads,
            table=None,
            level=level,
            tools=tools,
        )
        if typing.failed:
            said = "; ".join(f"{t} failed: {why}" for t, why in typing.failed.items())
            raise RunRefused(f"{lib_label(key)}: {said}; no slice 3f table written")
        scratch = work / ".slice3f-scratch"
        try:
            ours = calls_at(
                gvcf,
                [*TAG_ALLELES, *UGT1A1_28],
                fasta=fasta,
                haploid=haploid,
                scratch=scratch,
            )
        except (ToolFailed, ValueError, OSError) as e:
            raise _slice3f_refused(key, e) from e
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        try:
            truth_shas = [sha256_file(truth_vcf), sha256_file(bench)]
            roh_step = _cached(
                work,
                "slice3f_roh",
                SLICE3F_STEP_VERSION,
                [steps["roh_sites"].key, IMAGES["truvari"], sha256_file(af)]
                + [*truth_shas, json.dumps(SPIKE)]
                + ([] if key == "pcrfree" else [correct_copies_sha256()]),
                "spiked.txt",
                build_roh,
            )
            blood_step = _cached(
                work,
                "slice3f_truth_rbceq2",
                SLICE3F_STEP_VERSION,
                [IMAGES["rbceq2"], *truth_shas, sha256_file(POSITIONS_FILE)]
                + [
                    fasta_sha,
                    str(male),
                    json.dumps(BLOODGROUP_DEPENDS, sort_keys=True),
                ],
                "truth_bloodgroups.json",
                build_truth_blood,
            )
            mito_step = _cached(
                work,
                "slice3f_mito",
                SLICE3F_STEP_VERSION,
                [cram_key, IMAGES["mutserve"], IMAGES["haplocheck"], fasta_sha],
                "mito.json",
                build_mito,
            )
        except _STEP_ERRORS as e:
            raise _slice3f_refused(key, e) from e
    try:
        theirs = truth_calls(truth_vcf, bench, TAG_ALLELES, haploid=haploid)
        truth_genotype = truth_ta(truth_vcf)
        ta_inside = bool(
            StrataIndex.from_beds({"b": bench}).flags_at(*UGT1A1_28[0][:2])
        )
    except (ValueError, OSError) as e:
        raise _slice3f_refused(key, e) from e
    tags: dict[str, bool | None] = {}
    for rsid, allele in zip(TAG_RSIDS, TAG_ALLELES, strict=True):
        mine, known = ours[allele], theirs[allele]
        graded = mine.dosage != "unknown" and level(allele, mine) != Level.LOW
        tags[rsid] = (
            None
            if known.dosage == "unknown" or not graded
            else mine.dosage == known.dosage
        )
    dq25 = TAG_ALLELES[TAG_RSIDS.index(DQ25_TAG)]
    known = theirs[dq25].dosage
    return _Library(
        typing=typing,
        tags=tags,
        rs2187668={
            "ours": ours[dq25].dosage,
            "truth": None if known == "unknown" else known,
        },
        ours_ta=ours_ta(ours),
        truth_ta=truth_genotype,
        ta_inside=ta_inside,
        roh={
            name: read_runs(roh_step.path.with_name(f"{name}.txt"))
            for name in ("ours", "truth", "spiked")
        },
        truth_blood={
            s: None if c is None else SystemCall(**c)
            for s, c in json.loads(blood_step.path.read_text()).items()
        },
        mito=_read_mito(mito_step.path),
    )


def _mixtures(
    lab: Path,
    records: Mapping[str, _Record],
    *,
    fasta_sha: str,
    tools: Mapping[str, str],
    threads: int,
) -> dict[float, tuple[list[MitoCall], Haplocheck]]:
    """HG003's chrM reads subsampled (`samtools view -s`) and merged into HG002 PCR-free's
    at each `GRID` level that HG003's reads allow, each keyed by `mixture_level` of the two
    parts' counts (`samtools view -c`), with mutserve and haplocheck on each; cached under
    HG002 PCR-free's folder."""
    fasta = grch38_fasta(lab)
    hg002, hg003 = records["pcrfree"], records["hg003"]

    def build(out: Path) -> None:
        n2 = _chrm_reads(hg002[0], fasta, out / "hg002.chrM.bam")
        n3 = _chrm_reads(hg003[0], fasta, out / "hg003.chrM.bam")
        found: list[dict[str, Any]] = []
        for target in GRID:
            fraction = target * n2 / ((1 - target) * n3) if n3 else 1.0
            if fraction >= 1:
                continue
            name = f"mix{target:.2f}"
            sub, mix = out / f"{name}.hg003.bam", out / f"{name}.bam"
            seed = f"{MIXTURE_SEED}{f'{fraction:.6f}'[1:]}"
            _samtools(
                "view", "-b", "-s", seed, "-o", str(sub), str(out / "hg003.chrM.bam")
            )
            taken = int(_samtools("view", "-c", str(sub)).strip())
            _samtools(
                "merge", "-f", "-o", str(mix), str(out / "hg002.chrM.bam"), str(sub)
            )
            _samtools("index", str(mix))
            level = mixture_level(n2, taken)
            calls = call_chrm(mix, fasta, out / name, threads=threads)
            flag = run_haplocheck(out / name / VCF_FILE, out / name / "haplocheck")
            found.append(
                {
                    "target": target,
                    "level": level,
                    "n_hg002": n2,
                    "n_hg003": taken,
                    "mito": json.loads(_mito_json(calls, flag)),
                }
            )
            for path in (sub, mix, Path(f"{mix}.bai")):
                path.unlink(missing_ok=True)
        for name in ("hg002", "hg003"):
            for path in (out / f"{name}.chrM.bam", out / f"{name}.chrM.bam.bai"):
                path.unlink(missing_ok=True)
        (out / "mixtures.json").write_text(json.dumps(found))

    work = library_dir(lab, "pcrfree")
    with run_lock(work):
        try:
            step = _cached(
                work,
                "slice3f_mixtures",
                SLICE3F_STEP_VERSION,
                [hg002[1], hg003[1], IMAGES["mutserve"], IMAGES["haplocheck"]]
                + [fasta_sha, tools["samtools"], str(MIXTURE_SEED), json.dumps(GRID)],
                "mixtures.json",
                build,
            )
        except _STEP_ERRORS as e:
            raise _slice3f_refused("chrM mixtures", e) from e
    return {
        float(m["level"]): (
            [MitoCall(**c) for c in m["mito"]["calls"]],
            Haplocheck(**m["mito"]["haplocheck"]),
        )
        for m in json.loads(step.path.read_text())
    }


def _phase(
    lab: Path, record: _Record, truth: tuple[Path, Path], *, fasta_sha: str
) -> tuple[int, int]:
    """WhatsHap on HG002 PCR-free at its heterozygous PharmCAT positions, chained through
    each PharmCAT gene's span (`phasing.phase_sites`, as the run phases them), against
    HG002's phased truth (`compare_phase`)."""
    cram, cram_key, _, gvcf = record
    fasta = grch38_fasta(lab)
    vcf = gvcf.with_name(f"{GENOME_SAMPLE}.vcf.gz")
    positions = positions_path(lab)

    def build(out: Path) -> None:
        alleles = pharmcat_positions(positions)
        sites = sorted(
            het_snvs(vcf, [(a[0], a[1], a[1]) for a in alleles]) & {*alleles}
        )
        regions = [(c, s, e) for _, c, s, e in gene_regions(positions)]
        phased = phase_sites(cram, fasta, vcf, sites, regions)
        rows = [[list(site), ps, hap] for site, (ps, hap) in phased.items()]
        (out / "phased.json").write_text(json.dumps(rows))

    work = library_dir(lab, "pcrfree")
    with run_lock(work):
        try:
            step = _cached(
                work,
                "slice3f_phase",
                SLICE3F_STEP_VERSION,
                [cram_key, calls_key(gvcf), sha256_file(positions)]
                + [IMAGES["whatshap"], fasta_sha],
                "phased.json",
                build,
            )
        except _STEP_ERRORS as e:
            raise _slice3f_refused("phase", e) from e
    phased: dict[Allele, tuple[int, int]] = {
        (s[0], int(s[1]), s[2], s[3]): (int(ps), int(hap))
        for s, ps, hap in json.loads(step.path.read_text())
    }
    return compare_phase(phased, truth_haplotypes(truth[0], truth[1], list(phased)))


def slice3f_table(
    measured: Mapping[str, _Library],
    mixtures: Mapping[float, tuple[list[MitoCall], Haplocheck]],
    phase: tuple[int, int],
    grading: Grading,
) -> Slice3fTable:
    """The measures over the trio's results (`genokeeper.calibrate.slice3f`)."""
    libs = list(SLICE3F_LIBRARIES)
    typing = {lib: measured[lib].typing for lib in libs}
    report: dict[str, Any] = {"population": f"{TRIO_POPULATION} AF file (ROH)"}
    hla_calls = {
        lib: (typing[lib].optitype or {}, typing[lib].t1k or {}) for lib in libs
    }
    withheld = measure_hla(
        hla_calls, load_hg002_truth(), {lib: measured[lib].tags for lib in libs}
    )
    report["hla"] = {
        lib: {"optitype": o, "t1k": t} for lib, (o, t) in hla_calls.items()
    }
    report["dq_tags"] = {lib: measured[lib].tags for lib in libs}
    report["rs2187668"] = {lib: measured[lib].rs2187668 for lib in libs}
    ugt1a1 = measure_ugt1a1(
        {lib: measured[lib].ours_ta for lib in libs},
        {lib: measured[lib].truth_ta for lib in libs},
        {lib: measured[lib].ta_inside for lib in libs},
    )
    report["ugt1a1"] = {
        lib: {
            "ours": measured[lib].ours_ta,
            "truth": measured[lib].truth_ta,
            "inside": measured[lib].ta_inside,
        }
        for lib in libs
    }
    report["phase"] = {"agree": phase[0], "disagree": phase[1]}
    # chrM: HG003's homoplasmies HG002 lacks are what each mixture must show.
    unmixed = {lib: measured[lib].mito for lib in libs}
    homoplasmy = grading.mito.min_homoplasmy
    hg002_at = {c.pos for c in unmixed["pcrfree"][0]}
    expected = frozenset(
        c.pos
        for c in unmixed["hg003"][0]
        if c.level >= homoplasmy and c.pos not in hg002_at
    )
    false_levels = [
        c.level for calls, _ in unmixed.values() for c in calls if c.level < homoplasmy
    ]
    floor = measure_chrm_floor(
        {level: calls for level, (calls, _) in mixtures.items()}, expected, false_levels
    )
    if floor is None:
        top = max(false_levels, default=0.0)
        withheld[CHRM_KEY] = (
            f"not covered: false heteroplasmy at {top:.2f} in calibration"
            if top >= GRID[-1]
            else "not covered: no HG003 homoplasmy HG002 lacks to measure in calibration"
            if not expected
            else "not covered: no chrM mixture called within its level in calibration"
        )
    else:
        missed = measure_haplocheck(
            {level: found for level, (_, found) in mixtures.items()},
            max(floor, grading.mito.heteroplasmy_floor),
        )
        if missed is not None:
            withheld[CHRM_KEY], floor = missed, None
    report["chrm"] = {
        "expected_positions": len(expected),
        "highest_false_heteroplasmy": max(false_levels, default=None),
        "mixtures": {
            f"{level:.4f}": {"step": grid_step(level), "flagged": found.contaminated}
            for level, (_, found) in mixtures.items()
        },
        "unmixed_flagged": sum(found.contaminated for _, found in unmixed.values()),
    }
    ratios: dict[str, float] = {}
    for lib in libs:
        call = typing[lib].hhv6
        ratios[lib] = math.inf if call is None or math.isnan(call.ratio) else call.ratio
    ceiling, why = measure_hhv6(ratios, grading.hhv6)
    background: float | None = ceiling
    if why is not None:
        withheld[HHV6_KEY], background = why, None
    report["hhv6_ratios"] = ratios
    roh = {lib: measured[lib].roh for lib in libs}
    if (
        why := measure_roh(
            {lib: r["ours"] for lib, r in roh.items()},
            {lib: r["truth"] for lib, r in roh.items()},
            {lib: r["spiked"] for lib, r in roh.items()},
            SPIKE,
            grading.roh,
        )
    ) is not None:
        withheld[ROH_KEY] = why
    chrom, start, end = SPIKE
    report["roh"] = {
        lib: {name: len(runs) for name, runs in r.items()}
        | {
            "truth_runs_in_spike": sum(
                t.chrom == chrom and t.start <= end and t.end >= start
                for t in r["truth"]
            )
        }
        for lib, r in roh.items()
    }
    shares: dict[str, float] = {}
    for lib in libs:
        share = typing[lib].archaic
        shares[lib] = math.inf if share is None else share.share
    if (why := measure_archaic(shares, PUBLISHED_ARCHAIC)) is not None:
        withheld[ARCHAIC_KEY] = why
    report["neanderthal_share"] = shares
    ours_blood = {
        lib: _systems(typing[lib].bloodgroups, typing[lib].bloodgroup_uncovered)
        for lib in libs
    }
    blood = measure_bloodgroups(
        ours_blood, {lib: measured[lib].truth_blood for lib in libs}
    )
    withheld |= blood
    report["bloodgroups_withheld"] = sorted(blood)
    report["ancestry"] = consistency_ancestry(
        {lib: typing[lib].admixture or {} for lib in libs}
    )
    return Slice3fTable(
        withheld=withheld,
        ugt1a1_trusted=ugt1a1,
        phase_trusted=measure_phase(*phase),
        chrm_floor=floor,
        hhv6_background=background,
        report=report,
    )


def build_slice3f(
    lab: Path, *, grading: Grading | None = None, threads: int | None = None
) -> Path:
    """Measure the trio's three called libraries (`slice3f.LIBRARIES`; HG002's PCR library
    is never read) through the run's typing steps (`pipeline.typing_steps`), each under its
    calibration folder and run lock, and write the slice 3f table with this machine's
    versions (`slice3f_versions`, the function the run compares with); its folder.

    Beside the typing steps: our tag-SNP and UGT1A1*28 calls against each truth (graded by
    `deepvariant_level`), RBCeq2 on each truth's genotypes at its positions (a no-call
    outside the benchmark, so not measured there), `bcftools roh` on our sites, the truth's
    and our sites spiked with `SPIKE`, all restricted to the benchmark (on the EUR AF file,
    `TRIO_POPULATION`), mutserve and haplocheck on each library and on HG003's chrM reads
    mixed into HG002 PCR-free's (`_mixtures`), and WhatsHap at PharmCAT's positions on HG002
    PCR-free against its phased truth (`_phase`). Every step is cached. `grading` defaults
    to config/grading.toml, `threads` to this machine's CPUs.

    Raises `RunRefused` before anything runs for a library not aligned or not called
    (`read_aligned`'s and `read_called`'s own refusal, naming the command), or a missing
    truth or typing file; after, for Docker not running or a tool failing. Nothing is
    written then."""
    grading = grading if grading is not None else load_grading()
    threads = threads if threads is not None else (os.cpu_count() or 4)
    records: dict[str, _Record] = {
        key: (*read_aligned(lab, key), read_called(lab, key))
        for key in SLICE3F_LIBRARIES
    }
    truth = slice3f_truths(lab)
    for key, (vcf, bench) in truth.items():
        flag = {"pcrfree": "--with-hg002", "hg003": "--with-calibration"}.get(
            key, "--with-hg004"
        )
        for what, path in (
            ("truth", vcf),
            ("truth's index", Path(f"{vcf}.tbi")),
            ("benchmark", bench),
        ):
            if not path.is_file():
                raise RunRefused(
                    f"missing input: {lib_label(key)}'s {what} at {path}; run"
                    f" scripts/bootstrap.sh {flag}"
                )
    missing = [p for p in typing_files(lab).required() if not p.is_file()]
    if missing:
        raise RunRefused(
            f"missing input: slice 3f's typing files ({missing[0]}, {len(missing)} in"
            " all); run scripts/bootstrap.sh --with-typing"
        )
    try:
        versions = slice3f_versions(lab, grading)
    except FileNotFoundError as e:
        raise RunRefused(
            f"missing input: {' '.join(str(e).split())}; run scripts/bootstrap.sh"
        ) from e
    tools = genome_tool_versions()
    common: dict[str, Any] = {
        "grading": grading,
        "threads": threads,
        "tools": tools,
        "fasta_sha": versions["fasta"],
    }
    measured = {
        key: _measure_library(lab, key, records[key], truth[key], **common)
        for key in SLICE3F_LIBRARIES
    }
    mixtures = _mixtures(
        lab, records, fasta_sha=versions["fasta"], tools=tools, threads=threads
    )
    phase = _phase(
        lab, records["pcrfree"], truth["pcrfree"], fasta_sha=versions["fasta"]
    )
    table = slice3f_table(measured, mixtures, phase, grading)
    return write_table(
        lab, "slice3f", lambda staging: write_slice3f(staging, table), versions
    )


def small_variants_versions(lab: Path) -> dict[str, str]:
    """The small-variant table's versions: this machine's pipeline's seven."""
    return calibration_versions(
        lab, tools=genome_tool_versions(), identities=hg003_dir(lab) / IDENTITIES
    )


def hg002_cells_versions(lab: Path) -> dict[str, str]:
    """HG002's cells' versions: this machine's pipeline's seven, the key of HG002
    PCR-free's alignment (the gVCF they score), and what cuts and scores their regions
    (`lab_regions`): `[par]`, HG003's benchmark BED and v5.0q's truth VCF and BED."""
    record = read_record(lab, "pcrfree")
    identities = library_dir(lab, "pcrfree") / IDENTITIES
    truth, v5 = hg002_truth(lab)
    return {
        **calibration_versions(
            lab, tools=genome_tool_versions(), identities=identities
        ),
        "pcrfree": "" if record is None else record.key,
        "par": json.dumps(load_grading().par, sort_keys=True),
        "hg003_benchmark": sha256_file(layout(lab)["truth"] / BENCH_BED),
        "v5_truth": sha256_cached(truth, identities),
        "v5_benchmark": sha256_file(v5),
    }


def _fasta_sha(lab: Path) -> str:
    return sha256_cached(grch38_fasta(lab), hg002_dir(lab) / IDENTITIES)


def _paths(lab: Path) -> Paths:
    return replace(load_paths(), lab=lab)


def _small_variants(lab: Path, threads: int) -> None:
    table, _ = build_table(_paths(lab), threads=threads)
    print(
        f"calibration table: {table_dir(lab)},"
        f" {len(table.measured())} genotype cells measured"
    )


def _stage3(lab: Path, threads: int) -> None:
    table = build_stage3(_paths(lab), threads=threads)
    print(
        f"stage 3 calibration table: {stage3_dir(lab)},"
        f" {len(table.withheld)} withheld, {len(table.caps)} capped"
    )


def _deletions(lab: Path, threads: int) -> None:
    table = build_deletions(_paths(lab), threads=threads)
    measured = ", ".join(
        f"{band} {terms.calls} calls" for band, terms in table.bands.items()
    )
    print(f"deletion calibration table: {deletions_dir(lab)}, {measured}")


def _slice3f(lab: Path, threads: int) -> None:
    folder = build_slice3f(lab, threads=threads)
    print(f"slice 3f calibration table: {folder}")


type Versions = Callable[[Path], dict[str, str]]
type Builder = Callable[[Path, int], None]
# Each registered table's versions and builder (`tables.TABLES` reaches them here).
ROWS: dict[str, tuple[Versions, Builder]] = {
    "small_variants": (small_variants_versions, _small_variants),
    "hg002_cells": (hg002_cells_versions, build_hg002_cells),
    "stage3": (
        lambda lab: stage3_versions(lab, _fasta_sha(lab), load_grading()),
        _stage3,
    ),
    "deletions": (
        lambda lab: deletion_versions(lab, _fasta_sha(lab), load_grading()),
        _deletions,
    ),
    "foundation": (foundation_versions, build_foundation),
    "cnv": (lambda lab: cnv_versions(lab, _fasta_sha(lab), load_grading()), build_cnv),
    "insertions": (
        lambda lab: insertion_versions(lab, _fasta_sha(lab), load_grading()),
        build_insertions,
    ),
    "paralogs": (
        lambda lab: paralog_versions(lab, _fasta_sha(lab), load_grading()),
        build_paralogs,
    ),
    "slice3f": (slice3f_versions, _slice3f),
}


def _calibrate_command(args: argparse.Namespace) -> int:
    threads = args.threads if args.threads is not None else (os.cpu_count() or 4)
    return calibrate(
        load_paths().lab, only=args.only, report=args.report, threads=threads
    )


def register(subparsers: Subparsers) -> None:
    """Add the `calibrate` sub-command."""
    g = subparsers.add_parser(
        "calibrate",
        help="rebuild every stale calibration table, each before those that need it",
    )
    g.add_argument(
        "--only",
        choices=list(TABLES),
        metavar="NAME",
        help=f"that table alone, refused while one it needs is stale ({', '.join(TABLES)})",
    )
    g.add_argument("--threads", type=int, help="default: this machine's CPUs")
    g.add_argument(
        "--report",
        action="store_true",
        help="print each built table's report (with --only, a current one's too)",
    )
    g.set_defaults(handler=_calibrate_command)
