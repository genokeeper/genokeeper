"""`genokeeper run` end to end on a mini GRCh38 lab, its containers faked: every claim the
release keeps has a call row and the accounting balances, each row type reaches the store,
a candidate pair is phased, a QC failure stops before calling, a run is atomic, a rerun
reuses every cached artifact, stage 3 runs only with a current calibration table, and every
refusal comes before any step."""

import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal, NoReturn

import duckdb
import pysam
import pytest
from fixtures.clinvar_mini import write_clinvar
from fixtures.contract import (
    CONTRACT,
    COVERAGE_CONTRACT,
    PGX_CONTRACT,
    SLICE3F_FINDING_CONTRACT,
    key_paths,
)
from fixtures.genome_mini import (
    CALIBRATED_ERRORS,
    CHRM_LENGTH,
    CLASSIFIED_LOF_POS,
    LARGE_REGION,
    MINI_QC,
    PHARMGKB_DATE,
    PREDICTED_POS,
    STAGE3_WITHHELD,
    FakeDeepVariant,
    FakeDeletions,
    FakeHaplocheck,
    FakeHaplogrep,
    FakeMutserve,
    FakeNextflow,
    FakePgx,
    FakeSlice3d,
    FakeStage3,
    FakeTyping,
    FakeVep,
    FakeVerifyBamID,
    FakeWhatsHap,
    GenomeLab,
    write_calibration,
    write_cnv,
    write_deletions,
    write_delivered,
    write_lab,
    write_paralogs,
    write_pgx,
    write_prs,
    write_slice3d,
    write_slice3f,
    write_stage3,
    write_sv,
)
from fixtures.genome_mini import bases as genome_bases
from fixtures.genome_mini import block as genome_block
from fixtures.genome_mini import variant as genome_variant
from fixtures.grch38_mini import CONTIGS as GRCH38_CONTIGS
from fixtures.grch38_mini import READ
from fixtures.gwas_mini import EFO, Assoc, write_gwas
from fixtures.pgsc_mini import PINNED
from fixtures.pharmgkb_mini import Ann, write_annotations

from genokeeper import (
    checksum,
    cli,
    haplocheck,
    haplogroup,
    insertions,
    mito,
    paralogs,
    phasing,
    pipeline,
    qc,
    refresh,
    smn,
    tools,
    yhaplogroup,
)
from genokeeper import coverage as gene_coverage
from genokeeper.annotate import vep
from genokeeper.calibrate import (
    BENCH_BED,
    TRUTH_VCF,
    hg003_dir,
    layout,
    rtg_path,
    table_dir,
    tables,
)
from genokeeper.calibrate import build as calibrate_build
from genokeeper.calibrate import foundation as calibrate_foundation
from genokeeper.calibrate.deletions import DeletionTable, deletions_dir
from genokeeper.calibrate.hg002_cells import (
    HAPLOID,
    NOT_MEASURED,
    ONLY,
    SHARED,
    cells_dir,
    open_cells,
    write_cells,
)
from genokeeper.calibrate.records import Record, reference_sha256, write_record
from genokeeper.calibrate.repeat_rules import X_RULE_DID_NOT_HOLD
from genokeeper.calibrate.stage3 import (
    FALSE_CALL,
    YLEAF_KEY,
    YLEAF_WITHHELD,
    stage3_dir,
    two_x_key,
)
from genokeeper.calibrate.table import (
    VERSION_KEYS,
    CalibrationTable,
    CellKey,
    Counts,
    jeffreys,
)
from genokeeper.calibrate.tables import NOT_CALIBRATED
from genokeeper.calibrate.tables import write_table as write_calibration_table
from genokeeper.call import deepvariant
from genokeeper.call.deepvariant import (
    CALLING_FLAGS,
    CHECKPOINT_FILES,
    checkpoint_path,
    deepvariant_args,
)
from genokeeper.carrier import genotype
from genokeeper.checksum import sha256_cached, sha256_file
from genokeeper.cnv import callers as cnv_callers
from genokeeper.cnv.screen import LARGE_FILTERED
from genokeeper.config import MASKED_REFERENCE, Grading, QcThresholds, load_grading
from genokeeper.containers import IMAGES, DockerUnavailable
from genokeeper.deletions import callers as deletion_callers
from genokeeper.deletions import truvari as deletion_truvari
from genokeeper.deletions.exons import GBA1_HOMOLOGY, Exon, load_spans
from genokeeper.deletions.files import truth_files
from genokeeper.deletions.screen import (
    COMMON,
    NO_EXON,
    SAME_POSITION,
    SHORT,
    CallClass,
    warnings,
)
from genokeeper.evidence.clingen import DOSAGE_FILE
from genokeeper.evidence.clinvar import INHERITANCE_REASON
from genokeeper.evidence.frequency import FafTable, faf_folder
from genokeeper.evidence.gencc import GENCC_FILE
from genokeeper.evidence.gwas import DOWNLOAD as GWAS_DOWNLOAD
from genokeeper.evidence.pgx import CPIC_FILES, GENE_UNRESOLVED
from genokeeper.evidence.pharmgkb import DOWNLOADS as PHARMGKB_DOWNLOADS
from genokeeper.genelists import DPYD_NOTE
from genokeeper.genotype import clonal
from genokeeper.genotype.clonal import Span
from genokeeper.genotype.reliability import AF_CLONAL
from genokeeper.history import dates, snapshots
from genokeeper.identity import analyser_identity
from genokeeper.paralogs import LOCI
from genokeeper.pgx import cyp2d6, hla, pharmcat
from genokeeper.pgx.hla import t1k_index
from genokeeper.pgx.hla_typing import CONFIRM_DQ
from genokeeper.pgx.pharmcat import PHARMCAT_VERSION, positions_path
from genokeeper.phasing import IN_TRANS_REASON
from genokeeper.pipeline import STORE_NAME, RunRefused
from genokeeper.pipeline import pgx_inputs as read_pgx_inputs
from genokeeper.prs import pgsc
from genokeeper.prs.scores import load_not_covered, load_scores
from genokeeper.query import Filters, query, show
from genokeeper.repeats import expansionhunter
from genokeeper.repeats.expansionhunter import EH_VERSION
from genokeeper.repeats.loci import Rule
from genokeeper.repeats.not_covered import NOT_COVERED
from genokeeper.sample_identity import IDENTITY_TOOL
from genokeeper.sources import panel as panel_source
from genokeeper.sources import prs as prs_source
from genokeeper.sources.bloodgroups import HPA_REASON, MNS_REASON, RH_REASON
from genokeeper.sources.deletions import CONSEQUENCE_WARNING
from genokeeper.sources.dosage import STATEMENT_ID as DOSAGE_STATEMENT
from genokeeper.sources.hhv6 import BELOW as HHV6_BELOW
from genokeeper.sources.hhv6 import BETWEEN as HHV6_BETWEEN
from genokeeper.sources.hla import KIR_REASON
from genokeeper.sources.mitomap import load_mitomap
from genokeeper.sources.paralogs import GBA1_CONTRADICTED
from genokeeper.sources.paralogs import NOT_COVERED as PARALOG_NOT_COVERED
from genokeeper.sources.repeats import (
    CFTR_UNREAD,
    NORMAL_TWO_X,
    X_DOSAGE,
    X_NOT_CALIBRATED,
)
from genokeeper.sources.smn import TWO_COPIES_REASON
from genokeeper.sources.traits import TRAIT_ALLELES
from genokeeper.store.db import TABLES, Level, Validity, connect
from genokeeper.yhaplogroup import YLEAF_VERSION

pytestmark = pytest.mark.tools

TWO_STARS = "criteria_provided,_multiple_submitters,_no_conflicts"
CLINVAR_RELEASE = "20261012"  # the ##fileDate the miniature ClinVar files carry


def _alt(ref: str) -> str:
    return "A" if ref != "A" else "G"


def _row(
    con: duckdb.DuckDBPyConnection, run_id: str, variation_id: int
) -> dict[str, Any]:
    """The statement `clinvar:<id>` of the run with its call's columns."""
    cursor = con.execute(
        "SELECT s.kind, s.topic, s.reason, s.validity, c.gt, c.dosage, c.level,"
        " c.row_type, c.warnings, c.flags, c.ad, c.p_wrong, c.p_wrong_lo, c.p_wrong_hi,"
        " v.rsid"
        " FROM statement s LEFT JOIN call c ON c.run_id = s.run_id AND c.vid = s.vid"
        " LEFT JOIN variant v ON v.run_id = s.run_id AND v.vid = s.vid"
        " WHERE s.run_id = ? AND s.statement_id = ?",
        [run_id, f"clinvar:{variation_id}"],
    )
    row = cursor.fetchone()
    assert row is not None, f"no statement clinvar:{variation_id}"
    return dict(zip([d[0] for d in cursor.description], row, strict=True))


def _inputs_json(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, Any]:
    row = con.execute("SELECT inputs FROM run WHERE run_id = ?", [run_id]).fetchone()
    assert row is not None
    return json.loads(row[0])


def _tools(g: Genome, run_id: str) -> dict[str, Any]:
    """`run.tools` of the run, as the store holds it."""
    with connect(g.store) as con:
        row = con.execute("SELECT tools FROM run WHERE run_id = ?", [run_id]).fetchone()
    assert row is not None
    return json.loads(row[0])


def _step(g: Genome, run_id: str, name: str) -> dict[str, Any]:
    """The step's record (key, path, reused) in the run's inputs."""
    with connect(g.store) as con:
        return _inputs_json(con, run_id)["steps"][name]


def _statements(g: Genome, run_id: str) -> dict[str, dict[str, Any]]:
    """Every statement of the run by id, as its columns."""
    with connect(g.store) as con:
        cursor = con.execute("SELECT * FROM statement WHERE run_id = ?", [run_id])
        names = [d[0] for d in cursor.description]
        return {row[1]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()}


def _accounting(g: Genome, run_id: str, tool: str) -> SimpleNamespace:
    """The accounting row of `tool` (its counts, `dropped` as a dict)."""
    columns = "input_rows, findings, watch, non_findings, unassessed, dropped"
    with connect(g.store) as con:
        row = con.execute(
            f"SELECT {columns} FROM accounting WHERE run_id = ? AND tool = ?",
            [run_id, tool],
        ).fetchone()
    assert row is not None
    names = columns.split(", ")
    return SimpleNamespace(
        **dict(zip(names, row, strict=True)) | {"dropped": json.loads(row[5])}
    )


def _new_clinvar_release(g: Genome, date: str) -> Path:
    """A second ClinVar release of the mini lab: the first with its file date changed and
    one more record (a pathogenic BRCA1 claim where the sample has no call)."""
    with gzip.open(g.clinvar38, "rt") as handle:
        lines = handle.read().splitlines()
    head = [ln for ln in lines if ln.startswith("#")]
    head = [
        f"##fileDate={date[:4]}-{date[4:6]}-{date[6:]}"
        if ln.startswith("##fileDate")
        else ln
        for ln in head
    ]
    records = [ln for ln in lines if not ln.startswith("#")]
    ref = genome_bases(g.lab.fasta, "chr1", 11000)
    info = f"CLNSIG=Pathogenic;CLNREVSTAT={TWO_STARS};GENEINFO=BRCA1:672"
    records.append(f"1\t11000\t8099\t{ref}\t{_alt(ref)}\t.\t.\t{info}")
    # The mini ClinVar holds contigs 1 and MT only, in that order.
    records.sort(key=lambda ln: (ln.split("\t")[0] != "1", int(ln.split("\t")[1])))
    path = g.clinvar38.with_name(f"clinvar38-{date}.vcf.gz")
    plain = path.with_suffix("")
    plain.write_text("\n".join([*head, *records]) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    plain.unlink()
    return path


def _counts(store: Path) -> dict[str, tuple[int, ...] | None]:
    with connect(store) as con:
        return {
            table: con.execute(f"SELECT count(*) FROM {table}").fetchone()
            for table in TABLES
        }


# --- the store's checks and the step cache ---------------------------------------------------


def _check_store(tmp_path: Path, rows: list[tuple[str, str | None, str | None]]) -> Any:
    """A store holding one run's statements (kind, vid) and the calls `CALLS` gives."""
    con = connect(tmp_path / "findings.duckdb")
    run_id = pipeline.start_run(con, {}, {})
    for vid, dosage in {"V_UNKNOWN": "unknown", "V_HET": "1", "V_ZERO": "0"}.items():
        con.execute(
            "INSERT INTO call (run_id, vid, dosage, level, row_type)"
            " VALUES (?, ?, ?, 'low', 'PASS')",
            [run_id, vid, dosage],
        )
    for n, (kind, vid, _note) in enumerate(rows):
        con.execute(
            "INSERT INTO statement (run_id, statement_id, kind, vid, gene, validity,"
            " source_level, source, source_version, topic)"
            " VALUES (?, ?, ?, ?, 'G', 'strong', 'l', 'clinvar', 'v', 'risks')",
            [run_id, f"clinvar:{n}", kind, vid],
        )
    return con, run_id


# Consistent: a finding on a known call, a gene-level finding, unassessed on an unknown
# call or on none, and a non-finding with no call.
CONSISTENT = [
    ("finding", "V_HET", None),
    ("finding", None, None),
    ("unassessed", "V_UNKNOWN", None),
    ("unassessed", "V_NONE", None),
    ("non_finding", "V_NONE", None),
    ("non_finding", "V_ZERO", None),
]


MISMATCH = "variant findings on an unknown, missing or zero-copy call"


def test_check_calls_passes_a_run_whose_kinds_agree_with_their_calls(
    tmp_path: Path,
) -> None:
    con, run_id = _check_store(tmp_path, CONSISTENT)
    pipeline.check_calls(con, run_id)


@pytest.mark.parametrize(
    ("kind", "vid", "counts"),
    [
        ("finding", "V_UNKNOWN", f"1 {MISMATCH} and 0"),
        ("finding", "V_NONE", f"1 {MISMATCH} and 0"),
        ("finding", "V_ZERO", f"1 {MISMATCH} and 0"),
        ("unassessed", "V_HET", f"0 {MISMATCH} and 1"),
        ("unassessed", "V_ZERO", f"0 {MISMATCH} and 1"),
    ],
)
def test_check_calls_raises_when_a_kind_contradicts_its_call(
    tmp_path: Path, kind: str, vid: str, counts: str
) -> None:
    con, run_id = _check_store(tmp_path, [*CONSISTENT, (kind, vid, None)])
    with pytest.raises(pipeline.CallMismatch, match=counts):
        pipeline.check_calls(con, run_id)


def test_a_cached_artefact_that_no_longer_resolves_is_rebuilt(tmp_path: Path) -> None:
    targets = [tmp_path / "old.txt"]
    targets[0].write_text("x")

    def build(out: Path) -> None:
        (out / "link.txt").symlink_to(targets[-1])

    def cached() -> Any:
        return pipeline._cached(tmp_path / "w", "step", "1", ["a"], "link.txt", build)

    first = cached()
    assert cached().reused is True
    targets[0].rename(tmp_path / "new.txt")
    targets.append(tmp_path / "new.txt")
    second = cached()
    assert (first.reused, second.reused, second.key) == (False, False, first.key)
    assert second.path.read_text() == "x"
    assert cached().reused is True


def test_a_step_another_run_finished_first_is_reused(tmp_path: Path) -> None:
    work = tmp_path / "w"
    final = work / "step" / pipeline._key("step", "1", "a")

    def build(out: Path) -> None:
        # Another process renames its complete directory into place meanwhile.
        (out / "out.txt").write_text("mine")
        final.mkdir()
        (final / "out.txt").write_text("theirs")

    step = pipeline._cached(work, "step", "1", ["a"], "out.txt", build)
    assert (step.reused, step.path) == (True, final / "out.txt")
    assert step.path.read_text() == "theirs"
    assert [p.name for p in (work / "step").iterdir()] == [final.name]


def _staging(work: Path) -> list[Path]:
    """What a killed run leaves: each staging folder, under its step or the CRAM's folder."""
    staged = [
        work / "call" / ".building-dead",
        work / ".pharmcat-input-dead",
        work / ".pgsc-input-dead",
        *(work / "align" / "k" / f".{p}-dead" for p in ("qc", "sex", "phase")),
    ]
    for folder in staged:
        folder.mkdir(parents=True)
        (folder / "part").write_text("x")
    return staged


def test_a_run_removes_a_killed_runs_staging_and_keeps_what_it_finished(
    tmp_path: Path,
) -> None:
    work = tmp_path / "w"
    staged = _staging(work)
    (work / "align" / "k" / "aligned.cram").write_text("done")
    with pipeline.run_lock(work):
        pass
    assert [f for f in staged if f.exists()] == []
    assert (work / "align" / "k" / "aligned.cram").read_text() == "done"


def test_a_live_runs_staging_survives_another_runs_start(tmp_path: Path) -> None:
    work = tmp_path / "w"
    with pipeline.run_lock(work):
        staged = _staging(work)
        with (
            pytest.raises(RunRefused, match="another run is under way"),
            pipeline.run_lock(work),
        ):
            pass
        assert all(f.exists() for f in staged)


# --- `genokeeper run` on the mini GRCh38 lab --------------------------------------------------

# (variation id, position on the mini chr1, gene, CLNSIG): the claims of the mini ClinVar.
BRCA1_HET = (8001, 1000, "BRCA1:672", "Pathogenic")
CGD_RECESSIVE = (8002, 2000, "GENEX:1", "Pathogenic")
NO_INHERITANCE = (8003, 3000, "GENEY:2", "Pathogenic")
REFCALL = (8004, 4000, "BRCA1:672", "Pathogenic")
NOCALL = (8005, 5000, "GENEY:2", "Pathogenic")
NO_RECORD = (8006, 6000, "BRCA1:672", "Likely_pathogenic")
BENIGN_38 = (8007, 7000, "BRCA1:672", "Benign")
VUS_UNLISTED = (8008, 8000, "GENEY:2", "Uncertain_significance")
WRONG_REF = (8009, 9000, "BRCA1:672", "Pathogenic")
DRUG = (8010, 10000, "CYP2C19:1557", "drug_response")
# A second het in the CGD-recessive gene: with CGD_RECESSIVE, a candidate pair to phase.
SECOND_RECESSIVE = (8011, 2100, "GENEX:1", "Likely_pathogenic")
# Inside the gVCF's reference block (BLOCK): two reference copies, from the block's GQ.
IN_BLOCK = (8013, 20050, "BRCA1:672", "Pathogenic")
BLOCK = (20000, 20100)
GENOME_CLAIMS = (
    BRCA1_HET,
    CGD_RECESSIVE,
    NO_INHERITANCE,
    REFCALL,
    NOCALL,
    NO_RECORD,
    BENIGN_38,
    VUS_UNLISTED,
    WRONG_REF,
    DRUG,
    SECOND_RECESSIVE,
    IN_BLOCK,
)
# Deletions of 60 bases written with explicit bases: one inside the reference
# block, one from a stretch the gVCF has no row for across the no-call. Structural variants,
# which DeepVariant's gVCF cannot genotype, so no claim: neither "not a carrier" from the
# block nor "unknown" from the missing rows.
LONG_DELETIONS = (
    (8015, BLOCK[0] + 10, "BRCA1:672", "Pathogenic"),
    (8016, NOCALL[1] - 10, "GENEY:2", "Pathogenic"),
)
# On chrM (ClinVar's MT), called from mutserve: a homoplasmic variant, and a position
# mutserve covers without a variant, which is the reference.
CHRM_CLAIM = (8012, 3000, "MT-TL1:4567", "Pathogenic")
CHRM_REFERENCE = (8014, 4000, "MT-TL1:4567", "Pathogenic")
# Slice 1b. PharmCAT's positions (chrom, pos, gene), each written into its input VCF from
# our call there: the carried drug_response allele, a position in the reference block, the
# no-call, and G6PD on the male's chrX, where no record is one unknown allele.
PGX_POSITIONS = (
    ("chr1", DRUG[1], "CYP2C19"),
    ("chr1", 20060, "CYP2C19"),
    ("chr1", NOCALL[1], "CYP2C19"),
    ("chrX", 100, "G6PD"),
)
PGX_GENOTYPES = ["1/1", "0/0", "./.", "."]
# The rsIDs the mini ClinVar's records carry, where PharmGKB's annotations are located.
RSIDS = {BRCA1_HET: "8001", REFCALL: "8004", NOCALL: "8005", DRUG: "8010"}
# The outside calls PharmCAT is given, from the fake callers: CYP2D6 agreed on the activity
# score alone (Aldy offers two diplotypes), CYP2B6 from Aldy, the HLA genes where OptiType
# and T1K agree (OptiType's diplotype), MT-RNR1 from mutserve (no CPIC variant).
OUTSIDE_LINES = (
    "CYP2D6\t\t\t1.0\n"
    "CYP2B6\t*1/*6\n"
    "HLA-A\t*01:01/*01:01\n"
    "HLA-B\t*08:01/*57:01\n"
    "MT-RNR1\tReference\n"
)


# Slice 1d. The scores' variants (each chosen score's file holds both: the het BRCA1 site
# and a position in the reference block) and the panel's ancestry positions (one where the
# gVCF has no record, and the block's again): the target VCF holds each once, from our call.
PRS_SCORED = (BRCA1_HET[1], 20070)
PRS_ANCESTRY = (NO_RECORD[1], 20070)
PRS_GENOTYPES = ["0/1", "./.", "0/0"]
# The GWAS Catalog's associations, located at those ancestry positions (the panel's
# alleles there: the reference and `_alt`'s A). At the block's position a risk allele that
# is the reference (two copies: a finding) and one that is the ALT (none); at the position
# with no record the ALT (unknown). Then rows that are no statement: a weaker study of the
# first, a SNP the panel lacks, and a p-value above 5e-8.
GWAS_RELEASE = "20261003"
GWAS_ROWS = [
    Assoc(20070, "C"),
    Assoc(
        20070,
        "A",
        effect="0.02",
        ci="[0.01-0.03] unit decrease",
        trait="body height",
        uris=f"{EFO}EFO_0004339",
        gene="HMGA2",
    ),
    Assoc(6000, "A", trait="LDL cholesterol", uris=f"{EFO}EFO_0004611", gene="APOB"),
    Assoc(20070, "C", p="1E-9", study="GCST000002"),
    Assoc(BRCA1_HET[1], "A"),
    Assoc(20070, "C", p="1E-5"),
]


@dataclass(frozen=True)
class Genome:
    lab: GenomeLab
    clinvar38: Path
    deepvariant: FakeDeepVariant
    verifybamid: FakeVerifyBamID
    mutserve: FakeMutserve
    whatshap: FakeWhatsHap
    pgx: FakePgx
    table: CalibrationTable
    nextflow: FakeNextflow
    haplogrep: FakeHaplogrep
    vep: FakeVep
    stage3: FakeStage3
    deletions: FakeDeletions
    deletion_table: DeletionTable
    slice3d: FakeSlice3d
    haplocheck: FakeHaplocheck
    typing: FakeTyping

    @property
    def store(self) -> Path:
        return self.lab.paths.reports / STORE_NAME


@pytest.fixture
def genome(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Genome:
    """The mini GRCh38 lab with calibration tables of its own pipeline (stage 2b's, stage
    3's and the deletion screen's), its ClinVar release, slice 1b's inputs, the `sv` step's
    files, and the containers faked."""
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")
    # Trait sites off the mini contigs or past their ends are left out (read unknown).
    monkeypatch.setattr(
        pipeline,
        "TRAIT_ALLELES",
        tuple(a for a in TRAIT_ALLELES if a[1] <= GRCH38_CONTIGS.get(a[0], 0)),
    )
    lab = write_lab(tmp_path, cgd=[("GENEX", "AR")])
    table = write_calibration(lab, tmp_path / "calibration-scratch")
    write_stage3(lab)
    write_sv(lab)
    deletion_table = write_deletions(lab)
    write_slice3d(lab)

    def base(pos: int) -> str:
        return genome_bases(lab.fasta, "chr1", pos)

    records = []
    for claim in GENOME_CLAIMS:
        variation_id, pos, gene, clnsig = claim
        # WRONG_REF names a REF the reference does not have.
        ref = _alt(base(pos)) if (variation_id, pos) == WRONG_REF[:2] else base(pos)
        info = {
            "CLNSIG": clnsig,
            "CLNREVSTAT": TWO_STARS,
            "GENEINFO": gene,
            "CLNDN": "Some_disease",
        }
        if claim in RSIDS:
            info["RS"] = RSIDS[claim]
        records.append(("1", pos, variation_id, ref, _alt(ref), info))
    for variation_id, pos, gene, clnsig in LONG_DELETIONS:
        info = {"CLNSIG": clnsig, "CLNREVSTAT": TWO_STARS, "GENEINFO": gene}
        span = genome_bases(lab.fasta, "chr1", pos, 61)
        records.append(("1", pos, variation_id, span, span[0], info))
    chrm_ref = genome_bases(lab.fasta, "chrM", CHRM_CLAIM[1])
    for variation_id, pos, gene, clnsig in (CHRM_CLAIM, CHRM_REFERENCE):
        ref = genome_bases(lab.fasta, "chrM", pos)
        info = {"CLNSIG": clnsig, "CLNREVSTAT": TWO_STARS, "GENEINFO": gene}
        records.append(("MT", pos, variation_id, ref, _alt(ref), info))
    clinvar38 = write_clinvar(tmp_path / "clinvar38.vcf.gz", "GRCh38", records)

    def at(claim: tuple[int, int, str, str], **kw: Any) -> str:
        pos = claim[1]
        return genome_variant(lab.fasta, "chr1", pos, [_alt(base(pos))], **kw)

    gvcf = [
        at(BRCA1_HET, gt="0/1", gq=50, ad=[15, 15]),
        at(CGD_RECESSIVE, gt="0/1", gq=45, ad=[14, 16]),
        at(SECOND_RECESSIVE, gt="0/1", gq=45, ad=[15, 15]),
        at(NO_INHERITANCE, gt="0/1", gq=45, ad=[14, 16]),
        at(REFCALL, gt="0/0", gq=30, ad=[25, 3], filt="RefCall"),
        at(NOCALL, gt="./.", gq=5, ad=[2, 2], filt="NoCall"),
        at(DRUG, gt="1/1", gq=40, ad=[0, 28]),
        genome_block(lab.fasta, "chr1", *BLOCK, gq=50, min_dp=20),
        # ClinVar classifies it in no record: slice 1e's VEP annotates it (`FakeVep`).
        genome_variant(
            lab.fasta, "chr1", PREDICTED_POS, [_alt(base(PREDICTED_POS))],
            gt="0/1", gq=45, ad=[15, 15],
        ),
    ]  # fmt: skip
    fake_dv = FakeDeepVariant(gvcf, dict(GRCH38_CONTIGS))
    fake_vb = FakeVerifyBamID()
    fake_mutserve = FakeMutserve([(CHRM_CLAIM[1], chrm_ref, _alt(chrm_ref), 0.99)])
    fake_whatshap = FakeWhatsHap("in_trans")
    monkeypatch.setattr(deepvariant, "run_container", fake_dv)
    monkeypatch.setattr(deepvariant, "docker_memory", lambda: 64 * 2**30)
    monkeypatch.setattr(qc, "run_container", fake_vb)
    monkeypatch.setattr(mito, "run_container", fake_mutserve)
    monkeypatch.setattr(phasing, "run_container", fake_whatshap)

    def labels(claim: tuple[int, int, str, str]) -> tuple[str, ...]:
        ref = base(claim[1])
        return (ref + ref, ref + _alt(ref), _alt(ref) * 2)

    def position(chrom: str, pos: int, gene: str) -> tuple[str, int, str, str, str]:
        ref = genome_bases(lab.fasta, chrom, pos)
        return (chrom, pos, ref, _alt(ref), gene)

    write_pgx(
        lab,
        [position(*p) for p in PGX_POSITIONS],
        [
            Ann("701", "rs8001", "TNF", "2B", labels(BRCA1_HET), drugs="etanercept"),
            Ann("702", "rs8004", "ADRB2", "2A", labels(REFCALL), drugs="salmeterol"),
            Ann("703", "rs8005", "MTHFR", "1B", labels(NOCALL), drugs="methotrexate"),
            Ann("704", "rs8010", "CYP2C19", "1B", labels(DRUG)),
            Ann("705", "rs8010", "XRCC1", "1A", labels(DRUG)),
            Ann("706", "CYP2A6*1, CYP2A6*2", "CYP2A6", "1B", ("*1/*1",)),
            Ann("707", "rs99999", "EGFR", "2B"),
        ],
    )
    fake_pgx = FakePgx()
    fake_slice3d = FakeSlice3d()
    for module in (hla, pharmcat):
        monkeypatch.setattr(module, "run_container", fake_pgx)
    monkeypatch.setattr(cyp2d6, "run_container", fake_slice3d.over(fake_pgx))
    write_prs(lab, PRS_SCORED, PRS_ANCESTRY)
    fake_nextflow = FakeNextflow()
    monkeypatch.setattr(pgsc, "run_tool", fake_nextflow)
    monkeypatch.setattr(pgsc, "docker_memory", lambda: 64 * 2**30)
    write_gwas(lab.paths.lab / "gwas" / GWAS_RELEASE / GWAS_DOWNLOAD, GWAS_ROWS)
    fake_haplogrep = FakeHaplogrep()
    monkeypatch.setattr(haplogroup, "run_container", fake_haplogrep)
    fake_vep = FakeVep()
    monkeypatch.setattr(vep, "run_container", fake_vep)
    fake_stage3 = FakeStage3()
    for module in (expansionhunter, yhaplogroup):
        monkeypatch.setattr(module, "run_container", fake_stage3)
    monkeypatch.setattr(smn, "run_container", fake_slice3d.over(fake_stage3))
    monkeypatch.setattr(expansionhunter, "run_tool", fake_stage3.tool)
    fake_deletions = FakeDeletions()
    for module in (deletion_callers, deletion_truvari):
        monkeypatch.setattr(module, "run_container", fake_slice3d.over(fake_deletions))
    for module in (cnv_callers, insertions, paralogs):
        monkeypatch.setattr(module, "run_container", fake_slice3d)
    for module in (insertions, paralogs):
        monkeypatch.setattr(module, "run_tool", fake_slice3d.tool)
    fake_haplocheck = FakeHaplocheck()
    monkeypatch.setattr(haplocheck, "run_container", fake_haplocheck)
    fake_typing = FakeTyping(lab.fasta)
    fake_typing.install(monkeypatch)
    return Genome(
        lab,
        clinvar38,
        fake_dv,
        fake_vb,
        fake_mutserve,
        fake_whatshap,
        fake_pgx,
        table,
        fake_nextflow,
        fake_haplogrep,
        fake_vep,
        fake_stage3,
        fake_deletions,
        deletion_table,
        fake_slice3d,
        fake_haplocheck,
        fake_typing,
    )


def _grading(qc_thresholds: QcThresholds = MINI_QC) -> Grading:
    return replace(load_grading(), qc=qc_thresholds)


def _run_genome(
    g: Genome,
    *,
    qc_thresholds: QcThresholds = MINI_QC,
    clinvar38: Path | None = None,
    sex: Literal["male", "female"] | None = "male",
    grading: Grading | None = None,
    **kw: Any,
) -> str:
    return pipeline.run_genome(
        g.lab.paths,
        fastq1=g.lab.fastq1,
        fastq2=g.lab.fastq2,
        clinvar38=clinvar38 or g.clinvar38,
        sex=sex,
        grading=grading if grading is not None else _grading(qc_thresholds),
        gene_validity=g.lab.gene_validity,
        cgd=g.lab.cgd,
        threads=2,
        **kw,
    )


def test_a_run_opens_no_network_connection(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """README, "Privacy and network": `genokeeper run` needs no network once the lab is built.
    Every name lookup and socket connection from this process fails the test. Bound: a
    subprocess opens its own sockets; the containers have none (`--network none`, tested in
    test_containers.py), and the native tools (bwa, samtools) open none."""
    opened: list[str] = []

    def refuse(*args: object, **kwargs: object) -> NoReturn:
        opened.append(repr(args[:2]))
        raise OSError("network refused by the test")

    monkeypatch.setattr(socket, "getaddrinfo", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)
    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(socket.socket, "connect_ex", refuse)
    assert _run_genome(genome)
    assert opened == []


def test_run_records_the_analyser(genome: Genome) -> None:
    run_id = _run_genome(genome)
    assert _tools(genome, run_id)["analyser"] == analyser_identity()


def test_after_a_run_on_the_newest_files_a_refresh_needs_no_run(
    genome: Genome,
) -> None:
    # `refresh.run_needed` reads `run.inputs` as `run_genome` records them, not as a fake
    # writes them: right after a run on the lab's newest files, no grading input differs.
    lab = genome.lab.paths.lab
    newest = lab / "clinvar" / "20261011" / "clinvar_GRCh38.vcf.gz"
    newest.parent.mkdir(parents=True)
    shutil.copy(genome.clinvar38, newest)
    shutil.copy(f"{genome.clinvar38}.tbi", f"{newest}.tbi")
    _run_genome(genome, clinvar38=newest)
    latest = snapshots(genome.lab.paths.reports / "snapshots")[-1]
    assert not refresh.run_needed(lab, latest, analyser_identity())
    dosage = lab / "clingen-dosage" / "29991231" / DOSAGE_FILE
    dosage.parent.mkdir()
    dosage.write_text(
        next((lab / "clingen-dosage").glob(f"*/{DOSAGE_FILE}")).read_text() + "\n"
    )
    assert refresh.run_needed(lab, latest, analyser_identity())


def test_run_genome_end_to_end(genome: Genome) -> None:
    run_id = _run_genome(genome)
    with connect(genome.store) as con:

        def statement(claim: tuple[int, int, str, str]) -> dict[str, Any]:
            return _row(con, run_id, claim[0])

        accounting = con.execute(
            "SELECT tool, input_rows, findings, watch, non_findings, unassessed, dropped"
            " FROM accounting WHERE run_id = ? AND tool = 'clinvar'",
            [run_id],
        ).fetchall()
        without_call = con.execute(
            "SELECT count(*) FROM statement s LEFT JOIN call c"
            " ON c.run_id = s.run_id AND c.vid = s.vid"
            " WHERE s.run_id = ? AND s.vid IS NOT NULL AND c.vid IS NULL",
            [run_id],
        ).fetchone()
        row_types = dict(
            con.execute(
                "SELECT row_type, count(*) FROM call WHERE run_id = ? GROUP BY row_type",
                [run_id],
            ).fetchall()
        )
        inputs = _inputs_json(con, run_id)
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
        kept = (
            *GENOME_CLAIMS[:6],
            DRUG,
            SECOND_RECESSIVE,
            IN_BLOCK,
            CHRM_CLAIM,
            CHRM_REFERENCE,
        )
        rows = {c[0]: statement(c) for c in kept}

    # Every claim the release keeps has a call row, every row type among them; the
    # accounting balances.
    assert without_call == (0,)
    assert row_types == {
        # Slice 1e's predicted stop-gain among them, and slice 3d's insertion: Scramble's
        # HG002 Alu at chr22:50623930 lies in ARSA's exon.
        "PASS": 8,
        "RefCall": 1,
        "NoCall": 1,
        "none": 1,
        "ref_block": 3,  # a GWAS association's in the block
    }
    ((tool, input_rows, findings, watch, non_findings, unassessed, dropped),) = (
        accounting
    )
    reasons = {k: v for k, v in json.loads(dropped).items() if v}
    assert tool == "clinvar"
    # The carried drug_response claim is a finding since slice 1b.
    assert (input_rows, findings, watch, non_findings, unassessed) == (16, 6, 0, 3, 2)
    assert reasons == {
        "benign": 1,
        "vus_outside_lists": 1,
        "ref_disagrees_with_reference": 1,
        "structural_variant": 2,
    }
    assert input_rows == findings + watch + non_findings + unassessed + sum(
        reasons.values()
    )

    het = rows[BRCA1_HET[0]]
    assert (het["kind"], het["topic"], het["dosage"], het["row_type"]) == (
        "finding",
        "risks",
        "1",
        "PASS",
    )
    # The CGD-recessive gene's two hets, phased in trans by reads: a finding each.
    for pair in (CGD_RECESSIVE, SECOND_RECESSIVE):
        assert rows[pair[0]]["kind"] == "finding"
        assert IN_TRANS_REASON in rows[pair[0]]["reason"]
    assert genome.whatshap.sites == [("chr1", 2000), ("chr1", 2100)]
    assert inputs["phased_pairs"] == {"in_trans": 1, "in_cis": 0, "unphased": 0}
    # chrM from mutserve, never from DeepVariant: one copy at its level, medium; a covered
    # position without a variant is the reference.
    chrm = rows[CHRM_CLAIM[0]]
    assert (chrm["kind"], chrm["dosage"], chrm["row_type"], chrm["level"]) == (
        "finding",
        "1",
        "PASS",
        "medium",
    )
    chrm_reference = rows[CHRM_REFERENCE[0]]
    assert (chrm_reference["kind"], chrm_reference["dosage"]) == ("non_finding", "0")
    assert chrm_reference["row_type"] == "ref_block"
    assert genome.mutserve.runs == 1
    # A claim inside a reference block: two reference copies at the block's GQ.
    in_block = rows[IN_BLOCK[0]]
    assert (in_block["kind"], in_block["gt"], in_block["dosage"]) == (
        "non_finding",
        "0/0",
        "0",
    )
    assert (in_block["row_type"], in_block["level"]) == ("ref_block", "high")
    # p_wrong is the rate measured on HG003 for the call's cell, with its interval: here
    # one error in a million units in every cell.
    errors, total = CALIBRATED_ERRORS
    lo, hi = jeffreys(errors, total)
    for claim, cell in (
        (BRCA1_HET, CellKey("variant", "snv", "easy", 50)),
        (IN_BLOCK, CellKey("ref_block", "snv", "easy", 50)),
        (REFCALL, CellKey("RefCall", "snv", "easy", 30)),
    ):
        row = rows[claim[0]]
        assert genome.table.rate_for(cell).rate == errors / total
        assert (row["p_wrong"], row["p_wrong_lo"], row["p_wrong_hi"]) == pytest.approx(
            (errors / total, lo, hi)
        )
        assert "outside calibrated regions" not in row["warnings"]
    assert rows[NO_INHERITANCE[0]]["kind"] == "finding"
    assert INHERITANCE_REASON in rows[NO_INHERITANCE[0]]["reason"]
    assert (rows[REFCALL[0]]["kind"], rows[REFCALL[0]]["dosage"]) == (
        "non_finding",
        "0",
    )
    assert rows[NOCALL[0]]["kind"] == "unassessed"
    assert rows[NO_RECORD[0]]["kind"] == "unassessed"
    assert rows[NO_RECORD[0]]["row_type"] == "none"
    # A carried drug_response claim (slice 1b): a finding under drugs.
    assert (rows[DRUG[0]]["kind"], rows[DRUG[0]]["topic"]) == ("finding", "drugs")

    # QC was recorded, with the grading inputs.
    assert inputs["qc"]["metrics"]["mean_depth"] > 0
    assert inputs["sex"] == {"value": "male", "decided_by": "given"}
    assert inputs["grading"] == {"given": True}
    assert inputs["gene_validity"] == {
        "path": str(genome.lab.gene_validity),
        "sha256": sha256_file(genome.lab.gene_validity),
        "date": "20261018",
        "present": True,
    }
    assert inputs["cgd"] == {
        "path": str(genome.lab.cgd),
        "sha256": sha256_file(genome.lab.cgd),
        "present": True,
    }
    manifest = faf_folder(genome.lab.paths.lab) / "manifest.json"
    assert (inputs["faf"]["path"], inputs["faf"]["sha256"]) == (
        str(manifest),
        sha256_file(manifest),
    )
    assert genome.verifybamid.calls == 1 and genome.deepvariant.calls == 1
    assert set(inputs["steps"]) >= {
        "align",
        "qc",
        "call",
        "claims",
        "genotypes",
        "mito",
        "phases",
        "pgx_alleles",
        "pgx_genotypes",
        "bcyrius",
        "aldy",
        "hla_reads",
        "optitype",
        "t1k",
        "pharmcat",
        "prs_genotypes",
        "pgsc",
        "haplogroup",
        "gwas_alleles",
        "gwas_genotypes",
        "repeats",
        "smn",
        "yleaf",
        "delly",
        "dysgu",
        "deletions_union",
        "coverage",
        "dosage",
        *SLICE3D_STEPS,
    }
    assert not any(step["reused"] for step in inputs["steps"].values())
    assert tools_row is not None
    tools_json = json.loads(tools_row[0])
    assert tools_json["clinvar"] == CLINVAR_RELEASE
    assert tools_json["gnomad_faf"] == (
        f"{inputs['faf']['release']}, floor {inputs['faf']['floor']}"
    )
    # The run records the six versions the table is checked against, and the table.
    assert {key: tools_json[key] for key in VERSION_KEYS} == genome.table.versions
    folder = table_dir(genome.lab.paths.lab)
    bench = layout(genome.lab.paths.lab)["truth"] / BENCH_BED
    assert inputs["calibration"] == {
        "path": str(folder),
        "manifest_sha256": sha256_file(folder / "manifest.json"),
        "benchmark": {"path": str(bench), "sha256": sha256_file(bench)},
        "hg002_cells": {"reason": NOT_MEASURED},
    }
    # What a refresh compares with the lab: the run recorded each lab table as found.
    lab = genome.lab.paths.lab
    assert "small_variants" in inputs["tables"]
    for name in refresh.LAB_TABLES:
        assert refresh._recorded(inputs, name) == refresh._on_disk(lab, name)  # pyright: ignore[reportPrivateUsage]

    # A rerun reuses every step: nothing is aligned, measured or called again.
    again = _run_genome(genome)
    assert again != run_id
    with connect(genome.store) as con:
        rerun = _inputs_json(con, again)
    assert all(step["reused"] for step in rerun["steps"].values())
    assert {"coverage", "dosage", "repeats"} <= set(rerun["steps"])
    assert genome.verifybamid.calls == 1 and genome.deepvariant.calls == 1
    assert genome.mutserve.runs == 1 and len(genome.whatshap.sites) == 2
    assert genome.pgx.runs == dict.fromkeys(
        ("bcyrius", "aldy", "optitype", "t1k", "pharmcat"), 1
    ) | {"aldy": 2}  # CYP2D6 and CYP2B6
    assert genome.nextflow.runs == 1
    assert genome.haplogrep.runs == 1
    # No HG002 allele is above its normal range, so REViewer draws nothing.
    # Parascopy's depth is the shared step's (`FakeSlice3d`), run once for SMN1 and the
    # paralog loci.
    assert genome.stage3.runs == dict.fromkeys(
        ("expansionhunter", "smncaller", "parascopy", "yleaf"), 1
    )
    assert genome.slice3d.runs == dict.fromkeys(
        (
            *("parascopy depth", "cnvpytor", "insurveyor", "scramble", "truvari"),
            *("parascopy cn", "gauchian", "kilda", "aldy"),
        ),
        1,
    ) | {"delly cnv": 3}  # one per primary chromosome, as 3b's DELLY
    # Each caller one container per primary chromosome of the mini FASTA (chr1, chrX, chrY).
    assert genome.deletions.runs == {"delly": 3, "dysgu": 3, "truvari": 1}
    # Dysgu's shards reach its cached step: every shard's records, kept on its chromosome.
    dysgu = Path(inputs["steps"]["dysgu"]["path"])
    with gzip.open(dysgu, "rt") as handle:
        chroms = [x.split("\t", 1)[0] for x in handle if not x.startswith("#")]
    assert chroms and set(chroms) == {"chr1", "chrX", "chrY"}


def test_run_writes_predicted_findings(genome: Genome) -> None:
    """Slice 1e in every run: VEP is given each PASS allele the user carries, offline in
    its image behind dna:annotate (SpliceAI on the listed genes), and grading keeps those
    no ClinVar record classifies germline; the stop-gain it predicts on BRCA1's MANE transcript (haploinsufficient,
    expressed, no NMD escape) is a finding under risks, validity limited, "predicted, not
    classified", with its call; the accounting counts every variant VEP was given; the
    feed puts it after the ClinVar findings."""
    run_id = _run_genome(genome)
    predicted = f"GRCh38:chr1:{PREDICTED_POS}:"
    classified = f"GRCh38:chr1:{CLASSIFIED_LOF_POS}:"
    # VEP is given the classified allele too: ClinVar's filter is at grading.
    assert sum(vid.startswith(predicted) for vid in genome.vep.sites) == 1
    assert sum(vid.startswith(classified) for vid in genome.vep.sites) == 1
    vep_run, *spliceai = genome.vep.runs
    assert vep_run == ("vep", "dna:annotate")
    assert spliceai and set(spliceai) == {("spliceai", "dna:spliceai")}
    with connect(genome.store) as con:
        (row,) = con.execute(
            "SELECT s.kind, s.topic, s.validity, s.source, s.reason, s.actionable,"
            " c.dosage FROM statement s JOIN call c ON c.run_id = s.run_id"
            " AND c.vid = s.vid WHERE s.run_id = ? AND s.statement_id LIKE 'predicted:%'",
            [run_id],
        ).fetchall()
        accounting = con.execute(
            "SELECT input_rows, findings, dropped FROM accounting"
            " WHERE run_id = ? AND tool = 'vep'",
            [run_id],
        ).fetchone()
        inputs = _inputs_json(con, run_id)
        loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.LIMITED)
        order = [f["sources"][0]["source"] for f in query(con, None, loosest).findings]
    kind, topic, validity, source, reason, actionable, dosage = row
    assert (kind, topic, validity, source, actionable, dosage) == (
        "finding",
        "risks",
        "limited",
        "vep",
        False,
        "1",
    )
    assert reason.startswith("predicted, not classified") and "geneticist" in reason
    assert accounting is not None and accounting[:2] == (1, 1)
    assert not any(json.loads(accounting[2]).values())
    assert order.index("vep") > max(i for i, s in enumerate(order) if s == "clinvar")
    assert {"vep", "clingen_dosage"} <= set(inputs)
    assert {"clinvar_vids", "vep", "predicted_genotypes"} <= set(inputs["steps"])


def test_new_clinvar_release_reuses_vep(genome: Genome) -> None:
    first = _run_genome(genome)
    second = _run_genome(genome, clinvar38=_new_clinvar_release(genome, "20261101"))
    assert _tools(genome, second)["clinvar"] == "20261101"
    assert _step(genome, second, "vep")["reused"] is True
    assert _step(genome, first, "vep")["key"] == _step(genome, second, "vep")["key"]
    # The release's filter still runs at grading: its new records are another ClinVar's.
    assert (
        _step(genome, first, "clinvar_vids")["key"]
        != (_step(genome, second, "clinvar_vids")["key"])
    )


def test_clinvar_classified_allele_is_never_predicted(genome: Genome) -> None:
    run_id = _run_genome(genome)
    ref = genome_bases(genome.lab.fasta, "chr1", CLASSIFIED_LOF_POS)
    classified_lof_vid = f"GRCh38:chr1:{CLASSIFIED_LOF_POS}:{ref}:{_alt(ref)}"
    assert classified_lof_vid in genome.vep.sites
    ids = _statements(genome, run_id)
    assert not any(i.startswith(f"predicted:{classified_lof_vid}") for i in ids)
    # VEP was given the classified alleles too; grading neither predicts nor counts them.
    with gzip.open(_step(genome, run_id, "clinvar_vids")["path"], "rt") as handle:
        classified = set(handle.read().split())
    assert classified_lof_vid in classified
    acc = _accounting(genome, run_id, "vep")
    assert acc.input_rows == len(set(genome.vep.sites) - classified) == 1


def test_a_chrm_allele_vep_predicts_is_not_counted_no_predicted_effect(
    genome: Genome,
) -> None:
    # DeepVariant's chrM call, a stop-gain to the fake VEP: mutserve's, never VEP's.
    pos = CLASSIFIED_LOF_POS
    ref = genome_bases(genome.lab.fasta, "chrM", pos)
    genome.deepvariant.contigs["chrM"] = CHRM_LENGTH
    genome.deepvariant.records.append(
        genome_variant(
            genome.lab.fasta, "chrM", pos, [_alt(ref)], gt="1/1", gq=40, ad=[0, 30]
        )
    )
    run_id = _run_genome(genome)
    assert f"GRCh38:chrM:{pos}:{ref}:{_alt(ref)}" in genome.vep.sites
    acc = _accounting(genome, run_id, "vep")
    assert (acc.input_rows, acc.findings, acc.dropped["no_predicted_effect"]) == (
        1,
        1,
        0,
    )


def test_missing_vep_files_refused_before_any_step(genome: Genome) -> None:
    vep.vep_files(genome.lab.paths.lab).pext_max.unlink()
    _refused_before_any_step(genome, "no VEP 116 file.*bootstrap.sh --with-vep")


def test_run_states_the_haplogroup(genome: Genome) -> None:
    """Every run states the mitochondrial haplogroup from mutserve's homoplasmic calls,
    with Haplogrep3's version, tree and image in `run.tools`; when too much of chrM is
    uncovered, Haplogrep3 is not run and the haplogroup is unassessed with its reason."""
    run_id = _run_genome(genome)
    chrm_ref = genome_bases(genome.lab.fasta, "chrM", CHRM_CLAIM[1])
    # The homoplasmic call (0.99) is the VCF's one row.
    assert genome.haplogrep.vcf.splitlines()[-1] == (
        f"chrM\t{CHRM_CLAIM[1]}\t.\t{chrm_ref}\t{_alt(chrm_ref)}\t.\tPASS\t.\tGT\t1"
    )
    with connect(genome.store) as con:
        row = con.execute(
            "SELECT kind, topic, validity, claim, reason FROM statement"
            " WHERE run_id = ? AND statement_id = 'haplogrep3:haplogroup'",
            [run_id],
        ).fetchone()
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
    assert row is not None and row[:3] == ("finding", "ancestry", "strong")
    assert json.loads(row[3])["haplogroup"] == "J1c" and row[4] is None
    assert tools_row is not None
    tools_json = json.loads(tools_row[0])
    assert (tools_json["haplogrep3"], tools_json["haplogrep3_tree"]) == (
        haplogroup.HAPLOGREP_VERSION,
        haplogroup.TREE,
    )
    assert tools_json["haplogrep3_image"] == IMAGES["haplogrep3"]

    genome.mutserve.uncovered = range(1, 201)
    shutil.rmtree(genome.lab.paths.reports / "work" / "mito")
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        row = con.execute(
            "SELECT kind, reason FROM statement"
            " WHERE run_id = ? AND statement_id = 'haplogrep3:haplogroup'",
            [run_id],
        ).fetchone()
    assert row == (
        "unassessed",
        "200 of 16569 chrM positions below 10 reads; no haplogroup stated",
    )
    assert genome.haplogrep.runs == 1


STAGE3_TOOLS = ("expansionhunter", "smn", "yleaf")


def _stage3(
    con: duckdb.DuckDBPyConnection, run_id: str
) -> tuple[dict[str, dict[str, Any]], dict[str, tuple[int, ...]]]:
    """The run's stage-3 statements by statement id (with their gene calls), and each
    stage-3 tool's accounting (input_rows, findings, watch, non_findings, unassessed)."""
    cursor = con.execute(
        "SELECT s.statement_id, s.source, s.kind, s.reason, s.topic, s.actionable,"
        " s.claim, s.validity, s.exclusion, g.phenotype, g.level FROM statement s"
        " LEFT JOIN gene_call g"
        " ON g.run_id = s.run_id AND g.gene = s.gene AND s.source <> 'yleaf'"
        " WHERE s.run_id = ? AND s.source IN (?, ?, ?)",
        [run_id, *STAGE3_TOOLS],
    )
    names = [d[0] for d in cursor.description]
    rows = {row[0]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()}
    accounting = {
        tool: tuple(counts)
        for tool, *counts in con.execute(
            "SELECT tool, input_rows, findings, watch, non_findings, unassessed"
            " FROM accounting WHERE run_id = ? AND tool IN (?, ?, ?)",
            [run_id, *STAGE3_TOOLS],
        ).fetchall()
    }
    return rows, accounting


# The X-linked loci, unassessed by the X-dosage rule while the run's [dosage] is not
# calibrated (the shipped grading's `calibrated_on` is empty), once the run's own checks
# pass.
X_LINKED = {"expansionhunter:FMR1", "expansionhunter:AFF2"}
# One X copy (a depth ratio of 0.5) or two, the run's X and Y depth ratios.
ONE_X = (0.5, 0.5)
TWO_X = (1.0, 0.0)
# HG002's slice-3e output (`hg002.eh3e.*`, male, one X) as the run grades it: Task 4's
# record, with the flanking condition off for two distinct autosomal alleles (CNBP 15/16 and
# FGF14 17/19 normal). AFF2 needs a calibrated one-X dosage.
ELEVEN_NORMAL = {"HTT", "DMPK", "C9orf72", "ATXN1", "ATXN3", "ATXN7", "ATN1", "FXN"}
ELEVEN = ELEVEN_NORMAL | {"ATXN2", "CACNA1A", "AR"}
NEW_NORMAL = {"TCF4", "PABPN1", "CNBP", "ATXN10", "TBP", "FGF14", "AFF2"}
SAME_LENGTH = (
    "not assessed (the two alleles are the same length: a longer one may be hidden);"
    " needs PCR sizing"
)
NEW_UNASSESSED = {
    # 31 CGG on one spanning read.
    "FMR1": "not assessed (an allele has fewer than 3 spanning reads); needs a targeted"
    " PCR test",
    "NOP56": "not assessed (an allele has fewer than 3 spanning reads); needs"
    " repeat-primed PCR",
    **dict.fromkeys(("PPP2R2B", "JPH3", "ZFHX3", "THAP11"), SAME_LENGTH),
}
EXCLUSION_NORMAL = {"CNBP", "ATXN10", "TBP", "FGF14", "AFF2", "RFC1"}


def _dosage_at(
    monkeypatch: pytest.MonkeyPatch, ratios: tuple[float, float], **kw: Any
) -> Grading:
    """A grading whose [dosage] is calibrated (`kw` replacing its fields), the run's X and
    Y depth ratios `ratios`."""
    monkeypatch.setattr(pipeline, "measure_ratios", lambda *_a, **_k: ratios)
    g = _grading()
    fields: dict[str, Any] = {"calibrated_on": "2026-10-30 test"} | kw
    return replace(g, dosage=replace(g.dosage, **fields))


def _eh_line(vcf: str, varid: str, alt: str, sample: str) -> str:
    """`vcf` with `varid`'s ALT and sample column replaced."""
    lines = vcf.splitlines()
    (i,) = [n for n, line in enumerate(lines) if f";VARID={varid};" in line]
    fields = lines[i].split("\t")
    fields[4], fields[9] = alt, sample
    lines[i] = "\t".join(fields)
    return "\n".join(lines) + "\n"


# The 22 loci written as not covered, each with its own reason whatever the run.
NOT_COVERED_IDS = {f"expansionhunter:{entry.gene}" for entry in NOT_COVERED}


def _assert_stage3_unassessed(genome: Genome, run_id: str, reason: str) -> None:
    """Every stage-3 statement unassessed with `reason`, FMR1 and AFF2 included, but the
    not-covered loci's, with their own; no stage-3 container ran; each tool's accounting
    balances."""
    with connect(genome.store) as con:
        rows, accounting = _stage3(con, run_id)
        inputs = _inputs_json(con, run_id)
    repeats = [r for r in rows.values() if r["source"] == "expansionhunter"]
    assert len(repeats) == 48
    assert {r["kind"] for r in rows.values()} == {"unassessed"}
    assert NOT_COVERED_IDS <= set(rows)
    assert all(rows[i]["reason"].startswith("not covered:") for i in NOT_COVERED_IDS)
    graded = set(rows) - NOT_COVERED_IDS
    assert {"expansionhunter:RFC1", "expansionhunter:CFTR", *X_LINKED} <= graded
    assert {rows[i]["reason"] for i in graded} == {reason}
    assert {i for i, r in rows.items() if r["source"] != "expansionhunter"} == {
        "smn:SMN1",
        "yleaf:haplogroup",
    }
    assert genome.stage3.runs == {}
    assert accounting == {
        "expansionhunter": (48, 0, 0, 0, 48),
        "smn": (1, 0, 0, 0, 1),
        "yleaf": (1, 0, 0, 0, 1),
    }
    assert not {"repeats", "smn", "yleaf"} & set(inputs["steps"])
    assert inputs["stage3"] == {"table": None, "reason": reason}


def test_uncalibrated_stage3_unassessed(genome: Genome) -> None:
    """With no stage-3 table no stage-3 tool runs: every locus, SMN1 and the Y haplogroup
    are unassessed with the calibrate command, and the rest of the run is written."""
    shutil.rmtree(stage3_dir(genome.lab.paths.lab))
    run_id = _run_genome(genome)
    _assert_stage3_unassessed(genome, run_id, NOT_CALIBRATED)
    assert NOT_CALIBRATED.endswith("run `genokeeper calibrate`")


def test_missing_parascopy_data_leaves_stage3_unassessed(genome: Genome) -> None:
    """A lab without the bootstrap's Parascopy data runs: stages 1-2 are written and every
    stage-3 gene is unassessed with a reason naming the bootstrap."""
    root = genome.lab.paths.lab / "parascopy"
    shutil.rmtree(root)
    run_id = _run_genome(genome)
    reason = (
        f"stage 3 not run: {root}: no Parascopy data md5s; run scripts/bootstrap.sh"
    )
    _assert_stage3_unassessed(genome, run_id, reason)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("expansionhunter", "genomelab/expansionhunter:5.0.0-other"),
        ("exclusion", '{"min_spanning": 4}'),
    ],
)
def test_stale_stage3_table_unassessed(genome: Genome, key: str, value: str) -> None:
    """A stage-3 table measured with another ExpansionHunter image, or under another
    [exclusion], is refused as no table is, its reason naming the key that differs."""
    other: dict[str, Any] = {key: value}  # one version replaced
    write_stage3(genome.lab, **other)
    run_id = _run_genome(genome)
    reason = f"stale: {key} changed; run `genokeeper calibrate`"
    _assert_stage3_unassessed(genome, run_id, reason)


def test_calibrated_stage3_statements(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With a current table and a calibrated one-X dosage the three tools run on the run's
    CRAM and grade HG002's slice-3e output: eight of the eleven normal, ATXN2 and CACNA1A
    watch items (their intervals reach past the normal range), AR withheld by the table;
    of the thirteen new loci seven normal (AFF2 "normal, … on the one X"), FMR1 and NOP56
    on too few spanning reads and four with two alleles of one length unassessed; RFC1 and
    CFTR normal; the 22 not-covered loci unassessed; SMN1 two copies, reduced risk; the Y
    haplogroup stated. No repeat statement is actionable, none a watch but the eleven's,
    TCF4's or PABPN1's; the accounting balances; the run records the table's versions,
    ExpansionHunter's build and the five images."""
    run_id = _run_genome(genome, grading=_dosage_at(monkeypatch, ONE_X))
    with connect(genome.store) as con:
        rows, accounting = _stage3(con, run_id)
        inputs = _inputs_json(con, run_id)
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
    repeats = {
        i.removeprefix("expansionhunter:"): r
        for i, r in rows.items()
        if r["source"] == "expansionhunter"
    }
    kinds = {g: r["kind"] for g, r in repeats.items()}
    normal = ELEVEN_NORMAL | NEW_NORMAL | {"RFC1", "CFTR"}
    assert {g for g, k in kinds.items() if k == "non_finding"} == normal
    assert {g for g, k in kinds.items() if k == "watch"} == {"ATXN2", "CACNA1A"}
    assert {g for g, k in kinds.items() if k == "watch"} <= ELEVEN | {"TCF4", "PABPN1"}
    assert repeats["AR"]["kind"] == "unassessed" and repeats["AR"]["reason"] == (
        FALSE_CALL
    )
    assert {g: repeats[g]["reason"] for g in NEW_UNASSESSED} == NEW_UNASSESSED
    assert {g for g, r in repeats.items() if r["exclusion"]} == EXCLUSION_NORMAL
    aff2 = repeats["AFF2"]
    assert aff2["reason"].startswith("normal, 19 GCC on the one X")
    assert repeats["CNBP"]["reason"] == (
        "normal: 15 and 16 CAGG, each allele read across in full"
    )
    assert json.loads(repeats["RFC1"]["claim"])["rfc1"] == {
        "call": "normal",
        "description": "no expanded allele",
    }
    assert repeats["CFTR"]["reason"] == "no 5T allele: T 7/7, TG 10/12"
    # A finding at a non-standard locus, RFC1 or CFTR would be a possible one at validity
    # limited; HG002 has none.
    loci = load_grading().repeats.loci
    standard = {g for g, locus in loci.items() if locus.rule is Rule.STANDARD}
    for gene, row in repeats.items():
        if row["kind"] == "finding" and gene not in standard:
            assert row["validity"] == "limited" and "possible" in row["reason"]
    assert {f"expansionhunter:{g}" for g in repeats} >= NOT_COVERED_IDS
    assert {rows[i]["kind"] for i in NOT_COVERED_IDS} == {"unassessed"}
    assert repeats["HTT"]["phenotype"] == "normal"
    smn1 = rows["smn:SMN1"]
    assert (smn1["kind"], smn1["reason"], smn1["phenotype"]) == (
        "non_finding",
        TWO_COPIES_REASON,
        "reduced risk",
    )
    y = rows["yleaf:haplogroup"]
    assert (y["kind"], y["topic"]) == ("finding", "ancestry")
    assert json.loads(y["claim"])["haplogroup"] == "J-C134621"
    assert not any(r["actionable"] for r in rows.values())
    # 48 statements: 24 loci, RFC1, CFTR and 22 not covered. Watch: ATXN2, CACNA1A.
    # Normal: 8 of the eleven, the 7 new, RFC1 and CFTR (17). Unassessed: AR, FMR1, NOP56,
    # PPP2R2B, JPH3, ZFHX3, THAP11 and the 22 (29).
    assert accounting == {
        "expansionhunter": (48, 0, 2, 17, 29),
        "smn": (1, 0, 0, 1, 0),
        "yleaf": (1, 1, 0, 0, 0),
    }
    assert "reviewer" not in genome.stage3.runs
    versions = inputs["stage3"]["table"]
    assert inputs["stage3"]["reason"] is None
    assert versions["expansionhunter"] == IMAGES["expansionhunter"]
    assert tools_row is not None
    tools_json = json.loads(tools_row[0])
    assert tools_json["expansionhunter"] == EH_VERSION
    for name in ("expansionhunter", "reviewer", "smncaller", "parascopy", "yleaf"):
        assert tools_json[f"{name}_image"] == IMAGES[name]


def test_rfc1_failure_unassessed_rest_written(genome: Genome, tmp_path: Path) -> None:
    """ExpansionHunter's JSON without RFC1's motif analysis leaves RFC1 unassessed with the
    failure, keeps the repeats step, and the 24 loci and CFTR are graded as ever."""
    data = json.loads(genome.stage3.json.read_text())
    del data["LocusResults"]["RFC1"]["Variants"]["RFC1"]["RFC1MotifAnalysis"]
    genome.stage3.json = tmp_path / "eh.json"
    genome.stage3.json.write_text(json.dumps(data))
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        rows, accounting = _stage3(con, run_id)
        inputs = _inputs_json(con, run_id)
    rfc1 = rows["expansionhunter:RFC1"]
    assert rfc1["kind"] == "unassessed"
    assert rfc1["reason"].startswith("rfc1 failed: ")
    assert "no RFC1 motif analysis" in rfc1["reason"]
    assert "repeats" in inputs["steps"]
    assert rows["expansionhunter:HTT"]["kind"] == "non_finding"
    assert rows["expansionhunter:CNBP"]["kind"] == "non_finding"
    assert rows["expansionhunter:CFTR"]["kind"] == "non_finding"
    # The calibrated run's counts, under the shipped [dosage]: AFF2 unassessed by the X
    # dosage, RFC1 by the failure.
    assert accounting["expansionhunter"] == (48, 0, 2, 15, 31)


def test_female_run_fmr1_two_x(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A female run with two X copies grades FMR1 by the two-X rule: two distinct spanned
    alleles are normal; the table withholding `FMR1:two_x` leaves it unassessed with that
    reason."""
    grading = _dosage_at(monkeypatch, TWO_X)
    two = _eh_line(
        genome.stage3.vcf.read_text(),
        "FMR1",
        "<STR29>,<STR31>",
        "1/2:SPANNING/SPANNING:29/31:29-29/31-31:8/8:10/10:0/0:11.189189",
    )
    genome.stage3.vcf = tmp_path / "eh.vcf"
    genome.stage3.vcf.write_text(two)
    run_id = _run_genome(genome, sex="female", grading=grading)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    fmr1 = rows["expansionhunter:FMR1"]
    assert (fmr1["kind"], fmr1["exclusion"]) == ("non_finding", True)
    assert fmr1["reason"] == NORMAL_TWO_X.format(a=29, b=31, motif="CGG", gene="FMR1")
    write_stage3(
        genome.lab, withheld={**STAGE3_WITHHELD, two_x_key("FMR1"): X_RULE_DID_NOT_HOLD}
    )
    run_id = _run_genome(genome, sex="female", grading=grading)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    fmr1 = rows["expansionhunter:FMR1"]
    assert (fmr1["kind"], fmr1["reason"]) == ("unassessed", X_RULE_DID_NOT_HOLD)


def test_xxy_run_fmr1_unassessed(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Two X copies with a Y present (XXY): two distinct spanned FMR1 alleles are seen, but
    the closing design leaves an XXY sample unassessed, naming the measured dosage."""
    grading = _dosage_at(monkeypatch, (1.0, 0.5))
    two = _eh_line(
        genome.stage3.vcf.read_text(),
        "FMR1",
        "<STR29>,<STR31>",
        "1/2:SPANNING/SPANNING:29/31:29-29/31-31:8/8:10/10:0/0:11.189189",
    )
    genome.stage3.vcf = tmp_path / "eh.vcf"
    genome.stage3.vcf.write_text(two)
    run_id = _run_genome(genome, sex="female", grading=grading)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    fmr1 = rows["expansionhunter:FMR1"]
    assert (fmr1["kind"], fmr1["reason"]) == (
        "unassessed",
        X_DOSAGE.format(copies="2 X and a Y"),
    )


def test_uncalibrated_dosage_fmr1_unassessed(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A measured one-X dosage under a [dosage] whose `calibrated_on` is empty is not used:
    FMR1 and AFF2 are unassessed by the X-dosage rule, naming `genokeeper calibrate`."""
    grading = _dosage_at(monkeypatch, ONE_X, calibrated_on="")
    run_id = _run_genome(genome, grading=grading)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    assert {(rows[i]["kind"], rows[i]["reason"]) for i in X_LINKED} == {
        ("unassessed", X_NOT_CALIBRATED)
    }


def test_run_passes_min_spanning(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`[exclusion] min_spanning` reaches the source: at 50, under a table written with it,
    every non-standard locus is unassessed on too few spanning reads, and CFTR unread."""
    grading = _dosage_at(monkeypatch, ONE_X)
    grading = replace(grading, exclusion=replace(grading.exclusion, min_spanning=50))
    write_stage3(genome.lab, grading=grading)
    run_id = _run_genome(genome, grading=grading)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    loci = grading.repeats.loci
    other = {g for g, locus in loci.items() if locus.rule is not Rule.STANDARD}
    assert len(other) == 11
    for gene in other:
        row = rows[f"expansionhunter:{gene}"]
        assert row["kind"] == "unassessed", gene
        assert "an allele has fewer than 50 spanning reads" in row["reason"], gene
    assert rows["expansionhunter:CFTR"]["reason"] == CFTR_UNREAD


def test_repeats_step_keyed_on_what_its_build_reads(genome: Genome) -> None:
    """A [repeat_loci] edit that changes no call (a note, the evidence, a disease
    threshold) reuses the repeats step; one its build reads (`normal_max` picks the loci
    REViewer draws) runs it again."""
    _run_genome(genome)
    loci = _grading().repeats.loci
    htt = loci["HTT"]
    worded = replace(htt, note="reworded", evidence="other", pathogenic_min=41)
    graded = replace(
        _grading(), repeats=replace(_grading().repeats, loci={**loci, "HTT": worded})
    )
    write_stage3(genome.lab, grading=graded)  # the stage-3 table measures the thresholds
    run_id = _run_genome(genome, grading=graded)
    with connect(genome.store) as con:
        reused = _inputs_json(con, run_id)
    assert reused["steps"]["repeats"]["reused"]
    drawn = replace(htt, normal_max=htt.normal_max + 1)
    redrawn = replace(
        _grading(), repeats=replace(_grading().repeats, loci={**loci, "HTT": drawn})
    )
    write_stage3(genome.lab, grading=redrawn)  # the stage-3 table measures the thresholds
    run_id = _run_genome(genome, grading=redrawn)
    with connect(genome.store) as con:
        again = _inputs_json(con, run_id)
    assert not again["steps"]["repeats"]["reused"]


def test_expanded_allele_drawn_by_reviewer(genome: Genome, tmp_path: Path) -> None:
    """An HTT allele of 45 CAG (HG002's 24 grown past the full-penetrance threshold) is a
    finding, screening only, with the image REViewer drew of that locus alone."""
    vcf = FakeStage3().vcf.read_text()
    grown = vcf.replace("<STR24>", "<STR45>").replace(
        "17/24:17-17/24-24", "17/45:17-17/45-45"
    )
    assert grown.count("45-45") == 1
    genome.stage3.vcf = tmp_path / "eh.vcf"
    genome.stage3.vcf.write_text(grown)
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        rows, _accounting = _stage3(con, run_id)
    htt = rows["expansionhunter:HTT"]
    assert (htt["kind"], htt["actionable"]) == ("finding", False)
    assert "screening result, needs confirmation" in htt["reason"]
    image = json.loads(htt["claim"])["image"]
    assert image is not None and Path(image).name == "rev.HTT.svg"
    assert genome.stage3.runs["reviewer"] == 1


def test_low_reliability_gene_level_non_finding_is_probably_not(
    genome: Genome, tmp_path: Path
) -> None:
    """AR called with two alleles on this male's X is a non-finding at low reliability,
    which `genotype --gene AR` answers as probably not a carrier, as a site at low is."""
    vcf = FakeStage3().vcf.read_text()
    two = vcf.replace("\tT\t<STR27>\t", "\tT\t<STR20>,<STR22>\t").replace(
        "1:SPANNING:27:27-27:3:6:0:",
        "1/2:SPANNING/SPANNING:20/22:20-20/22-22:3/3:6/6:0/0:",
    )
    assert two.count("SPANNING/SPANNING:20/22") == 1
    write_stage3(genome.lab, withheld={})
    genome.stage3.vcf = tmp_path / "eh.vcf"
    genome.stage3.vcf.write_text(two)
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        found = genotype(con, run_id, gene="AR")
    (ar,) = [s for s in found.sites if s["marker"]["grch38"] is None]
    assert ar["kind"] == "non_finding"
    assert ar["answer"] == "probably_not_a_carrier"


def test_stage3_statements_answer_genotype_and_show(genome: Genome) -> None:
    """A stage-3 gene's statement answers `genotype --gene` and `show`, whatever its kind:
    HTT normal, SMN1 two copies (and no "main mechanism not assessed" line), ATXN2 a
    borderline watch."""
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        found = {g: genotype(con, run_id, gene=g) for g in ("HTT", "SMN1", "ATXN2")}
        _run, shown = show(con, run_id, "smn:SMN1")
    level = {
        g: [s for s in r.sites if s["marker"]["grch38"] is None]
        for g, r in found.items()
    }
    assert found["HTT"].assessed
    (htt,) = level["HTT"]
    assert (htt["kind"], htt["answer"]) == ("non_finding", "not_a_carrier")
    assert htt["not_a_finding_because"] == "normal"
    (smn1,) = level["SMN1"]
    assert (smn1["kind"], smn1["claim"]["smn1"]) == ("non_finding", 2)
    # Two copies leave a 2+0 carrier possible: never a plain "not a carrier".
    assert smn1["answer"] == "probably_not_a_carrier"
    assert smn1["not_a_finding_because"] == TWO_COPIES_REASON
    smn1_gene = found["SMN1"].gene
    assert (
        smn1_gene is not None and smn1_gene["dominant_mechanism_not_assessed"] is None
    )
    assert shown is not None
    assert (shown["kind"], shown["claim"]["smn1"]) == ("non_finding", 2)
    (atxn2,) = level["ATXN2"]
    assert (atxn2["kind"], atxn2["answer"]) == ("watch", "unknown")
    assert "borderline, indeterminate" in atxn2["not_a_finding_because"]


def test_fmr1_unassessed_at_the_default_filters(genome: Genome) -> None:
    """FMR1's unassessed statement is in `query`'s default output, with its reason."""
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        result = query(con, run_id, Filters())
    fmr1 = result.unassessed["by_gene"]["FMR1"]
    assert fmr1["reasons"] == {f"{GENE_UNRESOLVED}: {X_NOT_CALIBRATED}": 1}


def test_withheld_y_haplogroup_at_the_default_filters(genome: Genome) -> None:
    """A Y haplogroup the table withholds is unassessed in `query`'s default output, with
    the table's reason, as FMR1's unassessed statement is."""
    write_stage3(genome.lab, withheld={YLEAF_KEY: YLEAF_WITHHELD})
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        result = query(con, run_id, Filters())
    # Beside it the dosage statement, not calibrated under the shipped [dosage], and slice
    # 3f's HHV-6 and ROH statements without a slice 3f table (the same reason).
    assert result.unassessed["by_gene"][""]["reasons"] == {
        YLEAF_WITHHELD: 1,
        NOT_CALIBRATED: 3,
    }


def test_stage3_tool_failure_unassessed_rest_written(genome: Genome) -> None:
    """SMNCopyNumberCaller failing leaves SMN1 unassessed with the failure, keeps nothing of
    its step, and the repeat loci and the Y haplogroup are graded as ever."""
    genome.stage3.fail.add("smncaller")
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        rows, accounting = _stage3(con, run_id)
        inputs = _inputs_json(con, run_id)
    smn1 = rows["smn:SMN1"]
    assert smn1["kind"] == "unassessed"
    assert smn1["reason"] == "smn failed: exit 1: smncaller crashed"
    assert rows["expansionhunter:HTT"]["kind"] == "non_finding"
    assert rows["yleaf:haplogroup"]["kind"] == "finding"
    assert accounting["smn"] == (1, 0, 0, 0, 1)
    assert "smn" not in inputs["steps"] and {"repeats", "yleaf"} <= set(inputs["steps"])
    genome.stage3.fail.clear()
    _run_genome(genome)
    assert genome.stage3.runs["smncaller"] == 2


def test_shared_depth_fault_unassessed_everywhere(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One fault in the Parascopy depth step 3a and 3d share leaves both parts unassessed
    naming it, whichever reached it first, and never aborts the run."""

    def missing(*_a: Any, **_k: Any) -> None:
        raise tools.ToolMissing("samtools")

    monkeypatch.setattr(pipeline, "parascopy_depth_step", missing)
    run_id = _run_genome(genome)
    statements, _, _ = _slice3d(genome, run_id)
    with connect(genome.store) as con:
        rows, _ = _stage3(con, run_id)
    smn, strc = rows["smn:SMN1"], statements["paralogs:STRC"]
    assert smn["kind"] == strc["kind"] == "unassessed"
    assert smn["reason"].startswith("smn failed: samtools is not on PATH")
    assert strc["reason"].startswith("parascopy_depth failed: samtools is not on PATH")


# --- stage 3b: deletions --------------------------------------------------------------------

DELETION_TOOLS = {"delly", "dysgu", "truvari"}
SLICE3D_TOOLS = {"cnv", "insertions", "paralogs", "lpa"}
LOOSEST_FILTERS = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
# Every slice-3d step a calibrated run makes.
SLICE3D_STEPS = (
    *("parascopy_depth", "dysgu_ins_dup", "dellycnv", "cnvpytor", "insurveyor"),
    *("scramble", "insertions_union", "parascopy_loci", "gauchian", "kilda"),
    "aldy_cyp2a6",
)
# Hand-added exons in the middle of two fixture deletions no other union call overlaps
# (tests/fixtures/hg002/hg002.callers.*): Dysgu's 2,030 bp het one at chr20:2014856 that
# DELLY confirms (depth ratio 0.47), and Dysgu's 174 bp het one at chr20:17481850 alone.
EXONS = {
    "BRCA1": [Exon("BRCA1", "ENST00000357654.9", 7, "chr20", 2015800, 2015900)],
    "CFTR": [Exon("CFTR", "ENST00000003084.11", 3, "chr20", 17481900, 17481950)],
}
FIXTURE_UNION = 57 + 43 + 30  # tp-base, fn and fp of hg002.callers.*
FIXTURE_SHORT = 3  # under 50 bp: Dysgu's 49 bp records (SVLEN 50) in tp-base and fn


def _deletion_rows(
    g: Genome, run_id: str
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]], dict[str, Any]]:
    """The run's deletion statements by id, its deletion call rows by vid, and its inputs."""
    with connect(g.store) as con:
        cursor = con.execute(
            "SELECT * FROM statement WHERE run_id = ? AND source = 'deletions'",
            [run_id],
        )
        names = [d[0] for d in cursor.description]
        statements = {
            row[1]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()
        }
        cursor = con.execute(
            "SELECT * FROM call WHERE run_id = ? AND vid LIKE '%:<DEL>'", [run_id]
        )
        names = [d[0] for d in cursor.description]
        calls = {
            row[names.index("vid")]: dict(zip(names, row, strict=True))
            for row in cursor.fetchall()
        }
        inputs = _inputs_json(con, run_id)
    return statements, calls, inputs


def _assert_deletions_unassessed(genome: Genome, run_id: str, reason: str) -> None:
    """One `deletions:screen` statement, unassessed with `reason`; no union matched; the
    accounting balanced; the run's `tools` without the three deletion tools."""
    statements, calls, inputs = _deletion_rows(genome, run_id)
    assert list(statements) == ["deletions:screen"]
    screen = statements["deletions:screen"]
    assert (screen["kind"], screen["reason"], screen["gene"]) == (
        "unassessed",
        reason,
        None,
    )
    assert calls == {}
    assert "truvari" not in genome.deletions.runs
    assert "deletions_union" not in inputs["steps"]
    accounting = _accounting(genome, run_id, "deletions")
    assert (accounting.input_rows, accounting.unassessed, accounting.findings) == (
        1,
        1,
        0,
    )
    assert not DELETION_TOOLS & set(_tools(genome, run_id))
    assert inputs["deletions"]["reason"] == reason


def test_uncalibrated_deletions_unassessed(genome: Genome) -> None:
    shutil.rmtree(deletions_dir(genome.lab.paths.lab))
    run_id = _run_genome(genome)
    _assert_deletions_unassessed(genome, run_id, NOT_CALIBRATED)
    assert NOT_CALIBRATED.endswith("run `genokeeper calibrate`")
    assert genome.deletions.runs == {}
    assert _deletion_rows(genome, run_id)[2]["deletions"]["table"] is None


def test_stale_deletion_table_unassessed(genome: Genome) -> None:
    write_deletions(genome.lab, dysgu="genomelab/dysgu:1.9.0-other")
    run_id = _run_genome(genome)
    reason = "stale: dysgu changed; run `genokeeper calibrate`"
    _assert_deletions_unassessed(genome, run_id, reason)
    assert genome.deletions.runs == {}


def test_sv_files_missing_unassessed(genome: Genome) -> None:
    shutil.rmtree(genome.lab.paths.lab / "sv")
    run_id = _run_genome(genome)
    _assert_deletions_unassessed(
        genome, run_id, "gnomAD-SV missing: run scripts/bootstrap.sh --with-sv"
    )
    assert genome.deletions.runs == {}


def test_screen_runs_without_hg002_truth(genome: Genome) -> None:
    """Calibrated, then HG002's v5.0q SV truth removed: the run still screens."""
    shutil.rmtree(truth_files(genome.lab.paths.lab)[0].parent, ignore_errors=True)
    screen = _deletion_rows(genome, _run_genome(genome))[0]["deletions:screen"]
    assert screen["kind"] == "non_finding"
    assert screen["reason"].startswith(f"screened: {FIXTURE_UNION} deletions called")


def test_female_run_shards_no_chry(genome: Genome) -> None:
    """A female run gives the callers no chrY shard, and the screen still runs."""
    run_id = _run_genome(genome, sex="female")
    assert genome.deletions.shards == {
        "delly": {"chr1", "chrX"},
        "dysgu": {"chr1", "chrX"},
    }
    screen = _deletion_rows(genome, run_id)[0]["deletions:screen"]
    assert screen["kind"] == "non_finding"


def test_caller_failure_unassessed(genome: Genome) -> None:
    """DELLY failing leaves the screen unassessed with its failure: the table was measured
    with both callers' classes, so Dysgu's calls alone are not shown under them."""
    genome.deletions.fail.add("delly")
    run_id = _run_genome(genome)
    _assert_deletions_unassessed(genome, run_id, "delly failed: exit 1: delly crashed")
    steps = _deletion_rows(genome, run_id)[2]["steps"]
    assert "dysgu" in steps and "delly" not in steps


def test_calibrated_deletions_statements(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With a current table, the sv files and two hand-added exons, the callers run on the
    run's CRAM: the both-callers deletion over BRCA1's exon and the Dysgu-only one over
    CFTR's are findings at validity limited, each call row with its class's chance of being
    false from the table and its warnings, hidden by `query`'s defaults; the rest dropped
    by reason; the screen's own statement a non-finding; the accounting balanced."""
    monkeypatch.setattr(pipeline, "load_exons", lambda: EXONS)
    # The run hands the table's own cautions (BR-002) to each call's warnings.
    monkeypatch.setattr(
        DeletionTable, "cautions", lambda _, d, c: [f"caution {c.band}"]
    )
    run_id = _run_genome(genome)
    statements, calls, inputs = _deletion_rows(genome, run_id)
    brca1 = statements["deletions:BRCA1:chr20:2014856-2016886"]
    cftr = statements["deletions:CFTR:chr20:17481850-17482024"]
    assert set(statements) == {
        "deletions:screen",
        *(s["statement_id"] for s in (brca1, cftr)),
    }
    for statement, topic in ((brca1, "risks"), (cftr, "carrier")):
        assert (statement["kind"], statement["validity"], statement["topic"]) == (
            "finding",
            "limited",
            topic,
        )
        assert statement["actionable"] is False
    both = json.loads(brca1["claim"])["class"]
    alone = json.loads(cftr["claim"])["class"]
    assert both == {
        "band": "long",
        "evidence": "both+depth",
        "repeat": False,
        "masked": False,
    }
    assert alone == {
        "band": "short",
        "evidence": "dysgu",
        "repeat": False,
        "masked": False,
    }
    for statement, cls in ((brca1, both), (cftr, alone)):
        call_class = CallClass(**cls)
        call = calls[statement["vid"]]
        rate = genome.deletion_table.p_wrong(call_class)
        assert rate is not None
        assert (
            call["p_wrong"],
            call["p_wrong_lo"],
            call["p_wrong_hi"],
        ) == pytest.approx(rate)
        assert call["warnings"] == [
            *warnings(call_class),
            f"caution {call_class.band}",
            CONSEQUENCE_WARNING,
        ]
        assert call["level"] == "low"
    assert set(calls) == {brca1["vid"], cftr["vid"]}
    screen = statements["deletions:screen"]
    assert (screen["kind"], screen["validity"], screen["reason"]) == (
        "non_finding",
        "limited",
        (
            f"screened: {FIXTURE_UNION} deletions called, 2 shown; a miss is possible:"
            " on HG002's libraries the screen found 9 of 10 (pcr), 9 of 10 (pcrfree)"
            " rare long and 9 of 10 (pcr), 9 of 10 (pcrfree) rare short truth deletions"
        ),
    )
    assert genome.deletion_table.recall() == {
        band: {"pcr": (9, 10), "pcrfree": (9, 10)} for band in ("long", "short")
    }
    accounting = _accounting(genome, run_id, "deletions")
    assert accounting.dropped == {
        SHORT: FIXTURE_SHORT,
        NO_EXON: FIXTURE_UNION - FIXTURE_SHORT - 2,
        COMMON: 0,
        SAME_POSITION: 0,
    }
    assert (accounting.findings, accounting.non_findings, accounting.unassessed) == (
        2,
        1,
        0,
    )
    assert accounting.input_rows == 3 + FIXTURE_UNION - 2
    with connect(genome.store) as con:
        default = query(con, run_id, Filters())
        loosest = query(
            con,
            run_id,
            Filters(min_reliability=Level.LOW, min_validity=Validity.LIMITED),
        )
    sources = [{s["source"] for s in f["sources"]} for f in default.findings]
    assert not any("deletions" in found for found in sources)
    listed = [
        f
        for f in loosest.findings
        if any(s["source"] == "deletions" for s in f["sources"])
    ]
    assert len(listed) == 2
    tools_json = _tools(genome, run_id)
    assert {tools_json[name] for name in DELETION_TOOLS} == {"2.7.0", "1.9.0", "5.5.0"}
    assert inputs["deletions"] == {
        "table": genome.deletion_table.versions,
        "reason": None,
    }


# --- slice 3d: copy number, insertions, paralog loci ----------------------------------------

# Hand-added exons, on chr20 where no fixture event lies unless said: CFTR's under Dysgu's and
# INSurVeyor's 228 bp insertion at chr20:3053197 (tests/fixtures/hg002/hg002.ins.tp-base),
# under a depth-only loss and under a Scramble-only mobile element; DPYD's under a depth loss.
EXONS_3D = {
    **EXONS,
    "CFTR": [
        *EXONS["CFTR"],
        Exon("CFTR", "ENST00000003084.11", 12, "chr20", 3053150, 3053250),
        Exon("CFTR", "ENST00000003084.11", 14, "chr20", 55001000, 55001200),
        Exon("CFTR", "ENST00000003084.11", 16, "chr20", 57000000, 57000300),
    ],
    "DPYD": [Exon("DPYD", "ENST00000370192.8", 5, "chr20", 56001000, 56001200)],
}
SLICE3D_IMAGES = ("cnvpytor", "insurveyor", "scramble", "gauchian", "kilda")
SLICE3D_KEYS = {"delly cnv", *SLICE3D_IMAGES, *(f"{n}_image" for n in SLICE3D_IMAGES)}
# Every slice-3d container a calibrated run starts but the depth step 3a shares.
SLICE3D_CONTAINERS = {
    *("delly cnv", "cnvpytor", "insurveyor", "scramble", "truvari"),
    *("parascopy cn", "gauchian", "kilda", "aldy"),
}


def _cnv_record(chrom: str, start: int, end: int, cn: int, vid: str) -> str:
    """A `delly cnv` PASS record over 0-based `[start, end)` at copy number `cn`."""
    return (
        f"{chrom}\t{start + 1}\t{vid}\tN\t<CNV>\t60\tPASS"
        f"\tIMPRECISE;SVTYPE=CNV;END={end}\tGT:CN:FT\t./.:{cn}:PASS"
    )


def _mei(chrom: str, at: int) -> str:
    """A row of Scramble's MEI table: an Alu at 0-based `at`."""
    return f"{chrom}:{at}\talu\tPlus\t7\t130.0\t100\t100\tGCC\tleft\t2\t64" + (
        "\tNone Found" * 4 + "\tNA"
    )


def _slice3d(
    g: Genome, run_id: str
) -> tuple[dict[str, dict[str, Any]], dict[str, SimpleNamespace], dict[str, Any]]:
    """The run's slice-3d statements by id, each tool's accounting, and the run's inputs."""
    statements = {
        i: s for i, s in _statements(g, run_id).items() if s["source"] in SLICE3D_TOOLS
    }
    accounting = {tool: _accounting(g, run_id, tool) for tool in SLICE3D_TOOLS}
    for tool, a in accounting.items():
        counted = a.findings + a.watch + a.non_findings + a.unassessed
        assert a.input_rows == counted + sum(a.dropped.values()), tool
    with connect(g.store) as con:
        inputs = _inputs_json(con, run_id)
    return statements, accounting, inputs


def _assert_part_unassessed(
    statements: Mapping[str, Mapping[str, Any]], part: str, reason: str
) -> None:
    """Slice 3d's `part` (cnv, insertions, paralogs with lpa) unassessed with `reason`."""
    ids = {
        "cnv": ["cnv:screen"],
        "insertions": ["insertions:screen"],
        "paralogs": [*(f"paralogs:{n}" for n in LOCI), "lpa:LPA"],
    }[part]
    for i in ids:
        # RHD's reason carries its note whatever it says (`sources.paralogs.RHD_NOTE`).
        assert statements[i]["kind"] == "unassessed"
        assert statements[i]["reason"].startswith(reason), i
    others = [s for s in statements.values() if s["source"] == ids[0].split(":")[0]]
    if part != "paralogs":
        assert [s["statement_id"] for s in others] == ids


def test_uncalibrated_slice3d_unassessed(genome: Genome) -> None:
    """With no slice-3d table no slice-3d step runs: the depth and insertion screens, every
    paralog locus and Lp(a) are unassessed with the calibrate command, each tool's
    accounting balances, and the run records none of 3d's tools."""
    for name in ("cnv", "insertions", "paralogs"):
        shutil.rmtree(tables.folder(genome.lab.paths.lab, name))
    run_id = _run_genome(genome)
    statements, accounting, inputs = _slice3d(genome, run_id)
    for part in ("cnv", "insertions", "paralogs"):
        _assert_part_unassessed(statements, part, NOT_CALIBRATED)
    not_covered = {f"paralogs:{gene}" for gene in PARALOG_NOT_COVERED}
    assert {i for i, s in statements.items() if s["kind"] != "unassessed"} == set()
    assert set(statements) == {
        "cnv:screen",
        "insertions:screen",
        "lpa:LPA",
        *(f"paralogs:{n}" for n in LOCI),
        *not_covered,
    }
    assert accounting["cnv"].input_rows == accounting["insertions"].input_rows == 1
    # Only the Parascopy depth SMN1 reads: 3a's, under its own table.
    assert genome.slice3d.runs == {"parascopy depth": 1}
    assert not set(SLICE3D_STEPS) - {"parascopy_depth"} & set(inputs["steps"])
    assert not SLICE3D_KEYS & set(_tools(genome, run_id))
    assert inputs["slice3d"]["cnv"] == {
        "table": None,
        "reason": NOT_CALIBRATED,
        "ran": False,
    }


def test_deletion_screen_independent_of_cnv(genome: Genome) -> None:
    """3b's deletions never wait on the `cnv` table (C6): without it they are graded as
    ever; with it and `delly cnv` failing, the depth screen names the failure and 3b's
    deletions are graded as ever."""
    shutil.rmtree(tables.folder(genome.lab.paths.lab, "cnv"))
    run_id = _run_genome(genome)
    statements, _, _ = _slice3d(genome, run_id)
    _assert_part_unassessed(statements, "cnv", NOT_CALIBRATED)
    screen = _deletion_rows(genome, run_id)[0]["deletions:screen"]
    assert screen["kind"] == "non_finding"
    assert "delly cnv" not in genome.slice3d.runs
    write_cnv(genome.lab)
    genome.slice3d.fail.add("delly cnv")
    run_id = _run_genome(genome)
    statements, _, inputs = _slice3d(genome, run_id)
    reason = statements["cnv:screen"]["reason"]
    assert reason == "dellycnv failed: exit 1: delly crashed"
    _assert_part_unassessed(statements, "cnv", reason)
    screen = _deletion_rows(genome, run_id)[0]["deletions:screen"]
    assert screen["kind"] == "non_finding"
    assert "dellycnv" not in inputs["steps"] and "deletions_union" in inputs["steps"]


def test_joined_depth_support(genome: Genome, monkeypatch: pytest.MonkeyPatch) -> None:
    """A depth loss over a shown 3b deletion (HG002's own `delly cnv` call at chr20:2014856,
    CN 1) is no statement of its own: 3b's deletion says read depth also supports it, its
    class and chance of being false unchanged."""
    monkeypatch.setattr(pipeline, "load_exons", lambda: EXONS)
    run_id = _run_genome(genome)
    deletions, calls, _ = _deletion_rows(genome, run_id)
    brca1 = deletions["deletions:BRCA1:chr20:2014856-2016886"]
    assert "; read depth also supports it (DELLY cnv); " in brca1["reason"]
    cls = json.loads(brca1["claim"])["class"]
    assert cls["evidence"] == "both+depth"
    rate = genome.deletion_table.p_wrong(CallClass(**cls))
    assert rate is not None
    call = calls[brca1["vid"]]
    assert (call["p_wrong"], call["p_wrong_lo"], call["p_wrong_hi"]) == (
        pytest.approx(rate)
    )
    statements, accounting, _ = _slice3d(genome, run_id)
    assert not any("2014856" in i for i in statements)
    assert accounting["cnv"].dropped["a 3b deletion"] == 1


def test_parascopy_depth_once(genome: Genome) -> None:
    """SMN1 and the paralog loci read one Parascopy depth of the CRAM: one container."""
    run_id = _run_genome(genome)
    assert genome.slice3d.runs["parascopy depth"] == 1
    assert genome.stage3.runs["parascopy"] == 1  # SMN1's cn-using
    assert genome.slice3d.runs["parascopy cn"] == 1  # the six paralog models
    statements, _, inputs = _slice3d(genome, run_id)
    assert statements["paralogs:STRC"]["kind"] == "non_finding"
    assert _statements(genome, run_id)["smn:SMN1"]["kind"] == "non_finding"
    assert not inputs["steps"]["parascopy_depth"]["reused"]


def test_dysgu_raw_missing_fails_insertions(genome: Genome) -> None:
    """3b's Dysgu step folder without its raw calls: the depth and insertion screens, which
    both read them, are unassessed naming it; the paralog loci graded as calibrated."""
    first = _run_genome(genome)
    with connect(genome.store) as con:
        steps = _inputs_json(con, first)["steps"]
    (Path(steps["dysgu"]["path"]).parent / "dysgu.vcf.gz").unlink()
    shutil.rmtree(Path(steps["dysgu_ins_dup"]["path"]).parent)
    run_id = _run_genome(genome)
    statements, _, _ = _slice3d(genome, run_id)
    reason = statements["cnv:screen"]["reason"]
    assert reason.startswith("dysgu_ins_dup failed: ")
    assert reason.endswith(
        "Dysgu's raw calls are gone: remove its step folder and rerun"
    )
    _assert_part_unassessed(statements, "cnv", reason)
    _assert_part_unassessed(statements, "insertions", reason)
    assert statements["paralogs:STRC"]["kind"] == "non_finding"
    assert statements["lpa:LPA"]["kind"] == "finding"


def test_stale_cnv_table_only_cnv_unassessed(genome: Genome) -> None:
    """A `cnv` table measured with another DELLY image leaves the depth screen unassessed
    naming the key; insertions, the paralog loci and 3b's deletions as calibrated."""
    write_cnv(genome.lab, delly="genomelab/delly:2.7.0-other")
    run_id = _run_genome(genome)
    statements, _, inputs = _slice3d(genome, run_id)
    reason = "stale: delly changed; run `genokeeper calibrate`"
    _assert_part_unassessed(statements, "cnv", reason)
    assert statements["insertions:screen"]["kind"] == "non_finding"
    assert statements["paralogs:STRC"]["kind"] == "non_finding"
    assert _deletion_rows(genome, run_id)[0]["deletions:screen"]["kind"] == (
        "non_finding"
    )
    assert "delly cnv" not in genome.slice3d.runs
    assert inputs["slice3d"]["cnv"]["reason"] == reason
    assert inputs["slice3d"]["insertions"]["ran"]


def test_calibrated_slice3d_statements(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With current tables and hand-added exons: a depth-only loss, an insertion labelled
    by a Scramble Alu and a Scramble-only Alu are findings; a depth loss over 3b's BRCA1
    deletion only adds its depth note; GJB6-D13S1830 is a finding on GJB2; a DPYD loss is a
    drugs finding; STRC two copies a non-finding; the genes not covered unassessed. Each
    finding is hidden by `query`'s defaults and listed at the loosest filters, and the run
    records 3d's tools."""
    monkeypatch.setattr(pipeline, "load_exons", lambda: EXONS_3D)
    genome.slice3d.dellycnv += [
        _cnv_record("chr20", 55000000, 55003000, 1, "CNV_CFTR"),
        _cnv_record("chr20", 56000000, 56003000, 1, "CNV_DPYD"),
        _cnv_record("chr13", 20223037, 20531806, 1, "CNV_GJB6"),
    ]
    genome.slice3d.meis += [_mei("chr20", 3053300), _mei("chr20", 57000100)]
    run_id = _run_genome(genome)
    statements, accounting, _ = _slice3d(genome, run_id)
    depth = statements["cnv:CFTR:chr20:55000000-55003000"]
    dpyd = statements["cnv:DPYD:chr20:56000000-56003000"]
    labelled = statements["insertions:CFTR:chr20:3053197-3053198"]
    alone = statements["insertions:CFTR:chr20:57000100-57000101"]
    named = statements["cnv:named:GJB6-D13S1830"]
    findings = [depth, dpyd, labelled, alone, named]
    for s in findings:
        assert (s["kind"], s["validity"], s["actionable"]) == (
            "finding",
            "limited",
            False,
        )
    assert (depth["topic"], dpyd["topic"], named["topic"]) == (
        "carrier",
        "drugs",
        "carrier",
    )
    assert DPYD_NOTE in dpyd["reason"]
    assert named["gene"] == "GJB2"
    assert labelled["reason"].startswith("possible Alu insertion")
    assert json.loads(alone["claim"])["class"]["band"] == "mei"
    assert not any("2014856" in i for i in statements)
    brca1 = _deletion_rows(genome, run_id)[0]["deletions:BRCA1:chr20:2014856-2016886"]
    assert "read depth also supports it" in brca1["reason"]
    assert statements["paralogs:STRC"]["kind"] == "non_finding"
    for gene, (topic, reason) in PARALOG_NOT_COVERED.items():
        s = statements[f"paralogs:{gene}"]
        assert (s["kind"], s["topic"], s["reason"]) == ("unassessed", topic, reason)
    assert not any(s["kind"] == "watch" for s in statements.values())
    assert accounting["cnv"].findings == 3
    with connect(genome.store) as con:
        default = query(con, run_id, Filters())
        loosest = query(
            con,
            run_id,
            Filters(min_reliability=Level.LOW, min_validity=Validity.LIMITED),
        )
    ids = {s["statement_id"] for s in findings}
    assert not ids & {f["id"] for f in default.findings}
    assert ids <= {f["id"] for f in loosest.findings}
    # An information copy number is no pharmacogene: no drug, no recommendation.
    (cyp2a6,) = [f for f in loosest.findings if f["id"] == "paralogs:CYP2A6"]
    assert "drugs" not in cyp2a6["claim"] and cyp2a6["actionable"] is False
    assert cyp2a6["gene_call"]["diplotype"] == "CYP2A6 2"
    tools_json = _tools(genome, run_id)
    assert SLICE3D_KEYS <= set(tools_json)
    assert tools_json["delly cnv"] == "2.7.0"


def test_cnvpytor_dropped_not_run(genome: Genome) -> None:
    """A `cnv` table that dropped CNVpytor: no CNVpytor container, every depth call DELLY's
    alone, no large event, and the `cnv:screen` claim says why."""
    write_cnv(genome.lab, cnvpytor="dropped")
    chrom, start, _ = LARGE_REGION
    genome.slice3d.dellycnv.append(_cnv_record(chrom, start, start + 600_000, 1, "BIG"))
    genome.slice3d.cnvpytor.append(
        f"deletion\t{chrom}:{start + 1}-{start + 600_000}\t600000\t0.5\t0.0"
    )
    run_id = _run_genome(genome)
    assert "cnvpytor" not in genome.slice3d.runs
    statements, accounting, inputs = _slice3d(genome, run_id)
    claim = json.loads(statements["cnv:screen"]["claim"])
    assert claim["notes"] == [pipeline.CNVPYTOR_DROPPED]
    assert not any(i.startswith("cnv:large:") for i in statements)
    assert accounting["cnv"].dropped[LARGE_FILTERED] >= 1
    assert "cnvpytor" not in inputs["steps"] and "cnvpytor" not in _tools(
        genome, run_id
    )
    with connect(genome.store) as con:
        callers = [
            json.loads(row[0])
            for row in con.execute(
                "SELECT claim->'callers' FROM statement WHERE run_id = ?"
                " AND source = 'cnv' AND statement_id <> 'cnv:screen'",
                [run_id],
            ).fetchall()
        ]
    assert all(set(c) == {"dellycnv"} for c in callers)


def test_large_event_not_calibrated(
    genome: Genome, capsys: pytest.CaptureFixture[str]
) -> None:
    """A 600 kb loss both depth callers make over a ClinGen region is shown, not
    calibrated: its call's chance of being false not measured, its level low."""
    chrom, start, _ = LARGE_REGION
    end = start + 600_000
    genome.slice3d.dellycnv.append(_cnv_record(chrom, start, end, 1, "BIG"))
    genome.slice3d.cnvpytor.append(
        f"deletion\t{chrom}:{start + 1}-{end}\t600000\t0.5\t0.0"
    )
    run_id = _run_genome(genome)
    statements, _, _ = _slice3d(genome, run_id)
    large = statements[f"cnv:large:{chrom}:{start}-{end}"]
    assert (large["kind"], large["topic"]) == ("finding", "risks")
    with connect(genome.store) as con:
        row = con.execute(
            "SELECT p_wrong, level FROM call WHERE run_id = ? AND vid = ?",
            [run_id, large["vid"]],
        ).fetchone()
    assert row == (None, "low")
    db = str(genome.store)
    args = [
        "query",
        "--db",
        db,
        "--min-validity",
        "limited",
        "--min-reliability",
        "low",
    ]
    assert cli.main([*args, "--text", "cnv:large", "--json"]) == 0
    (found,) = json.loads(capsys.readouterr().out)["findings"]
    assert found["reliability"]["p_wrong"] is None
    assert cli.main([*args, "--text", "cnv:large"]) == 0
    assert "p_wrong n/a (not measured)" in capsys.readouterr().out


def test_clingen_regions_missing_withholds_large(genome: Genome) -> None:
    """Without ClinGen's region lists the large-event rule is unavailable: a 600 kb loss
    both depth callers make is not shown, the `cnv:screen` claim says why, and the rest of
    the depth screen runs."""
    shutil.rmtree(genome.lab.paths.lab / "clingen-regions")
    chrom, start, _ = LARGE_REGION
    end = start + 600_000
    genome.slice3d.dellycnv.append(_cnv_record(chrom, start, end, 1, "BIG"))
    genome.slice3d.cnvpytor.append(
        f"deletion\t{chrom}:{start + 1}-{end}\t600000\t0.5\t0.0"
    )
    run_id = _run_genome(genome)
    statements, accounting, inputs = _slice3d(genome, run_id)
    screen = statements["cnv:screen"]
    assert screen["kind"] == "non_finding"
    assert json.loads(screen["claim"])["large_withheld"] == dict.fromkeys(
        ("loss_large", "gain_large"), pipeline.REGIONS_MISSING
    )
    assert pipeline.REGIONS_MISSING == (
        "large events not screened: ClinGen region lists missing; run"
        " scripts/bootstrap.sh"
    )
    assert not any(i.startswith("cnv:large:") for i in statements)
    # Over no region and no dosage-sensitive gene, it is no large event: no exon either.
    assert accounting["cnv"].dropped["no listed exon"] >= 1
    assert inputs["slice3d"]["clingen_regions"] == {"reason": pipeline.REGIONS_MISSING}


@pytest.mark.parametrize("tool", ["insurveyor", "scramble"])
def test_insertion_caller_failure_unassessed(genome: Genome, tool: str) -> None:
    """INSurVeyor or Scramble exiting non-zero leaves the insertion screen unassessed with
    its failure; every other part as calibrated."""
    genome.slice3d.fail.add(tool)
    run_id = _run_genome(genome)
    statements, _, inputs = _slice3d(genome, run_id)
    reason = f"{tool} failed: exit 1: {tool} crashed"
    _assert_part_unassessed(statements, "insertions", reason)
    assert statements["cnv:screen"]["kind"] == "non_finding"
    assert statements["paralogs:STRC"]["kind"] == "non_finding"
    assert statements["lpa:LPA"]["kind"] == "finding"
    assert tool not in inputs["steps"] and "insertions_union" not in inputs["steps"]
    assert "insurveyor" not in _tools(genome, run_id)


def test_gba1_contradiction_applied(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A ClinVar statement on a call carrying its allele inside GBA1's exons 9-11 (here moved
    over the mini chr1's first claims) where Gauchian is wild type is unassessed, its call's
    genotype unknown and every statement at its vid with it; a reference call there stays
    as it was; without Gauchian nothing changes. The store stays consistent: the run checks
    its calls and `query` reads why each is unassessed."""
    spans = {**load_spans(), GBA1_HOMOLOGY[0]: ("chr1", 900, 4100)}
    monkeypatch.setattr(pipeline, "load_spans", lambda: spans)
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        het = _row(con, run_id, BRCA1_HET[0])
        ref = _row(con, run_id, REFCALL[0])
        listed = query(con, run_id, LOOSEST_FILTERS, list_unassessed=True).to_json()
        at_vid = con.execute(
            "SELECT DISTINCT kind FROM statement WHERE run_id = ? AND vid = ("
            "SELECT vid FROM statement WHERE run_id = ? AND statement_id = ?)",
            [run_id, run_id, f"clinvar:{BRCA1_HET[0]}"],
        ).fetchall()
    assert (het["kind"], het["reason"], het["dosage"]) == (
        "unassessed",
        GBA1_CONTRADICTED,
        "unknown",
    )
    assert het["warnings"][0] == f"genotype unknown: {GBA1_CONTRADICTED}"
    assert at_vid == [("unassessed",)]
    assert (ref["kind"], ref["dosage"]) == ("non_finding", "0")
    (site,) = [s for s in listed["unassessed_sites"] if s["id"] == "clinvar:8001"]
    assert site["validity"]["reason"] == GBA1_CONTRADICTED
    # Gauchian's call contradicts nothing where the calibration withheld GBA1, nor
    # without the table.
    write_paralogs(genome.lab, withheld={"GBA1": "not covered: GBA1 not wild type"})
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        het = _row(con, run_id, BRCA1_HET[0])
    assert (het["kind"], het["dosage"]) == ("finding", "1")
    shutil.rmtree(tables.folder(genome.lab.paths.lab, "paralogs"))
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        het = _row(con, run_id, BRCA1_HET[0])
    assert (het["kind"], het["dosage"]) == ("finding", "1")


def test_female_no_haploid(genome: Genome) -> None:
    """A female run gives `delly cnv` no `--haploid-chr` and no chrY shard."""
    run_id = _run_genome(genome, sex="female")
    assert genome.slice3d.haploid == {None}
    assert genome.slice3d.shards == {"chr1", "chrX"}
    statements, _, _ = _slice3d(genome, run_id)
    assert statements["cnv:screen"]["kind"] == "non_finding"


def _gene_calls(con: duckdb.DuckDBPyConnection, run_id: str) -> dict[str, Any]:
    cursor = con.execute(
        "SELECT g.gene, g.diplotype, g.callers, g.agreed, g.level, g.warnings, s.kind,"
        " s.topic, s.validity, s.actionable"
        " FROM gene_call g JOIN statement s ON s.run_id = g.run_id"
        " AND s.statement_id = 'pharmcat:' || g.gene WHERE g.run_id = ?",
        [run_id],
    )
    names = [d[0] for d in cursor.description]
    return {row[0]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()}


def test_run_writes_gene_calls(genome: Genome) -> None:
    """Every run calls the pharmacogenes: PharmCAT on our calls at each of its positions
    (graded: a low or unknown call is a no-call) with the outside callers' calls, and each
    gene result is stored as a gene call with its gene-level statement under drugs.
    PharmGKB's annotations in genes PharmCAT gives no drug are located at their ClinVar
    allele and graded on our call there. PharmCAT's HTML report stays under reports/."""
    run_id = _run_genome(genome)
    # PharmCAT's input: every position, from our calls; the outside calls as assembled.
    assert [line.split("\t")[-1] for line in genome.pgx.vcf] == PGX_GENOTYPES
    assert genome.pgx.outside == OUTSIDE_LINES
    with connect(genome.store) as con:
        calls = _gene_calls(con, run_id)
        accounting = {
            row[0]: row[1:]
            for row in con.execute(
                "SELECT tool, input_rows, findings, watch, non_findings, unassessed,"
                " dropped FROM accounting WHERE run_id = ?",
                [run_id],
            ).fetchall()
        }
        pharmgkb = dict(
            con.execute(
                "SELECT statement_id, kind FROM statement"
                " WHERE run_id = ? AND source = 'pharmgkb'",
                [run_id],
            ).fetchall()
        )
        claims = [
            row[0]
            for row in con.execute(
                "SELECT claim FROM statement WHERE run_id = ? AND topic = 'drugs'",
                [run_id],
            ).fetchall()
        ]
        inputs = _inputs_json(con, run_id)
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
    assert set(calls) == {"CYP2C19", "CYP2D6", "DPYD", "HLA-A", "HLA-B", "F2"}
    # F2: a gene PharmCAT calls with no drug (unannotated), PharmCAT alone.
    assert (calls["F2"]["kind"], calls["F2"]["level"]) == ("finding", "medium")
    cyp2d6_call = calls["CYP2D6"]
    assert (cyp2d6_call["agreed"], cyp2d6_call["level"]) == (True, "high")
    assert set(json.loads(cyp2d6_call["callers"])) == {"BCyrius", "Aldy"}
    assert (calls["HLA-B"]["agreed"], calls["HLA-B"]["level"]) == (True, "high")
    assert calls["DPYD"]["level"] == "medium"  # PharmCAT alone, nothing missing
    assert calls["CYP2C19"]["level"] == "low"  # a defining position missing
    assert (calls["HLA-A"]["kind"], calls["HLA-A"]["level"]) == ("unassessed", "low")
    assert {c["topic"] for c in calls.values()} == {"drugs"}
    assert calls["CYP2C19"]["actionable"] is True  # CPIC A: clopidogrel
    # The accounting balances for each pharmacogenomic source.
    pharmcat_counts = accounting["pharmcat"]
    assert pharmcat_counts[:5] == (7, 5, 0, 0, 2)  # HLA-A and the G6PD gap
    pgkb = accounting["pharmgkb"]
    reasons = {k: v for k, v in json.loads(pgkb[5]).items() if v}
    assert pgkb[:5] == (7, 1, 0, 1, 1)
    assert reasons == {
        "level_1a": 1,
        "pharmcat_gene": 1,
        "not_an_rsid": 1,
        "rsid_not_in_clinvar": 1,
    }
    assert pharmgkb == {
        "pharmgkb:701": "finding",
        "pharmgkb:702": "non_finding",
        "pharmgkb:703": "unassessed",
    }
    # No dosing text in any stored drug claim.
    assert claims and not [c for c in claims if pharmcat.DOSING.search(c)]
    # The HTML report stays beside the JSON one, under reports/.
    report = Path(inputs["steps"]["pharmcat"]["path"])
    assert report.with_suffix(".html").is_file()
    assert report.is_relative_to(genome.lab.paths.reports)
    # The run records its slice-1b inputs and versions.
    lab = genome.lab.paths.lab
    pgx_inputs = inputs["pgx"]
    assert pgx_inputs["positions"] == {
        "path": str(positions_path(lab)),
        "sha256": sha256_file(positions_path(lab)),
    }
    assert pgx_inputs["t1k_index"]["path"] == str(t1k_index(lab))
    assert set(pgx_inputs["cpic"]) == set(CPIC_FILES.values())
    assert pgx_inputs["pharmgkb"]["created"] == PHARMGKB_DATE
    assert tools_row is not None
    tools_json = json.loads(tools_row[0])
    assert tools_json["pharmcat"] == PHARMCAT_VERSION
    assert tools_json["pharmgkb_clinical_annotations"] == PHARMGKB_DATE
    assert tools_json["pharmcat_image"] == IMAGES["pharmcat"]


@pytest.mark.parametrize(
    ("caller", "named", "genes"),
    [
        ("bcyrius", "BCyrius", {"CYP2D6"}),
        ("aldy", "Aldy", {"CYP2D6"}),  # and CYP2B6, which the report does not list
        ("t1k", "T1K", {"HLA-A", "HLA-B"}),
    ],
)
def test_a_failing_caller_leaves_its_genes_unassessed(
    genome: Genome, caller: str, named: str, genes: set[str]
) -> None:
    """A caller exiting non-zero on a sample's CRAM, run after run, leaves the genes it serves unassessed gene calls naming the failure; the run stores
    ClinVar's findings and the other genes, and keeps nothing of the failed step, so the
    rerun runs the caller again."""
    genome.pgx.fail[caller] = 1
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        calls = _gene_calls(con, run_id)
        (clinvar,) = con.execute(
            "SELECT count(*) FROM statement WHERE run_id = ? AND source = 'clinvar'",
            [run_id],
        ).fetchone() or (0,)
    for gene in genes:
        assert (calls[gene]["kind"], calls[gene]["level"]) == ("unassessed", "low")
        failure = f"{named} failed (exit 1: {caller} crashed)"
        assert failure in calls[gene]["warnings"][-1]
        assert json.loads(calls[gene]["callers"]) == {named: "failed"}
        assert gene not in (genome.pgx.outside or "")
    assert calls["CYP2C19"]["kind"] == "finding" and clinvar > 0
    assert caller != "aldy" or "CYP2B6" not in (genome.pgx.outside or "")
    _run_genome(genome)
    assert genome.pgx.runs[caller] == 2


@pytest.mark.parametrize("failure", ["pharmcat", "image"])
def test_pharmcat_failing_or_an_image_missing_refuses_in_one_line(
    genome: Genome, failure: str
) -> None:
    if failure == "pharmcat":
        genome.pgx.fail["pharmcat"] = 1
    else:
        genome.pgx.missing.add("optitype")
    with pytest.raises(pipeline.RunRefused) as refused:
        _run_genome(genome)
    assert len(str(refused.value).splitlines()) == 1
    assert ("pharmcat crashed" if failure == "pharmcat" else "No such image") in str(
        refused.value
    )


@pytest.mark.parametrize(
    ("missing", "named"),
    [
        ("positions", "PharmCAT 3.4.0 positions"),
        ("cpic", "CPIC"),
        ("pharmgkb", "ClinPGx clinical annotations"),
        ("t1k", "T1K HLA index"),
    ],
)
def test_missing_pgx_inputs_refused_before_any_step(
    genome: Genome, missing: str, named: str
) -> None:
    """A slice-1b input the bootstrap supplies refuses the run in one line naming it and
    the bootstrap, before any step."""
    lab = genome.lab.paths.lab
    if missing == "positions":
        positions_path(lab).unlink()
    elif missing == "cpic":
        shutil.rmtree(lab / "cpic")
    elif missing == "pharmgkb":
        shutil.rmtree(lab / "pharmgkb")
    else:
        t1k_index(lab).unlink()
    _refused_before_any_step(
        genome, rf"missing input: .*{named}.*scripts/bootstrap\.sh"
    )


def test_pgx_inputs_take_the_newest_folder_holding_either_download(
    genome: Genome,
) -> None:
    """The ClinPGx download a run reads is the newest dated folder's, the live
    `summaryAnnotations.zip` or the frozen `clinicalAnnotations.zip` an older bootstrap
    left; in one folder the live file wins."""
    lab = genome.lab.paths.lab
    pharmgkb = lab / "pharmgkb"
    live = pharmgkb / PHARMGKB_DATE / "summaryAnnotations.zip"
    assert read_pgx_inputs(lab).pharmgkb == live
    write_annotations(
        pharmgkb / "20240101" / "clinicalAnnotations.zip",
        [],
        created="2024-01-01",
        layout="clinical",
    )
    assert read_pgx_inputs(lab).pharmgkb == live
    frozen = write_annotations(
        pharmgkb / "20270101" / "clinicalAnnotations.zip",
        [],
        created="2027-01-01",
        layout="clinical",
    )
    found = read_pgx_inputs(lab)
    assert (found.pharmgkb, found.annotations.created) == (frozen, "20270101")
    newer = write_annotations(
        pharmgkb / "20270101" / "summaryAnnotations.zip", [], created="2027-01-01"
    )
    assert read_pgx_inputs(lab).pharmgkb == newer


def test_run_writes_research_associations(genome: Genome) -> None:
    """Every run states the GWAS Catalog's single-SNP associations on our calls at their
    panel alleles: a reference risk allele on a hom-ref call is a finding with two copies
    (`check_calls` passes it), the ALT there a non-finding, an unknown call unassessed;
    each limited, the topic from its trait; the rest counted; the release recorded."""
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        rows = {
            row[0]: row[1:]
            for row in con.execute(
                "SELECT s.statement_id, s.kind, s.gene, s.topic, s.validity, s.claim,"
                " c.row_type FROM statement s LEFT JOIN call c"
                " ON c.run_id = s.run_id AND c.vid = s.vid"
                " WHERE s.run_id = ? AND s.source = 'gwas_catalog'",
                [run_id],
            ).fetchall()
        }
        accounting = con.execute(
            "SELECT input_rows, findings, non_findings, unassessed, dropped"
            " FROM accounting WHERE run_id = ? AND tool = 'gwas_catalog'",
            [run_id],
        ).fetchone()
        inputs = _inputs_json(con, run_id)
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
    found = {
        sid: (kind, gene, topic, validity, json.loads(claim)["risk_allele_copies"], rt)
        for sid, (kind, gene, topic, validity, claim, rt) in rows.items()
    }
    assert found == {
        "gwas:rs20070-C:MONDO_0005148": (
            *("finding", "PPARG", "risks", "limited", 2, "ref_block"),
        ),
        "gwas:rs20070-A:EFO_0004339": (
            *("non_finding", "HMGA2", "research", "limited", 0, "ref_block"),
        ),
        "gwas:rs6000-A:EFO_0004611": (
            *("unassessed", "APOB", "wellness", "limited", None, "none"),
        ),
    }
    assert accounting is not None and accounting[:4] == (len(GWAS_ROWS), 1, 1, 1)
    assert {k: v for k, v in json.loads(accounting[4]).items() if v} == {
        "duplicate_association": 1,
        "not_in_panel": 1,
        "p_value_above_5e-8": 1,
    }
    catalog = genome.lab.paths.lab / "gwas" / GWAS_RELEASE / GWAS_DOWNLOAD
    variants = pgsc.panel_pvar(genome.lab.paths.lab)
    assert inputs["gwas"] == {
        "path": str(catalog),
        "sha256": sha256_file(catalog),
        "release": GWAS_RELEASE,
        "panel_variants": {"path": str(variants), "sha256": sha256_file(variants)},
    }
    assert tools_row is not None
    assert json.loads(tools_row[0])["gwas_catalog"] == GWAS_RELEASE


@pytest.mark.parametrize("missing", ["GWAS Catalog associations", "panel variants"])
def test_missing_gwas_inputs_refused_before_any_step(
    genome: Genome, missing: str
) -> None:
    if missing == "panel variants":
        pgsc.panel_pvar(genome.lab.paths.lab).unlink()
    else:
        shutil.rmtree(genome.lab.paths.lab / "gwas")
    _refused_before_any_step(
        genome, rf"missing input: no .*{missing} .*scripts/bootstrap\.sh"
    )


def test_run_phases_and_states_an_apoe_double_heterozygote(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """APOE from the run's graded ClinVar claims at its two SNPs (here two of the mini
    release's hets stand in for rs429358 C and rs7412 T): a double heterozygote joins the
    pairs phased from the reads, and its phase (in trans) makes it ε2/ε4, with a gene
    call."""
    e4, e2 = (
        ("chr1", pos, ref, _alt(ref))
        for pos in (BRCA1_HET[1], NO_INHERITANCE[1])
        for ref in [genome_bases(genome.lab.fasta, "chr1", pos)]
    )
    for module in (pipeline, pgsc):
        monkeypatch.setattr(module, "APOE_E4", e4)
        monkeypatch.setattr(module, "APOE_E2", e2)
    monkeypatch.setattr(prs_source, "_SNPS", {e4: "rs429358", e2: "rs7412"})
    run_id = _run_genome(genome)
    assert ("chr1", BRCA1_HET[1]) in genome.whatshap.sites
    assert ("chr1", NO_INHERITANCE[1]) in genome.whatshap.sites
    with connect(genome.store) as con:
        row = con.execute(
            "SELECT s.kind, s.claim, g.diplotype FROM statement s LEFT JOIN gene_call g"
            " ON g.run_id = s.run_id AND g.gene = s.gene"
            " WHERE s.run_id = ? AND s.statement_id = 'apoe:APOE'",
            [run_id],
        ).fetchone()
    assert row is not None
    assert (row[0], json.loads(row[1])["diplotype"], row[2]) == (
        "finding",
        "ε2/ε4",
        "ε2/ε4",
    )


def test_run_writes_polygenic_scores(genome: Genome) -> None:
    """Every run scores the chosen PGS with pgsc_calc on our calls: its target VCF holds
    each score and ancestry position once, reference and no-call rows explicit; each score
    is a finding with its percentile, one below 0.75 overlap unassessed with its rate; the
    APOE result needs ClinVar's claims at its two SNPs (the mini release has none)."""
    genome.nextflow.rates["PGS000903"] = 0.6
    run_id = _run_genome(genome)
    assert [line.split("\t")[-1] for line in genome.nextflow.vcf] == PRS_GENOTYPES
    with connect(genome.store) as con:
        rows = {
            row[0]: row[1:]
            for row in con.execute(
                "SELECT source_id, kind, validity, reason, effect FROM statement"
                " WHERE run_id = ? AND source = 'pgsc_calc'",
                [run_id],
            ).fetchall()
        }
        accounting = con.execute(
            "SELECT input_rows, findings, unassessed, dropped FROM accounting"
            " WHERE run_id = ? AND tool = 'pgsc_calc'",
            [run_id],
        ).fetchone()
        inputs = _inputs_json(con, run_id)
        tools_row = con.execute(
            "SELECT tools FROM run WHERE run_id = ?", [run_id]
        ).fetchone()
    scores = load_scores()
    not_covered = load_not_covered()
    assert set(rows) == {s.pgs_id for s in scores} | {n.slug for n in not_covered}
    kind, validity, reason, effect = rows["PGS000903"]
    assert (kind, reason) == ("unassessed", "score overlap 0.60, below 0.75")
    kind, validity, reason, effect = rows["PGS000018"]
    assert (kind, validity) == ("finding", "moderate")  # an EUR evaluation, EUR closest
    assert json.loads(effect)["score"]["percentile"] == 90.0
    assert accounting is not None
    assert accounting[:3] == (
        len(scores) + 1 + len(not_covered),
        len(scores) - 1,
        1 + len(not_covered),
    )
    assert json.loads(accounting[3]) == {"apoe_no_clinvar_claim": 1}
    assert tools_row is not None
    tools_json = json.loads(tools_row[0])
    assert (
        tools_json["pgsc_calc"] == f"{pgsc.PGSC_CALC_VERSION} ({pgsc.PGSC_CALC_COMMIT})"
    )
    assert tools_json["nextflow"] == pgsc.NEXTFLOW_VERSION
    assert tools_json["pgsc_plink2_image"] == IMAGES["pgsc_plink2"]
    lab = genome.lab.paths.lab
    assert inputs["prs"]["scores"]["PGS000018"] == {
        "path": str(pgsc.scoring_file(lab, "PGS000018")),
        "sha256": sha256_file(pgsc.scoring_file(lab, "PGS000018")),
    }
    # pgsc_calc's launch dir, work dir and results lie in its step, under reports/work.
    step = Path(inputs["steps"]["pgsc"]["path"])
    assert step.is_relative_to(genome.lab.paths.reports / "work" / "pgsc")
    # A task launched from an unpinned image refuses the run in one line; nothing of the
    # step is kept.
    genome.nextflow.image = None
    shutil.rmtree(step.parent)
    with pytest.raises(pipeline.RunRefused, match="not one of the pinned") as refused:
        _run_genome(genome)
    assert len(str(refused.value).splitlines()) == 1
    assert not step.parent.exists()
    # A failed Nextflow run keeps its stdout (with its log) under the run's work dir, in
    # one place each failure replaces, and the one-line refusal names it.
    genome.nextflow.image, genome.nextflow.fail = PINNED["plink2"], "ERROR ~ Missing"
    kept = genome.lab.paths.reports / "work" / "pgsc-failed"
    kept.mkdir()
    (kept / "stale").touch()
    with pytest.raises(pipeline.RunRefused, match="exit status 1") as refused:
        _run_genome(genome)
    assert str(refused.value).endswith(str(kept / "nextflow.stdout"))
    assert len(str(refused.value).splitlines()) == 1
    assert (kept / "nextflow.stdout").read_text() == "ERROR ~ Missing"
    assert {p.name for p in kept.iterdir()} == {".nextflow.log", "nextflow.stdout"}
    assert not step.parent.exists()


@pytest.mark.parametrize(
    "missing",
    ["scoring file", "reference panel", "ancestry positions", "Nextflow"],
)
def test_missing_prs_inputs_refused_before_any_step(
    genome: Genome, missing: str
) -> None:
    """A slice-1d input the bootstrap supplies refuses the run in one line naming it and
    the bootstrap, before any step."""
    lab = genome.lab.paths.lab
    {
        "scoring file": pgsc.scoring_file(lab, "PGS000018"),
        "reference panel": pgsc.panel_path(lab),
        "ancestry positions": pgsc.ancestry_pvar(lab),
        "Nextflow": pgsc.nextflow_path(lab),
    }[missing].unlink()
    _refused_before_any_step(
        genome, rf"missing input: .*{missing}.*scripts/bootstrap\.sh"
    )


def test_run_exports_its_snapshot_and_query_carries_the_contract(
    genome: Genome,
) -> None:
    run_id = _run_genome(genome)
    snapshot = genome.lab.paths.reports / "snapshots" / run_id
    assert sorted(p.name for p in snapshot.iterdir()) == sorted(
        f"{t}.parquet"
        for t in (*TABLES, "meta")  # meta: the store's schema version
    )
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    # The run's own snapshot is one `history` reads: every finding carries its dates.
    found = dates(genome.lab.paths.reports / "snapshots")
    with connect(genome.store) as con:
        result = query(con, None, loosest, list_unassessed=True, dates=found).to_json()
        # A research association comes back only to a narrowed query: here by its id.
        gwas = query(con, None, replace(loosest, text="gwas:")).to_json()
        result["findings"] += gwas["findings"]
        result["unassessed"]["total"] += gwas["unassessed"]["total"]
    # ClinVar's two, PharmGKB's no-call, PharmCAT's HLA-A without a call and its G6PD gap,
    # a GWAS association at the site with no record (6); AR withheld by the stage-3
    # table; 3c's sex-chromosome dosage statement under the shipped, not yet calibrated
    # [dosage], one more than the slice-3e plan counted; seven of the thirteen slice-3e
    # loci (FMR1 and AFF2 by the X dosage, which that [dosage] leaves unmeasured; NOP56
    # on too few spanning reads; PPP2R2B, JPH3, ZFHX3 and THAP11 with two alleles of one
    # length); the 22 not-covered loci: 6 + 1 + 1 + 7 + 22; each condition no polygenic
    # score covers; slice 3d's genes not covered (`sources.paralogs.NOT_COVERED`); and
    # slice 3f's 78 without a slice 3f table: 15 HLA (its three disease associations are
    # a research layer), HHV-6, ROH, 48 blood groups, the two gated ancestry statements
    # and the eleven traits (no site of theirs on the mini contigs).
    uncovered = len(load_not_covered()) + len(PARALOG_NOT_COVERED)
    assert result["unassessed"]["total"] == 37 + 78 + uncovered
    # No source gives an absolute effect.
    effect_keys = {p for p in CONTRACT if p.startswith("effect.absolute.")}
    finding_only = {"finding"}  # `show`'s wrapper key, not part of `query`
    # Every source's keys appear: ClinVar's, PharmCAT's gene results and PharmGKB's.
    # SMN1's two copies here are a non-finding, which `query` does not return (`show` does);
    # `test_skill_contract` pins its claim keys.
    smn_keys = {"claim.smn1", "claim.smn2", "claim.g27134tg", "claim.callers"}
    # Slice 3d's statements (an Alu over ARSA's exon, the paralog loci, Lp(a)) carry claim
    # keys the contract does not pin, as 3b's deletions' do not: the skill words them.
    keyed = {
        **result,
        **{
            part: [
                f
                for f in result[part]
                if not {s["source"] for s in f["sources"]} & SLICE3D_TOOLS
            ]
            for part in ("findings", "unassessed_sites")
        },
    }
    # A query without `--gene` prints no coverage.
    # Slice 3f's statements here are unassessed but the ungated similarity (R23): none
    # else writes the keys only a result has.
    similarity = {"claim.population", "claim.low_confidence", "claim.probabilities"}
    expected = CONTRACT - COVERAGE_CONTRACT - effect_keys
    expected -= SLICE3F_FINDING_CONTRACT - similarity
    expected -= finding_only | (smn_keys - {"claim.callers"})
    assert key_paths(keyed) == expected
    assert PGX_CONTRACT <= key_paths(result)
    assert {f["sources"][0]["version"] for f in keyed["findings"]} == {
        CLINVAR_RELEASE,
        PHARMCAT_VERSION,
        PHARMGKB_DATE,
        pgsc.PGSC_CALC_VERSION,
        f"pgsc_calc {pgsc.PGSC_CALC_VERSION}",  # ancestry:similarity, ungated (R23)
        haplogroup.HAPLOGREP_VERSION,
        GWAS_RELEASE,
        vep.VEP_VERSION,
        YLEAF_VERSION,
    }
    # The run's versions say it was graded with the FAF table (the skill's marker).
    assert result["versions"]["gnomad_faf"].startswith("4.1, floor ")


def test_qc_failure_stops_before_calling(genome: Genome) -> None:
    failing = replace(MINI_QC, min_mean_depth=1000.0)
    with pytest.raises(RunRefused, match=r"QC gate: mean_depth .* < 1000"):
        _run_genome(genome, qc_thresholds=failing)
    assert genome.deepvariant.calls == 0
    assert not genome.store.exists()
    # The finished steps stay cached: a rerun with another threshold aligns nothing again.
    assert len(list((genome.lab.paths.reports / "work" / "align").iterdir())) == 1
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        steps = _inputs_json(con, run_id)["steps"]
    assert steps["align"]["reused"] and steps["qc"]["reused"]
    assert genome.verifybamid.calls == 1


def test_calls_outside_the_gate_refuse_after_calling(genome: Genome) -> None:
    """Ti/Tv and het/hom are judged once DeepVariant has called: a refusal there stores
    nothing and runs no step after it."""
    failing = replace(MINI_QC, min_titv=1e9)
    with pytest.raises(RunRefused, match=r"QC gate: titv .* < 1e\+09"):
        _run_genome(genome, qc_thresholds=failing)
    assert genome.deepvariant.calls == 1
    assert (genome.mutserve.runs, genome.whatshap.sites) == (0, [])
    assert not genome.store.exists()


def test_missing_inputs_refused_before_any_step(genome: Genome) -> None:
    """The lifecycle's missing inputs refuse in one line naming the bootstrap."""
    checkpoint = checkpoint_path(genome.lab.paths.lab)
    (checkpoint.parent / "example_info.json").unlink()
    with pytest.raises(RunRefused, match="T7 checkpoint.*bootstrap"):
        _run_genome(genome)
    assert not (genome.lab.paths.reports / "work" / "align").exists()


def test_run_refuses_without_faf_table(genome: Genome) -> None:
    shutil.rmtree(faf_folder(genome.lab.paths.lab))
    with pytest.raises(RunRefused, match="--with-gnomad") as refused:
        _run_genome(genome)
    assert str(faf_folder(genome.lab.paths.lab)) in str(refused.value)
    assert "scripts/bootstrap.sh" in str(refused.value)
    # Refused before any step ran.
    assert not (genome.lab.paths.reports / "work" / "align").exists()
    assert not genome.store.exists()


def _refused_before_any_step(genome: Genome, match: str) -> None:
    with pytest.raises(RunRefused, match=match):
        _run_genome(genome)
    assert not (genome.lab.paths.reports / "work" / "align").exists()
    assert not genome.store.exists()


@pytest.mark.parametrize("key", VERSION_KEYS)
def test_table_from_other_versions_refused(genome: Genome, key: str) -> None:
    """A table measured under another pipeline (another bwa, samtools, image, checkpoint,
    make_examples flags or reference) is refused, naming what differs and the command
    that measures it again."""
    folder = table_dir(genome.lab.paths.lab)
    table = CalibrationTable.read(folder)
    replace_versions = {**table.versions, key: f"other-{key}"}
    write_calibration_table(
        genome.lab.paths.lab,
        "small_variants",
        CalibrationTable(cells=table.cells, versions=replace_versions).write,
        replace_versions,
    )
    _refused_before_any_step(
        genome, rf"stale: {key} changed; run `genokeeper calibrate`"
    )


def test_a_replaced_checkpoint_refuses_the_table(genome: Genome) -> None:
    """The run computes its versions from what it runs, never from the table: a
    checkpoint file replaced after calibrating changes the run's checkpoint md5."""
    weights = checkpoint_path(genome.lab.paths.lab).parent / CHECKPOINT_FILES[2]
    weights.write_text("retrained\n")
    _refused_before_any_step(genome, "checkpoint_md5.*genokeeper calibrate")


def test_run_refuses_without_calibration(genome: Genome) -> None:
    """Without a table the run aligns and measures QC, recording the sample's bases
    (what `genokeeper calibrate` downsamples HG003 to), and refuses before calling."""
    shutil.rmtree(table_dir(genome.lab.paths.lab))
    with pytest.raises(
        RunRefused,
        match=r"80000 bases.* is not calibrated: run `genokeeper calibrate`$",
    ):
        _run_genome(genome)
    assert genome.deepvariant.calls == 0 and not genome.store.exists()
    assert pipeline.recorded_bases(genome.lab.paths.reports) == 2 * 400 * READ
    write_calibration(genome.lab, genome.lab.paths.repo / "calibration-scratch")
    (layout(genome.lab.paths.lab)["truth"] / BENCH_BED).unlink()
    shutil.rmtree(genome.lab.paths.reports / "work" / "align")
    _refused_before_any_step(genome, "benchmark.*--with-calibration")


def test_open_cells_stale_is_not_a_refusal(genome: Genome) -> None:
    """HG002's cells measured under another reference are not used, and the run goes on
    recording why; current ones are opened beside HG003's table, haploid part first."""
    lab = genome.lab.paths.lab
    versions = genome.table.versions
    regions = {
        HAPLOID: [("chrX", 0, 1000)],
        ONLY: [("chr1", 0, 1000)],
        SHARED: [("chr1", 0, 10)],
    }
    stale = CalibrationTable(
        cells=genome.table.cells, versions={**versions, "reference_md5": "other"}
    )
    write_cells(cells_dir(lab), dict.fromkeys(regions, stale), regions)
    extra, reason = open_cells(lab, versions)
    assert extra == () and reason is not None and "reference_md5" in reason
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        inputs = _inputs_json(con, run_id)
    assert inputs["calibration"]["hg002_cells"] == {"reason": reason}
    write_cells(cells_dir(lab), dict.fromkeys(regions, genome.table), regions)
    extra, reason = open_cells(lab, versions)
    assert reason is None and [part.haploid for part in extra] == [True, False]
    assert extra[0].regions.flags_at("chrX", 500) and not extra[1].regions.flags_at(
        "chrX", 500
    )
    calibration = pipeline.open_calibration(lab, versions)
    assert calibration is not None and calibration.extra == extra
    assert calibration.record["hg002_cells"] == {
        part: {
            "path": str(cells_dir(lab) / part),
            "manifest_sha256": sha256_file(cells_dir(lab) / part / "manifest.json"),
            "regions_sha256": sha256_file(cells_dir(lab) / part / "regions.bed"),
        }
        for part in (HAPLOID, ONLY)
    }


def test_hg002_cells_stale_when_their_regions_change(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The cells' regions come from `[par]`, HG003's benchmark and v5.0q's truth: an edit
    to any of them stales the `hg002_cells` table."""
    lab = genome.lab.paths.lab
    vcf, v5 = calibrate_build.hg002_truth(lab)
    vcf.parent.mkdir(parents=True, exist_ok=True)
    for path in (vcf, v5):
        path.write_text("chr1\t0\t10\n")
    hg003 = layout(lab)["truth"] / BENCH_BED
    g = load_grading()
    edits = {
        "par": lambda: monkeypatch.setattr(
            calibrate_build,
            "load_grading",
            lambda: replace(g, par={**g.par, "chrX": [(1, 2)]}),
        ),
        "hg003_benchmark": lambda: hg003.write_text(hg003.read_text() + "chr2\t0\t5\n"),
        "v5_truth": lambda: vcf.write_text("chr1\t0\t11\n"),
        "v5_benchmark": lambda: v5.write_text("chr1\t0\t12\n"),
    }
    for key, edit in edits.items():
        versions = calibrate_build.hg002_cells_versions(lab)
        write_calibration_table(lab, "hg002_cells", lambda _: None, versions)
        assert isinstance(tables.open_table(lab, "hg002_cells"), tables.Current)
        edit()
        assert tables.open_table(lab, "hg002_cells") == tables.Stale(
            f"stale: {key} changed; run `genokeeper calibrate`"
        )


def test_calibrate_refuses_hg002_cells_without_pcrfree_calls(
    genome: Genome,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """With the small-variant table current and no HG002 PCR-free calls, `genokeeper
    calibrate` skips HG002's cells with one line naming what is missing and the command
    that supplies it, and exits 2."""
    lab = genome.lab.paths.lab
    monkeypatch.setattr(
        tables,
        "TABLES",
        {name: tables.TABLES[name] for name in ("small_variants", "hg002_cells")},
    )
    assert tables.calibrate(lab, only=None, report=False, threads=2) == 2
    captured = capsys.readouterr()
    assert captured.out.splitlines() == ["small_variants: current"]
    (line,) = captured.err.strip().splitlines()
    assert line.startswith(
        "genokeeper calibrate: hg002_cells: HG002 cells not measured: "
    )
    assert line.endswith("`python -m genokeeper.calibrate.hg002 align pcrfree`")
    assert not cells_dir(lab).exists()


def test_calibration_versions_come_from_the_pipeline(genome: Genome) -> None:
    """Each of the six versions is read from what the run runs: the tools on PATH, the
    image digest and DeepVariant flags in the code, and the checkpoint's and reference's
    files on disk. The shard count is not among the flags: it follows the machine."""
    lab = genome.lab.paths.lab
    tool_versions = pipeline.genome_tool_versions()
    versions = pipeline.calibration_versions(
        lab, tools=tool_versions, identities=genome.lab.paths.repo / "ids.json"
    )
    checkpoint = checkpoint_path(lab)

    def md5(path: Path) -> str:
        return hashlib.md5(path.read_bytes()).hexdigest()

    assert versions == {
        "bwa": tool_versions["bwa"],
        "samtools": tool_versions["samtools"],
        "deepvariant_image": IMAGES["deepvariant"],
        "checkpoint_md5": ",".join(
            f"{name}:{md5(checkpoint.parent / name)}" for name in CHECKPOINT_FILES
        ),
        "make_examples_flags": " ".join(CALLING_FLAGS),
        "reference_md5": md5(genome.lab.fasta),
    }
    for shards in (1, 8):
        args = deepvariant_args(
            Path("x.cram"),
            genome.lab.fasta,
            Path("/out"),
            checkpoint=checkpoint,
            male=True,
            par_bed=Path("par.bed"),
            threads=shards,
        )
        assert set(CALLING_FLAGS) <= set(args)
    assert "num_shards" not in versions["make_examples_flags"]


def test_calibrate_writes_a_table_the_run_accepts(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """`genokeeper calibrate` downsamples HG003 to the bases the user's run recorded,
    runs it through the run's steps, checks its inferred sex is male, scores it and
    writes the table with the versions the run computes: a run then accepts it. Its
    downsampling, sex and scoring (rtg) are faked; their own tests run them for real."""
    lab = genome.lab.paths.lab
    shutil.rmtree(table_dir(lab))
    where = layout(lab)
    (where["truth"] / TRUTH_VCF).write_text("placeholder\n")
    where["sdf"].mkdir()
    rtg = rtg_path(lab)
    rtg.parent.mkdir(parents=True)
    rtg.write_text("placeholder\n")
    where["fastq"].mkdir(parents=True)  # HG003's full FASTQ, as the bootstrap puts it
    full = [
        Path(shutil.copy(fastq, where["fastq"] / fastq.name))
        for fastq in (genome.lab.fastq1, genome.lab.fastq2)
    ]
    scored: list[Path] = []
    targets: list[int] = []
    cell = CellKey("variant", "snv", "easy", 50)

    def count_cells(gvcf: Path, *_args: object, **_kwargs: object) -> object:
        scored.append(gvcf)
        return {"chr1": {("genotype", cell): Counts(2, 1000, 3e4)}}, {"snp": 0.9971}

    def downsample(*pair: Path, target_bases: int) -> tuple[Path, ...]:
        targets.append(target_bases)
        pair[2].mkdir(parents=True)
        return tuple(Path(shutil.copy(f, pair[2] / f.name)) for f in pair[:2])

    def failed(*_pair: Path, target_bases: int) -> tuple[Path, ...]:
        raise ValueError(f"cannot downsample 1 bases to {target_bases}")

    def refused(match: str) -> None:
        assert (
            cli.main(["calibrate", "--only", "small_variants", "--threads", "2"]) == 2
        )
        err = capsys.readouterr().err
        assert err.startswith("genokeeper calibrate: ") and re.search(match, err), err
        assert len(err.strip().splitlines()) == 1 and not table_dir(lab).exists()

    monkeypatch.setattr(calibrate_build, "count_cells", count_cells)
    monkeypatch.setattr(calibrate_build, "downsample", failed)
    monkeypatch.setattr(calibrate_build, "load_paths", lambda: genome.lab.paths)
    monkeypatch.setattr(calibrate_build, "load_grading", _grading)
    refused("no run of your FASTQ has recorded its bases.*genokeeper run")
    with pytest.raises(
        RunRefused, match="is not calibrated: run `genokeeper calibrate`"
    ):
        _run_genome(genome)  # records the user's bases, then refuses
    refused("cannot downsample")
    assert all(f.is_file() for f in full)  # a failed downsample deletes nothing
    monkeypatch.setattr(calibrate_build, "downsample", downsample)
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: "female")
    refused("hg003's inferred sex is female, not male")
    # Both downsampled files landed: the full FASTQ is deleted, never read again.
    assert not any(f.exists() for f in full)
    assert genome.deepvariant.calls == 0
    shutil.rmtree(hg003_dir(lab) / "sex")
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: "male")
    assert (
        cli.main(
            ["calibrate", "--only", "small_variants", "--threads", "2", "--report"]
        )
        == 0
    )
    assert targets == [2 * 400 * READ]  # once: a rerun reuses the downsampled pair
    out = capsys.readouterr().out
    assert f"calibration table: {table_dir(lab)}" in out
    assert "held out, odd chromosomes -> even chromosomes" in out
    assert "SNP F1 0.9971 (the part stops below 0.995)\n" in out
    # HG003 went through the run's own call step, under the calibration folder.
    (gvcf,) = scored
    assert gvcf.is_relative_to(lab / "calibration" / "hg003")
    assert genome.deepvariant.calls == 1
    table = CalibrationTable.read(table_dir(lab))
    assert table.rate_for(cell).rate == 2 / 1000
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        row = _row(con, run_id, BRCA1_HET[0])
    assert row["p_wrong"] == 2 / 1000
    # Missing inputs refuse in one line, and write no table.
    shutil.rmtree(table_dir(lab))
    (where["truth"] / TRUTH_VCF).unlink()
    refused("--with-calibration")


def _hg003_record(genome: Genome, reference: str | None = None) -> Path:
    """The calibration inputs and an `hg003` record called on `reference` (default: this
    one); returns its gVCF."""
    lab = genome.lab.paths.lab
    shutil.rmtree(table_dir(lab))
    where = layout(lab)
    (where["truth"] / TRUTH_VCF).write_text("placeholder\n")
    where["sdf"].mkdir()
    rtg = rtg_path(lab)
    rtg.parent.mkdir(parents=True)
    rtg.write_text("placeholder\n")
    cram = hg003_dir(lab) / "align" / "k" / "sample.cram"
    gvcf = hg003_dir(lab) / "call" / "c" / "sample.g.vcf.gz"
    for p in (cram, gvcf):
        p.parent.mkdir(parents=True)
        p.write_bytes(b"")
    write_record(
        lab,
        "hg003",
        Record(
            cram=str(cram),
            key="k",
            metrics={"total_bases": 8e4},
            gvcf=str(gvcf),
            reference_sha256=reference or reference_sha256(lab),
            accepted=True,
            snp_f1=None,
            replaced=None,
        ),
    )
    return gvcf


def test_hg003_on_another_reference_refuses(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An `hg003` record on another reference belongs to `realign hg003`: `build_table`
    refuses naming it and aligns no pair, even one at the downsampled names."""
    _hg003_record(genome, reference="0" * 64)
    down = layout(genome.lab.paths.lab)["downsampled"]
    down.mkdir(parents=True, exist_ok=True)
    for name in ("x_1.fq.gz", "x_2.fq.gz"):
        (down / name).write_bytes(b"")

    def never(*_a: object, **_k: object) -> object:
        raise AssertionError("HG003 aligned outside realign")

    monkeypatch.setattr(calibrate_build, "run_sample", never)
    with pytest.raises(
        RunRefused,
        match=r"^HG003 is aligned on another reference; run"
        r" `python -m genokeeper\.calibrate\.hg002 realign hg003`$",
    ):
        calibrate_build.build_table(genome.lab.paths, grading=_grading(), threads=2)
    assert not table_dir(genome.lab.paths.lab).exists()


def test_hg003_scored_from_its_record(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A current `hg003` record (HG003 re-aligned from its CRAM and called) is scored from
    its gVCF: `build_table` runs no `run_sample` and reads no FASTQ."""
    lab = genome.lab.paths.lab
    gvcf = _hg003_record(genome)
    scored: list[Path] = []
    cell = CellKey("variant", "snv", "easy", 50)

    def count_cells(gvcf: Path, *_args: object, **_kwargs: object) -> object:
        scored.append(gvcf)
        return {"chr1": {("genotype", cell): Counts(2, 1000, 3e4)}}, {"snp": 0.9971}

    def never(*_a: object, **_k: object) -> object:
        raise AssertionError("HG003 aligned or downsampled again")

    monkeypatch.setattr(calibrate_build, "count_cells", count_cells)
    monkeypatch.setattr(calibrate_build, "run_sample", never)
    monkeypatch.setattr(calibrate_build, "downsample", never)
    table, _ = calibrate_build.build_table(
        genome.lab.paths, grading=_grading(), threads=2
    )
    assert scored == [gvcf]
    assert CalibrationTable.read(table_dir(lab)).rate_for(cell).rate == 2 / 1000
    assert table.rate_for(cell).rate == 2 / 1000


@pytest.mark.parametrize(
    "content",
    [
        None,
        '"GENE SYMBOL","MOI"\n"OTC","XL"\n',
        # A field over the csv module's limit raises csv.Error, not ValueError.
        f'"FILE CREATED: 2026-10-18"\n"{"x" * 200_000}"\n',
    ],
    ids=["absent", "no classification column", "csv error"],
)
def test_missing_or_unreadable_gene_validity_refuses_before_any_step(
    genome: Genome, content: str | None
) -> None:
    if content is None:
        genome.lab.gene_validity.unlink()
    else:
        genome.lab.gene_validity.write_text(content)
    with pytest.raises(RunRefused, match="gene-validity.csv"):
        _run_genome(genome)
    assert not (genome.lab.paths.reports / "work" / "align").exists()
    assert not genome.store.exists()


def test_missing_tool_refuses(
    genome: Genome, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PATH", str(tmp_path / "no-tools"))
    with pytest.raises(RunRefused, match="bcftools"):
        _run_genome(genome)


def test_a_call_mismatch_rolls_the_run_back(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A writer that broke the one-call-per-vid rule: every call the source writes is
    # unknown while its statements keep their grades.
    call_row = panel_source.call_row

    def unknown(run_id: str, vid: str, called: Any, reliability: Any) -> Any:
        return call_row(run_id, vid, replace(called, dosage="unknown"), reliability)

    monkeypatch.setattr(panel_source, "call_row", unknown)
    with pytest.raises(pipeline.CallMismatch, match=MISMATCH):
        _run_genome(genome)
    assert _counts(genome.store) == dict.fromkeys(TABLES, (0,))


def test_interrupted_run_leaves_no_rows(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    def boom(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("injected after the rows")

    # Runs after every row of the run is written.
    monkeypatch.setattr(pipeline, "check_calls", boom)
    with pytest.raises(RuntimeError, match="injected"):
        _run_genome(genome)
    assert _counts(genome.store) == dict.fromkeys(TABLES, (0,))
    assert not (genome.lab.paths.reports / "snapshots").exists()


def test_a_run_exports_an_earlier_runs_snapshot_a_failed_export_left_out(
    genome: Genome,
) -> None:
    snapshots = genome.lab.paths.reports / "snapshots"
    first = _run_genome(genome)
    shutil.rmtree(snapshots / first)  # as an export that failed after the commit
    second = _run_genome(genome)
    assert {p.name for p in snapshots.iterdir()} == {first, second}


def test_a_held_lock_refuses_and_a_run_removes_a_killed_runs_rows(
    genome: Genome,
) -> None:
    # A SIGKILL mid-flush leaves the user's rows as JSONL; the next run, holding the lock,
    # removes them, and a refused run (another holds the lock) touches nothing.
    work = genome.lab.paths.reports / "work"
    stale = work / ".store-rows-killed"
    stale.mkdir(parents=True)
    (stale / "call.jsonl").write_text('{"dosage": "1"}\n')
    with (work / ".run.lock").open("w") as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RunRefused, match="another run is under way"):
            _run_genome(genome)
    assert (stale / "call.jsonl").exists()
    assert not genome.store.exists()
    _run_genome(genome)  # released: the next run goes ahead
    assert not stale.exists()
    assert list(work.glob(".store-rows-*")) == []


@pytest.mark.parametrize("offset", [1, -1])
def test_a_store_of_another_schema_version_refuses_before_any_step(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, offset: int
) -> None:
    with connect(genome.store) as con:
        con.execute("UPDATE meta SET schema_version = schema_version + ?", [offset])
    before = genome.store.read_bytes()

    def no_step(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("a step ran before the store was checked")

    monkeypatch.setattr(pipeline, "_check_genome_tools", no_step)
    monkeypatch.setattr(pipeline, "_cached", no_step)
    with pytest.raises(RunRefused, match="schema version"):
        _run_genome(genome)
    assert not any((genome.lab.paths.reports / "work").glob("*/*"))
    assert genome.store.read_bytes() == before


# --- the command line ---------------------------------------------------------------------


def _cli(genome: Genome, monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """`genokeeper run`'s arguments, its paths and grading the mini lab's; the FASTQ pair,
    ClinGen, CGD and the FAF table are found where the bootstrap puts them."""
    monkeypatch.setattr(pipeline, "load_paths", lambda: genome.lab.paths)
    monkeypatch.setattr(pipeline, "load_grading", _grading)
    return ["run", "--clinvar38", str(genome.clinvar38), "--sex", "male"] + [
        "--threads",
        "2",
    ]


def test_cli_prints_the_run_and_its_accounting(
    genome: Genome,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = _cli(genome, monkeypatch)
    umask = os.umask(0o022)
    os.umask(umask)
    assert cli.main(args) == 0
    assert os.umask(umask) == umask  # restored for the caller
    out = capsys.readouterr().out
    with connect(genome.store) as con:
        row = con.execute("SELECT run_id FROM run").fetchone()
        counts = con.execute(
            "SELECT input_rows, findings, watch, non_findings, unassessed FROM accounting"
            " WHERE tool = 'clinvar'"
        ).fetchone()
    assert row is not None and counts is not None
    with connect(genome.store) as con:
        inputs = _inputs_json(con, row[0])
    assert f"run {row[0]}" in out
    # The unassessed count sits beside the others.
    assert "tool\tinput_rows\tfindings\twatch\tnon_findings\tunassessed\tdropped" in out
    assert f"clinvar\t{'\t'.join(map(str, counts))}\t" in out
    assert "\npharmcat\t" in out and "\npharmgkb\t" in out
    report = Path(inputs["steps"]["pharmcat"]["path"]).with_suffix(".html")
    assert out.count(f"PharmCAT report: {report}\n") == 1
    # The report's path holds a step key (hex), not the user's data.
    out = out.replace(str(report), "")
    assert "QC: " in out and "mean_depth" in out
    faf = FafTable.open(faf_folder(genome.lab.paths.lab)).describe()
    assert out.count(f"population frequency: {faf}\n") == 1
    # Everything the run wrote is the user's alone.
    written = [genome.lab.paths.reports, *genome.lab.paths.reports.rglob("*")]
    assert {oct(p.stat().st_mode & 0o077) for p in written if not p.is_symlink()} == {
        "0o0"
    }
    # Counts only: no genotype or position of the user's.
    assert "0/1" not in out
    assert str(BRCA1_HET[1]) not in out


@pytest.mark.parametrize(
    ("what", "message"),
    [
        ("tool", "bcftools is not on PATH"),
        ("gene validity", "gene-validity.csv"),
        ("faf table", "--with-gnomad"),
        ("calibration table", "genokeeper calibrate"),
        ("pharmgkb", "ClinPGx clinical annotations"),
        ("older store", "move it aside and rerun `genokeeper run`"),
    ],
)
def test_cli_exits_2(
    genome: Genome,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    what: str,
    message: str,
) -> None:
    args = _cli(genome, monkeypatch)
    before = None
    if what == "tool":
        monkeypatch.setenv("PATH", str(tmp_path / "no-tools"))
    elif what == "gene validity":
        genome.lab.gene_validity.unlink()
    elif what == "faf table":
        (faf_folder(genome.lab.paths.lab) / "faf-chr22.parquet").unlink()
    elif what == "calibration table":
        (table_dir(genome.lab.paths.lab) / "cells.tsv").unlink()
    elif what == "pharmgkb":
        shutil.rmtree(genome.lab.paths.lab / "pharmgkb")
    else:
        with connect(genome.store) as con:
            con.execute("UPDATE meta SET schema_version = schema_version - 1")
        before = genome.store.read_bytes()
    assert cli.main(args) == 2
    err = capsys.readouterr().err
    assert err.startswith("genokeeper run: ")
    assert message in err
    assert len(err.strip().splitlines()) == 1
    if before is None:
        assert not genome.store.exists()
    else:
        assert genome.store.read_bytes() == before  # refused untouched


def test_a_store_another_process_holds_refuses_in_one_line(
    genome: Genome,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    args = _cli(genome, monkeypatch)
    genome.store.parent.mkdir(parents=True)
    duckdb.connect(str(genome.store)).close()
    hold = "import duckdb, sys, time; c = duckdb.connect(sys.argv[1]); print('held', flush=True); time.sleep(60)"
    holder = subprocess.Popen(
        [sys.executable, "-c", hold, str(genome.store)],
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        assert holder.stdout is not None and holder.stdout.readline() == "held\n"
        assert cli.main(args) == 2
    finally:
        holder.kill()
        holder.wait()
    err = capsys.readouterr().err
    assert err.startswith("genokeeper run: ")
    assert "in use" in err
    assert len(err.strip().splitlines()) == 1


def test_inferred_sex_is_worded_as_the_x_copy_count(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """The run's sex is the X copy count's, never chrY's: in the run record and in
    `--sex`'s help."""
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: "male")
    run_id = _run_genome(genome, sex=None)
    with connect(genome.store) as con:
        decided_by = _inputs_json(con, run_id)["sex"]["decided_by"]
    assert decided_by.startswith("X copy count") and "chrY" not in decided_by
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["run", "--help"])
    shown = " ".join(capsys.readouterr().out.split())
    start = shown.rindex("--sex")  # its help, after the usage line
    sex = shown[start : shown.index("--threads", start)]
    assert "X copy count" in sex and "chrY" not in sex


def test_run_command_takes_the_one_fastq_pair(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "x_1.fq.gz").write_bytes(b"")
    with pytest.raises(RunRefused, match="expected two FASTQ.*--fastq1"):
        pipeline.fastq_pair(raw)
    (raw / "x_2.fq.gz").write_bytes(b"")
    assert pipeline.fastq_pair(raw) == (raw / "x_1.fq.gz", raw / "x_2.fq.gz")
    args = cli.build_parser().parse_args(["run", "--sex", "male", "--threads", "4"])
    assert args.handler is pipeline._genome_command
    assert (args.sex, args.threads, args.fastq1) == ("male", 4, None)


def test_fastq_hash_cached_by_inode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A FASTQ is identified by its full sha256, computed again only when its path, size,
    mtime or inode changed."""
    fastq = tmp_path / "reads_1.fq.gz"
    fastq.write_bytes(b"@r1\nACGT\n+\nIIII\n")
    cache = tmp_path / "work" / "identities.json"
    real = checksum.sha256_file
    hashed: list[Path] = []

    def counting(path: Path) -> str:
        hashed.append(path)
        return real(path)

    monkeypatch.setattr(checksum, "sha256_file", counting)
    first = sha256_cached(fastq, cache)
    assert first == real(fastq) and len(hashed) == 1
    assert sha256_cached(fastq, cache) == first and len(hashed) == 1

    # Another file of the same size and mtime moved in over it: a new inode.
    before = fastq.stat()
    other = tmp_path / "other.fq.gz"
    other.write_bytes(b"@r1\nTTTT\n+\nIIII\n")
    os.utime(other, ns=(before.st_atime_ns, before.st_mtime_ns))
    os.replace(other, fastq)
    assert fastq.stat().st_ino != before.st_ino
    assert fastq.stat().st_mtime_ns == before.st_mtime_ns
    second = sha256_cached(fastq, cache)
    assert second == real(fastq) != first and len(hashed) == 2

    # Rewritten in place: the same inode, another mtime.
    fastq.write_bytes(b"@r1\nGGGG\n+\nIIII\n")
    os.utime(fastq, ns=(before.st_atime_ns, before.st_mtime_ns + 10**9))
    assert sha256_cached(fastq, cache) == real(fastq) and len(hashed) == 3
    # The cache holds file metadata and digests only.
    entry = json.loads(cache.read_text())[str(fastq.resolve())]
    assert set(entry) == {"size", "mtime_ns", "inode", "sha256"}


# --- slice 3c: identity, coverage, dosage, the allele-fraction guard -------------------------

# The one GWAS SNP of the mini lab both sides can hold: rs20070 in the reference block (two
# reference copies, not A/T or C/G); rs6000 has no call.
IDENTITY_POS = 20070
NOT_DELIVERED = "no delivered genotype file in data/raw"


def _deliver(g: Genome, genotype: str, name: str = "genome.txt") -> Path:
    return write_delivered(
        g.lab.paths.raw / name, [(f"rs{IDENTITY_POS}", "1", IDENTITY_POS, genotype)]
    )


def _reference_pair(g: Genome) -> tuple[str, str]:
    ref = genome_bases(g.lab.fasta, "chr1", IDENTITY_POS)
    return ref, _alt(ref)


def _min_compared(n: int) -> Grading:
    g = _grading()
    return replace(g, identity=replace(g.identity, min_compared=n))


def test_identity_refuses_a_mismatched_delivery(genome: Genome) -> None:
    """A delivered file whose genotypes oppose the run's calls refuses before the claims
    step, and nothing reaches the store."""
    _, alt = _reference_pair(genome)
    delivered = _deliver(genome, alt * 2)
    with pytest.raises(RunRefused, match="does not match the delivered genotypes"):
        _run_genome(genome, grading=_min_compared(1), delivered=delivered)
    assert not (genome.lab.paths.reports / "work" / "claims").exists()
    assert not genome.store.exists()


def test_identity_recorded_when_matched(genome: Genome) -> None:
    ref, _ = _reference_pair(genome)
    delivered = _deliver(genome, ref * 2)
    run_id = _run_genome(genome, grading=_min_compared(1), delivered=delivered)
    assert _tools(genome, run_id)["identity"] == (
        f"{IDENTITY_TOOL}: 1 SNPs, concordance 1.0000"
    )
    with connect(genome.store) as con:
        identity = _inputs_json(con, run_id)["identity"]
    skipped = identity.pop("skipped")
    assert sum(skipped.values()) == 1  # rs6000
    assert identity == {
        "path": str(delivered),
        "sha256": sha256_file(delivered),
        "verdict": "matched",
        "compared": 1,
        "concordant": 1,
        "reason": None,
    }


def test_identity_not_checked_without_a_file(genome: Genome) -> None:
    run_id = _run_genome(genome)
    assert _tools(genome, run_id)["identity"] == (
        f"{IDENTITY_TOOL}: not checked: {NOT_DELIVERED}"
    )
    with connect(genome.store) as con:
        identity = _inputs_json(con, run_id)["identity"]
    assert identity == {"verdict": "not checked", "reason": NOT_DELIVERED}


def test_identity_not_checked_below_min_compared(genome: Genome) -> None:
    """A matching file at the shipped `min_compared`, which the mini lab cannot reach: not
    checked, and the run completes."""
    assert load_grading().identity.min_compared == 1000
    ref, _ = _reference_pair(genome)
    run_id = _run_genome(genome, delivered=_deliver(genome, ref * 2))
    assert _tools(genome, run_id)["identity"].startswith(
        f"{IDENTITY_TOOL}: not checked: too few common SNPs"
    )


def test_two_delivered_files_refuse(
    genome: Genome, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    ref, _ = _reference_pair(genome)
    _deliver(genome, ref * 2, "one.txt")
    _deliver(genome, ref * 2, "two.txt")
    monkeypatch.setattr(pipeline, "load_paths", lambda: genome.lab.paths)
    args = cli.build_parser().parse_args(["run"])
    assert args.genotypes is None
    assert pipeline._genome(args) == 2  # pyright: ignore[reportPrivateUsage]
    err = capsys.readouterr().err
    assert "--genotypes" in err and err.count("\n") == 1
    assert not genome.store.exists()


def test_genotypes_names_the_delivered_file(genome: Genome) -> None:
    args = cli.build_parser().parse_args(["run", "--genotypes", "/x/genome.txt"])
    assert args.genotypes == Path("/x/genome.txt")


def test_foundation_steps_write_coverage_and_dosage(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One gene's CDS on mini chr1 [1000, 2000) is one coverage row of 1,004 bases (its
    2 bp margins); without a foundation table the dosage statement is unassessed."""
    one = {"GENEC": [Exon("GENEC", "NM_000001.1", 1, "chr1", 1000, 2000)]}
    monkeypatch.setattr(gene_coverage, "load_cds", lambda: one)
    run_id = _run_genome(genome)  # no HG002 record and no foundation table
    with connect(genome.store) as con:
        rows = con.execute(
            "SELECT gene, bases, owner, reliable FROM coverage WHERE run_id = ?",
            [run_id],
        ).fetchall()
        inputs = _inputs_json(con, run_id)
    assert rows == [("GENEC", 1004, None, False)]  # about 3x: below 20x
    dosage = _statements(genome, run_id)[DOSAGE_STATEMENT]
    assert (dosage["kind"], dosage["reason"]) == ("unassessed", NOT_CALIBRATED)
    assert inputs["coverage"] == {} and inputs["dosage"] == {}
    assert not inputs["steps"]["coverage"]["reused"]


def test_male_x_linked_gene_measured_per_copy(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A male run measures a chrX non-PAR gene at `min_depth_per_copy` x 1 and an
    autosomal one at x 2: one mosdepth call per expected copy count."""
    one = {
        "GENEC": [Exon("GENEC", "NM_000001.1", 1, "chr1", 1000, 2000)],
        "GENEXL": [Exon("GENEXL", "NM_000003.1", 1, "chrX", 1000, 2000)],
    }
    monkeypatch.setattr(gene_coverage, "load_cds", lambda: one)
    measured: dict[str, int | None] = {}
    real = pipeline.mosdepth

    def spy(
        cram: Path, fasta: Path, regions: Path, prefix: Path, **kw: Any
    ) -> dict[str, tuple[float, float]]:
        for line in regions.read_text().splitlines():
            measured[line.split("\t")[3]] = kw["thresholds"]
        return real(cram, fasta, regions, prefix, **kw)

    monkeypatch.setattr(pipeline, "mosdepth", spy)
    _run_genome(genome, sex="male")
    per_copy = load_grading().coverage.min_depth_per_copy
    assert measured == {"GENEC": per_copy * 2, "GENEXL": per_copy}


def test_dosage_follows_the_foundation_table(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A current `foundation` table grades the dosage; after `[dosage] tolerance` changes,
    which stales the table, the statement is unassessed with the table's own reason."""
    lab = genome.lab.paths.lab
    versions = calibrate_foundation.foundation_versions(lab)
    write_calibration_table(lab, "foundation", lambda _: None, versions)
    dosage = _statements(genome, _run_genome(genome))[DOSAGE_STATEMENT]
    assert dosage["kind"] != "unassessed"
    g = load_grading()
    edited = replace(g, dosage=replace(g.dosage, tolerance=0.14))
    monkeypatch.setattr(calibrate_foundation, "load_grading", lambda: edited)
    run_id = _run_genome(genome, grading=replace(edited, qc=MINI_QC))
    dosage = _statements(genome, run_id)[DOSAGE_STATEMENT]
    assert (dosage["kind"], dosage["reason"]) == (
        "unassessed",
        "stale: dosage changed; run `genokeeper calibrate`",
    )


def test_a_failed_dosage_step_leaves_the_statement_unassessed(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    def failing(*_a: Any, **_k: Any) -> tuple[float, float]:
        raise tools.ToolFailed(1, ["mosdepth"], "", "boom")

    monkeypatch.setattr(pipeline, "measure_ratios", failing)
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        inputs = _inputs_json(con, run_id)
    assert "dosage" not in inputs["steps"]
    assert inputs["dosage"] == {"failed": "exit 1: boom"}
    dosage = _statements(genome, run_id)[DOSAGE_STATEMENT]
    assert dosage["kind"] == "unassessed" and "boom" in dosage["reason"]


def test_guarded_call_in_the_store(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A ClinVar-claimed heterozygote at VAF 0.2 inside a clonal gene's span: level low,
    the clonal warning."""
    monkeypatch.setattr(
        clonal,
        "load_spans",
        lambda: {
            gene: Span(gene, "NM_000002.1", "chr1", 900, 1100)
            for gene in load_grading().af_guard.clonal_genes
        },
    )
    pos = BRCA1_HET[1]
    low = genome_variant(
        genome.lab.fasta,
        "chr1",
        pos,
        [_alt(genome_bases(genome.lab.fasta, "chr1", pos))],
        gt="0/1",
        gq=50,
        ad=[24, 6],
    )
    records = genome.deepvariant.records
    records[:] = [low if r.split("\t")[1] == str(pos) else r for r in records]
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        row = _row(con, run_id, BRCA1_HET[0])
    assert row["level"] == "low"
    assert AF_CLONAL in row["warnings"]


def test_tools_record_the_reference(genome: Genome) -> None:
    run_id = _run_genome(genome)
    assert _tools(genome, run_id)["reference"] == MASKED_REFERENCE


def test_sex_step_key_holds_the_dosage_tolerance(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`sex_from_ratios` reads `[dosage] tolerance`: a changed tolerance rebuilds the sex
    step."""
    built: list[float] = []

    def infer(*_a: Any, grading: Grading, **_k: Any) -> str:
        built.append(grading.dosage.tolerance)
        return "male"

    monkeypatch.setattr(pipeline, "infer_sex_from_cram", infer)
    work = genome.lab.paths.reports / "work"
    aligned = pipeline._Step(key="k", path=work / "x.cram", reused=False)  # pyright: ignore[reportPrivateUsage]
    g = _grading()

    def step(grading: Grading) -> bool:
        return pipeline.sex_step(
            work,
            aligned=aligned,
            fasta=genome.lab.fasta,
            lab=genome.lab.paths.lab,
            versions={"mosdepth": "0.3"},
            grading=grading,
            threads=1,
        ).reused

    assert not step(g) and step(g)
    wider = replace(g, dosage=replace(g.dosage, tolerance=g.dosage.tolerance + 0.05))
    assert not step(wider)
    assert built == [g.dosage.tolerance, wider.dosage.tolerance]


# --- slice 3f: the typing sources, chrM, GenCC ----------------------------------------------

SLICE3F_TOOLS = (
    "hla",
    "mitomap",
    "hhv6",
    "roh",
    "rbceq2",
    "ancestry",
    "traits",
    "pharmgkb_research",
)
# The statements a slice 3f table gates, by prefix (`ancestry:similarity` aside, ruling
# R23); those of them that are never assessed carry their own reason.
GATED = ("hla:", "hhv6:", "roh:", "bloodgroup:", "ancestry:a", "ancestry:n")
NEVER_ASSESSED = {
    "hla:KIR": KIR_REASON,
    "bloodgroup:RH": RH_REASON,
    "bloodgroup:MNS": MNS_REASON,
    "bloodgroup:HPA": HPA_REASON,
}
GENCC_DATE = "20261025"
NO_GENCC = {"present": False, "reason": "no GenCC download: run scripts/bootstrap.sh"}
# The six DQ tag SNPs moved into the gVCF's reference block (the mini lab has no chr6):
# each called 0/0 there.
MINI_TAGS = tuple(range(BLOCK[0] + 11, BLOCK[0] + 17))


def _write_gencc(g: Genome, rows: list[tuple[str, str, str]]) -> Path:
    """A GenCC export of (gene, disease, classification title) rows in the lab."""
    path = g.lab.paths.lab / "gencc" / GENCC_DATE / GENCC_FILE
    path.parent.mkdir(parents=True)
    head = "gene_symbol\tdisease_curie\tclassification_title\tsubmitted_run_date\n"
    path.write_text(head + "".join(f"{r}\t{d}\t{t}\t2026-10-25\n" for r, d, t in rows))
    return path


def _mitomap(pos: int) -> Any:
    return next(v for v in load_mitomap() if v.pos == pos)


def _chrm_call(g: Genome, level: float, reference: float | None = None) -> None:
    """mutserve's calls: the ClinVar chrM claim's allele at `level` and, with `reference`,
    the second chrM claim's allele at that level; the mito step built again."""
    calls = []
    for claim, at in ((CHRM_CLAIM, level), (CHRM_REFERENCE, reference)):
        if at is not None:
            ref = genome_bases(g.lab.fasta, "chrM", claim[1])
            calls.append((claim[1], ref, _alt(ref), at))
    g.mutserve.calls = calls
    shutil.rmtree(g.lab.paths.reports / "work" / "mito", ignore_errors=True)


def _call_of(g: Genome, run_id: str, vid: str) -> dict[str, Any]:
    with connect(g.store) as con:
        cursor = con.execute(
            "SELECT dosage, vaf, warnings FROM call WHERE run_id = ? AND vid = ?",
            [run_id, vid],
        )
        row = cursor.fetchone()
        assert row is not None, f"no call at {vid}"
        return dict(zip([d[0] for d in cursor.description], row, strict=True))


def _flagged(rows: dict[str, dict[str, Any]]) -> set[str]:
    return {sid for sid, r in rows.items() if r["exclusion"]}


def test_uncalibrated_slice3f_unassessed(genome: Genome) -> None:
    """Without a slice 3f table the run finishes and no typing tool runs: every statement
    the table gates is unassessed with 3c's reason (KIR, RH, MNS and the platelet antigens
    with their own), PharmCAT is given no phase and UGT1A1*28 is not trusted. The ungated
    ones are written as without 3f's calibration: the traits, MITOMAP's and every chrM
    call at the config's floor, the scores' deciles and Z_norm2, GenCC's validity, the
    level-3 annotations as research, the G6PD gap; nothing is flagged an exclusion."""
    m591, m616 = _mitomap(591), _mitomap(616)
    genome.mutserve.calls = [
        *genome.mutserve.calls,
        (m591.pos, m591.ref, m591.alt, 0.06),
        (m616.pos, m616.ref, m616.alt, 0.15),
    ]
    _write_gencc(genome, [("GENEY", "MONDO:0000002", "Limited")])

    def labels(pos: int) -> tuple[str, ...]:
        ref = genome_bases(genome.lab.fasta, "chr1", pos)
        return (ref + ref, ref + _alt(ref), _alt(ref) * 2)

    write_annotations(
        genome.lab.paths.lab / "pharmgkb" / "20250706" / PHARMGKB_DOWNLOADS[0],
        [
            Ann("801", "rs8001", "TNF", "3", labels(BRCA1_HET[1]), drugs="etanercept"),
            Ann("802", "rs8010", "CYP2C19", "1B", labels(DRUG[1])),
        ],
    )
    run_id = _run_genome(genome)
    rows = _statements(genome, run_id)
    with connect(genome.store) as con:
        inputs = _inputs_json(con, run_id)

    gated = {sid: r for sid, r in rows.items() if sid.startswith(GATED)}
    assert len(gated) == 15 + 3 + 1 + 1 + 45 + 3 + 2
    assert rows["ancestry:similarity"]["reason"] != NOT_CALIBRATED
    for sid, row in gated.items():
        assert (sid, row["kind"]) == (sid, "unassessed")
        assert (sid, row["reason"]) == (sid, NEVER_ASSESSED.get(sid, NOT_CALIBRATED))
    assert genome.typing.runs == {}
    assert inputs["slice3f"] == {"table": None, "reason": NOT_CALIBRATED}
    assert inputs["pgx"]["phase"] == NOT_CALIBRATED
    assert inputs["pgx"]["ugt1a1_28_trusted"] is False
    assert all(line.split("\t")[8] == "GT" for line in genome.pgx.vcf)
    assert not {"hhv6", "roh", "bloodgroups", "admixture", "archaic"} & set(
        inputs["steps"]
    )

    # Ungated: the traits on our genotypes (every trait site is off the mini contigs or
    # past their ends).
    traits = {sid: r for sid, r in rows.items() if r["source"] == "traits"}
    assert set(traits) == {
        f"trait:{n}"
        for n in ("lactase", "ADH1B", "ACTN3", "TAS2R38", "CYP1A2")
        + ("FUT2", "ABCC11", "MC1R", "HERC2", "AMPD1", "CCR5")
    }
    assert all(r["reason"] != NOT_CALIBRATED for r in traits.values())
    # AMPD1's chr1 site, left out of the step, is read unknown (a missing allele).
    ampd1 = traits["trait:AMPD1"]
    assert (ampd1["kind"], ampd1["reason"]) == (
        "unassessed",
        "genotype unknown at rs17602729",
    )
    # MITOMAP and chrM at the config's floor, 0.10.
    below, above = rows["mitomap:591C>T"], rows["mitomap:616T>C"]
    assert (below["kind"], above["kind"]) == ("unassessed", "finding")
    warned = _call_of(genome, run_id, below["vid"])["warnings"]
    assert "genotype unknown: heteroplasmy at level 0.06, below the floor 0.1" in warned
    assert inputs["haplocheck"]["status"] == "not contaminated"
    # The scores' deciles and Z_norm2.
    scores = [
        r
        for r in rows.values()
        if r["source"] == "pgsc_calc" and r["kind"] == "finding"
    ]
    assert scores
    for score in scores:
        effect = json.loads(score["effect"])["score"]
        assert effect["decile"] == 10 and "z_norm2" in effect
    # GenCC's validity for a gene ClinGen does not curate, named in the reason.
    assert "gene validity Limited (GenCC)" in rows["clinvar:8003"]["reason"]
    assert inputs["gencc"]["present"] is True
    assert inputs["gencc"]["sha256"] == sha256_file(
        genome.lab.paths.lab / "gencc" / GENCC_DATE / GENCC_FILE
    )
    # The level-3 annotation is research, the 1B one the clinical layer's.
    research = rows["pharmgkb_research:801"]
    assert (research["kind"], research["topic"], research["validity"]) == (
        "finding",
        "research",
        "limited",
    )
    assert "pharmgkb:801" not in rows and "pharmgkb_research:802" not in rows
    assert _accounting(genome, run_id, "pharmgkb_research").input_rows == 1
    assert rows["pharmcat:G6PD-gap"]["kind"] == "unassessed"
    assert _flagged(
        {s: r for s, r in rows.items() if r["source"] in SLICE3F_TOOLS}
    ) == (set())


def test_haplocheck_flag_demotes_without_table(genome: Genome) -> None:
    """haplocheck runs in every run's mito step and its flag demotes without a table: a
    ClinVar chrM claim at level 0.40 is unassessed with the contamination reason, a 0.95
    call is unaffected."""
    _chrm_call(genome, 0.40, reference=0.95)
    genome.haplocheck.level = 0.08
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        demoted = _row(con, run_id, CHRM_CLAIM[0])
        kept = _row(con, run_id, CHRM_REFERENCE[0])
        inputs = _inputs_json(con, run_id)
    assert (demoted["kind"], demoted["dosage"]) == ("unassessed", "unknown")
    assert (
        "genotype unknown: haplocheck: possible contamination at level 0.08"
        in demoted["warnings"]
    )
    assert (kept["kind"], kept["dosage"]) == ("finding", "1")
    assert inputs["haplocheck"] == {
        "status": "contaminated",
        "level": 0.08,
        "major": "H5a7",
        "minor": "H5a7",
        "reason": None,
    }
    assert inputs["slice3f"]["table"] is None
    assert "haplocheck" in inputs["steps"] and genome.haplocheck.runs == 1


def test_slice3f_tool_failure_unassessed(genome: Genome) -> None:
    """haplocheck failing leaves its line in `run.inputs` and every chrM heteroplasmy
    unknown with it (ruling R20); RBCeq2 failing leaves each blood-group statement unassessed with its line
    while the other typing tools' statements are written; a rerun runs both again."""
    write_slice3f(genome.lab)
    _chrm_call(genome, 0.40)
    genome.haplocheck.fail = True
    genome.typing.fail.add("rbceq2")
    run_id = _run_genome(genome)
    rows = _statements(genome, run_id)
    with connect(genome.store) as con:
        chrm = _row(con, run_id, CHRM_CLAIM[0])
        inputs = _inputs_json(con, run_id)
    assert inputs["haplocheck"] == {
        "status": "failed",
        "level": None,
        "major": None,
        "minor": None,
        "reason": "haplocheck failed: exit 1: haplocheck crashed",
    }
    assert (chrm["kind"], chrm["dosage"]) == ("unassessed", "unknown")
    assert (
        "genotype unknown: haplocheck failed: exit 1: haplocheck crashed"
        in chrm["warnings"]
    )
    blood = {s: r for s, r in rows.items() if s.startswith("bloodgroup:")}
    for sid, row in blood.items():
        assert (sid, row["kind"]) == (sid, "unassessed")
        assert (sid, row["reason"]) == (
            sid,
            NEVER_ASSESSED.get(sid, "RBCeq2 failed: exit 1: rbceq2 crashed"),
        )
    assert rows["roh:summary"]["kind"] == "finding"
    assert rows["ancestry:admixture"]["kind"] == "finding"
    assert "bloodgroups" not in inputs["steps"] and "haplocheck" not in inputs["steps"]
    genome.haplocheck.fail = False
    genome.typing.fail.clear()
    _run_genome(genome)
    assert genome.typing.runs["rbceq2"] == 2 and genome.haplocheck.runs == 2


def test_haplocheck_image_missing_is_a_failure(
    genome: Genome, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ruling R21: on a lab bootstrapped before slice 3f, haplocheck's image is missing;
    the run finishes with the chrM heteroplasmies unknown, naming it. Docker not running
    still refuses the run."""

    def no_image(*_: object, **__: object) -> None:
        raise DockerUnavailable(f"Docker could not start {IMAGES['haplocheck']}")

    monkeypatch.setattr(haplocheck, "run_container", no_image)
    monkeypatch.setattr(pipeline, "image_missing", lambda name: name == "haplocheck")
    _chrm_call(genome, 0.40)
    run_id = _run_genome(genome)
    with connect(genome.store) as con:
        chrm = _row(con, run_id, CHRM_CLAIM[0])
    assert (chrm["kind"], chrm["dosage"]) == ("unassessed", "unknown")
    assert any(
        "haplocheck failed: Docker could not start" in w for w in chrm["warnings"]
    )
    monkeypatch.setattr(pipeline, "image_missing", lambda name: False)
    _chrm_call(genome, 0.30)
    with pytest.raises(RunRefused, match="Docker could not start"):
        _run_genome(genome)


def test_run_without_gencc_uses_clingen_only(genome: Genome) -> None:
    """Without a GenCC download the run goes on with ClinGen's validity alone and records
    why; a refresh right after it needs no run, and a GenCC file appearing is a new
    grading input."""
    lab = genome.lab.paths.lab
    newest = lab / "clinvar" / "20261011" / "clinvar_GRCh38.vcf.gz"
    newest.parent.mkdir(parents=True)
    shutil.copy(genome.clinvar38, newest)
    shutil.copy(f"{genome.clinvar38}.tbi", f"{newest}.tbi")
    run_id = _run_genome(genome, clinvar38=newest)
    rows = _statements(genome, run_id)
    with connect(genome.store) as con:
        inputs = _inputs_json(con, run_id)
    assert inputs["gencc"] == NO_GENCC
    reason = rows["clinvar:8003"]["reason"]
    assert "no ClinGen gene validity" in reason and "GenCC" not in reason
    latest = snapshots(genome.lab.paths.reports / "snapshots")[-1]
    assert not refresh.run_needed(lab, latest, analyser_identity())
    _write_gencc(genome, [("GENEY", "MONDO:0000002", "Limited")])
    assert refresh.run_needed(lab, latest, analyser_identity())


def test_chrm_floor_is_max_of_config_and_table(genome: Genome) -> None:
    """The run's chrM floor is the larger of the config's (0.10) and the table's: the
    table can raise it, never lower it; a table withholding heteroplasmies makes them
    unknown with its reason; without a table the config's floor stands."""

    def graded(level: float, **table: Any) -> dict[str, Any]:
        if table:
            write_slice3f(genome.lab, **table)
        else:
            shutil.rmtree(tables.folder(genome.lab.paths.lab, "slice3f"))
        _chrm_call(genome, level)
        run_id = _run_genome(genome)
        with connect(genome.store) as con:
            return _row(con, run_id, CHRM_CLAIM[0])

    lowered = graded(0.06, chrm_floor=0.05)
    assert lowered["kind"] == "unassessed"
    assert any("below the floor 0.1" in w for w in lowered["warnings"])
    raised = graded(0.15, chrm_floor=0.20)
    assert raised["kind"] == "unassessed"
    assert any("below the floor 0.2" in w for w in raised["warnings"])
    why = "not covered: haplocheck missed a 0.05 mixture in calibration"
    withheld = graded(0.40, chrm_floor=None, withheld={"chrM heteroplasmy": why})
    assert withheld["kind"] == "unassessed"
    assert f"genotype unknown: {why}" in withheld["warnings"]
    config = graded(0.15)
    assert (config["kind"], config["dosage"]) == ("finding", "1")


def test_hhv6_reads_table_background(genome: Genome) -> None:
    """HHV-6 is classified against the table's background ceiling: with 0.03, a ratio of
    0.02 is a non-finding and 0.2 (between the ceiling and the band) unassessed."""
    write_slice3f(genome.lab, hhv6_background=0.03)
    kinds = []
    for ratio in (0.02, 0.2):
        genome.typing.ratio = ratio
        rows = _statements(genome, _run_genome(genome))
        kinds.append((rows["hhv6:ici"]["kind"], rows["hhv6:ici"]["reason"]))
    assert kinds == [("non_finding", HHV6_BELOW), ("unassessed", HHV6_BETWEEN)]
    assert genome.typing.runs["hhv6"] == 1


def test_rule_out_stored_flagged(
    genome: Genome, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """HLA calls that exclude DQ2 and DQ8 (T1K's DQA1 and DQB1, and the six tag SNPs
    called 0/0) store `hla:coeliac` as a non-finding flagged `exclusion` with the claim's
    `confirm`; no other slice 3f statement is flagged. HLA is typed once: the typing steps
    reuse the outside callers' MHC reads, OptiType and T1K. A table that trusts phase and
    UGT1A1*28 gives PharmCAT phase sets (`GT:PS`) and the trusted allele."""
    write_slice3f(genome.lab, phase_trusted=True, ugt1a1_trusted=True)
    t1k = tmp_path / "t1k_dq.tsv"
    t1k.write_text(
        genome.pgx.t1k.read_text()
        + "HLA-DQA1\t2\tHLA-DQA1*01:01:01\t0.5\t60\tHLA-DQA1*03:01:01\t0.5\t60\n"
        + "HLA-DQB1\t2\tHLA-DQB1*05:01:01\t0.5\t60\tHLA-DQB1*06:02:01\t0.5\t60\n"
    )
    genome.pgx.t1k = t1k
    tags = tuple(
        ("chr1", pos, genome_bases(genome.lab.fasta, "chr1", pos), "A")
        for pos in MINI_TAGS
    )
    monkeypatch.setattr(
        pipeline, "TAG_ALLELES", tuple((c, p, r, _alt(r)) for c, p, r, _ in tags)
    )
    run_id = _run_genome(genome)
    rows = _statements(genome, run_id)
    with connect(genome.store) as con:
        inputs = _inputs_json(con, run_id)
    coeliac = rows["hla:coeliac"]
    assert (coeliac["kind"], coeliac["exclusion"]) == ("non_finding", True)
    assert json.loads(coeliac["claim"])["confirm"] == CONFIRM_DQ
    assert _flagged(
        {s: r for s, r in rows.items() if r["source"] in SLICE3F_TOOLS}
    ) == {"hla:coeliac"}
    assert genome.pgx.runs["optitype"] == genome.pgx.runs["t1k"] == 1
    assert {line.split("\t")[8] for line in genome.pgx.vcf} == {"GT:PS"}
    assert (inputs["pgx"]["phase"], inputs["pgx"]["ugt1a1_28_trusted"]) == (None, True)
    assert inputs["slice3f"]["reason"] is None and inputs["slice3f"]["table"]
    assert not any(
        inputs["steps"][s]["reused"] for s in ("hla_reads", "optitype", "t1k")
    )
