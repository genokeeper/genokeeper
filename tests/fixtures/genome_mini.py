"""A GRCh38 lab for `genokeeper run`, built in a test's tmp_path: the mini reference of
`grch38_mini` plus a random chrM (bwa-indexed) where the bootstrap puts the GRCh38 FASTA, the
GIAB GRCh38 BEDs and the derived region BEDs, placeholder T7 checkpoint files, a ClinGen
export, a CGD table, an empty gnomAD FAF table, and DNBSEQ-T7-named paired FASTQ under `raw/`.

The containers a run starts are replaced by fakes: `FakeDeepVariant` writes a VCF and a gVCF
of given records (in the shapes DeepVariant 1.10 writes, `gvcf_mini`), `FakeVerifyBamID`
writes VerifyBamID2's `.selfSM`, `FakeMutserve` mutserve's two tables, `FakeHaplogrep`
Haplogrep3's haplogroup, `FakeWhatsHap` a phased VCF, `FakeStage3` what ExpansionHunter,
REViewer, SMNCopyNumberCaller, Parascopy and Yleaf wrote on HG002 (`fixtures/hg002/`), with
`write_stage3` the stage-3 table, and `FakeDeletions` what DELLY, Dysgu and Truvari wrote
there, with `write_sv` the `sv` step's files and `write_deletions` the deletion table,
`FakeSlice3d` what slice 3d's callers wrote there, with `write_slice3d` its lab files and
three tables, and slice 3f's `FakeHaplocheck` haplocheck's table, `FakeTyping` the typing
tools (at their wrappers in `genokeeper.pipeline`), with `write_slice3f` the slice 3f table.
Nothing here is committed as a genomic file.
"""

import gzip
import itertools
import json
import random
import shutil
import threading
from collections.abc import Callable, Mapping, Sequence
from compression import zstd
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam
import pytest

from fixtures.clingen_mini import write_dosage, write_gene_validity, write_recurrent
from fixtures.clingen_mini import write_regions as write_clingen_region_list
from fixtures.faf_mini import write_faf_table
from fixtures.grch38_mini import (
    CONTIGS,
    bwa_index,
    write_fastq_pair,
    write_reference,
    write_regions,
)
from fixtures.gvcf_mini import HEADER
from fixtures.pgsc_mini import PINNED, write_results, write_task
from fixtures.pharmgkb_mini import Ann, write_annotations
from genokeeper import pipeline
from genokeeper.annotate.vep import (
    CACHE_INFO,
    FIELDS,
    VEP_INPUT,
    VEP_OUTPUT,
    vep_files,
)
from genokeeper.calibrate import BENCH_BED, TRUTH_VCF, hg004_truth, layout, tables
from genokeeper.calibrate.classes import ClassTable, write_class_table
from genokeeper.calibrate.classes import measure as class_measure
from genokeeper.calibrate.classes import summary as class_summary
from genokeeper.calibrate.deletions import (
    DeletionTable,
    Geometry,
    deletion_versions,
    measure,
    summary,
)
from genokeeper.calibrate.deletions import write_table as write_deletion_table
from genokeeper.calibrate.records import Record, library_dir, write_record
from genokeeper.calibrate.slice3d import CNV_BANDS as CNV_TABLE_BANDS
from genokeeper.calibrate.slice3d import (
    CNV_FEATURES,
    INSERTION_BANDS,
    INSERTION_FEATURES,
    ParalogTable,
    cnv_versions,
    insertion_versions,
    paralog_versions,
    write_paralog_table,
)
from genokeeper.calibrate.slice3f import LIBRARIES as SLICE3F_LIBRARIES
from genokeeper.calibrate.slice3f import Slice3fTable, slice3f_versions
from genokeeper.calibrate.slice3f import write_slice3f as write_table3f
from genokeeper.calibrate.stage3 import (
    FALSE_CALL,
    Stage3Table,
    hg002_truth,
    stage3_versions,
    write_table,
)
from genokeeper.calibrate.table import (
    ALLELE_CLASSES,
    GENOTYPE,
    GQ_BINS,
    REGIONS,
    ROW_TYPES,
    CalibrationTable,
    CellKey,
    Counts,
)
from genokeeper.call.deepvariant import CHECKPOINT_FILES, checkpoint_path
from genokeeper.checksum import sha256_file
from genokeeper.config import Grading, Paths, QcThresholds, grch38_fasta, load_grading
from genokeeper.containers import DockerUnavailable
from genokeeper.deletions.exons import load_spans, locus_span
from genokeeper.deletions.files import sv_files
from genokeeper.deletions.screen import CallClass
from genokeeper.evidence.clingen import DOSAGE_FILE, RECURRENT_FILE, REGIONS_FILE
from genokeeper.evidence.frequency import faf_folder
from genokeeper.evidence.pgx import CPIC_FILES
from genokeeper.evidence.pharmgkb import DOWNLOADS
from genokeeper.masking import masked_bed
from genokeeper.paralogs import LOCI
from genokeeper.pgx.hla import t1k_index
from genokeeper.pgx.pharmcat import positions_path
from genokeeper.pipeline import DYSGU_INS
from genokeeper.prs.pgsc import (
    NF_SCHEMA_VERSION,
    ancestry_pvar,
    java_home,
    nextflow_path,
    nxf_home,
    panel_path,
    panel_pvar,
    pgsc_calc_dir,
    scoring_file,
)
from genokeeper.prs.scores import load_scores
from genokeeper.regions import QC_AUTOSOMES, SEX_REGIONS, giab38_dir
from genokeeper.sample_identity import DELIVERED_HEADER
from genokeeper.tools import ToolFailed
from genokeeper.typing_files import typing_files

CHRM_LENGTH = 16569
GENE_VALIDITY_DATE = "20261018"
CGD_DATE = "20250221"
# The CGD download's header (research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz, read
# 2026-10-20): twelve columns, the first commented.
CGD_COLUMNS = (
    "#GENE",
    "HGNC ID",
    "ENTREZ GENE ID",
    "CONDITION",
    "INHERITANCE",
    "AGE GROUP",
    "ALLELIC CONDITIONS",
    "MANIFESTATION CATEGORIES",
    "INTERVENTION CATEGORIES",
    "COMMENTS",
    "INTERVENTION/RATIONALE",
    "REFERENCES",
)
SAMPLE = "sample"
# VerifyBamID2's output columns (its README, "Output").
SELFSM_HEADER = (
    "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\tFREE_RH"
    "\tFREE_RA\tCHIPMIX\tCHIPLK1\tCHIPLK0\tCHIP_RH\tCHIP_RA\tDPREF\tRDPHET\tRDPALT"
)


# A QC gate the mini sample (about 3x of chr1, 800 reads) passes: the shipped thresholds are
# for a 30x genome; each test that wants a refusal sets its own.
MINI_QC = QcThresholds(
    min_mean_depth=1.0,
    depth_mapq=20,
    min_fraction_20x=0.0,
    min_mapped_fraction=0.9,
    max_duplicate_fraction=0.5,
    min_median_insert=100.0,
    max_freemix=0.03,
    min_titv=0.0,
    max_titv=1e9,
    max_het_hom=1e9,
)
# The calibration table `write_calibration` writes: one genotype error in a million units
# in every cell.
CALIBRATED_ERRORS = (1, 10**6)


def lab_paths(root: Path) -> Paths:
    """`Paths` with the repository, raw data, reports and lab all under `root`, the lab
    holding an empty gnomAD FAF table (a test writes its rows with `write_faf_table`)."""
    paths = Paths(
        repo=root, raw=root / "raw", reports=root / "reports", lab=root / "lab"
    )
    write_faf_table(faf_folder(paths.lab))
    return paths


def write_cgd(path: Path, rows: Sequence[tuple[str, str]]) -> Path:
    """A CGD table of (gene, inheritance) rows, gzipped as the download is."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(CGD_COLUMNS)]
    for i, (gene, inheritance) in enumerate(rows, start=1):
        fields = [gene, str(i), str(1000 + i), "Some condition", inheritance]
        lines.append("\t".join(fields + ["N/A"] * (len(CGD_COLUMNS) - len(fields))))
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    return path


@dataclass(frozen=True)
class GenomeLab:
    """The inputs `run_genome` reads, all under one root."""

    paths: Paths
    fasta: Path
    fastq1: Path
    fastq2: Path
    gene_validity: Path
    cgd: Path


def write_lab(
    root: Path, *, cgd: Sequence[tuple[str, str]] = (), pairs: int = 400
) -> GenomeLab:
    """The lab and the raw FASTQ under `root`; `cgd` rows go into the CGD table."""
    paths = lab_paths(root)
    fasta = write_reference(grch38_fasta(paths.lab))
    # The mitochondrial genome too, random as the rest, at rCRS's length.
    rng = random.Random(16569)
    chrm = "".join(rng.choice("ACGT") for _ in range(CHRM_LENGTH))
    with fasta.open("a") as f:
        f.write(
            ">chrM\n"
            + "".join(f"{chrm[i : i + 60]}\n" for i in range(0, CHRM_LENGTH, 60))
        )
    Path(f"{fasta}.fai").unlink()
    pysam.faidx(str(fasta))
    bwa_index(fasta)
    masked_bed(paths.lab).write_text("")  # the mini reference masks nothing
    giab = giab38_dir(paths.lab)
    autosomes = ("chr1", 0, CONTIGS["chr1"], "autosomes")
    write_regions(giab / QC_AUTOSOMES, [autosomes])
    write_regions(
        giab / SEX_REGIONS,
        [autosomes, ("chrY", 0, CONTIGS["chrY"], "xdr"), ("chrX", 0, 12000, "xnonpar")],
    )
    for name in load_grading().strata.values():
        # One interval far from every claim of the tests: no claim carries a flag.
        (giab / name).write_bytes(gzip.compress(b"chr1\t29000\t29100\n"))
    model = checkpoint_path(paths.lab).parent
    model.mkdir(parents=True)
    for name in CHECKPOINT_FILES:
        (model / name).write_text("placeholder\n")
    validity = write_gene_validity(
        paths.lab / "clingen" / GENE_VALIDITY_DATE / "gene-validity.csv",
        [("BRCA1", "Definitive")],
    )
    table = write_cgd(paths.lab / "cgd" / CGD_DATE / "CGD.txt.gz", cgd)
    write_vep(paths.lab)
    write_parascopy(paths.lab)
    paths.raw.mkdir(parents=True)
    fastq1, fastq2 = write_fastq_pair(paths.raw, fasta, pairs=pairs)
    return GenomeLab(
        paths=paths,
        fasta=fasta,
        fastq1=fastq1,
        fastq2=fastq2,
        gene_validity=validity,
        cgd=table,
    )


def write_delivered(path: Path, rows: Sequence[tuple[str, str, int, str]]) -> Path:
    """A delivered genotype file in 23andMe's format at `path`: a comment, the column line
    (`DELIVERED_HEADER`), then (rsid, chromosome, position, genotype) rows."""
    lines = ["# a delivered genotype file, made by the test", DELIVERED_HEADER]
    lines += [
        f"{rsid}\t{chrom}\t{pos}\t{genotype}" for rsid, chrom, pos, genotype in rows
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


# Slice 1e: a variant no ClinVar record classifies, which the fake VEP annotates as a stop-gain
# on BRCA1's MANE transcript, and its gene's base-level pext.
PREDICTED_POS = 25000
# A variant ClinVar does classify (the mini ClinVar's BRCA1 het, `BRCA1_HET`): VEP is given
# it all the same and the fake annotates it as a stop-gain too, which grading must not
# predict.
CLASSIFIED_LOF_POS = 1000
BRCA1_ID = "ENSG00000012048"
DOSAGE_DATE = "20261023"
# UCSC's NMD Escape track rows for BRCA1's MANE transcript, as bigBedToBed writes them from
# nmdEscMane.bb (read 2026-10-23): rule 1 twice, rules 4 and 3.
_MOUSEOVER = "<b>Rule {}</b> - {} (ENST00000357654.9)"
BRCA1_NMD_ROWS = [
    f"chr17\t{start}\t{end}\tBRCA1\t0\t-\t{start}\t{end}\t{color}\t"
    f"{_MOUSEOVER.format(rule, text)}\tENST00000357654.9\tNM_007294.4"
    for start, end, color, rule, text in (
        (
            43045677,
            43045802,
            "255,0,0",
            1,
            "CDS within 50bp of the last splice junction",
        ),
        (
            43047642,
            43047692,
            "255,0,0",
            1,
            "CDS within 50bp of the last splice junction",
        ),
        (43091434, 43094860, "255,215,0", 4, "PTC in long exon (>400 nt)"),
        (43115759, 43115779, "139,0,0", 3, "first 100bp of coding nucleotides"),
    )
]


def write_vep(lab: Path) -> None:
    """The bootstrap's slice-1e files: the VEP cache's info.txt and placeholders for the
    score files VEP reads (it is faked), a pext table of BRCA1 at `CLASSIFIED_LOF_POS` and `PREDICTED_POS`, an empty
    NMD Escape track, and ClinGen's dosage list rating BRCA1 haploinsufficient."""
    files = vep_files(lab)
    info = files.folder / CACHE_INFO
    info.parent.mkdir(parents=True)
    info.write_text("species\thomo_sapiens\n")
    for path in (files.revel, files.alphamissense, files.splicevault):
        path.write_text("placeholder\n")
        path.with_name(f"{path.name}.tbi").write_text("placeholder\n")
    plain = files.pext.with_suffix("")
    plain.write_text(
        "".join(
            f"chr1\t{pos - 1}\t{pos}\t{BRCA1_ID}\t0.9\n"
            for pos in (CLASSIFIED_LOF_POS, PREDICTED_POS)
        )
    )
    pysam.tabix_index(str(plain), preset="bed", force=True)
    files.pext_max.write_text(f"{BRCA1_ID}\t0.9\n")
    files.nmd.write_text("".join(f"{row}\n" for row in BRCA1_NMD_ROWS))
    write_dosage(lab / "clingen-dosage" / DOSAGE_DATE / DOSAGE_FILE, [("BRCA1", "3")])


class FakeVep:
    """Stands in for `run_container` in `genokeeper.annotate.vep`: VEP writes a stop-gain on
    BRCA1's MANE transcript for the variants at `PREDICTED_POS` and `CLASSIFIED_LOF_POS` and
    an intergenic row for any other, SpliceAI a score of 0 for every site; `sites` keeps the
    vids VEP was given."""

    def __init__(self) -> None:
        self.sites: list[str] = []
        self.runs: list[tuple[str, str | None]] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        self.runs.append((name, label))
        work = _mount(mounts, "/work")
        if name == "spliceai":
            given, written = (args[args.index(flag) + 1] for flag in ("-I", "-O"))
            rows = []
            for line in (work / Path(given).name).read_text().splitlines():
                if not line.startswith("#"):
                    f = line.split("\t")
                    line = "\t".join([*f[:7], f"SpliceAI={f[4]}|BRCA1|0|0|0|0|0|0|0|0"])
                rows.append(line)
            (work / Path(written).name).write_text("\n".join(rows) + "\n")
            return
        assert name == "vep"
        rows = []
        for line in (work / VEP_INPUT).read_text().splitlines():
            if line.startswith("#"):
                continue
            vid = line.split("\t")[2]
            self.sites.append(vid)
            empty = ["-"] * (len(FIELDS) - 6)
            if any(f":{pos}:" in vid for pos in (PREDICTED_POS, CLASSIFIED_LOF_POS)):
                row = [vid, BRCA1_ID, "BRCA1", "ENST00000357654", "Transcript"]
                row += ["stop_gained", "NM_007294.4", "-", "10/23", "1234/5592"]
                row += ["ENST00000357654.9:c.1234C>T"]
                row += ["ENSP00000350283.3:p.Gln412Ter", "-", "-", "-"]
            else:
                row = [vid, "-", "-", "-", "Transcript", "intergenic_variant", *empty]
            rows.append("\t".join(row))
        text = "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n"
        text += "#" + "\t".join(FIELDS) + "\n" + "".join(f"{r}\n" for r in rows)
        (work / VEP_OUTPUT).write_bytes(gzip.compress(text.encode()))


def write_parascopy(lab: Path) -> None:
    """Placeholders for the bootstrap's Parascopy data (Parascopy is faked), with the md5s
    the bootstrap leaves beside them."""
    root = lab / "parascopy"
    for name in ("homology_table/GRCh38.bed.gz", "models_GRCh38_1KGP/EUR/SMN1.gz"):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_bytes(b"placeholder\n")
    (root / "homology_table" / "GRCh38.md5").write_text(f"{'0' * 32}  GRCh38.bed.gz\n")
    for tarball in ("GRCh38_v1.7.tar.gz", "models_GRCh38_1KGP_v1.7.tar.gz"):
        (root / f".unpacked.{tarball}").write_text(f"{'1' * 32}\n")


# A table withholding AR (40 repeats on the HG002 fixtures) alone: `measure` there also
# withholds ATXN2 and CACNA1A, kept here so the run grades their watch items.
STAGE3_WITHHELD = {"AR": FALSE_CALL}


def write_stage3(
    lab: GenomeLab,
    withheld: Mapping[str, str] = STAGE3_WITHHELD,
    *,
    grading: Grading | None = None,
    **versions: str,
) -> Stage3Table:
    """A stage-3 table measured under the lab's own versions (`stage3_versions`, under
    `grading`, config/grading.toml by default), any of `versions` replacing them,
    withholding `withheld`."""
    grading = grading if grading is not None else load_grading()
    table = Stage3Table(
        versions=stage3_versions(lab.paths.lab, sha256_file(lab.fasta), grading)
        | versions,
        withheld=dict(withheld),
        caps={},
        flank={"pcrfree": 1.0, "pcr": 0.5},
        measured={},
    )
    tables.write_table(
        lab.paths.lab, "stage3", lambda s: write_table(s, table), table.versions
    )
    return table


def write_sv(lab: GenomeLab) -> None:
    """The `sv` step's gnomAD-SV (one common deletion on chr20, away from every fixture
    call), bgzipped and indexed."""
    gnomad = sv_files(lab.paths.lab).gnomad
    gnomad.parent.mkdir(parents=True, exist_ok=True)
    plain = gnomad.with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Alleles">\n'
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr20\t1000\tcommon\tN\t<DEL>\t999\tPASS\tEND=3000;SVTYPE=DEL;AN=5000;AF=0.5\n"
    )
    pysam.tabix_index(str(plain), preset="vcf", force=True)


# The classes HG002 shows in `write_deletions`'s table: (class, length, confirmed) x count.
_DELETION_CALLS = [
    (CallClass("long", "both+depth", False, False), 2000, False, 4),
    (CallClass("long", "both+depth", False, False), 2000, True, 300),
    (CallClass("long", "dysgu", False, True), 3000, False, 6),
    (CallClass("long", "dysgu", False, True), 3000, True, 40),
    (CallClass("short", "dysgu", True, False), 150, False, 30),
    (CallClass("short", "dysgu", False, False), 150, True, 200),
    (CallClass("short", "both", False, False), 300, True, 400),
]


def write_deletions(lab: GenomeLab, **versions: str) -> DeletionTable:
    """A deletion table measured under the lab's own versions (`deletion_versions`), any of
    `versions` replacing them, on `_DELETION_CALLS` and a recall of 9 in 10 per band, each
    library's `summary` as the measured counts."""
    calls = [(c, n, ok) for c, n, ok, count in _DELETION_CALLS for _ in range(count)]
    truth = [(2000, True)] * 9 + [(2000, False)] + [(200, True)] * 9 + [(200, False)]
    grading = load_grading()
    bands, likelihoods = measure(
        {"pcrfree": calls, "pcr": calls},
        {"pcrfree": truth, "pcr": truth},
        Geometry(
            tuple(("chr1", n * 10**5, n * 10**5 + 297) for n in range(1, 4013)),
            2_500_000_000,
        ),
        {"long": 0.063, "short": 0.008},
        grading.deletions,
    )
    table = DeletionTable(
        versions=deletion_versions(lab.paths.lab, sha256_file(lab.fasta), grading)
        | versions,
        bands=bands,
        likelihoods=likelihoods,
        measured={
            lib: summary(calls, truth, grading.deletions) for lib in ("pcrfree", "pcr")
        },
    )
    tables.write_table(
        lab.paths.lab,
        "deletions",
        lambda s: write_deletion_table(s, table),
        table.versions,
    )
    return table


def write_calibration(lab: GenomeLab, scratch: Path) -> CalibrationTable:
    """HG003's benchmark BED over the whole mini chr1, and a calibration table measured on
    the lab's own pipeline (`pipeline.calibration_versions`, its digests cached under
    `scratch`): `CALIBRATED_ERRORS` in every cell."""
    bench = layout(lab.paths.lab)["truth"] / BENCH_BED
    bench.parent.mkdir(parents=True, exist_ok=True)
    bench.write_text(f"chr1\t0\t{CONTIGS['chr1']}\n")
    errors, total = CALIBRATED_ERRORS
    counts = {
        (GENOTYPE, CellKey(row_type, allele_class, region, gq)): Counts(
            errors, total, 30.0 * total
        )
        for row_type in ROW_TYPES
        for allele_class in ALLELE_CLASSES
        for region in REGIONS
        for gq in GQ_BINS
    }
    versions = pipeline.calibration_versions(
        lab.paths.lab,
        tools=pipeline.genome_tool_versions(),
        identities=scratch / "identities.json",
    )
    table = CalibrationTable.from_counts(counts, versions)
    tables.write_table(lab.paths.lab, "small_variants", table.write, versions)
    return table


def bases(fasta: Path, contig: str, pos: int, length: int = 1) -> str:
    """The reference bases at 1-based `pos`."""
    with pysam.FastaFile(str(fasta)) as f:
        return f.fetch(contig, pos - 1, pos - 1 + length).upper()


def variant(
    fasta: Path,
    contig: str,
    pos: int,
    alts: list[str],
    *,
    gt: str,
    gq: int,
    ad: list[int],
    filt: str = "PASS",
    ref_len: int = 1,
) -> str:
    """A DeepVariant variant record on `fasta`: ALTs end with `<*>`, AD lists REF and `alts`
    (the `<*>` 0 is added), VAF follows from AD."""
    depth = sum(ad)
    vaf = ",".join(f"{n / depth:.6g}" if depth else "0" for n in ad[1:])
    n = len(alts) + 2
    pl = ",".join(["0"] + ["990"] * (n * (n + 1) // 2 - 1))
    return (
        f"{contig}\t{pos}\t.\t{bases(fasta, contig, pos, ref_len)}"
        f"\t{','.join([*alts, '<*>'])}\t40.0\t{filt}\t.\tGT:GQ:DP:AD:VAF:PL"
        f"\t{gt}:{gq}:{depth}:{','.join(map(str, [*ad, 0]))}:{vaf},0:{pl}"
    )


def block(
    fasta: Path, contig: str, start: int, end: int, *, gq: int = 50, min_dp: int = 17
) -> str:
    """A reference block over `start..end` (1-based, inclusive)."""
    pl = "0,0,0" if min_dp == 0 else f"0,{gq + 1},{10 * gq + 9}"
    return (
        f"{contig}\t{start}\t.\t{bases(fasta, contig, start)}\t<*>\t0\t.\tEND={end}"
        f"\tGT:GQ:MIN_DP:PL\t0/0:{gq}:{min_dp}:{pl}"
    )


def _write_vcf(path: Path, records: list[str], contigs: dict[str, int]) -> None:
    order = list(contigs)
    rows = sorted(
        records, key=lambda r: (order.index(r.split("\t")[0]), int(r.split("\t")[1]))
    )
    header = HEADER + "".join(
        f"##contig=<ID={c},length={n}>\n" for c, n in contigs.items()
    )
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + SAMPLE + "\n"
    plain = path.with_suffix("")
    plain.write_text(header + columns + "".join(f"{r}\n" for r in rows))
    pysam.tabix_index(str(plain), preset="vcf", force=True)


def _mount(mounts: Sequence[tuple[Path, str, bool]], target: str) -> Path:
    return next(host for host, inside, _ in mounts if inside == target)


class FakeDeepVariant:
    """Stands in for `run_container` in `genokeeper.call.deepvariant`: writes `records` as the
    gVCF, its non-block records as the VCF, both indexed, where the arguments name them;
    `args` keeps each call's arguments."""

    def __init__(self, records: list[str], contigs: dict[str, int]) -> None:
        self.records = records
        self.contigs = contigs
        self.calls = 0
        self.args: list[list[str]] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("deepvariant", "dna:call")
        self.calls += 1
        self.args.append(list(args))
        out = _mount(mounts, "/out")
        named = dict(a.removeprefix("--").split("=", 1) for a in args if "=" in a)
        gvcf = out / named["output_gvcf"].removeprefix("/out/")
        vcf = out / named["output_vcf"].removeprefix("/out/")
        _write_vcf(gvcf, self.records, self.contigs)
        _write_vcf(vcf, [r for r in self.records if "<*>\t0\t." not in r], self.contigs)


class FakeVerifyBamID:
    """Stands in for `run_container` in `genokeeper.qc`: writes `vb.selfSM` with FREEMIX
    `value` into the directory mounted at /out."""

    def __init__(self, value: float = 0.001) -> None:
        self.value = value
        self.calls = 0

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("verifybamid2", "dna:qc")
        self.calls += 1
        row = ["mini", "NA", "NA", "9000", "90000", "30", str(self.value), *["NA"] * 12]
        text = f"{SELFSM_HEADER}\n" + "\t".join(row) + "\n"
        (_mount(mounts, "/out") / "vb.selfSM").write_text(text)


# mutserve 2.0.3's two tables, their columns as it writes them (tests/test_mito.py).
MUTSERVE_COLUMNS = (
    "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMajorBase\tMajorLevel\tMinorBase"
    "\tMinorLevel\tCoverage\tCoverageFWD\tCoverageREV\tType\tMeanBaseQuality"
)
RAW_COLUMNS = (
    "SAMPLE\tPOS\tREF\tTOP-FWD\tMINOR-FWD\tTOP-REV\tMINOR-REV\tCOV-FWD\tCOV-REV"
    "\tCOV-TOTAL\tTYPE\tLEVEL"
)


class FakeMutserve:
    """Stands in for `run_container` in `genokeeper.mito`: writes `calls` (pos, ref, alt,
    level) at 2000x and every chrM position covered 2000x, but `uncovered`, not at all."""

    def __init__(self, calls: Sequence[tuple[int, str, str, float]]) -> None:
        self.calls = calls
        self.uncovered: Sequence[int] = ()
        self.runs = 0

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("mutserve", "dna:mito")
        self.runs += 1
        work = _mount(mounts, "/work")
        rows = [
            f"chrM.bam\tPASS\t{pos}\t{ref}\t{alt}\t{level}\t{alt}\t\t-\t\t2000\t1000"
            f"\t1000\t{1 if level >= 0.9 else 2}\t40"
            for pos, ref, alt, level in self.calls
        ]
        (work / "mutserve.txt").write_text("\n".join([MUTSERVE_COLUMNS, *rows]) + "\n")
        # haplocheck's input; the fake reads nothing of it.
        vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        (work / "mutserve.vcf.gz").write_bytes(gzip.compress(vcf.encode()))
        raw = [
            f"chrM.bam\t{pos}\tN\tN\t-\t-\t-\t1000\t1000\t2000\t0\t0.0"
            for pos in range(1, CHRM_LENGTH + 1)
            if pos not in self.uncovered
        ]
        (work / "mutserve_raw.txt").write_text("\n".join([RAW_COLUMNS, *raw]) + "\n")


# Haplogrep3 3.3.2's output, verbatim, on a synthetic chrM VCF of J1c's mutations in the
# image docker/haplogrep3.Dockerfile builds (2026-10-21).
HAPLOGREP_OUTPUT = (
    '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n'
    '"SAMPLE"\t"J1c"\t"1"\t"1.0000"\t"1-16569"\n'
)


class FakeHaplogrep:
    """Stands in for `run_container` in `genokeeper.haplogroup`: checks the command and
    writes Haplogrep3's output into the mounted output directory; `vcf` keeps the chrM VCF
    it was handed."""

    def __init__(self) -> None:
        self.args: list[str] = []
        self.vcf = ""
        self.runs = 0

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label, platform) == ("haplogrep3", "dna:haplogroup", None)
        self.runs += 1
        self.args = list(args)
        given = {target: (host, read_only) for host, target, read_only in mounts}
        assert given["/in"][1] and not given["/out"][1]
        self.vcf = (given["/in"][0] / "chrM.vcf").read_text()
        (given["/out"][0] / "haplogroups.txt").write_text(HAPLOGREP_OUTPUT)


class FakeWhatsHap:
    """Stands in for `run_container` in `genokeeper.phasing`: puts every site it is handed
    in one phase set, on alternating haplotypes (`in_trans`) or on one (`in_cis`), as
    WhatsHap 2.8 writes them (tests/test_phasing.py)."""

    def __init__(self, phase: str = "in_trans") -> None:
        self.phase = phase
        self.sites: list[tuple[str, int]] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("whatshap", "dna:phase")
        work = _mount(mounts, "/out")
        lines = (work / "sites.vcf").read_text().splitlines()
        header = [line for line in lines if line.startswith("#")]
        header.insert(
            -1, '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
        )
        out: list[str] = []
        for i, line in enumerate(line for line in lines if not line.startswith("#")):
            fields = line.split("\t")
            self.sites.append((fields[0], int(fields[1])))
            gt = "1|0" if self.phase == "in_trans" and i % 2 else "0|1"
            out.append("\t".join([*fields[:8], "GT:PS", f"{gt}:1"]))
        (work / "phased.vcf").write_text("\n".join(header + out) + "\n")


HG002 = Path(__file__).parent / "hg002"
# The depth per base `FakeStage3` gives every flank base (samtools bedcov).
STAGE3_FLANK_DEPTH = 10


class FakeStage3:
    """Stands in for `run_container` in `genokeeper.repeats.expansionhunter`, `genokeeper.smn`
    and `genokeeper.yhaplogroup`, and for `run_tool` in the first (samtools): each tool writes
    what it wrote on HG002 (`fixtures/hg002/`) where its arguments name it, REViewer an image
    per locus, and samtools bedcov `STAGE3_FLANK_DEPTH` per base. `runs` counts each
    container's runs; a container in `fail` exits 1; ExpansionHunter writes `vcf` and
    `json` (slice 3e's run on HG002, `hg002.eh3e.*`)."""

    def __init__(self) -> None:
        self.runs: dict[str, int] = {}
        self.fail: set[str] = set()
        self.vcf = HG002 / "hg002.eh3e.vcf.txt"
        self.json = HG002 / "hg002.eh3e.json"

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert label == {"smncaller": "dna:smn", "parascopy": "dna:smn"}.get(
            name, "dna:yleaf" if name == "yleaf" else "dna:repeats"
        )
        self.runs[name] = self.runs.get(name, 0) + 1
        if name in self.fail:
            raise ToolFailed(1, ["docker", name], "", f"{name} crashed\n")
        out = _mount(mounts, "/out")
        if name == "expansionhunter":
            shutil.copyfile(self.vcf, out / "eh.vcf")
            shutil.copyfile(self.json, out / "eh.json")
            (out / "eh_realigned.bam").write_bytes(b"")
        elif name == "reviewer":
            for locus in args[list(args).index("--locus") + 1].split(","):
                shutil.copyfile(
                    HG002 / "hg002.reviewer.HTT.svg", out / f"rev.{locus}.svg"
                )
        elif name == "smncaller":
            shutil.copyfile(HG002 / "hg002.smncaller.json", out / "smn.json")
        elif name == "parascopy":
            if args[1] == "cn-using":
                (out / "cn").mkdir()
                paralog = (HG002 / "hg002.parascopy.paralog.bed.txt").read_bytes()
                (out / "cn" / "res.paralog.bed.gz").write_bytes(gzip.compress(paralog))
        elif name == "yleaf":
            sample = Path(args[list(args).index("-cram") + 1]).stem
            folder = out / "yleaf"
            (folder / sample).mkdir(parents=True)
            shutil.copyfile(
                HG002 / "hg002.yleaf.hg_prediction.txt", folder / "hg_prediction.hg"
            )
            shutil.copyfile(
                HG002 / "hg002.yleaf.out.txt", folder / sample / f"{sample}.out"
            )
        else:
            raise AssertionError(f"unexpected container {name}")

    def tool(
        self,
        args: list[str],
        *,
        env: dict[str, str] | None = None,
        cwd: Path | None = None,
    ) -> str:
        if args[:2] == ["samtools", "sort"]:
            Path(args[args.index("-o") + 1]).write_bytes(b"")
        if args[:2] == ["samtools", "bedcov"]:
            rows = Path(args[-2]).read_text().splitlines()
            return "".join(
                f"{row}\t{STAGE3_FLANK_DEPTH * (int(row.split()[2]) - int(row.split()[1]))}\n"
                for row in rows
            )
        return ""


class FakeDeletions:
    """Stands in for `run_container` in `genokeeper.deletions.callers` and
    `genokeeper.deletions.truvari`: DELLY's template copy writes an exclude list, its first
    chromosome shard what DELLY wrote on HG002's chr20 (`fixtures/hg002/`) and every other
    shard that file's header alone; Dysgu writes what it wrote there on every shard, each
    record moved onto the shard's chromosome; Truvari between the
    callers writes the `hg002.callers.*` folder, and against the truth confirms every call
    and finds every truth deletion. `runs` counts each tool's containers (DELLY's
    shards, not its copy), `shards` names each caller's shards; a tool in `fail` exits 1."""

    def __init__(self) -> None:
        self.runs: dict[str, int] = {}
        self.shards: dict[str, set[str]] = {}
        self.fail: set[str] = set()
        self._lock = threading.Lock()
        self._first: str | None = None

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        out = _mount(mounts, "/out")
        if name == "delly" and args[0] == "cp":
            (out / "exclude.template.tsv").write_text("chrM\n")
            return
        assert label == f"dna:{name}"
        with self._lock:
            self.runs[name] = self.runs.get(name, 0) + 1
            if name in ("delly", "dysgu"):
                self.shards.setdefault(name, set()).add(out.name)
            first = self._first = self._first or out.name
        if name in self.fail:
            raise ToolFailed(1, ["docker", name], "", f"{name} crashed\n")
        if name == "delly":
            text = (HG002 / "hg002.delly.vcf.txt").read_text()
            if out.name != first:
                text = "".join(x for x in text.splitlines(True) if x.startswith("#"))
            (out / "delly.bcf").write_text(text)
        elif name == "dysgu":
            # Its records moved onto the shard's chromosome, so `bcftools view -t` keeps them.
            lines = (HG002 / "hg002.dysgu.vcf.txt").read_text().splitlines(True)
            header = [x for x in lines if x.startswith("#")]
            records = [x.split("\t", 1)[1] for x in lines if not x.startswith("#")]
            (out / "dysgu.vcf").write_text(
                "".join(header) + "".join(f"{out.name}\t{x}" for x in records)
            )
        elif name == "truvari" and "--includebed" in args:
            # Against the truth: every call confirmed, every truth deletion found.
            folder = out / "truth"
            folder.mkdir()
            calls = (
                _mount(mounts, "/calls") / Path(args[list(args).index("-c") + 1]).name
            )
            truth = (
                _mount(mounts, "/truth") / Path(args[list(args).index("-b") + 1]).name
            )
            shutil.copyfile(calls, folder / "tp-comp.vcf.gz")
            shutil.copyfile(truth, folder / "tp-base.vcf.gz")
            with gzip.open(calls, "rt") as f:
                header = "".join(x for x in f if x.startswith("#"))
            for part in ("fn", "fp"):
                (folder / f"{part}.vcf.gz").write_bytes(gzip.compress(header.encode()))
        elif name == "truvari":
            folder = out / list(args)[list(args).index("-o") + 1].removeprefix("/out/")
            folder.mkdir()
            for part in ("tp-base", "tp-comp", "fn", "fp"):
                fixture = HG002 / f"hg002.callers.{part}.vcf.txt"
                (folder / f"{part}.vcf.gz").write_bytes(
                    gzip.compress(fixture.read_bytes())
                )
        else:
            raise AssertionError(f"unexpected container {name}")


# --- slice 3d: copy number, insertions, paralog loci --------------------------------------

REGIONS_DATE = "20261030"
# A ClinGen region a large event may overlap, on chr20 away from every fixture event.
LARGE_REGION = ("chr20", 63_000_000, 63_700_000)


def _cnv_calls() -> list[tuple[CallClass, int, bool]]:
    """Per measured depth band and evidence, 20 true calls and 2 false ones."""
    return [
        (CallClass(band, evidence, False, False, False), 3000, ok)
        for band in CNV_TABLE_BANDS
        for evidence in CNV_FEATURES["evidence"]
        for ok in [True] * 20 + [False] * 2
    ]


def _insertion_calls() -> list[tuple[CallClass, int, bool]]:
    """Per insertion band and evidence, 20 true calls and 2 false ones."""
    return [
        (CallClass(band, evidence, False, False), 300, ok)
        for band in INSERTION_BANDS
        for evidence in INSERTION_FEATURES["evidence"]
        for ok in [True] * 20 + [False] * 2
    ]


def _class_table(
    lab: GenomeLab,
    name: str,
    versions: Mapping[str, str],
    calls: list[tuple[CallClass, int, bool]],
    bands: Sequence[str],
    features: Mapping[str, tuple[str, ...]],
    extra: Mapping[str, object],
) -> ClassTable:
    """A class table `name` of `calls` on two libraries, a recall of 9 in 10 per band."""
    truth = {
        lib: [(band, ok) for band in bands for ok in [True] * 9 + [False]]
        for lib in ("pcrfree", "pcr")
    }
    geometry = Geometry(
        tuple(("chr1", n * 10**5, n * 10**5 + 297) for n in range(1, 4013)),
        2_500_000_000,
    )
    measured, likelihoods = class_measure(
        {"pcrfree": calls, "pcr": calls},
        truth,
        dict.fromkeys(bands, geometry),
        dict.fromkeys(bands, 0.01),
        features,
        band=lambda c: c.band,
    )
    table = ClassTable(
        versions=dict(versions),
        bands=measured,
        likelihoods=likelihoods,
        measured={
            lib: class_summary(
                calls, truth[lib], features, bands, band=lambda c: c.band, key=str
            )
            for lib in truth
        },
        features=dict(features),
        extra=dict(extra),
    )
    write_class_table(lab.paths.lab, name, table)
    return table


def write_paralog_models(lab: Path) -> None:
    """Placeholders for the Parascopy models of `paralogs.LOCI`, each a `main` line over its
    locus's span; STRC's also a duplication whose copy leaves CATSPER2's span whole."""
    spans = load_spans()
    folder = lab / "parascopy" / "models_GRCh38_1KGP" / "EUR"
    for locus in LOCI.values():
        chrom, start, end = locus_span(spans, locus.genes)
        lines = [f"main\t{chrom}:{start + 1}-{end}\t*\tname={locus.model}\n"]
        if locus.locus == "STRC":
            lines.append(
                f"duplication\t{chrom}:{start + 1}-{end}"
                f"\t{chrom}:{start + 100_001}-{end + 100_000}:+\tix=0\n"
            )
        (folder / f"{locus.model}.gz").write_bytes(
            gzip.compress("".join(lines).encode())
        )


def write_clingen_regions(lab: Path) -> None:
    """ClinGen's region lists, the bootstrap's `clingen_regions`: one HI-3 region
    (`LARGE_REGION`) and one recurrent row there."""
    folder = lab / "clingen-regions" / REGIONS_DATE
    chrom, start, end = LARGE_REGION
    write_clingen_region_list(
        folder / REGIONS_FILE,
        [("20q13.33 region", f"{chrom}:{start + 1}-{end}", "3", "0")],
    )
    write_recurrent(
        folder / RECURRENT_FILE, [(chrom, start, end, "20q13.33", "255,111,0")]
    )


def write_cnv(
    lab: GenomeLab,
    *,
    cnvpytor: str = "kept",
    large_withheld: Mapping[str, str] | None = None,
    **versions: str,
) -> ClassTable:
    """A `cnv` table measured under the lab's own versions (`cnv_versions`), any of
    `versions` replacing them, recording CNVpytor `cnvpytor` (kept or dropped) and the large
    bands `large_withheld`."""
    grading = load_grading()
    own = cnv_versions(lab.paths.lab, sha256_file(lab.fasta), grading) | versions
    extra = {
        "cnvpytor": {"decision": cnvpytor},
        "large": {"calls": {}, "withheld": dict(large_withheld or {})},
    }
    return _class_table(
        lab, "cnv", own, _cnv_calls(), CNV_TABLE_BANDS, CNV_FEATURES, extra
    )


def write_insertions(lab: GenomeLab, **versions: str) -> ClassTable:
    """An `insertions` table measured under the lab's own versions, any of `versions`
    replacing them."""
    own = insertion_versions(lab.paths.lab, sha256_file(lab.fasta), load_grading())
    return _class_table(
        lab,
        "insertions",
        own | versions,
        _insertion_calls(),
        INSERTION_BANDS,
        INSERTION_FEATURES,
        {},
    )


def write_paralogs(
    lab: GenomeLab, withheld: Mapping[str, str] | None = None, **versions: str
) -> ParalogTable:
    """A `paralogs` table measured under the lab's own versions, any of `versions` replacing
    them, withholding `withheld` (locus -> reason; none by default)."""
    own = paralog_versions(lab.paths.lab, sha256_file(lab.fasta), load_grading())
    table = ParalogTable(
        versions=own | versions, withheld=dict(withheld or {}), measured={}
    )
    write_paralog_table(lab.paths.lab, table)
    return table


def write_slice3d(lab: GenomeLab) -> None:
    """Slice 3d's lab files (the paralog models, ClinGen's region lists) and its three
    current tables."""
    write_paralog_models(lab.paths.lab)
    write_clingen_regions(lab.paths.lab)
    write_cnv(lab)
    write_insertions(lab)
    write_paralogs(lab)


# The labels slice 3d's containers run behind.
SLICE3D_LABELS = frozenset(
    {
        "dna:cnv",
        "dna:cnvpytor",
        "dna:insurveyor",
        "dna:scramble",
        "dna:paralogs",
        "dna:gauchian",
        "dna:kilda",
    }
)
type RunContainer = Callable[..., None]


def _records(text: str, extra: Sequence[str]) -> str:
    """A VCF's header and its records with `extra` among them, sorted by the header's
    contig order and position (as `bcftools index` needs)."""
    lines = text.splitlines(True)
    header = [x for x in lines if x.startswith("#")]
    order = {
        x.split("ID=")[1].split(",")[0]: i
        for i, x in enumerate(header)
        if x.startswith("##contig=")
    }
    records = [x for x in lines if not x.startswith("#")]
    records += [f"{x.rstrip(chr(10))}\n" for x in extra]
    records.sort(key=lambda x: (order[x.split("\t")[0]], int(x.split("\t")[1])))
    return "".join(header) + "".join(records)


class FakeSlice3d:
    """Stands in for slice 3d's `run_container` calls (`cnv.callers`, `insertions`,
    `paralogs`, `truvari`'s union,
    `smn`'s Parascopy and `pgx.cyp2d6`'s Aldy: `over` routes a call here by its label, or
    Parascopy's depth and Truvari between Dysgu's and INSurVeyor's insertions), and for
    `run_tool` in `insertions` and `paralogs` (`tool`): each writes what it wrote on HG002
    (`fixtures/hg002/`, Task 6's runs), with `dellycnv`, `cnvpytor`, `meis` and `gauchian`
    adding records or replacing a result. `runs` counts each container (`delly cnv`,
    `parascopy depth` and `parascopy cn` named apart), `haploid` holds `delly cnv`'s
    `--haploid-chr` value (None without one), `excluded` the bare contig names on its
    exclude list; a name in `fail` exits 1."""

    def __init__(self) -> None:
        self.runs: dict[str, int] = {}
        self.fail: set[str] = set()
        self.haploid: set[str | None] = set()
        self.excluded: set[str] = set()
        self.dellycnv: list[str] = []
        self.cnvpytor: list[str] = []
        self.meis: list[str] = []
        self.gauchian: tuple[str, str] | None = None
        self._lock = threading.Lock()

    @staticmethod
    def handles(name: str, args: Sequence[str], label: str | None) -> bool:
        """Whether a container call is slice 3d's (or the depth step 3a and 3d share)."""
        return (
            label in SLICE3D_LABELS
            or (name == "parascopy" and args[1] == "depth")
            or (name == "truvari" and f"/base/{DYSGU_INS}" in args)
        )

    def over(self, other: RunContainer) -> RunContainer:
        """A `run_container` sending slice 3d's calls here and every other to `other`."""

        def run(
            name: str,
            args: Sequence[str],
            *,
            mounts: Sequence[tuple[Path, str, bool]],
            platform: str | None = None,
            label: str | None = None,
        ) -> None:
            target = self if self.handles(name, args, label) else other
            target(name, args, mounts=mounts, platform=platform, label=label)

        return run

    def _key(self, name: str, args: Sequence[str]) -> str:
        if name == "delly":
            return "delly cnv"
        if name == "parascopy":
            return f"parascopy {'depth' if args[1] == 'depth' else 'cn'}"
        return name

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        key = self._key(name, args)
        out = _mount(mounts, "/aldy-out" if name == "aldy" else "/out")
        with self._lock:
            self.runs[key] = self.runs.get(key, 0) + 1
            if key == "delly cnv":
                given = list(args)
                at = given.index("--haploid-chr") if "--haploid-chr" in given else None
                self.haploid.add(None if at is None else given[at + 1])
                lines = (out / "exclude.tsv").read_text().splitlines()
                self.excluded = {line for line in lines if "\t" not in line}
        if key in self.fail:
            raise ToolFailed(1, ["docker", name], "", f"{name} crashed\n")
        if key == "delly cnv":
            text = (HG002 / "hg002.dellycnv.vcf.txt").read_text()
            (out / "delly.bcf").write_text(_records(text, self.dellycnv))
        elif name == "cnvpytor":
            rows = (HG002 / "hg002.cnvpytor.calls.txt").read_text()
            (out / "calls.tsv").write_text(
                rows + "".join(f"{r}\n" for r in self.cnvpytor)
            )
        elif name == "insurveyor":
            text = (HG002 / "hg002.insurveyor.vcf.txt").read_text()
            (out / "work" / "out.pass.vcf.gz").write_text(text)
        elif name == "scramble":
            rows = (HG002 / "hg002.scramble.meis.txt").read_text()
            (out / "scramble_MEIs.txt").write_text(
                rows + "".join(f"{r}\n" for r in self.meis)
            )
        elif key == "parascopy depth":
            (out / "depth").mkdir()
            (out / "depth" / "depth.csv").write_text("# depth\n")
        elif key == "parascopy cn":
            (out / "cn").mkdir()
            loci = (HG002 / "hg002.parascopy.loci.bed.txt").read_bytes()
            (out / "cn" / "res.paralog.bed.gz").write_bytes(gzip.compress(loci))
        elif name == "gauchian":
            tsv, js = self.gauchian or (
                (HG002 / "hg002.gauchian.tsv.txt").read_text(),
                (HG002 / "hg002.gauchian.json.txt").read_text(),
            )
            (out / "gauchian.tsv").write_text(tsv)
            (out / "gauchian.json").write_text(js)
        elif name == "kilda":
            (out / "kilda-output").mkdir()
            shutil.copyfile(
                HG002 / "hg002.kilda.txt", out / "kilda-output" / "kilda_kiv2.CNs"
            )
        elif name == "aldy":
            shutil.copyfile(HG002 / "hg002.aldy.CYP2A6.txt", out / "CYP2A6.aldy")
        elif name == "truvari":
            folder = out / "both"
            folder.mkdir()
            for part in ("tp-base", "tp-comp", "fn", "fp"):
                fixture = HG002 / f"hg002.ins.{part}.vcf.txt"
                (folder / f"{part}.vcf.gz").write_bytes(
                    gzip.compress(fixture.read_bytes())
                )
        else:
            raise AssertionError(f"unexpected container {name}")

    def tool(
        self,
        args: list[str],
        *,
        env: dict[str, str] | None = None,
        cwd: Path | None = None,
    ) -> str:
        """samtools natively: `view` writes an empty BAM, `bedcov` gives
        `STAGE3_FLANK_DEPTH` per base."""
        if args[:2] == ["samtools", "view"]:
            Path(args[args.index("-o") + 1]).write_bytes(b"")
        if args[:2] == ["samtools", "bedcov"]:
            rows = Path(args[-2]).read_text().splitlines()
            return "".join(
                f"{row}\t{STAGE3_FLANK_DEPTH * (int(row.split()[2]) - int(row.split()[1]))}\n"
                for row in rows
            )
        return ""


# --- slice 1b: the pharmacogenomic inputs and callers ------------------------------------

CPIC_DATE = "20261020"
PHARMGKB_DATE = "20250705"
FIXTURES = Path(__file__).parent
# CPIC's data as the run reads it: the pairs' levels and MT-RNR1's allele names.
CPIC_PAIRS = [
    ("CYP2C19", "clopidogrel", "A"),
    ("CYP2D6", "codeine", "A"),
    ("HLA-B", "abacavir", "A"),
    ("DPYD", "fluorouracil", "A"),
]
MTRNR1_ALLELES = ["Reference", "m.1555A>G"]


def write_pgx(
    lab: GenomeLab,
    positions: Sequence[tuple[str, int, str, str, str]],
    annotations: Sequence[Ann],
) -> None:
    """The bootstrap's slice-1b files in the lab: PharmCAT's positions (chrom, pos, REF,
    ALT, gene), CPIC's pairs and MT-RNR1 alleles, the ClinPGx clinical
    annotations download of `annotations`, and a T1K index (a placeholder: T1K is faked)."""
    home = lab.paths.lab
    rows = [
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tPX={gene}\tGT\t0/0"
        for chrom, pos, ref, alt, gene in positions
    ]
    path = positions_path(home)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        gzip.compress(
            (
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n"
                + "".join(f"{r}\n" for r in rows)
            ).encode()
        )
    )
    cpic = home / "cpic" / CPIC_DATE
    cpic.mkdir(parents=True)
    pairs = [dict(zip(("genesymbol", "drugname", "cpiclevel"), p)) for p in CPIC_PAIRS]
    (cpic / CPIC_FILES["pairs"]).write_text(json.dumps(pairs))
    names = [{"name": name} for name in MTRNR1_ALLELES]
    (cpic / CPIC_FILES["MT-RNR1"]).write_text(json.dumps(names))
    write_annotations(
        home / "pharmgkb" / PHARMGKB_DATE / DOWNLOADS[0], list(annotations)
    )
    index = t1k_index(home)
    index.parent.mkdir(parents=True)
    index.write_text(">HLA00001 A*01:01:01:01\nACGT\n")


class FakePgx:
    """Stands in for `run_container` in `genokeeper.pgx`: each outside caller writes its
    result where its arguments name it (BCyrius `*1/*4`, Aldy its runs on public samples in
    `fixtures/`, OptiType and T1K theirs), and PharmCAT keeps its input VCF and outside
    calls (`vcf`, `outside`) and writes the published example report trimmed in
    `fixtures/pharmcat_report_mini.json`, with an HTML report beside it. `runs` counts each
    container's runs."""

    LABEL = "dna:pgx"

    def __init__(self) -> None:
        self.runs: dict[str, int] = {}
        self.vcf: list[str] = []
        self.outside: str | None = None
        # T1K's genotype table (slice 3f's tests type DQ from their own).
        self.t1k = FIXTURES / "t1k_genotype.tsv"
        # A container that exits with this status, and images Docker does not have.
        self.fail: dict[str, int] = {}
        self.missing: set[str] = set()

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert label == self.LABEL
        self.runs[name] = self.runs.get(name, 0) + 1
        if name in self.missing:
            raise DockerUnavailable(f"Docker could not start {name} (No such image)")
        if name in self.fail:
            raise ToolFailed(self.fail[name], ["docker", name], "", f"{name} crashed\n")
        if name == "bcyrius":
            (_mount(mounts, "/bcyrius-out") / "bcyrius.tsv").write_text(
                "Sample\tGenotype\tFilter\tCall method\tPreferred diplotype"
                "\tOther candidates\tActivity score\tPredicted phenotype\n"
                "bcyrius\t*1/*4\tPASS\tDirect\t\t\t1.0\tIntermediate Metabolizer\n"
            )
        elif name == "aldy":
            gene = args[list(args).index("-g") + 1]
            fixture = {
                "CYP2D6": "aldy_cyp2d6_na10860.aldy",
                "CYP2B6": "aldy_cyp2b6.aldy",
            }
            out = _mount(mounts, "/aldy-out") / f"{gene}.aldy"
            out.write_text((FIXTURES / fixture[gene]).read_text())
        elif name == "optitype":
            out = _mount(mounts, "/optitype-out") / "optitype_result.tsv"
            out.write_text((FIXTURES / "optitype_na11995_result.tsv").read_text())
        elif name == "t1k":
            out = _mount(mounts, "/t1k-out") / "t1k_genotype.tsv"
            out.write_text(self.t1k.read_text())
        elif name == "pharmcat":
            given = dict(itertools.pairwise(args))
            vcf = _mount(mounts, "/pharmcat-in") / Path(given["-vcf"]).name
            self.vcf = [
                line
                for line in vcf.read_text().splitlines()
                if not line.startswith("#")
            ]
            if "-po" in given:
                po = _mount(mounts, "/pharmcat-outside") / Path(given["-po"]).name
                self.outside = po.read_text(encoding="utf-8")
            out = _mount(mounts, "/pharmcat-out")
            report = (FIXTURES / "pharmcat_report_mini.json").read_text()
            (out / "pharmcat.report.json").write_text(report)
            (out / "pharmcat.report.html").write_text("<html>report</html>\n")
        else:
            raise AssertionError(f"unexpected container {name}")


def write_prs(lab: GenomeLab, scored: Sequence[int], ancestry: Sequence[int]) -> None:
    """The bootstrap's slice-1d files in the lab: every chosen score's scoring file, its
    variants at the `scored` mini chr1 positions (effect allele the ALT `_alt` gives), the
    panel's ancestry positions at `ancestry`, and placeholders for the panel, Nextflow, its
    plugin, pgsc_calc and Java (Nextflow is faked)."""
    home = lab.paths.lab

    def alleles(pos: int) -> tuple[str, str]:
        ref = bases(lab.fasta, "chr1", pos)
        return ref, "A" if ref != "A" else "G"

    columns = "rsID\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n"
    rows = "".join(
        f"rs{p}\t{alleles(p)[1]}\t{alleles(p)[0]}\t0.1\t1\t{p}\n" for p in scored
    )
    for score in load_scores():
        path = scoring_file(home, score.pgs_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            gzip.compress(f"#pgs_id={score.pgs_id}\n{columns}{rows}".encode())
        )
    ancestry_pvar(home).parent.mkdir(parents=True)
    for pvar in (ancestry_pvar(home), panel_pvar(home)):  # the panel's: the same here
        pvar.write_bytes(
            zstd.compress(
                (
                    "#CHROM\tPOS\tID\tREF\tALT\n"
                    + "".join(
                        f"1\t{p}\tv{p}\t{'\t'.join(alleles(p))}\n" for p in ancestry
                    )
                ).encode()
            )
        )
    panel_path(home).write_bytes(b"panel placeholder\n")
    for path in (
        nextflow_path(home),
        pgsc_calc_dir(home) / "main.nf",
        nxf_home(home) / "plugins" / f"nf-schema-{NF_SCHEMA_VERSION}" / "MANIFEST.MF",
        java_home(home) / "bin" / "java",
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder\n")


class FakeNextflow:
    """Stands in for `run_tool` in `genokeeper.prs.pgsc`: Nextflow running pgsc_calc keeps
    the target VCF's rows (`vcf`), leaves one task launched with the pinned image
    (`image`, None for an unpinned one) and writes the results at each score's match rate
    (`rates`, 0.9 by default); `fail`, when set, is what a failing run prints on stdout,
    beside its `.nextflow.log`. `runs` counts its runs."""

    def __init__(self) -> None:
        self.runs = 0
        self.fail: str | None = None
        self.vcf: list[str] = []
        self.rates: dict[str, float] = {}
        self.image: str | None = PINNED["plink2"]

    def __call__(
        self,
        args: list[str],
        *,
        env: dict[str, str] | None = None,
        cwd: Path | None = None,
    ) -> str:
        assert env is not None and env["NXF_OFFLINE"] == "true" and cwd is not None
        self.runs += 1
        if self.fail is not None:
            (cwd / ".nextflow.log").write_text("log\n")
            raise ToolFailed(1, args, self.fail, "")
        given = dict(itertools.pairwise(args))
        sheet = Path(given["--input"]).read_text().splitlines()[1]
        target = Path(f"{sheet.split(',')[1]}.vcf.gz")
        with gzip.open(target, "rt") as handle:
            self.vcf = [line for line in handle.read().splitlines() if line[0] != "#"]
        image = self.image or PINNED["plink2"].split("@")[0] + ":2.00a5.10"
        write_task(Path(given["-w"]), image)
        rates = {s.pgs_id: self.rates.get(s.pgs_id, 0.9) for s in load_scores()}
        write_results(Path(given["--outdir"]), rates)
        return ""


# --- slice 3f: haplocheck, the typing tools and the slice 3f table ---------------------------


class FakeHaplocheck:
    """Stands in for `run_container` in `genokeeper.haplocheck`: writes haplocheck's `--raw`
    table, HG002's (`fixtures/hg002/hg002.haplocheck.txt`: not contaminated) or, with
    `level` set, contaminated at that level; `fail` makes it exit 1. `runs` counts its
    runs."""

    def __init__(self) -> None:
        self.level: float | None = None
        self.fail = False
        self.runs = 0

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("haplocheck", "dna:mito")
        self.runs += 1
        if self.fail:
            raise ToolFailed(1, ["docker", name], "", "haplocheck crashed\n")
        header, row = (HG002 / "hg002.haplocheck.txt").read_text().splitlines()[:2]
        if self.level is not None:
            fields = row.split("\t")
            fields[1:3] = ['"YES"', f'"{self.level}"']
            row = "\t".join(fields)
        (_mount(mounts, "/out") / "haplocheck.raw.txt").write_text(f"{header}\n{row}\n")


# RBCeq2's one position in the mini lab: inside the gVCF's reference block, ABO's.
RBCEQ2_POS = 20080
HMMIX_DECODED = "decoded.diploid.txt"


class FakeTyping:
    """Stands in for slice 3f's typing tools at their wrappers in `genokeeper.pipeline`'s
    namespace (as `tests/test_calibrate_slice3f.py` fakes them), each writing what its
    reader takes: HHV-6 at `ratio` times the haploid human depth, one 6 Mb run of
    homozygosity, RBCeq2's HG002 tables on its one mini position (`RBCEQ2_POS`),
    fastmixture's EUR projection and hmmix's HG002 segments. `runs` counts each tool's
    runs (`hhv6`, `roh`, `rbceq2`, `fastmixture`, `hmmix`); a tool in `fail` exits 1."""

    def __init__(self, fasta: Path) -> None:
        self.fasta = fasta
        self.ratio = 0.0
        self.runs: dict[str, int] = {}
        self.fail: set[str] = set()

    def _ran(self, tool: str) -> None:
        self.runs[tool] = self.runs.get(tool, 0) + 1
        if tool in self.fail:
            raise ToolFailed(1, ["docker", tool], "", f"{tool} crashed\n")

    @staticmethod
    def _write(path: Path, text: str) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def install(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Replace the wrappers in `genokeeper.pipeline` through `monkeypatch`."""
        from genokeeper.bloodgroups import Rbceq2Positions
        from genokeeper.hhv6 import classify

        write = self._write
        abo = ("chr1", RBCEQ2_POS, bases(self.fasta, "chr1", RBCEQ2_POS), ".")
        ref = bases(self.fasta, "chr1", 20070)
        snp = ("chr1", 20070, ref, "A" if ref != "A" else "G")

        def unmapped(_cram: Path, _fasta: Path, out: Path) -> tuple[Path, Path]:
            self._ran("hhv6")
            return write(out / "u_1.fq", "@r\n"), write(out / "u_2.fq", "@r\n")

        def coverage(_bam: Path, out: Path) -> Path:
            return write(out, "#rname\tstartpos\tendpos\tmeandepth\n")

        def to_bed(vcf: Path, out: Path) -> Path:
            """plink2's `--make-bed`, which exits 7 on a row with two ALTs."""
            with gzip.open(vcf, "rt") as handle:
                if any("," in row.split("\t")[4] for row in handle if row[0] != "#"):
                    raise RuntimeError(
                        "target.bim cannot contain multiallelic variants"
                    )
            return write(out, "bed")

        def hhv6(_depths: Any, haploid: float, **kw: Any) -> Any:
            return classify({"NC_000898.1": self.ratio * haploid}, haploid, **kw)

        def roh(_vcf: Path, _af: Path, out: Path, **_kw: Any) -> Path:
            self._ran("roh")
            return write(out, "RG\tsample\tchr1\t100\t6000100\t6000001\t50\t80.0\n")

        def rbceq2(_vcf: Path, out: Path) -> Path:
            self._ran("rbceq2")
            out.mkdir(parents=True, exist_ok=True)
            for suffix in (
                "_geno.tsv",
                "_pheno_alphanumeric.tsv",
                "_pheno_numeric.tsv",
            ):
                shutil.copy(
                    HG002 / f"hg002.rbceq2{suffix}.txt", out / f"rbceq2{suffix}"
                )
            return out

        def projection(_bed: Path, _p: Path, _k: int, _out: Path) -> dict[str, float]:
            self._ran("fastmixture")
            return dict.fromkeys(("AFR", "CSA", "EAS", "MID"), 0.0) | {"EUR": 1.0}

        def hmmix(_vcf: Path, _files: Path, out: Path) -> Path:
            self._ran("hmmix")
            out.mkdir(parents=True, exist_ok=True)
            shutil.copy(HG002 / "hg002.hmmix.decoded.txt", out / HMMIX_DECODED)
            return out / HMMIX_DECODED

        for name, fake in {
            "extract_unmapped": unmapped,
            "align_viral": lambda _a, _b, _f, out, threads: write(out, "bam"),
            "viral_coverage": coverage,
            "classify_hhv6": hhv6,
            "sites_vcf": lambda _g, _af, _fa, out, **_kw: write(out, "sites"),
            "run_roh": roh,
            "rbceq2_positions": lambda _fasta: Rbceq2Positions(
                alleles=(abo,), systems={abo[:2]: frozenset({"ABO"})}, unchecked={}
            ),
            "run_rbceq2": rbceq2,
            "pruned_alleles": lambda _bim: [
                snp,
                (*snp[:3], "T" if ref != "T" else "C"),
            ],
            "to_bed": to_bed,
            "read_groups": lambda _path: ("AFR", "CSA", "EAS", "EUR", "MID"),
            "run_projection": projection,
            "run_hmmix": hmmix,
        }.items():
            monkeypatch.setattr(pipeline, name, fake)


def write_slice3f(
    lab: GenomeLab,
    *,
    withheld: Mapping[str, str] | None = None,
    ugt1a1_trusted: bool = False,
    phase_trusted: bool = False,
    chrm_floor: float | None = 0.05,
    hhv6_background: float | None = 0.03,
) -> Slice3fTable:
    """A slice 3f table measured under the lab's own versions (`slice3f_versions`): first
    a placeholder for every typing file and truth it hashes that the lab lacks, and the
    three called libraries' records (HG002 PCR-free, HG003, HG004) on the lab's FASTA."""
    home = lab.paths.lab
    for path in [
        *typing_files(home).required(),
        *hg002_truth(home),
        layout(home)["truth"] / TRUTH_VCF,
        layout(home)["truth"] / BENCH_BED,
        *hg004_truth(home),
    ]:
        if not path.is_file():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"{path.name}\n")
    for key in SLICE3F_LIBRARIES:
        folder = library_dir(home, key)
        cram, gvcf = folder / "lib.cram", folder / "call" / f"c-{key}" / "s.g.vcf.gz"
        gvcf.parent.mkdir(parents=True, exist_ok=True)
        cram.write_text("cram\n")
        gvcf.write_text("gvcf\n")
        record = Record(
            cram=str(cram),
            key=f"align-{key}",
            metrics={"mean_depth": 30.0},
            gvcf=str(gvcf),
            reference_sha256=sha256_file(lab.fasta),
            accepted=True,
            snp_f1=0.995 if key == "hg004" else None,
            replaced=None,
        )
        write_record(home, key, record)
    table = Slice3fTable(
        withheld=dict(withheld or {}),
        ugt1a1_trusted=ugt1a1_trusted,
        phase_trusted=phase_trusted,
        chrm_floor=chrm_floor,
        hhv6_background=hhv6_background,
    )
    versions = slice3f_versions(home, load_grading())
    tables.write_table(home, "slice3f", lambda s: write_table3f(s, table), versions)
    return table
