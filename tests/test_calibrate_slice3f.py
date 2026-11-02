"""Slice 3f's calibration table (`genokeeper.calibrate.slice3f`): each measure at its edge,
the versions it is measured under, and the table written whole and read back. Every input is
synthetic: the calibration itself runs on the trio's data, never in a test."""

import json
import shutil
from collections.abc import Callable, Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any

import pysam
import pytest

from genokeeper import cli, pipeline
from genokeeper.ancestry.archaic import ArchaicShare
from genokeeper.bloodgroups import Rbceq2Positions, SystemCall
from genokeeper.calibrate import (
    BENCH_BED,
    TRUTH_VCF,
    build,
    hg004_truth,
    layout,
    records,
    slice3f,
    tables,
)
from genokeeper.calibrate.records import Record, library_dir, write_record
from genokeeper.calibrate.slice3f import (
    ARCHAIC_KEY,
    CHRM_KEY,
    DQ_KEY,
    HHV6_KEY,
    LIBRARIES,
    ROH_KEY,
    SPIKE,
    Slice3fTable,
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
    open_slice3f,
    slice3f_versions,
    write_slice3f,
)
from genokeeper.calibrate.stage3 import hg002_truth
from genokeeper.calibrate.tables import Current, open_table, write_table
from genokeeper.checksum import sha256_file
from genokeeper.config import (
    Hhv6Thresholds,
    RohThresholds,
    grch38_fasta,
    load_grading,
    load_paths,
)
from genokeeper.containers import IMAGES, DockerUnavailable
from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.haplocheck import Haplocheck
from genokeeper.hhv6 import Hhv6Call
from genokeeper.mito import MitoCall
from genokeeper.pgx.hla import HlaCall, t1k_index
from genokeeper.pgx.hla_typing import TAG_RSIDS, load_hg002_truth
from genokeeper.pgx.pharmcat import UGT1A1_28, positions_path
from genokeeper.pipeline import RunRefused, TypingResults, _Step
from genokeeper.roh import Run
from genokeeper.store.db import Kind, Level, RowType
from genokeeper.tools import ToolFailed
from genokeeper.typing_files import typing_files

TRUTH = load_hg002_truth()
GOOD_TAGS = dict.fromkeys(TAG_RSIDS, True)


# --- HLA ----------------------------------------------------------------------------------


def _hg002_calls() -> tuple[dict[str, HlaCall], dict[str, HlaCall]]:
    """OptiType's and T1K's calls on HG002 matching its published types (maternal first);
    T1K's maternal DQA1 is the full-sequence name, one of the truth's two."""
    optitype = {gene: TRUTH[gene] for gene in ("HLA-A", "HLA-B", "HLA-C")}
    t1k = {**optitype, "HLA-DQA1": ("*01:05", "*03:01")}
    t1k |= {gene: TRUTH[gene] for gene in ("HLA-DQB1", "HLA-DRB1")}
    return optitype, t1k


def _parent(side: int, other: str) -> tuple[dict[str, HlaCall], dict[str, HlaCall]]:
    """A parent's calls: carrying HG002's allele from that side (0 maternal, 1 paternal),
    with `other` as its second allele at every gene."""
    optitype, t1k = _hg002_calls()

    def carry(calls: dict[str, HlaCall]) -> dict[str, HlaCall]:
        return {gene: (call[side], other) for gene, call in calls.items()}

    return carry(optitype), carry(t1k)


def _trio() -> dict[str, tuple[dict[str, HlaCall], dict[str, HlaCall]]]:
    return {
        "pcrfree": _hg002_calls(),
        "hg003": _parent(1, "*02:01"),
        "hg004": _parent(0, "*02:01"),
    }


def _tags(**libraries: Mapping[str, bool | None]) -> dict[str, dict[str, bool | None]]:
    return {lib: dict(libraries.get(lib, GOOD_TAGS)) for lib in LIBRARIES}


def test_hla_withheld_on_truth_mismatch() -> None:
    """HG002's calls matching the published types withhold nothing; OptiType's HLA-B off
    by one allele withholds HLA-B alone, and T1K's DQB1 off withholds DQB1 (class II read by
    T1K alone)."""
    assert measure_hla(_trio(), TRUTH, _tags()) == {}
    calls = _trio()
    optitype, t1k = calls["pcrfree"]
    calls["pcrfree"] = ({**optitype, "HLA-B": ("*35:01", "*38:01")}, t1k)
    withheld = measure_hla(calls, TRUTH, _tags())
    assert list(withheld) == ["HLA-B"]
    assert withheld["HLA-B"].startswith("not covered: OptiType")
    calls = _trio()
    optitype, t1k = calls["pcrfree"]
    calls["pcrfree"] = (optitype, {**t1k, "HLA-DQB1": ("*05:01", "*03:01")})
    assert list(measure_hla(calls, TRUTH, _tags())) == ["HLA-DQB1"]
    # Equivalent names that span the truth and another allele do not match.
    calls["pcrfree"] = (optitype, {**t1k, "HLA-DQB1": ("*05:01", "*03:02,*03:01")})
    assert list(measure_hla(calls, TRUTH, _tags())) == ["HLA-DQB1"]


def test_hla_trio_inconsistency_withholds() -> None:
    """HG004 not carrying HG002's maternal HLA-A (by either caller) withholds HLA-A, even
    with HG002's own call right; a parent's gene left uncalled withholds too."""
    calls = _trio()
    optitype, t1k = calls["hg004"]
    calls["hg004"] = ({**optitype, "HLA-A": ("*02:01", "*03:01")}, t1k)
    withheld = measure_hla(calls, TRUTH, _tags())
    assert list(withheld) == ["HLA-A"]
    assert "HG004" in withheld["HLA-A"]
    calls = _trio()
    optitype, t1k = calls["hg003"]
    calls["hg003"] = (optitype, {g: c for g, c in t1k.items() if g != "HLA-DQA1"})
    assert list(measure_hla(calls, TRUTH, _tags())) == ["HLA-DQA1"]


def test_hla_dq_tags_graded_only_inside_a_benchmark() -> None:
    """R2: a tag SNP called differently from a truth withholds DQ; a tag no library's
    benchmark covers with a matching call (rs2187668 and rs4639334 in v5.0q's gap) withholds
    DQ naming it, unless another library's benchmark grades it."""
    gap = {**GOOD_TAGS, "rs2187668": None, "rs4639334": None}
    graded_by_parents = _tags(pcrfree=gap)
    assert measure_hla(_trio(), TRUTH, graded_by_parents) == {}
    nowhere = _tags(pcrfree=gap, hg003=gap, hg004={**GOOD_TAGS, "rs2187668": None})
    withheld = measure_hla(_trio(), TRUTH, nowhere)
    assert list(withheld) == [DQ_KEY]
    assert "rs2187668" in withheld[DQ_KEY]
    assert "rs4639334" not in withheld[DQ_KEY]
    wrong = measure_hla(_trio(), TRUTH, _tags(hg003={**GOOD_TAGS, "rs7454108": False}))
    assert list(wrong) == [DQ_KEY]
    assert "rs7454108" in wrong[DQ_KEY]


# --- UGT1A1*28 and phase ------------------------------------------------------------------


def test_ugt1a1_outside_benchmark_not_trusted() -> None:
    """Review Focus 4: every library agreeing with a truth that carries TA7 is trusted, but
    not when one library's site lies outside its benchmark."""
    ours = {"pcrfree": "TA6/TA6", "hg003": "TA7/TA6", "hg004": "TA6/TA6"}
    truth = {"pcrfree": "TA6/TA6", "hg003": "TA6/TA7", "hg004": "TA6/TA6"}
    inside = dict.fromkeys(LIBRARIES, True)
    assert measure_ugt1a1(ours, truth, inside) is True
    assert measure_ugt1a1(ours, truth, {**inside, "hg004": False}) is False
    assert measure_ugt1a1({**ours, "hg003": None}, truth, inside) is False
    assert measure_ugt1a1({**ours, "hg003": "TA6/TA6"}, truth, inside) is False


def test_ugt1a1_reference_only_truth_not_trusted() -> None:
    """Every library and truth TA6/TA6, all agreeing: a measure that saw no TA7 proves
    nothing about it."""
    same = dict.fromkeys(LIBRARIES, "TA6/TA6")
    assert measure_ugt1a1(same, same, dict.fromkeys(LIBRARIES, True)) is False


def test_phase_needs_twenty_pairs() -> None:
    assert measure_phase(20, 0) is True
    assert measure_phase(19, 0) is False
    assert measure_phase(50, 1) is False


# --- chrM ---------------------------------------------------------------------------------

EXPECTED = frozenset({100, 200})


def _mixed(level: float, called: float | None = None) -> list[MitoCall]:
    """A mixture's calls at the expected positions, each at `called` (default `level`)."""
    at = level if called is None else called
    return [MitoCall(pos, "A", "G", at, 2000) for pos in sorted(EXPECTED)]


def _grid(missing_at: float = 0.02) -> dict[float, list[MitoCall]]:
    """Mixtures at levels near the grid; the one at `missing_at` misses position 200."""
    levels = (0.021, 0.049, 0.098, 0.203)
    found = {level: _mixed(level) for level in levels}
    low = next(level for level in levels if grid_step(level) == missing_at)
    found[low] = _mixed(low)[:1]
    return found


def test_chrm_floor_raised_by_false_heteroplasmy() -> None:
    """Every expected position called from the 0.05 mixture gives 0.05; a false
    heteroplasmy at 0.06 in an unmixed library raises it to the step above, 0.10. A call
    outside half to twice its mixture's level does not count."""
    assert measure_chrm_floor(_grid(), EXPECTED, []) == 0.05
    assert measure_chrm_floor(_grid(), EXPECTED, [0.03, 0.06]) == 0.10
    off = _grid()
    off[0.049] = _mixed(0.049, called=0.11)
    assert measure_chrm_floor(off, EXPECTED, []) == 0.10


def test_chrm_floor_none_above_grid() -> None:
    """A false heteroplasmy at 0.25 (or at 0.20) leaves no floor; so does no mixture
    passing."""
    assert measure_chrm_floor(_grid(), EXPECTED, [0.25]) is None
    assert measure_chrm_floor(_grid(), EXPECTED, [0.20]) is None
    none_pass = {level: calls[:1] for level, calls in _grid().items()}
    assert measure_chrm_floor(none_pass, EXPECTED, []) is None
    assert measure_chrm_floor(_grid(), frozenset(), []) is None  # nothing to measure


def test_mixture_level_from_counts() -> None:
    assert mixture_level(9_000, 1_000) == pytest.approx(0.10)


def test_grid_step_of_level() -> None:
    assert grid_step(0.093) == 0.10
    assert grid_step(0.07) == 0.05
    assert grid_step(0.012) is None


def _flag(contaminated: bool) -> Haplocheck:
    return Haplocheck(contaminated, 0.1 if contaminated else None, "H1", "H2")


def test_haplocheck_missed_mixture_withholds() -> None:
    """At floor 0.10 the 0.10 mixture unflagged withholds the heteroplasmies; only the 0.05
    mixture unflagged (below the floor) withholds nothing."""
    missed = {0.049: _flag(False), 0.098: _flag(False), 0.203: _flag(True)}
    reason = measure_haplocheck(missed, 0.10)
    assert reason == "not covered: haplocheck missed a 0.10 mixture in calibration"
    low_only = {0.049: _flag(False), 0.098: _flag(True), 0.203: _flag(True)}
    assert measure_haplocheck(low_only, 0.10) is None


# --- HHV-6, ROH and the archaic share ------------------------------------------------------


def test_hhv6_background_ceiling() -> None:
    """The highest background ratio times `background_factor`: 0.004–0.01 gives 0.03; one
    library at 0.2 gives 0.6, at or above `min_ratio`, which withholds HHV-6."""
    grading = Hhv6Thresholds()
    ceiling, why = measure_hhv6(
        {"pcrfree": 0.004, "hg003": 0.01, "hg004": 0.006}, grading
    )
    assert (ceiling, why) == (pytest.approx(0.03), None)
    ceiling, why = measure_hhv6(
        {"pcrfree": 0.004, "hg003": 0.2, "hg004": 0.006}, grading
    )
    assert ceiling == pytest.approx(0.6)
    assert why is not None and why.startswith("not covered:")


MB = 1_000_000
ROH = RohThresholds()
CHROM, START, END = SPIKE


def _run(start: int, length: int, chrom: str = "chr5") -> Run:
    return Run(chrom, start, start + length - 1, length, 80.0)


def _spiked() -> dict[str, list[Run]]:
    return {lib: [_run(START, END - START + 1, CHROM)] for lib in LIBRARIES}


def _roh(ours: Mapping[str, list[Run]], truth: Mapping[str, list[Run]]) -> str | None:
    return measure_roh(ours, truth, _spiked(), SPIKE, ROH)


def test_roh_false_long_run() -> None:
    """A run of ours of `min_run` or more that no truth run of 1 Mb or more overlaps by half
    its length withholds ROH; one a truth run covers by half keeps it."""
    truth = {lib: [_run(10 * MB, 4 * MB)] for lib in LIBRARIES}
    covered = {lib: [_run(9 * MB, 6 * MB)] for lib in LIBRARIES}
    assert _roh(covered, truth) is None
    false = {**covered, "hg004": [_run(40 * MB, 6 * MB)]}
    why = _roh(false, truth)
    assert why is not None and "hg004" in why.lower() and "false" in why


def test_roh_missed_truth_run() -> None:
    """A truth run of 6 Mb with no run of ours withholds ROH."""
    truth = {lib: [_run(10 * MB, 6 * MB)] for lib in LIBRARIES}
    ours = {lib: [_run(10 * MB, 6 * MB)] for lib in LIBRARIES}
    assert _roh(ours, truth) is None
    why = _roh({**ours, "hg003": []}, truth)
    assert why is not None and "missed" in why


def test_roh_spike_missed() -> None:
    """The spiked calls' longest run covering 80% of the 12 Mb spike withholds ROH; one
    covering it whole keeps it."""
    nothing: dict[str, list[Run]] = {lib: [] for lib in LIBRARIES}
    assert measure_roh(nothing, nothing, _spiked(), SPIKE, ROH) is None
    short = {**_spiked(), "pcrfree": [_run(START - 2_400_000, 12 * MB, CHROM)]}
    why = measure_roh(nothing, nothing, short, SPIKE, ROH)
    assert why is not None and "spike" in why
    assert ROH_KEY == "ROH"


def test_archaic_outside_published_band() -> None:
    """Published 0.015: 0.016 is kept, 0.030 and 0.007 are each withheld."""
    assert measure_archaic({"pcrfree": 0.016, "hg003": 0.015}, 0.015) is None
    for off in (0.030, 0.007):
        why = measure_archaic({"pcrfree": 0.016, "hg004": off}, 0.015)
        assert why is not None and "hg004" in why.lower()
    assert ARCHAIC_KEY == "neanderthal"


# --- blood groups ---------------------------------------------------------------------------


def test_bloodgroup_system_mismatch_withheld() -> None:
    """A system any library calls differently from RBCeq2 on its truth is withheld; one
    measured nowhere (None: a no-call position on either side) is withheld too, and so is
    one neither side typed anywhere; agreeing systems are not."""
    a, o = SystemCall("ABO*A1.01/ABO*O.01.01", "A"), SystemCall("ABO*O/ABO*O", "O")
    fy, untyped = SystemCall("FY*02/FY*02", "Fy(a-b+)"), SystemCall(None, None)
    ours = {lib: {"ABO": a, "FY": fy, "KEL": None, "LU": untyped} for lib in LIBRARIES}
    truth = {lib: {"ABO": a, "FY": fy, "KEL": None, "LU": untyped} for lib in LIBRARIES}
    truth["hg003"] = {**truth["hg003"], "ABO": o}
    withheld = measure_bloodgroups(ours, truth)
    assert sorted(withheld) == ["ABO", "KEL", "LU"]
    assert "HG003" in withheld["ABO"] and withheld["KEL"].startswith("not covered")
    assert withheld["LU"] == "not covered: LU not typed in calibration"


# --- the table and its versions -------------------------------------------------------------


def make_lab(lab: Path) -> Path:
    """A lab with every input `slice3f_versions` reads, each a tiny file, and the three
    called libraries' records (HG002 PCR-free, HG003, HG004) on its FASTA."""
    fasta = grch38_fasta(lab)
    files = [
        fasta,
        t1k_index(lab),
        positions_path(lab),
        *typing_files(lab).required(),
        *hg002_truth(lab),
        layout(lab)["truth"] / TRUTH_VCF,
        layout(lab)["truth"] / BENCH_BED,
        *hg004_truth(lab),
    ]
    for path in files:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{path.name}\n")
    for key in LIBRARIES:
        folder = library_dir(lab, key)
        cram, gvcf = (
            folder / "lib.cram",
            folder / "call" / f"c-{key}" / "sample.g.vcf.gz",
        )
        gvcf.parent.mkdir(parents=True, exist_ok=True)
        cram.write_text("cram\n")
        gvcf.write_text("gvcf\n")
        write_record(
            lab,
            key,
            Record(
                cram=str(cram),
                key=f"align-{key}",
                metrics={"mean_depth": 30.0},
                gvcf=str(gvcf),
                reference_sha256=sha256_file(fasta),
                accepted=True,
                snp_f1=0.995 if key == "hg004" else None,
                replaced=None,
            ),
        )
    return lab


def test_versions_cover_every_resource_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The versions name the images, the FASTA, each library's CRAM and calls keys, each
    truth file and the [hhv6], [roh] and [mito] floor values; a new file under
    `resources/hla/` or a changed `[hhv6] background_factor` changes them."""
    lab = make_lab(tmp_path / "lab")
    resources = tmp_path / "resources"
    for name in slice3f.RESOURCE_DIRS:
        (resources / name).mkdir(parents=True)
        (resources / name / "a.tsv").write_text(f"{name}\n")
    monkeypatch.setattr(slice3f, "RESOURCES", resources)
    grading = load_grading()
    first = slice3f_versions(lab, grading)
    assert first["pcrfree:cram"] == "align-pcrfree"
    assert first["hg004:calls"] == "c-hg004"
    assert "pcr:cram" not in first
    assert {"optitype", "t1k", "hmmix", "fasta", "[hhv6]", "[roh]"} <= set(first)
    assert "resources/hla/a.tsv" in first
    (resources / "hla" / "b.tsv").write_text("new\n")
    assert slice3f_versions(lab, grading) != first
    (resources / "hla" / "b.tsv").unlink()
    assert slice3f_versions(lab, grading) == first
    more = replace(grading, hhv6=Hhv6Thresholds(background_factor=4.0))
    assert slice3f_versions(lab, more) != first
    floor = replace(grading, mito=replace(grading.mito, heteroplasmy_floor=0.2))
    assert slice3f_versions(lab, floor) != first
    homoplasmy = replace(grading, mito=replace(grading.mito, min_homoplasmy=0.8))
    assert slice3f_versions(lab, homoplasmy) != first
    assert slice3f_versions(lab, replace(grading, p_low=0.05)) != first


def test_run_reads_no_calibration_file(tmp_path: Path) -> None:
    """Ruling R22: once the table is written, a library's CRAM and gVCF and a truth file
    removed leave the run's table current (record keys, recorded truth digests); the
    calibration's own check still finds it stale."""
    lab = make_lab(tmp_path / "lab")
    grading = load_grading()
    table = Slice3fTable({}, True, True, 0.05, 0.03)
    versions = slice3f_versions(lab, grading)
    write_table(lab, "slice3f", lambda s: write_slice3f(s, table), versions)
    for path in library_dir(lab, "hg003").rglob("*"):
        if path.suffix in (".cram", ".gz"):
            path.unlink()
    slice3f.truths(lab)["hg004"][0].unlink()
    assert open_slice3f(lab, grading) == (versions, table, None)
    assert not isinstance(open_table(lab, "slice3f"), Current)


def test_table_round_trip(tmp_path: Path) -> None:
    """Written through `tables.write_table` and read back unchanged; other versions make
    it stale."""
    lab = tmp_path / "lab"
    table = Slice3fTable(
        withheld={CHRM_KEY: "not covered: x", HHV6_KEY: "not covered: y"},
        ugt1a1_trusted=True,
        phase_trusted=False,
        chrm_floor=None,
        hhv6_background=None,
        report={"phase": {"agree": 3, "disagree": 0}},
    )
    assert Slice3fTable.from_json(table.to_json()) == table
    versions = {"fasta": "f", "t1k": "t"}
    write_table(lab, "slice3f", lambda s: write_slice3f(s, table), versions)
    found = open_table(lab, "slice3f", versions=versions)
    assert isinstance(found, Current)
    assert Slice3fTable.from_json((found.folder / "table.json").read_text()) == table
    stale = open_table(lab, "slice3f", versions={**versions, "fasta": "g"})
    assert not isinstance(stale, Current) and stale.reason.startswith("stale: fasta")


# --- Task 16: the typing steps, build_slice3f and the registration ---------------------------

HG002_FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
HMMIX_HEADER = (
    "chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tadmixpopvariants"
    "\tAltaiNeandertal\tChagyrskaya-Phalanx\tDenisova\tVindija33.19\n"
)
ABO_ALLELE = ("chr9", 10, "C", "CA")


def _ref_call(_allele: object = None) -> Called:
    return Called("0/0", "0", RowType.PASS, 50, 30, None, None, ())


def _rbceq2_tables(_vcf: Path, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    for suffix in ("_geno.tsv", "_pheno_alphanumeric.tsv", "_pheno_numeric.tsv"):
        shutil.copy(
            HG002_FIXTURES / f"hg002.rbceq2{suffix}.txt", out / f"rbceq2{suffix}"
        )
    return out


def _fake_typing_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every tool `typing_steps` runs, faked at its wrapper in the pipeline's namespace:
    each writes what its reader takes."""

    def write(path: Path, text: str) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def reads(_cram: Path, _fasta: Path, out: Path) -> tuple[Path, Path]:
        return write(out / "unmapped_1.fq", "@r\n"), write(
            out / "unmapped_2.fq", "@r\n"
        )

    def hla_reads(_cram: Path, _fasta: Path, out: Path) -> tuple[Path, Path]:
        return write(out / "hla_1.fq", "@r\n"), write(out / "hla_2.fq", "@r\n")

    def optitype(_fq1: Path, _fq2: Path, out: Path) -> dict[str, HlaCall]:
        write(
            out / "optitype_result.tsv",
            "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
            "0\tA*01:01\tA*26:01\tB*35:08\tB*38:01\tC*04:01\tC*12:03\t9\t9.0\n",
        )
        return {}

    def t1k(_fq1: Path, _fq2: Path, out: Path, *, index: Path) -> dict[str, HlaCall]:
        write(
            out / "t1k_genotype.tsv",
            "HLA-DQB1\t2\tHLA-DQB1*05:01:01\t9\t60\tHLA-DQB1*03:02:01\t9\t60\n",
        )
        return {}

    def coverage(_bam: Path, out: Path) -> Path:
        return write(
            out,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\n"
            "NC_001664.4\t8090\t151288\t0\t0\t0\t0.0\nNC_000898.1\t8794\t153321\t1\t1\t0\t0.15\n",
        )

    def hmmix(_vcf: Path, _files: Path, out: Path) -> Path:
        return write(
            out / "decoded.diploid.txt",
            HMMIX_HEADER
            + "chr1\t0\t1000000\t1000000\tArchaic\t0.9\t9\t9\t9\t1\t1\t1\n",
        )

    for name, fake in {
        "extract_hla_reads": hla_reads,
        "run_optitype": optitype,
        "run_t1k": t1k,
        "extract_unmapped": reads,
        "align_viral": lambda _a, _b, _f, out, threads: write(out, "bam"),
        "viral_coverage": coverage,
        "sites_vcf": lambda _g, _af, _fa, out, **_kw: write(out, "sites"),
        "run_roh": lambda _v, _af, out, **_kw: write(
            out, "RG\tsample\tchr1\t100\t6000100\t6000001\t50\t80.0\n"
        ),
        "rbceq2_positions": lambda _fasta: Rbceq2Positions(
            alleles=(ABO_ALLELE,),
            systems={ABO_ALLELE[:2]: frozenset({"ABO"})},
            unchecked={"FY": ("chr1:159205735:193bp>T",)},
        ),
        "calls_at": lambda _g, alleles, **_kw: {a: _ref_call() for a in alleles},
        "run_rbceq2": _rbceq2_tables,
        "pruned_alleles": lambda _bim: [("chr1", 100, "A", "G")],
        "to_bed": lambda _vcf, out: write(out, "bed"),
        "read_groups": lambda _path: ("AFR", "CSA", "EAS", "EUR", "MID"),
        "run_projection": lambda _b, _p, _k, _o: (
            dict.fromkeys(("AFR", "CSA", "EAS", "MID"), 0.0) | {"EUR": 1.0}
        ),
        "run_hmmix": hmmix,
    }.items():
        monkeypatch.setattr(pipeline, name, fake)


def _typing(lab: Path, work: Path) -> tuple[TypingResults, dict[str, _Step]]:
    gvcf = lab / "calls" / "sample.g.vcf.gz"
    gvcf.parent.mkdir(parents=True, exist_ok=True)
    for path in (gvcf, gvcf.with_name("sample.vcf.gz")):
        path.write_text("calls\n")
    return pipeline.typing_steps(
        work,
        cram=lab / "lib.cram",
        cram_key="align-key",
        fasta=grch38_fasta(lab),
        fasta_sha="fasta-sha",
        vcf=gvcf.with_name("sample.vcf.gz"),
        gvcf=gvcf,
        calls_key="calls-key",
        mean_depth=30.0,
        population="EUR",
        male=True,
        lab=lab,
        grading=load_grading(),
        threads=2,
        table=None,
        level=lambda _a, _c: Level.HIGH,
        tools={"samtools": "1.22", "bwa": "0.7.19"},
    )


def test_typing_step_keys(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Each tool is one cached step, its results parsed; changing RBCeq2's image changes
    the blood-group step's key alone, and a rebuilt panel reference (fastmixture's P file)
    the admixture step's alone. The HLA steps are keyed as `outside_steps` keys them."""
    lab = make_lab(tmp_path / "lab")
    _fake_typing_tools(monkeypatch)
    work = tmp_path / "work"
    results, steps = _typing(lab, work)
    assert results.failed == {}
    assert set(steps) == {
        "hla_reads",
        "optitype",
        "t1k",
        "hhv6",
        "roh_sites",
        "roh",
        "bloodgroups",
        "admixture",
        "archaic",
    }
    assert results.optitype is not None and results.optitype["HLA-B"] == (
        "*35:08",
        "*38:01",
    )
    assert results.t1k == {"HLA-DQB1": ("*05:01", "*03:02")}
    assert results.hhv6 is not None and results.hhv6.ratio == pytest.approx(0.01)
    assert results.roh is not None and results.roh[0].length == 6_000_001
    assert results.bloodgroups is not None and results.bloodgroups["ABO"].genotype
    assert results.bloodgroup_unchecked == {"FY": ("chr1:159205735:193bp>T",)}
    assert results.admixture is not None and results.admixture["EUR"] == 1.0
    assert results.archaic is not None and results.archaic.neanderthal_bp == 1_000_000
    keys = {name: step.key for name, step in steps.items()}
    reads_key = pipeline._key(  # pyright: ignore[reportPrivateUsage]
        "hla_reads",
        pipeline.HLA_READS_STEP_VERSION,
        "align-key",
        "1.22",
        "fasta-sha",
        json.dumps(pipeline.MHC_REGION),
    )
    assert keys["hla_reads"] == reads_key
    monkeypatch.setitem(IMAGES, "rbceq2", "genomelab/rbceq2:rebuilt")
    _, again = _typing(lab, work)
    changed = {name for name, step in again.items() if step.key != keys[name]}
    assert changed == {"bloodgroups"}
    monkeypatch.setitem(IMAGES, "rbceq2", IMAGES["rbceq2"])
    typing_files(lab).reference_p().write_text("a rebuilt panel reference\n")
    _, rebuilt = _typing(lab, work)
    changed = {name for name, step in rebuilt.items() if step.key != again[name].key}
    assert changed == {"admixture"}


def test_typing_step_failure_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A tool that fails is in `failed` with its line and its field None; the others run.
    Docker not running refuses."""
    lab = make_lab(tmp_path / "lab")
    _fake_typing_tools(monkeypatch)

    def broken(_vcf: Path, _out: Path) -> Path:
        raise ToolFailed(3, ["rbceq2"], "", "bad VCF\n")

    monkeypatch.setattr(pipeline, "run_rbceq2", broken)
    results, steps = _typing(lab, tmp_path / "work")
    assert results.failed == {"RBCeq2": "exit 3: bad VCF"}
    assert results.bloodgroups is None and "bloodgroups" not in steps
    assert results.admixture is not None

    def down(*_args: object, **_kwargs: object) -> None:
        raise DockerUnavailable("Docker is not running")

    monkeypatch.setattr(pipeline, "run_hmmix", down)
    with pytest.raises(RunRefused, match="Docker is not running"):
        _typing(lab, tmp_path / "work2")


def test_slice3f_registered(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`TABLES["slice3f"]` needs nothing; `genokeeper calibrate --only slice3f` builds it
    once (its builder faked), with no flag of 3f's own."""
    assert tables.TABLES["slice3f"].needs == ()
    seen: list[tuple[Path, int | None]] = []

    def fake(lab: Path, *, grading: object = None, threads: int | None = None) -> Path:
        seen.append((lab, threads))
        return lab

    lab = tmp_path / "lab"
    monkeypatch.setattr(build, "build_slice3f", fake)
    monkeypatch.setattr(build, "load_paths", lambda: replace(load_paths(), lab=lab))
    assert cli.main(["calibrate", "--only", "slice3f", "--threads", "2"]) == 0
    assert seen == [(lab, 2)]


def _good_typing(key: str) -> TypingResults:
    optitype, t1k = _trio()[key]
    return TypingResults(
        optitype=optitype,
        t1k=t1k,
        hhv6=Hhv6Call("HHV-6B", 0.01, Kind.NON_FINDING),
        roh=[],
        bloodgroups={"ABO": SystemCall("ABO*O.01.01/ABO*O.01.01", "O")},
        bloodgroup_uncovered={},
        bloodgroup_unchecked={},
        admixture={"EUR": 0.9, "MID": 0.1},
        archaic=ArchaicShare(60_000_000, 2_000_000, 0.0118),
        failed={},
    )


def _fake_build(
    lab: Path, monkeypatch: pytest.MonkeyPatch, typed: list[str]
) -> dict[str, list[str]]:
    """The tools `build_slice3f` runs, faked in its namespace, and its library reads
    spied on; returns the keys each read was asked for."""
    read: dict[str, list[str]] = {"aligned": [], "called": []}

    def spy[T](kind: str, real: Callable[[Path, str], T]) -> Callable[[Path, str], T]:
        def wrapped(lab: Path, key: str) -> T:
            read[kind].append(key)
            return real(lab, key)

        return wrapped

    monkeypatch.setattr(build, "read_aligned", spy("aligned", records.read_aligned))
    monkeypatch.setattr(build, "read_called", spy("called", records.read_called))

    def typing_steps(
        work: Path, **kwargs: Any
    ) -> tuple[TypingResults, dict[str, _Step]]:
        key = next(k for k in LIBRARIES if library_dir(lab, k) == work)
        typed.append(key)
        sites = work / "sites.vcf.gz"
        sites.write_text("sites\n")
        return _good_typing(key), {"roh_sites": _Step(f"sites-{key}", sites, False)}

    def write(path: Path, text: str = "x\n") -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    hg003_homoplasmy = MitoCall(16_000, "C", "T", 0.99, 2000)
    noise = MitoCall(300, "A", "G", 0.012, 2000)

    def chrm(source: Path, _fasta: Path, out: Path, **_kw: object) -> list[MitoCall]:
        out.mkdir(parents=True, exist_ok=True)
        if "hg003" in str(source):
            return [hg003_homoplasmy]
        if source.name.startswith("mix"):
            level = float(source.name.removeprefix("mix").removesuffix(".bam"))
            return [replace(hg003_homoplasmy, level=level)]
        return [noise]

    def samtools(args: list[str]) -> str:
        if args[1:3] == ["view", "-c"]:
            name = Path(args[3]).name
            if name.startswith("hg002"):
                return "9000"
            if name.startswith("hg003"):
                return "10000"
            target = float(name.removeprefix("mix").split(".hg003")[0])
            return str(round(9000 * target / (1 - target)))
        return ""

    fakes: dict[str, Any] = {
        "typing_steps": typing_steps,
        "genome_tool_versions": lambda: {"samtools": "1.22", "bwa": "0.7.19"},
        "calls_at": lambda _g, alleles, **_kw: {a: _ref_call() for a in alleles},
        "truth_calls": lambda _v, _b, alleles, **_kw: {a: _ref_call() for a in alleles},
        "truth_ta": lambda _vcf: "TA6/TA6",
        "truth_sites_vcf": lambda _t, _b, _af, out: write(out),
        "spike": lambda _v, _region, out: write(out),
        "run_roh": lambda _v, _af, out, **_kw: write(out, ""),
        "rbceq2_positions": lambda _fasta: Rbceq2Positions((ABO_ALLELE,), {}, {}),
        "run_rbceq2": _rbceq2_tables,
        "call_chrm": chrm,
        "run_haplocheck": lambda vcf, _out: Haplocheck(
            "mix" in str(vcf), None, "H1", None
        ),
        "run_tool": samtools,
        "het_snvs": lambda _vcf, _regions: set(),
        "pharmcat_positions": lambda _path: [],
        "gene_regions": lambda _path: [],
        "phase_sites": lambda *_args: {},
        "truth_haplotypes": lambda _v, _b, _sites: {},
    }
    for name, fake in fakes.items():
        monkeypatch.setattr(build, name, fake)
    return read


def test_build_measures_three_libraries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With the tools faked, `build_slice3f` reads HG002 PCR-free, HG003 and HG004 (never
    HG002's PCR library), runs the typing steps on each, and writes a table the run's own
    versions open."""
    lab = make_lab(tmp_path / "lab")
    for vcf, bench in slice3f.truths(lab).values():
        Path(f"{vcf}.tbi").write_text("index\n")
        bench.write_text("chr2\t0\t242000000\nchr6\t0\t171000000\n")
    typed: list[str] = []
    read = _fake_build(lab, monkeypatch, typed)
    folder = build.build_slice3f(lab, grading=load_grading(), threads=2)
    assert typed == list(LIBRARIES)
    assert set(read["aligned"]) == set(read["called"]) == set(LIBRARIES)
    _, table, why = open_slice3f(lab)
    assert why is None and table is not None
    assert folder == tables.folder(lab, "slice3f")
    assert (folder / "report.txt").is_file()
    # HG002's HLA and every tag match; TA6 alone is never trusted; no phase pair.
    assert not any(k.startswith("HLA-") or k == DQ_KEY for k in table.withheld)
    assert (table.ugt1a1_trusted, table.phase_trusted) == (False, False)
    assert table.hhv6_background == pytest.approx(0.03)
    assert table.chrm_floor == 0.02
    assert table.report["population"].startswith("EUR")


def test_build_slice3f_refuses_missing_library(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A library 3c has not called is refused with 3c's own line naming the command that
    supplies it, before any tool runs; no table is written."""
    lab = make_lab(tmp_path / "lab")
    typed: list[str] = []
    _fake_build(lab, monkeypatch, typed)
    gvcf = records.read_record(lab, "hg004")
    assert gvcf is not None and gvcf.gvcf is not None
    Path(gvcf.gvcf).unlink()
    with pytest.raises(
        RunRefused,
        match=(
            r"^HG004 is not called; move \S+/hg004/pcrfree/aligned\.json aside, then run "
            r"`python -m genokeeper\.calibrate\.hg002 align hg004`$"
        ),
    ):
        build.build_slice3f(lab, grading=load_grading(), threads=2)
    assert typed == []
    assert not tables.folder(lab, "slice3f").exists()


def test_truth_readers(tmp_path: Path) -> None:
    """The calibration's truth side: an allele's copies inside the benchmark (none where
    the truth calls nothing, unknown outside), the UGT1A1 TA repeat from a truth's own
    spelling (C>CAT is TA7) and from our calls at PharmCAT's three ALTs, the phased
    truth's haplotypes, and cis/trans compared between neighbours of one phase set."""
    text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr2>\n##contig=<ID=chr6>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003\n"
        "chr2\t233760233\t.\tC\tCAT\t.\tPASS\t.\tGT\t0|1\n"
        "chr6\t32638107\t.\tC\tT\t.\tPASS\t.\tGT\t1|0\n"
        "chr6\t32713706\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\n"
    )
    plain = tmp_path / "truth.vcf"
    plain.write_text(text)
    vcf = Path(pysam.tabix_index(str(plain), preset="vcf", force=True))
    bed = tmp_path / "bench.bed"
    bed.write_text("chr2\t0\t242000000\nchr6\t32600000\t32700000\n")
    dq25, dq8, absent = (
        ("chr6", 32638107, "C", "T"),
        ("chr6", 32713706, "T", "C"),
        ("chr6", 32650000, "A", "G"),
    )
    found = build.truth_calls(
        vcf, bed, [dq25, dq8, absent], haploid=lambda _c, _p: False
    )
    assert [found[a].dosage for a in (dq25, dq8, absent)] == ["1", "unknown", "0"]
    assert build.truth_ta(vcf) == "TA6/TA7"

    def at(dosages: tuple[str, str, str]) -> dict[Allele, Called]:
        return {
            a: replace(_ref_call(), dosage=d)
            for a, d in zip(UGT1A1_28, dosages, strict=True)
        }

    assert build.ours_ta(at(("0", "1", "0"))) == "TA6/TA7"
    assert build.ours_ta(at(("2", "1", "0"))) is None
    assert build.ours_ta(at(("unknown", "0", "0"))) is None
    assert build.truth_haplotypes(vcf, bed, [dq25, dq8]) == {dq25: 0}
    a, b, c = (("chr6", p, "A", "G") for p in (100, 200, 300))
    ours = {a: (1, 0), b: (1, 1), c: (2, 0)}
    assert build.compare_phase(ours, {a: 0, b: 1, c: 0}) == (1, 0)
    assert build.compare_phase(ours, {a: 0, b: 0, c: 0}) == (0, 1)
