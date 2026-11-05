"""Slice 3c's foundation table: the `[dosage]` and `[af_guard]` bounds measured on the trio.

The four controls (`CONTROLS`: HG002's two libraries, HG003 and HG004, each aligned on the
masked reference) give their X and Y depth ratios (`sex.measure_ratios`) and, where called,
their chrX heterozygous share (`sex.x_het_share`, from the VCF beside the gVCF; HG002's PCR
library is aligned only). HG002 PCR-free's true heterozygous SNVs (`measure.true_het_vafs`,
over v5.0q's benchmark minus the mask and the correct copies) give the allele fractions a
germline heterozygote shows. From them `propose_dosage` and `propose_af_bound` propose the
bounds; the table holds the numbers (`report.json`) and the report (`report.txt`) with the
TOML lines Task 15 writes
into `config/grading.toml` in a commit citing it. Nothing here edits the config: the bounds
are config, set from the report (3a's pattern).

`tolerance` is never proposed: it is a design value (half the smallest mosaic share the
screen reports), which the controls only check: twice the largest deviation of a control's
ratio from its expected copy number must stay below it, or the proposal refuses (four cell
lines from one batch cannot set a floor for a saliva library from another lab). HG004 is read
only once accepted (`records.read_aligned` refuses an HG004 below `HG004_MIN_SNP_F1`, the
label swap's gate), before anything runs.

All arithmetic is in `decimal.Decimal` built from each value's `repr`: binary floats would put
2 x 0.03 a hair above 0.06 and round it up to 0.07.
"""

import json
import shutil
import tempfile
from bisect import bisect_left
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from datetime import date
from decimal import ROUND_CEILING, Decimal
from pathlib import Path

from genokeeper.calibrate import RTG_VERSION, layout, rtg_path
from genokeeper.calibrate.measure import true_het_vafs
from genokeeper.calibrate.records import (
    library_dir,
    read_aligned,
    read_called,
    read_record,
    reference_sha256,
)
from genokeeper.calibrate.stage3 import hg002_truth
from genokeeper.calibrate.tables import REPORT, folder, write_table
from genokeeper.checksum import sha256_file
from genokeeper.config import load_grading
from genokeeper.masking import (
    correct_copies_sha256,
    masked_bed,
    without_false_duplications,
)
from genokeeper.pipeline import GENOME_SAMPLE, RunRefused, genome_tool_versions
from genokeeper.regions import SEX_REGIONS, giab38_dir, read_named
from genokeeper.sex import Dosage, SexRefused, measure_ratios, x_het_share

CONTROLS = {"pcrfree": "male", "pcr": "male", "hg003": "male", "hg004": "female"}
# HG002's PCR (UNIVERSAL) library is aligned only (the plan's "Called libraries").
UNCALLED = frozenset({"pcr"})
# The design's "falls about 0.6% of the time": the share of true heterozygotes the
# allele-fraction bound may warn.
TRUE_HET_SHARE = 0.006
# Expected X and Y ratios over the autosomes per copy count.
_EXPECTED = {
    "male": (Decimal("0.5"), Decimal("0.5")),
    "female": (Decimal("1.0"), Decimal(0)),
}
_MIN_GAP = Decimal("0.1")  # HG004's share above the largest male's
_Y_FLOOR = Decimal("0.02")
_STEP = Decimal("0.01")


def _d(value: float) -> Decimal:
    return Decimal(repr(value))


def _up(value: Decimal) -> float:
    return float(value.quantize(_STEP, rounding=ROUND_CEILING))


def propose_dosage(
    measured: Mapping[str, Dosage], tolerance: float
) -> dict[str, float]:
    """`x_het_min`, halfway between the largest male share and HG004's, and
    `y_absent_max`, twice HG004's Y ratio and at least 0.02, each rounded up to 0.01;
    `ValueError` when twice a control's deviation from its expected ratio reaches
    `tolerance`, or HG004's share is not above the largest male's by 0.1."""
    limit = _d(tolerance)
    for name, d in measured.items():
        x_expected, y_expected = _EXPECTED[CONTROLS[name]]
        for what, value, expected in (
            ("X", d.x_ratio, x_expected),
            ("Y", d.y_ratio, y_expected),
        ):
            deviation = abs(_d(value) - expected)
            if 2 * deviation >= limit:
                raise ValueError(
                    f"{name}: {what} ratio {value} is {deviation} from {expected};"
                    f" twice that reaches the tolerance {tolerance}"
                )
    female = measured.get("hg004")
    males = [
        _d(d.x_het_share)
        for name, d in measured.items()
        if CONTROLS[name] == "male" and d.x_het_share is not None
    ]
    if female is None or female.x_het_share is None or not males:
        raise ValueError("hg004 and a male control need a chrX heterozygous share")
    top, her = max(males), _d(female.x_het_share)
    if not her > top + _MIN_GAP:
        raise ValueError(
            f"hg004's chrX heterozygous share {her} is not above the largest male's"
            f" {top} by {_MIN_GAP}"
        )
    return {
        "x_het_min": _up((top + her) / 2),
        "y_absent_max": _up(max(_Y_FLOOR, 2 * _d(female.y_ratio))),
    }


def propose_af_bound(vafs: Sequence[float]) -> float:
    """The largest multiple of 0.01 below which at most `TRUE_HET_SHARE` of `vafs` lie."""
    if not vafs:
        raise ValueError("no true heterozygous calls to set the allele-fraction bound")
    ordered = sorted(vafs)
    allowed = _d(TRUE_HET_SHARE) * len(ordered)
    best = 0
    for k in range(101):
        if bisect_left(ordered, k / 100) <= allowed:
            best = k
    return best / 100


def share_below(vafs: Sequence[float], bound: float) -> float:
    """The share of `vafs` below `bound`."""
    return sum(v < bound for v in vafs) / len(vafs) if vafs else 0.0


def foundation_versions(lab: Path) -> dict[str, str]:
    """The versions the foundation table is measured under: the FASTA's sha256, each
    control record's key, the sex regions BED's sha256, `[sex] mapq`, the `[dosage]`
    values the measurement reads, mosdepth's and rtg's versions, and the correct copies'
    sha256 (the true heterozygotes are read without them)."""
    grading = load_grading()
    d = grading.dosage
    keys = {}
    for control in CONTROLS:
        record = read_record(lab, control)
        keys[f"{control}_key"] = "" if record is None else record.key
    return {
        "fasta": reference_sha256(lab),
        **keys,
        "sex_regions": sha256_file(giab38_dir(lab) / SEX_REGIONS),
        "sex_mapq": str(grading.sex.mapq),
        "dosage": json.dumps(
            {
                "tolerance": d.tolerance,
                "het_low": d.het_low,
                "het_high": d.het_high,
                "min_x_sites": d.min_x_sites,
            },
            sort_keys=True,
        ),
        "mosdepth": genome_tool_versions()["mosdepth"],
        "rtg": RTG_VERSION,
        "correct_copies": correct_copies_sha256(),
    }


def _report(
    measured: Mapping[str, Dosage],
    f1: float | None,
    vafs: Sequence[float],
    dosage: Mapping[str, float],
    bound: float,
) -> str:
    today = date.today().isoformat()
    lines = ["control\tsex\tX ratio\tY ratio\tchrX het share\tchrX sites"]
    for name, d in measured.items():
        share = "not called" if d.x_het_share is None else f"{d.x_het_share:.4f}"
        lines.append(
            f"{name}\t{CONTROLS[name]}\t{d.x_ratio:.4f}\t{d.y_ratio:.4f}\t{share}"
            f"\t{d.x_sites}"
        )
    shown = "not scored" if f1 is None else f"{f1:.4f}"
    lines += [
        f"HG004 SNP F1 against its v4.2.1 truth: {shown}",
        (
            f"HG002 PCR-free true heterozygous SNVs: {len(vafs)};"
            f" below 0.30: {share_below(vafs, 0.30):.3%}"
        ),
        "",
        "proposed for config/grading.toml:",
        "[dosage]",
        f"x_het_min = {dosage['x_het_min']}",
        f"y_absent_max = {dosage['y_absent_max']}",
        f'calibrated_on = "{today} HG002 pcrfree/pcr, HG003, HG004"',
        "[af_guard]",
        f"min_het_vaf = {bound}",
        f'calibrated_on = "{today} HG002 pcrfree"',
    ]
    return "\n".join(lines) + "\n"


def build_foundation(lab: Path, threads: int) -> None:
    """Measure the four controls and HG002 PCR-free's true heterozygotes and write the
    foundation table (`report.json`, `report.txt`) with `foundation_versions`. Raises
    `RunRefused` on an unaccepted HG004 (first, before anything runs), a control not
    aligned or called on this reference, a missing input, or controls that cannot set the
    bounds; nothing is written then."""
    read_aligned(lab, "hg004")  # refuses an HG004 below its SNP F1 gate
    truth, bench = hg002_truth(lab)
    where = layout(lab)
    regions = giab38_dir(lab) / SEX_REGIONS
    for what, path, command in (
        ("GIAB's HG002 v5.0q truth", truth, "scripts/bootstrap.sh --with-hg002"),
        ("its benchmark BED", bench, "scripts/bootstrap.sh --with-hg002"),
        ("rtg's SDF", where["sdf"], "scripts/bootstrap.sh --with-masked-reference"),
        ("rtg-tools", rtg_path(lab), "scripts/bootstrap.sh --with-calibration"),
        (
            "the mask's BED",
            masked_bed(lab),
            "scripts/bootstrap.sh --with-masked-reference",
        ),
        ("the sex regions BED", regions, "scripts/bootstrap.sh"),
    ):
        if not path.exists():
            raise RunRefused(f"missing input: {what} at {path}; run {command}")
    grading = load_grading()
    xnonpar = read_named(regions, "xnonpar")
    measured: dict[str, Dosage] = {}
    for control in CONTROLS:
        cram, _, _ = read_aligned(lab, control)
        try:
            x, y = measure_ratios(cram, grading=grading, lab=lab, threads=threads)
        except (FileNotFoundError, SexRefused) as e:
            raise RunRefused(f"{control}: {' '.join(str(e).split())}") from e
        share, sites = None, 0
        if control not in UNCALLED:
            vcf = read_called(lab, control).with_name(f"{GENOME_SAMPLE}.vcf.gz")
            share, sites = x_het_share(vcf, xnonpar=xnonpar, grading=grading)
        measured[control] = Dosage(x, y, share, sites)
    gvcf = read_called(lab, "pcrfree")
    work = Path(
        tempfile.mkdtemp(prefix=".foundation-", dir=library_dir(lab, "pcrfree"))
    )
    try:
        vafs = true_het_vafs(
            gvcf,
            truth,
            without_false_duplications(lab, bench, work / "v5.0q.masked.bed"),
            sdf=where["sdf"],
            rtg=rtg_path(lab),
            threads=threads,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    try:
        dosage = propose_dosage(measured, grading.dosage.tolerance)
        bound = propose_af_bound(vafs)
    except ValueError as e:
        raise RunRefused(f"the controls set no bounds: {e}") from e
    hg004 = read_record(lab, "hg004")
    f1 = None if hg004 is None else hg004.snp_f1
    document = {
        "controls": {name: asdict(d) for name, d in measured.items()},
        "hg004_snp_f1": f1,
        "true_hets": len(vafs),
        "true_hets_below_0.30": share_below(vafs, 0.30),
        "dosage": dosage,
        "af_guard": {"min_het_vaf": bound},
    }
    text = _report(measured, f1, vafs, dosage, bound)

    def write(staging: Path) -> None:
        (staging / "report.json").write_text(json.dumps(document, indent=1))
        (staging / REPORT).write_text(text)

    write_table(lab, "foundation", write, foundation_versions(lab))
    print(f"foundation table: {folder(lab, 'foundation')}")
