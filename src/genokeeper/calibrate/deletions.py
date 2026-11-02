"""Slice 3b's calibration: each deletion class's chance of being false, measured on HG002's
two libraries (`pcrfree`, `pcr`) against GIAB's v5.0q structural-variant truth, into the table
every run grades shown deletions with (design decisions 1 and 5, 2026-10-27; option D).

HG002 has no true deletion over the listed genes' exons to count, so a class's chance is not
a precision in the shown set: it is `genokeeper.calibrate.classes`'s model (the one every
structural table shares) over 3b's size bands and `FEATURE_VALUES`, with gnomAD-SV's rare
deletions over the targets (`expected_true_carried`) as each band's prior.

The table (`DeletionTable`, `<lab>/calibration/deletions/`) records the versions it was
measured under (`deletion_versions`). A run computes its own and, when the table is missing or
was made with others, runs no deletion step and states the screen unassessed with the reason
`open_deletions` gives. The check's anchor comes from the run's own computation (the image
names written in the code, the sha256 of the exon resource, the region BEDs, gnomAD-SV and
the FASTA, the thresholds and the Truvari options), which a stale table cannot move because
the run never reads its versions from the table. This module imports nothing from `pipeline`,
so `pipeline` may import it.
"""

import hashlib
import json
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from importlib.resources import files
from pathlib import Path
from typing import Any

from genokeeper.calibrate import classes
from genokeeper.calibrate.classes import (
    FORMAT,
    FORMAT_KEY,
    OUTSIDE_VERSION,
    BandTerms,
    ClassTable,
    Features,
    Geometry,
    Likelihoods,
)
from genokeeper.calibrate.tables import folder
from genokeeper.checksum import sha256_file
from genokeeper.config import DeletionThresholds, Grading
from genokeeper.containers import IMAGES
from genokeeper.deletions.callers import caller_filter
from genokeeper.deletions.exons import RESOURCE
from genokeeper.deletions.files import sv_files
from genokeeper.deletions.screen import CallClass, PopulationSv, Screen, key
from genokeeper.deletions.truvari import AGAINST_TRUTH, BETWEEN_CALLERS
from genokeeper.regions import Interval, giab38_dir, merge
from genokeeper.tools import run_tool

# The run's three deletion steps (`pipeline.deletion_steps`): bump when a step's output
# changes. Here, not in `pipeline`, because the table records them.
DELLY_STEP_VERSION = "1"
DYSGU_STEP_VERSION = "2"  # 2: one run per primary chromosome
DELETIONS_UNION_STEP_VERSION = "1"
DELETION_IMAGES = ("delly", "dysgu", "truvari")
BANDS = ("long", "short")
# Every value each class feature takes (`screen.classify`), as `feature_value` names it.
FEATURE_VALUES: dict[str, tuple[str, ...]] = {
    "evidence": ("both", "both+depth", "dysgu", "delly", "delly+depth"),
    "repeat": ("true", "false"),
    "masked": ("true", "false"),
    "junction": ("true", "false"),  # always false for a deletion
}
# What `read_population` keeps of gnomAD-SV besides `AN`: PASS records of the type asked
# for (design decision 2: deletions).
_PASS_TYPE = 'FILTER="PASS" && INFO/SVTYPE="{svtype}"'


def deletions_dir(lab: Path) -> Path:
    """The deletion table `genokeeper calibrate` writes (`tables.folder`)."""
    return folder(lab, "deletions")


@dataclass(frozen=True)
class DeletionTable(ClassTable):
    """3b's class table: its schema is `FEATURE_VALUES`."""

    features: Features = field(default_factory=lambda: dict(FEATURE_VALUES))


def band_of(length: int, t: DeletionThresholds) -> str | None:
    """A deletion's band by its length, as `screen.classify` gives it; None under
    `min_length`."""
    if length < t.min_length:
        return None
    return "long" if length >= t.long_from else "short"


def _band(c: CallClass) -> str:
    return c.band


def measure(
    calls: Mapping[str, Sequence[tuple[CallClass, int, bool]]],
    truth: Mapping[str, Sequence[tuple[int, bool]]],
    g: Geometry,
    priors: Mapping[str, float],
    t: DeletionThresholds,
) -> tuple[dict[str, BandTerms], dict[str, Likelihoods]]:
    """`classes.measure` over `BANDS` (one geometry for both) and `FEATURE_VALUES`, a truth
    deletion's band from its length (`band_of`: under `min_length`, none)."""
    return classes.measure(
        calls,
        {lib: [(band_of(n, t), ok) for n, ok in rows] for lib, rows in truth.items()},
        dict.fromkeys(BANDS, g),
        priors,
        FEATURE_VALUES,
        band=_band,
    )


def class_of(name: str) -> CallClass:
    """The class `screen.key` named."""
    fields = name.split(":")
    if (
        len(fields) != 5
        or fields[2] not in ("repeat", "norepeat")
        or fields[3] not in ("masked", "unmasked")
        or fields[4] not in ("junction", "nojunction")
    ):
        raise ValueError(f"not a class key: {name!r}")
    band, evidence, repeat, masked, junction = fields
    return CallClass(
        band, evidence, repeat == "repeat", masked == "masked", junction == "junction"
    )


def read_population(
    vcf: Path,
    spans: Iterable[Interval],
    min_an: int,
    work: Path,
    *,
    svtype: str = "DEL",
) -> list[PopulationSv]:
    """gnomAD-SV's PASS `SVTYPE=<svtype>` records with `AN` at least `min_an` overlapping any
    0-based half-open span, by `bcftools query -R` (natively, regions written under `work`;
    htslib places a symbolic record by its `END`, and `-R` keeps every record overlapping a
    region, not only those starting in one), with the overall `AF`; once each, in file order.
    An insertion is its point (`end = pos + 1`, whatever gnomAD's `END`) with `|SVLEN|`."""
    rows = sorted({(c, s + 1, e) for c, s, e in spans if e > s})
    if not rows:
        return []
    insertion = svtype == "INS"
    columns = "%CHROM\t%POS\t%ID\t%INFO/END\t%INFO/AF"
    if insertion:
        columns += "\t%INFO/SVLEN"
    with tempfile.TemporaryDirectory(dir=work, prefix=".gnomad-regions-") as folder:
        regions = Path(folder) / "regions.tsv"
        regions.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))
        text = run_tool(
            [
                "bcftools", "query", "-R", str(regions),
                "-i", f"{_PASS_TYPE.format(svtype=svtype)} && INFO/AN>={min_an}",
                "-f", f"{columns}\n", str(vcf),
            ]
        )  # fmt: skip
    found: dict[tuple[str, int, str], PopulationSv] = {}
    for line in text.splitlines():
        chrom, pos, vid, end, af, *svlen = line.split("\t")
        if af == ".":
            continue
        start = int(pos)
        length = abs(int(svlen[0])) if insertion and svlen[0] != "." else None
        found.setdefault(
            (chrom, start, vid),
            PopulationSv(
                vid,
                chrom,
                start,
                start + 1 if insertion else int(end),
                float(af.split(",")[0]),
                svtype=svtype,
                svlen=length,
            ),
        )
    return list(found.values())


def expected_true_carried(
    lab: Path,
    grading: Grading,
    s: Screen,
    band: str,
    *,
    svtype: str = "DEL",
    band_of: Callable[[int, DeletionThresholds], str | None] = band_of,
) -> float:
    """Σ 2·AF over gnomAD-SV's PASS records of `svtype` in `band` (by `band_of` of their
    length) with `AN` at least `min_gnomad_an` and `AF` at most `max_gnomad_af` whose bases
    touch a target of `s`: the rare true ones over the targets a genome carries. An
    insertion's length is its `SVLEN` (none: no band), its bases its point `[pos, pos + 1)`."""
    t = grading.deletions
    spans = [(x.chrom, x.start, x.end) for ts in s.targets.values() for x in ts]
    by_chrom: dict[str, list[Interval]] = {}
    for span in merge(spans):
        by_chrom.setdefault(span[0], []).append(span)
    total, gnomad = 0.0, sv_files(lab).gnomad
    population = read_population(
        gnomad, spans, t.min_gnomad_an, gnomad.parent, svtype=svtype
    )
    for p in population:
        touches = any(
            st < p.end and p.pos < en for _, st, en in by_chrom.get(p.chrom, ())
        )
        length = p.svlen if svtype == "INS" else p.end - p.pos
        if (
            touches
            and p.af <= t.max_gnomad_af
            and length is not None
            and band_of(length, t) == band
        ):
            total += 2 * p.af
    return total


def _exons_sha256() -> str:
    data = files("genokeeper").joinpath("resources", "deletions", RESOURCE).read_bytes()
    return hashlib.sha256(data).hexdigest()


def deletion_versions(lab: Path, fasta_sha: str, grading: Grading) -> dict[str, str]:
    """The versions a deletion table is measured under: the three images, the sha256 of the
    exon resource, of the mask's two BEDs, the repeats' BED and gnomAD-SV, `[deletions]` as
    JSON without `likely_below` (a wording rule the table never reads), both Truvari option
    sets, `deletions_only`'s per-caller filters (`caller_filter`, as the run's step keys
    hold them), the three step versions, the FASTA's sha256, the class table's format
    (`classes.FORMAT_KEY`: a table of another format is stale) and the form of the `extra`
    `ClassTable.cautions` reads (`OUTSIDE_VERSION`); not HG002's truth, which the run never
    reads. `FileNotFoundError` when a file is absent."""
    giab, t = giab38_dir(lab), grading.deletions
    return {
        **{name: IMAGES[name] for name in DELETION_IMAGES},
        "exons": _exons_sha256(),
        "segdup": sha256_file(giab / grading.strata["segdup"]),
        "low_mappability": sha256_file(giab / grading.strata["low_mappability"]),
        "repeats": sha256_file(giab / grading.strata["low_complexity"]),
        "gnomad_sv": sha256_file(sv_files(lab).gnomad),
        "deletions": json.dumps(
            {k: v for k, v in asdict(t).items() if k != "likely_below"}, sort_keys=True
        ),
        "between_callers": " ".join(BETWEEN_CALLERS),
        "against_truth": " ".join(AGAINST_TRUTH),
        "delly_filter": caller_filter("delly", t.min_length),
        "dysgu_filter": caller_filter("dysgu", t.min_length),
        "delly_step": DELLY_STEP_VERSION,
        "dysgu_step": DYSGU_STEP_VERSION,
        "deletions_union_step": DELETIONS_UNION_STEP_VERSION,
        "fasta": fasta_sha,
        FORMAT_KEY: str(FORMAT),
        "extra": OUTSIDE_VERSION,
    }


# The geometry and `tables.write_table`'s writer of a deletion table, which build_deletions
# imports from here.
geometry = classes.geometry
write_table = classes.write_classes


def open_deletions(
    lab: Path, versions: Mapping[str, str]
) -> tuple[DeletionTable | None, str | None]:
    """The lab's deletion table when it was measured under `versions` (the run's own,
    `deletion_versions`); otherwise None and why (`classes.open_class_table`)."""
    return classes.open_class_table(lab, "deletions", versions, DeletionTable)


def report(table: DeletionTable) -> str:
    """`classes.report` over 3b's class keys (`class_of`)."""
    return classes.report(table, class_of)


def summary(
    calls: Sequence[tuple[CallClass, int, bool]],
    truth: Sequence[tuple[int, bool]],
    t: DeletionThresholds,
) -> dict[str, Any]:
    """One library's counts per band (`classes.summary`), a truth deletion's band from its
    length."""
    return classes.summary(
        calls,
        [(band_of(n, t), ok) for n, ok in truth],
        FEATURE_VALUES,
        BANDS,
        band=_band,
        key=key,
    )
