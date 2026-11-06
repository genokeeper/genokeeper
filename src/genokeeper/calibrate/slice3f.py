"""Slice 3f's calibration table: what the trio's called libraries (HG002 PCR-free, HG003,
HG004; HG002's PCR library is aligned only, so never measured here) withhold or trust, into
the table every run grades slice 3f with (closing design, "Slice 3f"; plan Task 15).

The measures are pure: `genokeeper.calibrate.build.build_slice3f` runs the tools on the trio
and hands each its per-library results.

- HLA (`measure_hla`): a class I gene is withheld when OptiType's or T1K's call on HG002
  PCR-free is not HG002's published type (`resources/hla/hg002_truth.tsv`), or a parent does
  not carry the allele HG002 got from it (no HG003 or HG004 types are published, so that is
  the trio check there is); a class II gene by T1K alone. `DQ` is withheld when a tag SNP is
  called differently from a truth, or (R2) when no library's benchmark covers a tag with a
  matching call: an ungraded tag never backs the coeliac rule-out.
- UGT1A1*28 (`measure_ugt1a1`): trusted only when every library's site is inside its
  benchmark and its call agrees with its truth, and some truth carries TA7.
- Phase (`measure_phase`): trusted with no cis/trans disagreement against HG002's phased
  truth over at least `MIN_PHASE_PAIRS` pairs.
- chrM (`measure_chrm_floor`, `measure_haplocheck`): HG003's chrM reads mixed into HG002
  PCR-free's at levels near `GRID`, each keyed by its level from the two parts' read counts
  (`mixture_level`); the floor is the lowest grid step from which every mixture calls
  HG003's homoplasmies within half to twice its level, raised above the highest false
  heteroplasmy of the unmixed libraries. haplocheck missing a mixture at or above the run's
  floor withholds the heteroplasmies.
- HHV-6 (`measure_hhv6`): the background ceiling is the highest ratio times
  `[hhv6] background_factor`; one reaching `min_ratio` withholds HHV-6.
- ROH (`measure_roh`): our runs against the truth's in both directions, and an artificial
  isodisomy (`SPIKE`) the calls must find.
- The Neanderthal share (`measure_archaic`, R5): each library within ±30% of the published
  mean (`PUBLISHED_ARCHAIC`).
- Blood groups (`measure_bloodgroups`, R8: systems only): a system RBCeq2 types differently
  on our calls and on the truth's, or measured on no library, is withheld.

`withheld` is keyed as each source reads it: an HLA gene or `DQ_KEY`, a blood-group system,
`CHRM_KEY`, `HHV6_KEY`, `ROH_KEY`, `ARCHAIC_KEY`. The versions (`slice3f_versions`) cover
every input a measure reads; `tables.open_table` compares them with the run's own, and its
anchor is that computation, which a table measured under others cannot move.
"""

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from genokeeper.bloodgroups import SystemCall
from genokeeper.calibrate import BENCH_BED, TRUTH_VCF, hg004_truth, layout, tables
from genokeeper.calibrate.records import hg002_dir, read_record
from genokeeper.calibrate.stage3 import hg002_truth
from genokeeper.calibrate.tables import RUN, Current, open_table
from genokeeper.checksum import sha256_cached, sha256_file
from genokeeper.config import (
    Grading,
    Hhv6Thresholds,
    RohThresholds,
    grch38_fasta,
    load_grading,
)
from genokeeper.containers import IMAGES
from genokeeper.haplocheck import Haplocheck
from genokeeper.mito import MitoCall
from genokeeper.pgx.hla import CLASS_I, OPTITYPE, T1K, HlaCall, shown, t1k_index
from genokeeper.pgx.pharmcat import positions_path
from genokeeper.roh import Run
from genokeeper.typing_files import typing_files

# HG002 PCR-free, HG003, HG004 (`calibrate.records.RECORDED`; `pcr` is aligned only).
LIBRARIES: tuple[str, ...] = ("pcrfree", "hg003", "hg004")
SLICE3F_IMAGES = (
    "optitype",
    "t1k",
    "whatshap",
    "mutserve",
    "haplocheck",
    "rbceq2",
    "fastmixture",
    "hmmix",
    "truvari",
)
# Bump when a measure's arithmetic changes: a table measured by other code is stale.
MEASURE_VERSION = "3"
NAME = "slice3f"
TABLE = "table.json"
FORMAT = 1
RESOURCES = Path(__file__).resolve().parents[1] / "resources"
RESOURCE_DIRS = ("hla", "mito", "bloodgroups", "traits")
IDENTITIES = "slice3f.identities.json"

DQ_KEY = "DQ"
CHRM_KEY = "chrM heteroplasmy"
HHV6_KEY = "HHV-6"
ROH_KEY = "ROH"
ARCHAIC_KEY = "neanderthal"

# Which parent each of HG002's truth alleles comes from (maternal first).
PARENTS = (("hg004", "HG004", "maternal"), ("hg003", "HG003", "paternal"))
MIN_PHASE_PAIRS = 20
GRID = (0.02, 0.05, 0.10, 0.20)
# A level more than this factor from every grid step has no step.
GRID_TOLERANCE = 1.5
# The spike: 12 Mb of chr2's long arm, away from the centromere.
SPIKE = ("chr2", 50_000_000, 62_000_000)
SPIKE_COVER = 0.9
# A run counts as matching another when one of at least this length covers half of it.
MATCHING_RUN = 1_000_000
# Coll Macià et al. 2026's Neanderthal share of NA20810 by hmmix (R5), and the band.
PUBLISHED_ARCHAIC = 0.0118
ARCHAIC_BAND = 0.3


# --- HLA ----------------------------------------------------------------------------------


def _matches(called: str | None, truth: str) -> bool:
    """A called allele (equivalent names joined by `,`) matches a truth allele (names
    joined by `,`, each a match) when every name it could be is one of the truth's."""
    return called is not None and set(called.split(",")) <= set(truth.split(","))


def measure_hla(
    calls: Mapping[str, tuple[Mapping[str, HlaCall], Mapping[str, HlaCall]]],
    truth: Mapping[str, HlaCall],
    tags: Mapping[str, Mapping[str, bool | None]],
) -> dict[str, str]:
    """The HLA genes and `DQ_KEY` withheld, with the reason. `calls`: library ->
    (OptiType's, T1K's calls); `truth`: HG002's published types, maternal first; `tags`:
    library -> tag rsid -> True (inside the library's benchmark, our call matching its
    truth), False (inside, differing), None (outside, or not called: not graded there)."""
    withheld: dict[str, str] = {}
    for gene, published in truth.items():
        # A truth allele missing matches nothing, so the gene is withheld.
        maternal, paternal = (a or "?" for a in published)
        callers = [(OPTITYPE, 0), (T1K, 1)] if gene in CLASS_I else [(T1K, 1)]
        for caller, i in callers:
            call = calls["pcrfree"][i].get(gene)
            right = call is not None and (
                (_matches(call[0], maternal) and _matches(call[1], paternal))
                or (_matches(call[0], paternal) and _matches(call[1], maternal))
            )
            if not right:
                withheld.setdefault(
                    gene,
                    f"not covered: {caller} called HG002's {gene} {shown(call)}, not its"
                    f" published {maternal}/{paternal}, in calibration",
                )
                continue
            for (parent, name, side), allele in zip(
                PARENTS, (maternal, paternal), strict=True
            ):
                found = calls[parent][i].get(gene)
                if found is None or not any(_matches(a, allele) for a in found):
                    withheld.setdefault(
                        gene,
                        f"not covered: {caller} on {name} ({shown(found)}) does not"
                        f" carry HG002's {side} {gene}{allele} in calibration",
                    )
    rsids = sorted({rsid for found in tags.values() for rsid in found})
    differ = sorted(
        {rsid for found in tags.values() for rsid, ok in found.items() if ok is False}
    )
    ungraded = [r for r in rsids if not any(found.get(r) for found in tags.values())]
    if differ:
        withheld[DQ_KEY] = (
            f"not covered: DQ tag SNP {', '.join(differ)} called differently from the"
            " truth in calibration"
        )
    elif ungraded:
        withheld[DQ_KEY] = (
            f"not covered: DQ tag SNP {', '.join(ungraded)} not graded: no calibration"
            " benchmark covers it with a matching call"
        )
    return withheld


# --- UGT1A1*28 and phase ------------------------------------------------------------------


def ta_genotype(alleles: Sequence[int]) -> str:
    """The TA repeat's genotype from each allele's repeat count, sorted: `TA6/TA7`."""
    return "/".join(f"TA{n}" for n in sorted(alleles))


def _ordered(genotype: str | None) -> str | None:
    return None if genotype is None else "/".join(sorted(genotype.split("/")))


def measure_ugt1a1(
    ours: Mapping[str, str | None],
    truth: Mapping[str, str | None],
    inside: Mapping[str, bool],
) -> bool:
    """Whether rs3064744's call is trusted: every library inside its benchmark, its call
    (`ta_genotype`; None a no-call) the truth's, and some truth carrying TA7 (a measure
    that saw only reference genotypes proves nothing about TA7)."""
    for lib in LIBRARIES:
        call = _ordered(ours.get(lib))
        if (
            not inside.get(lib, False)
            or call is None
            or call != _ordered(truth.get(lib))
        ):
            return False
    return any("TA7" in (truth.get(lib) or "").split("/") for lib in LIBRARIES)


def measure_phase(agree: int, disagree: int) -> bool:
    """Whether read-backed phase is trusted: no cis/trans disagreement with HG002's phased
    truth over at least `MIN_PHASE_PAIRS` compared pairs."""
    return disagree == 0 and agree >= MIN_PHASE_PAIRS


# --- chrM ---------------------------------------------------------------------------------


def mixture_level(n_hg002: int, n_hg003: int) -> float:
    """A mixture's expected HG003 level from the two parts' chrM read counts."""
    return n_hg003 / (n_hg002 + n_hg003)


def grid_step(level: float) -> float | None:
    """The `GRID` step nearest `level` on a log scale; None when the level is more than
    `GRID_TOLERANCE` times from every step."""
    if level <= 0:
        return None
    step = min(GRID, key=lambda s: abs(math.log(level / s)))
    return step if abs(math.log(level / step)) <= math.log(GRID_TOLERANCE) else None


def _passes(level: float, calls: Sequence[MitoCall], expected: frozenset[int]) -> bool:
    found = {c.pos: c.level for c in calls}
    return all(
        pos in found and level / 2 <= found[pos] <= level * 2 for pos in expected
    )


def measure_chrm_floor(
    mixtures: Mapping[float, Sequence[MitoCall]],
    expected: frozenset[int],
    false_levels: Sequence[float],
) -> float | None:
    """The chrM heteroplasmy floor: the lowest grid step from which every mixture (keyed
    by its level, `mixture_level`; one with no `grid_step` is not counted) calls every
    `expected` position within half to twice its level, raised to the grid step above the
    highest of `false_levels` (heteroplasmies in the unmixed libraries). None when a false
    level reaches the top step, no mixture passes, or nothing is `expected`."""
    passing: dict[float, bool] = {}
    for level, calls in mixtures.items():
        step = grid_step(level)
        if step is not None:
            passing[step] = passing.get(step, True) and _passes(level, calls, expected)
    lowest = next(
        (
            step
            for step in sorted(passing)
            if all(ok for s, ok in passing.items() if s >= step)
        ),
        None,
    )
    top = max(false_levels, default=0.0)
    if lowest is None or top >= GRID[-1] or not expected:
        return None
    above = min(s for s in GRID if s > top)
    return max(lowest, above)


def measure_haplocheck(
    mixtures: Mapping[float, Haplocheck], floor: float
) -> str | None:
    """The heteroplasmies' withholding reason when haplocheck does not flag a mixture at or
    above `floor` (the run's); None otherwise. A flagged unmixed library only demotes, so it
    is reported (`report`), never withheld on."""
    missed = sorted(
        step
        for level, found in mixtures.items()
        if (step := grid_step(level) or level) >= floor and not found.contaminated
    )
    if not missed:
        return None
    return f"not covered: haplocheck missed a {missed[0]:.2f} mixture in calibration"


# --- HHV-6, ROH, the archaic share, blood groups ------------------------------------------


def measure_hhv6(
    ratios: Mapping[str, float], grading: Hhv6Thresholds
) -> tuple[float, str | None]:
    """The background ceiling (the highest ratio times `background_factor`), and a
    withholding reason when it reaches `min_ratio`."""
    lib, highest = max(ratios.items(), key=lambda item: item[1])
    ceiling = highest * grading.background_factor
    if ceiling < grading.min_ratio:
        return ceiling, None
    return ceiling, (
        f"not covered: HHV-6 background ratio {highest:.3g} on {lib.upper()} gives a"
        f" non-finding ceiling of {ceiling:.3g}, at or above the band's {grading.min_ratio}"
    )


def _overlap(a: Run, chrom: str, start: int, end: int) -> int:
    if a.chrom != chrom:
        return 0
    return max(0, min(a.end, end) - max(a.start, start) + 1)


def _covered(run: Run, others: Sequence[Run]) -> bool:
    """Whether a run of at least `MATCHING_RUN` among `others` overlaps half of `run`."""
    span = run.end - run.start + 1
    return any(
        o.length >= MATCHING_RUN
        and 2 * _overlap(o, run.chrom, run.start, run.end) >= span
        for o in others
    )


def measure_roh(
    ours: Mapping[str, Sequence[Run]],
    truth: Mapping[str, Sequence[Run]],
    spiked: Mapping[str, Sequence[Run]],
    spike: tuple[str, int, int],
    grading: RohThresholds,
) -> str | None:
    """ROH's withholding reason: a run of ours of `min_run` or more that no truth run of
    1 Mb or more overlaps by half its length (false), a truth run of `min_run` or more no
    run of ours of 1 Mb or more overlaps by half (missed), or a library's spiked calls with
    no single run of `isodisomy_min` or more covering 90% of `spike`; None otherwise."""
    for lib, runs in ours.items():
        for r in runs:
            if r.length >= grading.min_run and not _covered(r, truth.get(lib, ())):
                return (
                    f"not covered: a false run of homozygosity on {lib.upper()}"
                    f" ({r.chrom}:{r.start}-{r.end}) in calibration"
                )
    for lib, runs in truth.items():
        for t in runs:
            if t.length >= grading.min_run and not _covered(t, ours.get(lib, ())):
                return (
                    f"not covered: a run of homozygosity missed on {lib.upper()}"
                    f" ({t.chrom}:{t.start}-{t.end}) in calibration"
                )
    chrom, start, end = spike
    need = SPIKE_COVER * (end - start + 1)
    for lib, runs in spiked.items():
        if not any(
            r.length >= grading.isodisomy_min and _overlap(r, chrom, start, end) >= need
            for r in runs
        ):
            return (
                f"not covered: the artificial isodisomy (spike {chrom}:{start}-{end})"
                f" not found on {lib.upper()} in calibration"
            )
    return None


def measure_archaic(shares: Mapping[str, float], published: float) -> str | None:
    """The Neanderthal share's withholding reason when a library's share lies outside
    ±30% of the `published` mean; None otherwise."""
    for lib, share in shares.items():
        if abs(share - published) > ARCHAIC_BAND * published:
            return (
                f"not covered: Neanderthal share {share:.4f} on {lib.upper()} outside"
                f" ±30% of the published {published} in calibration"
            )
    return None


def measure_bloodgroups(
    ours: Mapping[str, Mapping[str, SystemCall | None]],
    from_truth: Mapping[str, Mapping[str, SystemCall | None]],
) -> dict[str, str]:
    """The blood-group systems withheld: library -> system -> RBCeq2's call on our calls
    and on the library's truth (None: not measurable there, a no-call position). A system
    whose genotypes differ on any library, measured on none, or typed by neither side on
    any, is withheld."""
    names = sorted(
        {s for found in (*ours.values(), *from_truth.values()) for s in found}
    )
    withheld: dict[str, str] = {}
    for system in names:
        measured = typed = 0
        for lib in LIBRARIES:
            mine = ours.get(lib, {}).get(system)
            theirs = from_truth.get(lib, {}).get(system)
            if mine is None or theirs is None:
                continue
            measured += 1
            typed += mine.genotype is not None or theirs.genotype is not None
            if mine.genotype != theirs.genotype:
                withheld[system] = (
                    f"not covered: RBCeq2 typed {system} {mine.genotype or 'no call'} on"
                    f" {lib.upper()}'s calls, {theirs.genotype or 'no call'} on its truth,"
                    " in calibration"
                )
                break
        else:
            if not measured:
                withheld[system] = (
                    f"not covered: {system} is measurable on no calibration library"
                )
            elif not typed:
                withheld[system] = f"not covered: {system} not typed in calibration"
    return withheld


def consistency_ancestry(
    admixture: Mapping[str, Mapping[str, float]],
) -> dict[str, Any]:
    """For the report only: each library's admixture proportions and, with all three,
    HG002's minus its parents' mean per group (a child's proportions are about the mean
    of its parents')."""
    found: dict[str, Any] = {"proportions": {k: dict(v) for k, v in admixture.items()}}
    if all(lib in admixture for lib in LIBRARIES):
        child, father, mother = (admixture[lib] for lib in LIBRARIES)
        found["child_minus_parents"] = {
            g: round(child[g] - (father.get(g, 0.0) + mother.get(g, 0.0)) / 2, 6)
            for g in child
        }
    return found


# --- the table ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Slice3fTable:
    """What the trio withholds (key -> reason, keyed as each source reads it), whether
    UGT1A1*28 and phase are trusted, the measured chrM floor (None with `CHRM_KEY`
    withheld) and HHV-6's background ceiling (None with `HHV6_KEY` withheld), and the
    report's numbers."""

    withheld: dict[str, str]
    ugt1a1_trusted: bool
    phase_trusted: bool
    chrm_floor: float | None
    hhv6_background: float | None
    report: dict[str, Any] = field(default_factory=dict[str, Any])

    def to_json(self) -> str:
        return json.dumps(
            {"format": FORMAT, **asdict(self)}, indent=1, sort_keys=True, default=str
        )

    @staticmethod
    def from_json(text: str) -> Slice3fTable:
        """The table `to_json` wrote; `ValueError` for another format."""
        data = json.loads(text)
        if data.pop("format", None) != FORMAT:
            raise ValueError(f"slice 3f table format is not {FORMAT}")
        return Slice3fTable(**data)


def write_slice3f(folder: Path, table: Slice3fTable) -> None:
    """Write `table` and its report into `folder`, the staging folder
    `tables.write_table` swaps in whole."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / TABLE).write_text(table.to_json())
    (folder / tables.REPORT).write_text(report(table) + "\n")


def open_slice3f(
    lab: Path, grading: Grading | None = None
) -> tuple[dict[str, str] | None, Slice3fTable | None, str | None]:
    """The lab's table with the versions it was checked against (`slice3f_versions`, the
    truths' digests as recorded: a run reads no calibration CRAM, gVCF or truth, ruling
    R22), or None for both and why: no table (`NOT_CALIBRATED`), an input those versions
    hash missing (`slice 3f not run: ...`, a typing file naming the bootstrap), or the
    table stale (`tables.open_table`'s reason) or unreadable."""
    path = tables.folder(lab, NAME) / tables.VERSIONS
    if not path.is_file():
        return None, None, tables.NOT_CALIBRATED
    try:
        versions = slice3f_versions(lab, grading, recorded=json.loads(path.read_text()))
    except (OSError, RuntimeError, ValueError) as e:  # RuntimeError: RunRefused
        return None, None, f"slice 3f not run: {' '.join(str(e).split())}"
    found = open_table(lab, NAME, versions=versions)
    if not isinstance(found, Current):
        return None, None, found.reason
    try:
        return (
            versions,
            Slice3fTable.from_json((found.folder / TABLE).read_text()),
            None,
        )
    except (OSError, ValueError, KeyError, TypeError) as e:
        return (
            None,
            None,
            f"slice 3f table unreadable ({' '.join(str(e).split())}); {RUN}",
        )


def report(table: Slice3fTable) -> str:
    """The table as text: what is withheld or trusted, then the report's numbers."""
    lines = ["withheld:"]
    lines += [f"  {key}: {why}" for key, why in sorted(table.withheld.items())] or [
        "  nothing"
    ]
    lines += [
        f"UGT1A1*28 trusted: {table.ugt1a1_trusted}",
        f"phase trusted: {table.phase_trusted}",
        f"chrM floor: {table.chrm_floor}",
        f"HHV-6 background ceiling: {table.hhv6_background}",
        "measured:",
        json.dumps(table.report, indent=1, sort_keys=True, default=str),
    ]
    return "\n".join(lines)


# --- versions -----------------------------------------------------------------------------


def truths(lab: Path) -> dict[str, tuple[Path, Path]]:
    """Each library's truth VCF and benchmark BED: HG002's v5.0q, HG003's and HG004's
    v4.2.1."""
    hg003 = layout(lab)["truth"]
    return {
        "pcrfree": hg002_truth(lab),
        "hg003": (hg003 / TRUTH_VCF, hg003 / BENCH_BED),
        "hg004": hg004_truth(lab),
    }


def calls_key(gvcf: Path) -> str:
    """A called library's calls key: its call step's key, the folder `pipeline._cached`
    names the gVCF's directory after."""
    return gvcf.parent.name


def resource_shas(root: Path | None = None) -> dict[str, str]:
    """The sha256 of every file under `resources/<dir>` for each of `RESOURCE_DIRS`,
    walked, keyed `resources/<path>`: a file added later joins by itself."""
    root = RESOURCES if root is None else root
    found: dict[str, str] = {}
    for name in RESOURCE_DIRS:
        for path in sorted((root / name).rglob("*")):
            if path.is_file() and "__pycache__" not in path.parts:
                key = f"resources/{path.relative_to(root).as_posix()}"
                found[key] = sha256_file(path)
    return found


def slice3f_versions(
    lab: Path,
    grading: Grading | None = None,
    *,
    recorded: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """The versions the table is measured under: the images, `MEASURE_VERSION`, the
    FASTA's, T1K's index's and PharmCAT's positions' sha256, every resource file's
    (`resource_shas`), the typing files' (`TypingFiles.required`), each library's align
    and call step keys from its record (never its CRAM or gVCF: "" without one), each
    truth VCF's and benchmark BED's sha256 (the run takes them from the table's
    `recorded` versions, ruling R22), the values of `[hhv6]`, `[roh]`, `[mito]` and
    `p_low`. `FileNotFoundError` for a missing file."""
    grading = grading if grading is not None else load_grading()
    identities = hg002_dir(lab).parent / IDENTITIES

    def sha(path: Path) -> str:
        if not path.is_file():
            raise FileNotFoundError(f"{path} is missing")
        return sha256_cached(path, identities)

    versions = {name: IMAGES[name] for name in SLICE3F_IMAGES}
    versions["measures"] = MEASURE_VERSION
    versions["fasta"] = sha(grch38_fasta(lab))
    versions["t1k_index"] = sha(t1k_index(lab))
    versions["pharmcat_positions"] = sha(positions_path(lab))
    versions |= resource_shas()
    versions |= {
        f"typing:{path.relative_to(lab).as_posix()}": sha(path)
        for path in typing_files(lab).required()
    }
    for lib in LIBRARIES:
        record = read_record(lab, lib)
        gvcf = None if record is None else record.gvcf
        versions[f"{lib}:cram"] = "" if record is None else record.key
        versions[f"{lib}:calls"] = "" if gvcf is None else calls_key(Path(gvcf))
    for lib, files in truths(lab).items():
        for key, path in zip((f"{lib}:truth", f"{lib}:benchmark"), files, strict=True):
            versions[key] = sha(path) if recorded is None else recorded.get(key, "")
    versions["[hhv6]"] = json.dumps(asdict(grading.hhv6), sort_keys=True)
    versions["[roh]"] = json.dumps(asdict(grading.roh), sort_keys=True)
    versions["[mito]"] = json.dumps(asdict(grading.mito), sort_keys=True)
    versions["p_low"] = str(grading.p_low)
    return versions
