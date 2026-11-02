"""Slice 3c: small-variant cells measured on HG002 where HG003's truth does not reach.

HG003's v4.2.1 truth covers chr1-22 and leaves out much of the genome GIAB's HG002 v5.0q
benchmark now holds. HG002 PCR-free's own gVCF (`read_called(lab, "pcrfree")`, called by the
run's pipeline) is scored as HG003's is (`measure.count_cells`) in three parts
(`cell_regions`), each the v5.0q benchmark minus the masked GRC false duplications:

- `haploid`: chrX and chrY outside the PARs, where a male's calls are haploid;
- `only`: the rest outside HG003's benchmark (a male's PARs among it, diploid);
- `shared`: inside HG003's benchmark, never used to grade, only to compare: a cell where
  HG002 scores better than HG003's interval allows may be one where HG002 was in
  DeepVariant's training (`training_overlap`), which the report names.

Each part is a `CalibrationTable` folder with its `regions.bed` under `cells_dir`, written
whole or not at all. A run opens the haploid and only parts (`open_cells`) beside HG003's
table, and `reliability.grade_called` tries them, by ploidy, outside HG003's benchmark; a
part missing or measured under other versions is not used, and the run says why, never
refusing. The version check's anchor comes from the versions the run computes itself
(`pipeline.calibration_versions`), which a part measured under another pipeline cannot move
because the run never reads its own versions from the part.
"""

import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

from genokeeper.calibrate import BENCH_BED, layout, table_dir
from genokeeper.calibrate import tables as registry
from genokeeper.calibrate.measure import measure
from genokeeper.calibrate.stage3 import hg002_truth
from genokeeper.calibrate.table import GENOTYPE, MANIFEST, CalibrationTable, CellKey
from genokeeper.config import load_grading
from genokeeper.genotype.reliability import ExtraCells, StrataIndex
from genokeeper.masking import masked_bed
from genokeeper.regions import Interval, merge, read_bed, subtract

HAPLOID = "haploid"
ONLY = "only"
SHARED = "shared"
PARTS = (HAPLOID, ONLY, SHARED)
REGIONS = "regions.bed"
NOT_MEASURED = "HG002 cells not measured"
_SEX_CONTIGS = ("chrX", "chrY")


def cells_dir(lab: Path) -> Path:
    """Where HG002's cells live (`tables.folder`): one folder per part."""
    return registry.folder(lab, "hg002_cells")


def _intersect(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    return subtract(a, subtract(a, b))


def cell_regions(
    v5: Sequence[Interval],
    hg003: Sequence[Interval],
    masked: Sequence[Interval],
    par: Mapping[str, list[tuple[int, int]]],
) -> dict[str, list[Interval]]:
    """The three parts' regions (see the module); `par` is grading.toml's [par], 1-based
    inclusive."""
    pars = [
        (contig, start - 1, end)
        for contig, spans in par.items()
        for start, end in spans
    ]
    sex = [i for i in v5 if i[0] in _SEX_CONTIGS]
    haploid = subtract(sex, [*pars, *masked])
    return {
        HAPLOID: haploid,
        ONLY: subtract(v5, [*hg003, *haploid, *masked]),
        SHARED: subtract(_intersect(merge(v5), merge(hg003)), masked),
    }


def training_overlap(
    shared: CalibrationTable, hg003: CalibrationTable
) -> list[CellKey]:
    """The genotype cells where HG002's rate in the shared part is below the lower end of
    HG003's interval for the same cell, sorted."""
    listed: list[CellKey] = []
    for key in shared.measured(errors=GENOTYPE):
        try:
            theirs = hg003.rate_for(key)
        except KeyError:
            continue
        if shared.rate_for(key).rate < theirs.lo:
            listed.append(key)
    return listed


def lab_regions(lab: Path) -> dict[str, list[Interval]]:
    """The three parts' regions from the lab's v5.0q benchmark, HG003's and the mask."""
    _, v5 = hg002_truth(lab)
    hg003 = layout(lab)["truth"] / BENCH_BED
    return cell_regions(
        read_bed(v5), read_bed(hg003), read_bed(masked_bed(lab)), load_grading().par
    )


def _write_bed(path: Path, intervals: Sequence[Interval]) -> Path:
    path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals))
    return path


def build_cells(
    lab: Path,
    *,
    gvcf: Path,
    sdf: Path,
    strata: StrataIndex,
    versions: Mapping[str, str],
    rtg: Path,
    threads: int,
) -> tuple[dict[str, CalibrationTable], list[CellKey]]:
    """Score `gvcf` (HG002 PCR-free's) against v5.0q's truth in each part's regions, into
    one table per part with the pipeline's `versions`; and the shared cells HG002 may have
    been in training for, against the lab's HG003 table when there is one."""
    truth, _ = hg002_truth(lab)
    regions = lab_regions(lab)
    tables: dict[str, CalibrationTable] = {}
    work = Path(tempfile.mkdtemp(prefix=".hg002-cells-", dir=gvcf.parent))
    try:
        for part in PARTS:
            bed = _write_bed(work / f"{part}.bed", regions[part])
            tables[part] = measure(
                gvcf,
                truth,
                bed,
                sdf=sdf,
                strata=strata,
                versions=versions,
                rtg=rtg,
                threads=threads,
            )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    hg003 = table_dir(lab)
    overlap = (
        training_overlap(tables[SHARED], CalibrationTable.read(hg003))
        if (hg003 / MANIFEST).is_file()
        else []
    )
    return tables, overlap


def write_cells(
    folder: Path,
    tables: Mapping[str, CalibrationTable],
    regions: Mapping[str, Sequence[Interval]],
) -> None:
    """Each part of `tables` as a `CalibrationTable` folder under `folder`, with its
    `regions.bed`; `build_hg002_cells` writes them through `tables.write_table`, which swaps
    the folder in whole: a reader sees the old parts or the new ones, never half of them."""
    for part, table in tables.items():
        table.write(folder / part)
        _write_bed(folder / part / REGIONS, regions[part])


def open_cells(
    lab: Path, versions: Mapping[str, str]
) -> tuple[tuple[ExtraCells, ...], str | None]:
    """The haploid and only parts, when both were measured under `versions` (the run's
    own, `pipeline.calibration_versions`); otherwise `()` and why. Never a refusal: the
    run grades without them."""
    folder = cells_dir(lab)
    parts = ((HAPLOID, True), (ONLY, False))
    if not all((folder / part / MANIFEST).is_file() for part, _ in parts):
        return (), NOT_MEASURED
    extra: list[ExtraCells] = []
    differ: set[str] = set()
    for part, haploid in parts:
        try:
            table = CalibrationTable.read(folder / part)
            regions = StrataIndex.from_beds({part: folder / part / REGIONS})
        except (OSError, ValueError, KeyError) as e:
            return (), f"HG002 cells unreadable ({' '.join(str(e).split())})"
        differ |= {
            key
            for key in {*table.versions, *versions}
            if table.versions.get(key) != versions.get(key)
        }
        extra.append(ExtraCells(table=table, regions=regions, haploid=haploid))
    if differ:
        return (
            (),
            f"HG002 cells measured under other versions: {', '.join(sorted(differ))}",
        )
    return tuple(extra), None


def _cell(key: CellKey) -> str:
    return f"{key.row_type}/{key.allele_class}/{key.region}/GQ{key.gq_bin}"


def cells_report(
    tables: Mapping[str, CalibrationTable], overlap: Sequence[CellKey]
) -> str:
    """The report's line: the genotype cells measured in the haploid and only parts, and
    the shared cells HG002 may have been in training for."""
    listed = ", ".join(_cell(key) for key in overlap) or "none"
    return (
        f"HG002 cells: haploid {len(tables[HAPLOID].measured())},"
        f" v5.0q-only {len(tables[ONLY].measured())}; HG002 may have been in training:"
        f" {listed}"
    )
