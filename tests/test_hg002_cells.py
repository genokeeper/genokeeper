"""HG002's small-variant cells where HG003's truth does not reach (slice 3c): the regions
each part is measured in, and the cells HG002 may have scored better on for having been in
DeepVariant's training."""

from pathlib import Path

import pytest

from genokeeper.calibrate.build import build_hg002_cells
from genokeeper.calibrate.hg002_cells import (
    HAPLOID,
    ONLY,
    SHARED,
    cell_regions,
    cells_dir,
    training_overlap,
)
from genokeeper.calibrate.table import (
    GENOTYPE,
    VERSION_KEYS,
    CalibrationTable,
    CellKey,
    Counts,
)
from genokeeper.pipeline import RunRefused

VERSIONS = {key: f"{key}-value" for key in VERSION_KEYS}


def test_cell_regions() -> None:
    """haploid: v5.0q on chrX/chrY outside the PARs; only: v5.0q outside HG003's benchmark
    and the haploid part (a male's PAR is diploid); shared: v5.0q inside HG003's; the mask
    out of all three."""
    regions = cell_regions(
        [("chr1", 0, 100), ("chrX", 0, 3_000_000)],
        [("chr1", 0, 50)],
        [("chr1", 90, 95)],
        {"chrX": [(10001, 2781479)]},  # 1-based, inclusive
    )
    assert regions == {
        HAPLOID: [("chrX", 0, 10000), ("chrX", 2781479, 3_000_000)],
        ONLY: [("chr1", 50, 90), ("chr1", 95, 100), ("chrX", 10000, 2781479)],
        SHARED: [("chr1", 0, 50)],
    }


def _one_cell(errors: int, total: int) -> CalibrationTable:
    key = CellKey("variant", "snv", "easy", 50)
    return CalibrationTable.from_counts(
        {(GENOTYPE, key): Counts(errors, total, 30.0 * total)}, VERSIONS
    )


def test_training_overlap() -> None:
    """A shared cell where HG002's genotype rate is below HG003's interval is listed:
    HG002 may have been in DeepVariant's training. Equal rates are not."""
    hg003 = _one_cell(50, 10_000)
    key = CellKey("variant", "snv", "easy", 50)
    assert hg003.rate_for(key).lo > 1 / 10_000
    assert training_overlap(_one_cell(1, 10_000), hg003) == [key]
    assert training_overlap(_one_cell(50, 10_000), hg003) == []


def test_build_without_pcrfree_refuses(tmp_path: Path) -> None:
    """Without HG002 PCR-free called on this reference the cells are not measured: the
    builder refuses naming what is missing and the command that supplies it, and no
    folder is written."""
    lab = tmp_path / "lab"
    with pytest.raises(
        RunRefused, match=r"^HG002 cells not measured: .*not aligned.*align pcrfree`$"
    ):
        build_hg002_cells(lab, 2)
    assert not cells_dir(lab).exists()
