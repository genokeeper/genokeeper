"""A `RunContext` for tests of sources that grade without the variant evidence."""

from pathlib import Path

from fixtures.faf_mini import write_faf_table
from genokeeper.config import load_grading
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import load_gene_lists
from genokeeper.sources import RunContext


def run_context(tmp_path: Path, *, male: bool = True) -> RunContext:
    return RunContext(
        run_id="run1",
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=male,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
    )
