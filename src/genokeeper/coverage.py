"""Per-gene coverage and callability over MANE CDS (slice 3c).

The genes (`coverage_genes`): ACMG SF v3.3 and carrier Tier 3, the stage-3 genes (SMN1,
SMN2 and the `[repeat_loci]` genes) and the reference-limited `REFERENCE_LIMITED`, whose
copies the masked GRCh38 sets apart. Their CDS rows of a MANE Select or Plus Clinical
transcript on the primary assembly are cut from MANE v1.5's GTF into
`resources/coverage/cds.GRCh38.tsv` (3b's columns, `deletions.exons.COLUMNS`):

    python -m genokeeper.coverage cut GTF    # refuses any GTF but the pinned MANE v1.5

A gene's callable share is the share of its CDS +/- `[coverage] cds_margin` covered at
MAPQ >= `mapq` by `min_depth_per_copy` reads per expected copy (`expected_copies`: one for a
male's chrX outside the PARs and chrY, else two), from one `qc.mosdepth` call per copy
group; the shares inside GIAB's difficult strata and outside the calibrated benchmark are
reported beside it, never folded in.
"""

import bisect
import hashlib
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from importlib.resources import files
from pathlib import Path

from genokeeper.config import load_grading
from genokeeper.deletions.exons import (
    COLUMNS,
    MANE_RELEASE,
    MANE_SHA256,
    MANE_URL,
    Exon,
    cut_exons,
    read_resource,
)
from genokeeper.genelists import LISTED, GeneList, load_gene_lists
from genokeeper.regions import Interval, merge, subtract
from genokeeper.repeats.loci import RepeatLocus

# Genes the GRC false duplications on chr21 made unreadable in the unmasked GRCh38 (named by
# the acceptance; KCNE1, CRYAA and U2AF1 are on neither gene list).
REFERENCE_LIMITED = ("CBS", "KCNE1", "CRYAA", "U2AF1")
RESOURCE = "cds.GRCh38.tsv"
_SMN = ("SMN1", "SMN2")


def coverage_genes(
    lists: Mapping[str, GeneList], loci: Mapping[str, RepeatLocus]
) -> set[str]:
    """The `LISTED` lists' genes (Tier 4 enters only the carrier gate), the repeat loci's
    genes, SMN1 and SMN2, and `REFERENCE_LIMITED`."""
    genes = {gene for name in LISTED for gene in lists[name].genes}
    return genes | set(loci) | set(_SMN) | set(REFERENCE_LIMITED)


def load_cds() -> dict[str, list[Exon]]:
    """The packaged resource: gene -> its CDS blocks, in the file's order."""
    return read_resource(
        files("genokeeper").joinpath("resources", "coverage", RESOURCE).read_text()
    )


def write_resource(gtf: Path, out: Path) -> int:
    """Cut the coverage genes' MANE CDS rows from `gtf` into `out`; refuse any file but the
    pinned release (`MANE_SHA256`). Returns the number of rows written."""
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    if digest != MANE_SHA256:
        raise ValueError(
            f"{gtf} has sha256 {digest}, not {MANE_RELEASE}'s {MANE_SHA256}"
        )
    genes = coverage_genes(load_gene_lists(), load_grading().repeats.loci)
    rows = cut_exons(gtf, genes, feature="CDS")
    lines = [
        f"# source: {MANE_URL}",
        f"# release: {MANE_RELEASE}",
        f"# sha256: {digest}",
        (
            "# cut: every CDS row of a MANE_Select or MANE_Plus_Clinical transcript whose "
            "gene_name is in the ACMG SF v3.3 or ACMG carrier Tier 3 list, a [repeat_loci] "
            "gene, SMN1, SMN2, CBS, KCNE1, CRYAA or U2AF1, on chr1-22, chrX or chrY; start "
            "is 0-based half-open (BED), number is the transcript's own exon number"
        ),
        "\t".join(COLUMNS),
    ]
    lines += [
        f"{e.gene}\t{e.transcript}\t{e.number}\t{e.chrom}\t{e.start}\t{e.end}\t{tag}"
        for e, tag in rows
    ]
    out.write_text("\n".join(lines) + "\n")
    return len(rows)


@dataclass(frozen=True)
class CoverageRow:
    """One gene's coverage in a run (the store's `coverage` row without its run)."""

    gene: str
    bases: int
    callable: float
    mean_depth: float
    difficult: float
    outside_benchmark: float
    owner: str | None
    reliable: bool


def cds_bed(
    cds: Mapping[str, Sequence[Exon]],
    *,
    margin: int,
    contigs: Mapping[str, int],
    out: Path,
) -> dict[str, list[Interval]]:
    """Per gene its CDS +/- `margin`, merged and clipped to `contigs` (name -> length; a
    gene with nothing left is dropped), written to `out` as a BED named by gene."""
    genes: dict[str, list[Interval]] = {}
    for gene, blocks in sorted(cds.items()):
        widened = [
            (b.chrom, max(b.start - margin, 0), min(b.end + margin, contigs[b.chrom]))
            for b in blocks
            if b.chrom in contigs
        ]
        kept = merge(i for i in widened if i[1] < i[2])
        if kept:
            genes[gene] = kept
    rows = sorted((i, gene) for gene, intervals in genes.items() for i in intervals)
    out.write_text("".join(f"{c}\t{s}\t{e}\t{gene}\n" for (c, s, e), gene in rows))
    return genes


def expected_copies(
    genes: Mapping[str, list[Interval]],
    *,
    male: bool,
    par: Mapping[str, list[tuple[int, int]]],
) -> dict[str, int]:
    """1 for a male's gene on chrY or with bases on chrX outside `par` (grading's 1-based
    inclusive PARs), else 2."""
    pars = [(chrom, s - 1, e) for chrom, spans in par.items() for s, e in spans]

    def one_copy(intervals: list[Interval]) -> bool:
        if any(chrom == "chrY" for chrom, _, _ in intervals):
            return True
        on_x = [i for i in intervals if i[0] == "chrX"]
        return bool(subtract(on_x, pars))

    return {gene: 1 if male and one_copy(i) else 2 for gene, i in genes.items()}


class _Index:
    """Merged regions per chromosome, for the bases of an interval set inside them."""

    def __init__(self, regions: Iterable[Interval]) -> None:
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        for chrom, start, end in merge(regions):
            self._starts.setdefault(chrom, []).append(start)
            self._ends.setdefault(chrom, []).append(end)

    def bases_inside(self, intervals: Iterable[Interval]) -> int:
        total = 0
        for chrom, start, end in merge(intervals):
            starts, ends = self._starts.get(chrom, []), self._ends.get(chrom, [])
            i = bisect.bisect_right(ends, start)
            while i < len(starts) and starts[i] < end:
                total += min(end, ends[i]) - max(start, starts[i])
                i += 1
        return total


def _bases(intervals: Iterable[Interval]) -> int:
    return sum(end - start for _, start, end in merge(intervals))


def _share(intervals: Sequence[Interval], index: _Index) -> float:
    bases = _bases(intervals)
    return index.bases_inside(intervals) / bases if bases else 0.0


def share_inside(intervals: Sequence[Interval], regions: Iterable[Interval]) -> float:
    """The share of `intervals`' bases inside `regions`; 0.0 for no bases."""
    return _share(intervals, _Index(regions))


def coverage_rows(
    depth: Mapping[str, tuple[float, float]],
    genes: Mapping[str, list[Interval]],
    *,
    difficult: Iterable[Interval],
    benchmark: Iterable[Interval],
    owners: Mapping[str, str],
    min_callable: float,
) -> list[CoverageRow]:
    """One row per gene of `genes`, sorted: `depth` is `qc.mosdepth`'s (mean, share at the
    gene's own per-copy threshold); `reliable` when that share reaches `min_callable`."""
    hard, calibrated = _Index(difficult), _Index(benchmark)
    rows: list[CoverageRow] = []
    for gene, intervals in sorted(genes.items()):
        mean, share = depth[gene]
        rows.append(
            CoverageRow(
                gene=gene,
                bases=_bases(intervals),
                callable=share,
                mean_depth=mean,
                difficult=_share(intervals, hard),
                outside_benchmark=1.0 - _share(intervals, calibrated),
                owner=owners.get(gene),
                reliable=share >= min_callable,
            )
        )
    return rows


def coverage_fields(run_id: str, row: CoverageRow) -> dict[str, object]:
    """The store's `coverage` row."""
    return {"run_id": run_id, **asdict(row)}


def main(argv: Sequence[str] | None = None) -> int:
    """`cut GTF` writes the packaged resource from MANE's GTF."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) == 2 and args[0] == "cut":
        out = Path(__file__).resolve().parent / "resources" / "coverage" / RESOURCE
        out.parent.mkdir(parents=True, exist_ok=True)
        print(f"{write_resource(Path(args[1]), out)} CDS rows -> {out}")
        return 0
    print("usage: python -m genokeeper.coverage cut GTF", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
