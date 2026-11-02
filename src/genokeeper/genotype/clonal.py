"""The clonal genes' spans, which the allele-fraction guard reads (slice 3c).

Clonal haematopoiesis puts a somatic variant in a share of the blood cells, so in the genes
it mutates most (`[af_guard] clonal_genes`) a heterozygous call at a low allele fraction may
be clonal, not inherited (`reliability.guard_allele_fraction`). Each gene's span runs from
the first base of its MANE Select transcript's first exon to the last base of its last,
cut from NCBI's MANE GTF by `python -m genokeeper.genotype.clonal cut GTF` into
`resources/af_guard/clonal_spans.GRCh38.tsv` (header: `# source:`, `# release:`,
`# sha256:`), after checking the GTF is the pinned release.

The pin: the GTF's anchor comes from `MANE_SHA256` (`deletions.exons`), which another
release or a file cut short cannot move because the comparison runs on the bytes the cut
reads.
"""

import hashlib
import sys
from bisect import bisect_right
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from genokeeper.config import load_grading
from genokeeper.deletions.exons import MANE_RELEASE, MANE_SHA256, MANE_URL, cut_exons
from genokeeper.variants.key import ucsc_contig

RESOURCE = "clonal_spans.GRCh38.tsv"
COLUMNS = ("gene", "transcript", "chrom", "start", "end")


@dataclass(frozen=True)
class Span:
    """A gene's MANE Select span; `start` and `end` are 0-based half-open, as BED."""

    gene: str
    transcript: str
    chrom: str
    start: int
    end: int


def cut_spans(gtf: Path, genes: Collection[str]) -> list[Span]:
    """Per gene of `genes`, the span of its MANE Select exons on a primary chromosome,
    sorted by gene."""
    exons: dict[str, list[tuple[str, str, int, int]]] = {}
    for exon, tag in cut_exons(gtf, genes):
        if tag == "MANE_Select":
            row = (exon.transcript, exon.chrom, exon.start, exon.end)
            exons.setdefault(exon.gene, []).append(row)
    spans: list[Span] = []
    for gene in sorted(exons):
        rows = exons[gene]
        transcripts = {(t, c) for t, c, _, _ in rows}
        if len(transcripts) != 1:
            raise ValueError(f"{gene}: MANE Select on {sorted(transcripts)}, not one")
        ((transcript, chrom),) = transcripts
        start = min(s for _, _, s, _ in rows)
        end = max(e for _, _, _, e in rows)
        spans.append(Span(gene, transcript, chrom, start, end))
    return spans


def write_resource(gtf: Path, out: Path) -> int:
    """Cut the clonal genes' spans from `gtf` into `out`; refuse any file but the pinned
    release, and a gene the GTF has no MANE Select transcript for. Returns the number of
    spans written."""
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    if digest != MANE_SHA256:
        raise ValueError(
            f"{gtf} has sha256 {digest}, not {MANE_RELEASE}'s {MANE_SHA256}"
        )
    genes = load_grading().af_guard.clonal_genes
    spans = cut_spans(gtf, genes)
    missing = sorted(set(genes) - {s.gene for s in spans})
    if missing:
        raise ValueError(f"{gtf}: no MANE Select transcript for {', '.join(missing)}")
    lines = [
        f"# source: {MANE_URL}",
        f"# release: {MANE_RELEASE}",
        f"# sha256: {digest}",
        (
            "# cut: per [af_guard] clonal_genes gene, the span of its MANE_Select "
            "transcript's exon rows, first exon start to last exon end; start is 0-based "
            "half-open (BED)"
        ),
        "\t".join(COLUMNS),
    ]
    lines += [f"{s.gene}\t{s.transcript}\t{s.chrom}\t{s.start}\t{s.end}" for s in spans]
    out.write_text("\n".join(lines) + "\n")
    return len(spans)


def load_spans() -> dict[str, Span]:
    """The packaged resource: gene -> its span."""
    text = files("genokeeper").joinpath("resources", "af_guard", RESOURCE).read_text()
    spans: dict[str, Span] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if tuple(fields) == COLUMNS:
            continue
        gene, transcript, chrom, start, end = fields
        spans[gene] = Span(gene, transcript, chrom, int(start), int(end))
    return spans


@dataclass(frozen=True)
class ClonalIndex:
    """The listed clonal genes' spans: contig -> sorted, merged starts and ends."""

    intervals: dict[str, tuple[list[int], list[int]]]

    @classmethod
    def load(cls, genes: Collection[str]) -> ClonalIndex:
        """The packaged spans of `genes`; a gene the resource lacks raises `ValueError`
        (a gene named in grading.toml but never cut would be guarded nowhere)."""
        spans = load_spans()
        missing = sorted(set(genes) - set(spans))
        if missing:
            raise ValueError(
                f"no clonal span for {', '.join(missing)}; rerun "
                "`python -m genokeeper.genotype.clonal cut <MANE GTF>`"
            )
        raw: dict[str, list[tuple[int, int]]] = {}
        for gene in genes:
            span = spans[gene]
            raw.setdefault(span.chrom, []).append((span.start, span.end))
        intervals: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, pairs in raw.items():
            starts: list[int] = []
            ends: list[int] = []
            for start, end in sorted(pairs):
                if starts and start <= ends[-1]:
                    ends[-1] = max(ends[-1], end)
                else:
                    starts.append(start)
                    ends.append(end)
            intervals[chrom] = (starts, ends)
        return cls(intervals)

    def at(self, chrom: str, pos: int) -> bool:
        """Whether the 1-based position `pos` on `chrom` lies in a listed gene's span."""
        starts, ends = self.intervals.get(ucsc_contig(chrom), ([], []))
        i = bisect_right(starts, pos - 1) - 1
        return i >= 0 and pos - 1 < ends[i]


def main(argv: Sequence[str] | None = None) -> int:
    """`cut GTF` writes the packaged resource from MANE's GTF."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) == 2 and args[0] == "cut":
        out = (
            Path(__file__).resolve().parent.parent / "resources" / "af_guard" / RESOURCE
        )
        out.parent.mkdir(parents=True, exist_ok=True)
        print(f"{write_resource(Path(args[1]), out)} spans -> {out}")
        return 0
    print("usage: python -m genokeeper.genotype.clonal cut GTF", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
