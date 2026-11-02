"""The MANE exons of the ACMG SF, carrier and screen genes, the targets the deletion screen
reads, and the MANE spans of slice 3d's paralog loci.

`cut_exons` takes every exon of a MANE Select or MANE Plus Clinical transcript of the listed
genes on the primary assembly (`PRIMARY`) from the MANE GTF: no caller shard runs on a fix
or alt contig, so an exon there could never be screened (MANE v1.5 places SLC37A4's only
transcript on `chr11_KN196481v1_fix`, so SLC37A4 is not covered);
`python -m genokeeper.deletions.exons cut GTF` writes them to
`resources/deletions/exons.GRCh38.tsv` (header: `# source:`, `# release:`, `# sha256:`),
after checking the GTF is the pinned release. `load_exons` reads that resource back. Exon numbers
are each transcript's own, so a gene with a Select and a Plus Clinical transcript has two
numberings. `targets` widens each gene's exons by a margin and merges those that touch, and
numbers each target by its MANE Select exons, or by the Plus Clinical ones where it has none
(`Target.plus_clinical`).

`cut_spans` takes each gene's MANE Select transcript span on the primary assembly (MANE
v1.5 also places a C4A transcript on `chr6_GL000254v2_alt`, which must not replace chr6's);
`python -m genokeeper.deletions.exons spans GTF` writes the spans of `PARALOG_SPAN_GENES`, and
`GBA1_ex9-11` (GBA1's exons 9 to 11, the part GBAP1 copies), to
`resources/paralogs/genes.GRCh38.tsv`, which `load_spans` reads back.

The pin: the GTF's anchor comes from `MANE_SHA256` written here, which another release or a
file cut short cannot move because the comparison runs on the bytes the cut reads.
"""

import gzip
import hashlib
import re
import sys
from collections.abc import Collection, Iterator, Mapping, Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from genokeeper.genelists import load_gene_lists, load_screen_genes

MANE_RELEASE = "MANE v1.5"
MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
    "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
)
MANE_SHA256 = "71d4b3c89c9d948683bf0db5d81b00d9ae9f3b943177b74cbd87e68f08e34d66"
MANE_TAGS = ("MANE_Select", "MANE_Plus_Clinical")
RESOURCE = "exons.GRCh38.tsv"
COLUMNS = ("gene", "transcript", "number", "chrom", "start", "end", "tag")
SCREEN_GENES_NOTE = "and the screen genes of genelists/screen_genes.tsv (DPYD, CHEK2)"
SPANS = "genes.GRCh38.tsv"
SPAN_COLUMNS = ("gene", "chrom", "start", "end")
PARALOG_SPAN_GENES = (
    "STRC", "CATSPER2", "RHD", "RHCE", "UGT2B17", "UGT2B15", "AMY1A", "AMY1B",
    "AMY1C", "GSTM1", "C4A", "C4B", "CYP2A6", "CYP2A7", "GBA1", "LPA", "CYP2D6", "SMN1",
    "SMN2",
)  # fmt: skip
# GBA1's exons 9 to 11 (MANE Select numbering), the part of the gene GBAP1 copies.
GBA1_HOMOLOGY = ("GBA1_ex9-11", "GBA1", (9, 10, 11))
# The CYP2D locus Cyrius owns: CYP2D6's MANE v1.5 Select span, chr22:[42,126,498,
# 42,130,810) 0-based, from 5 kb before it to 25 kb after it, which holds CYP2D7 and CYP2D8P
# (NCBI Gene 1564 and 1568, GRCh38 chr22:42,139,576-42,144,483 and 42,149,883-42,155,016,
# 1-based, read from E-utilities esummary on 2026-10-30).
OWNED_CYP2D = ("chr22", 42_121_498, 42_155_810)
# The primary chromosomes: the exons `cut_exons` keeps, and the shards the callers run
# (`pipeline.deletion_steps`), one container each; a deletion never spans two chromosomes.
PRIMARY = (*(f"chr{n}" for n in range(1, 23)), "chrX", "chrY")

_GENE_NAME = re.compile(r'gene_name "([^"]*)"')
_TRANSCRIPT = re.compile(r'transcript_id "([^"]*)"')
_NUMBER = re.compile(r"exon_number (\d+)")
_TAG = re.compile(r'tag "([^"]*)"')


@dataclass(frozen=True)
class Exon:
    """One exon of a MANE transcript; `start` and `end` are 0-based half-open, as BED;
    `select` False for a MANE Plus Clinical transcript's."""

    gene: str
    transcript: str
    number: int
    chrom: str
    start: int
    end: int
    select: bool = True


@dataclass(frozen=True)
class Target:
    """A region of a gene the screen reads: its exons widened by a margin and merged;
    `exons` are the numbers merged into it, the Plus Clinical transcript's where
    `plus_clinical` (it holds no Select exon), which may name another exon in Select's."""

    chrom: str
    start: int
    end: int
    exons: tuple[int, ...]
    plus_clinical: bool = False


def _mane_rows(
    gtf: Path, genes: Collection[str], feature: str = "exon"
) -> Iterator[tuple[Exon, str]]:
    """Each `feature` row (`exon`, or `CDS` for slice 3c's coverage) of a MANE transcript of
    a listed gene, with its MANE tag."""
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != feature:
                continue
            attributes = fields[8]
            name = _GENE_NAME.search(attributes)
            tags = [t for t in _TAG.findall(attributes) if t in MANE_TAGS]
            transcript = _TRANSCRIPT.search(attributes)
            number = _NUMBER.search(attributes)
            if not (name and tags and transcript and number) or name[1] not in genes:
                continue
            if fields[0] not in PRIMARY:
                continue
            exon = Exon(
                gene=name[1],
                transcript=transcript[1],
                number=int(number[1]),
                chrom=fields[0],
                start=int(fields[3]) - 1,
                end=int(fields[4]),
            )
            yield exon, tags[0]


def _order(item: tuple[Exon, str]) -> tuple[str, str, int, int, str, int]:
    e = item[0]
    return e.gene, e.chrom, e.start, e.end, e.transcript, e.number


def cut_exons(
    gtf: Path, genes: Collection[str], feature: str = "exon"
) -> list[tuple[Exon, str]]:
    """Every `feature` row (default `exon`) on a `PRIMARY` chromosome of a transcript tagged
    MANE_Select or MANE_Plus_Clinical whose `gene_name` is in `genes`, with that tag, sorted
    by gene, chromosome name, then start."""
    return sorted(_mane_rows(gtf, genes, feature), key=_order)


def listed_genes() -> set[str]:
    """The ACMG SF and carrier Tier 3 genes and the screen genes (192)."""
    lists = load_gene_lists()
    return (
        set(lists["acmg_sf"].genes)
        | set(lists["carrier_tier3"].genes)
        | set(load_screen_genes())
    )


def cut_spans(gtf: Path, genes: Collection[str]) -> dict[str, tuple[str, int, int]]:
    """Each listed gene's MANE Select `transcript` row on a `PRIMARY` chromosome: gene ->
    (chrom, start, end), 0-based half-open."""
    spans: dict[str, tuple[str, int, int]] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            if fields[0] not in PRIMARY or "MANE_Select" not in _TAG.findall(fields[8]):
                continue
            name = _GENE_NAME.search(fields[8])
            if name and name[1] in genes:
                spans[name[1]] = (fields[0], int(fields[3]) - 1, int(fields[4]))
    return spans


def _read_spans(text: str) -> dict[str, tuple[str, int, int]]:
    spans: dict[str, tuple[str, int, int]] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if tuple(fields) == SPAN_COLUMNS:
            continue
        gene, chrom, start, end = fields
        spans[gene] = (chrom, int(start), int(end))
    return spans


def load_spans() -> dict[str, tuple[str, int, int]]:
    """The packaged spans: gene (or `GBA1_ex9-11`) -> (chrom, start, end), 0-based."""
    path = files("genokeeper").joinpath("resources", "paralogs", SPANS)
    return _read_spans(path.read_text())


def locus_span(
    spans: Mapping[str, tuple[str, int, int]], genes: Sequence[str]
) -> tuple[str, int, int]:
    """The union of the genes' spans; refuse genes on more than one chromosome."""
    chroms = {spans[g][0] for g in genes}
    if len(chroms) != 1:
        raise ValueError(f"{', '.join(genes)}: spans on chromosomes {sorted(chroms)}")
    return (
        chroms.pop(),
        min(spans[g][1] for g in genes),
        max(spans[g][2] for g in genes),
    )


def load_exons() -> dict[str, list[Exon]]:
    """The packaged resource: gene -> its exons, in the file's order."""
    return read_resource(
        files("genokeeper").joinpath("resources", "deletions", RESOURCE).read_text()
    )


def read_resource(text: str) -> dict[str, list[Exon]]:
    """A resource in `COLUMNS` (this module's exons, slice 3c's CDS): gene -> its rows, in
    the file's order."""
    exons: dict[str, list[Exon]] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "gene":
            continue
        gene, transcript, number, chrom, start, end, tag = fields
        exons.setdefault(gene, []).append(
            Exon(
                gene,
                transcript,
                int(number),
                chrom,
                int(start),
                int(end),
                select=tag == "MANE_Select",
            )
        )
    return exons


def targets(
    exons: Mapping[str, Sequence[Exon]], margin: int
) -> dict[str, list[Target]]:
    """Per gene, its exons widened by `margin` on each side and merged where they overlap
    or touch, keeping the exon numbers merged into each: the MANE Select exons' numbers, or
    the Plus Clinical ones' where the target holds no Select exon, so one exon the two
    transcripts number differently reads as one."""
    result: dict[str, list[Target]] = {}
    for gene, gene_exons in exons.items():
        merged: list[tuple[str, int, int, list[Exon]]] = []
        for e in sorted(gene_exons, key=lambda x: (x.chrom, x.start, x.end)):
            start, end = max(e.start - margin, 0), e.end + margin
            if merged and merged[-1][0] == e.chrom and start <= merged[-1][2]:
                chrom, first, last, held = merged[-1]
                merged[-1] = (chrom, first, max(last, end), [*held, e])
            else:
                merged.append((e.chrom, start, end, [e]))
        result[gene] = [_target(c, s, e, held) for c, s, e, held in merged]
    return result


def _target(chrom: str, start: int, end: int, held: Sequence[Exon]) -> Target:
    """The target numbered by its MANE Select exons, else by all of them (Plus Clinical's)."""
    select = [e for e in held if e.select]
    numbers = tuple(sorted({e.number for e in select or held}))
    return Target(chrom, start, end, numbers, plus_clinical=not select)


def _pinned(gtf: Path) -> str:
    """The GTF's sha256; refuse any file but the pinned release."""
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    if digest != MANE_SHA256:
        raise ValueError(
            f"{gtf} has sha256 {digest}, not {MANE_RELEASE}'s {MANE_SHA256}"
        )
    return digest


def _header(digest: str, cut: str) -> list[str]:
    return [
        f"# source: {MANE_URL}",
        f"# release: {MANE_RELEASE}",
        f"# sha256: {digest}",
        f"# cut: {cut}",
    ]


def write_resource(gtf: Path, out: Path) -> int:
    """Cut the listed genes' MANE exons from `gtf` into `out`; refuse any file but the
    pinned release. Returns the number of exons written."""
    digest = _pinned(gtf)
    rows = cut_exons(gtf, listed_genes())
    lines = [
        *_header(
            digest,
            "every exon row of a MANE_Select or MANE_Plus_Clinical transcript whose "
            "gene_name is in the ACMG SF v3.3 or ACMG carrier Tier 3 list, "
            f"{SCREEN_GENES_NOTE}; start is 0-based half-open (BED), number is the "
            "transcript's own exon number",
        ),
        "\t".join(COLUMNS),
    ]
    lines += [
        f"{e.gene}\t{e.transcript}\t{e.number}\t{e.chrom}\t{e.start}\t{e.end}\t{tag}"
        for e, tag in rows
    ]
    out.write_text("\n".join(lines) + "\n")
    return len(rows)


def write_spans(gtf: Path, out: Path) -> int:
    """Write the MANE Select spans of `PARALOG_SPAN_GENES` and `GBA1_HOMOLOGY`'s exons to
    `out`; refuse any file but the pinned release or a gene without a span. Returns the
    number of rows written."""
    digest = _pinned(gtf)
    spans = cut_spans(gtf, PARALOG_SPAN_GENES)
    missing = sorted(set(PARALOG_SPAN_GENES) - set(spans))
    if missing:
        raise ValueError(f"{gtf}: no MANE Select span for {', '.join(missing)}")
    name, gene, numbers = GBA1_HOMOLOGY
    exons = [
        e
        for e, tag in cut_exons(gtf, {gene})
        if tag == "MANE_Select" and e.number in numbers
    ]
    if {e.number for e in exons} != set(numbers):
        raise ValueError(f"{gtf}: {gene} lacks MANE Select exons {numbers}")
    spans[name] = (
        exons[0].chrom,
        min(e.start for e in exons),
        max(e.end for e in exons),
    )
    lines = [
        *_header(
            digest,
            "each gene's MANE_Select transcript row on chr1-22, chrX or chrY; "
            f"{name} is the union of {gene}'s MANE_Select exons "
            f"{numbers[0]}-{numbers[-1]}; start is 0-based half-open (BED)",
        ),
        "\t".join(SPAN_COLUMNS),
    ]
    for g in (*PARALOG_SPAN_GENES, name):
        chrom, start, end = spans[g]
        lines.append(f"{g}\t{chrom}\t{start}\t{end}")
    out.write_text("\n".join(lines) + "\n")
    return len(spans)


def main(argv: Sequence[str] | None = None) -> int:
    """`cut GTF` writes the exon resource from MANE's GTF, `spans GTF` the spans one."""
    args = list(sys.argv[1:] if argv is None else argv)
    resources = Path(__file__).resolve().parent.parent / "resources"
    if len(args) == 2 and args[0] == "cut":
        out = resources / "deletions" / RESOURCE
        out.parent.mkdir(parents=True, exist_ok=True)
        print(f"{write_resource(Path(args[1]), out)} exons -> {out}")
        return 0
    if len(args) == 2 and args[0] == "spans":
        out = resources / "paralogs" / SPANS
        print(f"{write_spans(Path(args[1]), out)} spans -> {out}")
        return 0
    print("usage: python -m genokeeper.deletions.exons cut|spans GTF", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
