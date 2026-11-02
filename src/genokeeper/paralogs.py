"""The paralog loci's callers (slice 3d): Parascopy 1.19.0's loci, Gauchian 1.0.2 for GBA1,
KILDA 1.2.3 for Lp(a)'s KIV-2 repeat and Aldy 4.8.3 for CYP2A6; the wrappers and readers the
paralog sources grade (`sources.paralogs`).

Each tool runs offline in its image (`containers.IMAGES`) behind the core gate, on the CRAM,
its index and the FASTA mounted read-only and one output directory read-write. Parascopy's
loci run through 3a's own steps (`smn.parascopy_depth` once per CRAM, then
`smn.parascopy_cn` with every model of `LOCI` in one `cn-using`, behind `dna:paralogs`), so
this module only reads them: a locus is the rows of the paralog table overlapping the union
of its genes' MANE spans (`exons.locus_span`), each row's paralog-specific copy number
(psCN), or the aggregate over the gene and its paralog (agCN) where only that is called
(AMY1, the C4 total, GBA1). STRC's own depth is read natively (samtools `bedcov`) over the
part of CATSPER2 no STRC duplication covers (`unique_part`, `depth_ratio`).
"""

import json
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from genokeeper.cnv.screen import read_model
from genokeeper.containers import run_container
from genokeeper.deletions.exons import GBA1_HOMOLOGY
from genokeeper.pgx.cyp2d6 import _ALLELE, run_aldy
from genokeeper.pgx.outside import CallerCall
from genokeeper.regions import Interval, subtract
from genokeeper.smn import SAMPLE, _inputs, known, paralog_rows
from genokeeper.tools import run_tool

GAUCHIAN_VERSION = "1.0.2"
KILDA_VERSION = "1.2.3"
PARALOGS_LABEL = "dna:paralogs"
GAUCHIAN_LABEL = "dna:gauchian"
KILDA_LABEL = "dna:kilda"
# bp each side of LPA's MANE span whose reads KILDA counts (Task 6's commands).
KILDA_FLANK = 5000
# Inside KILDA's image (docker/kilda.Dockerfile): its script and its GRCh38 k-mer lists.
_KILDA = "/opt/kilda/bin/kilda.py"
_KMERS = "/opt/kilda/data/kmers_GRCh38_LPA_k31"
_GAUCHIAN_PREFIX = "gauchian"
# Parascopy's paralog-table columns (`smn.paralog_rows`) per measure: filter, value, quality.
_COLUMNS = {"pscn": (4, 5, 6), "agcn": (7, 8, 9)}


@dataclass(frozen=True)
class Locus:
    """A paralog locus: its name, the rows of `load_spans` whose spans `locus_span` joins
    (genes' MANE spans, or a part of one such as `GBA1_ex9-11`), the Parascopy model (its file stem), its role (`carrier` or `information`), the copy number
    read (`pscn` or `agcn`) and the unique neighbour whose depth confirms a loss (STRC's
    CATSPER2), if any."""

    locus: str
    genes: tuple[str, ...]
    model: str
    role: str
    measure: str
    unique: str | None = None


# The six loci (2026-10-30 pre-flight): UGT2B17 is not covered (Parascopy's UGT2B15 model
# covers 21% of its MANE span), nor GSTM1, NCF1 or SORD (A2, A3).
LOCI: dict[str, Locus] = {
    locus.locus: locus
    for locus in (
        Locus("STRC", ("STRC",), "STRC", "carrier", "pscn", "CATSPER2"),
        Locus("RHD", ("RHD",), "RHCE", "information", "pscn"),
        Locus("AMY1", ("AMY1A", "AMY1B", "AMY1C"), "AMY1C", "information", "agcn"),
        # The total C4 copy number (R2): C4A and C4B cannot be told apart reliably.
        Locus("C4", ("C4A", "C4B"), "C4A", "information", "agcn"),
        Locus("CYP2A6", ("CYP2A6",), "CYP2A7", "information", "pscn"),
        # GBA1's GBAP1-homologous exons 9-11, the part Gauchian's CN(GBA+GBAP1) counts: the
        # rest of GBA1 has no GBAP1 copy, so its rows carry agCN 2 (HG002: 4 and 2 over the
        # whole span).
        Locus("GBA1", (GBA1_HOMOLOGY[0],), "GBA", "carrier", "agcn"),
    )
}


@dataclass(frozen=True)
class PsCall:
    """A locus's copy number from Parascopy's rows over its span: the one value they all
    give (None: no row, a row unknown, rows that disagree, or a row not PASS), the lowest
    quality among them (None: no row or one unknown), whether every row passed its filter,
    and how many rows there were."""

    cn: int | None
    qual: float | None
    passed: bool
    rows: int


def read_locus(bed: Path, locus: Locus, span: tuple[str, int, int]) -> PsCall:
    """`locus`'s call in Parascopy's paralog table `bed`: its rows overlapping `span`
    (0-based, `locus_span(load_spans(), locus.genes)`), read at `locus.measure`."""
    chrom, start, end = span
    filt, value, quality = _COLUMNS[locus.measure]
    rows = [
        r
        for r in paralog_rows(bed)
        if r[0] == chrom and int(r[1]) < end and int(r[2]) > start
    ]
    passed = bool(rows) and all(r[filt] == "PASS" for r in rows)
    values = {known(r[value]) for r in rows}
    quals = [known(r[quality]) for r in rows]
    one = values.pop() if len(values) == 1 else None
    return PsCall(
        cn=int(one) if one is not None and passed else None,
        qual=None if not rows or None in quals else min(float(q or 0) for q in quals),
        passed=passed,
        rows=len(rows),
    )


def unique_part(span: tuple[str, int, int], model: Path) -> list[Interval]:
    """`span` minus both copies of each of the Parascopy `model`'s duplications
    (`cnv.screen.read_model`): the bases there only one gene holds."""
    _, copies = read_model(model)
    return subtract([span], copies)


def depth_ratio(
    cram: Path, fasta: Path, parts: Sequence[Interval], mean_depth: float
) -> float | None:
    """The mean depth of `cram` over `parts` (samtools `bedcov` natively, its default
    filters: no unmapped, secondary, QC-failed or duplicate read) over `mean_depth`, the
    run's; None when `parts` hold no base."""
    length = sum(end - start for _, start, end in parts)
    if length == 0:
        return None
    with tempfile.TemporaryDirectory(prefix=".bedcov-") as folder:
        bed = Path(folder) / "parts.bed"
        bed.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in parts))
        out = run_tool(
            ["samtools", "bedcov", "--reference", str(fasta), str(bed), str(cram)]
        )
    bases = sum(int(line.split("\t")[3]) for line in out.splitlines() if line)
    return bases / length / mean_depth


# --- Gauchian -------------------------------------------------------------------------------


@dataclass(frozen=True)
class GauchianCall:
    """Gauchian's call of GBA1: the total GBA+GBAP1 copy number, whether it calls a carrier
    or a biallelic GBAP1-like variant in exons 9-11, the variants it names there and the
    other unphased ones, the reads supporting the variant it called (the least
    supported, if several; None when it called none) and the deletion breakpoint it reports
    in GBA (None: `N/A`, none). None: a `None`, `?` or `N/A` cell."""

    total_cn: int | None
    carrier: bool | None
    biallelic: bool | None
    variants: tuple[str, ...]
    other: tuple[str, ...]
    reads: int | None
    breakpoint: str | None = None


def run_gauchian(cram: Path, fasta: Path, out_dir: Path, *, threads: int) -> Path:
    """Gauchian on `cram` (GRCh38) behind `dna:gauchian`, its manifest and outputs in
    `out_dir`; returns its TSV, `gauchian.tsv` (its JSON beside it, `gauchian.json`)."""
    (out_dir / "manifest.txt").write_text(f"/cram/{cram.name}\n")
    run_container(
        "gauchian",
        [
            "gauchian", "--manifest", "/out/manifest.txt", "--genome", "38",
            "--reference", f"/ref/{fasta.name}", "--prefix", _GAUCHIAN_PREFIX,
            "--outDir", "/out", "--threads", str(threads),
        ],
        mounts=_inputs(cram, fasta, out_dir),
        label=GAUCHIAN_LABEL,
    )  # fmt: skip
    return out_dir / f"{_GAUCHIAN_PREFIX}.tsv"


def _cell(value: str) -> str | None:
    return None if value.strip() in ("None", "?", "N/A", "") else value.strip()


def _flag(value: str) -> bool | None:
    cell = _cell(value)
    return None if cell is None else cell == "True"


def _names(value: str | None) -> tuple[str, ...]:
    cell = None if value is None else _cell(value)
    return (
        () if cell is None else tuple(v.strip() for v in cell.split(",") if v.strip())
    )


def read_gauchian(tsv: Path, json_path: Path) -> GauchianCall:
    """Gauchian's one sample from its TSV (the README's columns) and, for the reads behind
    the variant it called, its JSON: `variant_raw_count` gives each checked variant's
    `<alt>,<ref>` read counts (alt first: gauchian/depth_calling/gene_caller.py line 115 at
    v1.0.2), and `reads` is the least alt count among the variants its
    `GBAP1_like_variant_exon9_11` names (None when it names none or a count is missing)."""
    header, *rows = (line.split("\t") for line in tsv.read_text().splitlines() if line)
    if len(rows) != 1:
        raise ValueError(f"{tsv}: {len(rows)} Gauchian rows, not one")
    cell = dict(zip(header, rows[0], strict=True))
    total = _cell(cell["CN(GBA+GBAP1)"])
    (sample,) = json.loads(json_path.read_text()).values()
    called = _names(sample.get("GBAP1_like_variant_exon9_11"))
    counts = sample.get("variant_raw_count", {})
    alts = [int(counts[v].split(",")[0]) if v in counts else None for v in called]
    return GauchianCall(
        total_cn=None if total is None else int(total),
        carrier=_flag(cell["is_carrier(GBAP1-like_variant_exon9-11)"]),
        biallelic=_flag(cell["is_biallelic(GBAP1-like_variant_exon9-11)"]),
        variants=_names(cell["GBAP1-like_variant_exon9-11"]),
        other=_names(cell["other_unphased_variants"]),
        reads=None if not alts or None in alts else min(a or 0 for a in alts),
        breakpoint=_cell(cell["deletion_breakpoint_in_GBA"]),
    )


def wild_type(g: GauchianCall) -> bool:
    """Gauchian's wild type: four copies, no carrier, no biallelic call, no variant and no
    deletion breakpoint."""
    return (
        g.total_cn == 4
        and g.carrier is False
        and g.biallelic is False
        and not g.variants
        and not g.other
        and g.breakpoint is None
    )


# --- KILDA ----------------------------------------------------------------------------------


def run_kilda(
    cram: Path,
    fasta: Path,
    out_dir: Path,
    *,
    threads: int,
    span: tuple[str, int, int],
) -> Path:
    """KILDA's KIV-2 copy number on `cram` in one container behind `dna:kilda`: the reads
    over LPA's MANE `span` (0-based) +/- `KILDA_FLANK` and the unmapped ones, as pairs, their
    31-mers counted by jellyfish, then `kilda.py` against its GRCh38 k-mer lists; the
    temporaries removed. Returns `out_dir/kilda-output/kilda_kiv2.CNs`."""
    chrom, start, end = span
    region = f"{chrom}:{max(start - KILDA_FLANK, 0) + 1}-{end + KILDA_FLANK}"
    ref, reads = f"/ref/{fasta.name}", f"/cram/{cram.name}"
    (out_dir / "counts_list.tsv").write_text(f"{SAMPLE}\t/out/s.counts.tsv\n")
    temporaries = (
        "lpa.bam", "unmapped.bam", "both.bam", "r1.fq", "r2.fq", "s.jf", "s.counts.tsv",
    )  # fmt: skip
    pairs = " ".join(
        [
            "samtools collate -O -u /out/both.bam /out/tmp |",
            "samtools fastq -0 /dev/null",
        ]
        + ["-s /dev/null -1 /out/r1.fq -2 /out/r2.fq -"]
    )
    count = f"jellyfish count -m 31 -C -s 100M -t {threads} -o /out/s.jf"
    kilda = " ".join(
        [f"python {_KILDA} -c /out/counts_list.tsv"]
        + ["-k KIV2_hg38_kmers_6copies_specific.tsv"]
        + ["-l LPA_hg38_kmers_1copies_specific.tsv -o /out/kilda-output"]
    )
    script = " && ".join(
        [
            "set -o pipefail",
            "rm -rf /out/kilda-output",
            f"samtools view -b -T {ref} -o /out/lpa.bam {reads} {region}",
            f"samtools view -b -f 4 -T {ref} -o /out/unmapped.bam {reads} '*'",
            "samtools cat -o /out/both.bam /out/lpa.bam /out/unmapped.bam",
            pairs,
            f"{count} /out/r1.fq /out/r2.fq",
            "jellyfish dump -c -t /out/s.jf > /out/s.counts.tsv",
            f"cd {_KMERS}",
            kilda,
            "rm -f " + " ".join(f"/out/{name}" for name in temporaries),
        ]
    )
    run_container(
        "kilda",
        ["bash", "-c", script],
        mounts=_inputs(cram, fasta, out_dir),
        label=KILDA_LABEL,
    )
    return out_dir / "kilda-output" / "kilda_kiv2.CNs"


def read_kilda(path: Path) -> float | None:
    """The one sample's `KIV2_CN` in KILDA's table (columns `ID KIV2_CN quantile`); None
    for `NA`."""
    header, *rows = (line.split("\t") for line in path.read_text().splitlines() if line)
    if len(rows) != 1:
        raise ValueError(f"{path}: {len(rows)} KILDA rows, not one")
    value = dict(zip(header, rows[0], strict=True))["KIV2_CN"]
    return None if value in ("NA", "") else float(value)


# --- Aldy CYP2A6 ----------------------------------------------------------------------------

# CYP2A6*4: the whole-gene deletion.
_DELETION_STAR = "4"


def run_aldy_cyp2a6(cram: Path, fasta: Path, out_dir: Path) -> CallerCall:
    """Aldy on CYP2A6 (`cyp2d6.run_aldy`) behind `dna:paralogs`."""
    return run_aldy(cram, fasta, out_dir, "CYP2A6", label=PARALOGS_LABEL)


def aldy_copies(diplotype: str) -> int | None:
    """CYP2A6's copy number in an Aldy diplotype: per allele 0 for `*4` (the deletion), `N`
    for `xN`, else 1, summed over both haplotypes and each `+` part; None for alternatives
    or an allele it cannot read (CYP2D6's star-allele pattern, `cyp2d6._ALLELE`; a copy
    bin such as `x≥3` included)."""
    haplotypes = diplotype.split("/")
    if " or " in diplotype or len(haplotypes) != 2:
        return None
    total = 0
    for haplotype in haplotypes:
        for part in haplotype.split("+"):
            match = _ALLELE.match(part.strip())
            copies = "1" if match is None else match["copies"] or "1"
            if match is None or not copies.isdigit():
                return None
            copies = int(copies)
            total += 0 if match["star"] == _DELETION_STAR else copies
    return total
