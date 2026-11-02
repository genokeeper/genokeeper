"""Blood groups by RBCeq2 2.4.4 on our own calls (slice 3f, "Blood groups").

RBCeq2's input is a VCF of our calls at the positions of its ISBT allele table (`db.tsv`,
inside the pinned wheel; DB 2.5.1), written by PharmCAT's writer
(`pharmcat.write_pharmcat_vcf`, no second writer): one row per position, reference and
no-call rows explicit, a male's chrX outside the PARs one allele. The positions come from a
packaged cut of `db.tsv`'s `GRCh38` column (`POSITIONS_FILE`, its header naming db.tsv's
sha256; cut once, by hand, for RBCeq2 2.4.4, ruling R23): each `pos_REF_ALT` item an allele,
each `pos_ref` item a reference-only position whose REF is the FASTA's base and ALT `.`
(dropped where a variant item shares the position, whose row says the same). Structural
items (`…_del_53kb`, `…_INS_6431`), `unknown…`, `na` and empty items are skipped; an odd item
is read by its prefix (`44254581_ref_`, `25284544_G_C_no_phenotype`).

The cut leaves out three groups, never reported: the platelet antigens (R8: RBCeq2 2.4.4's
HPA rows hold GRCh37 positions in its GRCh38 column), RH (`RHD`, `RHCE`: copy number short
reads do not give) and MNS (`GYPA`, `GYPB` and their hybrids `GYP`: hybrid genes the truth
genomes do not carry). RBCeq2 runs without `--HPAs`, `--RH` or `--phased`.

Two rulings decide what a system's call is worth before RBCeq2 is read:
- R9, per allele (R14): `rbceq2_positions` checks each allele's REF against the GRCh38
  FASTA; a mismatching allele is not written (calls_at refuses a REF the FASTA does not
  have) and is named under its system (`unchecked`), whose statement says so; the system
  stays shown.
- R10: RBCeq2 reads an explicit `./.` as the reference and makes a confident call from it
  (`./.` at ABO c.261delG gave A1/A1), so a system with any no-call among its written
  positions is unassessed whatever RBCeq2 says (`uncovered`, read from the VCF as written).

`run_rbceq2` runs it offline in its image, behind the core gate as `dna:bloodgroups`;
`read_rbceq2` reads its genotype and alphanumeric phenotype tables (one row per sample, an
index header `UUID: <uuid>` that changes per run, so the sample's row is read by position).
"""

import csv
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import pysam

from genokeeper.containers import PLATFORMS, run_container
from genokeeper.genotype.gvcf import Allele

RBCEQ2_VERSION = "2.4.4"
DB_VERSION = "2.5.1"
SOURCE_VERSION = f"{RBCEQ2_VERSION} (DB {DB_VERSION})"
POSITIONS_FILE = (
    Path(__file__).resolve().parent
    / "resources"
    / "bloodgroups"
    / f"rbceq2_{RBCEQ2_VERSION}_positions.tsv"
)
LABEL = "dna:bloodgroups"
OUT_NAME = "rbceq2"
GENO_SUFFIX = "_geno.tsv"
PHENO_SUFFIX = "_pheno_alphanumeric.tsv"
# RBCeq2's genotype for a system it could not call.
UNDETERMINED = "Undetermined/Undetermined"
# MNS's genes, never reported (see the module docstring).
MNS_GENES = frozenset({"GYPA", "GYPB", "GYP"})


@dataclass(frozen=True)
class DbItem:
    """One position of `db.tsv`'s `GRCh38` column: its system (RBCeq2's name, `KLF1` read
    as `KLF` as RBCeq2 reads it), and REF and ALT (None for a reference-only item)."""

    system: str
    chrom: str
    pos: int
    ref: str | None
    alt: str | None


def _order(item: DbItem) -> tuple[int, int, str, str, str]:
    number = item.chrom.removeprefix("chr")
    rank = int(number) if number.isdigit() else 23 if number == "X" else 24
    return (rank, item.pos, item.ref or "", item.alt or "", item.system)


@cache
def load_positions(path: Path = POSITIONS_FILE) -> tuple[DbItem, ...]:
    """The packaged items."""
    with path.open(newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    return tuple(
        DbItem(
            row["system"],
            row["chrom"],
            int(row["pos"]),
            None if row["ref"] == "." else row["ref"],
            None if row["alt"] == "." else row["alt"],
        )
        for row in csv.DictReader(lines, delimiter="\t")
    )


def systems(items: Sequence[DbItem] | None = None) -> tuple[str, ...]:
    """The systems reported, sorted: those with a packaged position."""
    return tuple(sorted({i.system for i in (items or load_positions())}))


@dataclass(frozen=True)
class Rbceq2Positions:
    """What RBCeq2's input is written at: the alleles (`calls_at` and
    `write_pharmcat_vcf` take them; a reference-only position's ALT is `.`), each written
    position's systems, and each system's alleles left out because their REF disagrees
    with the FASTA (`allele_name`s)."""

    alleles: tuple[Allele, ...]
    systems: Mapping[tuple[str, int], frozenset[str]]
    unchecked: Mapping[str, tuple[str, ...]]


def allele_name(item: DbItem) -> str:
    """`chrom:pos:REF>ALT`, a REF or ALT longer than ten bases written as its length."""

    def short(bases: str | None) -> str:
        bases = bases or "."
        return bases if len(bases) <= 10 else f"{len(bases)}bp"

    return f"{item.chrom}:{item.pos}:{short(item.ref)}>{short(item.alt)}"


def rbceq2_positions(
    fasta: Path, items: Sequence[DbItem] | None = None
) -> Rbceq2Positions:
    """The packaged items (or `items`) against the GRCh38 `fasta`: a variant's REF must be
    the FASTA's bases there, else the allele is left out and named under its system
    (`unchecked`); a reference-only position takes the FASTA's base as REF and `.` as ALT, and is
    left out where a variant allele is written at the same position."""
    items = load_positions() if items is None else items
    alleles: dict[Allele, set[str]] = defaultdict(set)
    unchecked: dict[str, list[str]] = defaultdict(list)
    with pysam.FastaFile(str(fasta)) as reference:

        def bases(chrom: str, pos: int, length: int) -> str:
            return reference.fetch(chrom, pos - 1, pos - 1 + length).upper()

        for i in items:
            if i.ref is None or i.alt is None:
                continue
            if bases(i.chrom, i.pos, len(i.ref)) != i.ref:
                unchecked[i.system].append(allele_name(i))
                continue
            alleles[(i.chrom, i.pos, i.ref, i.alt)].add(i.system)
        variant_at = {(chrom, pos) for chrom, pos, _, _ in alleles}
        for i in items:
            if i.ref is not None or (i.chrom, i.pos) in variant_at:
                continue
            alleles[(i.chrom, i.pos, bases(i.chrom, i.pos, 1), ".")].add(i.system)
    at: dict[tuple[str, int], set[str]] = defaultdict(set)
    for (chrom, pos, _, _), names in alleles.items():
        at[(chrom, pos)] |= names
    order = sorted(alleles, key=lambda a: _order(DbItem("", a[0], a[1], a[2], a[3])))
    return Rbceq2Positions(
        alleles=tuple(order),
        systems={key: frozenset(names) for key, names in at.items()},
        unchecked={system: tuple(names) for system, names in unchecked.items()},
    )


def uncovered(vcf: Path, positions: Rbceq2Positions) -> dict[str, tuple[int, ...]]:
    """Each system with a no-call (a `.` in GT) at a position of RBCeq2's input `vcf` as
    written, with those positions (R10)."""
    found: dict[str, set[int]] = defaultdict(set)
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        chrom, pos, gt = fields[0], int(fields[1]), fields[9].split(":")[0]
        if "." in gt:
            for system in positions.systems.get((chrom, pos), ()):
                found[system].add(pos)
    return {system: tuple(sorted(found[system])) for system in sorted(found)}


@dataclass(frozen=True)
class SystemCall:
    """RBCeq2's call of one system: the genotype (None when it made none) and the
    alphanumeric phenotype (None when empty)."""

    genotype: str | None
    phenotype: str | None


def run_rbceq2(vcf: Path, out_dir: Path) -> Path:
    """RBCeq2 on `vcf` (written by `write_pharmcat_vcf` at `rbceq2_positions`), its tables
    written into `out_dir` as `OUT_NAME<suffix>`; `out_dir`. `RuntimeError` when it exits 0
    without them, and what `run_container` raises."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "rbceq2",
        [
            "rbceq2", "--vcf", f"/in/{vcf.name}", "--out", f"/out/{OUT_NAME}",
            "--reference_genome", "GRCh38", "--processes", "1",
        ],
        mounts=[(vcf.parent, "/in", True), (out_dir, "/out", False)],
        platform=PLATFORMS["rbceq2"],
        label=LABEL,
    )  # fmt: skip
    for suffix in (GENO_SUFFIX, PHENO_SUFFIX):
        if not (out_dir / f"{OUT_NAME}{suffix}").is_file():
            raise RuntimeError(f"RBCeq2 exited 0 without {OUT_NAME}{suffix}")
    return out_dir


def _one_sample(path: Path) -> dict[str, str]:
    """The table's one sample row, column -> value (the index column left out)."""
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows or not rows[0] or not rows[0][0].startswith("UUID:"):
        raise ValueError(f"{path}: not an RBCeq2 table")
    if len(rows) != 2:
        raise ValueError(f"{path}: need one sample, found {len(rows) - 1}")
    header, values = rows
    return dict(zip(header[1:], values[1:], strict=True))


def read_rbceq2(out_dir: Path) -> dict[str, SystemCall]:
    """Each system of RBCeq2's tables in `out_dir`: its genotype and phenotype."""
    geno = _one_sample(out_dir / f"{OUT_NAME}{GENO_SUFFIX}")
    pheno = _one_sample(out_dir / f"{OUT_NAME}{PHENO_SUFFIX}")
    return {
        system: SystemCall(
            genotype=None if value in ("", UNDETERMINED) else value,
            phenotype=pheno.get(system) or None,
        )
        for system, value in geno.items()
    }
