"""Sample identity (slice 3c): the run's genotypes against the delivered 23andMe-format file.

The delivered file (`find_delivered`) is the one `*.txt` in data/raw, not `*.vc.txt`, whose
header carries 23andMe's column line: every platform position, hom-ref included (the
`.vc.txt` lists variant sites only and would bias the comparison). `compare` joins it by rsID
to the run's GWAS panel SNPs, located on GRCh38 with their alleles and called by the run;
autosomal biallelic SNVs only, a delivered genotype compared as an unordered pair on the
plus strand or, for a SNP that is not A/T or C/G, complemented (GRCh37 and GRCh38 plus
strands differ in a few inverted regions). `judge` applies `[identity]`.

The provider derived the delivered genotypes from the same reads, so the check catches a FASTQ
mix-up, not a lab's sample swap (`IDENTITY_TOOL`).
"""

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from genokeeper.config import IdentityThresholds

DELIVERED_HEADER = "# rsid\tchromosome\tposition\tgenotype"
IDENTITY_TOOL = "FASTQ mix-up check against the delivered genotypes"

NO_CALL = "no call"
PALINDROMIC = "A/T or C/G"
NOT_SNV = "not a biallelic SNV"
OTHER_ALLELE = "allele outside REF/ALT"
UNKNOWN = "our genotype unknown"
NOT_DELIVERED = "not in the delivered file"

_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_PALINDROMES = (frozenset("AT"), frozenset("CG"))
_AUTOSOMES = frozenset(str(n) for n in range(1, 23))


def _has_header(path: Path) -> bool:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("#"):
                return False  # the first data row ends the header
            if line.rstrip("\r\n") == DELIVERED_HEADER:
                return True
    return False


def find_delivered(raw: Path) -> Path | None:
    """The one `*.txt` in `raw`, not `*.vc.txt`, whose `#` header lines (up to the first
    data row) include `DELIVERED_HEADER`; any other `.txt` is ignored. None when there is
    none; `ValueError` for two or more (an ambiguous input the run refuses, naming
    `--genotypes`)."""
    if not raw.is_dir():
        return None
    found = [
        path
        for path in sorted(raw.glob("*.txt"))
        if path.is_file() and not path.name.endswith(".vc.txt") and _has_header(path)
    ]
    if len(found) > 1:
        names = ", ".join(path.name for path in found)
        raise ValueError(f"two or more delivered genotype files in {raw}: {names}")
    return found[0] if found else None


def read_delivered(path: Path) -> dict[str, str]:
    """rsID -> genotype (`"AG"`) from the 23andMe rows `rsid<TAB>chromosome<TAB>position<TAB>
    genotype`: `#` lines skipped, only `rs` ids on chromosomes 1-22 kept."""
    genotypes: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                continue
            rsid, chromosome, _, genotype = fields[:4]
            if rsid.startswith("rs") and chromosome in _AUTOSOMES:
                genotypes[rsid] = genotype.strip().upper()
    return genotypes


@dataclass(frozen=True)
class Identity:
    """The SNPs compared, those that agree, and the rest counted by skip reason."""

    compared: int
    concordant: int
    skipped: dict[str, int]

    @property
    def concordance(self) -> float:
        return self.concordant / self.compared if self.compared else 0.0


def _oriented(theirs: str, ref: str, alt: str) -> str:
    """The delivered genotype on the plus strand, or complemented where it does not fit."""
    return theirs if set(theirs) <= {ref, alt} else theirs.translate(_COMPLEMENT)


def _skip_reason(theirs: str, ref: str, alt: str, dosage: str) -> str | None:
    if len(ref) != 1 or len(alt) != 1 or not {ref, alt} <= _BASES or ref == alt:
        return NOT_SNV
    if frozenset((ref, alt)) in _PALINDROMES:
        return PALINDROMIC
    if not theirs:
        return NOT_DELIVERED
    if theirs == "--":
        return NO_CALL
    if len(theirs) != 2 or not set(theirs) <= _BASES:
        return NOT_SNV  # indel codes (`DI`, `II`) and anything else
    if dosage not in ("0", "1", "2"):
        return UNKNOWN
    if not set(_oriented(theirs, ref, alt)) <= {ref, alt}:
        return OTHER_ALLELE
    return None


def compare(
    delivered: Mapping[str, str], ours: Iterable[tuple[str, str, str, str]]
) -> Identity:
    """`ours` is (rsid, REF, ALT, dosage), one per rsID (the first kept); our genotype is
    REF x (2 - dosage) + ALT x dosage, compared with the delivered one as an unordered pair,
    on the plus strand or complemented."""
    seen: set[str] = set()
    skipped: Counter[str] = Counter()
    compared = concordant = 0
    for rsid, ref, alt, dosage in ours:
        if rsid in seen:
            continue
        seen.add(rsid)
        theirs = delivered.get(rsid, "")
        reason = _skip_reason(theirs, ref, alt, dosage)
        if reason is not None:
            skipped[reason] += 1
            continue
        copies = int(dosage)
        compared += 1
        mine = ref * (2 - copies) + alt * copies
        concordant += sorted(mine) == sorted(_oriented(theirs, ref, alt))
    return Identity(compared=compared, concordant=concordant, skipped=dict(skipped))


class Verdict(StrEnum):
    MATCHED = "matched"
    MISMATCH = "mismatch"
    NOT_CHECKED = "not checked"


def judge(identity: Identity, t: IdentityThresholds) -> tuple[Verdict, str]:
    """Below `t.min_compared` SNPs, not checked (the run goes on); else below
    `t.min_concordance`, a mismatch (the run refuses); else matched."""
    if identity.compared < t.min_compared:
        return Verdict.NOT_CHECKED, (
            "too few common SNPs to check the sample's identity:"
            f" {identity.compared:,} compared"
        )
    if identity.concordance < t.min_concordance:
        return Verdict.MISMATCH, (
            "the FASTQ does not match the delivered genotypes: concordance"
            f" {identity.concordance:.3f} at {identity.compared:,} SNPs"
        )
    return (
        Verdict.MATCHED,
        f"{identity.compared:,} SNPs, concordance {identity.concordance:.4f}",
    )
