"""The genotype a VCF row we write for an outside tool carries (PharmCAT's input, slice 1b;
pgsc_calc's, slice 1d): one position, from our calls at each of its ALTs and their grades.

- each position's genotype comes from the dosage of each of its ALTs: `0/1`, `1/2`, `1/1`,
  and `0/0` only when no ALT has a copy and nothing else was called there;
- a position with any ALT's call unknown (or missing) is `./.`; so is one with any ALT's
  call graded low or not graded, so the tool reads it missing rather than resting a result
  on a weak call (the slice 1b ruling); and so is one where another allele the row does not
  list was called: the rest of the genotype is that allele, not the reference (builder's
  call). A call's `OTHER_ALLELE_WARNING` alone does not say so, because the call layer
  gives it to every ALT the record does not call, a sibling ALT of the position included;
  an unlisted allele is one the record's GT carries beyond the listed ALTs' copies. Its
  anchor comes from that GT, which the call layer copies from the gVCF record unchanged, so
  neither a sibling's warning nor the row's ALT list can move it;
- on a haploid locus (a male's chrX outside the PARs) the genotype is one allele (`0`, `1`,
  `.`); a heterozygote there is `.` (the call layer keeps it as written, an error for one
  copy).

A GWAS risk allele's copies (`risk_copies`, slice 1d; ruling R26) follow the same call: an
ALT's are that ALT's dosage (0 when another allele is called); a reference risk allele's
are the ploidy less that dosage (a reference block or RefCall: the ploidy), counted only
when every non-reference allele of the GT is the located ALT. Beside another allele (a SNP,
or an indel that keeps the reference base on its haplotype), or for a risk allele that is
neither the REF nor the ALT, they are `NOT_COUNTABLE`; unknown only when the call is.
"""

from collections.abc import Sequence
from typing import Final, Literal

from genokeeper.genotype.gvcf import OTHER_ALLELE_WARNING, Called
from genokeeper.store.db import Level

NOT_COUNTABLE: Final = "not_countable"


def _alt_copies(gt: str | None) -> int:
    """The ALT copies a record's GT calls, any ALT (`0/1` 1, `1/2` 2); unknown counts
    as more than any position's."""
    if gt is None:
        return 3
    return sum(a not in ("0", ".") for a in gt.replace("|", "/").split("/"))


def site_genotype(
    called: Sequence[Called | None], levels: Sequence[Level | None], *, haploid: bool
) -> str:
    """One position's GT from the calls of each of its ALTs and their levels, in ALT
    order."""
    nocall = "." if haploid else "./."
    if any(c is None or c.dosage == "unknown" for c in called):
        return nocall
    if any(level is None or level == Level.LOW for level in levels):
        return nocall
    known = [c for c in called if c is not None]
    if haploid and any(
        c.gt is not None and len(c.gt.replace("|", "/").split("/")) > 1 for c in known
    ):
        return nocall
    ploidy = 1 if haploid else 2
    copies = [int(c.dosage) for c in known]
    total = sum(copies)
    if total > ploidy:
        return nocall
    if total < ploidy and any(
        (OTHER_ALLELE_WARNING in c.warnings or c.gt is not None)
        and _alt_copies(c.gt) > total
        for c in known
    ):
        return nocall  # the rest of the genotype is an ALT this row does not list
    indices = [0] * (ploidy - total) + [
        i for i, n in enumerate(copies, start=1) for _ in range(n)
    ]
    return "/".join(str(i) for i in indices)


def risk_copies(
    vid: str, risk_allele: str, gt: str | None, dosage: str | None
) -> int | Literal["not_countable"] | None:
    """The copies of `risk_allele` (exactly the REF or the ALT of `vid`, with or without
    its `GRCh38:` prefix) the call at `vid` carries, from its GT and dosage as the call
    layer writes them (None for no call row); None when the genotype is unknown,
    `NOT_COUNTABLE` for any other allele or a reference risk allele beside another."""
    ref, alt = vid.split(":")[-2:]
    if risk_allele not in (ref, alt):
        return NOT_COUNTABLE
    if dosage is None or dosage == "unknown" or gt is None:
        return None
    if risk_allele == alt:
        return int(dosage)
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) - alleles.count("0") != int(dosage):
        return NOT_COUNTABLE
    return len(alleles) - int(dosage)
