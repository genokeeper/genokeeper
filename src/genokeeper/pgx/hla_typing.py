"""HLA alleles beyond PharmCAT's, and the coeliac DQ types (slice 3f, "HLA").

`load_hla_alleles` reads `resources/hla/alleles.tsv`: the class I alleles PharmGKB grades at
2A or 2B that PharmCAT does not read (`hla.PGX_HLA` is refused there), HLA-B*27 as risk
information, and the research associations (HLA-C*06:02 with psoriasis, HLA-B*51 with
Behçet's disease, HLA-DQB1*06:02 with narcolepsy).

The coeliac rule-out asks two readings for the same answer (closing design, "Coeliac"):

- T1K's DQA1 and DQB1 (`dq_from_t1k`): DQ2.5 is DQA1*05 with DQB1*02:01, DQ2.2 DQA1*02:01
  with DQB1*02:02, DQ8 DQB1*03:02; the half-heterodimer is DQA1*05 or DQB1*02 (ESPGHAN
  2020). Unphased, a type's copies are the fewer of its two halves'. Each gene needs two kept
  alleles, and an allele whose equivalent names fall in different classes of that scheme
  (`DQB1*03:02,*03:01`) gives no status: a risk allele may hide in it.
- Monsuur et al. 2008's six tag SNPs (`dq_from_tags`, `resources/hla/dq_tags.tsv`, the tag
  allele the ALT in each): DQ2.5, DQ8 and DQ7 copies are their tag's dosage, DQ7 (DQ7.5,
  carrying DQA1*05) the half-heterodimer; DQ2.2 is the T,C,A haplotype at rs2395182,
  rs7775228 and rs4713586 (T,C,G is DQ4). The dosages are unphased, so DQ2.2's copies are
  counted over every phasing the dosages allow, and are no status unless all agree.

`coeliac`: "practically excluded" only when both see no DQ2, DQ8 or half-heterodimer; a
finding only when both see one and agree on the full types carried; anything else, a type
disagreement included, unassessed.
"""

import csv
import itertools
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from genokeeper.genotype.gvcf import Allele
from genokeeper.pgx.hla import PGX_HLA, HlaCall
from genokeeper.store.db import Kind, Topic, Validity

RESOURCES = Path(__file__).resolve().parents[1] / "resources" / "hla"
ALLELES_FILE = RESOURCES / "alleles.tsv"
DQ_TAGS_FILE = RESOURCES / "dq_tags.tsv"
HG002_TRUTH_FILE = RESOURCES / "hg002_truth.tsv"

# Monsuur et al. 2008's tags, by what each tells.
DQ25_TAG = "rs2187668"
DQ22_FIRST_TAG = "rs2395182"
DQ22_SECOND_TAG = "rs7775228"
DQ4_TAG = "rs4713586"
DQ8_TAG = "rs7454108"
DQ7_TAG = "rs4639334"

DRUG_NOTE = "screening result, needs confirmation by HLA typing in a laboratory"
B27_NOTE = (
    "only meaningful with symptoms (inflammatory back pain, uveitis); about 8% of"
    " Europeans carry it"
)
ASSOCIATION_NOTE = "an association: most carriers never develop the disease"
EXCLUDED = "coeliac disease practically excluded: no HLA-DQ2 or -DQ8"
TAGS_EUROPEAN = "the tag SNPs are validated in Europeans"
CONFIRM_DQ = "HLA-DQ typing in a laboratory"
COMMON = ": common (about 30–40% of Europeans), not a diagnosis"
HALF_ONLY = (
    "possibly carries half of HLA-DQ2 (DQA1*05 or DQB1*02 alone): coeliac disease is rare"
    " with it alone; not a diagnosis"
)
T1K_DQ_UNTYPED = (
    "T1K did not type DQA1 and DQB1 at two fields (an allele dropped, or equivalent"
    " alleles spanning a risk allele)"
)
TAGS_UNTYPED = (
    "a DQ tag SNP not called, or the unphased genotypes leave DQ2.2's copies open"
)
DQ_DISAGREE = "T1K's class II types and the DQ tag SNPs disagree"
# Each topic's reason for a carried allele; B*27 is the one risks row.
_NOTES = {
    Topic.DRUGS: DRUG_NOTE,
    Topic.RISKS: B27_NOTE,
    Topic.RESEARCH: ASSOCIATION_NOTE,
}


@dataclass(frozen=True)
class HlaAllele:
    """One row of `alleles.tsv`: the allele (two fields, or a one-field group such as
    `*27`), where a carrier's statement goes, its condition and drugs, its source level and
    citation, and the reason a carrier's statement gives (`note`, from the topic)."""

    gene: str
    allele: str
    topic: Topic
    validity: Validity
    condition: str
    drugs: tuple[str, ...]
    level: str
    note: str
    source: str


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


@cache
def load_hla_alleles() -> tuple[HlaAllele, ...]:
    """`alleles.tsv`'s rows; one PharmCAT reads itself (`PGX_HLA`) is refused."""
    found: list[HlaAllele] = []
    for row in _rows(ALLELES_FILE):
        gene, allele, topic = row["gene"], row["allele"], Topic(row["topic"])
        if allele in PGX_HLA.get(gene, ()):
            raise ValueError(f"{ALLELES_FILE}: {gene}{allele} is PharmCAT's")
        if topic not in _NOTES:
            raise ValueError(f"{ALLELES_FILE}: {gene}{allele}: topic {topic}")
        found.append(
            HlaAllele(
                gene=gene,
                allele=allele,
                topic=topic,
                validity=Validity(row["validity"]),
                condition=row["condition"],
                drugs=tuple(d for d in row["drugs"].split(";") if d),
                level=row["level"],
                note=_NOTES[topic],
                source=row["source"],
            )
        )
    return tuple(found)


def load_hg002_truth() -> dict[str, HlaCall]:
    """HG002's published types (`hg002_truth.tsv`), maternal first; names joined by `,` are
    each a match."""
    return {r["gene"]: (r["maternal"], r["paternal"]) for r in _rows(HG002_TRUTH_FILE)}


def _load_tags() -> dict[str, Allele]:
    tags = {
        r["rsid"]: (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
        for r in _rows(DQ_TAGS_FILE)
    }
    expected = {DQ25_TAG, DQ22_FIRST_TAG, DQ22_SECOND_TAG, DQ4_TAG, DQ8_TAG, DQ7_TAG}
    if set(tags) != expected:
        raise ValueError(f"{DQ_TAGS_FILE}: need exactly {sorted(expected)}")
    return tags


_TAGS = _load_tags()
TAG_RSIDS: tuple[str, ...] = tuple(_TAGS)
TAG_ALLELES: tuple[Allele, ...] = tuple(_TAGS.values())


@dataclass(frozen=True)
class DqStatus:
    """Copies of DQ2.5, DQ2.2 and DQ8, and whether half of DQ2 is carried (DQA1*05 or
    DQB1*02 from T1K; the DQ7 tag from the SNPs)."""

    dq25: int
    dq22: int
    dq8: int
    half: bool

    @property
    def full(self) -> frozenset[str]:
        """The full types carried."""
        copies = {"DQ2.5": self.dq25, "DQ2.2": self.dq22, "DQ8": self.dq8}
        return frozenset(name for name, n in copies.items() if n)

    @property
    def carries(self) -> bool:
        """Any DQ2, DQ8 or half-heterodimer."""
        return bool(self.full) or self.half


def _dqa1(name: str) -> str:
    if name.startswith("*05:"):
        return "05"
    return "02:01" if name == "*02:01" else "other"


def _dqb1(name: str) -> str:
    if name in ("*02:01", "*02:02", "*03:02"):
        return name[1:]
    return "02" if name.startswith("*02:") else "other"


def _classes(call: HlaCall | None, of: Callable[[str], str]) -> list[str] | None:
    """Each allele's class under `of` (`_dqa1`, `_dqb1`); None when an allele is missing
    or its equivalent names fall in different classes."""
    if call is None or None in call:
        return None
    found: list[str] = []
    for allele in call:
        names = {of(n) for n in (allele or "").split(",")}
        if len(names) != 1:
            return None
        found.append(names.pop())
    return found


def dq_from_t1k(t1k: Mapping[str, HlaCall]) -> DqStatus | None:
    """The DQ status of T1K's DQA1 and DQB1 (see the module docstring); None without two
    kept alleles in each, or when an allele's equivalents span classes."""
    a = _classes(t1k.get("HLA-DQA1"), _dqa1)
    b = _classes(t1k.get("HLA-DQB1"), _dqb1)
    if a is None or b is None:
        return None
    return DqStatus(
        dq25=min(a.count("05"), b.count("02:01")),
        dq22=min(a.count("02:01"), b.count("02:02")),
        dq8=b.count("03:02"),
        half="05" in a or any(x.startswith("02") for x in b),
    )


def _phasings(dosage: int) -> list[tuple[int, int]]:
    return [(0, 0)] if dosage == 0 else [(1, 1)] if dosage == 2 else [(1, 0), (0, 1)]


def _dq22_copies(first: int, second: int, dq4: int) -> int | None:
    """Chromosomes with the T,C,A haplotype, when every phasing of the three unphased
    dosages gives the same number; else None."""
    counts = {
        sum(a[h] == 1 and b[h] == 1 and c[h] == 0 for h in (0, 1))
        for a, b, c in itertools.product(
            _phasings(first), _phasings(second), _phasings(dq4)
        )
    }
    return counts.pop() if len(counts) == 1 else None


def dq_from_tags(dosages: Mapping[str, str]) -> DqStatus | None:
    """The DQ status of the six tag SNPs' dosages (rsid -> `0`, `1`, `2` or `unknown`) by
    Monsuur's rules; None when any is not known or DQ2.2's copies are left open."""
    d: dict[str, int] = {}
    for rsid in TAG_RSIDS:
        value = dosages.get(rsid, "unknown")
        if value not in ("0", "1", "2"):
            return None
        d[rsid] = int(value)
    dq22 = _dq22_copies(d[DQ22_FIRST_TAG], d[DQ22_SECOND_TAG], d[DQ4_TAG])
    if dq22 is None:
        return None
    return DqStatus(dq25=d[DQ25_TAG], dq22=dq22, dq8=d[DQ8_TAG], half=d[DQ7_TAG] > 0)


def carried_reason(t1k: DqStatus) -> str:
    """A carrier's wording from T1K's status: the full types, else half of DQ2."""
    dq2, dq8 = bool(t1k.dq25 or t1k.dq22), bool(t1k.dq8)
    if dq2 and dq8:
        return f"possibly carries HLA-DQ2 and -DQ8{COMMON}"
    if dq2 or dq8:
        return f"possibly carries HLA-DQ{'2' if dq2 else '8'}{COMMON}"
    if t1k.half:
        return HALF_ONLY
    raise ValueError("carries no DQ2, DQ8 or half-heterodimer")


def coeliac(t1k: DqStatus | None, tags: DqStatus | None) -> tuple[Kind, str]:
    """The coeliac statement's kind and reason from the two readings."""
    if t1k is None:
        return Kind.UNASSESSED, T1K_DQ_UNTYPED
    if tags is None:
        return Kind.UNASSESSED, TAGS_UNTYPED
    if not t1k.carries and not tags.carries:
        return Kind.NON_FINDING, f"{EXCLUDED}; {TAGS_EUROPEAN}"
    if t1k.carries and tags.carries and t1k.full == tags.full:
        return Kind.FINDING, f"{carried_reason(t1k)}; {TAGS_EUROPEAN}"
    return Kind.UNASSESSED, DQ_DISAGREE
