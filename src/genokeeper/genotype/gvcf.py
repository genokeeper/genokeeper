"""The call layer: our own DeepVariant gVCF read at every allele a claim is about.

A claim's allele and the gVCF's alleles are compared after normalisation (`bcftools norm -f
<reference> -c e -m -both`, each allele on its own), so a homopolymer deletion DeepVariant
writes as `CAA>CA` beside `CAA>C` is the claim's `CA>C`. Records come from pysam's `fetch`,
which honours a block's `END`. For each claim (track design, stage 2a, "The call layer"):

- a variant record holding the claim's allele gives its genotype, GQ, DP, AD and VAF; the
  dosage counts that allele's copies (`1/2` with the claim's ALT 2 is 1);
- a RefCall at the position, or one whose REF spans it, gives `0/0` with the network's GQ;
- a reference block gives `0/0` with the block's GQ (binomial, capped at `BLOCK_GQ_CAP`) and
  its MIN_DP as depth, as its own row type (`RowType.REF_BLOCK`);
- a no-call (`./.`, a record's or a block's) stays unknown (`NOCALL_REASON`), and so does no
  coverage (a block with MIN_DP 0, which DeepVariant writes `0/0` with GQ 1, or no record at
  all; `NO_COVERAGE_REASON`);
- another allele called at the position, and not the claim's, gives dosage 0 with
  `OTHER_ALLELE_WARNING`: DeepVariant writes every candidate at a position into one record;
- a position inside a called deletion with no record of its own, or a called `*` allele where
  the claim's allele has no copy, is unknown (`INSIDE_DELETION_REASON`) (builder's call);
  so is a position inside another called multi-base record (`MULTI_BASE_REASON`).

On a haploid locus (`haploid(chrom, pos)`: a male's chrX or chrY outside the PARs) a
homozygous genotype is one copy: DeepVariant writes a hom call `1/1` and every reference
block and RefCall `0/0` there, which read as `1` and `0`. A heterozygote is kept as written;
`reliability.grade_called` makes it low.

Grading (levels, the indel caps, the depth warning) is
`reliability.grade_called`'s; this module reports what the gVCF says and nothing else.
"""

import tempfile
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.store.db import RowType
from genokeeper.tools import run_tool
from genokeeper.variants.key import ucsc_contig

type Allele = tuple[str, int, str, str]

OTHER_ALLELE_WARNING = "another allele called at this position"
INSIDE_DELETION_REASON = "genotype unknown: inside a called deletion"
INDEL_BLOCK_WARNING = "indel site in a reference block; reference inferred"
NOCALL_REASON = "genotype unknown: no-call"
NO_COVERAGE_REASON = "genotype unknown: no coverage"
MULTI_BASE_REASON = "genotype unknown: inside a multi-base record"

# A reference block's GQ is binomial and capped at 50 (DeepVariant's gVCF; track design 2a).
BLOCK_GQ_CAP = 50
# Normalisation only moves an allele left (or right within its own REF), so a record that
# normalises onto a claim starts at most this far right of it; DeepVariant writes its indels
# left-aligned, so this is a margin, not a repeat length.
NORM_WINDOW = 100

_SYMBOLIC = "<*>"
_STAR = "*"
_ROW_TYPES = {
    "PASS": RowType.PASS,
    "RefCall": RowType.REFCALL,
    "NoCall": RowType.NOCALL,
}


@dataclass(frozen=True)
class Called:
    """What the gVCF says of one claimed allele: its genotype as written (a haploid locus's
    homozygote as one index), the allele's dosage (`0`, `1`, `2` or `unknown`), the row it
    rests on, that row's GQ and depth (a block's MIN_DP), AD as written, the claim allele's
    VAF, and the warnings (an unknown's reason among them)."""

    gt: str | None
    dosage: str
    row_type: RowType
    gq: int | None
    dp: int | None
    ad: tuple[int, ...] | None
    vaf: float | None
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class _Record:
    """One gVCF row; `end` is its last reference position (a block's END), 1-based."""

    chrom: str
    pos: int
    end: int
    ref: str
    alts: tuple[str, ...]
    filt: str
    gt: tuple[int | None, ...]
    phased: bool
    gq: int | None
    dp: int | None
    min_dp: int | None
    ad: tuple[int, ...] | None
    vaf: tuple[float, ...] | None

    @property
    def block(self) -> bool:
        return self.alts == (_SYMBOLIC,)

    @property
    def nocall(self) -> bool:
        return not self.gt or None in self.gt or self.filt == "NoCall"

    def called(self) -> set[int]:
        """The ALT numbers its genotype calls (1-based)."""
        return {i for i in self.gt if i is not None and i > 0}

    def row_type(self) -> RowType:
        if self.block:
            return RowType.REF_BLOCK
        row_type = _ROW_TYPES.get(self.filt)
        if row_type is None:
            raise ValueError(
                f"unexpected FILTER {self.filt!r} at {self.chrom}:{self.pos}:"
                " expected PASS, RefCall or NoCall"
            )
        return row_type

    def gt_text(self, *, haploid: bool) -> str:
        indices = ["." if i is None else str(i) for i in self.gt]
        if haploid and "." not in indices and len(set(indices)) == 1:
            return indices[0]
        return ("|" if self.phased else "/").join(indices)


def _optional_int(value: object) -> int | None:
    return None if value is None else int(str(value))


def _tuple[T](value: object, kind: Callable[[str], T]) -> tuple[T, ...] | None:
    if not isinstance(value, tuple) or any(v is None for v in value):
        return None
    return tuple(kind(str(v)) for v in value)


def _record(rec: pysam.VariantRecord) -> _Record:
    sample = rec.samples[0]
    filters = [str(name) for name in rec.filter]
    return _Record(
        chrom=rec.chrom,
        pos=rec.pos,
        end=rec.stop,  # 0-based exclusive: the last base, 1-based
        ref=(rec.ref or "").upper(),
        alts=tuple(a.upper() if a != _SYMBOLIC else a for a in rec.alts or ()),
        filt=";".join(filters) if filters else ".",
        gt=tuple(sample.get("GT") or ()),
        phased=bool(sample.phased),
        gq=_optional_int(sample.get("GQ")),
        dp=_optional_int(sample.get("DP")),
        min_dp=_optional_int(sample.get("MIN_DP")),
        ad=_tuple(sample.get("AD"), int),
        vaf=_tuple(sample.get("VAF"), float),
    )


def normalise(
    alleles: list[Allele], *, fasta: Path, contigs: dict[str, int], scratch: Path
) -> list[Allele]:
    """Each allele left-aligned and trimmed against `fasta` by `bcftools norm`, in order.

    `-c e` stops at a REF the reference does not have: the alleles and the gVCF are on one
    reference, so a mismatch means the wrong FASTA, never a claim to skip.

    The comparison's anchor comes from the lab's GRCh38 FASTA, which a record written in
    another frame (`CAA>CA` for `CA>C`) cannot move, because both sides are rewritten
    against that one reference before they are compared.

    The alleles (the user's calls) are written under `scratch`, the run's work directory
    inside reports/ (deny-ruled, git-ignored), never the system temp, and removed after."""
    if not alleles:
        return []
    order = {name: i for i, name in enumerate(contigs)}
    lines = sorted(
        (order[chrom], pos, f"{chrom}\t{pos}\t{i}\t{ref}\t{alt}\t.\t.\t.\n")
        for i, (chrom, pos, ref, alt) in enumerate(alleles)
    )
    header = "##fileformat=VCFv4.2\n" + "".join(
        f"##contig=<ID={name},length={length}>\n" for name, length in contigs.items()
    )
    with tempfile.TemporaryDirectory(prefix=".norm-", dir=scratch) as work:
        raw = Path(work) / "alleles.vcf"
        raw.write_text(
            header
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join(line for _, _, line in lines)
        )
        out = run_tool(
            [
                "bcftools", "norm", "--no-version", "-f", str(fasta), "-c", "e",
                "-m", "-both", "-Ov", str(raw),
            ]
        )  # fmt: skip
    normalised: dict[int, Allele] = {}
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        chrom, pos, number, ref, alt = line.split("\t")[:5]
        normalised[int(number)] = (chrom, int(pos), ref.upper(), alt.upper())
    missing = [alleles[i] for i in range(len(alleles)) if i not in normalised]
    if missing:
        raise RuntimeError(f"bcftools norm returned no record for {missing[:3]}")
    return [normalised[i] for i in range(len(alleles))]


def _unknown(
    reason: str,
    row_type: RowType = RowType.NONE,
    *,
    gt: str | None = None,
    gq: int | None = None,
    dp: int | None = None,
    ad: tuple[int, ...] | None = None,
) -> Called:
    return Called(
        gt=gt, dosage="unknown", row_type=row_type, gq=gq, dp=dp, ad=ad, vaf=None,
        warnings=(reason,),
    )  # fmt: skip


def _from_record(rec: _Record, allele: int | None, *, haploid: bool) -> Called:
    """The claim's call from the record at its position; `allele` is the claim's ALT number
    in it, or None when the record does not hold it."""
    row_type = rec.row_type()
    gt = rec.gt_text(haploid=haploid)
    if rec.nocall:
        return _unknown(
            NOCALL_REASON, RowType.NOCALL, gt=gt, gq=rec.gq, dp=rec.dp, ad=rec.ad
        )
    called = rec.called()
    copies = sum(i == allele for i in rec.gt) if allele is not None else 0
    if haploid and len(set(rec.gt)) == 1:
        copies = min(copies, 1)
    vaf = rec.vaf[allele - 1] if allele is not None and rec.vaf is not None else None
    others = called - {allele}
    if copies == 0 and any(rec.alts[i - 1] == _STAR for i in others):
        return _unknown(
            INSIDE_DELETION_REASON, row_type, gt=gt, gq=rec.gq, dp=rec.dp, ad=rec.ad
        )
    warnings = (OTHER_ALLELE_WARNING,) if copies == 0 and others else ()
    return Called(
        gt=gt,
        dosage=str(copies),
        row_type=row_type,
        gq=rec.gq,
        dp=rec.dp,
        ad=rec.ad,
        vaf=vaf,
        warnings=warnings,
    )


def _deletes(rec: _Record) -> bool:
    """Whether a called allele of `rec` may remove a base its REF spans: a `*`, or any
    allele shorter than REF (read conservatively: which bases it removes is not worked
    out)."""
    alleles = [rec.alts[i - 1] for i in rec.called()]
    return any(a == _STAR or len(a) < len(rec.ref) for a in alleles)


def _call(
    claim: Allele,
    records: list[_Record],
    normal: dict[tuple[_Record, int], Allele],
    *,
    haploid: bool,
) -> Called:
    _, pos, ref, _ = claim
    variants = [r for r in records if not r.block]
    at = [
        r
        for r in variants
        if r.pos == pos
        or any(
            normal[(r, i)][1] == pos
            for i in range(1, len(r.alts) + 1)
            if (r, i) in normal
        )
    ]
    for rec in at:
        for i in range(1, len(rec.alts) + 1):
            if normal.get((rec, i)) == claim:
                return _from_record(rec, i, haploid=haploid)
    if at:
        # The claim's allele is in no record at its position: the first that calls
        # something (or fails to call) is the verdict on the locus.
        deciding = next((r for r in at if r.nocall or r.called()), at[0])
        return _from_record(deciding, None, haploid=haploid)
    spanning = [r for r in variants if r.pos < pos <= r.end]
    if any(_deletes(r) for r in spanning if not r.nocall):
        return _unknown(INSIDE_DELETION_REASON)
    if spanning:
        rec = spanning[0]
        if rec.nocall or rec.called():
            reason = NOCALL_REASON if rec.nocall else MULTI_BASE_REASON
            return _unknown(reason, rec.row_type(), gt=rec.gt_text(haploid=haploid))
        return _from_record(rec, None, haploid=haploid)  # a RefCall's span
    last = pos + max(len(ref), 1) - 1
    blocks = [r for r in records if r.block and r.pos <= last and pos <= r.end]
    if not any(b.pos <= pos <= b.end for b in blocks):
        return _unknown(NO_COVERAGE_REASON)
    min_dp = min((b.min_dp for b in blocks if b.min_dp is not None), default=None)
    # A no-call record on a base the claim deletes leaves it unknown, as a no-call block does.
    deleted = [r for r in variants if pos < r.pos <= last]
    nocall = next((r for r in blocks + deleted if r.nocall), None)
    if nocall is not None:
        return _unknown(
            NOCALL_REASON,
            RowType.NOCALL,
            gt=nocall.gt_text(haploid=haploid),
            gq=nocall.gq,
            dp=min_dp,
        )
    if min_dp is None or min_dp == 0:
        return _unknown(NO_COVERAGE_REASON, dp=min_dp)
    gqs = [b.gq for b in blocks if b.gq is not None]
    return Called(
        gt=blocks[0].gt_text(haploid=haploid),
        dosage="0",
        row_type=RowType.REF_BLOCK,
        gq=min(min(gqs), BLOCK_GQ_CAP) if gqs else None,
        dp=min_dp,
        ad=None,
        vaf=None,
        warnings=(),
    )


def calls_at(
    gvcf: Path,
    alleles: Iterable[Allele],
    *,
    fasta: Path,
    haploid: Callable[[str, int], bool],
    scratch: Path,
) -> dict[Allele, Called]:
    """The call at each `(chrom, pos, ref, alt)` allele, keyed by the allele as given.

    `gvcf` is DeepVariant's indexed gVCF on `fasta`; `haploid(chrom, pos)` says whether a
    locus has one copy (asked at the normalised position, in the gVCF's contig names);
    `scratch` takes `normalise`'s files."""
    wanted = list(dict.fromkeys(alleles))
    with pysam.VariantFile(str(gvcf)) as vcf:
        contigs = {
            str(name): int(c.length or 0) for name, c in vcf.header.contigs.items()
        }

        def contig(chrom: str) -> str | None:
            for name in (chrom, ucsc_contig(chrom)):
                if name in contigs:
                    return name
            return None

        known = [(a, contig(a[0])) for a in wanted]
        on_contig = [(a, (c, a[1], a[2].upper(), a[3].upper())) for a, c in known if c]
        claims = normalise(
            [raw for _, raw in on_contig], fasta=fasta, contigs=contigs, scratch=scratch
        )
        nearby: dict[Allele, list[_Record]] = {}
        for claim in claims:
            chrom, pos, ref, _ = claim
            start, stop = pos - 1, pos - 1 + max(len(ref), 1) + NORM_WINDOW
            nearby[claim] = [_record(rec) for rec in vcf.fetch(chrom, start, stop)]
    records = list(
        dict.fromkeys(r for rs in nearby.values() for r in rs if not r.block)
    )
    keys = [
        (r, i)
        for r in records
        for i, alt in enumerate(r.alts, start=1)
        if alt not in (_SYMBOLIC, _STAR)
    ]
    normed = normalise(
        [(r.chrom, r.pos, r.ref, r.alts[i - 1]) for r, i in keys],
        fasta=fasta,
        contigs=contigs,
        scratch=scratch,
    )
    normal = dict(zip(keys, normed, strict=True))
    out: dict[Allele, Called] = {}
    for allele, c in known:
        if c is None:
            out[allele] = _unknown(NO_COVERAGE_REASON)
    for (allele, _), claim in zip(on_contig, claims, strict=True):
        out[allele] = _call(
            claim, nearby[claim], normal, haploid=haploid(claim[0], claim[1])
        )
    return out
