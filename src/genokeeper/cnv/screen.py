"""The depth copy-number screen (slice 3d): which DELLY cnv and CNVpytor events are shown,
under which genes and kinds, which are large events, named loci or 3b deletions, and why the
others are dropped.

It is 3b's screen over depth events (`deletions.screen`): the same targets, mask, gnomAD-SV
rule (`is_common`, of the event's own type) and class model, under the depth callers
(`Screen.callers` = `DEPTH_CALLERS`). `screen_cnv` takes each event through, in order:
`SHORT` (under `[cnv] min_length`); `MASKED_REF` (an end, or half its length, in 3c's masked
intervals, C19); `OWNED` (the same in a region a locus-specific caller owns:
`owned_regions`); `COMMON`; a loss that is a shown 3b deletion (`JOINED`, recorded so the
deletion's statement says depth supports it, C6); a loss over a named locus (`named`); a
large event (`[cnv] large_from` or more), stated on its own (`LargeEvent`) only when both
depth callers agree, either no end lies in 3b's mask or it matches a ClinGen recurrent
region (whose LCR-mediated ends lie in segmental duplications by definition), and it
overlaps a ClinGen region or a dosage-sensitive gene, else `LARGE_FILTERED` (never an
ordinary event, whatever exons it covers); `NO_EXON`; then `SAME_POSITION` as 3b's. A gain's class says whether a Dysgu DUP
junction supports it, and its kind per gene (`gain_kind`) decides its placement.
"""

import gzip
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from importlib.resources import files
from pathlib import Path

from genokeeper.config import CnvThresholds
from genokeeper.deletions.exons import OWNED_CYP2D, Target
from genokeeper.deletions.screen import (
    CALLER_NAMES,
    COMMON,
    NO_EXON,
    SAME_POSITION,
    CallClass,
    Deletion,
    Screen,
    Shown,
    _Index,
    classify,
    exon_hits,
    is_common,
    one_per_position,
    population_match,
    reciprocal,
)
from genokeeper.evidence.clingen import DosageList
from genokeeper.regions import Interval

SHORT = "under 1 kb"
MASKED_REF = "in the reference's masked false duplications"
OWNED = "a locus-specific caller owns the region"
JOINED = "a 3b deletion"
LARGE_FILTERED = (
    "large, but one caller only, an end in a segmental duplication outside a recurrent"
    " CNV region, or over no ClinGen region or dosage-sensitive gene"
)
DROP_REASONS = (
    SHORT,
    MASKED_REF,
    OWNED,
    COMMON,
    JOINED,
    LARGE_FILTERED,
    NO_EXON,
    SAME_POSITION,
)
DEPTH_CALLERS = ("dellycnv", "cnvpytor")
# Each direction's three size bands (the plan's Decisions, "Depth bands"); no `:`, the
# class key's separator.
CNV_BANDS = tuple(
    f"{direction}_{size}"
    for direction in ("loss", "gain")
    for size in ("1-5kb", "5kb+", "large")
)
# Every value each class feature of a depth event takes (`classify_cnv`).
CNV_FEATURE_VALUES: dict[str, tuple[str, ...]] = {
    "evidence": ("both", *DEPTH_CALLERS),
    "repeat": ("false",),
    "masked": ("true", "false"),
    "junction": ("true", "false"),  # always false for a loss
}
# The reciprocal overlap at which a depth loss is a named locus, and a gain matches a DUP
# junction or a large event a ClinGen recurrent region (the plan's Decisions).
MATCH_OVERLAP = 0.5
NAMED = "named.tsv"
NAMED_COLUMNS = ("name", "chrom", "start", "end", "gene", "reason", "homozygous_reason")
SMN_GENES = ("SMN1", "SMN2")
_DIRECTION = {"DEL": "loss", "DUP": "gain"}


@dataclass(frozen=True)
class NamedCnv:
    """A named copy-number locus: its span (0-based half-open), the gene a loss there is
    stated on, and the reason that statement carries on one copy and on both."""

    name: str
    chrom: str
    start: int
    end: int
    gene: str
    reason: str
    homozygous_reason: str


@dataclass(frozen=True)
class LargeEvent:
    """A large event the screen shows on its own: the event, its class, and the ClinGen
    regions and dosage-sensitive genes it overlaps."""

    event: Deletion
    call_class: CallClass
    regions: tuple[str, ...]


@dataclass(frozen=True)
class CnvScreen:
    """What the screen made of the depth events: the shown ones (3b's shape, with each
    gene's kind), the large events, the named loci with the event over each, each shown 3b
    deletion's (chrom, pos, end) -> the depth callers that joined it, and the drop counts."""

    shown: list[Shown]
    large: list[LargeEvent]
    named: list[tuple[NamedCnv, Shown]]
    joined: dict[tuple[str, int, int], str]
    dropped: dict[str, int]


def load_named() -> list[NamedCnv]:
    """The packaged named loci (`resources/cnv/named.tsv`)."""
    text = files("genokeeper").joinpath("resources", "cnv", NAMED).read_text()
    rows: list[NamedCnv] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if tuple(fields) == NAMED_COLUMNS:
            continue
        name, chrom, start, end, gene, reason, homozygous = fields
        rows.append(
            NamedCnv(name, chrom, int(start), int(end), gene, reason, homozygous)
        )
    return rows


def _model_region(text: str) -> Interval:
    """Parascopy's `chr:start-end[:strand]` (1-based inclusive) as 0-based half-open."""
    chrom, span, *_ = text.split(":")
    start, end = span.split("-")
    return chrom, int(start) - 1, int(end)


def read_model(model: Path) -> tuple[Interval, list[Interval]]:
    """A Parascopy model's `main` region and both copies of each `duplication` line. A
    model with no `main` line is refused: its anchor comes from Parascopy's model format,
    which a truncated or other file cannot pass because a model always opens with its
    `main` region."""
    main: Interval | None = None
    copies: list[Interval] = []
    with gzip.open(model, "rt") as handle:
        for line in handle:
            kind, *rest = line.rstrip("\n").split("\t")
            if kind == "main":
                main = _model_region(rest[0])
            elif kind == "duplication":
                copies += [_model_region(rest[0]), _model_region(rest[1])]
    if main is None:
        raise ValueError(f"{model}: no Parascopy 'main' line")
    return main, copies


def owned_regions(
    spans: Mapping[str, tuple[str, int, int]], models: Sequence[Path]
) -> list[Interval]:
    """The regions a locus-specific caller owns: SMN1's and SMN2's spans (SMNCopyNumberCaller
    and Parascopy, 3a), each Parascopy model's `main` region and both copies of each
    `duplication` line (`read_model`), and the CYP2D locus (`OWNED_CYP2D`)."""
    regions: list[Interval] = [spans[gene] for gene in SMN_GENES]
    for model in models:
        main, copies = read_model(model)
        regions += [main, *copies]
    regions.append(OWNED_CYP2D)
    return regions


def merge_depth(
    dellycnv: Sequence[Deletion],
    cnvpytor: Sequence[Deletion] | None,
    *,
    min_overlap: float,
) -> list[Deletion]:
    """One event per DELLY event and the CNVpytor event of its type it overlaps most, at
    reciprocal overlap `min_overlap` or more: DELLY's coordinates and genotype, both
    callers; each CNVpytor event joins one at most; the rest alone, DELLY's first. With
    `cnvpytor` None (CNVpytor dropped by the calibration), DELLY's alone."""
    if cnvpytor is None:
        return list(dellycnv)
    index = _Index((e.chrom, e.pos, e.end, i) for i, e in enumerate(cnvpytor))
    used: set[int] = set()
    merged: list[Deletion] = []
    for d in dellycnv:
        span = (d.chrom, d.pos, d.end)
        overlaps = [
            (reciprocal(span, (e.chrom, e.pos, e.end)), i)
            for i in index.hits(d.chrom, d.pos, d.end)
            if i not in used and (e := cnvpytor[i]).svtype == d.svtype
        ]
        best = max(overlaps, default=None)
        if best is None or best[0] < min_overlap:
            merged.append(d)
            continue
        used.add(best[1])
        merged.append(replace(d, callers={**d.callers, **cnvpytor[best[1]].callers}))
    return merged + [e for i, e in enumerate(cnvpytor) if i not in used]


def cnv_band(d: Deletion, t: CnvThresholds) -> str:
    """A depth event's band by its direction and length: `<loss|gain>_<1-5kb|5kb+|large>`."""
    size = (
        "large"
        if d.length >= t.large_from
        else "5kb+"
        if d.length >= t.band_from
        else "1-5kb"
    )
    return f"{_DIRECTION[d.svtype]}_{size}"


def classify_cnv(
    d: Deletion, s: Screen, t: CnvThresholds, *, junction: bool
) -> CallClass:
    """A depth event's class: its band (`cnv_band`), its evidence from `s.callers` and an
    end in 3b's mask (`classify`), no repeat feature, and for a gain whether a DUP junction
    supports it (a loss with a junction is a 3b deletion)."""
    return replace(
        classify(d, s),
        band=cnv_band(d, t),
        repeat=False,
        junction=junction and d.svtype == "DUP",
    )


def gain_kind(d: Deletion, span: tuple[str, int, int], junction: bool) -> str:
    """A gain's kind in a gene of `span`: `tandem` (both ends inside and a junction),
    `intragenic` (both inside, no junction), `whole` (both outside), else `partial`."""
    _, start, end = span
    first, last = d.pos, d.end - 1
    if start <= first and last < end:
        return "tandem" if junction else "intragenic"
    if first < start and last >= end:
        return "whole"
    return "partial"


def gene_spans(
    targets: Mapping[str, Sequence[Target]],
) -> dict[str, tuple[str, int, int]]:
    """Each gene's extent: its first target's start to its last target's end (exons and
    their margins), for `gain_kind`."""
    return {
        gene: (ts[0].chrom, min(t.start for t in ts), max(t.end for t in ts))
        for gene, ts in targets.items()
        if ts
    }


def _interval_index(intervals: Iterable[Interval]) -> _Index[Interval]:
    return _Index((c, s, e, (c, s, e)) for c, s, e in intervals)


def _covered(index: _Index[Interval], chrom: str, start: int, end: int) -> int:
    """How many bases of `[start, end)` the intervals cover, overlaps counted once."""
    covered, reach = 0, start
    for _, s, e in sorted(index.hits(chrom, start, end), key=lambda x: x[1]):
        lo, hi = max(s, reach), min(e, end)
        if hi > lo:
            covered += hi - lo
            reach = hi
    return covered


def _within(d: Deletion, index: _Index[Interval]) -> bool:
    """An end base of the event, or half its length, inside the intervals."""
    ends = (d.pos, d.end - 1)
    if any(index.hits(d.chrom, b, b + 1) for b in ends):
        return True
    return 2 * _covered(index, d.chrom, d.pos, d.end) >= d.end - d.pos


def segdup_exons(
    targets: Mapping[str, Sequence[Target]],
    mask: Sequence[Interval],
    owned: Sequence[Interval],
) -> dict[str, tuple[int, ...]]:
    """Per gene, the exon numbers of its targets wholly inside the mask and outside every
    owned region: target exons in 3b's mask (segmental duplications or low mappability) no
    Parascopy model covers, which the depth screen names as not covered."""
    masked, claimed = _interval_index(mask), _interval_index(owned)
    found: dict[str, tuple[int, ...]] = {}
    for gene, ts in targets.items():
        numbers = sorted(
            {
                n
                for t in ts
                if _covered(masked, t.chrom, t.start, t.end) == t.end - t.start
                and not claimed.hits(t.chrom, t.start, t.end)
                for n in t.exons
            }
        )
        if numbers:
            found[gene] = tuple(numbers)
    return found


def _named(d: Deletion, named: Sequence[NamedCnv]) -> NamedCnv | None:
    """The named locus the loss is, at reciprocal overlap `MATCH_OVERLAP` or more."""
    span = (d.chrom, d.pos, d.end)
    matches = [
        (o, n)
        for n in named
        if (o := reciprocal(span, (n.chrom, n.start, n.end))) >= MATCH_OVERLAP
    ]
    return max(matches, key=lambda m: m[0])[1] if matches else None


def named_deletions(
    deletions: Sequence[Deletion], named: Sequence[NamedCnv], s: Screen
) -> list[tuple[NamedCnv, Shown]]:
    """3b's deletions that are a named locus, before 3b's `NO_EXON` drop (a named locus's
    gene need not be screened), each in its 3b class; common ones are not."""
    found: list[tuple[NamedCnv, Shown]] = []
    for d in deletions:
        row = _named(d, named)
        if row is not None and not is_common(d, s):
            found.append((row, Shown(d, classify(d, s), {}, population_match(d, s))))
    return found


def _names(callers: Iterable[str]) -> str:
    return " and ".join(CALLER_NAMES[c] for c in DEPTH_CALLERS if c in set(callers))


def screen_cnv(
    events: Sequence[Deletion],
    deletions: Sequence[Deletion],
    junctions: Sequence[Deletion],
    s: Screen,
    t: CnvThresholds,
    *,
    masked: Sequence[Interval],
    owned: Sequence[Interval],
    spans: Mapping[str, tuple[str, int, int]],
    dosage: DosageList,
    regions: Sequence[tuple[str, int, int, str]],
    recurrent: Sequence[tuple[str, int, int, str]],
    named: Sequence[NamedCnv],
) -> CnvScreen:
    """Screen the depth events (see the module docstring). `deletions` are 3b's shown
    deletions, `junctions` Dysgu's DUP records, `masked` 3c's masked intervals, `owned`
    `owned_regions`, `spans` each screened gene's span (`gene_spans`), `regions` and
    `recurrent` ClinGen's (`load_regions`, `load_recurrent`)."""
    d_t = s.thresholds
    dropped: dict[str, int] = dict.fromkeys(DROP_REASONS, 0)
    masked_at, owned_at = _interval_index(masked), _interval_index(owned)
    deletion_at = _Index((x.chrom, x.pos, x.end, x) for x in deletions)
    junction_at = _Index(
        (j.chrom, j.pos, j.end, j) for j in junctions if j.svtype == "DUP"
    )
    sensitive = dosage.haploinsufficient | dosage.triplosensitive
    region_at = _Index(
        [
            *regions,
            *((*dosage.spans[g], g) for g in sorted(sensitive) if g in dosage.spans),
        ]
    )
    recurrent_at = _Index((c, a, b, (c, a, b)) for c, a, b, _ in recurrent)
    joined: dict[tuple[str, int, int], set[str]] = {}
    large: list[LargeEvent] = []
    found: list[tuple[NamedCnv, Shown]] = []
    candidates: list[Shown] = []
    for d in events:
        span = (d.chrom, d.pos, d.end)
        if d.length < t.min_length:
            dropped[SHORT] += 1
            continue
        if _within(d, masked_at):
            dropped[MASKED_REF] += 1
            continue
        if _within(d, owned_at):
            dropped[OWNED] += 1
            continue
        if is_common(d, s):
            dropped[COMMON] += 1
            continue
        if d.svtype == "DEL":
            partners = [
                (reciprocal(span, (x.chrom, x.pos, x.end)), (x.chrom, x.pos, x.end))
                for x in deletion_at.hits(d.chrom, d.pos, d.end)
            ]
            best = max(partners, default=None)
            if best is not None and best[0] >= d_t.min_reciprocal_overlap:
                joined.setdefault(best[1], set()).update(d.callers)
                dropped[JOINED] += 1
                continue
            row = _named(d, named)
            if row is not None:
                c = classify_cnv(d, s, t, junction=False)
                found.append((row, Shown(d, c, {}, population_match(d, s))))
                continue
        junction = d.svtype == "DUP" and any(
            reciprocal(span, (j.chrom, j.pos, j.end)) >= MATCH_OVERLAP
            for j in junction_at.hits(d.chrom, d.pos, d.end)
        )
        c = classify_cnv(d, s, t, junction=junction)
        if d.length >= t.large_from:
            ends_masked = any(s._mask.covers(d.chrom, b) for b in (d.pos, d.end - 1))
            recurrent_match = any(
                reciprocal(span, r) >= MATCH_OVERLAP
                for r in recurrent_at.hits(d.chrom, d.pos, d.end)
            )
            over = tuple(dict.fromkeys(region_at.hits(d.chrom, d.pos, d.end)))
            if (
                c.evidence != "both"
                or (ends_masked and not recurrent_match)
                or not over
            ):
                dropped[LARGE_FILTERED] += 1
            else:
                large.append(LargeEvent(d, c, over))
            continue
        exons = exon_hits(d, s)
        if not exons:
            dropped[NO_EXON] += 1
            continue
        kinds = {
            gene: "loss" if d.svtype == "DEL" else gain_kind(d, spans[gene], junction)
            for gene in exons
        }
        candidates.append(Shown(d, c, exons, population_match(d, s), kinds=kinds))
    shown, same = one_per_position(candidates)
    dropped[SAME_POSITION] += same
    return CnvScreen(
        shown=shown,
        large=large,
        named=found,
        joined={at: _names(callers) for at, callers in joined.items()},
        dropped=dropped,
    )
