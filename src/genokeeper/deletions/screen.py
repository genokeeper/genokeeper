"""The deletion screen: which called deletions are shown, under which genes, in which class.

Coordinates are VCF's: `pos` is the base before the deletion and `end` its last deleted base
(both 1-based), so the deleted bases are `[pos, end)` 0-based and `end - pos` long. Only
three filters cut (`SHORT`, `NO_EXON`, `COMMON`, in that order; design decision 5); every
other feature that lowers confidence is a `CallClass`, whose chance of being false the
calibration measures and whose `warnings` the call carries. Interval lookups bisect a
per-chromosome list sorted by start (`_Index`).

Slice 3d runs the same screen over duplications and insertions (`Deletion.svtype`) and over
other pairs of callers (`Screen.callers`). An insertion is a point, `end = pos + 1`, which
touches a target as a base would; its length is its `svlen`, and it is the same as a
gnomAD-SV insertion within `InsertionThresholds.max_distance` of it (a size ratio as well
unless its length is only a lower bound, C20), not by overlap. A caller in `LABELS`
(Scramble) only labels a call another caller made, which keeps its class; a call only it
made is a class of its own (band `mei`, its length unknown, never cut as short).
"""

from bisect import bisect_left
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cached_property

from genokeeper.config import DeletionThresholds, InsertionThresholds
from genokeeper.deletions.exons import Target
from genokeeper.regions import Interval

SHORT = "under 50 bp"
NO_EXON = "no listed exon"
COMMON = "common in gnomAD-SV"
SAME_POSITION = "another deletion at the same position"
DROP_REASONS = (SHORT, NO_EXON, COMMON, SAME_POSITION)
CALLER_NAMES = {
    "delly": "DELLY",
    "dysgu": "Dysgu",
    "dellycnv": "DELLY cnv",
    "cnvpytor": "CNVpytor",
    "insurveyor": "INSurVeyor",
    "scramble": "Scramble",
}
# The depth line of a class whose evidence lacks `+depth` (`warnings`): in the long band
# depth should move and does not; under `long_from` a deletion often moves it too little.
NO_DEPTH = "read depth does not support it"
SHORT_DEPTH = "read depth does not show it, which a deletion this short often cannot"
# Callers that label a call rather than find one (slice 3d: Scramble names a mobile
# element's family, the plan's "Insertion class features"): beside the screen's callers they
# change no class; alone they make the `mei` band, evidence the caller.
LABELS = frozenset({"scramble"})


@dataclass(frozen=True)
class Deletion:
    """One deletion either caller made; `callers` maps `"dysgu"` and/or `"delly"` to that
    caller's id and QUAL, `depth_ratio` is DELLY's inside read count over its two flanks'
    together (None without a DELLY record). Of another `svtype` (`"DUP"`, `"INS"`), the same
    shape: an insertion's `end` is `pos + 1` and `svlen` its inserted length, which
    `lower_bound` marks as only a minimum."""

    chrom: str
    pos: int
    end: int
    ref: str
    gt: str
    callers: Mapping[str, str]
    depth_ratio: float | None
    svtype: str = "DEL"
    svlen: int | None = None
    lower_bound: bool = False

    @property
    def length(self) -> int:
        return self.end - self.pos if self.svlen is None else self.svlen


@dataclass(frozen=True)
class PopulationSv:
    """A gnomAD-SV record already past its PASS, `SVTYPE` and `AN` filters, in the same VCF
    coordinates as `Deletion` (an insertion's `svlen` its `SVLEN`)."""

    id: str
    chrom: str
    pos: int
    end: int
    af: float
    svtype: str = "DEL"
    svlen: int | None = None


@dataclass(frozen=True)
class CallClass:
    """A shown deletion's class: its size band, its evidence (which callers, and whether
    DELLY's depth supports it: one feature, as depth exists only where DELLY called), an end
    in a tandem repeat (short band only) and an end in a segdup or low-mappability region;
    for a depth gain, whether a DUP junction supports it (always False for a deletion)."""

    band: str
    evidence: str
    repeat: bool
    masked: bool
    junction: bool = False


class _Index[T]:
    """Intervals per chromosome sorted by start; `hits` bisects for those overlapping a span."""

    def __init__(self, rows: Iterable[tuple[str, int, int, T]]) -> None:
        by_chrom: dict[str, list[tuple[int, int, T]]] = {}
        for chrom, start, end, value in rows:
            by_chrom.setdefault(chrom, []).append((start, end, value))
        self._rows = {c: sorted(r, key=lambda x: x[:2]) for c, r in by_chrom.items()}
        self._starts = {c: [r[0] for r in rows] for c, rows in self._rows.items()}
        self._longest = {
            c: max(e - s for s, e, _ in rows) for c, rows in self._rows.items()
        }

    def hits(self, chrom: str, start: int, end: int) -> list[T]:
        """The values of the intervals overlapping 0-based `[start, end)`."""
        rows = self._rows.get(chrom)
        if not rows:
            return []
        starts = self._starts[chrom]
        # An interval starting at or before `start - longest` ends at or before `start`.
        lo = bisect_left(starts, start - self._longest[chrom] + 1)
        hi = bisect_left(starts, end)
        return [value for _, e, value in rows[lo:hi] if e > start]

    def covers(self, chrom: str, base: int) -> bool:
        return bool(self.hits(chrom, base, base + 1))


@dataclass(frozen=True)
class Screen:
    """What the screen reads: gene -> its targets (exons widened by the margin), chrom -> the
    mask (segdups and low mappability) and the repeats (GIAB's tandem repeats and
    homopolymers), gnomAD-SV's records of the calls' type, the thresholds, and the two
    callers whose union the calls are (`classify`)."""

    targets: Mapping[str, Sequence[Target]]
    mask: Mapping[str, Sequence[Interval]]
    repeats: Mapping[str, Sequence[Interval]]
    population: Sequence[PopulationSv]
    thresholds: DeletionThresholds
    callers: tuple[str, str] = ("dysgu", "delly")

    @cached_property
    def _targets(self) -> _Index[tuple[str, Target]]:
        return _Index(
            (t.chrom, t.start, t.end, (gene, t))
            for gene, ts in self.targets.items()
            for t in ts
        )

    @cached_property
    def _population(self) -> _Index[tuple[int, PopulationSv]]:
        return _Index(
            (p.chrom, p.pos, p.end, (i, p)) for i, p in enumerate(self.population)
        )

    @cached_property
    def _mask(self) -> _Index[None]:
        return _Index((c, s, e, None) for v in self.mask.values() for c, s, e in v)

    @cached_property
    def _repeats(self) -> _Index[None]:
        return _Index((c, s, e, None) for v in self.repeats.values() for c, s, e in v)


@dataclass(frozen=True)
class Shown:
    """A deletion the screen shows: its class, gene -> exon numbers hit (sorted genes), the
    highest-AF gnomAD-SV record it matches, rare or not, for a depth event gene -> its
    kind there (`loss`, `tandem`, `intragenic`, `whole`, `partial`; empty for 3b's), and
    the genes whose targets it touches only in the margin beside an exon (`flank_only`)."""

    deletion: Deletion
    call_class: CallClass
    exons: Mapping[str, tuple[int, ...]]
    population: PopulationSv | None
    kinds: Mapping[str, str] = field(default_factory=dict[str, str])
    flank: frozenset[str] = frozenset()


def key(c: CallClass) -> str:
    """A stable name of the class, for `INFO/CLASS`."""
    repeat = "repeat" if c.repeat else "norepeat"
    masked = "masked" if c.masked else "unmasked"
    junction = "junction" if c.junction else "nojunction"
    return f"{c.band}:{c.evidence}:{repeat}:{masked}:{junction}"


def warnings(c: CallClass, callers: tuple[str, str] = ("dysgu", "delly")) -> list[str]:
    """One line per feature of the class that lowers confidence; the depth lines only where
    DELLY is one of the `callers` (its depth ratio, 3b's one depth feature)."""
    lines: list[str] = []
    found = c.evidence.removesuffix("+depth")
    if found != "both":
        lines.append(f"found by one caller ({CALLER_NAMES[found]})")
    if "delly" in callers:
        if found not in ("both", "delly"):
            lines.append("read depth not measured")
        elif not c.evidence.endswith("+depth"):
            lines.append(SHORT_DEPTH if c.band == "short" else NO_DEPTH)
    if c.repeat:
        lines.append("an end in a tandem repeat: may be a repeat-length change")
    if c.masked:
        lines.append("an end in a segmental duplication or low-mappability region")
    return lines


def classify(d: Deletion, s: Screen) -> CallClass:
    """The deletion's class; an end base is its first deleted base or its last. Its
    evidence is `both` when both of `s.callers` called it, else the one caller's name, and
    `+depth` where DELLY called it and its depth ratio supports it; `LABELS` callers count
    only alone (band `mei`, evidence that caller)."""
    t = s.thresholds
    ends = (d.pos, d.end - 1)
    masked = any(s._mask.covers(d.chrom, b) for b in ends)
    named = set(d.callers)
    callers, pair = named - LABELS, set(s.callers)
    if not callers and len(named) == 1:
        return CallClass(band="mei", evidence=named.pop(), repeat=False, masked=masked)
    band = "long" if d.length >= t.long_from else "short"
    if callers == pair:
        evidence = "both"
    elif len(callers) == 1 and callers <= pair:
        evidence = callers.pop()
    else:
        raise ValueError(f"{d.chrom}:{d.pos}: unknown callers {sorted(callers)}")
    ratio = d.depth_ratio
    if "delly" in d.callers and ratio is not None and ratio < t.max_depth_ratio:
        evidence += "+depth"
    return CallClass(
        band=band,
        evidence=evidence,
        repeat=band == "short" and any(s._repeats.covers(d.chrom, b) for b in ends),
        masked=masked,
    )


def max_distance(d: Deletion, t: DeletionThresholds) -> int:
    """How far from the call a gnomAD-SV record may lie and still be it: `t.max_distance` for
    an insertion (`InsertionThresholds` only), 0 for any other type."""
    if d.svtype != "INS":
        return 0
    if not isinstance(t, InsertionThresholds):
        raise TypeError(f"{d.chrom}:{d.pos}: an insertion needs InsertionThresholds")
    return t.max_distance


def _reciprocal(d: Deletion, p: PopulationSv, distance: int = 0) -> float:
    """Reciprocal overlap for a deletion or duplication; for an insertion within `distance`
    of `p`, 1.0 when either length is unknown or only a minimum (position and type alone),
    else the size ratio. 0 for another chromosome or type."""
    if p.chrom != d.chrom or p.svtype != d.svtype:
        return 0.0
    if d.svtype == "INS":
        if abs(d.pos - p.pos) > distance:
            return 0.0
        if d.lower_bound or p.svlen is None:
            return 1.0
        return min(d.length, p.svlen) / max(d.length, p.svlen)
    overlap = min(d.end, p.end) - max(d.pos, p.pos)
    return max(overlap, 0) / max(d.length, p.end - p.pos)


def common_match(
    d: Deletion,
    population: Sequence[PopulationSv],
    *,
    min_overlap: float,
    distance: int = 0,
) -> PopulationSv | None:
    """The highest-AF record with reciprocal overlap at least `min_overlap` (an insertion's
    within `distance`, `_reciprocal`), if any."""
    matches = [p for p in population if _reciprocal(d, p, distance) >= min_overlap]
    return max(matches, key=lambda p: p.af, default=None)


def is_common(d: Deletion, s: Screen) -> bool:
    """The `COMMON` test alone: a match (`population_match`) above `max_gnomad_af`."""
    match = population_match(d, s)
    return match is not None and match.af > s.thresholds.max_gnomad_af


def is_short(d: Deletion, t: DeletionThresholds) -> bool:
    """The `SHORT` test: under `min_length`, unless its length is unknown (no `svlen`, only
    a bound: a Scramble mobile element)."""
    unknown = d.svlen is None and d.lower_bound
    return d.length < t.min_length and not unknown


def population_spans(calls: Sequence[Deletion], *, pad: int) -> list[Interval]:
    """Each call's `[pos - pad, end + pad)` (clipped at 0): where gnomAD-SV is read for it,
    `pad` being `max_distance` for insertions and 0 for any other type."""
    return [(d.chrom, max(d.pos - pad, 0), d.end + pad) for d in calls]


def rare(deletions: Sequence[Deletion], s: Screen) -> list[Deletion]:
    """The deletions `is_common` finds not common: the run's `COMMON` test, for the
    calibration's genome-wide call sets."""
    return [d for d in deletions if not is_common(d, s)]


def reciprocal(a: Interval, b: Interval) -> float:
    """The reciprocal overlap of two 0-based half-open spans: their shared bases over the
    longer's length; 0 on different chromosomes."""
    if a[0] != b[0]:
        return 0.0
    overlap = min(a[2], b[2]) - max(a[1], b[1])
    return max(overlap, 0) / max(a[2] - a[1], b[2] - b[1], 1)


def exon_hits(d: Deletion, s: Screen) -> dict[str, tuple[int, ...]]:
    """The genes whose targets the call touches, with the exon numbers it hits (sorted)."""
    hit: dict[str, set[int]] = {}
    for gene, t in s._targets.hits(d.chrom, d.pos, d.end):
        hit.setdefault(gene, set()).update(t.exons)
    return {gene: tuple(sorted(hit[gene])) for gene in sorted(hit)}


def flank_only(d: Deletion, s: Screen) -> frozenset[str]:
    """The genes whose targets the call touches only within `exon_margin` of their ends:
    a target is its exons widened by the margin and merged where they touch, so its bases
    past the margin at either end are exon or an intron of at most twice the margin (40 bp
    by default), shorter than any deletion the screen keeps (50 bp)."""
    m = s.thresholds.exon_margin
    hit: dict[str, bool] = {}
    for gene, t in s._targets.hits(d.chrom, d.pos, d.end):
        core = t.start + m < d.end and d.pos < t.end - m
        hit[gene] = hit.get(gene, False) or core
    return frozenset(gene for gene, core in hit.items() if not core)


def population_match(d: Deletion, s: Screen) -> PopulationSv | None:
    """`common_match` with the screen's own overlap and distance, over only the records
    (in the screen's order, so a tie in AF picks as the whole list would) that overlap the
    call (an insertion's, those within `max_distance`, as
    `population_spans` pads it): the same answer as over the whole list, as
    `min_reciprocal_overlap` is above 0 (`config.load_grading` refuses another), in
    O(overlapping records) per call, for the calibration's genome-wide call sets."""
    distance = max_distance(d, s.thresholds)
    ((chrom, start, end),) = population_spans([d], pad=distance)
    return common_match(
        d,
        [p for _, p in sorted(s._population.hits(chrom, start, end))],
        min_overlap=s.thresholds.min_reciprocal_overlap,
        distance=distance,
    )


def screen(d: Deletion, s: Screen) -> tuple[dict[str, tuple[int, ...]], str | None]:
    """The genes the deletion shows under, with the exon numbers it hits; or `{}` and the
    first drop reason (`is_short`, `is_common`)."""
    if is_short(d, s.thresholds):
        return {}, SHORT
    hit = exon_hits(d, s)
    if not hit:
        return {}, NO_EXON
    if is_common(d, s):
        return {}, COMMON
    return hit, None


def one_per_position(shown: Iterable[Shown]) -> tuple[list[Shown], int]:
    """Of the calls shown at one chrom, POS and type (one vid), the one with the most
    callers, then the longest, in the order first seen; and how many others were dropped
    (`SAME_POSITION`)."""
    best: dict[tuple[str, int, str], Shown] = {}
    same = 0
    for item in shown:
        d = item.deletion
        at = (d.chrom, d.pos, d.svtype)
        if at in best:
            same += 1
            kept = best[at].deletion
            if (len(d.callers), d.length) <= (len(kept.callers), kept.length):
                continue
        best[at] = item
    return list(best.values()), same


def screen_all(
    deletions: Sequence[Deletion], s: Screen
) -> tuple[list[Shown], dict[str, int]]:
    """`screen` and `classify` on each deletion, then `one_per_position`."""
    dropped: dict[str, int] = dict.fromkeys(DROP_REASONS, 0)
    candidates: list[Shown] = []
    for d in deletions:
        exons, reason = screen(d, s)
        if reason is not None:
            dropped[reason] += 1
            continue
        candidates.append(
            Shown(
                d,
                classify(d, s),
                exons,
                population_match(d, s),
                flank=flank_only(d, s),
            )
        )
    shown, same = one_per_position(candidates)
    dropped[SAME_POSITION] += same
    return shown, dropped
