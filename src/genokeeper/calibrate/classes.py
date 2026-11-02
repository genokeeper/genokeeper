"""The class model every structural calibration table shares (slice 3b's, generalised in 3d):
each class's chance of being false, measured on calibration libraries against a truth,
under a feature schema each table names (`Features`: feature -> its values).

Within each band the expected false calls over the targets (the band's false calls per
genome inside the benchmark times `h`, the chance a call of its length touches a target) and
the expected true ones (the band's prior, the rare true events a genome carries over the
targets, times the band's recall) are each moved by the class's features' likelihoods, taken
as independent (naive Bayes); every term's Jeffreys bound, pushed to the same side, makes the
class's interval, so a class measured on few calls gets a wide one (`ClassTable.p_wrong`).

A table (`ClassTable`) is written and opened through 3c's table module
(`genokeeper.calibrate.tables`): `write_class_table` swaps it in whole, `open_class_table`
returns it only when it was measured under the versions the run computed itself. Format 2
adds the schema (`features`) and a table's own fields (`extra`); a table every run must
re-measure after a format change names the format among its versions (`FORMAT_KEY`), so a
table written in another format is stale under 3c's own check, and `read_class_table`
refuses any other format besides. This module imports nothing from `pipeline`.
"""

import json
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Protocol

from genokeeper.calibrate.table import jeffreys
from genokeeper.calibrate.tables import REPORT, RUN, Current, open_table, write_table
from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import CallClass, Deletion, key
from genokeeper.regions import Interval, merge, subtract

TABLE = "table.json"
FORMAT = 2
# The versions key naming the table format a run reads: a table of another format is stale.
FORMAT_KEY = "table_format"
type Features = Mapping[str, tuple[str, ...]]
# `ClassTable.extra`'s key for each band's target parts outside the region the band was
# measured in (`outside_benchmark`), the versions value a table holding it is measured under
# (a table without it is stale, not read as all inside), and the two lines
# `ClassTable.cautions` adds where a chance of being false is extrapolated.
OUTSIDE_BENCHMARK = "outside_benchmark"
OUTSIDE_VERSION = "outside_benchmark per band"
UNSEEN = "no calibration call of this class: its chance of being false is extrapolated"
OUTSIDE = (
    "over a target outside the benchmark the chance of being false was measured in"
)


class Banded(Protocol):
    """A call class: any frozen dataclass with a band and the schema's features."""

    @property
    def band(self) -> str: ...


@dataclass(frozen=True)
class Share:
    """A share measured as `count` in `total`: smoothed point, Jeffreys 95% interval."""

    point: float
    lo: float
    hi: float


def share(count: int, total: int) -> Share:
    """`(count + 0.5) / (total + 1)`, never 0, with `jeffreys(count, total)`; nothing
    measured (`total` 0) is `Share(0.5, 0, 1)`, the widest interval."""
    if total == 0:
        return Share(0.5, 0.0, 1.0)
    return Share((count + 0.5) / (total + 1), *jeffreys(count, total))


@dataclass(frozen=True)
class BandTerms:
    """A band's terms: its false calls among its rare calls inside the benchmark, pooled
    over `libraries` (per genome: `point * calls / libraries`), its recall of the truth's
    rare events, `h` (None without a call) and the prior."""

    false_calls: Share
    calls: int
    libraries: int
    recall: Share
    h: float | None
    prior: float


@dataclass(frozen=True)
class Likelihoods:
    """Per feature and value (the table's `features`), its share of the band's false calls
    and of its true calls."""

    false: dict[str, dict[str, Share]]
    true: dict[str, dict[str, Share]]


def feature_value(c: object, feature: str) -> str:
    """The class's value of `feature` as a schema names it (a bool as `true`/`false`)."""
    value = getattr(c, feature)
    return str(value).lower() if isinstance(value, bool) else str(value)


def _ratio(false: float, true: float, empty: float) -> float:
    return empty if false + true == 0 else false / (false + true)


@dataclass(frozen=True)
class ClassTable:
    """The versions measured under, the terms and likelihoods per band, what was measured
    on each library (for the report), the feature schema and the table's own fields."""

    versions: dict[str, str]
    bands: dict[str, BandTerms]
    likelihoods: dict[str, Likelihoods]
    measured: dict[str, Any]
    features: Features
    extra: dict[str, Any] = field(default_factory=dict[str, Any])

    def expected(self, c: Banded, x: str) -> tuple[float, float]:
        """The class's expected false and true calls over the targets per genome, every
        term at its `x` ("point", "lo" or "hi"); the band must be measured."""
        band, lk = self.bands[c.band], self.likelihoods[c.band]
        if band.h is None:
            raise ValueError(f"band {c.band} not measured")
        values = [(f, feature_value(c, f)) for f in self.features]
        false = math.prod(getattr(lk.false[f][v], x) for f, v in values)
        true = math.prod(getattr(lk.true[f][v], x) for f, v in values)
        per_genome = getattr(band.false_calls, x) * band.calls / band.libraries
        return per_genome * band.h * false, band.prior * getattr(band.recall, x) * true

    def recall(self) -> dict[str, dict[str, tuple[int, int]]]:
        """Per band measured, each library's truth events the screen found and their number
        (`measured`, `summary`'s shape): the counts the band's recall is made of. Not summed:
        two libraries of one genome would count each truth event twice."""
        found: dict[str, dict[str, tuple[int, int]]] = {}
        for library, m in sorted(self.measured.items()):
            for band, counts in m.items():
                if isinstance(counts, dict) and {"found", "truth"} <= counts.keys():
                    found.setdefault(band, {})[library] = (
                        counts["found"],
                        counts["truth"],
                    )
        return found

    def p_wrong(self, c: Banded) -> tuple[float, float, float] | None:
        """The class's chance of being false `F / (F + T)` and its interval: every term's
        Jeffreys bound pushed to the same side (the false terms low against the true ones
        high, and the other way), a range of unknown coverage that leaves out the
        uncertainty in `h`, in the prior and in the features' correlation. None when its
        band has no call or no `h` (not measured)."""
        band = self.bands.get(c.band)
        if band is None or band.calls == 0 or band.h is None:
            return None
        point = _ratio(*self.expected(c, "point"), 1.0)
        f_lo, t_lo = self.expected(c, "lo")
        f_hi, t_hi = self.expected(c, "hi")
        return point, _ratio(f_lo, t_hi, 0.0), _ratio(f_hi, t_lo, 1.0)

    def cautions(self, d: Deletion, c: CallClass) -> list[str]:
        """`UNSEEN` when no library's calibration calls had the class `c` (`summary`'s
        `classes`), `OUTSIDE` when `d`'s bases touch a target part outside the region `c`'s
        band was measured in (`extra[OUTSIDE_BENCHMARK]`): where `p_wrong` is
        extrapolated, which it does not say itself."""
        seen = any(
            key(c) in m.get(c.band, {}).get("classes", ())
            for m in self.measured.values()
        )
        outside = any(
            chrom == d.chrom and start < d.end and d.pos < end
            for chrom, start, end in self.extra.get(OUTSIDE_BENCHMARK, {}).get(
                c.band, ()
            )
        )
        return [*([] if seen else [UNSEEN]), *([OUTSIDE] if outside else [])]


@dataclass(frozen=True)
class Geometry:
    """The targets merged over all genes (sorted, 0-based half-open) and the benchmark's
    bases."""

    targets: tuple[Interval, ...]
    region_bp: int

    def reach(self, length: int) -> int:
        """The bases where an event of `length` could start (its first base) and touch a
        target: the union over every target `[start, end)` of `[start - length, end)`,
        clipped at 0, so bases two targets share count once."""
        total, chrom, first, last = 0, None, 0, 0
        for c, s, e in self.targets:
            start = max(s - length, 0)
            if c == chrom and start <= last:
                last = max(last, e)
                continue
            total += last - first
            chrom, first, last = c, start, e
        return total + last - first


def geometry(
    targets: Mapping[str, Sequence[Target]], region: Iterable[Interval]
) -> Geometry:
    """The targets merged over all genes, and the benchmark BED's size (`region`)."""
    merged = merge((t.chrom, t.start, t.end) for ts in targets.values() for t in ts)
    return Geometry(
        targets=tuple(merged),
        region_bp=sum(e - s for _, s, e in merge(region)),
    )


def outside_benchmark(
    targets: Mapping[str, Sequence[Target]], regions: Mapping[str, Iterable[Interval]]
) -> dict[str, list[list[str | int]]]:
    """Per band, the targets' parts outside the region it was measured in (`regions`, the
    same its geometry reads), merged and sorted, as a table's `extra` stores them."""
    spans = [(t.chrom, t.start, t.end) for ts in targets.values() for t in ts]
    return {
        band: [list(x) for x in subtract(spans, region)]
        for band, region in regions.items()
    }


def measure(
    calls: Mapping[str, Sequence[tuple[Any, int, bool]]],
    truth: Mapping[str, Sequence[tuple[str | None, bool]]],
    geometries: Mapping[str, Geometry],
    priors: Mapping[str, float],
    features: Features,
    *,
    band: Callable[[Any], str],
) -> tuple[dict[str, BandTerms], dict[str, Likelihoods]]:
    """Per band of `geometries`, its terms and likelihoods from each library's rare calls
    inside the benchmark (class, length, confirmed by the truth; the band `band(class)`) and
    its rare truth events (band, found; a band None counts in no band's recall). `h` is the
    mean over the band's false calls (all its calls when none is false) of
    `Geometry.reach(L) / region_bp` on the band's own geometry, at most 1."""
    bands: dict[str, BandTerms] = {}
    likelihoods: dict[str, Likelihoods] = {}
    for name, g in geometries.items():
        rows = [
            (c, n, ok)
            for found in calls.values()
            for c, n, ok in found
            if band(c) == name
        ]
        false = [c for c, _, ok in rows if not ok]
        true = [c for c, _, ok in rows if ok]
        lengths = [n for _, n, ok in rows if not ok] or [n for _, n, _ in rows]
        h = (
            None
            if not lengths
            else sum(min(g.reach(n), g.region_bp) / g.region_bp for n in lengths)
            / len(lengths)
        )
        found = [ok for one in truth.values() for b, ok in one if b == name]
        bands[name] = BandTerms(
            false_calls=share(len(false), len(rows)),
            calls=len(rows),
            libraries=len(calls),
            recall=share(sum(found), len(found)),
            h=h,
            prior=priors[name],
        )

        def shares(classes: Sequence[object]) -> dict[str, dict[str, Share]]:
            return {
                f: {
                    v: share(
                        sum(feature_value(c, f) == v for c in classes), len(classes)
                    )
                    for v in values
                }
                for f, values in features.items()
            }

        likelihoods[name] = Likelihoods(false=shares(false), true=shares(true))
    return bands, likelihoods


def summary(
    calls: Sequence[tuple[Any, int, bool]],
    truth: Sequence[tuple[str | None, bool]],
    features: Features,
    bands: Sequence[str],
    *,
    band: Callable[[Any], str],
    key: Callable[[Any], str],
) -> dict[str, Any]:
    """One library's counts per band, for the table's `measured` and the report: its calls,
    false calls, truth events and those found, each feature value's false and true counts,
    and the classes seen (`key`)."""
    out: dict[str, Any] = {}
    for name in bands:
        rows = [(c, ok) for c, _, ok in calls if band(c) == name]
        found = [ok for b, ok in truth if b == name]
        out[name] = {
            "calls": len(rows),
            "false": sum(not ok for _, ok in rows),
            "truth": len(found),
            "found": sum(found),
            "features": {
                f: {
                    v: [
                        sum(feature_value(c, f) == v and not ok for c, ok in rows),
                        sum(feature_value(c, f) == v and ok for c, ok in rows),
                    ]
                    for v in values
                }
                for f, values in features.items()
            },
            "classes": sorted({key(c) for c, _ in rows}),
        }
    return out


def _share(s: Share) -> list[float]:
    return [s.point, s.lo, s.hi]


def write_classes(folder: Path, table: ClassTable) -> None:
    """Write `table` into `folder`, the staging folder `tables.write_table` swaps in whole."""
    document = {
        "format": FORMAT,
        "versions": table.versions,
        "features": {f: list(values) for f, values in table.features.items()},
        "bands": {
            band: {
                **asdict(terms),
                "false_calls": _share(terms.false_calls),
                "recall": _share(terms.recall),
            }
            for band, terms in table.bands.items()
        },
        "likelihoods": {
            band: {
                side: {
                    f: {v: _share(s) for v, s in values.items()}
                    for f, values in getattr(lk, side).items()
                }
                for side in ("false", "true")
            }
            for band, lk in table.likelihoods.items()
        },
        "measured": table.measured,
        "extra": table.extra,
    }
    text = json.dumps(document, indent=1, sort_keys=True)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / TABLE).write_text(text)


def write_class_table(
    lab: Path, name: str, table: ClassTable, *, report: str | None = None
) -> Path:
    """`table` (and `report`, when given) as the table `name`, swapped in whole under its
    versions by 3c's `tables.write_table`. Returns the folder."""

    def write(staging: Path) -> None:
        write_classes(staging, table)
        if report is not None:
            (staging / REPORT).write_text(report + "\n")

    return write_table(lab, name, write, table.versions)


def read_class_table[T: ClassTable](path: Path, kind: type[T] = ClassTable) -> T:
    """The table in the folder `path`, as `kind` (a `ClassTable` or a subclass);
    `ValueError` for another format."""
    document = json.loads((path / TABLE).read_text())
    if document.get("format") != FORMAT:
        raise ValueError(f"{path}: class table format {document.get('format')}")

    def sides(found: dict[str, Any], side: str) -> dict[str, dict[str, Share]]:
        return {
            f: {v: Share(*s) for v, s in values.items()}
            for f, values in found[side].items()
        }

    read = {
        "versions": dict(document["versions"]),
        "bands": {
            band: BandTerms(
                **{
                    **terms,
                    "false_calls": Share(*terms["false_calls"]),
                    "recall": Share(*terms["recall"]),
                }
            )
            for band, terms in document["bands"].items()
        },
        "likelihoods": {
            band: Likelihoods(false=sides(found, "false"), true=sides(found, "true"))
            for band, found in document["likelihoods"].items()
        },
        "measured": dict(document["measured"]),
        "features": {f: tuple(v) for f, v in document["features"].items()},
        "extra": dict(document["extra"]),
    }
    return kind(**{f.name: read[f.name] for f in fields(kind)})


def open_read[T](
    lab: Path, name: str, versions: Mapping[str, str], read: Callable[[Path], T]
) -> tuple[T | None, str | None]:
    """The lab's table `name`, read from its folder by `read`, when it was measured under
    `versions` (the run's own); otherwise None and why (`tables.open_table`'s reason, or the
    table unreadable)."""
    found = open_table(lab, name, versions=versions)
    if not isinstance(found, Current):
        return None, found.reason
    try:
        return read(found.folder), None
    except (OSError, ValueError, KeyError, TypeError) as e:
        return None, f"{name} table unreadable ({' '.join(str(e).split())}); {RUN}"


def open_class_table[T: ClassTable](
    lab: Path, name: str, versions: Mapping[str, str], kind: type[T] = ClassTable
) -> tuple[T | None, str | None]:
    """`open_read` of a class table, read as `kind`."""
    return open_read(lab, name, versions, lambda path: read_class_table(path, kind))


def _rate(found: tuple[float, float, float] | None) -> str:
    if found is None:
        return "not measured"
    return f"{found[0]:.3g} ({found[1]:.3g}-{found[2]:.3g})"


def report(table: ClassTable, parse: Callable[[str], Banded]) -> str:
    """The versions, then per band and library the counts and each feature value's false
    and true counts, the band's `h` and prior, and `p_wrong` with its interval for every
    class seen (each named as `summary`'s `key` names it, read back by `parse`)."""
    lines = ["versions measured under:"]
    lines += [f"  {key}: {value}" for key, value in sorted(table.versions.items())]
    for band, terms in table.bands.items():
        h = "none" if terms.h is None else f"{terms.h:.3g}"
        lines += ["", f"{band} band: h {h}, prior {terms.prior:.3g}"]
        seen: set[str] = set()
        for library, measured in table.measured.items():
            found = measured[band]
            lines.append(
                f"  {library}: {found['calls']} rare calls, {found['false']} false;"
                f" truth {found['found']} of {found['truth']} found"
            )
            for f, values in found["features"].items():
                cells = ", ".join(
                    f"{v} {counts[0]} false/{counts[1]} true"
                    for v, counts in values.items()
                )
                lines.append(f"    {f}: {cells}")
            seen |= set(found["classes"])
        for name in sorted(seen):
            lines.append(f"  {name}: p_wrong {_rate(table.p_wrong(parse(name)))}")
    return "\n".join(lines)
