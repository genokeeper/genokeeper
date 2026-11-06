"""The calibration table: genotype error per cell, measured on a GIAB sample, with the
versions of the pipeline it was measured on.

A cell is a row type (a variant record, a RefCall, a reference block) x the claim's allele
class (SNV, indel) x region (easy, or difficult: any GIAB v3.6 stratum the grading flags)
x GQ bin (0-9, 10-19, ..., >= 50). Each cell counts errors over units (variant records, or
reference positions for RefCalls and blocks) for two kinds of error: a wrong genotype, and a
wrong allele (zygosity ignored). The kinds coincide for reference rows: a truth variant at a
reference call is both. Depth is kept per cell as a covariate (mean DP, or MIN_DP per
position).

Within each (kind, row type, allele class, region), the rates are made non-increasing in GQ
by pool-adjacent-violators: a bin measured worse than the one below it is pooled with it,
and the pooled rate and its Jeffreys 95% interval go to every bin of the pool, while each
cell keeps its own counts. The interval is taken at the pool's effective size, its counts
divided by their between-chromosome dispersion (`dispersion`), so a rate that varies
between chromosomes more than sampling explains is held less precise; the rate is not
changed by it. A cell with no units is absent: `rate_for` raises `KeyError`.

The table is a folder: `cells.tsv` and `manifest.json`, which records the seven versions that
invalidate it (`VERSION_KEYS`), the sample's SNP and INDEL F1 (`f1`, which Task 21's stop
rule reads), and the cells file's sha256. A table whose manifest lacks a
version, or whose cells file does not hash to the manifest's sha256, is refused on read.
The versions' anchor comes from the manifest the measurement writes beside its cells, which
a table measured without them (by older code, or assembled by hand) cannot pass for a
complete one because the constructor, and so `read`, refuses any of the seven missing; the
run compares them with its own (stage 2b, Task 20).
"""

import hashlib
import json
import math
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path

GQ_BINS = (0, 10, 20, 30, 40, 50)  # lower edges; the last is >= 50
ROW_TYPES = ("variant", "RefCall", "ref_block")
ALLELE_CLASSES = ("snv", "indel")
REGIONS = ("easy", "difficult")
GENOTYPE = "genotype"
ALLELE = "allele"
KINDS = (GENOTYPE, ALLELE)
# What a table was measured on; a run whose pipeline differs in any of them cannot use it.
VERSION_KEYS = (
    "bwa",
    "samtools",
    "deepvariant_image",
    "checkpoint_md5",
    "make_examples_flags",
    "reference_md5",
    "correct_copies",
)
FORMAT = 1
CELLS = "cells.tsv"
MANIFEST = "manifest.json"
_COLUMNS = (
    "kind",
    "row_type",
    "allele_class",
    "region",
    "gq_bin",
    "errors",
    "total",
    "rate",
    "lo",
    "hi",
    "depth",
)


def gq_bin(gq: int) -> int:
    """The lower edge of the GQ bin `gq` falls in."""
    if gq < 0:
        raise ValueError(f"GQ {gq} is negative")
    return max(edge for edge in GQ_BINS if edge <= gq)


@dataclass(frozen=True, order=True)
class CellKey:
    row_type: str  # "variant" | "RefCall" | "ref_block"
    allele_class: str  # "snv" | "indel"
    region: str  # "easy" | "difficult"
    gq_bin: int

    def __post_init__(self) -> None:
        if (
            self.row_type not in ROW_TYPES
            or self.allele_class not in ALLELE_CLASSES
            or self.region not in REGIONS
            or self.gq_bin not in GQ_BINS
        ):
            raise ValueError(f"not a calibration cell: {self}")


@dataclass(frozen=True)
class CellRate:
    """A cell's own counts, the rate of the pool it falls in with that pool's Jeffreys 95%
    interval, and the cell's mean depth."""

    errors: int
    total: int
    rate: float
    lo: float
    hi: float  # Jeffreys 95%
    depth: float | None = None


@dataclass
class Counts:
    """Errors over units of one cell, with the depth summed over its units."""

    errors: int = 0
    total: int = 0
    depth_sum: float = 0.0

    def add(self, errors: int, units: int, depth: float) -> None:
        self.errors += errors
        self.total += units
        self.depth_sum += depth * units

    def merge(self, other: Counts) -> None:
        self.errors += other.errors
        self.total += other.total
        self.depth_sum += other.depth_sum


# --- the Jeffreys interval ---------------------------------------------------------------


def _betacf(a: float, b: float, x: float) -> float:
    """The continued fraction of the incomplete beta (modified Lentz; Numerical Recipes
    6.4). It converges in about sqrt(max(a, b)) terms where `betainc` uses it."""
    tiny = 1e-300
    c, d = 1.0, 1.0 - (a + b) * x / (a + 1.0)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 10_000_000):
        m2 = 2 * m
        for aa in (
            m * (b - m) * x / ((a - 1.0 + m2) * (a + m2)),
            -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1.0 + m2)),
        ):
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            h *= d * c
        if abs(d * c - 1.0) < 1e-15:
            return h
    raise ArithmeticError(f"incomplete beta did not converge for a={a}, b={b}, x={x}")


def _log_beta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def betainc(a: float, b: float, x: float) -> float:
    """The regularised incomplete beta I_x(a, b): the Beta(a, b) distribution function."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(a * math.log(x) + b * math.log1p(-x) - _log_beta(a, b))
    below = x < (a + 1.0) / (a + b + 2.0)
    if front == 0.0:  # far in a tail: the fraction would not change the answer
        return 0.0 if below else 1.0
    if below:
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _beta_ppf(p: float, a: float, b: float) -> float:
    """The Beta(a, b) quantile `p`: Newton steps on `betainc`, kept inside a bisection
    bracket, from the normal approximation."""
    mean = a / (a + b)
    sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1.0)))
    z = math.copysign(1.959963984540054, p - 0.5) if p != 0.5 else 0.0
    x = min(max(mean + z * sd, 1e-300), 1.0 - 1e-16)
    lo, hi = 0.0, 1.0
    for _ in range(200):
        f = betainc(a, b, x) - p
        if f > 0:
            hi = x
        else:
            lo = x
        if abs(f) < 1e-13 or hi - lo < 1e-15 * max(x, 1e-300):
            return x
        log_pdf = (a - 1.0) * math.log(x) + (b - 1.0) * math.log1p(-x) - _log_beta(a, b)
        pdf = math.exp(log_pdf) if log_pdf < 700 else math.inf
        step = x - f / pdf if pdf > 0 else math.nan
        x = step if lo < step < hi else (lo + hi) / 2.0
    return x


def jeffreys(errors: float, total: float) -> tuple[float, float]:
    """The Jeffreys 95% interval of `errors` in `total`: the 2.5% and 97.5% quantiles of
    Beta(errors + 1/2, total - errors + 1/2), with 0 below no errors and 1 above all.
    Effective (non-integer) counts are taken as they are."""
    if not 0 <= errors <= total or total == 0:
        raise ValueError(f"no interval for {errors} errors in {total}")
    a, b = errors + 0.5, total - errors + 0.5
    lo = 0.0 if errors == 0 else _beta_ppf(0.025, a, b)
    hi = 1.0 if errors == total else _beta_ppf(0.975, a, b)
    return lo, hi


def beta_binomial_interval(n: int, a: float, b: float) -> tuple[int, int]:
    """The central 95% interval of the errors in `n` units whose rate is Beta(a, b) (the
    beta-binomial): its 2.5% and 97.5% quantiles. Its terms are walked out from the mean
    by the ratio of neighbours, pmf(k + 1) / pmf(k) = (n - k)(k + a) / ((k + 1)(n - k -
    1 + b)), until they fall below 1e-16 of the mass so far: no factorial is formed, so
    any `n` works, and the walk spans about ten standard deviations each way."""
    start = min(n, round(n * a / (a + b)))
    weights, total = {start: 1.0}, 1.0
    for step in (-1, 1):
        k, w = start, 1.0
        while 0 <= k + step <= n:
            j = k if step == 1 else k - 1
            ratio = (n - j) * (j + a) / ((j + 1) * (n - j - 1 + b))
            w, k = (w * ratio if step == 1 else w / ratio), k + step
            if w < 1e-16 * total:
                break
            weights[k] = w
            total += w
    seen, lo = 0.0, -1
    for k in sorted(weights):
        seen += weights[k] / total
        lo = k if lo < 0 and seen >= 0.025 else lo
        if seen >= 0.975:
            return lo, k
    return lo, n


def dispersion(per_chrom: Iterable[tuple[int, int]]) -> float:
    """The extra-binomial dispersion of a cell's (errors, total) per chromosome: Pearson's
    X^2 against their pooled rate over its degrees of freedom (the quasi-binomial moment
    estimator, McCullagh and Nelder 4.5), never below 1, the binomial. It fits because a
    cell's units are counted the same way on every chromosome, while each chromosome's
    mix of sequence moves the cell's true rate: X^2 / df estimates by how much the
    variance exceeds sampling's, and dividing the counts by it gives the effective sample
    size (Kish's design effect) that the Jeffreys and beta-binomial intervals then use.
    A cell measured on fewer than two chromosomes, or with no errors or nothing but
    errors, shows no spread to measure: 1."""
    found = [(e, n) for e, n in per_chrom if n > 0]
    errors, total = sum(e for e, _ in found), sum(n for _, n in found)
    if len(found) < 2 or errors in (0, total):
        return 1.0
    p = errors / total
    x2 = sum((e - n * p) ** 2 / (n * p * (1 - p)) for e, n in found)
    return max(1.0, x2 / (len(found) - 1))


# --- monotone rates ----------------------------------------------------------------------


def _pooled(cells: list[tuple[int, Counts]]) -> list[tuple[int, int, list[int]]]:
    """Pool-adjacent-violators over `cells` sorted by GQ bin: (errors, total, GQ bins) of
    the pool each cell falls in, so that the pooled rates never rise with GQ."""
    pools: list[tuple[int, int, list[int]]] = []
    for gq, counts in cells:
        pools.append((counts.errors, counts.total, [gq]))
        while (
            len(pools) > 1 and pools[-1][0] * pools[-2][1] > pools[-2][0] * pools[-1][1]
        ):
            errors, total, bins = pools.pop()
            below = pools.pop()
            pools.append((below[0] + errors, below[1] + total, below[2] + bins))
    return [pool for pool in pools for _ in pool[2]]


# --- the table ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CalibrationTable:
    """Each error kind's cells (`GENOTYPE`, `ALLELE`) and the seven versions measured on."""

    cells: dict[str, dict[CellKey, CellRate]]
    versions: dict[str, str] = field(default_factory=dict)
    f1: dict[str, float] = field(default_factory=dict)  # "snp", "indel"

    def __post_init__(self) -> None:
        absent = [key for key in VERSION_KEYS if not self.versions.get(key)]
        if absent:
            raise ValueError(
                f"a calibration table needs the versions {', '.join(absent)}"
            )
        extra = sorted(set(self.versions) - set(VERSION_KEYS))
        if extra:
            raise ValueError(f"unknown calibration versions {', '.join(extra)}")
        unknown = sorted(set(self.cells) - set(KINDS))
        if unknown:
            raise ValueError(f"unknown error kinds {', '.join(unknown)}")

    @classmethod
    def from_counts(
        cls,
        counts: Mapping[tuple[str, CellKey], Counts],
        versions: Mapping[str, str],
        f1: Mapping[str, float] | None = None,
        *,
        by_chrom: Mapping[str, Mapping[tuple[str, CellKey], Counts]] | None = None,
    ) -> CalibrationTable:
        """The table of `counts` (kind, cell -> counts): monotone rates with intervals.
        `by_chrom` holds the per-chromosome counts `counts` sums, whose spread widens
        each pool's interval (`dispersion`); without it the interval is binomial."""
        chroms = list((by_chrom or {}).values())
        groups: dict[tuple[str, str, str, str], list[tuple[int, Counts]]] = {}
        for (kind, key), found in counts.items():
            if found.total > 0:
                group = (kind, key.row_type, key.allele_class, key.region)
                groups.setdefault(group, []).append((key.gq_bin, found))
        cells: dict[str, dict[CellKey, CellRate]] = {kind: {} for kind in KINDS}
        for (kind, row_type, allele_class, region), found in sorted(groups.items()):
            found.sort(key=lambda item: item[0])
            for (gq, own), (errors, total, bins) in zip(
                found, _pooled(found), strict=True
            ):
                pool = [
                    (kind, CellKey(row_type, allele_class, region, g)) for g in bins
                ]
                phi = dispersion(
                    (
                        sum(c[k].errors for k in pool if k in c),
                        sum(c[k].total for k in pool if k in c),
                    )
                    for c in chroms
                )
                lo, hi = jeffreys(errors / phi, total / phi)
                cells[kind][CellKey(row_type, allele_class, region, gq)] = CellRate(
                    errors=own.errors,
                    total=own.total,
                    rate=errors / total,
                    lo=lo,
                    hi=hi,
                    depth=own.depth_sum / own.total,
                )
        return cls(cells=cells, versions=dict(versions), f1=dict(f1 or {}))

    def rate_for(self, key: CellKey, *, errors: str = GENOTYPE) -> CellRate:
        """The cell's genotype error (`p_wrong`), or with `errors=ALLELE` its allele error;
        `KeyError` for a cell nothing was measured in."""
        try:
            return self.cells[errors][key]
        except KeyError:
            raise KeyError(f"no {errors} error measured for {key}") from None

    def measured(self, *, errors: str = GENOTYPE) -> list[CellKey]:
        """The measured cells of one error kind, sorted."""
        return sorted(self.cells.get(errors, {}))

    def _rows(self) -> Iterator[str]:
        yield "\t".join(_COLUMNS) + "\n"
        for kind in KINDS:
            for key, rate in sorted(self.cells.get(kind, {}).items()):
                values = (
                    kind,
                    key.row_type,
                    key.allele_class,
                    key.region,
                    str(key.gq_bin),
                    str(rate.errors),
                    str(rate.total),
                    repr(rate.rate),
                    repr(rate.lo),
                    repr(rate.hi),
                    "" if rate.depth is None else repr(rate.depth),
                )
                yield "\t".join(values) + "\n"

    def write(self, folder: Path) -> None:
        """Write the table's cells and manifest into `folder`; `tables.write_table` swaps
        the folder in whole."""
        text = "".join(self._rows())
        manifest = {
            "format": FORMAT,
            "versions": self.versions,
            "cells": CELLS,
            "cells_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "f1": self.f1,
        }
        folder.mkdir(parents=True, exist_ok=True)
        (folder / CELLS).write_text(text)
        (folder / MANIFEST).write_text(json.dumps(manifest, indent=1, sort_keys=True))

    @classmethod
    def read(cls, folder: Path) -> CalibrationTable:
        """The table in `folder`; `ValueError` for another format, a cells file whose
        sha256 is not the manifest's, or versions missing.

        The check's anchor comes from the manifest's sha256 of the cells file, written
        with it, which an edit of the cells file cannot move because the manifest is not
        rewritten with it; a deliberate edit of both is the user's own decision."""
        manifest = json.loads((folder / MANIFEST).read_text())
        if manifest.get("format") != FORMAT:
            raise ValueError(
                f"{folder}: calibration table format {manifest.get('format')}"
            )
        data = (folder / CELLS).read_bytes()
        if hashlib.sha256(data).hexdigest() != manifest.get("cells_sha256"):
            raise ValueError(f"{folder}: {CELLS} does not match the manifest's sha256")
        cells: dict[str, dict[CellKey, CellRate]] = {kind: {} for kind in KINDS}
        lines = data.decode().splitlines()
        if not lines or tuple(lines[0].split("\t")) != _COLUMNS:
            raise ValueError(f"{folder}: {CELLS} lacks its header")
        for line in lines[1:]:
            kind, row_type, allele_class, region, gq, *numbers = line.split("\t")
            errors, total, rate, lo, hi, depth = numbers
            cells.setdefault(kind, {})[
                CellKey(row_type, allele_class, region, int(gq))
            ] = CellRate(
                errors=int(errors),
                total=int(total),
                rate=float(rate),
                lo=float(lo),
                hi=float(hi),
                depth=float(depth) if depth else None,
            )
        versions = dict(manifest.get("versions", {}))
        return cls(cells=cells, versions=versions, f1=dict(manifest.get("f1", {})))
