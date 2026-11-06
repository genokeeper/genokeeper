"""Slice 3d's calibration, its pure part: the trio, HG002's truth, the gain truth's k-mer rule,
trio judging, spike-in recall, the CNVpytor decision, the paralog loci withheld, the large
bands pooled per genome, the versions the three tables (`cnv`, `insertions`, `paralogs`) are
measured under, and their report. The builders that run the tools are `calibrate.build`'s.

The depth and insertion tables are `classes.ClassTable`s over their own schemas
(`CNV_FEATURES`, `INSERTION_FEATURES`); the paralog table (`ParalogTable`) holds the loci
withheld and what each library gave. All three are written and opened through 3c's table
module (`calibrate.tables`), whose `open_table` is the one staleness check: its anchor is the
versions a run computes itself (`cnv_versions`, `insertion_versions`, `paralog_versions`: the
image names written in the code, the step versions, the resources' and the lab's files'
sha256, the thresholds and the FASTA's sha256), which a table measured under others cannot
move. ClinGen's lists are not among them: they decide what is shown, not how often a call is
false (amended spec). This module imports nothing from `pipeline`, so `pipeline` may import
it.
"""

import hashlib
import json
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from importlib.resources import files
from pathlib import Path
from typing import Any

from genokeeper.calibrate import classes, spikein
from genokeeper.calibrate.classes import FORMAT, FORMAT_KEY, OUTSIDE_VERSION, ClassTable
from genokeeper.calibrate.deletions import (
    DYSGU_STEP_VERSION,
    FEATURE_VALUES,
    class_of,
)
from genokeeper.calibrate.tables import REPORT, write_table
from genokeeper.checksum import sha256_file
from genokeeper.cnv.screen import (
    CNV_BANDS as DEPTH_BANDS,
)
from genokeeper.cnv.screen import (
    CNV_FEATURE_VALUES,
    LargeEvent,
    cnv_band,
)
from genokeeper.config import CnvThresholds, Grading, LpaThresholds
from genokeeper.containers import IMAGES
from genokeeper.deletions.exons import RESOURCE, SPANS
from genokeeper.deletions.files import sv_files
from genokeeper.deletions.screen import LABELS, CallClass, Deletion, _Index, reciprocal
from genokeeper.deletions.truvari import AGAINST_TRUTH, BETWEEN_CALLERS
from genokeeper.insertions import INSERTION_CALLERS, SCRAMBLE_MARGIN
from genokeeper.masking import masked_bed
from genokeeper.paralogs import LOCI, GauchianCall, wild_type
from genokeeper.regions import giab38_dir
from genokeeper.smn import parascopy_md5s
from genokeeper.tools import run_tool

# The calibration libraries (3c's `records.RECORDED`) and the genome each is; `pcr` is read
# at the CRAM level only.
TRIO = {"pcrfree": "HG002", "pcr": "HG002", "hg003": "HG003", "hg004": "HG004"}
_HG002 = tuple(k for k, sample in TRIO.items() if sample == "HG002")
_PARENTS = ("hg003", "hg004")
# The run's 3d steps: bump when a step's output changes. Here, not in `pipeline`, because
# the tables record them.
CNV_STEP_VERSION = "2"
CNVPYTOR_STEP_VERSION = "1"
INSURVEYOR_STEP_VERSION = "1"
SCRAMBLE_STEP_VERSION = "1"
INSERTIONS_UNION_STEP_VERSION = "1"
PARASCOPY_LOCI_STEP_VERSION = "1"
GAUCHIAN_STEP_VERSION = "1"
KILDA_STEP_VERSION = "1"
ALDY_CYP2A6_STEP_VERSION = "1"
DYSGU_INS_DUP_STEP_VERSION = "1"
PARASCOPY_DEPTH_STEP_VERSION = "1"
CNV_IMAGES = ("delly", "cnvpytor", "dysgu")
INSERTION_IMAGES = ("dysgu", "insurveyor", "scramble", "truvari")
PARALOG_IMAGES = ("parascopy", "gauchian", "kilda", "aldy")
# The depth classes' schema: `classify_cnv`'s features but `repeat`, which a depth event
# never has (the plan's Decisions, "Depth class features").
CNV_FEATURES: dict[str, tuple[str, ...]] = {
    f: values for f, values in CNV_FEATURE_VALUES.items() if f != "repeat"
}
# The measured depth bands; the two `*_large` bands are counted, never measured (Decisions,
# "Large band").
LARGE_BANDS = tuple(b for b in DEPTH_BANDS if b.endswith("_large"))
CNV_BANDS = tuple(b for b in DEPTH_BANDS if b not in LARGE_BANDS)
# The insertion classes' schema: 3b's, the evidence from the Truvari union or Scramble alone.
INSERTION_FEATURES: dict[str, tuple[str, ...]] = {
    "evidence": ("both", *INSERTION_CALLERS, *sorted(LABELS)),
    "repeat": FEATURE_VALUES["repeat"],
    "masked": FEATURE_VALUES["masked"],
}
INSERTION_BANDS = ("short", "long", "mei")
# The gain truth (Decisions, "Gain truth"): a v5.0q insertion of `GAIN_MIN_LENGTH` or more
# is a tandem duplication when `GAIN_CONTAINMENT` of its canonical `KMER`-mers occur in the
# reference within its length and `GAIN_PAD` of its position.
GAIN_MIN_LENGTH = 1000
GAIN_CONTAINMENT = 0.8
GAIN_PAD = 100
KMER = 31
TRUTH = "hg002_truth.tsv"
TRUTH_COLUMNS = (
    "locus", "sample", "haplotype_1", "haplotype_2", "total", "source", "method", "reason",
)  # fmt: skip
NO_ROW = "no row in the truth resource"
NO_CALL = "not covered: no call in calibration"
DISAGREE = "not covered: HG002's libraries disagree"
EXCEEDS = "not covered: HG002 exceeds its parents' copies"
PARALOGS = "paralogs"
PARALOG_FORMAT = 1
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _resource_sha256(*parts: str) -> str:
    data = files("genokeeper").joinpath("resources", *parts).read_bytes()
    return hashlib.sha256(data).hexdigest()


def _thresholds(section: Any) -> str:
    """A grading section as JSON, without `likely_below` (a wording rule no table reads)."""
    found = {k: v for k, v in asdict(section).items() if k != "likely_below"}
    return json.dumps(found, sort_keys=True)


# --- HG002's truth --------------------------------------------------------------------------


@dataclass(frozen=True)
class Truth:
    """A locus's copy number on one sample (`total`), or None and why there is none."""

    total: int | None
    reason: str | None


def read_truth(text: str) -> dict[tuple[str, str], Truth]:
    """The truth resource's rows, (locus, sample) -> its row: an empty `total` is no truth,
    with the row's `reason` (`ValueError` when it gives none)."""
    rows = [line.split("\t") for line in text.splitlines() if line and line[0] != "#"]
    if not rows or tuple(rows[0]) != TRUTH_COLUMNS:
        raise ValueError(f"truth resource: header is not {', '.join(TRUTH_COLUMNS)}")
    truth: dict[tuple[str, str], Truth] = {}
    for fields in rows[1:]:
        row = dict(zip(TRUTH_COLUMNS, fields + [""] * len(TRUTH_COLUMNS), strict=False))
        locus, sample, total, reason = (
            row[k].strip() for k in ("locus", "sample", "total", "reason")
        )
        if not total and not reason:
            raise ValueError(
                f"truth resource: {locus} {sample} has no total and no reason"
            )
        truth[locus, sample] = Truth(int(total), None) if total else Truth(None, reason)
    return truth


def load_truth() -> dict[tuple[str, str], Truth]:
    """The packaged truth (`resources/paralogs/hg002_truth.tsv`)."""
    return read_truth(
        files("genokeeper").joinpath("resources", "paralogs", TRUTH).read_text()
    )


def load_kiv2_truth() -> Truth:
    """HG002's KIV-2 count from the packaged truth."""
    return load_truth().get(("KIV2", "HG002"), Truth(None, NO_ROW))


# --- the gain truth -------------------------------------------------------------------------


def _kmers(sequence: str) -> set[str]:
    """The canonical `KMER`-mers of `sequence` holding only A, C, G and T."""
    s = sequence.upper()
    found: set[str] = set()
    for i in range(len(s) - KMER + 1):
        k = s[i : i + KMER]
        if k.strip("ACGT"):
            continue
        found.add(min(k, k.translate(_COMPLEMENT)[::-1]))
    return found


def _fetch(fasta: Path, regions: Sequence[str], work: Path) -> list[str]:
    """Each 1-based region's sequence, by `samtools faidx -r` natively."""
    listed = work / "regions.txt"
    listed.write_text("".join(f"{r}\n" for r in regions))
    text = run_tool(["samtools", "faidx", "-r", str(listed), str(fasta)])
    return ["".join(record.splitlines()[1:]) for record in text.split(">")[1:]]


def classify_gain_truth(
    ins: Sequence[tuple[Deletion, str]], fasta: Path
) -> list[Deletion]:
    """The truth insertions (each with its inserted sequence, the ALT without its anchor
    base) of `GAIN_MIN_LENGTH` or more whose canonical k-mers occur at `GAIN_CONTAINMENT`
    or more in the reference `[POS - L - GAIN_PAD, POS + L + GAIN_PAD]` (0-based
    `[pos - L - pad, pos + L + pad)`, read by `samtools faidx`): each a DUP `[POS, POS + L)`
    or `[POS - L, POS)`, the side holding more of those k-mers (after it on a tie)."""
    kept = [
        (d, seq)
        for d, seq in ins
        if d.svtype == "INS"
        and not d.lower_bound
        and d.svlen is not None
        and d.svlen >= GAIN_MIN_LENGTH
    ]
    if not kept:
        return []
    regions: list[str] = []
    for d, _ in kept:
        n = d.length
        regions += [
            f"{d.chrom}:{max(d.pos - n - GAIN_PAD, 0) + 1}-{d.pos}",
            f"{d.chrom}:{d.pos + 1}-{d.pos + n + GAIN_PAD}",
        ]
    with tempfile.TemporaryDirectory(prefix=".gain-truth-") as folder:
        windows = _fetch(fasta, regions, Path(folder))
    found: list[Deletion] = []
    for i, (d, seq) in enumerate(kept):
        kmers = _kmers(seq)
        left = kmers & _kmers(windows[2 * i])
        right = kmers & _kmers(windows[2 * i + 1])
        if not kmers or len(left | right) / len(kmers) < GAIN_CONTAINMENT:
            continue
        n = d.length
        start, end = (
            (d.pos, d.pos + n) if len(right) >= len(left) else (d.pos - n, d.pos)
        )
        found.append(
            replace(d, pos=start, end=end, svtype="DUP", svlen=None, lower_bound=False)
        )
    return found


# --- trio, spike-ins, CNVpytor --------------------------------------------------------------


def _span(d: Deletion) -> tuple[str, int, int]:
    return d.chrom, d.pos, d.end


def _matcher(
    events: Iterable[Deletion], min_overlap: float
) -> Callable[[Deletion], bool]:
    """Whether an event of `events` has the call's svtype at reciprocal overlap
    `min_overlap` or more."""
    index = _Index((e.chrom, e.pos, e.end, e) for e in events)

    def matches(d: Deletion) -> bool:
        return any(
            e.svtype == d.svtype and reciprocal(_span(d), _span(e)) >= min_overlap
            for e in index.hits(d.chrom, d.pos, d.end)
        )

    return matches


def trio_judged(
    hg002: Sequence[Deletion],
    parents: Sequence[Sequence[Deletion]],
    *,
    min_overlap: float,
    common: Callable[[Deletion], bool],
) -> list[tuple[Deletion, bool]]:
    """Each HG002 call outside the benchmark judged by the trio (Decisions, "False gains
    outside the benchmark"; C18): true when present (same svtype, reciprocal overlap
    `min_overlap`) in exactly one parent and not `common`; false when absent from both,
    present in both, or common."""
    present = [_matcher(p, min_overlap) for p in parents]
    return [
        (d, not common(d) and sum(match(d) for match in present) == 1) for d in hg002
    ]


def spike_recall(
    found: Sequence[Deletion],
    spikes: Sequence[Deletion],
    t: CnvThresholds,
    *,
    min_overlap: float,
) -> list[tuple[str, bool]]:
    """Each spike's band (`cnv_band`) and whether a found event (one that passed the
    screen) of its svtype matches it at reciprocal overlap `min_overlap`: the depth callers'
    recall only."""
    matches = _matcher(found, min_overlap)
    return [(cnv_band(s, t), matches(s)) for s in spikes]


def cnvpytor_kept(
    calls: Mapping[str, Sequence[tuple[CallClass, int, bool]]],
) -> tuple[bool, dict[str, Any]]:
    """Whether CNVpytor stays (Decisions, "CNVpytor stays"): pooled over libraries and
    bands, the false share of the calls both depth callers made (`evidence` `both`) below
    the false share of every call DELLY made (`both` or `dellycnv`), point estimates; with
    the counts, the record the `cnv` table keeps (`extra["cnvpytor"]`)."""
    rows = [(c.evidence, ok) for found in calls.values() for c, _, ok in found]

    def counts(evidence: tuple[str, ...]) -> dict[str, int]:
        judged = [ok for e, ok in rows if e in evidence]
        return {"calls": len(judged), "false": sum(not ok for ok in judged)}

    both, delly = counts(("both",)), counts(("both", "dellycnv"))
    kept = (
        both["calls"] > 0
        and delly["calls"] > 0
        and both["false"] / both["calls"] < delly["false"] / delly["calls"]
    )
    return kept, {
        "decision": "kept" if kept else "dropped",
        "both": both,
        "delly": delly,
    }


# --- large events ---------------------------------------------------------------------------


def pool_large(
    events_by_library: Mapping[str, Sequence[LargeEvent]], *, min_overlap: float
) -> dict[str, dict[str, int]]:
    """Per genome (`TRIO`) and large band, its large events (`screen_cnv`'s): HG002's two
    libraries are one genome, so an event matching one already counted on that genome (the
    same svtype at reciprocal overlap `min_overlap`) counts once."""
    pooled = {s: dict.fromkeys(LARGE_BANDS, 0) for s in dict.fromkeys(TRIO.values())}
    counted: dict[str, list[Deletion]] = {}
    for library, events in events_by_library.items():
        sample = TRIO[library]
        seen = counted.setdefault(sample, [])
        for e in events:
            if _matcher(seen, min_overlap)(e.event):
                continue
            seen.append(e.event)
            pooled[sample][e.call_class.band] += 1
    return pooled


def large_withheld(pooled: Mapping[str, Mapping[str, int]]) -> dict[str, str]:
    """The large bands with a call on two or more of the genomes, each with its reason (the
    `cnv` table's `extra["large"]`; Decisions, "Large band": three genomes measure no
    rate)."""
    withheld: dict[str, str] = {}
    for band in LARGE_BANDS:
        n = sum(counts[band] > 0 for counts in pooled.values())
        if n >= 2:
            direction = band.split("_")[0]
            withheld[band] = (
                f"not covered: large {direction} calls on {n} of {len(pooled)}"
                " calibration genomes"
            )
    return withheld


# --- paralogs -------------------------------------------------------------------------------


def _no_truth(row: Truth | None) -> str:
    reason = NO_ROW if row is None else row.reason
    return f"not calibrated: no HG002 truth: {reason}"


def _copies_withheld(
    copies: Mapping[str, int | None],
    truths: Mapping[str, Truth],
) -> str | None:
    """A locus's reason from its copies per library: a call missing, HG002's libraries
    apart, a library off its sample's truth (`truths`, sample -> row), HG002 above its
    parents' sum (Decisions, "Trio rule")."""
    if any(cn is None for cn in copies.values()):
        return NO_CALL
    if len({copies[k] for k in _HG002}) > 1:
        return DISAGREE
    for library, cn in copies.items():
        total = truths.get(TRIO[library], Truth(None, None)).total
        if total is not None and cn != total:
            return f"not covered: {library} gave {cn}, truth {total}"
    hg002, *_ = (copies[k] for k in _HG002)
    parents = [copies[k] for k in _PARENTS]
    if (
        hg002 is not None
        and None not in parents
        and hg002 > sum(p or 0 for p in parents)
    ):
        return EXCEEDS
    return None


def paralog_measure(
    per_sample: Mapping[str, Mapping[str, int | None]],
    truth: Mapping[tuple[str, str], Truth],
    gauchian: Mapping[str, GauchianCall | None],
    kiv2: Mapping[str, float | None],
    kiv2_truth: float | None,
    t: LpaThresholds,
) -> tuple[dict[str, str], dict[str, Any]]:
    """The loci of `paralogs.LOCI` and `LPA` the calibration withholds, each with its
    reason, and what each library of `TRIO` gave (its copies per locus, KILDA's count,
    whether Gauchian found GBA1 wild type). `per_sample` is library -> locus -> Parascopy's
    copy number (a library or locus absent: no call).

    A locus needs HG002's truth row (else "not calibrated"); GBA1 has none: its HG002
    truth is wild type (R2), so Gauchian must call both HG002 libraries wild type. LPA needs
    `kiv2_truth` (else the KIV2 row's reason) and each HG002 library's KILDA count within
    `t.max_kiv2_error` of it, as a fraction of it."""
    withheld: dict[str, str] = {}
    for locus in LOCI:
        copies = {k: per_sample.get(k, {}).get(locus) for k in TRIO}
        if locus == "GBA1":
            calls = [gauchian.get(k) for k in _HG002]
            reason = NO_CALL if None in calls else _copies_withheld(copies, {})
            off = [
                k
                for k, g in zip(_HG002, calls, strict=True)
                if g is not None and not wild_type(g)
            ]
            if reason in (None, EXCEEDS) and off:
                reason = f"not covered: GBA1 not wild type on HG002 {off[0]}"
        else:
            rows = {s: row for (name, s), row in truth.items() if name == locus}
            row = rows.get("HG002")
            reason = (
                _no_truth(row)
                if row is None or row.total is None
                else _copies_withheld(copies, rows)
            )
        if reason is not None:
            withheld[locus] = reason
    if kiv2_truth is None:
        withheld["LPA"] = _no_truth(truth.get(("KIV2", "HG002")))
    else:
        counts = [kiv2.get(k) for k in _HG002]
        if any(c is None for c in counts):
            withheld["LPA"] = NO_CALL
        else:
            d = max(abs((c or 0.0) - kiv2_truth) for c in counts)
            if d > t.max_kiv2_error * kiv2_truth:
                withheld["LPA"] = (
                    f"not covered: KILDA off HG002's count by {d:.1f}"
                    f" ({100 * d / kiv2_truth:.0f}%)"
                )
    measured: dict[str, Any] = {
        library: {
            "copies": {locus: per_sample.get(library, {}).get(locus) for locus in LOCI},
            "kiv2": kiv2.get(library),
            "gba1_wild_type": None
            if (g := gauchian.get(library)) is None
            else wild_type(g),
        }
        for library in TRIO
    }
    return withheld, measured


@dataclass(frozen=True)
class ParalogTable:
    """The versions measured under, the loci withheld (locus -> reason) and what each
    library gave (`paralog_measure`)."""

    versions: dict[str, str]
    withheld: dict[str, str]
    measured: dict[str, Any]


def write_paralog_table(
    lab: Path, table: ParalogTable, *, report: str | None = None
) -> Path:
    """`table` (and `report`) as the `paralogs` table, through `tables.write_table`."""

    def write(staging: Path) -> None:
        document = {"format": PARALOG_FORMAT, **asdict(table)}
        (staging / classes.TABLE).write_text(
            json.dumps(document, indent=1, sort_keys=True)
        )
        if report is not None:
            (staging / REPORT).write_text(report + "\n")

    return write_table(lab, PARALOGS, write, table.versions)


def read_paralog_table(path: Path) -> ParalogTable:
    """The paralog table in the folder `path`; `ValueError` for another format."""
    document = json.loads((path / classes.TABLE).read_text())
    if document.get("format") != PARALOG_FORMAT:
        raise ValueError(f"{path}: paralog table format {document.get('format')}")
    return ParalogTable(
        versions=dict(document["versions"]),
        withheld=dict(document["withheld"]),
        measured=dict(document["measured"]),
    )


def open_paralog_table(
    lab: Path, versions: Mapping[str, str]
) -> tuple[ParalogTable | None, str | None]:
    """The lab's paralog table when it was measured under `versions` (the run's own,
    `paralog_versions`); otherwise None and why (`classes.open_read`)."""
    return classes.open_read(lab, PARALOGS, versions, read_paralog_table)


# --- versions -------------------------------------------------------------------------------


def _models() -> str:
    return ",".join(sorted({locus.model for locus in LOCI.values()}))


def _mask(lab: Path, grading: Grading) -> dict[str, str]:
    giab = giab38_dir(lab)
    return {
        "segdup": sha256_file(giab / grading.strata["segdup"]),
        "low_mappability": sha256_file(giab / grading.strata["low_mappability"]),
    }


def cnv_versions(lab: Path, fasta_sha: str, grading: Grading) -> dict[str, str]:
    """The versions the `cnv` table is measured under: the depth callers' and Dysgu's images
    (its DUP records are the gains' junctions), the exon and span resources' sha256, the
    mask's two BEDs', 3c's masked intervals' and gnomAD-SV's, Parascopy's data md5s and the
    models whose regions the screen leaves to their callers, `[cnv]`, `[deletions]` (its
    margin, gnomAD-SV rule and overlap) and `[par]` (a male's loss or gain) as JSON, the
    spike-in rule's constants, the step versions, the FASTA's sha256, the table format and
    the form of its `extra` `ClassTable.cautions` reads. `FileNotFoundError` when a file is
    absent."""
    return {
        **{name: IMAGES[name] for name in CNV_IMAGES},
        "exons": _resource_sha256("deletions", RESOURCE),
        "spans": _resource_sha256("paralogs", SPANS),
        **_mask(lab, grading),
        "masked_intervals": sha256_file(masked_bed(lab)),
        "gnomad_sv": sha256_file(sv_files(lab).gnomad),
        "parascopy_data": parascopy_md5s(lab),
        "parascopy_models": _models(),
        "cnv": _thresholds(grading.cnv),
        "deletions": _thresholds(grading.deletions),
        "par": json.dumps(grading.par, sort_keys=True),
        "spikein": json.dumps(
            [spikein.SEED, spikein.N_PER, spikein.EXON_COUNTS, spikein.FLANK_CAP]
            + [spikein.MIN_SPIKE, spikein.SPACING]
        ),
        "cnv_step": CNV_STEP_VERSION,
        "cnvpytor_step": CNVPYTOR_STEP_VERSION,
        "dysgu_step": DYSGU_STEP_VERSION,
        "dysgu_ins_dup_step": DYSGU_INS_DUP_STEP_VERSION,
        "fasta": fasta_sha,
        FORMAT_KEY: str(FORMAT),
        "extra": OUTSIDE_VERSION,
    }


def insertion_versions(lab: Path, fasta_sha: str, grading: Grading) -> dict[str, str]:
    """The versions the `insertions` table is measured under: the four images, the exon
    resource's sha256, the mask's two BEDs', the repeats' and gnomAD-SV's, `[insertions]` as
    JSON, Scramble's exon margin, both Truvari option sets, the step versions, the FASTA's
    sha256, the table format and the form of its `extra` `ClassTable.cautions` reads.
    `FileNotFoundError` when a file is absent."""
    giab = giab38_dir(lab)
    return {
        **{name: IMAGES[name] for name in INSERTION_IMAGES},
        "exons": _resource_sha256("deletions", RESOURCE),
        **_mask(lab, grading),
        "repeats": sha256_file(giab / grading.strata["low_complexity"]),
        "gnomad_sv": sha256_file(sv_files(lab).gnomad),
        "insertions": _thresholds(grading.insertions),
        "scramble_margin": str(SCRAMBLE_MARGIN),
        "between_callers": " ".join(BETWEEN_CALLERS),
        "against_truth": " ".join(AGAINST_TRUTH),
        "dysgu_step": DYSGU_STEP_VERSION,
        "dysgu_ins_dup_step": DYSGU_INS_DUP_STEP_VERSION,
        "insurveyor_step": INSURVEYOR_STEP_VERSION,
        "scramble_step": SCRAMBLE_STEP_VERSION,
        "insertions_union_step": INSERTIONS_UNION_STEP_VERSION,
        "fasta": fasta_sha,
        FORMAT_KEY: str(FORMAT),
        "extra": OUTSIDE_VERSION,
    }


def paralog_versions(lab: Path, fasta_sha: str, grading: Grading) -> dict[str, str]:
    """The versions the `paralogs` table is measured under: the four images, the span
    resource's and the truth resource's sha256, Parascopy's data md5s and models,
    `[paralogs]` and `[lpa]` as JSON, the step versions, the FASTA's sha256 and the table
    format. `FileNotFoundError` when a file is absent."""
    return {
        **{name: IMAGES[name] for name in PARALOG_IMAGES},
        "spans": _resource_sha256("paralogs", SPANS),
        "truth": _resource_sha256("paralogs", TRUTH),
        "parascopy_data": parascopy_md5s(lab),
        "parascopy_models": _models(),
        "paralogs": _thresholds(grading.paralogs),
        "lpa": _thresholds(grading.lpa),
        "parascopy_depth_step": PARASCOPY_DEPTH_STEP_VERSION,
        "parascopy_loci_step": PARASCOPY_LOCI_STEP_VERSION,
        "gauchian_step": GAUCHIAN_STEP_VERSION,
        "kilda_step": KILDA_STEP_VERSION,
        "aldy_cyp2a6_step": ALDY_CYP2A6_STEP_VERSION,
        "fasta": fasta_sha,
        FORMAT_KEY: str(PARALOG_FORMAT),
    }


# --- the report -----------------------------------------------------------------------------


def _cnv_lines(table: ClassTable) -> list[str]:
    lines = ["cnv table (recall from spike-ins: the depth callers only)"]
    lines.append(classes.report(table, class_of))
    record = table.extra.get("cnvpytor")
    if record is not None:
        both, delly = record["both"], record["delly"]
        lines.append(
            f"CNVpytor {record['decision']}: calls both depth callers made"
            f" {both['false']} false of {both['calls']}; DELLY's calls {delly['false']}"
            f" false of {delly['calls']}"
        )
    large = table.extra.get("large")
    if large is not None:
        for sample, counts in large["calls"].items():
            cells = ", ".join(f"{band} {n}" for band, n in counts.items())
            lines.append(f"large calls on {sample}: {cells}")
        for band, reason in large["withheld"].items():
            lines.append(f"{band}: {reason}")
    return lines


def _paralog_lines(table: ParalogTable) -> list[str]:
    lines = ["paralogs table", "versions measured under:"]
    lines += [f"  {key}: {value}" for key, value in sorted(table.versions.items())]
    for library, found in table.measured.items():
        copies = ", ".join(f"{locus} {cn}" for locus, cn in found["copies"].items())
        lines.append(
            f"  {library}: {copies}; KIV-2 {found['kiv2']};"
            f" GBA1 wild type {found['gba1_wild_type']}"
        )
    lines += [
        f"  {locus}: {reason}" for locus, reason in sorted(table.withheld.items())
    ]
    if not table.withheld:
        lines.append("  no locus withheld")
    return lines


def report(
    cnv: ClassTable | None, insertions: ClassTable | None, paralogs: ParalogTable | None
) -> str:
    """The tables given: each table's counts per band and library, the CNVpytor decision
    with its counts, the large calls per genome and the bands withheld, every class's
    `p_wrong` (`classes.report`) and each withheld paralog locus."""
    parts: list[str] = []
    if cnv is not None:
        parts.append("\n".join(_cnv_lines(cnv)))
    if insertions is not None:
        parts.append("insertions table\n" + classes.report(insertions, class_of))
    if paralogs is not None:
        parts.append("\n".join(_paralog_lines(paralogs)))
    return "\n\n".join(parts)
