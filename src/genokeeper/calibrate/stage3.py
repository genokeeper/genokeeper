"""Slice 3a's calibration: HG002's two libraries (`pcrfree`, `pcr`) through the run's three
stage-3 steps, into the table every run grades stage 3 with (spec, slice 3a "Calibration").

HG002 is healthy, so a repeat call on either library the run would not grade normal (an
interval reaching past the normal range, `sources.repeats.graded_normal`) is a false call and
withholds its locus, as does a call the run would leave unassessed for low depth; an allele shorter than the read whose size is off GIAB's HG002
v5.0q truth (`truth_repeats`) by more than `MAX_SIZE_ERROR` repeats caps its locus at low; SMN
is withheld unless both callers give SMN1 = 2 and SMN2 = 2 on both libraries, and the Y
haplogroup unless both libraries derive L816 with the same predicted haplogroup (`measure`). Each library's C9orf72 flank ratio
is kept, so a run can tell which library its own looks like. Slice 3e adds the trio checks on
HG002, HG003 and HG004 (`genokeeper.calibrate.repeat_rules`), kept under `measured[RULES]`.

The table (`Stage3Table`, `<lab>/calibration/stage3/`) records the versions it was measured
under (`stage3_versions`). A run computes its own and, when the table is missing or was made
with others, runs no stage-3 tool and states each gene unassessed with the reason
`open_stage3` gives. The version check's anchor is the run's own computation (the image names
written in the code, the catalog resource's and the FASTA's sha256, the Parascopy data's md5s
on disk, the grading's `[repeats]`, `[repeat_loci]` but its wording fields and `[exclusion]`,
the trio checks' constants), which a stale table cannot move, because the run never reads its versions from the table. It stops an honest
mistake: a rebuilt image, another reference or an edited threshold graded with a calibration
it was never measured with.
"""

import json
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pysam

from genokeeper.calibrate.tables import RUN, Current, folder, open_table
from genokeeper.config import Grading
from genokeeper.containers import IMAGES
from genokeeper.repeats.expansionhunter import catalog_sha256
from genokeeper.repeats.loci import RepeatCall, RepeatLocus, classify
from genokeeper.repeats.rfc1 import RFC1_LOCUS
from genokeeper.smn import parascopy_md5s
from genokeeper.sources.repeats import CFTR_GENE, graded_normal
from genokeeper.store.db import Level
from genokeeper.yhaplogroup import derived_markers

if TYPE_CHECKING:
    from genokeeper.pipeline import Stage3Results

STAGE3_IMAGES = ("expansionhunter", "reviewer", "smncaller", "parascopy", "yleaf")
TABLE = "table.json"
# 3 (slice 3e, after its final review): no `homozygous_ok`, so a table the false-normal
# trials wrote (format 2), or slice 3a's (format 1), is not used. 4 (BR-008): a low-depth
# locus withheld, a flanking SNV no longer voiding a truth, and each allele's interval and
# reads recorded, so a table measured without them (format 3) is not used.
FORMAT = 4
FALSE_CALL = "not covered: false call in calibration"
NO_CALL = "not covered: no call in calibration"
LOW_DEPTH = "not covered: low depth in calibration"
SMN_WITHHELD = "not covered: SMN1 and SMN2 not 2 and 2 from both callers in calibration"
YLEAF_WITHHELD = (
    "not covered: Yleaf did not give one haplogroup deriving L816 on both libraries"
    " in calibration"
)
# The `[repeat_loci]` fields no measurement reads (the statement's wording and validity):
# an edit to one leaves the table current.
UNMEASURED = ("note", "topic", "evidence", "validity")
# The keys `withheld` gives SMN and the Y haplogroup under; a repeat locus's is its gene.
SMN_KEY = "SMN1"
YLEAF_KEY = "yleaf"
# The `measured` key holding `repeat_rules.measure_repeat_rules`'s measurements.
RULES = "rules"
# The spec's stop condition: a size off by more than this caps the locus at low.
MAX_SIZE_ERROR = 2
# HG002's published Y lineage is J-L816 (Rhie et al., Nature, August 2023).
L816 = "L816"


# GIAB's HG002 v5.0q GRCh38 small-variant truth, as the bootstrap (step_hg002) names it.
HG002_TRUTH = "HG002-v5.0q-GRCh38"
HG002_TRUTH_VCF = "HG002_GRCh38_v5.0q_smvar.vcf.gz"
HG002_BENCH_BED = "HG002_GRCh38_v5.0q_smvar.benchmark.bed"


def two_x_key(gene: str) -> str:
    """The withheld key of an X-linked locus's two-X branch."""
    return f"{gene}:two_x"


def hg002_truth(lab: Path) -> tuple[Path, Path]:
    """GIAB's HG002 v5.0q truth VCF and its benchmark BED in the lab."""
    folder = lab / "giab" / HG002_TRUTH
    return folder / HG002_TRUTH_VCF, folder / HG002_BENCH_BED


def stage3_dir(lab: Path) -> Path:
    """The stage-3 calibration table `genokeeper calibrate` writes (`tables.folder`)."""
    return folder(lab, "stage3")


@dataclass(frozen=True)
class Stage3Table:
    """The versions measured under, the genes withheld (key -> reason: a repeat gene,
    `<gene>:two_x`, `SMN_KEY`, `YLEAF_KEY`, RFC1, CFTR) and capped (gene -> level), each
    library's C9orf72 flank ratio (`pcrfree`, `pcr`), what was measured on each library
    (and, under `RULES`, on each trio sample), for the report."""

    versions: dict[str, str]
    withheld: dict[str, str]
    caps: dict[str, Level]
    flank: dict[str, float]
    measured: dict[str, Any]


def stage3_versions(lab: Path, fasta_sha: str, grading: Grading) -> dict[str, str]:
    """The versions a stage-3 table is measured under: the five images, the catalog's sha256,
    the Parascopy data's md5s (`FileNotFoundError` naming the bootstrap when absent), the
    FASTA's sha256, the grading's `[repeats]` with every `[repeat_loci]` row but its
    `UNMEASURED` fields and its `[exclusion]`, and the trio checks' constants
    (`repeat_rules.MENDEL_TOLERANCE`, `HG002_RFC1`), each as JSON."""
    # Imported here: `repeat_rules` imports this module.
    from genokeeper.calibrate import repeat_rules

    repeats = asdict(grading.repeats)
    for row in repeats["loci"].values():
        for name in UNMEASURED:
            del row[name]
    rules = {
        "mendel_tolerance": repeat_rules.MENDEL_TOLERANCE,
        "hg002_rfc1": list(repeat_rules.HG002_RFC1),
    }
    return {
        **{name: IMAGES[name] for name in STAGE3_IMAGES},
        "catalog": catalog_sha256(),
        "parascopy_data": parascopy_md5s(lab),
        "fasta": fasta_sha,
        "repeats": json.dumps(repeats, sort_keys=True),
        "exclusion": json.dumps(asdict(grading.exclusion), sort_keys=True),
        "rules": json.dumps(rules, sort_keys=True),
    }


def _in_benchmark(bench: Path, chrom: str, start: int, end: int) -> bool:
    """Whether 0-based `[start, end)` lies wholly in one interval of the BED `bench`."""
    with bench.open() as f:
        for line in f:
            fields = line.split("\t")
            if fields[0] == chrom and int(fields[1]) <= start and end <= int(fields[2]):
                return True
    return False


def truth_repeats(
    vcf: Path,
    bench: Path,
    *,
    region: tuple[str, int, int],
    motif: str,
    copies: int = 2,
) -> tuple[int, ...] | None:
    """HG002's repeat count on each of its `copies` haplotypes (1 for a male's chrX) over
    `region` (1-based, inclusive) in GIAB's v5.0q phased small-variant truth `vcf`: the
    reference count plus each haplotype's net length change in records touching the region,
    in motif units, in haplotype order. None unless the region lies wholly in the benchmark
    BED `bench`, and when a record there is filtered, uncalled, of another ploidy or an
    unphased heterozygote, or anchored in the flank and not whole motif units: no truth
    is guessed. A record in the flank that changes no length is skipped first."""
    chrom, start, end = region
    if not _in_benchmark(bench, chrom, start - 1, end):
        return None
    reference = (end - start + 1) // len(motif)
    deltas = [0] * copies
    # A VCF indel is anchored on the base before it, and left normalisation slides one in
    # the repeat into the flank by up to a motif length less one base where the flank
    # ends as a unit does (GIAB anchors HG002's FXN insertions two bases before FXN's GAA).
    # Such a record counts only when what it inserts or deletes is whole units of the
    # motif; one that also changes the flank (FXN's poly-A) gives no truth.
    first = start - len(motif)
    units = {motif[i:] + motif[:i] for i in range(len(motif))}
    with pysam.VariantFile(str(vcf)) as f:
        for rec in f.fetch(chrom, first - 1, end):
            ref = rec.ref or ""
            if rec.pos + len(ref) - 1 < first or rec.pos > end:
                continue
            # A record in the slop that ends before the region and changes no length (a
            # flanking SNV) changes no count, whatever its filter or genotype.
            if rec.pos + len(ref) - 1 < start and all(
                len(a) == len(ref) for a in rec.alleles or ()
            ):
                continue
            if set(rec.filter.keys()) - {"PASS"}:
                return None
            sample = rec.samples[0]
            gt = sample["GT"]
            if gt is None or len(gt) != copies or None in gt:
                return None
            if len(set(gt)) > 1 and not sample.phased:
                return None
            alleles = rec.alleles or ()
            for i, allele in enumerate(gt):
                if not allele:
                    continue
                alt = str(alleles[allele])
                if rec.pos < start - 1:
                    short, long = sorted((alt, ref), key=len)
                    change = long[len(short) :]
                    unit = change[: len(motif)]
                    if (
                        not long.startswith(short)
                        or unit not in units
                        or (change != unit * (len(change) // len(motif)))
                    ):
                        return None
                deltas[i] += len(alt) - len(ref)
    return tuple(reference + round(d / len(motif)) for d in deltas)


def _size_error(
    locus: RepeatLocus, call: RepeatCall, truth: tuple[int, ...] | None
) -> int | None:
    """The largest difference between a called allele shorter than the read and its truth,
    the two sorted and paired; None without a truth or with another allele count."""
    if truth is None or len(truth) != len(call.alleles):
        return None
    pairs = zip(sorted(a.repeats for a in call.alleles), sorted(truth), strict=True)
    errors = [abs(n - t) for n, t in pairs if n * len(locus.motif) < call.read_length]
    return max(errors, default=0)


def measure(
    results: Mapping[str, Stage3Results],
    truth: Mapping[str, tuple[int, ...] | None],
    loci: Mapping[str, RepeatLocus],
) -> tuple[dict[str, str], dict[str, Level], dict[str, Any]]:
    """What the libraries' `results` withhold (key -> reason) and cap (gene -> level), and
    what was measured on each library: per locus its alleles, their intervals, the reads
    behind each, their classes, the truth and the size error; both SMN callers' copy
    numbers; Yleaf's haplogroup and whether it derives L816; the flank ratio. A locus with
    no call or low depth on a library is withheld, and one whose allele count differs from
    a known truth capped at low."""
    withheld: dict[str, str] = {}
    caps: dict[str, Level] = {}
    measured: dict[str, Any] = {}
    smn_ok = l816_ok = True
    haplogroups: set[str | None] = set()
    for library, found in results.items():
        repeats: dict[str, Any] = {}
        for gene, locus in loci.items():
            known = truth.get(gene)
            call = (
                None if found.repeats is None else found.repeats.get(locus.variant_id)
            )
            if call is None or not call.alleles:
                withheld.setdefault(gene, NO_CALL)
                repeats[gene] = {"alleles": None, "classes": None, "truth": known}
                continue
            counts = [a.repeats for a in call.alleles]
            classes = [classify(locus, n) for n in counts]
            if call.low_depth:
                withheld.setdefault(gene, LOW_DEPTH)
            if not graded_normal(locus, call.alleles):
                withheld[gene] = FALSE_CALL
            error = _size_error(locus, call, known)
            if (error is not None and error > MAX_SIZE_ERROR) or (
                known is not None and len(known) != len(counts)
            ):
                caps[gene] = Level.LOW
            repeats[gene] = {
                "alleles": counts,
                "ci": [list(a.ci) for a in call.alleles],
                "reads": [a.spanning + a.flanking + a.inrepeat for a in call.alleles],
                "classes": [str(c) for c in classes],
                "truth": None if known is None else list(known),
                "size_error": error,
                "read_length": call.read_length,
            }
        caller, parascopy = found.caller, found.parascopy
        smn = {
            "smncaller": None if caller is None else [caller.smn1, caller.smn2],
            "parascopy": None
            if parascopy is None
            else [parascopy.smn1, parascopy.smn2],
        }
        smn_ok = smn_ok and smn["smncaller"] == [2, 2] and smn["parascopy"] == [2, 2]
        derived = found.yleaf_out is not None and L816 in derived_markers(
            found.yleaf_out
        )
        l816_ok = l816_ok and derived
        haplogroups.add(None if found.yleaf is None else found.yleaf.name)
        measured[library] = {
            "repeats": repeats,
            "smn": smn,
            "yleaf": {
                "haplogroup": None if found.yleaf is None else found.yleaf.name,
                "l816": derived,
            },
            "flank": found.flank,
        }
    if not smn_ok:
        withheld[SMN_KEY] = SMN_WITHHELD
    if not l816_ok or len(haplogroups) != 1:
        withheld[YLEAF_KEY] = YLEAF_WITHHELD
    return withheld, caps, measured


def write_table(folder: Path, table: Stage3Table) -> None:
    """Write `table` into `folder`, the staging folder `tables.write_table` swaps in whole."""
    document = {
        "format": FORMAT,
        "versions": table.versions,
        "withheld": table.withheld,
        "caps": {gene: str(level) for gene, level in table.caps.items()},
        "flank": table.flank,
        "measured": table.measured,
    }
    folder.mkdir(parents=True, exist_ok=True)
    (folder / TABLE).write_text(json.dumps(document, indent=1, sort_keys=True))


def read_table(folder: Path) -> Stage3Table:
    """The table in `folder`; `ValueError` for another format."""
    document = json.loads((folder / TABLE).read_text())
    if document.get("format") != FORMAT:
        raise ValueError(f"{folder}: stage 3 table format {document.get('format')}")
    return Stage3Table(
        versions=dict(document["versions"]),
        withheld=dict(document["withheld"]),
        caps={gene: Level(level) for gene, level in document["caps"].items()},
        flank={k: float(v) for k, v in document["flank"].items()},
        measured=dict(document["measured"]),
    )


def open_stage3(
    lab: Path, versions: Mapping[str, str]
) -> tuple[Stage3Table | None, str | None]:
    """The lab's stage-3 table when it was measured under `versions` (the run's own,
    `stage3_versions`); otherwise None and why (`tables.open_table`'s reason, or the
    table unreadable)."""
    found = open_table(lab, "stage3", versions=versions)
    if not isinstance(found, Current):
        return None, found.reason
    try:
        table = read_table(found.folder)
    except (OSError, ValueError, KeyError) as e:
        return None, f"stage 3 table unreadable ({' '.join(str(e).split())}); {RUN}"
    return table, None


def _shown(values: list[Any] | None) -> str:
    return (
        "?"
        if values is None
        else "/".join("?" if v is None else str(v) for v in values)
    )


def _rules_lines(table: Stage3Table, verdict: Callable[[str], str]) -> list[str]:
    """The repeat rules' lines: per locus and trio sample the alleles, kind, fewest spanning
    reads and (HG002 pcrfree) the Mendelian check; the two-X branches withheld; RFC1's
    calls; CFTR's counts."""
    rules: dict[str, Any] = table.measured.get(RULES, {})
    samples = list(rules)
    lines = [
        "",
        (
            "repeat rules (per sample: alleles, kind, fewest spanning reads, Mendelian;"
            " verdict):"
        ),
    ]

    def mendel(value: bool | None) -> str:
        return "?" if value is None else "yes" if value else "no"

    genes = list(rules[samples[0]]["repeats"]) if samples else []
    for gene in genes:
        cells = []
        for sample in samples:
            found = rules[sample]["repeats"][gene]
            cell = (
                f"{sample} {_shown(found['alleles'])} {found['kind'] or '?'},"
                f" {'?' if found['spanning'] is None else found['spanning']} spanning"
            )
            if "mendelian" in found:
                cell += f", Mendelian {mendel(found['mendelian'])}"
            cells.append(cell)
        lines.append(f"  {gene}\t" + "\t".join(cells) + f"\t{verdict(gene)}")
    lines += [
        f"  {key}\t{verdict(key)}" for key in sorted(table.withheld) if ":" in key
    ]
    rfc1 = []
    for sample in samples:
        found = rules[sample]["rfc1"]
        rfc1.append(
            f"{sample} ?"
            if found is None
            else f"{sample} {found['call']} ({found['description']})"
        )
    lines.append("  RFC1\t" + "\t".join(rfc1) + f"\t{verdict(RFC1_LOCUS)}")
    cftr = []
    for sample in samples:
        found = rules[sample]["cftr"]
        cftr.append(
            f"{sample} TG {_shown(found['tg'])} T {_shown(found['t'])}"
            f" {found['kind'] or '?'}"
        )
    lines.append("  CFTR\t" + "\t".join(cftr) + f"\t{verdict(CFTR_GENE)}")
    return lines


def report(table: Stage3Table) -> str:
    """Per locus and tool, each library's result against HG002's truth and the verdict,
    both libraries' flank ratios, then the repeat rules (`_rules_lines`)."""
    libraries = [key for key in table.measured if key != RULES]
    lines = ["versions measured under:"]
    lines += [f"  {key}: {value}" for key, value in sorted(table.versions.items())]

    def verdict(key: str) -> str:
        if key in table.withheld:
            return f"withheld: {table.withheld[key]}"
        if key in table.caps:
            return f"capped at {table.caps[key]}"
        return "kept"

    lines += ["", "repeat loci (per library: alleles, truth, size error; verdict):"]
    genes = list(table.measured[libraries[0]]["repeats"]) if libraries else []
    for gene in genes:
        cells = []
        for library in libraries:
            found = table.measured[library]["repeats"][gene]
            error = found.get("size_error")
            ci = found.get("ci")
            interval = (
                "" if ci is None else " CI " + "/".join(f"{a}-{b}" for a, b in ci)
            )
            cells.append(
                f"{library} {_shown(found['alleles'])}{interval}"
                f" (truth {_shown(found['truth'])}, off {'?' if error is None else error})"
            )
        lines.append(f"  {gene}\t" + "\t".join(cells) + f"\t{verdict(gene)}")
    # `[repeats] min_allele_reads` is derived from these: the fewest reads behind any of
    # HG002's alleles.
    fewest = []
    for library in libraries:
        reads = [
            (n, gene)
            for gene, found in table.measured[library]["repeats"].items()
            for n in found.get("reads") or ()
        ]
        if reads:
            n, gene = min(reads)
            fewest.append(f"{library} {n} ({gene})")
    if fewest:
        lines.append(f"  fewest reads behind an allele: {', '.join(fewest)}")
    smn = [
        f"{library} SMNCopyNumberCaller {_shown(table.measured[library]['smn']['smncaller'])},"
        f" Parascopy {_shown(table.measured[library]['smn']['parascopy'])}"
        for library in libraries
    ]
    lines.append("  SMN1/SMN2\t" + "\t".join(smn) + f"\t{verdict(SMN_KEY)}")
    yleaf = [
        f"{library} {table.measured[library]['yleaf']['haplogroup']}"
        f" (L816 {'derived' if table.measured[library]['yleaf']['l816'] else 'not derived'})"
        for library in libraries
    ]
    lines.append("  Y haplogroup\t" + "\t".join(yleaf) + f"\t{verdict(YLEAF_KEY)}")
    ratios = ", ".join(
        f"{library} {ratio:.3f}" for library, ratio in table.flank.items()
    )
    lines += ["", f"C9orf72 flank ratio: {ratios}"]
    lines += _rules_lines(table, verdict)
    return "\n".join(lines)
