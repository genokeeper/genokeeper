"""Repeat expansions as a source (`RepeatSource`; slice 3a, "Grading", repeat expansions;
slice 3e's rules).

Per locus one gene-level statement (`expansionhunter:<gene>`) and one `gene_call` row, each
graded by its locus's `Rule` (`grade`). Under the standard rule only an allele whose
interval's lower bound reaches the pathogenic minimum, with enough reads, is a finding,
worded "screening result, needs confirmation". The other rules are exclusion-only: "normal"
(flagged `exclusion`, with its confirming test) only when `excluded` finds every allele
read across in full, else unassessed, or at most a "possible" finding at validity limited;
the X-linked loci follow the run's X dosage and are never a finding. No statement is
actionable. A tool that did not run, a withheld gene, a missing call and low depth are
unassessed, never normal, at every locus before its rule.
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass, replace

from genokeeper.evidence.pgx import GENE_UNRESOLVED
from genokeeper.repeats.expansionhunter import EH_VERSION
from genokeeper.repeats.loci import (
    Allele,
    RepeatCall,
    RepeatClass,
    RepeatLocus,
    Rfc1Call,
    Rule,
    classify,
)
from genokeeper.repeats.not_covered import NOT_COVERED, NotCovered
from genokeeper.repeats.rfc1 import RFC1_LOCUS, RFC1_MOTIF, RFC1_PAIRS
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Level, Topic, Validity

REPEATS_TOOL = "expansionhunter"
CONFIRM = "screening result, needs confirmation"
FINDING_REASON = f"{CONFIRM} by repeat-primed PCR"
FINDING_REASONS = {"TCF4": f"{CONFIRM} by repeat-primed PCR and an eye examination"}
NORMAL_ONE_X = (
    "normal, {n} {motif} on the one X, read across in full; says nothing about anyone"
    " else's {gene}"
)
NORMAL_TWO_X = (
    "normal: {a} and {b} {motif}, each X's allele read across in full; says nothing about"
    " anyone else's {gene}"
)
NORMAL_BOTH = "normal: {a} and {b} {motif}, each allele read across in full"
X_DOSAGE = (
    "not assessed: the run's X dosage ({copies}) does not match the call; needs a"
    " targeted PCR test"
)
X_NOT_CALIBRATED = (
    "not assessed: the X dosage is not calibrated; run `genokeeper calibrate`"
)
POSSIBLE_ABOVE = (
    f"possible expansion; {CONFIRM}: a size estimate in the disease range, needs"
    " repeat-primed PCR"
)
BETWEEN = (
    "not assessed: an allele estimated above the normal range but below the disease"
    " range; needs repeat-primed PCR"
)
FEW_READS = (
    "not assessed: an allele estimated in the disease range on too few reads; needs"
    " repeat-primed PCR"
)
POSSIBLE_RECESSIVE = (
    f"possible recessive {{condition}}: both alleles carry {{n}} {{motif}}; {CONFIRM}"
    " by PCR sizing"
)
NOT_EXCLUDED = "not assessed ({why}); needs {test}"
NEEDS_PCR_X = "not assessed ({why}); needs a targeted PCR test"
NEEDS_LONG = "not assessed ({why}); needs long-range PCR or long reads"
NOT_SET = "not assessed: the exclusion rule's read threshold was not given"
SAME_LENGTH = "the two alleles are the same length: a longer one may be hidden"
# The test that would confirm a normal each exclusion-only rule flags.
CONFIRM_TEST = {
    Rule.X_LINKED: "a targeted PCR test",
    Rule.EXCLUSION: "PCR sizing",
    Rule.DISTINCT: "PCR sizing",
}
NOT_ASSESSED = "not assessed"
RFC1_NORMAL = (
    "a biallelic AAGGG or ACAGG expansion is not seen; rarer pathogenic motifs (AAAGG,"
    " CAGGG) are not tested"
)
RFC1_UNCLEAR = (
    "not assessed: an expanded allele's motif was not observed or not resolved; needs"
    " repeat-primed PCR or a Southern blot"
)
RFC1_NOT_RUN = "RFC1 motif analysis not run"
RFC1_CONFIRM = "repeat-primed PCR or a Southern blot"
RFC1_CONDITION = (
    "CANVAS (cerebellar ataxia, neuropathy and vestibular areflexia syndrome)"
)
BENIGN_MOTIF = "an expanded allele of a benign motif"
# The analysis's normals that saw every allele's motif, and their levels.
_RFC1_SEEN = {
    "no expanded allele": Level.HIGH,
    "1 expanded benign allele, 1 short reference allele": Level.MEDIUM,
}
_RFC1_UNCLEAR = {
    "expanded allele may exist but not observed",
    "2 expanded alleles may exist but not observed",
    "2 expanded alleles (possibly reference)",
}
CFTR_NO_5T = "no 5T allele: T {t}, TG {tg}"
CFTR_5T = (
    "possible variable-penetrance CFTR allele: 5T with TG{tg}{effect}; "
    f"{CONFIRM} by CFTR poly-T/TG testing"
)
CFTR_5T_ONE = (
    "; alone it causes no disease, in trans with a cystic-fibrosis variant it can cause a"
    " CFTR-related disorder or, in a man, absent vas deferens; the TG and T tracts are"
    " typed apart, so which TG goes with the 5T is not determined"
)
CFTR_5T_BOTH = (
    " on both copies; 5T on both copies beside a long TG can itself cause a CFTR-related"
    " disorder (absent vas deferens in a man)"
)
CFTR_TG11 = "information: 5T with TG11, a CFTR allele of low penetrance{both}"
CFTR_UNREAD = "not assessed: the CFTR TG/T tract was not read across in full"
CFTR_GENE = "CFTR"
CFTR_TG, CFTR_T = "CFTR_TG", "CFTR_T"
CFTR_CONDITION = (
    "CFTR-related disorders (variable-penetrance intron 9 poly-T/TG alleles)"
)
CFTR_SOURCE_LEVEL = "CFTR poly-T/TG (c.1210-34TG[n]T[n])"
NOT_COVERED_LEVEL = "not covered"
TWO_ALLELES_ON_X = "two alleles on a male X"
PCR_LIKE_WARNING = "library looks PCR-amplified; a miss is possible"
NO_SPANNING = "an allele without a spanning read"
BEYOND_FRAGMENT = (
    "allele longer than the fragment: the expansion is detected reliably, its exact size"
    " is not; a lab test (PCR or Southern blot) gives the size"
)

# What a gene's repeat expansion causes, so `--text <disease>` finds its statement.
CONDITIONS = {
    "HTT": "Huntington disease",
    "AR": "spinal and bulbar muscular atrophy (Kennedy disease)",
    "DMPK": "myotonic dystrophy type 1",
    "C9orf72": "frontotemporal dementia and amyotrophic lateral sclerosis",
    "ATXN1": "spinocerebellar ataxia type 1",
    "ATXN2": "spinocerebellar ataxia type 2",
    "ATXN3": "spinocerebellar ataxia type 3 (Machado-Joseph disease)",
    "ATXN7": "spinocerebellar ataxia type 7",
    "CACNA1A": "spinocerebellar ataxia type 6",
    "ATN1": "dentatorubral-pallidoluysian atrophy",
    "FXN": "Friedreich ataxia",
    "FMR1": "fragile X syndrome",
    "TCF4": "Fuchs endothelial corneal dystrophy",
    "PABPN1": "oculopharyngeal muscular dystrophy",
    "AFF2": "fragile X E intellectual disability (FRAXE)",
    "CNBP": "myotonic dystrophy type 2",
    "ATXN10": "spinocerebellar ataxia type 10",
    "NOP56": "spinocerebellar ataxia type 36",
    "PPP2R2B": "spinocerebellar ataxia type 12",
    "ZFHX3": "spinocerebellar ataxia type 4",
    "THAP11": "spinocerebellar ataxia type 51",
    "TBP": "spinocerebellar ataxia type 17",
    "FGF14": "spinocerebellar ataxia type 27B",
    "JPH3": "Huntington disease-like 2",
}

_VALIDITY = {"Definitive": Validity.STRONG, "Strong": Validity.MODERATE}
_RANK = {Level.LOW: 0, Level.MEDIUM: 1, Level.HIGH: 2}
_SEVERITY = {
    RepeatClass.NORMAL: 0,
    RepeatClass.INTERMEDIATE: 1,
    RepeatClass.GAP: 2,
    RepeatClass.REDUCED_PENETRANCE: 3,
    RepeatClass.PATHOGENIC: 4,
}


def _support(allele: Allele) -> int:
    return allele.spanning + allele.flanking + allele.inrepeat


def _length(locus: RepeatLocus, allele: Allele) -> int:
    return allele.repeats * len(locus.motif)


def _worst(locus: RepeatLocus, alleles: tuple[Allele, ...]) -> Allele:
    return max(alleles, key=lambda a: _SEVERITY[classify(locus, a.ci[1])])


def graded_normal(locus: RepeatLocus, alleles: tuple[Allele, ...]) -> bool:
    """Whether the run grades `alleles` normal (a non-finding): every interval's upper bound
    within the normal range. The stage-3 calibration withholds a locus HG002 fails it at."""
    return all(a.ci[1] <= locus.normal_max for a in alleles)


def _watch(locus: RepeatLocus, allele: Allele, min_reads: int) -> tuple[str, str]:
    """The watch's phenotype and the words naming its class for `allele`."""
    lo, hi = classify(locus, allele.ci[0]), classify(locus, allele.ci[1])
    if lo is not hi or hi is RepeatClass.GAP:
        return "borderline", "borderline, indeterminate"
    if hi is RepeatClass.PATHOGENIC:
        return "expanded, too few reads", (
            f"expanded range, but {_support(allele)} supporting reads (needs {min_reads})"
        )
    return str(hi), str(hi)


def _level(
    locus: RepeatLocus,
    call: RepeatCall,
    source: RepeatSource,
    start: Level | None = None,
) -> tuple[Level, list[str]]:
    """The gene call's reliability and its warnings; `start`, a grade's own level, stands
    in for the read-based one before the calibration's cap."""
    warnings: list[str] = []
    level = Level.HIGH
    if any(_length(locus, a) > source.insert for a in call.alleles):
        level = Level.MEDIUM
        warnings.append(BEYOND_FRAGMENT)
    if source.male and locus.chrom == "chrX" and len(call.alleles) == 2:
        level = Level.LOW
        warnings.append(TWO_ALLELES_ON_X)
    if level is Level.HIGH:
        if any(_length(locus, a) > call.read_length for a in call.alleles):
            level = Level.MEDIUM
            warnings.append("allele longer than the read: its size is an estimate")
        if not all(a.spanning > 0 for a in call.alleles):
            level = Level.MEDIUM
            warnings.append(NO_SPANNING)
    if start is not None:
        level = start
    cap = source.caps.get(locus.gene)
    if cap is not None and _RANK[cap] < _RANK[level]:
        level = cap
        warnings.append(f"capped at {cap} by the calibration")
    return level, warnings


@dataclass(frozen=True)
class Graded:
    """One call's grade: the statement's kind and reason, the gene call's phenotype, the
    level before the calibration's cap (None: `_level` decides from the reads), whether it
    is a normal flagged `exclusion` and the test that would confirm it, and the validity
    (None: the locus's)."""

    kind: Kind
    reason: str | None
    phenotype: str | None
    level: Level | None = None
    exclusion: bool = False
    confirm: str | None = None
    validity: Validity | None = None


def _unassessed(reason: str) -> Graded:
    return Graded(Kind.UNASSESSED, reason, NOT_ASSESSED)


def unread(call: RepeatCall, *, min_spanning: int, flanking: bool = True) -> str | None:
    """Why `call`'s alleles were not all read across in full, if they were not: an allele
    with fewer than `min_spanning` spanning reads, a read inside the repeat (the variant's
    count, `RepeatCall.inrepeat`), or (when `flanking`) a flanking read longer than every
    called allele."""
    alleles = call.alleles
    if any(a.spanning < min_spanning for a in alleles):
        return f"an allele has fewer than {min_spanning} spanning reads"
    if call.inrepeat > 0:
        return "reads inside the repeat: a longer allele may be hidden"
    if flanking and alleles and call.longest_flanking > max(a.repeats for a in alleles):
        return (
            "a flanking read is longer than every called allele: a longer allele may be"
            " hidden"
        )
    return None


def excluded(
    locus: RepeatLocus, call: RepeatCall, *, min_spanning: int, alleles: int
) -> str | None:
    """None when the exclusion rule holds for `call` (every allele read across in full,
    every interval's upper bound within the normal range, `alleles` alleles and, for two,
    of different lengths: the homozygous case is not measured); else the first unmet
    condition. The flanking-read condition guards a homozygous-looking or one-X call, not
    two distinct autosomal alleles, as the closing design's exclusion rule writes it: two
    distinct spanned alleles leave no room for a hidden one, while a flanking read past a
    one-X call's allele can be size mosaicism."""
    found = call.alleles
    homozygous = len(found) == 2 and found[0].repeats == found[1].repeats
    why = unread(
        call, min_spanning=min_spanning, flanking=homozygous or len(found) == 1
    )
    if why is not None:
        return why
    if not graded_normal(locus, found):
        return "an allele's upper bound is above the normal range"
    if len(found) != alleles:
        return f"{len(found)} alleles reported, {alleles} expected"
    if homozygous:
        return SAME_LENGTH
    return None


@dataclass(frozen=True)
class XDosage:
    """The run's X dosage under a calibrated [dosage] (slice 3c): its X copies (None:
    between copies), its X depth ratio, and whether a Y is present (a Y ratio above
    `[dosage] y_absent_max`)."""

    copies: int | None
    ratio: float
    y: bool


def _x_dosage(locus: RepeatLocus, *, male: bool, x: XDosage | None) -> Graded | None:
    """The X-dosage rule: an X-linked locus is graded only when the run's X dosage matches
    the call's ploidy (one X called male; two X and no Y called female); unassessed naming
    `genokeeper calibrate` when it is not calibrated, else naming the measured dosage."""
    if locus.rule is not Rule.X_LINKED:
        return None
    if x is None:
        return _unassessed(X_NOT_CALIBRATED)
    # An XXY sample called male shows one X's allele and hides the other, and one called
    # female is unassessed as the closing design writes it. Its anchor comes from the run's
    # measured X and Y depth ratios (slice 3c's dosage), which the sex ExpansionHunter was
    # run with cannot move because the ratios are read from depth, not from the flag.
    if (male and x.copies == 1) or (not male and x.copies == 2 and not x.y):
        return None
    copies = (
        f"X ratio {x.ratio:.2f}, between copies"
        if x.copies is None
        else f"{x.copies} X" + (" and a Y" if x.y else "")
    )
    return _unassessed(X_DOSAGE.format(copies=copies))


def _standard(locus: RepeatLocus, call: RepeatCall, min_allele_reads: int) -> Graded:
    """Slice 3a's grade: a finding from an allele reaching the disease minimum with enough
    reads, a possible recessive finding from both alleles at `recessive_at` or above, a
    normal, or a watch."""
    alleles = call.alleles
    if any(
        a.ci[0] >= locus.pathogenic_min and _support(a) >= min_allele_reads
        for a in alleles
    ):
        reason = FINDING_REASONS.get(locus.gene, FINDING_REASON)
        return Graded(Kind.FINDING, reason, "expanded (screening)")
    low = locus.recessive_at
    if (
        low is not None
        and len(alleles) == 2
        and all(low <= a.ci[0] and a.ci[1] < locus.pathogenic_min for a in alleles)
    ):
        reason = POSSIBLE_RECESSIVE.format(
            condition=CONDITIONS.get(locus.gene, locus.gene),
            n=" and ".join(dict.fromkeys(str(a.repeats) for a in alleles)),
            motif=locus.motif,
        )
        return Graded(
            Kind.FINDING,
            reason,
            "possible recessive (screening)",
            level=Level.LOW,
            validity=Validity.LIMITED,
        )
    if graded_normal(locus, alleles):
        return Graded(Kind.NON_FINDING, None, "normal")
    phenotype, words = _watch(locus, _worst(locus, alleles), min_allele_reads)
    reason = f"{CONFIRM}: {words}"
    return Graded(
        Kind.WATCH, f"{reason}; {locus.note}" if locus.note else reason, phenotype
    )


def grade(
    locus: RepeatLocus,
    call: RepeatCall,
    *,
    male: bool,
    x_dosage: XDosage | None,
    min_allele_reads: int,
    min_spanning: int | None,
) -> Graded:
    """One call's grade by its locus's rule, after the run's withheld and not-run checks;
    an X-linked locus takes the X-dosage rule first."""
    dosage = _x_dosage(locus, male=male, x=x_dosage)
    if dosage is not None:
        return dosage
    if locus.rule is Rule.STANDARD:
        return _standard(locus, call, min_allele_reads)
    if min_spanning is None:
        return _unassessed(NOT_SET)
    normal = Graded(
        Kind.NON_FINDING,
        None,
        "normal",
        exclusion=True,
        confirm=CONFIRM_TEST[locus.rule],
    )
    found = call.alleles
    if locus.rule is Rule.X_LINKED:
        why = excluded(locus, call, min_spanning=min_spanning, alleles=1 if male else 2)
        if why is not None:
            return _unassessed(NEEDS_PCR_X.format(why=why))
        words = (
            NORMAL_ONE_X.format(n=found[0].repeats, motif=locus.motif, gene=locus.gene)
            if male
            else NORMAL_TWO_X.format(
                a=found[0].repeats,
                b=found[1].repeats,
                motif=locus.motif,
                gene=locus.gene,
            )
        )
        return replace(normal, reason=words)
    why = excluded(locus, call, min_spanning=min_spanning, alleles=2)
    if why is None:
        words = NORMAL_BOTH.format(
            a=found[0].repeats, b=found[1].repeats, motif=locus.motif
        )
        return replace(normal, reason=words)
    if locus.rule is Rule.DISTINCT:
        return _unassessed(NEEDS_LONG.format(why=why))
    if any(
        a.ci[0] >= locus.pathogenic_min and _support(a) >= min_allele_reads
        for a in found
    ):
        return Graded(
            Kind.FINDING,
            POSSIBLE_ABOVE,
            "possible expansion (screening)",
            level=Level.LOW,
            validity=Validity.LIMITED,
        )
    if any(a.ci[0] >= locus.pathogenic_min for a in found):
        return _unassessed(FEW_READS)
    if any(a.ci[0] > locus.normal_max for a in found):
        return _unassessed(BETWEEN)
    # Two equal alleles need sizing; any other doubt, the test that sees a long allele.
    test = CONFIRM_TEST[locus.rule] if why == SAME_LENGTH else "repeat-primed PCR"
    return _unassessed(NOT_EXCLUDED.format(why=why, test=test))


def rfc1_grade(call: Rfc1Call) -> Graded:
    """RFC1's grade from its motif analysis's (label, description) pair: a pair outside
    `RFC1_PAIRS` is unassessed whatever its label; a normal flagged `exclusion` only when
    every allele's motif was seen, unassessed when an expanded allele's was not, and a
    "possible" finding at validity limited for a carrier or affected call; never a watch."""
    label, description = call.call, call.description
    if (label, description) not in RFC1_PAIRS:
        return _unassessed(f"not assessed: unknown RFC1 call {label!r} ({description})")
    if description in _RFC1_SEEN:
        return Graded(
            Kind.NON_FINDING,
            RFC1_NORMAL,
            label,
            level=_RFC1_SEEN[description],
            exclusion=True,
            confirm=RFC1_CONFIRM,
        )
    if description in _RFC1_UNCLEAR:
        return _unassessed(RFC1_UNCLEAR)
    return Graded(
        Kind.FINDING,
        f"possible {label}; {CONFIRM}: {description}; needs {RFC1_CONFIRM}",
        f"possible {label} (screening)",
        level=Level.LOW,
        validity=Validity.LIMITED,
    )


def cftr_grade(tg: RepeatCall, t: RepeatCall, *, min_spanning: int) -> Graded:
    """CFTR's intron 9 poly-T/TG: unassessed unless both tracts were read across in full;
    no 5T allele a non-finding; a 5T a finding at validity limited, "possible" with a TG
    allele of 12 or more and "information" (low penetrance) with every TG 11 or fewer.
    Never cystic-fibrosis carrier status."""
    if any(unread(c, min_spanning=min_spanning) for c in (tg, t)):
        return _unassessed(CFTR_UNREAD)
    ts = [a.repeats for a in t.alleles]
    tgs = [a.repeats for a in tg.alleles]
    if 5 not in ts:
        reason = CFTR_NO_5T.format(t="/".join(map(str, ts)), tg="/".join(map(str, tgs)))
        return Graded(Kind.NON_FINDING, reason, "no 5T")
    both = ts.count(5) == 2
    long = sorted({n for n in tgs if n >= 12})
    reason = (
        CFTR_5T.format(
            tg="/".join(map(str, long)),
            effect=CFTR_5T_BOTH if both else CFTR_5T_ONE,
        )
        if long
        else CFTR_TG11.format(both=" on both copies" if both else "")
    )
    return Graded(
        Kind.FINDING,
        reason,
        "5T (screening)",
        level=Level.LOW,
        validity=Validity.LIMITED,
    )


def _allele_fields(alleles: tuple[Allele, ...]) -> list[dict[str, object]]:
    return [
        {
            "repeats": a.repeats,
            "ci": list(a.ci),
            "spanning": a.spanning,
            "flanking": a.flanking,
            "inrepeat": a.inrepeat,
        }
        for a in alleles
    ]


def _callers(alleles: tuple[Allele, ...]) -> dict[str, str]:
    if not alleles:
        return {}
    counts = "/".join(str(a.repeats) for a in alleles)
    cis = "/".join(f"{a.ci[0]}-{a.ci[1]}" for a in alleles)
    reads = "; ".join(f"{a.spanning}/{a.flanking}/{a.inrepeat}" for a in alleles)
    return {f"ExpansionHunter {EH_VERSION}": (f"{counts}, CI {cis}, reads {reads}")}


@dataclass(frozen=True)
class RepeatSource:
    """ExpansionHunter's calls by catalog VariantId (None: not run, for `not_run`), the
    loci's thresholds, the genes the calibration withholds (gene -> reason) and caps
    (gene -> level), REViewer's images (gene -> SVG path), whether the library looks
    PCR-amplified, the fragment length (`insert`, bp), the sample's sex and the reads an
    allele needs; the exclusion rule's spanning reads (None: every non-standard locus
    unassessed), the run's X dosage (None: [dosage] not calibrated, or not measured)
    and RFC1's motif call (None:
    `rfc1_not_run` says why, or there was none)."""

    calls: Mapping[str, RepeatCall] | None
    not_run: str | None
    loci: Mapping[str, RepeatLocus]
    withheld: Mapping[str, str]
    caps: Mapping[str, Level]
    images: Mapping[str, str]
    pcr_like: bool
    insert: float
    male: bool
    min_allele_reads: int
    tool: str = REPEATS_TOOL
    min_spanning: int | None = None
    x_dosage: XDosage | None = None
    rfc1: Rfc1Call | None = None
    rfc1_not_run: str | None = RFC1_NOT_RUN

    def _unassessed(self, locus: RepeatLocus) -> str | None:
        """Why the locus cannot be graded, if it cannot."""
        call = None if self.calls is None else self.calls.get(locus.variant_id)
        if self.not_run is not None:
            return self.not_run
        if locus.gene in self.withheld:
            return self.withheld[locus.gene]
        if call is None or not call.alleles:
            return "no call at this locus"
        if call.low_depth:
            return "low depth at this locus"
        return None

    def collect(self, ctx: RunContext) -> Collected:
        """One statement and one gene call per locus, RFC1 and CFTR, and one statement per
        not-covered locus."""
        rows: list[Row] = []
        gene_calls: list[Fields] = []
        graded = [self._locus(ctx, locus) for locus in self.loci.values()]
        for statement, gene_call in [*graded, self._rfc1(ctx), self._cftr(ctx)]:
            rows.append(Row(None, None, statement))
            gene_calls.append(gene_call)
        rows += [Row(None, None, self._not_covered(ctx, n)) for n in NOT_COVERED]
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={"gene_call": gene_calls},
            input_rows=len(rows),
            dropped={},
        )

    def _graded(self, locus: RepeatLocus, call: RepeatCall | None) -> Graded | None:
        """The locus's grade, or None when 3a's checks (not run, withheld, no call, low
        depth) leave it unassessed: at every locus, the X-linked ones included, before its
        rule."""
        if self._unassessed(locus) is not None or call is None:
            return None
        return grade(
            locus,
            call,
            male=self.male,
            x_dosage=self.x_dosage,
            min_allele_reads=self.min_allele_reads,
            min_spanning=self.min_spanning,
        )

    def _locus(self, ctx: RunContext, locus: RepeatLocus) -> tuple[Fields, Fields]:
        call = None if self.calls is None else self.calls.get(locus.variant_id)
        alleles = () if call is None else call.alleles
        graded = self._graded(locus, call)
        if graded is None:
            why = self._unassessed(locus)
            assert why is not None
            graded = Graded(Kind.UNASSESSED, why, None)
        level = Level.LOW
        warnings = [f"{GENE_UNRESOLVED}: {graded.reason}"]
        assessed = graded.kind is not Kind.UNASSESSED
        if assessed:
            assert call is not None
            level, warnings = _level(locus, call, self, graded.level)
            if self.pcr_like and locus.gene == "C9orf72":
                warnings.append(PCR_LIKE_WARNING)
        validity = (
            graded.validity
            or locus.validity
            or _VALIDITY.get(locus.evidence, Validity.LIMITED)
        )
        claim = {
            "locus_id": locus.locus_id,
            "motif": locus.motif,
            "alleles": _allele_fields(alleles),
            "class": str(classify(locus, _worst(locus, alleles).repeats))
            if assessed
            else None,
            "thresholds": {
                "normal_max": locus.normal_max,
                "intermediate": locus.intermediate,
                "reduced_penetrance": locus.reduced_penetrance,
                "pathogenic_min": locus.pathogenic_min,
            },
            "image": self.images.get(locus.gene),
        }
        statement = self._statement(
            ctx,
            gene=locus.gene,
            condition=CONDITIONS.get(locus.gene),
            graded=graded,
            claim=claim,
            validity=validity,
            topic=locus.topic,
            source_level=locus.evidence,
            source_id=locus.locus_id,
        )
        diplotype = f"{locus.motif} " + "/".join(str(a.repeats) for a in alleles)
        gene_call = self._gene_call(
            ctx,
            locus.gene,
            graded,
            diplotype=diplotype if assessed else None,
            callers=_callers(alleles),
            level=level,
            warnings=warnings,
        )
        return statement, gene_call

    def _statement(
        self,
        ctx: RunContext,
        *,
        gene: str,
        condition: str | None,
        graded: Graded,
        claim: Mapping[str, object],
        validity: Validity,
        topic: Topic,
        source_level: str,
        source_id: str,
        rfc1: Rfc1Call | None = None,
    ) -> Fields:
        """The statement of one graded repeat claim: `claim` holds its locus's keys."""
        full = {
            "schema": CLAIM_SCHEMA,
            "type": "repeat_expansion",
            "gene": gene,
            "condition": condition,
            **claim,
            "confirm": graded.confirm if graded.exclusion else None,
            "rfc1": None
            if rfc1 is None
            else {"call": rfc1.call, "description": rfc1.description},
        }
        return {
            "run_id": ctx.run_id,
            "statement_id": f"{REPEATS_TOOL}:{gene}",
            "kind": str(graded.kind),
            "vid": None,
            "gene": gene,
            "claim": json.dumps(full),
            "validity": str(validity),
            "source_level": source_level,
            "reason": graded.reason,
            "effect": None,
            "actionable": False,
            "topic": str(topic),
            "source": REPEATS_TOOL,
            "source_id": source_id,
            "source_version": EH_VERSION,
            "points": None,
            "exclusion": graded.exclusion,
        }

    def _gene_call(
        self,
        ctx: RunContext,
        gene: str,
        graded: Graded,
        *,
        diplotype: str | None,
        callers: dict[str, str],
        level: Level,
        warnings: list[str],
    ) -> Fields:
        return {
            "run_id": ctx.run_id,
            "gene": gene,
            "diplotype": diplotype,
            "phenotype": graded.phenotype,
            "activity_score": None,
            "callers": json.dumps(callers),
            "agreed": None,
            "level": str(level),
            "warnings": warnings,
        }

    def _capped(self, gene: str, graded: Graded, warnings: list[str]) -> Level:
        """An RFC1 or CFTR grade's level (low when unassessed, high when not set), capped
        by the calibration."""
        if graded.kind is Kind.UNASSESSED:
            return Level.LOW
        level = graded.level or Level.HIGH
        cap = self.caps.get(gene)
        if cap is not None and _RANK[cap] < _RANK[level]:
            warnings.append(f"capped at {cap} by the calibration")
            return cap
        return level

    def _rfc1_unassessed(self) -> str | None:
        """Why RFC1 cannot be graded, if it cannot."""
        if self.not_run is not None:
            return self.not_run
        if self.rfc1_not_run is not None:
            return self.rfc1_not_run
        if RFC1_LOCUS in self.withheld:
            return self.withheld[RFC1_LOCUS]
        if self.rfc1 is None:
            return "no call at this locus"
        if self.rfc1.low_depth:
            return "low depth at this locus"
        return None

    def _rfc1(self, ctx: RunContext) -> tuple[Fields, Fields]:
        """RFC1's statement and gene call, graded by its motif call."""
        call = self.rfc1
        alleles = () if call is None else call.alleles
        why = self._rfc1_unassessed()
        graded = (
            rfc1_grade(call)
            if why is None and call is not None
            else Graded(Kind.UNASSESSED, why, None)
        )
        warnings = (
            [f"{GENE_UNRESOLVED}: {graded.reason}"]
            if graded.kind is Kind.UNASSESSED
            else []
        )
        if graded.kind is Kind.NON_FINDING and graded.level is Level.MEDIUM:
            warnings.append(BENIGN_MOTIF)
        level = self._capped(RFC1_LOCUS, graded, warnings)
        claim = {
            "locus_id": RFC1_LOCUS,
            "motif": RFC1_MOTIF,
            "alleles": _allele_fields(alleles),
            "class": None,
            "thresholds": None,
            "image": None,
        }
        statement = self._statement(
            ctx,
            gene=RFC1_LOCUS,
            condition=RFC1_CONDITION,
            graded=graded,
            claim=claim,
            # STRchive: Definitive.
            validity=graded.validity or Validity.STRONG,
            topic=Topic.RISKS,
            source_level="Definitive",
            source_id=RFC1_LOCUS,
            rfc1=call if why is None else None,
        )
        diplotype = f"{RFC1_MOTIF} " + "/".join(str(a.repeats) for a in alleles)
        gene_call = self._gene_call(
            ctx,
            RFC1_LOCUS,
            graded,
            diplotype=diplotype if graded.kind is not Kind.UNASSESSED else None,
            callers=_callers(alleles),
            level=level,
            warnings=warnings,
        )
        return statement, gene_call

    def _cftr(self, ctx: RunContext) -> tuple[Fields, Fields]:
        """CFTR's statement and gene call, graded on its TG and T tracts."""
        calls = {} if self.calls is None else self.calls
        tg, t = calls.get(CFTR_TG), calls.get(CFTR_T)
        if self.not_run is not None:
            graded = Graded(Kind.UNASSESSED, self.not_run, None)
        elif CFTR_GENE in self.withheld:
            graded = Graded(Kind.UNASSESSED, self.withheld[CFTR_GENE], None)
        elif tg is None or t is None or not tg.alleles or not t.alleles:
            graded = Graded(Kind.UNASSESSED, "no call at this locus", None)
        elif tg.low_depth or t.low_depth:
            graded = Graded(Kind.UNASSESSED, "low depth at this locus", None)
        elif self.min_spanning is None:
            graded = _unassessed(NOT_SET)
        else:
            graded = cftr_grade(tg, t, min_spanning=self.min_spanning)
        warnings = (
            [f"{GENE_UNRESOLVED}: {graded.reason}"]
            if graded.kind is Kind.UNASSESSED
            else []
        )
        level = self._capped(CFTR_GENE, graded, warnings)
        tgs = () if tg is None else tg.alleles
        ts = () if t is None else t.alleles
        claim = {
            "locus_id": CFTR_GENE,
            "motif": "TG/T",
            "alleles": {"tg": _allele_fields(tgs), "t": _allele_fields(ts)},
            "class": None,
            "thresholds": None,
            "image": None,
        }
        statement = self._statement(
            ctx,
            gene=CFTR_GENE,
            condition=CFTR_CONDITION,
            graded=graded,
            claim=claim,
            validity=Validity.LIMITED,
            topic=Topic.CARRIER,
            source_level=CFTR_SOURCE_LEVEL,
            source_id=CFTR_GENE,
        )
        diplotype = (
            "TG " + "/".join(str(a.repeats) for a in tgs)
            + " T " + "/".join(str(a.repeats) for a in ts)
        )  # fmt: skip
        tracts = [
            f"{name} {found}"
            for name, alleles in (("TG", tgs), ("T", ts))
            for found in _callers(alleles).values()
        ]
        callers = {f"ExpansionHunter {EH_VERSION}": "; ".join(tracts)} if tracts else {}
        gene_call = self._gene_call(
            ctx,
            CFTR_GENE,
            graded,
            diplotype=diplotype if graded.kind is not Kind.UNASSESSED else None,
            callers=callers,
            level=level,
            warnings=warnings,
        )
        return statement, gene_call

    def _not_covered(self, ctx: RunContext, locus: NotCovered) -> Fields:
        """A not-covered locus's statement: unassessed with its reason, no gene call."""
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "repeat_expansion",
            "gene": locus.gene,
            "condition": locus.condition,
        }
        return {
            "run_id": ctx.run_id,
            "statement_id": f"{REPEATS_TOOL}:{locus.gene}",
            "kind": str(Kind.UNASSESSED),
            "vid": None,
            "gene": locus.gene,
            "claim": json.dumps(claim),
            "validity": str(_VALIDITY.get(locus.evidence, Validity.LIMITED)),
            "source_level": NOT_COVERED_LEVEL,
            "reason": locus.reason,
            "effect": None,
            "actionable": False,
            "topic": str(locus.topic),
            "source": REPEATS_TOOL,
            "source_id": locus.gene,
            "source_version": EH_VERSION,
            "points": None,
            "exclusion": False,
        }
