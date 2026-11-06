"""Deletions of 50 bp or more in the listed genes' exons as a source (`DeletionSource`;
slice 3b).

One statement for the screen itself (`deletions:screen`, every run): unassessed with its
reason when no caller ran, else a non-finding at validity limited naming how many deletions
were called and shown and the screen's recall on HG002 (`DeletionTable.recall`). Then one
statement per shown deletion and gene, at validity limited and never actionable. A deletion
is a loss of function, so it goes through 1e's gene gate (`lof_gene_gate`): where that lets
it count, it is placed as a ClinVar pathogenic allele there would be (`placement`); on both
copies of a recessive gene outside the carrier lists (Tier 3 and 4, `carrier_genes`) and
the haploinsufficient genes it is under
research; elsewhere a non-finding, its reason saying why. Each is worded by its class's
chance of being false: "likely present" below `likely_below`, else "possible", always "needs confirmation by MLPA, a microarray, qPCR or
long reads"; a class whose band the table did not measure is withheld (unassessed, no vid,
`NOT_MEASURED`, its claim's `withheld` keeping the grade it would have). Each deletion's `variant` and `call` rows go on its first statement only.

Slice 3d's copy-number and insertion screens are the same source in their own words
(`Wording`: the tool every row names, the noun, the claim type, the test that confirms, the
ALT and the callers), with their own placement (`place`, which reads `Shown.kinds`), a
per-gene override of the gate (`placement_overrides`, A4: a DPYD loss is a drugs finding)
and the depth callers that also support a 3b deletion (`depth_support`, a reason line that
changes no class: 3b's DELLY ratio stays the one depth feature, C6). A wording's `qualifier`
names a call's kind in its "possible" reason (Scramble's mobile-element family: "possible Alu
insertion"), and `screen_claim` joins the screen statement's claim (what the screen does
not cover).
"""

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

from genokeeper.deletions.screen import CallClass, Deletion, Shown, warnings
from genokeeper.evidence.clinvar import placement
from genokeeper.evidence.predicted import BOTH_COPIES, GENETICIST_NOTE, lof_gene_gate
from genokeeper.genelists import carrier_genes, gene_mode
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA, WITHHELD, withheld_entry
from genokeeper.store.db import Kind, Level, RowType, Topic, Validity
from genokeeper.variants.key import make_vid

DELETIONS_TOOL = "deletions"
DELLY_VERSION = "2.7.0"
DYSGU_VERSION = "1.9.0"
TRUVARI_VERSION = "5.5.0"
SOURCE_LEVEL = f"DELLY {DELLY_VERSION} and Dysgu {DYSGU_VERSION}, matched by Truvari {TRUVARI_VERSION}"
SOURCE_VERSION = (
    f"DELLY {DELLY_VERSION}; Dysgu {DYSGU_VERSION}; Truvari {TRUVARI_VERSION}"
)
CONFIRM_BY = "needs confirmation by MLPA, a microarray, qPCR or long reads"
# Each call's warnings after its class's (`screen.warnings`): HG002's truth carries no
# zygosity the table measures, and the screen reads no transcript consequence.
ZYGOSITY_WARNING = "both copies called; zygosity is not calibrated"
CONSEQUENCE_WARNING = (
    "its consequence (whole or partial exon, in frame, last exon, decay escape) is not"
    " assessed"
)

# A gene's reason line where the call touches only the margin beside its exons
# (`Shown.flank`): placed as an exon hit, as the margin holds the splice sites.
FLANK_ONLY = "touches only the bases beside an exon, not the exon itself"

type Rate = tuple[float, float, float]
# Why a class whose band the table did not measure (`p_wrong` None) is withheld.
NOT_MEASURED = "not calibrated: no {} calls measured"


type Placed = tuple[Kind, Topic | None, str | None]

# What a withheld call's `value` holds of its claim.
_EVENT = ("chrom", "start", "end", "length", "exons", "callers", "depth_ratio")


@dataclass(frozen=True)
class Wording:
    """What a screen's rows say: the tool every row names (`source`, the statement ids'
    prefix), the noun its reasons use, the claim type, the test that confirms one,
    the ALT of its vids, its source level and version, the two callers whose evidence
    the class names (`screen.warnings`), what names a call's kind before the noun in its
    "possible" reason (None: nothing), and how far `Deletion.pos` lies before its record's
    POS (`pos_offset`: 0 for a padded symbolic record, 1 for a depth event, whose POS is its
    first base), so the vid and the haploid look-up name the record's own POS."""

    tool: str
    noun: str
    claim_type: str
    confirm_by: str
    alt: str
    source_level: str
    source_version: str
    callers: tuple[str, str]
    qualifier: Callable[[Deletion], str | None] | None = None
    pos_offset: int = 0

    def vid(self, d: Deletion) -> str:
        return make_vid(d.chrom, d.pos + self.pos_offset, d.ref, self.alt)

    @property
    def nouns(self) -> str:
        return f"{self.noun}es" if self.noun.endswith("s") else f"{self.noun}s"

    @property
    def likely(self) -> str:
        return f"likely present, {self.confirm_by}"

    @property
    def possible(self) -> str:
        return f"possible {self.noun}, {self.confirm_by}"

    def possible_for(self, d: Deletion) -> str:
        """`possible`, the call's kind before the noun where `qualifier` names one."""
        kind = None if self.qualifier is None else self.qualifier(d)
        noun = self.noun if kind is None else f"{kind} {self.noun}"
        return f"possible {noun}, {self.confirm_by}"


DELETION_WORDING = Wording(
    tool=DELETIONS_TOOL,
    noun="deletion",
    claim_type="deletion",
    confirm_by=CONFIRM_BY,
    alt="<DEL>",
    source_level=SOURCE_LEVEL,
    source_version=SOURCE_VERSION,
    callers=("dysgu", "delly"),
)
LIKELY_REASON = DELETION_WORDING.likely
POSSIBLE_REASON = DELETION_WORDING.possible


@dataclass(frozen=True)
class DeletionSource:
    """The screen's shown deletions, how many deletions the callers made (`called`), the
    screen's drop counts by reason, why the screen did not run (`not_run`), each class's
    chance of being false (rate, lo, hi; None where the table has none), ClinGen's
    haploinsufficient genes (`lof_gene_gate`), and per band the truth's rare deletions
    the screen found on HG002 and their number (`DeletionTable.recall`). Then the slice-3d
    hooks the module docstring names, each 3b's behaviour by default; `screen_row` False
    writes no `*:screen` row (a source that writes its own); `cautions` the lines a call's
    warnings carry after its class's where its rate is extrapolated
    (`DeletionTable.cautions`; none by default)."""

    shown: Sequence[Shown]
    called: int
    dropped: Mapping[str, int]
    not_run: str | None
    p_wrong: Callable[[CallClass], Rate | None]
    haploinsufficient: frozenset[str]
    recall: Mapping[str, Mapping[str, tuple[int, int]]]
    wording: Wording = DELETION_WORDING
    place: Callable[[RunContext, str, Shown], Placed] | None = None
    placement_overrides: Mapping[str, tuple[str, str]] = field(
        default_factory=dict[str, tuple[str, str]]
    )
    screen_row: bool = True
    depth_support: Mapping[tuple[str, int, int], str] = field(
        default_factory=dict[tuple[str, int, int], str]
    )
    screen_claim: Mapping[str, object] = field(default_factory=dict[str, object])
    cautions: Callable[[Shown], list[str]] = lambda _: []

    @property
    def tool(self) -> str:
        return self.wording.tool

    def _statement(
        self, ctx: RunContext, statement_id: str, **fields: object
    ) -> Fields:
        w = self.wording
        return {
            "run_id": ctx.run_id,
            "statement_id": statement_id,
            "effect": None,
            "actionable": False,
            "source_level": w.source_level,
            "source": w.tool,
            "source_version": w.source_version,
            "points": None,
            **fields,
        }

    def _screen_row(
        self,
        ctx: RunContext,
        *,
        shown: int | None = None,
        extra: Mapping[str, object] | None = None,
    ) -> Row:
        """The `<tool>:screen` statement; `shown` counts what was shown where that is more
        than `self.shown` (a source with rows of its own), and `extra` joins its claim."""
        count = len(self.shown) if shown is None else shown
        if self.not_run is not None:
            kind, reason, validity = Kind.UNASSESSED, self.not_run, Validity.STRONG
        else:
            kind, validity = Kind.NON_FINDING, Validity.LIMITED
            reason = (
                f"screened: {self.called} {self.wording.nouns} called,"
                f" {count} shown;"
                f" a miss is possible: {_recall(self.recall, self.wording.nouns)}"
            )
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": f"{self.wording.claim_type}_screen",
            "called": None if self.not_run is not None else self.called,
            "shown": None if self.not_run is not None else count,
            **(extra or {}),
        }
        statement = self._statement(
            ctx,
            f"{self.tool}:screen",
            kind=str(kind),
            vid=None,
            gene=None,
            claim=json.dumps(claim),
            validity=str(validity),
            reason=reason,
            topic=str(Topic.CARRIER),
            source_id="screen",
        )
        return Row(None, None, statement)

    def _call(
        self, ctx: RunContext, vid: str, shown: Shown, rate: Rate | None
    ) -> Fields:
        g = ctx.grading
        if rate is None or rate[0] > g.p_low:
            level = Level.LOW
        elif rate[0] <= g.p_high:
            level = Level.HIGH
        else:
            level = Level.MEDIUM
        d = shown.deletion
        return {
            "run_id": ctx.run_id,
            "vid": vid,
            "gt": d.gt,
            "dosage": _dosage(d.gt),
            "p_wrong": None if rate is None else rate[0],
            "p_wrong_lo": None if rate is None else rate[1],
            "p_wrong_hi": None if rate is None else rate[2],
            "level": str(level),
            "gq": None,
            "dp": None,
            "ad": None,
            "vaf": None,
            "row_type": str(RowType.PASS),
            "flags": [],
            "warnings": [
                *warnings(shown.call_class, self.wording.callers),
                *self.cautions(shown),
                *([ZYGOSITY_WARNING] if _dosage(d.gt) == "2" else []),
                CONSEQUENCE_WARNING,
            ],
        }

    def _place(self, ctx: RunContext, gene: str, shown: Shown) -> Placed:
        """The deletion's kind and topic in `gene`, and the line its reason carries for
        that: why the gene gate stopped it, or the gene's override note (A4), which places
        it as the override's class instead of the gate; None where the gate placed it."""
        d = shown.deletion
        inheritance = ctx.inheritance.get(gene)
        if gene in self.placement_overrides:
            cls, note = self.placement_overrides[gene]
            kind, topic, _ = placement(
                cls=cls,
                gene=gene,
                lists=dict(ctx.lists),
                dosage=_dosage(d.gt),
                inheritance=inheritance,
                male=ctx.male,
                haploid=ctx.haploid(d.chrom, d.pos + self.wording.pos_offset),
            )
            return kind, topic, note
        blocked = lof_gene_gate(
            gene,
            mode=gene_mode(gene, ctx.lists, inheritance),
            dosage=_dosage(d.gt),
            filtered=True,
            carrier=carrier_genes(ctx.lists),
            haploinsufficient=self.haploinsufficient,
        )
        if blocked == BOTH_COPIES:
            return Kind.FINDING, Topic.RESEARCH, blocked
        if blocked is not None:
            return Kind.NON_FINDING, None, blocked
        kind, topic, _ = placement(
            cls="pathogenic",
            gene=gene,
            lists=dict(ctx.lists),
            dosage=_dosage(d.gt),
            inheritance=inheritance,
            male=ctx.male,
            haploid=ctx.haploid(d.chrom, d.pos + self.wording.pos_offset),
        )
        return kind, topic, None

    def _gene_rows(self, ctx: RunContext, shown: Shown, seen: set[str]) -> list[Row]:
        d, c, w = shown.deletion, shown.call_class, self.wording
        vid = w.vid(d)
        rate = self.p_wrong(c)
        likely = rate is not None and rate[0] < ctx.grading.deletions.likely_below
        depth = self.depth_support.get((d.chrom, d.pos, d.end))
        place = self.place or self._place
        rows: list[Row] = []
        for gene, exons in shown.exons.items():
            kind, topic, line = place(ctx, gene, shown)
            reason = "; ".join(
                [
                    w.likely if likely else w.possible_for(d),
                    *([] if line is None else [line]),
                    *([FLANK_ONLY] if gene in shown.flank else []),
                    *(
                        []
                        if depth is None
                        else [f"read depth also supports it ({depth})"]
                    ),
                    "predicted, not classified",
                    GENETICIST_NOTE,
                ]
            )
            p = shown.population
            claim = {
                "schema": CLAIM_SCHEMA,
                "type": w.claim_type,
                "gene": gene,
                "chrom": d.chrom,
                "start": d.pos + 1,
                "end": d.end,
                "length": None if d.svlen is None and d.lower_bound else d.length,
                "exons": list(exons),
                "plus_clinical_exons": list(shown.plus_clinical.get(gene, ())),
                "callers": dict(d.callers),
                "depth_ratio": d.depth_ratio,
                "class": {
                    "band": c.band,
                    "evidence": c.evidence,
                    "repeat": c.repeat,
                    "masked": c.masked,
                },
                "gnomad_sv": None if p is None else {"id": p.id, "af": p.af},
                "false_call_rate": None if rate is None else list(rate),
            }
            locus = f"{d.chrom}:{d.pos}-{d.end}"
            if rate is None:
                why = NOT_MEASURED.format(c.band)
                claim[WITHHELD] = [
                    withheld_entry(
                        why,
                        {key: claim[key] for key in _EVENT},
                        would_be=kind,
                        significance=reason,
                        confirm=w.confirm_by,
                    )
                ]
                kind, reason = Kind.UNASSESSED, why
            statement = self._statement(
                ctx,
                f"{self.tool}:{gene}:{locus}",
                kind=str(kind),
                vid=None if rate is None else vid,
                gene=gene,
                claim=json.dumps(claim),
                validity=str(Validity.LIMITED),
                reason=reason,
                topic=None if topic is None else str(topic),
                source_id=locus,
            )
            first = rate is not None and vid not in seen
            seen.add(vid)
            variant: Fields = {"run_id": ctx.run_id, "vid": vid, "rsid": None}
            rows.append(
                Row(
                    variant if first else None,
                    self._call(ctx, vid, shown, rate) if first else None,
                    statement,
                )
            )
        return rows

    def collect(self, ctx: RunContext) -> Collected:
        """The screen's statement, then one per shown deletion and gene."""
        rows = (
            [self._screen_row(ctx, extra=self.screen_claim)] if self.screen_row else []
        )
        if self.not_run is not None:
            return Collected(self.tool, rows, {}, input_rows=len(rows), dropped={})
        seen: set[str] = set()
        for shown in self.shown:
            rows += self._gene_rows(ctx, shown, seen)
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=len(rows) + sum(self.dropped.values()),
            dropped=dict(self.dropped),
        )


def _recall(recall: Mapping[str, Mapping[str, tuple[int, int]]], nouns: str) -> str:
    if not recall:
        return "the screen's recall on HG002 is not recorded"
    bands = " and ".join(
        ", ".join(f"{f} of {n} ({library})" for library, (f, n) in per.items())
        + f" rare {band}"
        for band, per in recall.items()
    )
    return f"on HG002's libraries the screen found {bands} truth {nouns}"


def _dosage(gt: str) -> str:
    return "2" if gt == "1/1" else "1"
