"""Depth copy number as a source (`CnvSource`; slice 3d).

One statement for the screen (`cnv:screen`, 3b's shape): unassessed with its reason when the
depth callers did not run, else a non-finding at validity limited naming how many events
were called and shown, with what the depth screen does not cover, the target exons in 3b's
mask (segmental duplications or low mappability) no Parascopy model covers, and the
large-event bands the calibration withholds in its claim. Then the shown losses and gains through 3b's statements
(`DeletionSource`) in their own words (`LOSS_WORDING`, `GAIN_WORDING`), placed by each
gene's kind (`Shown.kinds`, `cnv.screen.gain_kind`): a loss or a tandem duplication goes
through 1e's gene gate as a loss of function, or the gene's override (A4: a DPYD loss is a
drugs finding); a whole-gene gain is a finding under risks in a triplosensitive gene, else a
non-finding "no known effect of an extra copy"; any other gain a finding under research,
"possible extra copy, effect unknown". Each large event is a statement of its own (no gene,
"needs a chromosomal microarray", shown not calibrated: its call's chance of being false is
not measured), unless the calibration withholds its band; each named locus one statement on
its gene, a carrier finding at its event's class's chance of being false. A shown event or
named locus whose band the table did not measure is withheld (`NOT_MEASURED`). No statement
is a watch item or actionable.
"""

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace

from genokeeper.cnv.screen import DEPTH_CALLERS, LargeEvent, NamedCnv
from genokeeper.deletions.screen import CallClass, Shown
from genokeeper.sources import Collected, Fields, Row, RunContext
from genokeeper.sources.deletions import (
    CONFIRM_BY,
    DELETION_WORDING,
    DELLY_VERSION,
    NOT_MEASURED,
    DeletionSource,
    Placed,
    Rate,
    Wording,
)
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

CNV_TOOL = "cnv"
CNVPYTOR_VERSION = "1.3.2"
SOURCE_LEVEL = f"DELLY {DELLY_VERSION} cnv and CNVpytor {CNVPYTOR_VERSION}"
SOURCE_VERSION = f"DELLY {DELLY_VERSION}; CNVpytor {CNVPYTOR_VERSION}"
LOSS_WORDING = Wording(
    tool=CNV_TOOL,
    noun="copy-number loss",
    claim_type="copy_number_loss",
    confirm_by=CONFIRM_BY,
    # DELLY cnv's own ALT: a loss's vid never takes a 3b deletion's at the same POS (C6).
    alt="<CNV>",
    source_level=SOURCE_LEVEL,
    source_version=SOURCE_VERSION,
    callers=DEPTH_CALLERS,
    pos_offset=1,
)
GAIN_WORDING = replace(
    LOSS_WORDING, noun="copy-number gain", claim_type="copy_number_gain", alt="<DUP>"
)
# The `cnv:screen` statement's words: losses and gains together.
SCREEN_WORDING = replace(
    LOSS_WORDING, noun="copy-number change", claim_type="copy_number"
)
LARGE_REASON = "screening result, needs a chromosomal microarray"
LARGE_NOT_CALIBRATED = "shown, not calibrated: no truth set holds a large CNV"
TANDEM_NOTE = "a tandem duplication inside the gene may disrupt it"
GAIN_DEPTH_ONLY = "possible extra copy, effect unknown"
WHOLE_GAIN = "no known effect of an extra copy"
# The spec's "Not covered, named" for depth copy number, in the `cnv:screen` claim.
NOT_COVERED_DEPTH = (
    "depth events under about 1 kb (most single-exon events)",
    "target exons in segmental duplications or low mappability without a Parascopy model",
    "GSTT1 (off the primary assembly)",
    "the location of an extra whole-gene copy",
)
# Drop reasons the source adds to the screen's.
LARGE_WITHHELD = "a large event in a withheld band"
NAMED_AGAIN = "a named locus already stated"


@dataclass(frozen=True)
class CnvSource:
    """The depth screen's shown events (`cnv.screen.CnvScreen`), its large events and
    named loci (3b's deletions over a named locus included, `named_deletions`, the first
    per locus stated), how many events the depth callers made (`called`), the drop counts,
    why the screen did not run (`not_run`), each class's chance of being false (rate, lo,
    hi; None where the table has none; a named 3b deletion's class is 3b's), ClinGen's
    haploinsufficient and triplosensitive genes, the per-gene placement overrides
    (`genelists.screen_overrides`, A4), the target exons in 3b's mask
    (`segdup_exons`), the large-event bands the calibration withholds, band -> its
    reason, and what the screen's claim also says of how it ran (`notes`: CNVpytor dropped
    by the calibration)."""

    shown: Sequence[Shown]
    large: Sequence[LargeEvent]
    named: Sequence[tuple[NamedCnv, Shown]]
    called: int
    dropped: Mapping[str, int]
    not_run: str | None
    p_wrong: Callable[[CallClass], Rate | None]
    haploinsufficient: frozenset[str]
    triplosensitive: frozenset[str]
    placement_overrides: Mapping[str, tuple[str, str]] = field(
        default_factory=dict[str, tuple[str, str]]
    )
    segdup_exons: Mapping[str, tuple[int, ...]] = field(
        default_factory=dict[str, tuple[int, ...]]
    )
    large_withheld: Mapping[str, str] = field(default_factory=dict[str, str])
    notes: Sequence[str] = ()

    @property
    def tool(self) -> str:
        return CNV_TOOL

    def _rows_by(self, wording: Wording) -> DeletionSource:
        """3b's source in `wording`, writing no screen row of its own, placed by kind."""
        return DeletionSource(
            shown=self.shown,
            called=self.called,
            dropped=self.dropped,
            not_run=self.not_run,
            p_wrong=self.p_wrong,
            haploinsufficient=self.haploinsufficient,
            recall={},
            wording=wording,
            place=self._place,
            placement_overrides=self.placement_overrides,
            screen_row=False,
        )

    def _place(self, ctx: RunContext, gene: str, shown: Shown) -> Placed:
        """The event's kind and topic in `gene` by its kind there (module docstring)."""
        kind = shown.kinds[gene]
        if kind in ("loss", "tandem"):
            # 3b's own placement: the gene's override, else the gene gate.
            placed, topic, line = self._rows_by(LOSS_WORDING)._place(ctx, gene, shown)
            if kind == "tandem":
                line = "; ".join([TANDEM_NOTE, *([] if line is None else [line])])
            return placed, topic, line
        if kind == "whole":
            if gene in self.triplosensitive:
                return Kind.FINDING, Topic.RISKS, None
            return Kind.NON_FINDING, None, WHOLE_GAIN
        return Kind.FINDING, Topic.RESEARCH, GAIN_DEPTH_ONLY

    def _own_row(
        self,
        ctx: RunContext,
        statement_id: str,
        shown: Shown,
        rate: Rate | None,
        seen: set[str],
        **fields: object,
    ) -> Row:
        """A statement of the source's own (a large event, a named locus) over `shown`'s
        event, with its `variant` and `call` rows where it is the first at its vid."""
        d = shown.deletion
        if not set(d.callers) <= set(DEPTH_CALLERS):
            # A named 3b deletion: 3b's callers, source and warnings, under this tool.
            wording = replace(DELETION_WORDING, tool=CNV_TOOL)
        else:
            wording = LOSS_WORDING if d.svtype == "DEL" else GAIN_WORDING
        rows = self._rows_by(wording)
        vid = wording.vid(d)
        withheld = fields["kind"] == str(Kind.UNASSESSED)
        locus = f"{d.chrom}:{d.pos}-{d.end}"
        statement = rows._statement(
            ctx,
            statement_id,
            vid=None if withheld else vid,
            validity=str(Validity.LIMITED),
            source_id=locus,
            **fields,
        )
        first = not withheld and vid not in seen
        seen.add(vid)
        variant: Fields = {"run_id": ctx.run_id, "vid": vid, "rsid": None}
        return Row(
            variant if first else None,
            rows._call(ctx, vid, shown, rate) if first else None,
            statement,
        )

    def _large_row(self, ctx: RunContext, event: LargeEvent, seen: set[str]) -> Row:
        d, c = event.event, event.call_class
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "large_copy_number",
            "chrom": d.chrom,
            "start": d.pos + 1,
            "end": d.end,
            "length": d.length,
            "direction": d.svtype,
            "callers": dict(d.callers),
            "class": {"band": c.band, "evidence": c.evidence, "masked": c.masked},
            "regions": list(event.regions),
        }
        return self._own_row(
            ctx,
            f"{CNV_TOOL}:large:{d.chrom}:{d.pos}-{d.end}",
            Shown(d, c, {}, None),
            None,
            seen,
            kind=str(Kind.FINDING),
            gene=None,
            claim=json.dumps(claim),
            reason=f"{LARGE_REASON}; {LARGE_NOT_CALIBRATED}",
            topic=str(Topic.RISKS),
        )

    def _named_row(
        self, ctx: RunContext, named: NamedCnv, shown: Shown, seen: set[str]
    ) -> Row:
        d, c = shown.deletion, shown.call_class
        rate = self.p_wrong(c)
        # Both copies: DFNB1 itself, a risk (controller ruling, 2026-10-30).
        reason = named.homozygous_reason if d.gt == "1/1" else named.reason
        kind = Kind.FINDING
        if rate is None:
            kind, reason = Kind.UNASSESSED, NOT_MEASURED.format(c.band)
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "named_copy_number",
            "name": named.name,
            "gene": named.gene,
            "chrom": d.chrom,
            "start": d.pos + 1,
            "end": d.end,
            "callers": dict(d.callers),
            "class": {"band": c.band, "evidence": c.evidence, "masked": c.masked},
            "false_call_rate": None if rate is None else list(rate),
        }
        return self._own_row(
            ctx,
            f"{CNV_TOOL}:named:{named.name}",
            shown,
            rate,
            seen,
            kind=str(kind),
            gene=named.gene,
            claim=json.dumps(claim),
            reason=reason,
            topic=str(Topic.RISKS if d.gt == "1/1" else Topic.CARRIER),
        )

    def collect(self, ctx: RunContext) -> Collected:
        """The screen's statement, then one per shown event and gene, per large event in
        a band not withheld, and per named locus."""
        stated: dict[str, tuple[NamedCnv, Shown]] = {}
        for named, item in self.named:
            stated.setdefault(named.name, (named, item))
        large = [e for e in self.large if e.call_class.band not in self.large_withheld]
        own = {
            LARGE_WITHHELD: len(self.large) - len(large),
            NAMED_AGAIN: len(self.named) - len(stated),
        }
        dropped = {**self.dropped, **{reason: n for reason, n in own.items() if n}}
        extra = {
            "large_withheld": dict(self.large_withheld),
            "segdup_exons": {g: list(n) for g, n in self.segdup_exons.items()},
            "not_covered": list(NOT_COVERED_DEPTH),
            "notes": list(self.notes),
        }
        screen = self._rows_by(SCREEN_WORDING)
        shown = len(self.shown) + len(large) + len(stated)
        rows = [screen._screen_row(ctx, shown=shown, extra=extra)]
        if self.not_run is not None:
            return Collected(self.tool, rows, {}, input_rows=len(rows), dropped={})
        seen: set[str] = set()
        losses, gains = self._rows_by(LOSS_WORDING), self._rows_by(GAIN_WORDING)
        for item in self.shown:
            source = gains if item.deletion.svtype == "DUP" else losses
            rows += source._gene_rows(ctx, item, seen)
        rows += [self._large_row(ctx, event, seen) for event in large]
        rows += [self._named_row(ctx, n, item, seen) for n, item in stated.values()]
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=len(rows) + sum(dropped.values()),
            dropped=dropped,
        )
