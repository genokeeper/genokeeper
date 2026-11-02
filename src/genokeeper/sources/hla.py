"""HLA typing as a source (`HlaSource`; slice 3f, "HLA"). Tool `hla`; nothing is actionable,
nothing a watch item, and only the coeliac rule-out is flagged `exclusion`.

- One statement per row of `resources/hla/alleles.tsv` (`hla:<gene><allele>`; a research
  association `hla:<gene><allele>:association`, ruling R3, so that C*06:02's drug statement
  and its psoriasis association differ; its source `hla_association` is a research layer,
  ruling R19). A class I allele needs OptiType and T1K to agree:
  both positive a finding, both negative a non-finding, anything else unassessed naming
  both (`hla.status`, which reads a one-field group such as `*27` too). A class II row reads
  T1K alone and says so (`CLASS_II_ONE_CALLER`).
- `hla:coeliac`: T1K's DQA1/DQB1 and the six DQ tag SNPs (`hla_typing.coeliac`); the
  rule-out is a non-finding flagged `exclusion` whose claim names the confirming test.
- `hla:table`: the class I genes both callers agree on and T1K's classical class II genes,
  a research finding.
- `hla:KIR`: always unassessed, not covered.

A gene in `withheld` (the calibration table's, keyed by gene, `DQ` for the tag SNPs) makes
its statements unassessed with the table's reason; `not_run` makes every statement but
`hla:KIR` unassessed with it. A caller that made no call for a gene is unassessed, never
negative.
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.pgx.hla import (
    CLASS_I,
    IMGT_HLA_VERSION,
    OPTITYPE,
    T1K,
    HlaCall,
    shown,
    status,
)
from genokeeper.pgx.hla_typing import (
    CONFIRM_DQ,
    TAGS_EUROPEAN,
    DqStatus,
    HlaAllele,
    coeliac,
    dq_from_t1k,
    dq_from_tags,
    load_hg002_truth,
    load_hla_alleles,
)
from genokeeper.sources import INFORMATION_NOTE, Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

HLA_TOOL = "hla"
HLA_ASSOCIATION_TOOL = "hla_association"
KIR_REASON = (
    "not covered: KIR typing has no use for an adult and the tools disagree"
    " (Castelli et al., HLA 2025)"
)
CLASS_II_ONE_CALLER = "class II typed by T1K alone"
CALLERS_DISAGREE = "the callers do not agree"
CLASS_II = (
    "HLA-DRB1",
    "HLA-DRB3",
    "HLA-DRB4",
    "HLA-DRB5",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DPA1",
    "HLA-DPB1",
)
# The keys of `withheld` that hold back the coeliac statement.
DQ_WITHHELD = ("HLA-DQA1", "HLA-DQB1", "DQ")
MONSUUR = "Monsuur et al. 2008, doi:10.1371/journal.pone.0002270"
VERSION = f"IPD-IMGT/HLA {IMGT_HLA_VERSION}"


def _alleles(call: HlaCall | None) -> list[str] | None:
    """A call's two alleles in order; None when it has none or one was dropped."""
    if call is None or call[0] is None or call[1] is None:
        return None
    return sorted((call[0], call[1]))


def _dq(found: DqStatus | None) -> dict[str, object] | None:
    if found is None:
        return None
    return {
        "dq25": found.dq25,
        "dq22": found.dq22,
        "dq8": found.dq8,
        "half": found.half,
    }


@dataclass(frozen=True)
class HlaSource:
    """OptiType's and T1K's calls (None: the caller gave none; T1K's from
    `read_t1k(..., genes=None)`), the six tag SNPs' dosages (rsid -> `0`, `1`, `2` or
    `unknown`), why HLA was not run, and the calibration's withheld genes."""

    optitype: Mapping[str, HlaCall] | None
    t1k: Mapping[str, HlaCall] | None
    tags: Mapping[str, str]
    not_run: str | None
    withheld: Mapping[str, str]
    tool: str = HLA_TOOL

    def _calls(self, gene: str) -> dict[str, HlaCall | None]:
        callers: Sequence[tuple[str, Mapping[str, HlaCall] | None]] = (
            [(T1K, self.t1k)]
            if gene not in CLASS_I
            else [(OPTITYPE, self.optitype), (T1K, self.t1k)]
        )
        return {name: None if m is None else m.get(gene) for name, m in callers}

    def _allele(self, ctx: RunContext, a: HlaAllele) -> Row:
        calls = self._calls(a.gene)
        statuses = {name: status(call, a.allele) for name, call in calls.items()}
        silent = [name for name, call in calls.items() if call is None]
        values = set(statuses.values())
        one_caller = a.gene not in CLASS_I
        kind: Kind
        if self.not_run is not None:
            kind, reason = Kind.UNASSESSED, self.not_run
        elif a.gene in self.withheld:
            kind, reason = Kind.UNASSESSED, self.withheld[a.gene]
        elif silent:
            kind, reason = Kind.UNASSESSED, f"{' and '.join(silent)} made no call"
        elif values == {"positive"}:
            kind = Kind.FINDING
            reason = f"{a.note}; {CLASS_II_ONE_CALLER}" if one_caller else a.note
        elif values == {"negative"}:
            kind, reason = Kind.NON_FINDING, CLASS_II_ONE_CALLER if one_caller else None
        else:
            read = ", ".join(f"{name} {s}" for name, s in statuses.items())
            reason = (
                f"{read}: an allele dropped or equivalent alleles"
                if one_caller
                else f"{CALLERS_DISAGREE} ({read})"
            )
            kind = Kind.UNASSESSED
        research = a.topic is Topic.RESEARCH
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "hla_association" if research else "hla_allele",
            "gene": a.gene,
            "allele": a.allele,
            "condition": a.condition,
            "drugs": list(a.drugs),
            "callers": {name: shown(call) for name, call in calls.items()},
            "status": statuses,
        }
        sid = f"{HLA_TOOL}:{a.gene}{a.allele}" + (":association" if research else "")
        return self._row(
            ctx,
            sid,
            kind=kind,
            gene=a.gene,
            claim=claim,
            validity=a.validity,
            level=a.level,
            reason=reason,
            topic=a.topic,
            source_id=a.source,
            source=HLA_ASSOCIATION_TOOL if research else HLA_TOOL,
        )

    def _coeliac(self, ctx: RunContext) -> Row:
        t1k = None if self.t1k is None else dq_from_t1k(self.t1k)
        tags = dq_from_tags(self.tags)
        withheld = next(
            (self.withheld[k] for k in DQ_WITHHELD if k in self.withheld), None
        )
        if self.not_run is not None:
            kind, reason = Kind.UNASSESSED, self.not_run
        elif withheld is not None:
            kind, reason = Kind.UNASSESSED, withheld
        elif self.t1k is None:
            kind, reason = Kind.UNASSESSED, f"{T1K} made no call"
        else:
            kind, reason = coeliac(t1k, tags)
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "hla_coeliac",
            "condition": "coeliac disease",
            "t1k": _dq(t1k),
            "tags": _dq(tags),
            "tag_dosages": dict(self.tags),
            "note": TAGS_EUROPEAN,
            "confirm": CONFIRM_DQ,
        }
        return self._row(
            ctx,
            f"{HLA_TOOL}:coeliac",
            kind=kind,
            gene="HLA-DQ",
            claim=claim,
            validity=Validity.LIMITED,
            level=f"{T1K} and the DQ tag SNPs",
            reason=reason,
            topic=Topic.RISKS,
            source_id=MONSUUR,
            exclusion=kind is Kind.NON_FINDING,
        )

    def _table(self, ctx: RunContext) -> Row:
        genes: list[dict[str, object]] = []
        for gene in CLASS_I + CLASS_II:
            calls = self._calls(gene)
            if all(call is None for call in calls.values()):
                continue
            entry: dict[str, object] = {"gene": gene}
            if gene in self.withheld:
                entry |= {"alleles": None, "callers": None, "note": self.withheld[gene]}
            elif gene in CLASS_I:
                o, t = calls[OPTITYPE], calls[T1K]
                agreed = _alleles(o)
                if agreed is not None and agreed == _alleles(t):
                    entry |= {"alleles": agreed, "callers": f"{OPTITYPE} and {T1K}"}
                else:
                    entry |= {
                        "alleles": None,
                        "callers": f"{OPTITYPE} {shown(o)}, {T1K} {shown(t)}",
                        "note": CALLERS_DISAGREE,
                    }
            else:
                t = calls[T1K]
                assert t is not None
                entry |= {
                    "alleles": sorted(a or "?" for a in t),
                    "callers": T1K,
                    "note": CLASS_II_ONE_CALLER,
                }
            genes.append(entry)
        if self.not_run is not None:
            kind, reason = Kind.UNASSESSED, self.not_run
        elif self.t1k is None:
            kind, reason = Kind.UNASSESSED, f"{T1K} made no call"
        else:
            # The genes no calibration measure covers (HG002's truth types none).
            other = [g["gene"] for g in genes if g["gene"] not in load_hg002_truth()]
            note = (
                f"; {', '.join(map(str, other))}: {INFORMATION_NOTE}" if other else ""
            )
            kind, reason = Kind.FINDING, CLASS_II_ONE_CALLER + note
        claim = {"schema": CLAIM_SCHEMA, "type": "hla_table", "genes": genes}
        return self._row(
            ctx,
            f"{HLA_TOOL}:table",
            kind=kind,
            gene=None,
            claim=claim,
            validity=Validity.LIMITED,
            level=f"{OPTITYPE} and {T1K}",
            reason=reason,
            topic=Topic.RESEARCH,
            source_id=VERSION,
        )

    def _kir(self, ctx: RunContext) -> Row:
        claim: dict[str, object] = {
            "schema": CLAIM_SCHEMA,
            "type": "hla_kir",
            "gene": "KIR",
        }
        return self._row(
            ctx,
            f"{HLA_TOOL}:KIR",
            kind=Kind.UNASSESSED,
            gene=None,
            claim=claim,
            validity=Validity.LIMITED,
            level="not covered",
            reason=KIR_REASON,
            topic=Topic.RESEARCH,
            source_id="KIR",
        )

    def _row(
        self,
        ctx: RunContext,
        statement_id: str,
        *,
        kind: Kind,
        gene: str | None,
        claim: dict[str, object],
        validity: Validity,
        level: str,
        reason: str | None,
        topic: Topic,
        source_id: str,
        exclusion: bool = False,
        source: str = HLA_TOOL,
    ) -> Row:
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": statement_id,
            "kind": str(kind),
            "vid": None,
            "gene": gene,
            "claim": json.dumps(claim, ensure_ascii=False),
            "validity": str(validity),
            "source_level": level,
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(topic),
            "source": source,
            "source_id": source_id,
            "source_version": VERSION,
            "points": None,
            "exclusion": exclusion,
        }
        return Row(None, None, statement)

    def collect(self, ctx: RunContext) -> Collected:
        """One statement per allele row, then the coeliac statement, the table and KIR."""
        rows = [self._allele(ctx, a) for a in load_hla_alleles()]
        rows += [self._coeliac(ctx), self._table(ctx), self._kir(ctx)]
        return Collected(
            tool=self.tool, rows=rows, extra={}, input_rows=len(rows), dropped={}
        )
