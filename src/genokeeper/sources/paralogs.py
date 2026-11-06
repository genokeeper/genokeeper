"""The paralog loci and Lp(a) as sources (`ParalogSource`, `LpaSource`; slice 3d).

`ParalogSource` writes per locus of `paralogs.LOCI` one statement (`paralogs:<locus>`) and
one `gene_call` (3a's SMN shape: the copy numbers, every caller's call, whether two
agreed, the level). STRC is the one paralog carrier locus (A3): one copy that passed at
`[paralogs] min_qual` with CATSPER2's unique depth inside `unique_loss_ratio` is a carrier
finding, "screening result, needs confirmation"; two or more a non-finding; anything else,
no copy included, unassessed naming why. RHD, AMY1, C4 (the total) and CYP2A6 are
information: a finding at validity limited under traits (CYP2A6 under drugs, and only when
Aldy gives the same copies), never actionable; RHD is a copy number, never a blood type
(design decision). GBA1 is graded from Gauchian in exactly four shapes (C21): wild type a non-finding;
a heterozygous carrier whose GBA+GBAP1 copies match Parascopy's agCN over exons 9-11, with
`gba_min_reads` reads behind its variant, a carrier finding; biallelic withheld; Parascopy
missing or disagreeing, or too few reads, unassessed; any other output unassessed "Gauchian
output not handled". A locus the calibration withholds, or every locus when the callers did
not run, is unassessed with that reason. Then one unassessed statement per gene the slice
names as not covered (`NOT_COVERED`). `gba1_contradicted` tells which ClinVar or predicted
statement in GBA1's exons 9-11 Gauchian contradicts (applied by the run, not here).

`LpaSource` writes one statement, `lpa:LPA`: KILDA's KIV-2 count as information under
wellness, a genetic estimate beside the blood test that measures what matters.

No statement here is a watch item or actionable.
"""

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass

from genokeeper.config import ParalogThresholds
from genokeeper.evidence.pgx import GENE_UNRESOLVED
from genokeeper.paralogs import (
    GAUCHIAN_VERSION,
    KILDA_VERSION,
    LOCI,
    GauchianCall,
    Locus,
    PsCall,
    wild_type,
)
from genokeeper.pgx.cyp2d6 import ALDY_VERSION
from genokeeper.smn import PARASCOPY_VERSION
from genokeeper.sources import Collected, Fields, Row, RunContext, present
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.sources.predicted import PREDICTED_TOOL
from genokeeper.sources.repeats import CONFIRM
from genokeeper.store.db import Kind, Level, Topic, Validity

PARALOGS_TOOL = "paralogs"
LPA_TOOL = "lpa"
PARASCOPY_NAME = f"Parascopy {PARASCOPY_VERSION}"
GAUCHIAN_NAME = f"Gauchian {GAUCHIAN_VERSION}"
ALDY_NAME = f"Aldy {ALDY_VERSION}"
KILDA_NAME = f"KILDA {KILDA_VERSION}"
SOURCE_LEVEL = f"{PARASCOPY_NAME}, {GAUCHIAN_NAME} and {ALDY_NAME}"
SOURCE_VERSION = f"{PARASCOPY_VERSION}; {GAUCHIAN_VERSION}; {ALDY_VERSION}"

STRC_REASON = (
    f"{CONFIRM} by MLPA or long reads: one STRC copy, carrier of DFNB16 hearing loss"
)
GBA_REASON = (
    f"{CONFIRM} by long-range PCR and sequencing: carrier for Gaucher disease; one GBA1"
    " variant also raises the risk of Parkinson's disease"
)
# The spec (Tayebi et al. 2025).
BIALLELIC_WITHHELD = "biallelic calls are withheld: disputed in an independent check"
# The spec (2025 ESC/EAS and 2026 ACC/AHA).
LPA_REASON = "a genetic estimate: a one-off Lp(a) blood test measures what matters"
RHD_NOTE = (
    "zero copies usually means RhD-negative in European ancestry; serology decides"
)
ZERO_REASON = "no copy in a healthy adult: a calling error"
UNHANDLED = "Gauchian output not handled"
GBA1_CONTRADICTED = (
    "GBA1 exons 9–11: Gauchian calls no variant at this position (it tests only its"
    " pseudogene-like sites); a pseudogene read may have mismapped"
)
_INFORMATION: dict[str, tuple[Topic, str]] = {
    "RHD": (Topic.TRAITS, RHD_NOTE),
    "AMY1": (Topic.TRAITS, "information: AMY1A, AMY1B and AMY1C copies together"),
    "C4": (
        Topic.TRAITS,
        "information: C4A and C4B copies together, which cannot be told apart reliably",
    ),
    "CYP2A6": (Topic.DRUGS, "no CPIC guideline; information"),
}
_CONDITIONS = {"STRC": "DFNB16 hearing loss", "GBA1": "Gaucher disease"}
_NCF1_SORD = (
    "their common carrier alleles, the GT deletion and c.757delG, change no copy number,"
    " so a copy-number rule cannot see them"
)
_NO_CALLER = "no open short-read caller"
_HBA = (
    "-α3.7 and -α4.2 not covered: no open caller; deletions removing both genes are in"
    " the deletion screen"
)
_PMS2 = (
    "exons 12-15 not covered: Parascopy cannot tell PMS2 from PMS2CL, and HG002 has no"
    " truth there"
)
# The genes the slice names as not covered (the spec's "Not covered, named", A2, A3 and the
# 2026-10-30 pre-flight), each an unassessed gene-level statement.
NOT_COVERED: dict[str, tuple[Topic, str]] = {
    **dict.fromkeys(("HBA1", "HBA2"), (Topic.CARRIER, _HBA)),
    **dict.fromkeys(("CYP21A2", "TNXB"), (Topic.CARRIER, _NO_CALLER)),
    "PMS2": (Topic.RISKS, _PMS2),
    "F8": (
        Topic.RISKS,
        "intron 1 and 22 inversions not covered: no short-read method",
    ),
    "OTOA": (
        Topic.CARRIER,
        "GRCh38 lacks a copy, so its model's reference copy number is wrong",
    ),
    "GSTT1": (Topic.TRAITS, "not on GRCh38's primary assembly"),
    # A2.
    "GSTM1": (
        Topic.TRAITS,
        "Parascopy's GSTM5 model covers 974 bp, 18%, of GSTM1's MANE span",
    ),
    **dict.fromkeys(("NCF1", "SORD"), (Topic.CARRIER, _NCF1_SORD)),
    "UGT2B17": (
        Topic.TRAITS,
        "Parascopy's UGT2B15 model covers 21% of UGT2B17's MANE span",
    ),
}
# Gauchian's named variants inside GBA1's exons 9-11 span (`GBA1_ex9-11`) -> their GRCh38
# chr1 position (1-based): the rows of gauchian/data/GBA_target_variant_38.txt and
# GBA_target_variant_homology_region_38.txt there, Gauchian v1.0.2 (commit
# 371a6d79599389e590038c4fecf13ec853abb8ab), read from the built image on 2026-10-30.
GAUCHIAN_SITES: dict[str, int] = {
    "c.1604G>A(p.Arg535His)": 155235002,
    "c.1603C>T(p.Arg535Cys)": 155235003,
    "c.1599G>A(p.Trp533Ter)": 155235007,
    "c.1549G>A(p.Gly517Ser)": 155235057,
    "c.1506-1G>A": 155235101,
    "c.1504C>T(p.Arg502Cys)": 155235196,
    "c.1503C>G(p.Asn501Lys)": 155235197,
    "c.1495G>C(p.Val499Leu)": 155235205,
    "c.1459G>A(p.Ala487Thr)": 155235241,
    "c.1448T>G(p.Leu483Arg)": 155235252,
    "c.1397T>G(p.Ile466Ser)": 155235303,
    "c.1388+2T>C": 155235679,
    "c.1361C>G(p.Pro454Arg)": 155235708,
    "c.1357C>T(p.Gln453Ter)": 155235712,
    "c.1348T>A(p.Phe450Ile)": 155235721,
    "c.1343A>T(p.Asp448Val)": 155235726,
    "c.1319C>T(p.Pro440Leu)": 155235750,
    "c.1312G>A(p.Asp438Asn)": 155235757,
    "c.1309G>T(p.Val437Phe)": 155235760,
    "c.1297G>T(p.Val433Leu)": 155235772,
    "c.1289C>T(p.Pro430Leu)": 155235780,
    "c.1279G>T(p.Glu427Ter)": 155235790,
    "c.1271T>C(p.Leu424Pro)": 155235798,
    "c.1259G>A(p.Trp420Ter)": 155235810,
    "c.1250G>A(p.Trp417Ter)": 155235819,
    "c.1246G>A(p.Gly416Ser)": 155235823,
    "c.1240G>T(p.Val414Leu)": 155235829,
    "c.1238A>C(p.His413Pro)": 155235831,
    "c.1228C>G(p.Leu410Val)": 155235841,
    "c.1226A>C(p.Asn409Thr)": 155235843,
    "c.1226A>G(p.Asn409Ser)": 155235843,
}
_GRADED = ("clinvar", PREDICTED_TOOL)

type Graded = tuple[Kind, str | None]


def _shown(value: object) -> str:
    return "?" if value is None else str(value)


def _usable(call: PsCall | None, locus: str, t: ParalogThresholds) -> int | str:
    """The locus's copy number from Parascopy's call, or why there is none."""
    why = None
    if call is None or call.rows == 0:
        why = f"no row over {locus}"
    elif not call.passed:
        why = "a row failed its filter"
    elif call.cn is None:
        why = f"its {call.rows} rows over {locus} disagree or one is unknown"
    elif call.qual is None or call.qual < t.min_qual:
        why = f"quality {_shown(call.qual)} below {t.min_qual:g}"
    if why is not None or call is None or call.cn is None:
        return f"{PARASCOPY_NAME} gave no usable {locus} copy number ({why})"
    return call.cn


def _fields(g: GauchianCall) -> str:
    """Gauchian's call, field by field, for an output no shape handles."""
    return (
        f"CN(GBA+GBAP1) {_shown(g.total_cn)}, carrier {_shown(g.carrier)},"
        f" biallelic {_shown(g.biallelic)}, variants {', '.join(g.variants) or 'none'},"
        f" other {', '.join(g.other) or 'none'}, breakpoint {g.breakpoint or 'none'}"
    )


@dataclass(frozen=True)
class ParalogSource:
    """Parascopy's call per locus of `LOCI` (a locus absent: no call), CATSPER2's unique
    depth over the run's mean (`paralogs.depth_ratio`; None: not measured), Gauchian's GBA1
    call, Aldy's CYP2A6 copies (`paralogs.aldy_copies`), the loci the calibration withholds
    (locus -> its reason) and why the callers did not run (`not_run`). The thresholds are
    the run's `[paralogs]`."""

    calls: Mapping[str, PsCall]
    unique_ratio: float | None
    gauchian: GauchianCall | None
    aldy_cyp2a6: int | None
    withheld: Mapping[str, str]
    not_run: str | None

    @property
    def tool(self) -> str:
        return PARALOGS_TOOL

    def _strc(self, t: ParalogThresholds) -> Graded:
        cn = _usable(self.calls.get("STRC"), "STRC", t)
        if isinstance(cn, str):
            return Kind.UNASSESSED, cn
        if cn == 0:
            return Kind.UNASSESSED, ZERO_REASON
        if cn >= 2:
            return Kind.NON_FINDING, None
        low, high = t.unique_loss_ratio
        ratio = self.unique_ratio
        if ratio is None:
            return (
                Kind.UNASSESSED,
                "one STRC copy, but CATSPER2's depth was not measured",
            )
        if not low <= ratio <= high:
            return Kind.UNASSESSED, (
                f"one STRC copy, but CATSPER2's unique depth is {ratio:.2f} of the mean,"
                f" outside {low:g}-{high:g}, where a one-copy loss puts it"
            )
        return Kind.FINDING, STRC_REASON

    def _information(self, locus: str, t: ParalogThresholds) -> Graded:
        cn = _usable(self.calls.get(locus), locus, t)
        note = _INFORMATION[locus][1]
        if isinstance(cn, str):
            return Kind.UNASSESSED, cn
        if locus == "CYP2A6" and self.aldy_cyp2a6 != cn:
            if self.aldy_cyp2a6 is None:
                return Kind.UNASSESSED, f"{ALDY_NAME} gave no CYP2A6 copy number"
            return (
                Kind.UNASSESSED,
                f"Parascopy and Aldy disagree ({cn} and {self.aldy_cyp2a6} copies)",
            )
        return Kind.FINDING, note

    def _gba1(self, t: ParalogThresholds) -> Graded:
        """Gauchian's call in the four shapes the module docstring names, else unhandled."""
        g = self.gauchian
        if g is None:
            return Kind.UNASSESSED, f"{GAUCHIAN_NAME} gave no GBA1 call"
        if wild_type(g):
            return Kind.NON_FINDING, None
        if g.biallelic is True:
            return Kind.UNASSESSED, BIALLELIC_WITHHELD
        carrier = g.carrier is True and g.biallelic is False and g.variants
        if not carrier or g.other or g.total_cn is None or g.breakpoint is not None:
            return Kind.UNASSESSED, f"{UNHANDLED}: {_fields(g)}"
        ag = _usable(self.calls.get("GBA1"), "GBA+GBAP1", t)
        if isinstance(ag, str):
            return Kind.UNASSESSED, ag
        if g.total_cn != ag:
            return Kind.UNASSESSED, (
                f"Gauchian and Parascopy disagree on GBA+GBAP1 copies"
                f" ({_shown(g.total_cn)} and {ag})"
            )
        if g.reads is None:
            return (
                Kind.UNASSESSED,
                "Gauchian's carrier call gives no count of supporting reads",
            )
        if g.reads < t.gba_min_reads:
            return Kind.UNASSESSED, (
                f"Gauchian's carrier call has {g.reads} supporting reads,"
                f" below {t.gba_min_reads}"
            )
        return Kind.FINDING, GBA_REASON

    def _grade(self, locus: Locus, t: ParalogThresholds) -> Graded:
        """The locus's kind and reason; every RHD statement carries `RHD_NOTE` (design decision)."""
        kind, reason = self._graded(locus, t)
        if locus.locus == "RHD" and reason != RHD_NOTE:
            reason = f"{reason}; {RHD_NOTE}"
        return kind, reason

    def _graded(self, locus: Locus, t: ParalogThresholds) -> Graded:
        if self.not_run is not None:
            return Kind.UNASSESSED, self.not_run
        if locus.locus in self.withheld:
            return Kind.UNASSESSED, self.withheld[locus.locus]
        if locus.locus == "STRC":
            return self._strc(t)
        if locus.locus == "GBA1":
            return self._gba1(t)
        return self._information(locus.locus, t)

    def _callers(self, locus: Locus) -> dict[str, str]:
        """Every caller's own call of the locus."""
        callers: dict[str, str] = {}
        call = self.calls.get(locus.locus)
        if call is not None:
            callers[PARASCOPY_NAME] = (
                f"{locus.measure} {_shown(call.cn)}, qual {_shown(call.qual)},"
                f" {'pass' if call.passed else 'fail'}, {call.rows} rows"
            )
        if locus.unique is not None and self.unique_ratio is not None:
            callers[f"{locus.unique} depth"] = f"{self.unique_ratio:.2f} of the mean"
        if locus.locus == "CYP2A6":
            callers[ALDY_NAME] = f"{_shown(self.aldy_cyp2a6)} copies"
        if locus.locus == "GBA1" and self.gauchian is not None:
            g = self.gauchian
            callers[GAUCHIAN_NAME] = f"{_fields(g)}, reads {_shown(g.reads)}"
        return callers

    def _copies(self, locus: Locus) -> int | None:
        """The copy number the statement shows: Gauchian's total for GBA1."""
        if locus.locus == "GBA1":
            return None if self.gauchian is None else self.gauchian.total_cn
        call = self.calls.get(locus.locus)
        return None if call is None else call.cn

    def _measured(self, locus: Locus) -> Fields | None:
        """Each caller's own call of the locus."""
        call = self.calls.get(locus.locus)
        g = self.gauchian if locus.locus == "GBA1" else None
        aldy = self.aldy_cyp2a6 if locus.locus == "CYP2A6" else None
        return present(
            {
                PARASCOPY_NAME: call and present({"copies": call.cn}),
                GAUCHIAN_NAME: g and present(asdict(g)),
                ALDY_NAME: present({"copies": aldy}),
            }
        )

    def _agreed(self, locus: Locus) -> bool | None:
        """Whether the two callers agree (None with one, or one without a value)."""
        call = self.calls.get(locus.locus)
        cn = None if call is None else call.cn
        other = None
        if locus.locus == "CYP2A6":
            other = self.aldy_cyp2a6
        elif locus.locus == "GBA1" and self.gauchian is not None:
            other = self.gauchian.total_cn
        return None if cn is None or other is None else cn == other

    def _rows(self, ctx: RunContext, locus: Locus) -> tuple[Row, Fields]:
        kind, reason = self._grade(locus, ctx.grading.paralogs)
        assessed = kind is not Kind.UNASSESSED
        copies = self._copies(locus) if assessed else None
        callers = self._callers(locus)
        carrier = locus.role == "carrier"
        topic = Topic.CARRIER if carrier else _INFORMATION[locus.locus][0]
        claim: dict[str, object] = {
            "schema": CLAIM_SCHEMA,
            "type": "paralog_copy_number",
            "locus": locus.locus,
            "genes": list(locus.genes),
            "measure": locus.measure,
            "copy_number": copies,
            "callers": callers,
        }
        if locus.locus in _CONDITIONS:
            claim["condition"] = _CONDITIONS[locus.locus]
        if locus.locus == "GBA1" and self.gauchian is not None:
            claim["variants"] = list(self.gauchian.variants)
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{PARALOGS_TOOL}:{locus.locus}",
            "kind": str(kind),
            "vid": None,
            "gene": locus.locus,
            "claim": json.dumps(claim),
            # Every 3d result is at validity limited: Task 16's "possible carrier".
            "validity": str(Validity.LIMITED),
            "source_level": SOURCE_LEVEL,
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(topic),
            "source": PARALOGS_TOOL,
            "source_id": locus.locus,
            "source_version": SOURCE_VERSION,
            "points": None,
        }
        phenotype = None
        if carrier and assessed:
            phenotype = (
                f"{_CONDITIONS[locus.locus]} carrier (screening)"
                if kind is Kind.FINDING
                else "not a carrier"
            )
        gene_call: Fields = {
            "run_id": ctx.run_id,
            "gene": locus.locus,
            "diplotype": f"{locus.locus} {copies}" if copies is not None else None,
            "phenotype": phenotype,
            "activity_score": None,
            "callers": json.dumps(callers),
            "agreed": self._agreed(locus),
            "level": str(Level.MEDIUM if assessed else Level.LOW),
            "warnings": []
            if assessed or reason is None
            else [f"{GENE_UNRESOLVED}: {reason}"],
        }
        return Row(None, None, statement, self._measured(locus)), gene_call

    def collect(self, ctx: RunContext) -> Collected:
        """Each locus's statement and gene call, then the not-covered genes'."""
        rows: list[Row] = []
        gene_calls: list[Fields] = []
        for locus in LOCI.values():
            row, gene_call = self._rows(ctx, locus)
            rows.append(row)
            gene_calls.append(gene_call)
        for gene, (topic, reason) in NOT_COVERED.items():
            claim = {
                "schema": CLAIM_SCHEMA,
                "type": "paralog_not_covered",
                "gene": gene,
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"{PARALOGS_TOOL}:{gene}",
                "kind": str(Kind.UNASSESSED),
                "vid": None,
                "gene": gene,
                "claim": json.dumps(claim),
                "validity": str(Validity.LIMITED),
                "source_level": SOURCE_LEVEL,
                "reason": reason,
                "effect": None,
                "actionable": False,
                "topic": str(topic),
                "source": PARALOGS_TOOL,
                "source_id": gene,
                "source_version": SOURCE_VERSION,
                "points": None,
            }
            rows.append(Row(None, None, statement))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={"gene_call": gene_calls},
            input_rows=len(rows),
            dropped={},
        )


def gba1_contradicted(
    statement: Fields, gauchian: GauchianCall | None, exons: tuple[str, int, int]
) -> bool:
    """Whether Gauchian contradicts a ClinVar or predicted statement in GBA1's exons 9-11
    (`exons`, 0-based, `GBA1_ex9-11`): Gauchian ran and is wild type, or names no variant
    at the statement's position (`GAUCHIAN_SITES`). A name Gauchian gives that the table
    cannot place contradicts nothing: the statement keeps its grade."""
    vid = statement.get("vid")
    if statement.get("source") not in _GRADED or gauchian is None:
        return False
    if not isinstance(vid, str):
        return False
    _, chrom, pos, *_ = vid.split(":")
    contig, start, end = exons
    if chrom != contig or not start < int(pos) <= end:
        return False
    if wild_type(gauchian):
        return True
    named = [GAUCHIAN_SITES.get(n) for n in (*gauchian.variants, *gauchian.other)]
    return None not in named and int(pos) not in named


@dataclass(frozen=True)
class LpaSource:
    """KILDA's KIV-2 copy number (None: it gave none), why the calibration withholds LPA
    (`withheld`) and why KILDA did not run (`not_run`)."""

    kiv2: float | None
    withheld: str | None
    not_run: str | None

    @property
    def tool(self) -> str:
        return LPA_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """The one `lpa:LPA` statement."""
        reason = self.not_run or self.withheld
        if reason is None and self.kiv2 is None:
            reason = f"{KILDA_NAME} gave no KIV-2 copy number"
        found = reason is None and self.kiv2 is not None
        kiv2 = round(self.kiv2, 1) if found and self.kiv2 is not None else None
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "lpa_kiv2",
            "gene": "LPA",
            "kiv2": kiv2,
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"{LPA_TOOL}:LPA",
            "kind": str(Kind.FINDING if found else Kind.UNASSESSED),
            "vid": None,
            "gene": "LPA",
            "claim": json.dumps(claim),
            "validity": str(Validity.LIMITED),
            "source_level": KILDA_NAME,
            "reason": LPA_REASON if found else reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.WELLNESS),
            "source": LPA_TOOL,
            "source_id": "LPA",
            "source_version": KILDA_VERSION,
            "points": None,
        }
        return Collected(
            tool=self.tool,
            rows=[Row(None, None, statement, present({"kiv2": self.kiv2}))],
            extra={},
            input_rows=1,
            dropped={},
        )
