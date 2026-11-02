"""Deletions as a source (`sources.deletions`): screened deletions built by hand become
statements, variants and calls, with the screen's own statement and a balanced account."""

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from pathlib import Path

import pytest
from fixtures.run_context import run_context

from genokeeper.deletions.screen import (
    NO_EXON,
    SHORT,
    CallClass,
    Deletion,
    PopulationSv,
    Shown,
)
from genokeeper.evidence.predicted import (
    BOTH_COPIES,
    NOT_HAPLOINSUFFICIENT,
    ONE_RECESSIVE_COPY,
)
from genokeeper.sources import Collected, RunContext
from genokeeper.sources.deletions import (
    CONFIRM_BY,
    CONSEQUENCE_WARNING,
    DELETION_WORDING,
    DELETIONS_TOOL,
    FLANK_ONLY,
    LIKELY_REASON,
    POSSIBLE_REASON,
    SOURCE_LEVEL,
    ZYGOSITY_WARNING,
    DeletionSource,
    Wording,
)
from genokeeper.store.db import Kind, Topic

LONG_BOTH = CallClass("long", "both+depth", False, False)
SHORT_DYSGU = CallClass("short", "dysgu", True, False)
CALLERS = {"dysgu": "1234 QUAL 45.0", "delly": "DEL00012 QUAL 300"}
# ClinGen rates both haploinsufficient (dosage score 3); PCSK9 and MUTYH it does not.
HAPLOINSUFFICIENT = frozenset({"BRCA1", "PALB2"})
# A deletion table's recall per band (`DeletionTable.recall`).
RECALL = {
    "long": {"pcr": (85, 197), "pcrfree": (85, 197)},
    "short": {"pcr": (1191, 7260), "pcrfree": (1325, 7260)},
}
type Rate = tuple[float, float, float] | None


def _shown(
    gene: str = "BRCA1",
    *,
    chrom: str = "chr17",
    pos: int = 43044000,
    end: int = 43046000,
    gt: str = "0/1",
    call_class: CallClass = LONG_BOTH,
    exons: Mapping[str, tuple[int, ...]] | None = None,
    population: PopulationSv | None = None,
) -> Shown:
    d = Deletion(chrom, pos, end, "A", gt, CALLERS, 0.4)
    return Shown(d, call_class, exons or {gene: (2, 3)}, population)


def _rate(value: float | None) -> Callable[[CallClass], Rate]:
    return lambda c: None if value is None else (value, value / 2, value * 2)


RATE_005 = _rate(0.05)


def _collect(
    tmp_path: Path,
    shown: Sequence[Shown],
    *,
    dropped: Mapping[str, int] | None = None,
    not_run: str | None = None,
    p_wrong: Callable[[CallClass], Rate] = RATE_005,
    haploinsufficient: frozenset[str] = HAPLOINSUFFICIENT,
    recall: Mapping[str, Mapping[str, tuple[int, int]]] = RECALL,
    inheritance: Mapping[str, str] | None = None,
) -> Collected:
    dropped = dropped or {}
    called = len(shown) + sum(dropped.values())
    source = DeletionSource(
        shown, called, dropped, not_run, p_wrong, haploinsufficient, recall
    )
    ctx = replace(run_context(tmp_path), inheritance=dict(inheritance or {}))
    return source.collect(ctx)


def _gene_rows(collected: Collected) -> list[dict[str, object]]:
    return [
        r.statement
        for r in collected.rows
        if r.statement["statement_id"] != "deletions:screen"
    ]


def test_finding_wording_and_not_actionable(tmp_path: Path) -> None:
    collected = _collect(tmp_path, [_shown()])
    (statement,) = _gene_rows(collected)
    assert statement["statement_id"] == "deletions:BRCA1:chr17:43044000-43046000"
    assert statement["vid"] == "GRCh38:chr17:43044000:A:<DEL>"
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["validity"] == "limited"
    assert statement["topic"] == "risks"
    assert statement["actionable"] is False
    assert statement["source_level"] == SOURCE_LEVEL
    assert statement["source_version"] == "DELLY 2.7.0; Dysgu 1.9.0; Truvari 5.5.0"
    reason = str(statement["reason"])
    assert reason == (
        f"{LIKELY_REASON}; predicted, not classified; discuss with a geneticist"
    )
    assert LIKELY_REASON == (
        "likely present, needs confirmation by MLPA, a microarray, qPCR or long reads"
    )
    # "Likely" below `likely_below` (0.1): a chance of being true above 0.9.
    for rate in (0.1, 0.3, 0.6):
        (possible,) = _gene_rows(_collect(tmp_path, [_shown()], p_wrong=_rate(rate)))
        assert str(possible["reason"]).startswith(POSSIBLE_REASON)


def test_carrier_gene_topic_carrier(tmp_path: Path) -> None:
    shown = _shown("CFTR", chrom="chr7", pos=117530000, end=117531000)
    (statement,) = _gene_rows(_collect(tmp_path, [shown]))
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == "carrier"


def test_het_in_recessive_gene_outside_carrier_lists_non_finding(
    tmp_path: Path,
) -> None:
    """GENEZ is recessive (CGD), on no list and not haploinsufficient: one copy is a
    non-finding and both copies go under research (1e's gene gate, `lof_gene_gate`).
    MUTYH, recessive in ACMG SF and a carrier Tier 4 gene, is carrier status for one."""
    ar = {"GENEZ": "AR"}
    het = _shown("GENEZ", chrom="chr1", pos=45330000, end=45332000)
    (statement,) = _gene_rows(_collect(tmp_path, [het], inheritance=ar))
    assert statement["kind"] == str(Kind.NON_FINDING)
    assert ONE_RECESSIVE_COPY in str(statement["reason"])
    mutyh = _shown("MUTYH", chrom="chr1", pos=45330000, end=45332000)
    (statement,) = _gene_rows(_collect(tmp_path, [mutyh]))
    assert (statement["kind"], statement["topic"]) == (str(Kind.FINDING), "carrier")
    hom = _shown("GENEZ", chrom="chr1", pos=45330000, end=45332000, gt="1/1")
    collected = _collect(tmp_path, [hom], inheritance=ar)
    (statement,) = _gene_rows(collected)
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["topic"] == "research"
    assert statement["validity"] == "limited"
    assert BOTH_COPIES in str(statement["reason"])
    call = collected.rows[1].call
    assert call is not None
    assert (call["gt"], call["dosage"]) == ("1/1", "2")


def test_not_haploinsufficient_non_finding(tmp_path: Path) -> None:
    """A deletion in an ACMG SF gene neither haploinsufficient nor Tier 3 (PCSK9, whose
    risk is a gain of function) is a non-finding that says why, at either dosage."""
    for gt in ("0/1", "1/1"):
        shown = _shown("PCSK9", chrom="chr1", pos=55039000, end=55041000, gt=gt)
        (statement,) = _gene_rows(_collect(tmp_path, [shown]))
        assert statement["kind"] == str(Kind.NON_FINDING)
        assert statement["topic"] is None
        assert statement["reason"] == (
            f"{LIKELY_REASON}; {NOT_HAPLOINSUFFICIENT}; predicted, not classified;"
            " discuss with a geneticist"
        )


def test_gene_gate_reads_haploinsufficient(tmp_path: Path) -> None:
    """BRCA1 is placed under risks only as a haploinsufficient gene; CFTR, in Tier 3, is
    placed whatever ClinGen rates."""
    (brca1,) = _gene_rows(_collect(tmp_path, [_shown()], haploinsufficient=frozenset()))
    assert (brca1["kind"], brca1["topic"]) == (str(Kind.NON_FINDING), None)
    assert NOT_HAPLOINSUFFICIENT in str(brca1["reason"])
    cftr = _shown("CFTR", chrom="chr7", pos=117530000, end=117531000)
    (statement,) = _gene_rows(_collect(tmp_path, [cftr], haploinsufficient=frozenset()))
    assert (statement["kind"], statement["topic"]) == (str(Kind.FINDING), "carrier")


def test_two_genes_one_call_row(tmp_path: Path) -> None:
    shown = _shown(exons={"BRCA1": (1,), "PALB2": (4,)})
    collected = _collect(tmp_path, [shown])
    genes = [r for r in collected.rows if r.statement["vid"] is not None]
    assert [r.statement["gene"] for r in genes] == ["BRCA1", "PALB2"]
    with_rows = [r for r in genes if r.variant is not None and r.call is not None]
    assert len(with_rows) == 1
    assert [r for r in genes if r.variant is None and r.call is None] == genes[1:]


def test_claim_lists_exons(tmp_path: Path) -> None:
    gnomad = PopulationSv(
        "gnomAD-SV_v3_DEL_chr17_9", "chr17", 43043900, 43046100, 0.001
    )
    shown = _shown(pos=43040000, end=43130000, exons={"BRCA1": tuple(range(1, 24))})
    shown = Shown(shown.deletion, shown.call_class, shown.exons, gnomad)
    (statement,) = _gene_rows(_collect(tmp_path, [shown]))
    claim = json.loads(str(statement["claim"]))
    assert claim["type"] == "deletion"
    assert claim["exons"] == list(range(1, 24))
    assert (claim["gene"], claim["chrom"]) == ("BRCA1", "chr17")
    assert (claim["start"], claim["end"], claim["length"]) == (
        43040001,
        43130000,
        90000,
    )
    assert claim["callers"] == CALLERS
    assert claim["depth_ratio"] == 0.4
    assert claim["class"] == {
        "band": "long",
        "evidence": "both+depth",
        "repeat": False,
        "masked": False,
    }
    assert claim["gnomad_sv"] == {"id": "gnomAD-SV_v3_DEL_chr17_9", "af": 0.001}
    assert claim["false_call_rate"] == [0.05, 0.025, 0.1]


@pytest.mark.parametrize(
    ("rate", "level"),
    [(0.05, "low"), (0.005, "medium"), (0.0005, "high")],
)
def test_level_from_p_wrong(tmp_path: Path, rate: float | None, level: str) -> None:
    collected = _collect(tmp_path, [_shown()], p_wrong=_rate(rate))
    call = collected.rows[1].call
    assert call is not None
    assert call["level"] == level
    assert call["row_type"] == "PASS"
    assert call["p_wrong"] == rate


def test_unmeasured_band_withheld(tmp_path: Path) -> None:
    """A class whose band the table did not measure (no call there: `p_wrong` None) is
    withheld: unassessed, naming the band, with no vid, variant or call."""
    collected = _collect(tmp_path, [_shown()], p_wrong=_rate(None))
    (statement,) = _gene_rows(collected)
    assert (statement["kind"], statement["vid"]) == (str(Kind.UNASSESSED), None)
    assert statement["reason"] == "not calibrated: no long calls measured"
    assert [(r.variant, r.call) for r in collected.rows] == [(None, None)] * 2


def test_class_sets_rate_and_warnings(tmp_path: Path) -> None:
    """Review Focus 6: a Dysgu-only 120 bp deletion with an end in a tandem repeat, over a
    CFTR exon, beside a deletion of another class."""
    rates = {LONG_BOTH: (0.1, 0.05, 0.2), SHORT_DYSGU: (0.7, 0.6, 0.8)}
    d = Deletion("chr7", 117530000, 117530120, "A", "0/1", {"dysgu": "7 QUAL 9"}, None)
    cftr = Shown(d, SHORT_DYSGU, {"CFTR": (4,)}, None)
    collected = _collect(tmp_path, [_shown(), cftr], p_wrong=rates.get)
    calls = [r.call for r in collected.rows if r.call is not None]
    assert [(c["p_wrong"], c["p_wrong_lo"], c["p_wrong_hi"]) for c in calls] == [
        rates[LONG_BOTH],
        rates[SHORT_DYSGU],
    ]
    assert calls[0]["warnings"] == [CONSEQUENCE_WARNING]
    assert calls[1]["warnings"] == [
        "found by one caller (Dysgu)",
        "read depth not measured",
        "an end in a tandem repeat: may be a repeat-length change",
        CONSEQUENCE_WARNING,
    ]
    assert str(_gene_rows(collected)[1]["reason"]).startswith(POSSIBLE_REASON)


def test_cautions_join_the_class_warnings(tmp_path: Path) -> None:
    """What the table says of the call itself (`cautions`: a class it never saw, a target
    outside its benchmark) follows the class's warnings."""
    source = DeletionSource(
        [_shown()], 1, {}, None, RATE_005, HAPLOINSUFFICIENT, RECALL,
        cautions=lambda shown: [f"caution {shown.call_class.band}"],
    )  # fmt: skip
    collected = source.collect(run_context(tmp_path))
    (call,) = [r.call for r in collected.rows if r.call is not None]
    assert call["warnings"] == ["caution long", CONSEQUENCE_WARNING]


def test_flank_only_gene_reason_says_so(tmp_path: Path) -> None:
    """A deletion that touches one gene's exon and only another's flank says so in the
    second gene's reason alone; its placement stays the gate's."""
    shown = replace(
        _shown(exons={"BRCA1": (2,), "PALB2": (1,)}), flank=frozenset({"PALB2"})
    )
    rows = {
        r["gene"]: str(r["reason"]) for r in _gene_rows(_collect(tmp_path, [shown]))
    }
    assert FLANK_ONLY in rows["PALB2"]
    assert FLANK_ONLY not in rows["BRCA1"]


def test_claim_names_plus_clinical_exons_apart(tmp_path: Path) -> None:
    """An exon only the Plus Clinical transcript has goes in `plus_clinical_exons`, never
    in `exons`, which are Select's numbers."""
    shown = replace(
        _shown("DYNC2H1", exons={"DYNC2H1": (64,)}),
        plus_clinical={"DYNC2H1": (64,)},
    )
    (row,) = _gene_rows(_collect(tmp_path, [shown]))
    claim = json.loads(str(row["claim"]))
    assert (claim["exons"], claim["plus_clinical_exons"]) == ([64], [64])
    (plain,) = _gene_rows(_collect(tmp_path, [_shown()]))
    assert json.loads(str(plain["claim"]))["plus_clinical_exons"] == []


def test_zygosity_and_consequence_warnings(tmp_path: Path) -> None:
    """A deletion called on both copies warns that zygosity is not calibrated; every one
    warns that its consequence is not assessed, last."""
    for gt, expected in (
        ("1/1", [ZYGOSITY_WARNING, CONSEQUENCE_WARNING]),
        ("0/1", [CONSEQUENCE_WARNING]),
    ):
        collected = _collect(tmp_path, [_shown(gt=gt)])
        (call,) = [r.call for r in collected.rows if r.call is not None]
        assert call["warnings"] == expected
    assert CONSEQUENCE_WARNING == (
        "its consequence (whole or partial exon, in frame, last exon, decay escape) is"
        " not assessed"
    )


def test_not_run_unassessed(tmp_path: Path) -> None:
    collected = _collect(tmp_path, [], not_run="no deletion table")
    (row,) = collected.rows
    s = row.statement
    assert s["statement_id"] == "deletions:screen"
    assert (s["kind"], s["validity"], s["topic"]) == ("unassessed", "strong", "carrier")
    assert s["gene"] is None
    assert s["reason"] == "no deletion table"
    assert s["source"] == DELETIONS_TOOL
    assert s["source_level"] == SOURCE_LEVEL
    assert s["source_version"] == "DELLY 2.7.0; Dysgu 1.9.0; Truvari 5.5.0"
    assert (collected.input_rows, collected.dropped) == (1, {})


def test_screened_none_shown_non_finding(tmp_path: Path) -> None:
    collected = _collect(tmp_path, [], dropped={SHORT: 2, NO_EXON: 1})
    (row,) = collected.rows
    assert row.statement["kind"] == str(Kind.NON_FINDING)
    assert row.statement["validity"] == "limited"
    assert row.statement["reason"] == (
        "screened: 3 deletions called, 0 shown; a miss is possible: on HG002's libraries"
        " the screen found 85 of 197 (pcr), 85 of 197 (pcrfree) rare long and 1191 of"
        " 7260 (pcr), 1325 of 7260 (pcrfree) rare short truth deletions"
    )
    (unrecorded,) = _collect(tmp_path, [], recall={}).rows
    assert str(unrecorded.statement["reason"]).endswith(
        "a miss is possible: the screen's recall on HG002 is not recorded"
    )
    assert collected.input_rows == 4
    assert collected.dropped == {SHORT: 2, NO_EXON: 1}


def test_accounting_balanced(tmp_path: Path) -> None:
    two = _shown(exons={"BRCA1": (1,), "PALB2": (4,)})
    one = _shown("CFTR", chrom="chr7", pos=117530000, end=117531000)
    collected = _collect(tmp_path, [two, one], dropped={SHORT: 1, NO_EXON: 2})
    assert len(collected.rows) == 4
    assert collected.input_rows == 7
    assert collected.rows[0].statement["reason"] == (
        "screened: 5 deletions called, 2 shown; a miss is possible: on HG002's libraries"
        " the screen found 85 of 197 (pcr), 85 of 197 (pcrfree) rare long and 1191 of"
        " 7260 (pcr), 1325 of 7260 (pcrfree) rare short truth deletions"
    )


# --- the wording and placement hooks (slice 3d) ---------------------------------------------


def _source(
    shown: Sequence[Shown], *, not_run: str | None = None, **options: object
) -> DeletionSource:
    return DeletionSource(
        shown, 3, {SHORT: 2}, not_run, RATE_005, HAPLOINSUFFICIENT,
        {"long": {"pcr": (85, 197)}}, **options,  # pyright: ignore[reportArgumentType]
    )  # fmt: skip


def test_deletion_wording_unchanged(tmp_path: Path) -> None:
    """3b's statements byte for byte under the defaults (taken from 94da76e's code)."""
    source = _source([_shown(exons={"BRCA1": (2, 3)})])
    assert source.wording is DELETION_WORDING
    collected = source.collect(run_context(tmp_path))
    assert (collected.tool, collected.input_rows) == ("deletions", 4)
    screen_row, gene_row = (r.statement for r in collected.rows)
    common = {
        "source": "deletions",
        "source_level": "DELLY 2.7.0 and Dysgu 1.9.0, matched by Truvari 5.5.0",
        "source_version": "DELLY 2.7.0; Dysgu 1.9.0; Truvari 5.5.0",
    }
    assert {k: screen_row[k] for k in (*common, "statement_id", "reason", "claim")} == {
        **common,
        "statement_id": "deletions:screen",
        "reason": "screened: 3 deletions called, 1 shown; a miss is possible: on HG002's"
        " libraries the screen found 85 of 197 (pcr) rare long truth deletions",
        "claim": '{"schema": "genomelab.claim/1", "type": "deletion_screen", "called": 3,'
        ' "shown": 1}',
    }
    assert {
        k: gene_row[k]
        for k in (*common, "statement_id", "reason", "claim", "vid", "source_id")
    } == {
        **common,
        "statement_id": "deletions:BRCA1:chr17:43044000-43046000",
        "reason": "likely present, needs confirmation by MLPA, a microarray, qPCR or long"
        " reads; predicted, not classified; discuss with a geneticist",
        "claim": '{"schema": "genomelab.claim/1", "type": "deletion", "gene": "BRCA1",'
        ' "chrom": "chr17", "start": 43044001, "end": 43046000, "length": 2000, "exons":'
        ' [2, 3], "plus_clinical_exons": [], "callers": {"dysgu": "1234 QUAL 45.0",'
        ' "delly": "DEL00012 QUAL 300"},'
        ' "depth_ratio": 0.4, "class": {"band": "long", "evidence": "both+depth",'
        ' "repeat": false, "masked": false}, "gnomad_sv": null, "false_call_rate":'
        " [0.05, 0.025, 0.1]}",
        "vid": "GRCh38:chr17:43044000:A:<DEL>",
        "source_id": "chr17:43044000-43046000",
    }
    call = collected.rows[1].call
    assert call is not None
    assert call["warnings"] == [CONSEQUENCE_WARNING]


CNV_LOSS = Wording(
    tool="cnv",
    noun="copy-number loss",
    claim_type="copy_number_loss",
    confirm_by=CONFIRM_BY,
    alt="<DEL>",
    source_level="DELLY 2.7.0 cnv and CNVpytor 1.3.2",
    source_version="DELLY 2.7.0; CNVpytor 1.3.2",
    callers=("dellycnv", "cnvpytor"),
)


def test_wording_sets_source_and_type(tmp_path: Path) -> None:
    """No row of another wording says `deletions` (C24); `screen_row=False` drops the
    `*:screen` row; a `place` given decides the kind, topic and placement line."""
    depth = CallClass("loss_1-5kb", "dellycnv", False, False)
    d = Deletion(
        "chr17", 43044000, 43046000, "N", "0/1", {"dellycnv": "CNV1 CN 1"}, None
    )
    shown = Shown(d, depth, {"BRCA1": (2, 3)}, None, kinds={"BRCA1": "loss"})
    collected = _source([shown], wording=CNV_LOSS).collect(run_context(tmp_path))
    assert collected.tool == "cnv"
    statements = [r.statement for r in collected.rows]
    assert [s["statement_id"] for s in statements] == [
        "cnv:screen",
        "cnv:BRCA1:chr17:43044000-43046000",
    ]
    assert {s["source"] for s in statements} == {"cnv"}
    assert {s["source_level"] for s in statements} == {CNV_LOSS.source_level}
    assert {s["source_version"] for s in statements} == {CNV_LOSS.source_version}
    claims = [json.loads(str(s["claim"])) for s in statements]
    assert [c["type"] for c in claims] == [
        "copy_number_loss_screen",
        "copy_number_loss",
    ]
    assert str(statements[0]["reason"]).startswith(
        "screened: 3 copy-number losses called, 1 shown;"
    )
    assert str(statements[0]["reason"]).endswith("rare long truth copy-number losses")
    assert statements[1]["vid"] == "GRCh38:chr17:43044000:N:<DEL>"
    assert str(statements[1]["reason"]).startswith(f"likely present, {CONFIRM_BY}")
    (possible,) = _gene_rows(
        DeletionSource(
            [shown],
            1,
            {},
            None,
            _rate(0.5),
            HAPLOINSUFFICIENT,
            {},
            wording=CNV_LOSS,
            screen_row=False,
        ).collect(run_context(tmp_path))
    )
    assert str(possible["reason"]).startswith(
        f"possible copy-number loss, {CONFIRM_BY}"
    )
    call = collected.rows[1].call
    assert call is not None
    assert call["warnings"] == ["found by one caller (DELLY cnv)", CONSEQUENCE_WARNING]

    seen: list[str] = []

    def place(
        ctx: RunContext, gene: str, s: Shown
    ) -> tuple[Kind, Topic | None, str | None]:
        seen.append(s.kinds[gene])
        return Kind.FINDING, Topic.RESEARCH, "possible extra copy, effect unknown"

    alone = _source([shown], wording=CNV_LOSS, screen_row=False, place=place)
    collected = alone.collect(run_context(tmp_path))
    (row,) = collected.rows
    assert row.statement["statement_id"] == "cnv:BRCA1:chr17:43044000-43046000"
    assert (row.statement["kind"], row.statement["topic"]) == ("finding", "research")
    assert "possible extra copy, effect unknown" in str(row.statement["reason"])
    assert seen == ["loss"]
    assert collected.input_rows == 1 + 2
    (unrun,) = [
        _source([], wording=CNV_LOSS, screen_row=False, not_run="x").collect(
            run_context(tmp_path)
        )
    ]
    assert (unrun.rows, unrun.input_rows) == ([], 0)


def test_depth_support_reason_only(tmp_path: Path) -> None:
    shown = _shown()
    d = shown.deletion
    plain = _source([shown]).collect(run_context(tmp_path))
    joined = _source(
        [shown], depth_support={(d.chrom, d.pos, d.end): "DELLY cnv"}
    ).collect(run_context(tmp_path))
    before, after = (_gene_rows(c)[0] for c in (plain, joined))
    assert after["reason"] == str(before["reason"]).replace(
        f"{LIKELY_REASON}; ",
        f"{LIKELY_REASON}; read depth also supports it (DELLY cnv); ",
    )
    assert "read depth also supports it (DELLY cnv)" in str(after["reason"])
    assert after["claim"] == before["claim"]
    calls = [c.rows[1].call for c in (plain, joined)]
    assert calls[0] == calls[1]


def test_placement_override(tmp_path: Path) -> None:
    """A4: a gene's override replaces the gate: a het DPYD deletion is a drugs finding."""
    dpyd = _shown("DPYD", chrom="chr1", pos=97500000, end=97502000)
    (gated,) = _gene_rows(_source([dpyd]).collect(run_context(tmp_path)))
    assert gated["kind"] == str(Kind.NON_FINDING)
    assert NOT_HAPLOINSUFFICIENT in str(gated["reason"])
    overridden = _source(
        [dpyd], placement_overrides={"DPYD": ("drug_response", "a note")}
    ).collect(run_context(tmp_path))
    (statement,) = _gene_rows(overridden)
    assert (statement["kind"], statement["topic"]) == ("finding", "drugs")
    gates = (NOT_HAPLOINSUFFICIENT, ONE_RECESSIVE_COPY, BOTH_COPIES)
    assert not any(g in str(statement["reason"]) for g in gates)
    assert statement["reason"] == (
        f"{LIKELY_REASON}; a note; predicted, not classified; discuss with a geneticist"
    )
