"""Depth copy number as a source (`sources.cnv`, slice 3d): shown losses and gains through
3b's statements in their own words and placements, large events and named loci as rows of
their own, one `cnv:screen` statement, and a balanced account."""

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import cast

from fixtures.run_context import run_context

from genokeeper.cnv.screen import LargeEvent, NamedCnv, load_named
from genokeeper.deletions.screen import CallClass, Deletion, Shown
from genokeeper.evidence.predicted import ONE_RECESSIVE_COPY
from genokeeper.genelists import DPYD_NOTE, screen_overrides
from genokeeper.sources import Collected
from genokeeper.sources.cnv import (
    CNV_TOOL,
    GAIN_DEPTH_ONLY,
    GAIN_WORDING,
    LARGE_NOT_CALIBRATED,
    LARGE_REASON,
    LOSS_WORDING,
    NOT_COVERED_DEPTH,
    TANDEM_NOTE,
    WHOLE_GAIN,
    CnvSource,
)
from genokeeper.store.db import Kind, Topic

type Rate = tuple[float, float, float] | None

BOTH = {"dellycnv": "CNV00001 CN 1", "cnvpytor": "chr7:1-2 e1 1e-9"}
LOSS = CallClass("loss_5kb+", "both", False, False)
GAIN = CallClass("gain_5kb+", "both", False, False, True)
LARGE = CallClass("loss_large", "both", False, False)


def _rate(c: CallClass) -> Rate:
    return None if c.band.endswith("large") else (0.05, 0.02, 0.1)


def _shown(
    gene: str,
    kind: str,
    *,
    chrom: str = "chr7",
    pos: int = 117_530_000,
    end: int = 117_540_000,
    gt: str = "0/1",
) -> Shown:
    svtype = "DEL" if kind == "loss" else "DUP"
    d = Deletion(chrom, pos, end, "N", gt, BOTH, None, svtype=svtype)
    c = LOSS if kind == "loss" else GAIN
    return Shown(d, c, {gene: (5, 6)}, None, kinds={gene: kind})


def _source(
    shown: Sequence[Shown] = (),
    *,
    large: Sequence[LargeEvent] = (),
    named: Sequence[tuple[NamedCnv, Shown]] = (),
    dropped: Mapping[str, int] | None = None,
    not_run: str | None = None,
    p_wrong: Callable[[CallClass], Rate] = _rate,
    triplosensitive: frozenset[str] = frozenset(),
    overrides: Mapping[str, tuple[str, str]] | None = None,
    large_withheld: Mapping[str, str] | None = None,
    cautions: Callable[[Shown], list[str]] = lambda _: [],
) -> CnvSource:
    dropped = dict(dropped or {})
    return CnvSource(
        shown=list(shown),
        large=list(large),
        named=list(named),
        called=len(shown) + len(large) + len(named) + sum(dropped.values()),
        dropped=dropped,
        not_run=not_run,
        p_wrong=p_wrong,
        haploinsufficient=frozenset({"BRCA1"}),
        triplosensitive=triplosensitive,
        placement_overrides=overrides or {},
        segdup_exons={"GENEA": (4,)},
        large_withheld=large_withheld or {},
        deletion_cautions=cautions,
    )


def _collect(
    source: CnvSource, tmp_path: Path, inheritance: Mapping[str, str] | None = None
) -> Collected:
    ctx = replace(run_context(tmp_path), inheritance=dict(inheritance or {}))
    return source.collect(ctx)


def _statements(collected: Collected) -> dict[str, dict[str, object]]:
    return {str(r.statement["statement_id"]): r.statement for r in collected.rows}


def _gene_row(collected: Collected) -> dict[str, object]:
    (row,) = [
        s for sid, s in _statements(collected).items() if sid != f"{CNV_TOOL}:screen"
    ]
    return row


def _parts(statement: Mapping[str, object]) -> list[str]:
    return str(statement["reason"]).split("; ")


def test_wordings() -> None:
    assert (LOSS_WORDING.tool, GAIN_WORDING.tool) == ("cnv", "cnv")
    assert LOSS_WORDING.nouns == "copy-number losses"
    assert GAIN_WORDING.nouns == "copy-number gains"
    assert (LOSS_WORDING.alt, GAIN_WORDING.alt) == ("<CNV>", "<DUP>")
    assert LOSS_WORDING.claim_type != GAIN_WORDING.claim_type
    assert LOSS_WORDING.callers == ("dellycnv", "cnvpytor")


def test_one_cnv_screen_row(tmp_path: Path) -> None:
    """A loss and a gain shown: exactly one `cnv:screen` statement, its claim naming what
    the depth screen does not cover and the target exons in segmental duplications; every
    row's source is `cnv`."""
    collected = _collect(
        _source([_shown("CFTR", "loss"), _shown("GENEA", "partial", pos=1_000)]),
        tmp_path,
    )
    statements = _statements(collected)
    screens = [sid for sid in statements if sid.endswith(":screen")]
    assert screens == ["cnv:screen"]
    screen = statements["cnv:screen"]
    assert screen["kind"] == str(Kind.NON_FINDING)
    assert str(screen["reason"]).startswith(
        "screened: 2 copy-number changes called, 2 shown"
    )
    claim = json.loads(str(screen["claim"]))
    assert claim["not_covered"] == list(NOT_COVERED_DEPTH)
    assert len(NOT_COVERED_DEPTH) == 4
    assert claim["segdup_exons"] == {"GENEA": [4]}
    assert {r.statement["source"] for r in collected.rows} == {CNV_TOOL}
    assert {s["vid"] for s in statements.values() if s["vid"]} == {
        "GRCh38:chr7:117530001:N:<CNV>",
        "GRCh38:chr7:1001:N:<DUP>",
    }
    assert not any(s["kind"] == str(Kind.WATCH) for s in statements.values())


def test_depth_vid_names_its_record(tmp_path: Path) -> None:
    """A depth event's vid names its record's POS (its first base, `pos + 1`) and REF, as a
    3b deletion's names its own, and the haploid look-up takes that position; its ALT
    (`<CNV>`, not `<DEL>`) keeps it apart from a 3b deletion at that POS it was not joined
    to (C6)."""
    seen: list[int] = []
    ctx = replace(run_context(tmp_path), haploid=lambda _, pos: bool(seen.append(pos)))
    collected = _source([_shown("CFTR", "loss")]).collect(ctx)
    assert _gene_row(collected)["vid"] == "GRCh38:chr7:117530001:N:<CNV>"
    assert set(seen) == {117_530_001}


def test_not_run_is_one_unassessed_screen_row(tmp_path: Path) -> None:
    collected = _collect(
        _source([_shown("CFTR", "loss")], not_run="not calibrated"), tmp_path
    )
    (row,) = collected.rows
    assert (row.statement["kind"], row.statement["reason"]) == (
        str(Kind.UNASSESSED),
        "not calibrated",
    )


def test_whole_gene_gain_with_inner_junction(tmp_path: Path) -> None:
    """Review Focus 2: a whole-gene gain (whatever junction lies inside) is a finding under
    risks in a triplosensitive gene, else a non-finding "no known effect of an extra
    copy"."""
    whole = _shown("MECP2", "whole", chrom="chrX", pos=154_000_000, end=154_200_000)
    finding = _gene_row(
        _collect(_source([whole], triplosensitive=frozenset({"MECP2"})), tmp_path)
    )
    assert (finding["kind"], finding["topic"]) == (str(Kind.FINDING), "risks")
    assert finding["validity"] == "limited" and finding["actionable"] is False
    other = _gene_row(_collect(_source([whole]), tmp_path))
    assert (other["kind"], other["topic"]) == (str(Kind.NON_FINDING), None)
    assert WHOLE_GAIN in _parts(other)
    assert WHOLE_GAIN == "no known effect of an extra copy"


def test_tandem_gain_through_gate(tmp_path: Path) -> None:
    """Both ends inside the gene with a junction: 1e's gene gate, as a loss of function. A
    het in a Tier 3 gene is a carrier finding; a het in a recessive gene outside the carrier
    lists a non-finding naming the gate."""
    carrier = _gene_row(_collect(_source([_shown("CFTR", "tandem")]), tmp_path))
    assert (carrier["kind"], carrier["topic"]) == (str(Kind.FINDING), "carrier")
    assert TANDEM_NOTE in _parts(carrier)
    gated = _gene_row(
        _collect(_source([_shown("GENEZ", "tandem")]), tmp_path, {"GENEZ": "AR"})
    )
    assert (gated["kind"], gated["topic"]) == (str(Kind.NON_FINDING), None)
    assert ONE_RECESSIVE_COPY in _parts(gated)


def test_depth_only_gain_research(tmp_path: Path) -> None:
    """Both ends inside and no junction, or one end inside: a finding under research at
    validity limited, "possible extra copy, effect unknown"; never a watch."""
    for kind in ("intragenic", "partial"):
        row = _gene_row(_collect(_source([_shown("CFTR", kind)]), tmp_path))
        assert (row["kind"], row["topic"], row["validity"]) == (
            str(Kind.FINDING),
            "research",
            "limited",
        ), kind
        assert GAIN_DEPTH_ONLY in _parts(row)
        assert GAIN_DEPTH_ONLY == "possible extra copy, effect unknown"


def test_dpyd_depth_loss_drugs(tmp_path: Path) -> None:
    """A het loss over DPYD exons is a drugs finding with DPYD_NOTE (A4), not the gate's
    one-copy non-finding."""
    loss = _shown("DPYD", "loss", chrom="chr1", pos=97_500_000, end=97_510_000)
    row = _gene_row(_collect(_source([loss], overrides=screen_overrides()), tmp_path))
    assert (row["kind"], row["topic"]) == (str(Kind.FINDING), str(Topic.DRUGS))
    assert DPYD_NOTE in _parts(row)
    gated = _gene_row(_collect(_source([loss]), tmp_path))
    assert gated["kind"] == str(Kind.NON_FINDING)


def _large(pos: int = 18_650_000, end: int = 21_150_000) -> LargeEvent:
    d = Deletion("chr22", pos, end, "N", "0/1", BOTH, None)
    return LargeEvent(d, LARGE, ("22q11.2 recurrent (DGS)", "TBX1"))


def test_large_event_needs_microarray(tmp_path: Path) -> None:
    """A large event is one statement of its own: gene None, a finding under risks at
    validity limited, "needs a chromosomal microarray", shown not calibrated, its call's
    chance of being false not measured (level low), never actionable. A band the table
    withholds gives no statement, and the screen's claim names it."""
    collected = _collect(_source(large=[_large()]), tmp_path)
    statements = _statements(collected)
    row = statements["cnv:large:chr22:18650000-21150000"]
    assert (row["gene"], row["kind"], row["topic"], row["validity"]) == (
        None,
        str(Kind.FINDING),
        "risks",
        "limited",
    )
    assert row["actionable"] is False
    assert _parts(row)[:2] == [LARGE_REASON, LARGE_NOT_CALIBRATED]
    assert LARGE_REASON == "screening result, needs a chromosomal microarray"
    claim = json.loads(str(row["claim"]))
    assert claim["regions"] == ["22q11.2 recurrent (DGS)", "TBX1"]
    (call,) = [r.call for r in collected.rows if r.call is not None]
    assert (call["p_wrong"], call["level"]) == (None, "low")
    assert call["vid"] == row["vid"] == "GRCh38:chr22:18650001:N:<CNV>"
    withheld = _collect(
        _source(large=[_large()], large_withheld={"loss_large": "on 2 of 3 genomes"}),
        tmp_path,
    )
    assert list(_statements(withheld)) == ["cnv:screen"]
    claim = json.loads(str(_statements(withheld)["cnv:screen"]["claim"]))
    assert claim["large_withheld"] == {"loss_large": "on 2 of 3 genomes"}


def test_named_gjb6_deletion(tmp_path: Path) -> None:
    """A loss over GJB6-D13S1830 is one `cnv:named:GJB6-D13S1830` statement on GJB2: a
    carrier finding at validity limited with the named row's reason, at its class's
    chance of being false; a 3b deletion at the same place adds no second statement."""
    (gjb6,) = load_named()
    depth = Deletion("chr13", 20_210_000, 20_520_000, "N", "0/1", BOTH, None)
    deletion = Deletion(
        "chr13", 20_210_000, 20_520_000, "A", "0/1", {"delly": "D"}, 0.4
    )
    named = [
        (gjb6, Shown(depth, LOSS, {}, None)),
        (
            gjb6,
            Shown(deletion, CallClass("long", "delly+depth", False, False), {}, None),
        ),
    ]
    collected = _collect(_source(named=named), tmp_path)
    statements = _statements(collected)
    assert sorted(statements) == ["cnv:named:GJB6-D13S1830", "cnv:screen"]
    row = statements["cnv:named:GJB6-D13S1830"]
    assert (row["gene"], row["topic"], row["validity"], row["kind"]) == (
        "GJB2",
        "carrier",
        "limited",
        str(Kind.FINDING),
    )
    assert row["reason"] == gjb6.reason
    assert gjb6.reason.startswith("possible GJB6-D13S1830 deletion")
    (call,) = [r.call for r in collected.rows if r.call is not None]
    assert call["p_wrong"] == 0.05
    assert collected.input_rows == len(collected.rows) + sum(collected.dropped.values())
    homozygous = replace(depth, gt="1/1")
    both = _collect(
        _source(named=[(gjb6, Shown(homozygous, LOSS, {}, None))]), tmp_path
    )
    row = _statements(both)["cnv:named:GJB6-D13S1830"]
    assert (row["topic"], row["kind"]) == ("risks", str(Kind.FINDING))
    assert row["reason"] == (
        "possible homozygous GJB6-D13S1830 deletion: causes DFNB1 hearing loss;"
        " needs confirmation by PCR across the breakpoint"
    )
    # A band the table did not measure withholds the named locus too.
    unmeasured = _source(named=named[:1], p_wrong=lambda c: None)
    row = _statements(_collect(unmeasured, tmp_path))["cnv:named:GJB6-D13S1830"]
    assert (row["kind"], row["vid"]) == (str(Kind.UNASSESSED), None)
    assert row["reason"] == "not calibrated: no loss_5kb+ calls measured"


def test_named_3b_deletion_carries_the_deletion_tables_cautions(tmp_path: Path) -> None:
    """A named locus over a 3b deletion shows the deletion table's chance of being false,
    so its call carries that table's cautions (a class it never saw, a target outside its
    benchmark); a named depth loss, rated by the cnv table, carries none."""
    (gjb6,) = load_named()
    deletion = Deletion(
        "chr13", 20_210_000, 20_520_000, "A", "0/1", {"delly": "D"}, 0.4
    )
    long_delly = CallClass("long", "delly+depth", False, False)

    def caution(shown: Shown) -> list[str]:
        return [f"caution {shown.call_class.band}"]

    named = _collect(
        _source(
            named=[(gjb6, Shown(deletion, long_delly, {}, None))], cautions=caution
        ),
        tmp_path,
    )
    (call,) = [r.call for r in named.rows if r.call is not None]
    assert "caution long" in cast(list[str], call["warnings"])
    depth = Deletion("chr13", 20_210_000, 20_520_000, "N", "0/1", BOTH, None)
    loss = _collect(
        _source(named=[(gjb6, Shown(depth, LOSS, {}, None))], cautions=caution),
        tmp_path,
    )
    (call,) = [r.call for r in loss.rows if r.call is not None]
    assert not any(w.startswith("caution") for w in cast(list[str], call["warnings"]))


def test_cnv_accounting_balanced(tmp_path: Path) -> None:
    """Input rows are the statements plus the drops; one `call` row per vid, on its first
    statement."""
    two_genes = Shown(
        Deletion("chr1", 1_000, 50_000, "N", "0/1", BOTH, None),
        LOSS,
        {"GENEA": (1,), "GENEB": (2,)},
        None,
        kinds={"GENEA": "loss", "GENEB": "loss"},
    )
    collected = _collect(
        _source(
            [two_genes, _shown("CFTR", "partial")],
            large=[_large()],
            dropped={"under 1 kb": 3, "no listed exon": 5},
        ),
        tmp_path,
    )
    assert collected.input_rows == len(collected.rows) + 8
    assert collected.dropped == {"under 1 kb": 3, "no listed exon": 5}
    calls = [str(r.call["vid"]) for r in collected.rows if r.call is not None]
    vids = {str(r.statement["vid"]) for r in collected.rows if r.statement["vid"]}
    assert sorted(calls) == sorted(vids) and len(calls) == len(set(calls)) == 3
