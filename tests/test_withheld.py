"""A value a calibration table withholds stays visible to the agent, and out of the graded
result (`sources.withheld`; the owner's ruling, 2026-11-06)."""

import json
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import duckdb
import pytest
from fixtures.run_context import run_context

from genokeeper import dashboard, pipeline
from genokeeper.ancestry.archaic import ArchaicShare
from genokeeper.bloodgroups import SystemCall
from genokeeper.carrier import Answer, genotype
from genokeeper.cnv.screen import LargeEvent, NamedCnv
from genokeeper.config import load_grading
from genokeeper.deletions.screen import CallClass, Deletion, Shown
from genokeeper.explain import LOOSEST, grounding
from genokeeper.hhv6 import Hhv6Call
from genokeeper.mito import MitoCall, mito_called
from genokeeper.paralogs import PsCall
from genokeeper.pgx.hla import read_optitype, read_t1k
from genokeeper.pgx.hla_typing import TAG_RSIDS
from genokeeper.prs.pgsc import PopSimilarity
from genokeeper.query import query, show
from genokeeper.repeats.loci import Allele, RepeatCall, Rfc1Call
from genokeeper.roh import RohSummary, Run
from genokeeper.smn import ParascopyCall, SmnCaller
from genokeeper.sources import Collected, Source, write_collected
from genokeeper.sources.ancestry import AncestrySource
from genokeeper.sources.bloodgroups import BloodGroupSource
from genokeeper.sources.cnv import CnvSource
from genokeeper.sources.deletions import DeletionSource
from genokeeper.sources.hhv6 import Hhv6Source
from genokeeper.sources.hla import HlaSource
from genokeeper.sources.mitomap import MitomapSource
from genokeeper.sources.paralogs import LpaSource, ParalogSource
from genokeeper.sources.repeats import RepeatSource, XDosage
from genokeeper.sources.roh import RohSource
from genokeeper.sources.smn import SmnSource
from genokeeper.sources.withheld import reveal
from genokeeper.sources.yleaf import YleafSource
from genokeeper.store.db import Kind, Level, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter
from genokeeper.yhaplogroup import YHaplogroup

REASON = "not covered: false call in calibration"
HG002 = Path(__file__).parent / "fixtures" / "hg002"
RUN = Run("chr11", 48_025_953, 64_604_776, 16_578_824, 72.2)
ONE_X = XDosage(1, 0.5, True)


def _repeats(variant: str, *alleles: Allele, **kw: Any) -> RepeatSource:
    call = RepeatCall(variant, alleles, 30.0, False, 150)
    fields: dict[str, Any] = {
        "calls": {variant: call},
        "not_run": None,
        "loci": load_grading().repeats.loci,
        "withheld": {variant: REASON},
        "caps": {},
        "images": {},
        "pcr_like": False,
        "insert": 400.0,
        "male": True,
        "min_allele_reads": 5,
        "min_spanning": 3,
        "x_dosage": ONE_X,
        "rfc1_not_run": None,
    }
    return RepeatSource(**(fields | kw))


def _spanned(repeats: int, spanning: int = 20) -> Allele:
    return Allele(repeats, (repeats, repeats), spanning, 0, 0)


def _hla() -> HlaSource:
    return HlaSource(
        read_optitype(HG002 / "hg002.optitype_result.tsv.txt"),
        read_t1k(HG002 / "hg002.t1k_genotype.tsv.txt", genes=None),
        {rsid: "0" for rsid in TAG_RSIDS},
        None,
        {"HLA-B": REASON, "DQ": REASON},
    )


def _paralog() -> ParalogSource:
    return ParalogSource(
        calls={"STRC": PsCall(1, 40.0, True, 1)},
        unique_ratio=0.5,
        gauchian=None,
        aldy_cyp2a6=None,
        withheld={"STRC": REASON},
        not_run=None,
    )


def _ancestry() -> AncestrySource:
    return AncestrySource(
        PopSimilarity("EUR", False, {"EUR": 0.9, "MID": 0.1}),
        {"EUR": 0.9, "MID": 0.1},
        ArchaicShare(57_768_000, 1_704_000, 0.011324),
        not_run={},
        withheld={"neanderthal": REASON},
    )


def _mitomap(withheld: str | None = REASON) -> MitomapSource:
    return MitomapSource(
        calls=[MitoCall(pos=3243, ref="A", alt="G", level=0.30, dp=3000)],
        coverage={pos: 3000 for pos in range(1, 16570)},
        clinvar_vids=frozenset(),
        floor=0.10,
        min_homoplasmy=0.9,
        min_dp=10,
        heteroplasmy_withheld=withheld,
    )


# name -> the source, the withheld statement, what it would be, a key of its value.
CASES: dict[str, tuple[Callable[[], Source], str, Kind, str, object]] = {
    "repeat in the pathogenic range": (
        lambda: _repeats("HTT", Allele(42, (41, 43), 6, 3, 0)),
        "expansionhunter:HTT",
        Kind.FINDING,
        "range",
        "pathogenic",
    ),
    "repeat in the normal range": (
        lambda: _repeats("HTT", _spanned(17), _spanned(19)),
        "expansionhunter:HTT",
        Kind.NON_FINDING,
        "range",
        "normal",
    ),
    "an exclusion-only normal": (
        lambda: _repeats("FMR1", _spanned(30)),
        "expansionhunter:FMR1",
        Kind.NON_FINDING,
        "range",
        "normal",
    ),
    "a size its own rule leaves unassessed": (
        lambda: _repeats("FMR1", _spanned(60)),
        "expansionhunter:FMR1",
        Kind.UNASSESSED,
        "range",
        "intermediate",
    ),
    "RFC1": (
        lambda: _repeats(
            "other",
            _spanned(1),
            withheld={"RFC1": REASON},
            rfc1=Rfc1Call(
                "affected",
                "2 expanded pathogenic alleles",
                (Allele(42, (41, 43), 6, 3, 0),),
                30.0,
                False,
            ),
        ),
        "expansionhunter:RFC1",
        Kind.FINDING,
        "rfc1",
        {"call": "affected", "description": "2 expanded pathogenic alleles"},
    ),
    "SMN1": (
        lambda: SmnSource(
            SmnCaller(1, 2, 0, "high"),
            ParascopyCall(1, 2, 40.0, True),
            None,
            REASON,
            20.0,
        ),
        "smn:SMN1",
        Kind.FINDING,
        "smn1",
        1,
    ),
    "the Y haplogroup": (
        lambda: YleafSource(
            YHaplogroup("J-C134621", "C134621", 1.0, 170393), None, REASON, 0.95, True
        ),
        "yleaf:haplogroup",
        Kind.FINDING,
        "haplogroup",
        "J-C134621",
    ),
    "a paralog locus": (
        _paralog,
        "paralogs:STRC",
        Kind.FINDING,
        "copy_number",
        1,
    ),
    "Lp(a)": (
        lambda: LpaSource(20.0, REASON, None),
        "lpa:LPA",
        Kind.FINDING,
        "kiv2",
        20.0,
    ),
    "an HLA allele": (_hla, "hla:HLA-B*27", Kind.NON_FINDING, "status", None),
    "the coeliac screen": (_hla, "hla:coeliac", Kind.UNASSESSED, "tags", None),
    "HHV-6": (
        lambda: Hhv6Source(Hhv6Call("HHV-6B", 1.0, Kind.FINDING), None, REASON),
        "hhv6:ici",
        Kind.FINDING,
        "ratio",
        1.0,
    ),
    "runs of homozygosity": (
        lambda: RohSource(RohSummary(RUN.length, 1, RUN, (RUN,)), None, REASON),
        "roh:summary",
        Kind.FINDING,
        "n_long",
        1,
    ),
    "a blood group": (
        lambda: BloodGroupSource(
            calls={"FY": SystemCall("FY*02/FY*02", "Fy(a-b+)")},
            not_run=None,
            withheld={"FY": REASON},
            uncovered={},
        ),
        "bloodgroup:FY",
        Kind.FINDING,
        "genotype",
        "FY*02/FY*02",
    ),
    "the Neanderthal share": (
        _ancestry,
        "ancestry:neanderthal",
        Kind.FINDING,
        "share",
        0.011324,
    ),
    "a MITOMAP heteroplasmy": (
        _mitomap,
        "mitomap:3243A>G",
        Kind.FINDING,
        "level",
        0.30,
    ),
}


def _case(name: str) -> tuple[Source, str, Kind, str, object]:
    make, statement, would_be, key, value = CASES[name]
    return make(), statement, would_be, key, value


def _row(collected: Collected, statement_id: str) -> Any:
    (row,) = (r for r in collected.rows if r.statement["statement_id"] == statement_id)
    return row


def _claim(row: Any) -> dict[str, Any]:
    return json.loads(str(row.statement["claim"]))


def _without_withheld(document: Any) -> Any:
    """`document` with every `withheld` key dropped and its run id blanked."""
    if isinstance(document, dict):
        kept = cast(dict[str, Any], document)
        return {
            k: "" if k == "run_id" else _without_withheld(v)
            for k, v in kept.items()
            if k != "withheld"
        }
    if isinstance(document, list):
        return [_without_withheld(v) for v in cast(list[Any], document)]
    return document


@pytest.mark.parametrize("name", CASES)
def test_only_the_claim_changes(tmp_path: Path, name: str) -> None:
    source, *_ = _case(name)
    ctx = run_context(tmp_path)
    shown, revealed = source.collect(ctx), reveal(source, ctx)
    assert revealed.extra == shown.extra
    assert (revealed.tool, revealed.input_rows, revealed.dropped) == (
        shown.tool,
        shown.input_rows,
        shown.dropped,
    )
    assert len(revealed.rows) == len(shown.rows)
    for after, before in zip(revealed.rows, shown.rows, strict=True):
        assert (after.variant, after.call) == (before.variant, before.call)
        assert {**after.statement, "claim": ""} == {**before.statement, "claim": ""}
        assert _without_withheld(_claim(after)) == _without_withheld(_claim(before))


@pytest.mark.parametrize("name", CASES)
def test_the_withheld_entry(tmp_path: Path, name: str) -> None:
    source, statement, would_be, key, value = _case(name)
    (entry,) = _claim(_row(reveal(source, run_context(tmp_path)), statement))[
        "withheld"
    ]
    assert entry["reliability"] == Level.LOW
    assert entry["reason"].startswith(REASON)
    assert entry["would_be"] == would_be
    if value is None:
        assert entry["value"][key]
    else:
        assert entry["value"][key] == value


@pytest.mark.parametrize(
    ("source", "statement"),
    [
        (_repeats("HTT", withheld={"HTT": REASON}), "expansionhunter:HTT"),
        (SmnSource(None, None, None, REASON, 20.0), "smn:SMN1"),
        (YleafSource(None, None, REASON, 0.95, True), "yleaf:haplogroup"),
        (LpaSource(None, REASON, None), "lpa:LPA"),
        (Hhv6Source(None, None, REASON), "hhv6:ici"),
        (RohSource(None, None, REASON), "roh:summary"),
        (
            BloodGroupSource(
                calls={}, not_run=None, withheld={"FY": REASON}, uncovered={}
            ),
            "bloodgroup:FY",
        ),
        (
            HlaSource(
                {}, {}, {rsid: "unknown" for rsid in TAG_RSIDS}, None, {"HLA-B": REASON}
            ),
            "hla:HLA-B*27",
        ),
        (
            AncestrySource(
                None, None, None, not_run={}, withheld={"neanderthal": REASON}
            ),
            "ancestry:neanderthal",
        ),
    ],
)
def test_no_value_nothing_changes(
    tmp_path: Path, source: Source, statement: str
) -> None:
    ctx = run_context(tmp_path)
    assert _row(reveal(source, ctx), statement) == _row(source.collect(ctx), statement)


def test_a_tool_not_run_changes_nothing(tmp_path: Path) -> None:
    source = SmnSource(
        SmnCaller(1, 2, 0, "high"),
        ParascopyCall(1, 2, 40.0, True),
        "not calibrated",
        REASON,
        20.0,
    )
    ctx = run_context(tmp_path)
    assert reveal(source, ctx) == source.collect(ctx)


def test_nothing_withheld_changes_nothing(tmp_path: Path) -> None:
    source = replace(_repeats("HTT", _spanned(42)), withheld={})
    ctx = run_context(tmp_path)
    assert reveal(source, ctx) == source.collect(ctx)


def _stored(
    tmp_path: Path, source: Source, *, revealed: bool
) -> duckdb.DuckDBPyConnection:
    """A store holding `source`'s rows, as the run writes them, or as it did before."""
    tmp_path.mkdir()
    con = connect(tmp_path / "findings.duckdb")
    run_id = start_run(con, {}, {})
    ctx = replace(run_context(tmp_path), run_id=run_id)
    with transaction(con), StoreWriter(con, tmp_path / "work") as writer:
        if revealed:
            pipeline.store_sources(writer, con, ctx, [source])
        else:
            write_collected(writer, con, run_id, [source.collect(ctx)])
    return con


@pytest.mark.parametrize("name", CASES)
def test_readers_of_the_graded_result_see_what_they_saw(
    tmp_path: Path, name: str
) -> None:
    source, statement, _, _, _ = _case(name)
    before = _stored(tmp_path / "before", source, revealed=False)
    after = _stored(tmp_path / "after", source, revealed=True)

    def read(con: duckdb.DuckDBPyConnection) -> Any:
        run = str(query(con, None, LOOSEST).run_id)
        accounting = con.execute(
            "SELECT * EXCLUDE (run_id) FROM accounting ORDER BY ALL"
        ).fetchall()
        return _without_withheld(
            {
                "query": query(con, None, LOOSEST, list_unassessed=True).to_json(),
                "excluded": dashboard._excluded(con, run),  # pyright: ignore[reportPrivateUsage]
                "screened": dashboard._screened(con, run),  # pyright: ignore[reportPrivateUsage]
                "accounting": accounting,
            }
        )

    assert read(after) == read(before)
    _, shown = show(after, None, statement)
    assert shown is not None
    assert shown["kind"] == Kind.UNASSESSED
    assert shown["claim"]["withheld"]


@pytest.mark.parametrize(
    ("alleles", "exclusion_gene", "would_be"),
    [
        ((Allele(42, (41, 43), 6, 3, 0),), "HTT", Kind.FINDING),
        ((_spanned(17), _spanned(19)), "HTT", Kind.NON_FINDING),
        ((_spanned(30),), "FMR1", Kind.NON_FINDING),
    ],
)
def test_a_withheld_repeat_is_never_a_result(
    tmp_path: Path,
    alleles: tuple[Allele, ...],
    exclusion_gene: str,
    would_be: Kind,
) -> None:
    """The safety line: in the pathogenic range no finding or alert, in the normal range no
    normal and no exclusion, and `genotype` answers unknown either way."""
    source = _repeats(exclusion_gene, *alleles)
    con = _stored(tmp_path / "store", source, revealed=True)
    run = str(query(con, None, LOOSEST).run_id)
    assert dashboard._excluded(con, run) == []  # pyright: ignore[reportPrivateUsage]
    screened = dashboard._screened(con, run)  # pyright: ignore[reportPrivateUsage]
    assert screened == {"watch": [], "normal": []}
    assert exclusion_gene not in {f["gene"] for f in query(con, None, LOOSEST).findings}
    (site, *_) = genotype(con, None, gene=exclusion_gene).sites
    assert site["answer"] == Answer.UNKNOWN
    claim = site["claim"]
    # What graded readers read stays as unassessed: no class, no confirmation.
    assert (claim["class"], claim["confirm"]) == (None, None)
    assert site["gene_call"]["diplotype"] is None
    assert site["gene_call"]["level"] == Level.LOW
    assert not con.execute("SELECT count(*) FROM statement WHERE exclusion").fetchall()[
        0
    ][0]
    (entry,) = claim["withheld"]
    assert entry["would_be"] == would_be


def test_the_query_counts_the_values_a_table_withheld(tmp_path: Path) -> None:
    source = _repeats("HTT", Allele(42, (41, 43), 6, 3, 0))
    con = _stored(tmp_path / "store", source, revealed=True)
    by_gene = query(con, None, LOOSEST).unassessed["by_gene"]
    assert by_gene["HTT"]["withheld"] == 1
    assert by_gene["ATXN1"]["withheld"] == 0


def test_a_withheld_hla_gene_stays_out_of_the_table_and_the_explanation(
    tmp_path: Path,
) -> None:
    con = _stored(tmp_path / "store", _hla(), revealed=True)
    _, table = show(con, None, "hla:table")
    assert table is not None
    assert table["kind"] == Kind.FINDING
    held = next(g for g in table["claim"]["genes"] if g["gene"] == "HLA-B")
    assert (held["alleles"], held["callers"], held["note"]) == (None, None, REASON)
    entry = next(e for e in table["claim"]["withheld"] if e["value"]["gene"] == "HLA-B")
    assert entry["reliability"] == Level.LOW
    assert entry["reason"] == REASON
    assert entry["value"]["callers"]["T1K"] != "no call"
    # An explanation is written from the claim: a withheld value must not become one.
    (own, *_) = grounding(table, None)
    assert "withheld" not in json.loads(own.text)


CLASS = CallClass("loss_5kb+", "both", False, False)


def _rate(_: CallClass) -> None:
    return None


def _shown() -> Shown:
    d = Deletion("chr7", 117_530_000, 117_540_000, "N", "0/1", {"delly": "D1"}, 0.5)
    return Shown(d, CLASS, {"CFTR": (5, 6)}, None, kinds={"CFTR": "loss"})


def test_a_deletion_of_a_band_not_measured_keeps_its_grade(tmp_path: Path) -> None:
    source = DeletionSource(
        shown=[_shown()],
        called=1,
        dropped={},
        not_run=None,
        p_wrong=_rate,
        haploinsufficient=frozenset(),
        recall={},
    )
    row = next(
        r
        for r in source.collect(run_context(tmp_path)).rows
        if r.statement["gene"] == "CFTR"
    )
    assert row.statement["kind"] == Kind.UNASSESSED
    (entry,) = _claim(row)["withheld"]
    assert entry["reason"] == row.statement["reason"]
    assert entry["would_be"] == Kind.FINDING
    assert entry["value"]["start"] == 117_530_001
    assert "needs confirmation" in entry["significance"]
    assert entry["confirm"]


def _cnv(**kw: Any) -> CnvSource:
    fields: dict[str, Any] = {
        "shown": [],
        "large": [],
        "named": [],
        "called": 1,
        "dropped": {},
        "not_run": None,
        "p_wrong": _rate,
        "haploinsufficient": frozenset(),
        "triplosensitive": frozenset(),
        "deletion_cautions": lambda _: [],
        "cautions": lambda _: [],
    }
    return CnvSource(**(fields | kw))


def test_a_large_event_of_a_withheld_band_stays_in_the_screen(tmp_path: Path) -> None:
    d = Deletion("chr7", 1_000, 2_000_000, "N", "0/1", {"dellycnv": "C1"}, None)
    large = LargeEvent(d, CallClass("loss_large", "both", False, False), ("GENE",))
    collected = _cnv(large=[large], large_withheld={"loss_large": REASON}).collect(
        run_context(tmp_path)
    )
    (row,) = collected.rows
    assert row.statement["statement_id"] == "cnv:screen"
    assert row.statement["kind"] == Kind.NON_FINDING
    claim = _claim(row)
    assert claim["shown"] == 0
    (entry,) = claim["withheld"]
    assert (entry["reason"], entry["reliability"]) == (REASON, Level.LOW)
    assert entry["would_be"] == Kind.FINDING
    assert entry["value"]["length"] == d.length
    assert entry["value"]["regions"] == ["GENE"]
    assert entry["confirm"] == "a chromosomal microarray"


def test_no_large_event_withheld_leaves_the_screen_as_it_was(tmp_path: Path) -> None:
    (row,) = _cnv().collect(run_context(tmp_path)).rows
    assert "withheld" not in _claim(row)


def test_a_named_locus_of_a_band_not_measured_keeps_its_grade(tmp_path: Path) -> None:
    named = NamedCnv("GJB6-D13S1830", "chr13", 1, 2, "GJB6", "one", "both")
    row = next(
        r
        for r in _cnv(named=[(named, _shown())]).collect(run_context(tmp_path)).rows
        if r.statement["statement_id"] == "cnv:named:GJB6-D13S1830"
    )
    assert row.statement["kind"] == Kind.UNASSESSED
    (entry,) = _claim(row)["withheld"]
    assert (entry["would_be"], entry["significance"]) == (Kind.FINDING, "one")
    assert entry["value"]["name"] == "GJB6-D13S1830"


def test_a_withheld_heteroplasmy_keeps_its_level_after_the_reason(
    tmp_path: Path,
) -> None:
    """Any chrM claim's, a ClinVar one's too: `query` prints warnings, and reads the first."""
    called = mito_called(
        ("chrM", 3243, "A", "G"),
        [MitoCall(pos=3243, ref="A", alt="G", level=0.30, dp=3000)],
        {},
        min_dp=10,
        min_homoplasmy=0.9,
        heteroplasmy_withheld=REASON,
    )
    assert called.dosage == "unknown"
    assert called.warnings == (
        f"genotype unknown: {REASON}",
        "heteroplasmy at level 0.3",
    )
