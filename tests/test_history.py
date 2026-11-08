"""Snapshot history: which statements changed between two runs, why, and their dates."""

import itertools
import time
from collections.abc import Iterator, Sequence
from datetime import date
from pathlib import Path
from typing import Any

import pytest
from fixtures import snapshots_mini as sm
from fixtures.snapshots_mini import call, file, gene_call, run_id, started, statement

from genokeeper import history
from genokeeper.dashboard import _day
from genokeeper.history import (
    COMPARED,
    REFERENCE_CAUSE,
    ChangeKind,
    Dates,
    dates,
    diff,
    kind_of,
    snapshots,
    states,
    trusted_ids,
)
from genokeeper.query import Filters, query
from genokeeper.store.db import Level, Validity, export_snapshot

ANALYSER = "a" * 64
TOOLS = {"clinvar": "20261011", "cpic": "20261001", "analyser": ANALYSER}
INPUTS: dict[str, Any] = {
    "fasta": file("fa"),
    "clinvar38": file("cv1"),
    "gene_validity": {**file("gv1"), "date": "20261001"},
    "pgx": {
        "positions": file("pos1"),
        "cpic": {"a.json": file("c1"), "b.json": file("c2")},
        "pharmgkb": {**file("pg1"), "created": "2026-09-01"},
    },
    "prs": {"scores": {"PGS000018": file("s1")}, "panel": file("p1")},
    "steps": {"align": {"key": "k1", "path": "/w", "reused": False}},
}
D1, D2, D3 = run_id(1), run_id(2), run_id(3)


def write(
    root: Path,
    rid: str,
    statements: list[dict[str, Any]],
    calls: Sequence[dict[str, Any]] = (),
    gene_calls: Sequence[dict[str, Any]] = (),
    *,
    tools: dict[str, Any] | None = None,
    inputs: dict[str, Any] | None = None,
) -> None:
    sm.write(
        root,
        rid,
        started(rid),
        TOOLS if tools is None else tools,
        INPUTS if inputs is None else inputs,
        statements,
        calls,
        gene_calls,
    )


def vid(pos: int) -> str:
    return f"GRCh38:chr17:{pos}:C:T"


def at(sid: str, pos: int, **kw: Any) -> dict[str, Any]:
    """A statement at its own position."""
    return statement(sid, vid=vid(pos), **kw)


def kinds(root: Path, base: str, run: str, f: Filters | None = None) -> dict[str, str]:
    f = f or Filters()
    before, after = trusted_ids(root, base, f), trusted_ids(root, run, f)
    return {
        c.statement_id: str(
            kind_of(c, f, c.statement_id in before, c.statement_id in after)
        )
        for c in diff(root, base, run)
    }


def _feed_rows() -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    """A GWAS row, an unassessed claim, a gene-level finding without a gene call, one with
    a low gene call, and findings at each level and validity."""
    statements = [
        at("clinvar:1", 1001),
        at("clinvar:2", 1002, validity="moderate"),
        at("clinvar:3", 1003, validity="limited"),
        at("clinvar:4", 1004, validity="conflicting"),
        at("clinvar:5", 1005),
        at("clinvar:6", 1006),
        at(
            "gwas:rs1-T:EFO_1",
            1007,
            gene="PPARG",
            rsid="rs1",
            source="gwas_catalog",
            topic="research",
        ),
        at("clinvar:8", 1008, kind="unassessed"),
        statement("pharmcat:CYP2C19", vid=None, gene="CYP2C19", source="pharmcat"),
        statement("pharmcat:CYP2D6", vid=None, gene="CYP2D6", source="pharmcat"),
        statement("pgsc_calc:PGS1", vid=None, gene=None, source="pgsc_calc"),
        at("clinvar:9", 1009, kind="watch", validity=None, topic=None),
    ]
    calls = [
        call(vid(1001)),
        call(vid(1002), level="medium"),
        call(vid(1003)),
        call(vid(1004)),
        call(vid(1005), level="medium"),
        call(vid(1006), level="low"),
        call(vid(1007)),
        call(
            vid(1008), dosage="unknown", level="low", warnings=["genotype unknown: x"]
        ),
        call(vid(1009)),
    ]
    return statements, calls, [gene_call("CYP2D6", level="low")]


def test_trusted_ids_match_query_for_every_slider_pair(tmp_path: Path) -> None:
    statements, calls, gene_calls = _feed_rows()
    con = sm.store(tmp_path / "genomelab.duckdb")
    sm.fill(con, D1, started(D1), TOOLS, INPUTS, statements, calls, gene_calls)
    root = tmp_path / "snapshots"
    export_snapshot(con, D1, root)
    for rel, val in itertools.product(Level, Validity):
        for gene in (None, "PPARG"):
            f = Filters(min_reliability=rel, min_validity=val, gene=gene)
            expected = {x["id"] for x in query(con, D1, f).findings}
            assert trusted_ids(root, D1, f) == expected, (rel, val, gene)
    assert trusted_ids(root, D1, Filters()) == {
        "clinvar:1",
        "clinvar:2",
        "clinvar:5",
        "pharmcat:CYP2C19",
        "pgsc_calc:PGS1",
    }


def test_each_change_kind(tmp_path: Path) -> None:
    root = tmp_path
    base = [
        at("clinvar:1", 1001, validity="limited"),
        at("clinvar:2", 1002),
        at("clinvar:3", 1003),
        at("clinvar:5", 1005, kind="watch", validity=None, topic=None),
        at("clinvar:6", 1006),
        at("clinvar:8", 1008, kind="non_finding"),
    ]
    run = [
        at("clinvar:1", 1001),
        at("clinvar:2", 1002, validity="conflicting"),
        at("clinvar:3", 1003, validity="moderate"),
        at("clinvar:4", 1004, kind="watch", validity=None, topic=None),
        at("clinvar:6", 1006, source_level="Pathogenic (3★)"),
        at("clinvar:7", 1007),
        at("clinvar:8", 1008, kind="non_finding", validity="limited"),
        at("clinvar:9", 1009, kind="non_finding"),
    ]
    calls = [call(vid(p)) for p in range(1001, 1010)]
    write(root, D1, base, calls)
    write(root, D2, run, calls)
    assert kinds(root, D1, D2) == {
        "clinvar:1": ChangeKind.ENTERED,
        "clinvar:2": ChangeKind.LEFT,
        "clinvar:3": ChangeKind.LEVEL,
        "clinvar:4": ChangeKind.NEW,
        "clinvar:5": ChangeKind.GONE,
        "clinvar:6": ChangeKind.DETAILS,
        "clinvar:7": ChangeKind.ENTERED,
    }
    (six,) = (c for c in diff(root, D1, D2) if c.statement_id == "clinvar:6")
    assert six.fields == ("source_level",)


def test_a_finding_that_becomes_a_non_finding_keeps_its_after_state(
    tmp_path: Path,
) -> None:
    """A reclassification to benign leaves the trusted tier with its new classification;
    a non-finding that appears or goes stays invisible."""
    calls = [call(vid(p)) for p in (1001, 1002, 1003, 1004)]
    watch = {"kind": "watch", "validity": None}
    benign = {"kind": "non_finding", "source_level": "Benign (2★)"}
    nothing = {"kind": "non_finding"}
    base = [
        at("clinvar:1", 1001),
        at("clinvar:2", 1002, **watch),
        at("clinvar:3", 1003, **nothing),
    ]
    write(tmp_path, D1, base, calls)
    run = [
        at("clinvar:1", 1001, **benign),
        at("clinvar:2", 1002, **{**watch, **benign}),
    ]
    write(tmp_path, D2, [*run, at("clinvar:4", 1004, **nothing)], calls)
    found = {c.statement_id: c for c in diff(tmp_path, D1, D2)}
    assert set(found) == {"clinvar:1", "clinvar:2"}
    left = found["clinvar:1"]
    assert left.after is not None
    assert (left.after.kind, left.after.source_level) == ("non_finding", "Benign (2★)")
    assert left.fields == ("kind", "source_level")
    assert kinds(tmp_path, D1, D2) == {
        "clinvar:1": ChangeKind.LEFT,
        "clinvar:2": ChangeKind.GONE,
    }


def test_research_associations_never_enter_or_leave(tmp_path: Path) -> None:
    row = {
        "gene": "PPARG",
        "rsid": "rs1",
        "source": "gwas_catalog",
        "topic": "research",
    }
    write(
        tmp_path,
        D1,
        [at("gwas:rs1-T:E", 1001, validity="limited", **row)],
        [call(vid(1001))],
    )
    write(tmp_path, D2, [at("gwas:rs1-T:E", 1001, **row)], [call(vid(1001))])
    named = Filters(gene="PPARG")
    # Named, the query keeps it once it passes, so the trust sets differ ...
    assert "gwas:rs1-T:E" in trusted_ids(tmp_path, D2, named)
    assert "gwas:rs1-T:E" not in trusted_ids(tmp_path, D1, named)
    # ... and still it neither enters nor leaves.
    assert kinds(tmp_path, D1, D2, named) == {"gwas:rs1-T:E": ChangeKind.LEVEL}
    assert kinds(tmp_path, D1, D2) == {"gwas:rs1-T:E": ChangeKind.LEVEL}


RESEARCH_ROW: dict[str, Any] = {
    "gene": "RGS4",
    "rsid": "rs951439",
    "validity": "limited",
    "source": "pharmgkb_research",
    "source_level": "ClinPGx level 3",
    "topic": "research",
}


def test_level_3_never_enters_change_feed(tmp_path: Path) -> None:
    """A new level-3 PharmGKB annotation is a research association: even named, and so
    trusted in the run, it is new, never entered."""
    write(tmp_path, D1, [at("clinvar:1", 1001)], [call(vid(1001))])
    write(
        tmp_path,
        D2,
        [
            at("clinvar:1", 1001),
            at("pharmgkb_research:655384602", 1002, **RESEARCH_ROW),
        ],
        [call(vid(1001)), call(vid(1002))],
    )
    named = Filters(
        gene="RGS4", min_reliability=Level.LOW, min_validity=Validity.CONFLICTING
    )
    assert "pharmgkb_research:655384602" in trusted_ids(tmp_path, D2, named)
    for f in (named, None):
        assert kinds(tmp_path, D1, D2, f) == {
            "pharmgkb_research:655384602": ChangeKind.NEW
        }


def test_version_only_change_is_no_change(tmp_path: Path) -> None:
    write(tmp_path, D1, [at("clinvar:1", 1001)], [call(vid(1001))])
    write(
        tmp_path,
        D2,
        [at("clinvar:1", 1001, source_version="20261018")],
        [call(vid(1001))],
        tools={**TOOLS, "clinvar": "20261018"},
        inputs={**INPUTS, "clinvar38": file("cv2")},
    )
    assert diff(tmp_path, D1, D2) == []


def test_clinvar_id_change_at_same_vid_is_one_statement(tmp_path: Path) -> None:
    write(
        tmp_path,
        D1,
        [at("clinvar:10", 1001), at("clinvar:20", 1002)],
        [call(vid(1001)), call(vid(1002))],
    )
    write(
        tmp_path,
        D2,
        [at("clinvar:11", 1001), at("clinvar:21", 1002, validity="moderate")],
        [call(vid(1001)), call(vid(1002))],
    )
    (change,) = diff(tmp_path, D1, D2)
    assert change.statement_id == "clinvar:21"
    assert change.before is not None and change.before.statement_id == "clinvar:20"
    assert change.fields == ("validity",)
    found = dates(tmp_path)
    assert found["clinvar:11"].first_found == date(2026, 10, 1)
    assert "clinvar:10" not in found


def test_three_dates_over_three_snapshots(tmp_path: Path) -> None:
    root = tmp_path
    d2, d3 = started(D2).date(), started(D3).date()
    calls = [call(vid(1001))]
    vus = {"kind": "watch", "validity": "limited", "source_level": "Uncertain (2★)"}
    write(root, D1, [at("clinvar:1", 1001, **vus)], calls)
    write(
        root, D2, [at("clinvar:1", 1001, source_level="Likely pathogenic (2★)")], calls
    )
    write(
        root, D3, [at("clinvar:1", 1001, source_level="Likely pathogenic (3★)")], calls
    )
    assert dates(root)["clinvar:1"] == Dates(
        first_found=d2, last_changed=d3, assessment_changed=d2, on_or_before=False
    )


def test_snapshot_without_compared_columns_is_skipped(tmp_path: Path) -> None:
    root = tmp_path
    old, older = run_id(28, month=9), run_id(27, month=9)
    rows, calls = [at("clinvar:1", 1001, validity="moderate")], [call(vid(1001))]
    write(root, older, rows, calls)
    write(root, old, rows, calls)
    write(root, D1, [at("clinvar:1", 1001)], calls)
    write(root, D2, [at("clinvar:1", 1001)], calls)
    # An older schema wrote `old` without `statement.reason`, `older` without gene calls.
    sm.store().execute(
        f"COPY (SELECT * EXCLUDE (reason) FROM read_parquet('{root / old / 'statement.parquet'}'))"
        f" TO '{root / old / 'statement.parquet'}' (FORMAT parquet)"
    )
    (root / older / "gene_call.parquet").unlink()
    assert [s.run_id for s in snapshots(root)] == [D1, D2]
    assert dates(root)["clinvar:1"] == Dates(
        first_found=date(2026, 10, 1),
        last_changed=date(2026, 10, 1),
        assessment_changed=date(2026, 10, 1),
        on_or_before=True,
    )


def test_causes_are_the_union_over_intermediate_pairs(tmp_path: Path) -> None:
    root = tmp_path
    calls = [call(vid(1001)), call(vid(1002))]
    write(root, D1, [at("clinvar:1", 1001), at("clinvar:2", 1002)], calls)
    write(
        root,
        D2,
        [at("clinvar:1", 1001, validity="moderate"), at("clinvar:2", 1002)],
        calls,
        tools={**TOOLS, "analyser": "b" * 64},
    )
    write(
        root,
        D3,
        [
            at("clinvar:1", 1001, validity="moderate", source_level="Pathogenic (3★)"),
            at("clinvar:2", 1002, reason="BA1"),
        ],
        calls,
        tools={**TOOLS, "analyser": "b" * 64, "clinvar": "20261018"},
        inputs={**INPUTS, "clinvar38": file("cv2")},
    )
    causes = {c.statement_id: c.causes for c in diff(root, D1, D3)}
    assert causes == {
        "clinvar:1": ("analyser", "source:clinvar38@20261018"),
        "clinvar:2": ("source:clinvar38@20261018",),
    }


def test_gencc_cause_in_change_feed(tmp_path: Path) -> None:
    """Two runs whose inputs differ only by GenCC: a changed ClinVar statement, and a
    changed predicted one, name GenCC's release as the cause."""
    calls = [call(vid(1001)), call(vid(1002))]
    gencc = {**file("gc1"), "date": "20261025"}
    write(
        tmp_path,
        D1,
        [at("clinvar:1", 1001), at("vep:2", 1002, source="vep")],
        calls,
        inputs={**INPUTS, "gencc": gencc},
    )
    write(
        tmp_path,
        D2,
        [
            at(
                "clinvar:1",
                1001,
                validity="limited",
                reason="gene validity Limited (GenCC)",
            ),
            at("vep:2", 1002, source="vep", validity="limited"),
        ],
        calls,
        inputs={**INPUTS, "gencc": {**file("gc2"), "date": "20261101"}},
    )
    causes = {c.statement_id: c.causes for c in diff(tmp_path, D1, D2)}
    assert causes == {
        "clinvar:1": ("source:gencc@20261101",),
        "vep:2": ("source:gencc@20261101",),
    }


def test_snapshot_lacking_the_analyser_differs_from_any(tmp_path: Path) -> None:
    calls = [call(vid(1001))]
    no_analyser = {k: v for k, v in TOOLS.items() if k != "analyser"}
    write(tmp_path, D1, [at("clinvar:1", 1001)], calls, tools=no_analyser)
    write(tmp_path, D2, [at("clinvar:1", 1001, reason="x")], calls)
    assert [c.causes for c in diff(tmp_path, D1, D2)] == [("analyser",)]


def test_nested_input_digest_names_cpic(tmp_path: Path) -> None:
    root = tmp_path
    row = {"vid": None, "gene": "CYP2C19", "source": "pharmcat", "topic": "drugs"}
    write(
        root,
        D1,
        [statement("pharmcat:CYP2C19", **row)],
        gene_calls=[gene_call("CYP2C19")],
    )
    pgx = INPUTS["pgx"]
    write(
        root,
        D2,
        [statement("pharmcat:CYP2C19", **row, validity="moderate")],
        gene_calls=[gene_call("CYP2C19")],
        tools={**TOOLS, "cpic": "20261020"},
        inputs={
            **INPUTS,
            "pgx": {
                **pgx,
                "cpic": {**pgx["cpic"], "b.json": file("c3")},
                "pharmgkb": file("pg2"),
            },
        },
    )
    assert [c.causes for c in diff(root, D1, D2)] == [("source:pgx.cpic@20261020",)]


def test_unknown_source_names_every_changed_input(tmp_path: Path) -> None:
    root = tmp_path
    row = {"vid": None, "gene": None, "source": "pgsc_calc"}
    write(root, D1, [statement("pgsc_calc:PGS1", **row)])
    write(
        root,
        D2,
        [statement("pgsc_calc:PGS1", **row, effect={"score": {"percentile": 90.0}})],
        tools={**TOOLS, "clinvar": "20261018"},
        inputs={
            **INPUTS,
            "clinvar38": file("cv2"),
            "prs": {**INPUTS["prs"], "scores": {"PGS000018": file("s2")}},
            "steps": {"align": {"key": "k1", "path": "/w", "reused": True}},
        },
    )
    ((first, second),) = [c.causes for c in diff(root, D1, D2)]
    assert first == "source:clinvar38@20261018"
    assert second.startswith("source:prs@") and len(second) == len("source:prs@") + 12


def test_reference_change_is_a_cause(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Two snapshots identical but for what the masked reference changes (the FASTA's md5
    and `inputs["fasta"]`) and one call's level, on a statement of a source
    `SOURCE_INPUTS` does not list: the one cause is the reference, not the FASTA input."""
    row: dict[str, Any] = {"vid": vid(1001), "gene": "CBS", "source": "expansionhunter"}

    def snapshots_in(root: Path) -> None:
        write(
            root,
            D1,
            [statement("expansionhunter:CBS", **row)],
            [call(vid(1001))],
            tools={**TOOLS, "reference_md5": "m1"},
        )
        write(
            root,
            D2,
            [statement("expansionhunter:CBS", **row)],
            [call(vid(1001), level="medium")],
            tools={**TOOLS, "reference_md5": "m2"},
            inputs={**INPUTS, "fasta": file("fb")},
        )

    snapshots_in(tmp_path / "a")
    (change,) = diff(tmp_path / "a", D1, D2)
    assert change.causes == (REFERENCE_CAUSE,) == ("reference",)
    # With `fasta` read as an input, the same pair names it too: the assertion can fail.
    monkeypatch.setattr(history, "_NOT_INPUTS", history._NOT_INPUTS - {"fasta"})  # pyright: ignore[reportPrivateUsage]
    snapshots_in(tmp_path / "b")
    (change,) = diff(tmp_path / "b", D1, D2)
    reference, fasta = change.causes
    assert reference == "reference" and fasta.startswith("source:fasta@")


def test_genotype_cause_when_only_the_call_moved(tmp_path: Path) -> None:
    write(tmp_path, D1, [at("clinvar:1", 1001)], [call(vid(1001))])
    write(tmp_path, D2, [at("clinvar:1", 1001)], [call(vid(1001), level="medium")])
    (change,) = diff(tmp_path, D1, D2)
    assert change.fields == ("level",) and change.causes == ("genotype",)


def test_states_carry_every_compared_field(tmp_path: Path) -> None:
    effect = {"relative": {"value": 2.0, "measure": "OR"}}
    write(
        tmp_path,
        D1,
        [at("clinvar:1", 1001, rsid="rs9", effect=effect)],
        [call(vid(1001), dosage="2")],
    )
    state = states(tmp_path, D1)["clinvar:1"]
    assert {f: getattr(state, f) for f in COMPARED} == {
        "kind": "finding",
        "validity": "strong",
        "source_level": "Pathogenic (2★)",
        "reason": None,
        "effect": '{"relative": {"measure": "OR", "value": 2.0}}',
        "actionable": False,
        "topic": "risks",
        "dosage": "2",
        "level": "high",
        "result": None,
    }
    assert (state.rsid, state.source, state.source_version) == (
        "rs9",
        "clinvar",
        "20261011",
    )


def test_trusted_ids_read_a_snapshot_without_newer_columns(tmp_path: Path) -> None:
    write(tmp_path, D1, [at("clinvar:1", 1001)], [call(vid(1001))])
    calls = tmp_path / D1 / "call.parquet"
    sm.store().execute(
        f"COPY (SELECT * EXCLUDE (p_wrong_lo, p_wrong_hi) FROM read_parquet('{calls}'))"
        f" TO '{calls}' (FORMAT parquet)"
    )
    assert trusted_ids(tmp_path, D1, Filters()) == {"clinvar:1"}


def _y(name: str) -> dict[str, Any]:
    return statement(
        "yleaf:haplogroup",
        vid=None,
        gene=None,
        source="yleaf",
        topic="ancestry",
        claim={"type": "y_haplogroup", "haplogroup": name},
    )


def test_a_gene_level_result_change_is_a_change(tmp_path: Path) -> None:
    """A diplotype and phenotype, or a haplogroup, that moved under the same validity and
    level is a change of the user's result, dated."""
    pgx = statement("pharmcat:CYP2C19", vid=None, gene="CYP2C19", source="pharmcat")
    poor = {
        **gene_call("CYP2C19"),
        "diplotype": "*2/*2",
        "phenotype": "Poor Metabolizer",
    }
    write(tmp_path, D1, [pgx, _y("R1b-M269")], gene_calls=[gene_call("CYP2C19")])
    write(tmp_path, D2, [pgx, _y("I2a-M423")], gene_calls=[poor])
    changed = {c.statement_id: c for c in diff(tmp_path, D1, D2)}
    assert set(changed) == {"pharmcat:CYP2C19", "yleaf:haplogroup"}
    pharmcat, y = changed["pharmcat:CYP2C19"], changed["yleaf:haplogroup"]
    assert pharmcat.fields == y.fields == ("result",)
    assert pharmcat.before is not None and pharmcat.after is not None
    assert pharmcat.before.result == "*1/*2 Intermediate Metabolizer"
    assert pharmcat.after.result == "*2/*2 Poor Metabolizer"
    assert y.after is not None and y.after.result == "I2a-M423"
    found = dates(tmp_path)
    for sid in changed:
        assert found[sid].last_changed == date(2026, 10, 2), sid


def test_two_tools_on_one_gene_change_apart(tmp_path: Path) -> None:
    """CFTR's poly-T/TG call moved under PharmCAT's unchanged one: the tracts' statement
    alone is a change, with its own result."""
    pgx = statement("pharmcat:CFTR", vid=None, gene="CFTR", source="pharmcat")
    tracts = statement(
        "expansionhunter:CFTR-polyT-TG",
        vid=None,
        gene="CFTR",
        source="expansionhunter",
    )
    before = {
        **gene_call("CFTR-polyT-TG"),
        "diplotype": "TG 11/12 T 7/9",
        "phenotype": "no 5T",
    }
    moved = {**before, "diplotype": "TG 11/12 T 5/9", "phenotype": "5T"}
    write(tmp_path, D1, [pgx, tracts], gene_calls=[gene_call("CFTR"), before])
    write(tmp_path, D2, [pgx, tracts], gene_calls=[gene_call("CFTR"), moved])
    (change,) = diff(tmp_path, D1, D2)
    assert change.statement_id == "expansionhunter:CFTR-polyT-TG"
    assert change.after is not None and change.after.result == "TG 11/12 T 5/9 5T"


def test_a_change_and_back_is_no_change(tmp_path: Path) -> None:
    calls = [call(vid(1001)), call(vid(1002))]
    write(tmp_path, D1, [at("clinvar:1", 1001), at("clinvar:2", 1002)], calls)
    write(tmp_path, D2, [at("clinvar:1", 1001, validity="moderate")], calls)
    write(tmp_path, D3, [at("clinvar:1", 1001), at("clinvar:2", 1002)], calls)
    assert [c.statement_id for c in diff(tmp_path, D1, D2)] == [
        "clinvar:1",
        "clinvar:2",
    ]
    assert diff(tmp_path, D1, D3) == []
    assert dates(tmp_path)["clinvar:2"].last_changed == date(2026, 10, 3)


def test_effect_keys_in_another_order_are_no_change(tmp_path: Path) -> None:
    calls = [call(vid(1001))]
    effect = {"relative": {"measure": "OR", "value": 2.0}, "absolute": None}
    write(tmp_path, D1, [at("clinvar:1", 1001, effect=effect)], calls)
    reordered = {"absolute": None, "relative": {"value": 2.0, "measure": "OR"}}
    write(tmp_path, D2, [at("clinvar:1", 1001, effect=reordered)], calls)
    assert diff(tmp_path, D1, D2) == []
    assert dates(tmp_path)["clinvar:1"].last_changed == date(2026, 10, 1)


@pytest.fixture
def utc_plus_3(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """A time zone east of UTC: UTC+3 (POSIX writes it Etc/GMT-3)."""
    monkeypatch.setenv("TZ", "Etc/GMT-3")
    time.tzset()
    yield
    monkeypatch.undo()
    time.tzset()


@pytest.mark.usefixtures("utc_plus_3")
def test_dates_are_the_local_calendar_day(tmp_path: Path) -> None:
    """A run started at 01:30 at UTC+3 (22:30 UTC the day before) is dated that day."""
    rid = "20261011T223000Z-aaaaaa"
    write(tmp_path, rid, [at("clinvar:1", 1001)], [call(vid(1001))])
    day = date(2026, 10, 12)
    assert dates(tmp_path)["clinvar:1"] == Dates(day, day, day, on_or_before=True)
    assert _day(snapshots(tmp_path)[0]) == "2026-10-12"
