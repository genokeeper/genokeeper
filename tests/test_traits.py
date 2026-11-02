"""Five trait statements (slice 3f, "Traits"): one per trait, not per rsid, from our graded
genotypes at `traits.tsv`'s alleles (GRCh38, forward strand, dbSNP; the pre-flight
record's "Traits" table)."""

import json
from dataclasses import replace
from pathlib import Path

from fixtures.run_context import run_context

from genokeeper import pipeline
from genokeeper.genotype.gvcf import NOCALL_REASON, Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.query import Filters, query
from genokeeper.sources import INFORMATION_NOTE, write_collected
from genokeeper.sources.traits import (
    RARE_TAS2R38,
    TRAIT_ALLELES,
    TRAITS_TOOL,
    Trait,
    TraitSource,
    load_traits,
)
from genokeeper.store.db import Level, RowType, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

RELIABLE = Reliability(p_wrong=1e-6, level=Level.HIGH, flags=(), warnings=())
LACTASE = {
    "rs4988235": ("chr2", 135851076, "G", "A"),
    "rs41525747": ("chr2", 135851073, "G", "C"),
    "rs41380347": ("chr2", 135851081, "A", "C"),
    "rs869051967": ("chr2", 135851175, "A", "C"),
    "rs145946881": ("chr2", 135851176, "C", "G"),
}
ADH1B: Allele = ("chr4", 99318162, "T", "C")  # effect allele T, the REF (His48)
ACTN3: Allele = ("chr11", 66560624, "C", "T")  # effect allele T (577X)
CYP1A2: Allele = ("chr15", 74749576, "C", "A")  # effect allele A (inducible)
# (rs713598, rs1726866, rs10246939); PAV is forward (G, G, C), so rs1726866's PAV base is
# the REF; the GRCh38 reference haplotype (C, G, T) is AAI.
TAS2R38: tuple[Allele, ...] = (
    ("chr7", 141973545, "C", "G"),
    ("chr7", 141972905, "G", "A"),
    ("chr7", 141972804, "T", "C"),
)


def _called(dosage: str) -> tuple[Called, Reliability]:
    if dosage == "unknown":
        called = Called(
            "./.", "unknown", RowType.NOCALL, None, None, None, None, (NOCALL_REASON,)
        )
    else:
        gt = {"0": "0/0", "1": "0/1", "2": "1/1"}[dosage]
        called = Called(gt, dosage, RowType.PASS, 40, 30, None, None, ())
    return called, RELIABLE


def _graded(dosages: dict[Allele, str]) -> dict[Allele, tuple[Called, Reliability]]:
    """Every trait allele at dosage 0 (ALT absent), except those given."""
    return {a: _called(dosages.get(a, "0")) for a in TRAIT_ALLELES}


def _statements(
    tmp_path: Path, dosages: dict[Allele, str]
) -> dict[str, dict[str, object]]:
    collected = TraitSource(graded=_graded(dosages)).collect(run_context(tmp_path))
    assert collected.tool == TRAITS_TOOL == "traits"
    assert collected.input_rows == len(collected.rows) == 5
    out: dict[str, dict[str, object]] = {}
    for row in collected.rows:
        s = row.statement
        assert (s["topic"], s["validity"], s["actionable"], s["exclusion"]) == (
            "traits",
            "moderate",
            False,
            False,
        )
        assert (row.variant, row.call, s["vid"]) == (None, None, None)
        out[str(s["statement_id"])] = s
    assert set(out) == {
        "trait:lactase",
        "trait:ADH1B",
        "trait:ACTN3",
        "trait:TAS2R38",
        "trait:CYP1A2",
    }
    return out


def test_traits_file_holds_the_preflight_alleles() -> None:
    traits = {t.name: t for t in load_traits()}
    assert set(traits) == {"lactase", "ADH1B", "ACTN3", "TAS2R38", "CYP1A2"}
    lactase = traits["lactase"]
    assert (lactase.gene, lactase.rule) == ("LCT", "any")
    assert lactase.alleles == tuple(LACTASE.values())
    assert lactase.effects == tuple(a[3] for a in LACTASE.values())
    assert (traits["ADH1B"].alleles, traits["ADH1B"].effects) == ((ADH1B,), ("T",))
    assert (traits["ACTN3"].alleles, traits["ACTN3"].effects) == ((ACTN3,), ("T",))
    assert (traits["CYP1A2"].alleles, traits["CYP1A2"].effects) == ((CYP1A2,), ("A",))
    assert (traits["TAS2R38"].alleles, traits["TAS2R38"].effects) == (
        TAS2R38,
        ("G", "G", "C"),
    )
    assert traits["TAS2R38"].rule == "tas2r38"
    assert lactase.wordings[0] == (
        "lactase persistence unlikely: dairy intolerance in adulthood is common"
    )
    assert lactase.wordings[1] == lactase.wordings[2] == "lactase persistence likely"
    assert lactase.note == (
        "validated in Europeans; other enhancer variants in Middle Eastern and African"
        " ancestry"
    )
    assert set(TRAIT_ALLELES) == {a for t in traits.values() for a in t.alleles}
    assert all(isinstance(t, Trait) for t in traits.values())


def test_lactase_is_one_or_statement(tmp_path: Path) -> None:
    """rs4988235 carried, the other enhancers not: one `trait:lactase` "likely", no
    statement per rsid; none carried, all called: "unlikely"; none carried, one unknown:
    unassessed."""
    likely = _statements(tmp_path, {LACTASE["rs4988235"]: "1"})["trait:lactase"]
    assert likely["kind"] == "finding"
    assert str(likely["reason"]).startswith("lactase persistence likely")
    assert "validated in Europeans" in str(likely["reason"])
    east_africa = _statements(tmp_path, {LACTASE["rs145946881"]: "1"})["trait:lactase"]
    assert str(east_africa["reason"]).startswith("lactase persistence likely")
    unlikely = _statements(tmp_path, {})["trait:lactase"]
    assert unlikely["kind"] == "finding"
    assert str(unlikely["reason"]).startswith("lactase persistence unlikely")
    unknown = _statements(tmp_path, {LACTASE["rs41380347"]: "unknown"})["trait:lactase"]
    assert unknown["kind"] == "unassessed"
    assert "rs41380347" in str(unknown["reason"])
    carried_beside_unknown = _statements(
        tmp_path, {LACTASE["rs4988235"]: "2", LACTASE["rs41380347"]: "unknown"}
    )["trait:lactase"]
    assert carried_beside_unknown["kind"] == "finding"


def test_tas2r38_diplotype(tmp_path: Path) -> None:
    """Three homozygous PAV: PAV/PAV; three heterozygous: PAV/AVI; three homozygous AVI:
    AVI/AVI; one heterozygous beside two homozygous: unassessed, a rare haplotype. The
    reference (C, G, T) is AAI, never read as a non-taster."""
    pav_pav = {TAS2R38[0]: "2", TAS2R38[1]: "0", TAS2R38[2]: "2"}
    het = {TAS2R38[0]: "1", TAS2R38[1]: "1", TAS2R38[2]: "1"}
    avi_avi = {TAS2R38[0]: "0", TAS2R38[1]: "2", TAS2R38[2]: "0"}
    rare = {TAS2R38[0]: "1", TAS2R38[1]: "0", TAS2R38[2]: "2"}
    reference = {TAS2R38[0]: "0", TAS2R38[1]: "0", TAS2R38[2]: "0"}
    for dosages, diplotype in (
        (pav_pav, "PAV/PAV"),
        (het, "PAV/AVI"),
        (avi_avi, "AVI/AVI"),
    ):
        s = _statements(tmp_path, dosages)["trait:TAS2R38"]
        assert s["kind"] == "finding"
        assert json.loads(str(s["claim"]))["result"] == diplotype
        assert str(s["reason"]).startswith(diplotype)
    for dosages in (rare, reference):
        s = _statements(tmp_path, dosages)["trait:TAS2R38"]
        assert (s["kind"], s["reason"]) == ("unassessed", RARE_TAS2R38)
    unknown = _statements(tmp_path, {**het, TAS2R38[2]: "unknown"})["trait:TAS2R38"]
    assert unknown["kind"] == "unassessed" and "rs10246939" in str(unknown["reason"])


def test_single_site_traits(tmp_path: Path) -> None:
    """ADH1B, ACTN3 and CYP1A2 read one site, by effect-allele copies: ADH1B's effect
    allele is the REF, so a reference call is two fast copies (His/His); an unknown
    dosage is unassessed."""
    ref = _statements(tmp_path, {})
    adh1b = json.loads(str(ref["trait:ADH1B"]["claim"]))
    assert (ref["trait:ADH1B"]["kind"], adh1b["result"]) == ("finding", 2)
    assert (
        json.loads(str(ref["trait:ACTN3"]["claim"]))["result"],
        ref["trait:ACTN3"]["kind"],
    ) == (
        0,
        "finding",
    )
    carried = _statements(tmp_path, {ADH1B: "2", ACTN3: "2", CYP1A2: "1"})
    assert json.loads(str(carried["trait:ADH1B"]["claim"]))["result"] == 0
    assert str(carried["trait:ACTN3"]["reason"]).startswith("XX")
    assert json.loads(str(carried["trait:CYP1A2"]["claim"]))["result"] == 1
    assert ref["trait:ADH1B"]["reason"] != carried["trait:ADH1B"]["reason"]
    unknown = _statements(
        tmp_path, {ADH1B: "unknown", ACTN3: "unknown", CYP1A2: "unknown"}
    )
    for sid in ("trait:ADH1B", "trait:ACTN3", "trait:CYP1A2"):
        assert unknown[sid]["kind"] == "unassessed"
        assert json.loads(str(unknown[sid]["claim"]))["result"] is None


def test_traits_are_information(tmp_path: Path) -> None:
    """Every trait statement is information: a finding (or unassessed) under traits at
    validity moderate, its reasons the trait's wording and `INFORMATION_NOTE`."""
    collected = TraitSource(graded=_graded({})).collect(run_context(tmp_path))
    for row in collected.rows:
        if row.statement["kind"] == "finding":
            assert INFORMATION_NOTE in str(row.statement["reason"])
    assert collected.extra == {}


def test_stored_traits_pass_check_calls(tmp_path: Path) -> None:
    """Ruling R24: gene-level trait statements on reference calls (all-REF; TAS2R38
    (2, 2, 1), a rare haplotype) pass the run's call invariant; the claim holds each site's
    genotype and level, and the lowest over the sites the result reads."""
    low = Reliability(p_wrong=0.1, level=Level.LOW, flags=(), warnings=())
    for dosages in ({}, {TAS2R38[0]: "2", TAS2R38[1]: "2", TAS2R38[2]: "1"}):
        graded = _graded(dosages)
        graded[LACTASE["rs145946881"]] = (graded[LACTASE["rs145946881"]][0], low)
        con = connect(tmp_path / f"s{len(dosages)}.duckdb")
        run_id = start_run(con, {}, {})
        collected = TraitSource(graded=graded).collect(
            replace(run_context(tmp_path), run_id=run_id)
        )
        with transaction(con):
            writer = StoreWriter(con, tmp_path / f"w{len(dosages)}")
            write_collected(writer, con, run_id, [collected])
            pipeline.check_calls(con, run_id)
        lactase = json.loads(str(collected.rows[0].statement["claim"]))
        assert (lactase["result"], lactase["reliability"]) == ("unlikely", "low")
        sites = lactase["sites"]
        assert (sites[0]["genotype"], sites[4]["level"]) == ("0/0", "low")
        assert query(con, run_id, Filters()).hidden["reliability_below:medium"] == 1
