"""Eleven trait statements (slice 3f's five, "Traits", and slice 3g's six single-gene
traits): one per trait, not per rsid, from our graded genotypes at `traits.tsv`'s alleles
(GRCh38, forward strand, dbSNP). Each new wording and caveat is written here as the closing
design's "B. Single-gene traits" table words it, never read back from the file."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from fixtures.run_context import run_context

from genokeeper import pipeline
from genokeeper.genotype.gvcf import NOCALL_REASON, Allele, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.query import Filters, query
from genokeeper.sources import INFORMATION_NOTE, write_collected
from genokeeper.sources.traits import (
    MC1R_WEAK_ONLY,
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
FUT2: Allele = ("chr19", 48703417, "G", "A")  # rs601338, effect A (W154X)
ABCC11: Allele = ("chr16", 48224287, "C", "T")  # rs17822931, effect T (dry earwax)
HERC2: Allele = ("chr15", 28120472, "A", "G")  # rs12913832, effect G (blue)
AMPD1: Allele = ("chr1", 114693436, "G", "A")  # rs17602729, effect A (Q12X)
# rs333, the 32-bp deletion left-aligned (dbSNP writes it right-shifted at 46373453).
CCR5: Allele = ("chr3", 46373452, "TACAGTCAGTATCAATTCTGGAAGAATTTCCAG", "T")
# MC1R's strong (R) sites, then its weak (r) ones; every effect allele is the ALT.
D84E: Allele = ("chr16", 89919510, "C", "A")  # rs1805006
R142H: Allele = ("chr16", 89919683, "G", "A")  # rs11547464
R151C: Allele = ("chr16", 89919709, "C", "T")  # rs1805007
I155T: Allele = ("chr16", 89919722, "T", "C")  # rs1110400
R160W: Allele = ("chr16", 89919736, "C", "T")  # rs1805008
D294H: Allele = ("chr16", 89920138, "G", "C")  # rs1805009
V60L: Allele = ("chr16", 89919436, "G", "T")  # rs1805005
V92M: Allele = ("chr16", 89919532, "G", "A")  # rs2228479
R163Q: Allele = ("chr16", 89919746, "G", "A")  # rs885479
MC1R_STRONG = (D84E, R142H, R151C, I155T, R160W, D294H)
MC1R_WEAK = (V60L, V92M, R163Q)
NEW_TRAITS = ("FUT2", "ABCC11", "MC1R", "HERC2", "AMPD1", "CCR5")
CAVEATS = {
    "FUT2": "the common European and African null allele only; the rule does not hold in"
    " East Asian ancestry",
    "ABCC11": "one SNP decides earwax type; odour is only a tendency",
    "MC1R": "hair colour depends on other genes too: in UK Biobank about 15% of people"
    " with two MC1R variants are red-haired; melanoma risk is in its polygenic score",
    "HERC2": "one SNP, for European ancestry; it cannot tell green from hazel",
    "AMPD1": "a common, mostly harmless variant; it predicts neither sports ability nor"
    " muscle disease",
    "CCR5": "not protection against HIV in general and no substitute for prevention",
}
MC1R_RED = (
    "red or strawberry-blond hair likely, with fair skin, freckling and easy sunburn"
)
MC1R_ONE = (
    "one red-hair variant: red hair is uncommon, fairer skin and sun sensitivity"
    " somewhat more likely"
)
MC1R_NONE = "no MC1R red-hair variant found"


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
    assert collected.input_rows == len(collected.rows) == 11
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
        *(f"trait:{name}" for name in NEW_TRAITS),
    }
    return out


def _shown(
    tmp_path: Path, dosages: dict[Allele, str], name: str
) -> tuple[object, str | None]:
    """The trait's claim result and its wording: the reason before the caveat, which must
    follow it with `INFORMATION_NOTE`; None for an unassessed statement."""
    s = _statements(tmp_path, dosages)[f"trait:{name}"]
    result = json.loads(str(s["claim"]))["result"]
    if s["kind"] == "unassessed":
        assert result is None
        return None, str(s["reason"])
    tail = f"; {CAVEATS[name]}; {INFORMATION_NOTE}"
    reason = str(s["reason"])
    assert reason.endswith(tail), reason
    return result, reason.removesuffix(tail)


def test_traits_file_holds_the_preflight_alleles() -> None:
    traits = {t.name: t for t in load_traits()}
    assert list(traits) == [
        "lactase",
        "ADH1B",
        "ACTN3",
        "TAS2R38",
        "CYP1A2",
        *NEW_TRAITS,
    ]
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
    for name, allele, effect in (
        ("FUT2", FUT2, "A"),
        ("ABCC11", ABCC11, "T"),
        ("HERC2", HERC2, "G"),
        ("AMPD1", AMPD1, "A"),
        ("CCR5", CCR5, "T"),
    ):
        trait = traits[name]
        assert (trait.gene, trait.rule) == (name, "single")
        assert (trait.alleles, trait.effects) == ((allele,), (effect,))
        assert trait.note == CAVEATS[name]
    mc1r = traits["MC1R"]
    assert (mc1r.gene, mc1r.rule, mc1r.note) == ("MC1R", "mc1r", CAVEATS["MC1R"])
    assert mc1r.alleles == (*MC1R_STRONG, *MC1R_WEAK)
    assert mc1r.rsids == (
        *("rs1805006", "rs11547464", "rs1805007", "rs1110400", "rs1805008"),
        *("rs1805009", "rs1805005", "rs2228479", "rs885479"),
    )
    assert mc1r.effects == tuple(a[3] for a in mc1r.alleles)
    assert mc1r.classes == ("strong",) * 6 + ("weak",) * 3
    assert all(c == "" for t in traits.values() if t.rule != "mc1r" for c in t.classes)
    assert CCR5 in TRAIT_ALLELES
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


def test_fut2_secretor(tmp_path: Path) -> None:
    non_secretor = (
        "likely non-secretor: no working FUT2 enzyme; non-secretors are rarely infected"
        " by the common norovirus GII.4 strains"
    )
    for dosage, shown in (
        ("2", (2, non_secretor)),
        ("1", (1, "likely secretor")),
        ("0", (0, "likely secretor")),
    ):
        assert _shown(tmp_path, {FUT2: dosage}, "FUT2") == shown
    assert _shown(tmp_path, {FUT2: "unknown"}, "FUT2") == (
        None,
        "genotype unknown at rs601338",
    )


def test_abcc11_earwax(tmp_path: Path) -> None:
    dry = "dry earwax likely; little typical underarm odour"
    for dosage, shown in (
        ("2", (2, dry)),
        ("1", (1, "wet earwax likely")),
        ("0", (0, "wet earwax likely")),
    ):
        assert _shown(tmp_path, {ABCC11: dosage}, "ABCC11") == shown
    assert _shown(tmp_path, {ABCC11: "unknown"}, "ABCC11") == (
        None,
        "genotype unknown at rs17822931",
    )


def test_herc2_eye_colour(tmp_path: Path) -> None:
    for dosage, shown in (
        ("2", (2, "blue or other light eyes likely (brown in about 3%)")),
        ("1", (1, "brown or intermediate (hazel, green) most likely")),
        ("0", (0, "brown eyes likely")),
    ):
        assert _shown(tmp_path, {HERC2: dosage}, "HERC2") == shown


def test_ampd1_three_genotypes(tmp_path: Path) -> None:
    no_enzyme = (
        "no working muscle AMPD1 enzyme: usually no symptoms, occasionally"
        " exercise-related cramps or tiredness"
    )
    for dosage, shown in (
        ("2", (2, no_enzyme)),
        ("1", (1, "one copy: no effect expected")),
        ("0", (0, "no AMPD1 Q12X variant")),
    ):
        assert _shown(tmp_path, {AMPD1: dosage}, "AMPD1") == shown


def test_ccr5_delta32_left_aligned(tmp_path: Path) -> None:
    """The deletion is left-aligned, as DeepVariant emits it (dbSNP writes it right-shifted
    at 46373453); its effect allele is the ALT, the deletion."""
    assert CCR5 in TRAIT_ALLELES
    resistant = (
        "two copies of CCR5Δ32: strongly resistant to HIV-1 strains that use CCR5,"
        " not to those that use CXCR4"
    )
    for dosage, shown in (
        ("2", (2, resistant)),
        ("1", (1, "one copy: may slow progression after HIV infection")),
        ("0", (0, "no CCR5Δ32")),
    ):
        assert _shown(tmp_path, {CCR5: dosage}, "CCR5") == shown


def test_mc1r_two_strong_alleles(tmp_path: Path) -> None:
    """R counts strong copies across sites (two heterozygous strong sites are two: strong
    loss-of-function alleles are not seen on one haplotype); weak copies never count
    toward R."""
    for dosages, shown in (
        ({R151C: "2"}, ("two strong", MC1R_RED)),
        ({R151C: "1", D294H: "1"}, ("two strong", MC1R_RED)),
        ({R151C: "1", V60L: "2"}, ("one strong", MC1R_ONE)),
        ({V60L: "1", V92M: "1"}, ("weak only", MC1R_WEAK_ONLY)),
        ({V60L: "2", R163Q: "2"}, ("weak only", MC1R_WEAK_ONLY)),
        ({}, ("none", MC1R_NONE)),
    ):
        assert _shown(tmp_path, dosages, "MC1R") == shown, dosages
    assert MC1R_WEAK_ONLY == "weak variants only: little effect expected"


def test_mc1r_reliability_reads_the_sites_its_result_reads(tmp_path: Path) -> None:
    """ "two strong" reads the carried strong sites, "one strong" every strong site,
    "weak only" every strong site and the carried weak ones, "none" every site."""
    low = Reliability(p_wrong=0.1, level=Level.LOW, flags=(), warnings=())
    for dosages, lowered, reliability in (
        ({R151C: "2"}, V60L, "high"),
        ({R151C: "2", V60L: "1"}, V60L, "high"),
        ({R151C: "2"}, D84E, "high"),
        ({R151C: "1"}, D84E, "low"),
        ({R151C: "1"}, V60L, "high"),
        ({V92M: "1"}, V60L, "high"),
        ({V92M: "1"}, D84E, "low"),
        ({}, V60L, "low"),
    ):
        graded = _graded(dosages)
        graded[lowered] = (graded[lowered][0], low)
        rows = TraitSource(graded=graded).collect(run_context(tmp_path)).rows
        mc1r = next(r.statement for r in rows if r.statement["source_id"] == "MC1R")
        assert json.loads(str(mc1r["claim"]))["reliability"] == reliability


def test_mc1r_unknown_strong_site(tmp_path: Path) -> None:
    """An unknown strong site leaves the result unassessed unless two strong copies are
    already known; an unknown weak site matters only when no strong and no weak copy is
    known."""
    assert _shown(tmp_path, {R151C: "1", R160W: "unknown"}, "MC1R") == (
        None,
        "genotype unknown at rs1805008",
    )
    assert _shown(tmp_path, {R151C: "1", D294H: "1", R160W: "unknown"}, "MC1R") == (
        "two strong",
        MC1R_RED,
    )
    assert _shown(tmp_path, {V60L: "unknown"}, "MC1R") == (
        None,
        "genotype unknown at rs1805005",
    )
    assert _shown(tmp_path, {V92M: "1", V60L: "unknown"}, "MC1R") == (
        "weak only",
        MC1R_WEAK_ONLY,
    )
    assert _shown(tmp_path, {V92M: "1", R160W: "unknown"}, "MC1R") == (
        None,
        "genotype unknown at rs1805008",
    )
    assert _shown(tmp_path, {R151C: "1", V60L: "unknown"}, "MC1R") == (
        "one strong",
        MC1R_ONE,
    )


def test_new_traits_are_information(tmp_path: Path) -> None:
    """Each of slice 3g's six statements: information under traits at validity moderate,
    never actionable, its reasons its wording, its caveat and `INFORMATION_NOTE`."""
    statements = _statements(tmp_path, {})
    for name in NEW_TRAITS:
        s = statements[f"trait:{name}"]
        assert (s["kind"], s["topic"], s["validity"], s["actionable"]) == (
            "finding",
            "traits",
            "moderate",
            False,
        )
        assert str(s["reason"]).endswith(f"; {CAVEATS[name]}; {INFORMATION_NOTE}")


def test_load_traits_refuses_a_misplaced_class(tmp_path: Path) -> None:
    """An `mc1r` site needs a class (strong or weak); any other rule's site has none."""
    header = (
        "trait\tgene\trule\trsid\tchrom\tpos\tref\talt\teffect\tclass\tnone\tone"
        "\ttwo\tnote\tsource\n"
    )
    for rule, cls in (("mc1r", ""), ("mc1r", "medium"), ("single", "strong")):
        path = tmp_path / f"{rule}-{cls}.tsv"
        path.write_text(
            header + f"X\tX\t{rule}\trs1\tchr1\t1\tA\tG\tG\t{cls}\ta\tb\tc\t\ts\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="class"):
            load_traits(path)
