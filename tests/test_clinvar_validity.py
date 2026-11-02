from typing import Any

import pytest

from genokeeper.config import load_grading
from genokeeper.evidence.clingen import Curation, gene_level
from genokeeper.evidence.clinvar import (
    ASSERTION_TERMS,
    BA1_REASON,
    BASE_CLASSES,
    BS1_REASON,
    CLNSIG_MAP,
    CONFLICT_REASON,
    GENOTYPE_REASON,
    GERMLINE_TERMS,
    HEMIZYGOUS_REASON,
    LEGACY_RISK_TERM,
    NOT_MATCHED_REASON,
    ValidityGrade,
    clnsig_class,
    grade_clinvar,
    parse_clnsig,
    stars,
)
from genokeeper.evidence.disease import match_curation
from genokeeper.evidence.frequency import load_ba1_exceptions
from genokeeper.genelists import load_gene_lists
from genokeeper.store.db import Kind, Topic, Validity

GRADING = load_grading()
LISTS = load_gene_lists()

TWO_STARS = "criteria_provided,_multiple_submitters,_no_conflicts"
ONE_STAR = "criteria_provided,_single_submitter"
CONFLICT_STATUS = "criteria_provided,_conflicting_classifications"
NO_STARS = "no_assertion_criteria_provided"
EXPERT_PANEL = "reviewed_by_expert_panel"
# `_grade`'s Definitive curation is gene-level: matched to none of the claim's diseases.
NOT_MATCHED = NOT_MATCHED_REASON


def _grade(**overrides: Any) -> ValidityGrade:
    """Pathogenic, 2 stars, in BRCA1 (SF, AD, Definitive), het, no frequency; then overrides."""
    kwargs: dict[str, Any] = {
        "clnsig": "Pathogenic",
        "clnrevstat": TWO_STARS,
        "gene": "BRCA1",
        "gene_validity": "Definitive",
        "lists": LISTS,
        "faf_grpmax": None,
        "dosage": "1",
        "grading": GRADING,
    }
    kwargs.update(overrides)
    return grade_clinvar(**kwargs)


@pytest.mark.parametrize(
    ("clnrevstat", "expected"),
    [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        (TWO_STARS, 2),
        (ONE_STAR, 1),
        (CONFLICT_STATUS, 1),
        ("criteria_provided,_conflicting_interpretations", 1),
        (NO_STARS, 0),
        ("no_classification_provided", 0),
        ("no_classifications_from_unflagged_records", 0),
        ("", 0),
    ],
)
def test_stars(clnrevstat: str, expected: int) -> None:
    assert stars(clnrevstat) == expected


@pytest.mark.parametrize(
    ("case", "overrides", "level", "kind", "topic", "reason"),
    [
        (
            "P 2 stars, Definitive gene",
            {},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "low penetrance is one term",
            {"clnsig": "Pathogenic,_low_penetrance"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "LP 1 star",
            {"clnsig": "Likely_pathogenic", "clnrevstat": ONE_STAR},
            Validity.MODERATE,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "Disputed gene",
            {"gene_validity": "Disputed"},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            f"gene validity Disputed; {NOT_MATCHED}",
        ),
        (
            "Refuted gene, any case",
            {"gene_validity": "REFUTED"},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            f"gene validity REFUTED; {NOT_MATCHED}",
        ),
        (
            "conflicting review status",
            {"clnrevstat": CONFLICT_STATUS},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "conflicting classifications",
        ),
        (
            "conflicting class",
            {
                "clnsig": "Conflicting_classifications_of_pathogenicity",
                "clnrevstat": CONFLICT_STATUS,
            },
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "conflicting classifications",
        ),
        (
            "frequency above BA1",
            {"faf_grpmax": 0.2},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "population frequency above BA1",
        ),
        (
            "frequency above the gene's BS1 (BRCA1 0.0001)",
            {"faf_grpmax": 0.0005},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "population frequency above BS1",
        ),
        (
            "same frequency under the fallback BS1 (TTR has no per-gene value)",
            {"faf_grpmax": 0.0005, "gene": "TTR"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "frequency at BS1 is not above it",
            {"faf_grpmax": 0.0001},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "4 stars stand on their own over a Limited gene",
            {"clnrevstat": "practice_guideline", "gene_validity": "Limited"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            None,
        ),
        (
            "no curation",
            {"gene_validity": None},
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            "no ClinGen gene validity",
        ),
        (
            "an uncurated word",
            {"gene_validity": "No Known Disease Relationship"},
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"gene validity No Known Disease Relationship; {NOT_MATCHED}",
        ),
        (
            "frequency conflict beats a Limited gene",
            {"gene_validity": "Limited", "faf_grpmax": 0.2},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "population frequency above BA1",
        ),
        (
            "0 stars",
            {"clnrevstat": NO_STARS},
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"no assertion criteria; {NOT_MATCHED}",
        ),
        (
            "VUS in an SF gene",
            {"clnsig": "Uncertain_significance"},
            None,
            Kind.WATCH,
            None,
            None,
        ),
        (
            "VUS in a Tier 3 gene",
            {"clnsig": "Uncertain_significance", "gene": "CFTR"},
            None,
            Kind.WATCH,
            None,
            None,
        ),
        (
            "VUS outside the lists",
            {"clnsig": "Uncertain_significance", "gene": "F5"},
            None,
            Kind.NON_FINDING,
            None,
            None,
        ),
        ("Benign", {"clnsig": "Benign"}, None, Kind.NON_FINDING, None, None),
        (
            "Benign/Likely_benign",
            {"clnsig": "Benign/Likely_benign"},
            None,
            Kind.NON_FINDING,
            None,
            None,
        ),
        (
            "dosage 0",
            {"dosage": "0"},
            Validity.STRONG,
            Kind.NON_FINDING,
            None,
            NOT_MATCHED,
        ),
        (
            "het P in an AR Tier 3 gene",
            {"gene": "CFTR"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.CARRIER,
            NOT_MATCHED,
        ),
        (
            "hom P in an AR Tier 3 gene",
            {"gene": "CFTR", "dosage": "2"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "het P in an AR ACMG SF gene in carrier Tier 4 (MUTYH)",
            {"gene": "MUTYH"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.CARRIER,
            NOT_MATCHED,
        ),
        (
            "hom P in an AR gene outside Tier 3",
            {"gene": "MUTYH", "dosage": "2"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "het P in an XL Tier 3 gene, female",
            {"gene": "DMD"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.CARRIER,
            NOT_MATCHED,
        ),
        (
            "hemizygous P in an XL Tier 3 gene, male",
            {"gene": "DMD", "male": True},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "unlisted gene, inheritance AR, het",
            {"gene": "F5", "inheritance": "AR"},
            Validity.STRONG,
            Kind.NON_FINDING,
            None,
            NOT_MATCHED,
        ),
        (
            "unlisted gene, inheritance AD",
            {"gene": "F5", "inheritance": "AD"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "lists win over inheritance",
            {"gene": "CFTR", "inheritance": "AD"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.CARRIER,
            NOT_MATCHED,
        ),
        (
            "unlisted gene, inheritance unknown",
            {"gene": "F5"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            f"{NOT_MATCHED}; inheritance unknown",
        ),
        (
            "frequency reason stays first",
            {"gene": "F5", "faf_grpmax": 0.03},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "population frequency above BS1; inheritance unknown",
        ),
        (
            "risk allele in an AR gene",
            {"clnsig": "risk_factor", "gene": "MUTYH"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "established risk allele, frequent, not conflicting",
            {"clnsig": "Established_risk_allele", "gene": "F5", "faf_grpmax": 0.03},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "P at an unknown genotype",
            {"dosage": "unknown"},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{NOT_MATCHED}; genotype unknown",
        ),
        (
            "P at an unknown genotype in a Tier 3 gene",
            {"dosage": "unknown", "gene": "CFTR", "clnrevstat": ONE_STAR},
            Validity.MODERATE,
            Kind.UNASSESSED,
            Topic.CARRIER,
            f"{NOT_MATCHED}; genotype unknown",
        ),
        (
            "P at an unknown genotype in an XL Tier 3 gene, male",
            {"dosage": "unknown", "gene": "DMD", "male": True},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{NOT_MATCHED}; genotype unknown",
        ),
        (
            "P at an unknown genotype in an AR gene outside the carrier lists",
            {"dosage": "unknown", "gene": "F5", "inheritance": "AR"},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{NOT_MATCHED}; genotype unknown",
        ),
        (
            "P at an unknown genotype with no inheritance known",
            {"dosage": "unknown", "gene": "NOVEL1", "gene_validity": "Definitive"},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{NOT_MATCHED}; inheritance unknown; genotype unknown",
        ),
        (
            "a risk allele at an unknown genotype",
            {"dosage": "unknown", "clnsig": "Established_risk_allele", "gene": "F5"},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{NOT_MATCHED}; genotype unknown",
        ),
        (
            "a conflicting claim at an unknown genotype",
            {"dosage": "unknown", "clnrevstat": CONFLICT_STATUS},
            Validity.CONFLICTING,
            Kind.UNASSESSED,
            Topic.RISKS,
            f"{CONFLICT_REASON}; genotype unknown",
        ),
        (
            "a listed VUS stays a watch at an unknown genotype",
            {"dosage": "unknown", "clnsig": "Uncertain_significance"},
            None,
            Kind.WATCH,
            None,
            None,
        ),
        (
            "benign at an unknown genotype stays a non-finding",
            {"dosage": "unknown", "clnsig": "Benign"},
            None,
            Kind.NON_FINDING,
            None,
            None,
        ),
        (
            "drug response: stars alone, under drugs",
            {"clnsig": "drug_response"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.DRUGS,
            None,
        ),
        (
            "P and B terms of one classification conflict",
            {"clnsig": "Pathogenic/Benign"},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            "conflicting classifications",
        ),
        (
            "P with a risk-factor part stays pathogenic",
            {"clnsig": "Pathogenic/Likely_pathogenic|risk_factor"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "P keeps its own stars after a risk-factor part",
            {"clnsig": "Pathogenic|risk_factor", "clnrevstat": ONE_STAR},
            Validity.MODERATE,
            Kind.FINDING,
            Topic.RISKS,
            NOT_MATCHED,
        ),
        (
            "a benign germline part hides a later risk factor",
            {"clnsig": "Benign|risk_factor"},
            None,
            Kind.NON_FINDING,
            None,
            None,
        ),
        (
            "VUS then risk factor: a risk allele capped at limited",
            {"clnsig": "Uncertain_significance|risk_factor"},
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"risk_factor has no review status; {NOT_MATCHED}",
        ),
        (
            "VUS then risk factor in an AR gene, any star count",
            {
                "clnsig": "Uncertain_significance|risk_factor",
                "gene": "MUTYH",
                "clnrevstat": "practice_guideline",
            },
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"risk_factor has no review status; {NOT_MATCHED}",
        ),
        (
            "the germline review status does not conflict the risk part",
            {
                "clnsig": "Uncertain_significance|risk_factor",
                "clnrevstat": CONFLICT_STATUS,
            },
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"risk_factor has no review status; {NOT_MATCHED}",
        ),
        (
            "VUS then risk factor in an uncurated gene keeps both reasons",
            {
                "clnsig": "Uncertain_significance|risk_factor",
                "gene": "F5",
                "gene_validity": None,
            },
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            "risk_factor has no review status; no ClinGen gene validity",
        ),
        (
            "VUS then risk factor in a Disputed gene conflicts",
            {
                "clnsig": "Uncertain_significance|risk_factor",
                "gene_validity": "Disputed",
            },
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            f"risk_factor has no review status; gene validity Disputed; {NOT_MATCHED}",
        ),
        (
            "empty germline part, then risk factor",
            {"clnsig": "|risk_factor", "gene": "F5"},
            Validity.LIMITED,
            Kind.FINDING,
            Topic.RISKS,
            f"risk_factor has no review status; {NOT_MATCHED}",
        ),
        (
            "uncertain risk allele is a VUS: watch in a listed gene",
            {"clnsig": "Uncertain_risk_allele"},
            None,
            Kind.WATCH,
            None,
            None,
        ),
        (
            "uncertain risk allele outside the lists",
            {"clnsig": "Uncertain_risk_allele", "gene": "F5"},
            None,
            Kind.NON_FINDING,
            None,
            None,
        ),
    ],
)
def test_grade_clinvar(
    case: str,
    overrides: dict[str, Any],
    level: Validity | None,
    kind: Kind,
    topic: Topic | None,
    reason: str | None,
) -> None:
    grade = _grade(**overrides)
    assert (grade.level, grade.kind, grade.topic, grade.reason) == (
        level,
        kind,
        topic,
        reason,
    ), case
    assert grade.alert is None, case


def test_tier4_clinvar_carrier() -> None:
    """A heterozygous ClinVar P in a carrier Tier 4 gene outside Tier 3 (ABCA4) is a
    carrier finding, as in Tier 3; two copies are under risks. Tier 4 enters the carrier
    gate only: a VUS there is no watch item. SYNE1, Tier 4 but AD in CGD, stays a dominant
    finding (`gene_mode` reads CGD before Tier 4)."""
    assert "ABCA4" not in LISTS["carrier_tier3"].genes
    het = _grade(gene="ABCA4")
    assert (het.kind, het.topic) == (Kind.FINDING, Topic.CARRIER)
    hom = _grade(gene="ABCA4", dosage="2")
    assert (hom.kind, hom.topic) == (Kind.FINDING, Topic.RISKS)
    vus = _grade(gene="ABCA4", clnsig="Uncertain_significance")
    assert (vus.kind, vus.topic) == (Kind.NON_FINDING, None)
    dominant = _grade(gene="SYNE1", inheritance="AD")
    assert (dominant.kind, dominant.topic) == (Kind.FINDING, Topic.RISKS)


def test_unmapped_value_is_a_limited_finding_with_an_alert() -> None:
    grade = _grade(clnsig="Weird_new_term", clnrevstat="practice_guideline")
    assert grade.level == Validity.LIMITED
    assert grade.alert == "unmapped CLNSIG: Weird_new_term"
    assert (grade.kind, grade.topic, grade.reason) == (Kind.FINDING, Topic.RISKS, None)


def test_unmapped_value_at_an_unknown_genotype_is_unassessed() -> None:
    grade = _grade(clnsig="Weird_new_term", dosage="unknown")
    assert (grade.kind, grade.topic) == (Kind.UNASSESSED, Topic.RISKS)
    assert grade.alert == "unmapped CLNSIG: Weird_new_term"
    assert grade.reason == "genotype unknown"


def test_unmapped_value_at_dosage_0_is_not_a_finding() -> None:
    grade = _grade(clnsig="Weird_new_term", dosage="0")
    assert grade.kind == Kind.NON_FINDING
    assert grade.alert == "unmapped CLNSIG: Weird_new_term"


def test_unmapped_part_keeps_the_mapped_class() -> None:
    grade = _grade(clnsig="Pathogenic|Weird_new_term")
    assert grade.level == Validity.LIMITED
    assert grade.alert == "unmapped CLNSIG: Pathogenic|Weird_new_term"
    assert (grade.kind, grade.topic) == (Kind.FINDING, Topic.RISKS)


def test_source_level_is_the_raw_value_and_stars() -> None:
    grade = _grade(clnsig="Pathogenic,_low_penetrance", clnrevstat=ONE_STAR)
    assert grade.source_level == "Pathogenic,_low_penetrance (1★)"


def test_clnsig_is_never_split_on_commas() -> None:
    assert clnsig_class("Pathogenic,_low_penetrance") == "pathogenic"
    assert clnsig_class("Likely_pathogenic,_low_penetrance") == "pathogenic"
    assert clnsig_class("Pathogenic/Likely_pathogenic") == "pathogenic"
    assert clnsig_class("Pathogenic/Pathogenic,_low_penetrance") == "pathogenic"
    assert clnsig_class("Uncertain_significance|drug_response") == "vus"
    assert clnsig_class("Uncertain_significance|risk_factor") == "risk"
    assert clnsig_class("Benign|risk_factor") == "benign"
    assert clnsig_class("Uncertain_risk_allele") == "vus"
    # "Pathogenic" and "_low_penetrance" are not looked up as parts.
    assert clnsig_class("Pathogenic,_weird") is None


@pytest.mark.parametrize(
    "value",
    [
        "Pathogenic,_low_penetrance",
        "Established_risk_allele",
        "risk_factor",
        "drug_response",
        "confers_sensitivity",
        "VUS-high",
        "VUS-mid",
        "VUS-low",
        "no_classifications_from_unflagged_records",
        "protective",
        "association",
        "Affects",
    ],
)
def test_tallied_values_are_mapped(value: str) -> None:
    assert value in CLNSIG_MAP


def test_every_key_maps_to_a_base_class() -> None:
    assert BASE_CLASSES == {
        "pathogenic",
        "benign",
        "vus",
        "risk",
        "drug_response",
        "other",
        "conflicting",
    }
    assert set(CLNSIG_MAP.values()) <= BASE_CLASSES
    assert all("/" not in key and "|" not in key for key in CLNSIG_MAP)


def test_refuses_an_unknown_dosage_value() -> None:
    with pytest.raises(ValueError, match="dosage"):
        _grade(dosage="3")


@pytest.mark.parametrize(
    ("case", "overrides", "level", "kind", "topic", "reason"),
    [
        (
            "one copy in an AR carrier gene (MUTYH)",
            {"gene": "MUTYH"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "two copies written on a haploid locus",
            {"gene": "MUTYH", "dosage": "2"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "an AR Tier 3 gene is never carrier status",
            {"gene": "GAA"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "an AD gene",
            {"gene": "BRCA1"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "an XL Tier 3 gene",
            {"gene": "DMD"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "an unlisted gene of unknown inheritance",
            {"gene": "F5"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "an unlisted gene the caller calls AR",
            {"gene": "F5", "inheritance": "AR"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            HEMIZYGOUS_REASON,
        ),
        (
            "a conflicting claim in an AR gene",
            {"gene": "MUTYH", "clnrevstat": CONFLICT_STATUS},
            Validity.CONFLICTING,
            Kind.FINDING,
            Topic.RISKS,
            f"{CONFLICT_REASON}; {HEMIZYGOUS_REASON}",
        ),
        (
            "an unknown genotype in an AR Tier 3 gene is risks, not carrier",
            {"gene": "GAA", "dosage": "unknown"},
            Validity.STRONG,
            Kind.UNASSESSED,
            Topic.RISKS,
            GENOTYPE_REASON,
        ),
        (
            "no copy",
            {"gene": "MUTYH", "dosage": "0"},
            Validity.STRONG,
            Kind.NON_FINDING,
            None,
            None,
        ),
        (
            "a VUS stays a watch",
            {"gene": "MUTYH", "clnsig": "Uncertain_significance"},
            None,
            Kind.WATCH,
            None,
            None,
        ),
        (
            "a risk allele keeps its own grading",
            {"gene": "MUTYH", "clnsig": "Established_risk_allele"},
            Validity.STRONG,
            Kind.FINDING,
            Topic.RISKS,
            None,
        ),
    ],
)
def test_a_haploid_locus_is_hemizygous_whatever_the_gene_label(
    case: str,
    overrides: dict[str, Any],
    level: Validity | None,
    kind: Kind,
    topic: Topic | None,
    reason: str | None,
) -> None:
    # A disease-matched curation: the reasons are the locus's alone.
    grade = _grade(male=True, haploid=True, validity_matched=True, **overrides)
    assert (grade.level, grade.kind, grade.topic, grade.reason) == (
        level,
        kind,
        topic,
        reason,
    ), case


@pytest.mark.parametrize(
    ("male", "gene", "kind", "topic"),
    [
        (True, "F5", Kind.NON_FINDING, None),
        (False, "F5", Kind.NON_FINDING, None),
        (True, "GAA", Kind.FINDING, Topic.CARRIER),
        (False, "GAA", Kind.FINDING, Topic.CARRIER),
    ],
)
def test_a_diploid_locus_keeps_the_genes_label(
    male: bool, gene: str, kind: Kind, topic: Topic | None
) -> None:
    # A female's X, or a male's PAR: one copy of an AR allele is a het.
    grade = _grade(
        gene=gene, inheritance="AR", male=male, haploid=False, validity_matched=True
    )
    assert (grade.kind, grade.topic, grade.reason) == (kind, topic, None)


def _unmapped(
    value: str,
) -> tuple[str, dict[str, Any], Kind, Validity, Topic, None, str]:
    """A table row for a value the whitelist refuses: a limited finding with its alert."""
    return (
        value,
        {},
        Kind.FINDING,
        Validity.LIMITED,
        Topic.RISKS,
        None,
        f"unmapped CLNSIG: {value}",
    )


RISK_PART = "risk_factor has no review status"


@pytest.mark.parametrize(
    ("clnsig", "overrides", "kind", "level", "topic", "reason", "alert"),
    [
        _unmapped(""),
        _unmapped("|"),
        _unmapped("Weird"),
        _unmapped("Uncertain_significance|Weird"),
        _unmapped("Pathogenic|Weird"),
        _unmapped("Uncertain_significance/risk_factor"),
        _unmapped("Benign/Uncertain_significance|risk_factor"),
        _unmapped("protective|risk_factor"),
        _unmapped("not_provided|risk_factor"),
        _unmapped("drug_response|risk_factor"),
        _unmapped("Pathogenic/risk_factor"),
        _unmapped("Pathogenic/Established_risk_allele"),
        _unmapped("Uncertain_significance/Likely_benign"),
        _unmapped("Pathogenic/Benign/Uncertain_significance"),
        _unmapped("|drug_response"),
        _unmapped("Pathogenic|"),
        _unmapped("Pathogenic/"),
        _unmapped("Pathogenic,_weird"),
        # A germline classification term after `|` is never ignored: the value is unmapped.
        _unmapped("Pathogenic|Benign"),
        _unmapped("Benign|Pathogenic"),
        _unmapped("Likely_benign|Likely_pathogenic"),
        _unmapped("Benign|Conflicting_classifications_of_pathogenicity"),
        _unmapped("Uncertain_significance|Pathogenic"),
        _unmapped("Uncertain_significance|Pathogenic|risk_factor"),
        _unmapped("Benign|Established_risk_allele"),
        _unmapped("Likely_benign|Established_risk_allele"),
        _unmapped("|Established_risk_allele"),
        _unmapped("Uncertain_significance|Uncertain_risk_allele"),
        _unmapped("Pathogenic|VUS-high"),
        (
            "Uncertain_risk_allele|risk_factor",
            {},
            Kind.FINDING,
            Validity.LIMITED,
            Topic.RISKS,
            RISK_PART,
            None,
        ),
        (
            "Weird",
            {"dosage": "0"},
            Kind.NON_FINDING,
            Validity.LIMITED,
            None,
            None,
            "unmapped CLNSIG: Weird",
        ),
        (
            "Weird",
            {"dosage": "unknown"},
            Kind.UNASSESSED,
            Validity.LIMITED,
            Topic.RISKS,
            "genotype unknown",
            "unmapped CLNSIG: Weird",
        ),
        (
            "Pathogenic/Benign",
            {},
            Kind.FINDING,
            Validity.CONFLICTING,
            Topic.RISKS,
            "conflicting classifications",
            None,
        ),
        (
            "Pathogenic/Likely_pathogenic",
            {},
            Kind.FINDING,
            Validity.STRONG,
            Topic.RISKS,
            None,
            None,
        ),
        (
            "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance",
            {},
            Kind.FINDING,
            Validity.STRONG,
            Topic.RISKS,
            None,
            None,
        ),
        ("Benign/Likely_benign", {}, Kind.NON_FINDING, None, None, None, None),
        (
            "Pathogenic|risk_factor",
            {},
            Kind.FINDING,
            Validity.STRONG,
            Topic.RISKS,
            None,
            None,
        ),
        (
            "Conflicting_classifications_of_pathogenicity|risk_factor",
            {},
            Kind.FINDING,
            Validity.CONFLICTING,
            Topic.RISKS,
            "conflicting classifications",
            None,
        ),
        (
            "Uncertain_significance|risk_factor",
            {},
            Kind.FINDING,
            Validity.LIMITED,
            Topic.RISKS,
            RISK_PART,
            None,
        ),
        (
            "Uncertain_significance|drug_response|risk_factor",
            {},
            Kind.FINDING,
            Validity.LIMITED,
            Topic.RISKS,
            RISK_PART,
            None,
        ),
        (
            "|risk_factor",
            {},
            Kind.FINDING,
            Validity.LIMITED,
            Topic.RISKS,
            RISK_PART,
            None,
        ),
        ("Benign|risk_factor", {}, Kind.NON_FINDING, None, None, None, None),
        (
            "Uncertain_significance|drug_response",
            {},
            Kind.WATCH,
            None,
            None,
            None,
            None,
        ),
        (
            "risk_factor",
            {"gene": "F5"},
            Kind.FINDING,
            Validity.STRONG,
            Topic.RISKS,
            None,
            None,
        ),
        (
            "Established_risk_allele/Likely_risk_allele",
            {"gene": "F5"},
            Kind.FINDING,
            Validity.STRONG,
            Topic.RISKS,
            None,
            None,
        ),
        ("drug_response", {}, Kind.FINDING, Validity.STRONG, Topic.DRUGS, None, None),
        ("not_provided", {}, Kind.NON_FINDING, None, None, None, None),
    ],
)
def test_clnsig_values(
    clnsig: str,
    overrides: dict[str, Any],
    kind: Kind,
    level: Validity | None,
    topic: Topic | None,
    reason: str | None,
    alert: str | None,
) -> None:
    # A disease-matched curation: the reasons are the CLNSIG value's alone.
    grade = _grade(clnsig=clnsig, validity_matched=True, **overrides)
    assert (grade.kind, grade.level, grade.topic, grade.reason, grade.alert) == (
        kind,
        level,
        topic,
        reason,
        alert,
    )


def test_only_a_legacy_risk_factor_after_a_bar_makes_a_risk_part() -> None:
    after = sorted(ASSERTION_TERMS | GERMLINE_TERMS | set(CLNSIG_MAP))
    for germline in ["", *sorted(GERMLINE_TERMS)]:
        for term in after:
            parsed = parse_clnsig(f"{germline}|{term}")
            if term not in ASSERTION_TERMS:
                assert parsed.cls is None, (germline, term)
            if parsed.risk_part:
                assert term == LEGACY_RISK_TERM, (germline, term)
                grade = _grade(clnsig=f"{germline}|{term}")
                assert grade.reason is not None
                assert grade.reason.startswith(RISK_PART)


def test_assertion_terms_hold_no_germline_term() -> None:
    assert ASSERTION_TERMS.isdisjoint(GERMLINE_TERMS)
    assert {CLNSIG_MAP[t] for t in ASSERTION_TERMS} == {
        "risk",
        "drug_response",
        "other",
    }
    assert [t for t in ASSERTION_TERMS if CLNSIG_MAP[t] == "risk"] == [LEGACY_RISK_TERM]


# ClinVar's GRCh38 records of variation IDs 9 (HFE C282Y) and 1900 (BTD D444H).
HFE_C282Y = "GRCh38:chr6:26092913:G:A"
BTD_D444H = "GRCh38:chr3:15645186:G:C"
EXCEPTIONS = load_ba1_exceptions()


def test_exception_list_holds_hfe_c282y_and_btd_d444h() -> None:
    assert {HFE_C282Y, BTD_D444H} <= EXCEPTIONS
    assert all(vid.startswith("GRCh38:chr") for vid in EXCEPTIONS)


def test_recessive_gene_outside_table_has_no_bs1() -> None:
    # CFTR: AR (ACMG Tier 3), no [bs1_by_gene] value.
    grade = _grade(gene="CFTR", faf_grpmax=0.02)
    assert grade.level == Validity.STRONG
    assert BS1_REASON not in (grade.reason or "")
    # The caller's inheritance counts too, for a gene on neither list.
    grade = _grade(gene="F5", inheritance="AR", faf_grpmax=0.02)
    assert grade.level == Validity.STRONG
    assert BS1_REASON not in (grade.reason or "")


def test_recessive_gene_outside_table_still_has_ba1() -> None:
    grade = _grade(gene="CFTR", faf_grpmax=0.06)
    assert (grade.level, grade.reason) == (Validity.CONFLICTING, BA1_REASON)


def test_dominant_gene_outside_table_uses_fallback() -> None:
    grade = _grade(gene="F5", inheritance="AD", faf_grpmax=0.02)
    assert (grade.level, grade.reason) == (Validity.CONFLICTING, BS1_REASON)


def test_xlinked_gene_uses_fallback() -> None:
    # DMD: XL (ACMG Tier 3), no [bs1_by_gene] value.
    grade = _grade(gene="DMD", faf_grpmax=0.02)
    assert (grade.level, grade.reason) == (Validity.CONFLICTING, BS1_REASON)


def test_unknown_inheritance_uses_fallback() -> None:
    grade = _grade(gene="F5", faf_grpmax=0.02)
    assert grade.level == Validity.CONFLICTING
    assert grade.reason == f"{BS1_REASON}; inheritance unknown"


def test_table_gene_uses_its_own_bs1() -> None:
    # GAA is AR, and its own BS1 (0.005) applies all the same.
    above = _grade(gene="GAA", faf_grpmax=0.006)
    assert (above.level, above.reason) == (Validity.CONFLICTING, BS1_REASON)
    below = _grade(gene="GAA", faf_grpmax=0.004)
    assert below.level == Validity.STRONG


def test_exception_exempt_from_ba1() -> None:
    kwargs: dict[str, Any] = {"gene": "HFE", "inheritance": "AR", "faf_grpmax": 0.06}
    grade = _grade(vid=HFE_C282Y, ba1_exceptions=EXCEPTIONS, **kwargs)
    assert grade.level == Validity.STRONG
    # Another allele at the same frequency is not exempt.
    other = _grade(vid="GRCh38:chr6:26092913:G:T", ba1_exceptions=EXCEPTIONS, **kwargs)
    assert (other.level, other.reason) == (Validity.CONFLICTING, BA1_REASON)


def test_exception_exempt_from_bs1_in_table_gene() -> None:
    # BTD's own BS1 is 0.005726.
    kwargs: dict[str, Any] = {"gene": "BTD", "faf_grpmax": 0.04}
    grade = _grade(vid=BTD_D444H, ba1_exceptions=EXCEPTIONS, **kwargs)
    assert grade.level == Validity.STRONG
    other = _grade(vid="GRCh38:chr3:15645186:G:T", ba1_exceptions=EXCEPTIONS, **kwargs)
    assert (other.level, other.reason) == (Validity.CONFLICTING, BS1_REASON)


def test_three_stars_strong_without_curation() -> None:
    for status in (EXPERT_PANEL, "practice_guideline"):
        grade = _grade(
            clnrevstat=status, gene="F5", inheritance="AD", gene_validity=None
        )
        assert (grade.level, grade.reason) == (Validity.STRONG, None)


def test_three_stars_unmatched_gene_negative_still_strong() -> None:
    grade = _grade(clnrevstat=EXPERT_PANEL, gene_validity="Disputed")
    assert grade.level == Validity.STRONG
    assert grade.reason is not None
    assert "gene validity Disputed" in grade.reason
    assert NOT_MATCHED_REASON in grade.reason


def test_three_stars_matched_negative_conflicting() -> None:
    grade = _grade(
        clnrevstat=EXPERT_PANEL, gene_validity="Refuted", validity_matched=True
    )
    assert (grade.level, grade.reason) == (
        Validity.CONFLICTING,
        "gene validity Refuted",
    )


def test_three_stars_conflict_and_frequency_still_conflicting() -> None:
    grade = _grade(clnrevstat=EXPERT_PANEL, faf_grpmax=0.2)
    assert (grade.level, grade.reason) == (Validity.CONFLICTING, BA1_REASON)
    grade = _grade(clnsig="Pathogenic/Benign", clnrevstat=EXPERT_PANEL)
    assert (grade.level, grade.reason) == (Validity.CONFLICTING, CONFLICT_REASON)


def test_two_stars_unmatched_uses_gene_level_with_reason() -> None:
    grade = _grade(gene_validity="Definitive")
    assert (grade.level, grade.reason) == (Validity.STRONG, NOT_MATCHED_REASON)
    matched = _grade(gene_validity="Definitive", validity_matched=True)
    assert (matched.level, matched.reason) == (Validity.STRONG, None)
    negative = _grade(gene_validity="Disputed")
    assert (negative.level, negative.reason) == (
        Validity.CONFLICTING,
        f"gene validity Disputed; {NOT_MATCHED_REASON}",
    )
    limited = _grade(gene_validity="Limited", validity_matched=True)
    assert (limited.level, limited.reason) == (
        Validity.LIMITED,
        "gene validity Limited",
    )


def test_two_stars_unmatched_limited_and_refuted_gene_is_refuted() -> None:
    found = [Curation("MONDO:1", "Limited"), Curation("MONDO:2", "Refuted")]
    match = match_curation("BRCA1", [frozenset()], {"BRCA1": found}, gene_level(found))
    grade = _grade(gene_validity=match.classification, validity_matched=match.matched)
    assert (grade.level, grade.reason) == (
        Validity.CONFLICTING,
        f"gene validity Refuted; {NOT_MATCHED_REASON}",
    )


def test_validity_source_named_in_every_validity_reason() -> None:
    """A validity GenCC gave (slice 3f) says so in every reason naming it: the negative
    one, at 3 stars too, and the one that leaves a claim limited."""
    gencc = {"validity_source": "GenCC"}
    negative = _grade(gene_validity="Disputed", **gencc)
    assert negative.reason == f"gene validity Disputed (GenCC); {NOT_MATCHED_REASON}"
    expert = _grade(gene_validity="Refuted", clnrevstat=EXPERT_PANEL, **gencc)
    assert (expert.level, expert.reason) == (
        Validity.STRONG,
        f"gene validity Refuted (GenCC); {NOT_MATCHED_REASON}",
    )
    limited = _grade(gene_validity="Limited", validity_matched=True, **gencc)
    assert (limited.level, limited.reason) == (
        Validity.LIMITED,
        "gene validity Limited (GenCC)",
    )


@pytest.mark.parametrize(
    ("clnrevstat", "validity"),
    [
        ("practice_guideline", Validity.STRONG),
        (EXPERT_PANEL, Validity.STRONG),
        (TWO_STARS, Validity.STRONG),
        (ONE_STAR, Validity.MODERATE),
        (NO_STARS, Validity.LIMITED),
    ],
)
@pytest.mark.parametrize("dosage", ["1", "2"])
def test_drug_response_star_grading(
    clnrevstat: str, validity: Validity, dosage: str
) -> None:
    """A carried drug_response claim is a finding under drugs, its validity from the stars
    alone: no ClinGen gate (a Refuted gene, or none) and no frequency filter."""
    for gene_validity in (None, "Refuted", "Definitive"):
        grade = _grade(
            clnsig="drug_response",
            clnrevstat=clnrevstat,
            gene="CYP2C19",
            gene_validity=gene_validity,
            faf_grpmax=0.4,
            dosage=dosage,
        )
        assert (grade.kind, grade.level, grade.topic, grade.reason) == (
            Kind.FINDING,
            validity,
            Topic.DRUGS,
            None,
        )
        assert grade.source_level == f"drug_response ({stars(clnrevstat)}★)"


def test_drug_response_dosage_zero_non_finding() -> None:
    """Not carried: a non-finding, graded all the same; unknown: unassessed under drugs."""
    absent = _grade(clnsig="drug_response", gene="CYP2C19", dosage="0")
    assert (absent.kind, absent.level, absent.topic) == (
        Kind.NON_FINDING,
        Validity.STRONG,
        None,
    )
    unknown = _grade(clnsig="drug_response", gene="CYP2C19", dosage="unknown")
    assert (unknown.kind, unknown.level, unknown.topic, unknown.reason) == (
        Kind.UNASSESSED,
        Validity.STRONG,
        Topic.DRUGS,
        GENOTYPE_REASON,
    )
    # A drug_response after a germline classification is not this claim's class.
    assert _grade(clnsig="Benign|drug_response").kind == Kind.NON_FINDING
