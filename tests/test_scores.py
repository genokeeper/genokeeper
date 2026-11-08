"""The polygenic score list (`genokeeper.prs.scores`): one score per condition of the track
design's slice 1d, each with its PGS Catalog evaluations, licence and topic."""

import re

import pytest

from genokeeper.prs.scores import (
    MIXED,
    POPULATIONS,
    NotCovered,
    load_not_covered,
    load_research_scores,
    load_scores,
    parse_not_covered,
    parse_research_scores,
    parse_scores,
    topic_for,
)
from genokeeper.store.db import Topic

# Track design, slice 1d "Scores" (design decision 2026-10-19): heart and vessels,
# metabolism, cancer, neuro.
CONDITIONS = (
    "coronary artery disease",
    "atrial fibrillation",
    "hypertension",
    "type 2 diabetes",
    "LDL cholesterol",
    "body mass index",
    "prostate cancer",
    "colorectal cancer",
    "melanoma",
    "Alzheimer's disease",
    "Parkinson's disease",
    # Closing design, "Slice 3g" (design decision 2026-10-30): every candidate condition with a
    # score that passes the rule read strictly.
    "ischaemic stroke",
    "heart failure",
    "abdominal aortic aneurysm",
    "peripheral artery disease",
    "venous thromboembolism",
    "primary open-angle glaucoma",
    "chronic kidney disease",
    "gout",
    "bone mineral density (a low percentile means a higher osteoporosis risk)",
    "triglycerides",
    "HDL cholesterol",
    "age-related macular degeneration",
    "lung cancer",
    "pancreatic cancer",
    "thyroid cancer",
    "type 1 diabetes",
    "rheumatoid arthritis",
    "systemic lupus erythematosus",
    "multiple sclerosis",
    "Crohn's disease",
    "ulcerative colitis",
    "hypothyroidism",
    "asthma",
    "atopic dermatitis",
    "COPD",
    "gallstones",
    "diverticular disease",
    "chronic lymphocytic leukaemia",
    "multiple myeloma",
    "testicular cancer",
    "depression",
    "schizophrenia",
    "height",
    "vitamin D",
)

# Closing design, "Slice 3g", "Weak predictive power": tier 3, or a weak qualifying European
# evaluation.
WEAK = {
    "PGS004885",
    "PGS000134",
    "PGS002724",
    "PGS001790",
    "PGS002378",
    "PGS019980",
    "PGS019973",
    "PGS002755",
}

# Closing design, "Slice 3g", "Not covered", in its order.
NOT_COVERED = (
    ("intraocular pressure", Topic.WELLNESS),
    ("urate", Topic.WELLNESS),
    ("apolipoprotein B", Topic.WELLNESS),
    ("lipoprotein(a)", Topic.WELLNESS),
    ("basal cell carcinoma", Topic.RISKS),
    ("psoriasis", Topic.RISKS),
    ("kidney stones", Topic.RISKS),
    ("benign prostatic hyperplasia", Topic.RISKS),
    ("migraine", Topic.RISKS),
    ("MASLD", Topic.RISKS),
    ("cataract", Topic.RISKS),
    ("hearing loss", Topic.RISKS),
    ("glioma", Topic.RISKS),
    ("coeliac disease", Topic.RISKS),
    ("anxiety", Topic.RISKS),
    ("insomnia", Topic.RISKS),
    ("bipolar disorder", Topic.RISKS),
    ("chronotype", Topic.TRAITS),
    ("neuroticism", Topic.TRAITS),
    ("risk tolerance", Topic.TRAITS),
    ("sleep duration", Topic.TRAITS),
    ("educational attainment", Topic.TRAITS),
    ("hair colour", Topic.TRAITS),
    ("FEV1", Topic.WELLNESS),
)

_SCORE = (
    '[[score]]\npgs_id = "PGS000001"\ncondition = "c"\ntrait_efo = "EFO_1"\n'
    'publication = "p"\ndevelopment_ancestry = "d"\nlicence = "l"\n'
    'evaluations = [["PPM000001", "EUR"]]\n'
)


def _not_covered(condition: str, topic: str = "risks") -> str:
    return f'[[not_covered]]\ncondition = "{condition}"\ntopic = "{topic}"\nwhy = "w"\n'


def test_every_condition_has_one_score() -> None:
    scores = load_scores()
    assert sorted(s.condition for s in scores) == sorted(CONDITIONS)
    assert len({s.pgs_id for s in scores}) == len(scores)
    assert all(re.fullmatch(r"PGS\d{6}", s.pgs_id) for s in scores)


def test_evaluation_ids_recorded() -> None:
    for score in load_scores():
        assert score.evaluations, score.pgs_id
        for ppm, ancestry in score.evaluations:
            assert re.fullmatch(r"PPM\d{6}", ppm), (score.pgs_id, ppm)
            assert ancestry in (*POPULATIONS, MIXED), (score.pgs_id, ancestry)
        assert score.evaluation_ancestry == {a for _, a in score.evaluations}
    # PGS Catalog's own ancestry words are not pgsc_calc's population codes: refused.
    text = (
        '[[score]]\npgs_id = "PGS000001"\ncondition = "c"\ntrait_efo = "EFO_1"\n'
        'publication = "p"\ndevelopment_ancestry = "d"\nlicence = "l"\n'
        'evaluations = [["PPM000001", "European"]]\n'
    )
    with pytest.raises(ValueError, match="European"):
        parse_scores(text)


def test_licence_recorded() -> None:
    for score in load_scores():
        assert score.licence.strip(), score.pgs_id
        assert score.publication.strip() and score.development_ancestry.strip()


def test_topic_mapping_falls_back_to_research() -> None:
    assert topic_for("EFO_9999999") is Topic.RESEARCH
    assert topic_for("MONDO_0005010") is Topic.RISKS  # coronary artery disorder
    # Every chosen score's trait is mapped: a disease under risks, a measurement under
    # wellness, height under traits.
    topics = {s.trait_efo: s.topic for s in load_scores()}
    assert set(topics.values()) == {Topic.RISKS, Topic.WELLNESS, Topic.TRAITS}
    assert topics["EFO_0004340"] is Topic.WELLNESS  # body mass index
    assert topics["OBA_VT0001253"] is Topic.TRAITS  # body height
    assert topics["EFO_0009270"] is Topic.WELLNESS  # heel bone mineral density


def test_weak_scores_are_the_spec_list() -> None:
    scores = load_scores()
    assert {s.pgs_id for s in scores if s.weak} == WEAK
    assert all(s.weak is False for s in scores if s.pgs_id not in WEAK)


def test_weak_must_be_boolean() -> None:
    assert parse_scores(_SCORE + "weak = true\n")[0].weak is True
    with pytest.raises(ValueError, match="PGS000001"):
        parse_scores(_SCORE + 'weak = "yes"\n')


def test_not_covered_list() -> None:
    entries = load_not_covered()
    assert [(e.condition, e.topic) for e in entries] == list(NOT_COVERED)
    assert all(e.why.strip() for e in entries)


def test_not_covered_never_names_a_chosen_condition() -> None:
    """A condition with a chosen score is never also not covered: refused, not written
    twice."""
    assert parse_not_covered(_SCORE + _not_covered("d")) == (
        NotCovered("d", Topic.RISKS, "w"),
    )
    with pytest.raises(ValueError, match="^c: a chosen score names it"):
        parse_not_covered(_SCORE + _not_covered("c"))
    crohn = _SCORE.replace('"c"', '"Crohn\'s disease"')
    with pytest.raises(ValueError, match="chosen score names it"):
        parse_not_covered(crohn + _not_covered("crohn's disease"))
    chosen = {s.condition for s in load_scores()}
    assert not chosen & {e.condition for e in load_not_covered()}


def test_not_covered_ids_are_unique_slugs() -> None:
    assert NotCovered("lipoprotein(a)", Topic.WELLNESS, "w").slug == "lipoprotein-a"
    assert NotCovered("Crohn's disease", Topic.RISKS, "w").slug == "crohn-s-disease"
    slugs = [e.slug for e in load_not_covered()]
    assert len(set(slugs)) == len(slugs) == 24
    assert all(re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", slug) for slug in slugs)
    with pytest.raises(ValueError, match="duplicate"):
        parse_not_covered(_not_covered("Hair colour") + _not_covered("hair-colour"))


def test_not_covered_unknown_topic_refused() -> None:
    with pytest.raises(ValueError, match="health"):
        parse_not_covered(_not_covered("d", topic="health"))


# --- the research layer: scores that do not pass the rule ------------------------------------

_RESEARCH = """
[[score]]
pgs_id = "PGS000001"
condition = "chosen"
trait_efo = "EFO_1"
publication = "p"
development_ancestry = "d"
evaluations = [["PPM1", "EUR"]]
licence = "l"

[[research_score]]
pgs_id = "PGS000002"
condition = "ADHD"
trait_efo = "EFO_2"
publication = "p"
development_ancestry = "d"
evaluations = [["PPM2", "AMR"]]
licence = "l"
why = "evaluated only on behaviour scales"

[[research_score]]
pgs_id = "PGS000003"
condition = "ADHD"
trait_efo = "EFO_2"
publication = "p"
development_ancestry = "d"
evaluations = [["PPM3", "EUR"]]
licence = "l"
why = "evaluated only on questionnaires"
"""


def test_a_condition_can_hold_several_research_scores() -> None:
    research = parse_research_scores(_RESEARCH)
    assert [(s.pgs_id, s.condition) for s in research] == [
        ("PGS000002", "ADHD"),
        ("PGS000003", "ADHD"),
    ]
    assert research[0].why == "evaluated only on behaviour scales"
    assert research[0].evaluations == (("PPM2", "AMR"),)


@pytest.mark.parametrize(
    ("text", "match"),
    [
        (re.sub(r'why = "evaluated only on behaviour scales"\n', "", _RESEARCH), "why"),
        (_RESEARCH.replace('"evaluated only on behaviour scales"', '"  "'), "why"),
        (_RESEARCH.replace("PGS000002", "PGS000001"), "chosen"),
        (_RESEARCH.replace("PGS000003", "PGS000002"), "duplicate"),
        (_RESEARCH.replace('"PPM2", "AMR"', '"PPM2", "Hispanic"'), "ancestry"),
    ],
)
def test_a_wrong_research_score_is_refused(text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        parse_research_scores(text)


def test_the_research_list_is_separate_from_the_chosen_one() -> None:
    assert [s.pgs_id for s in parse_scores(_RESEARCH)] == ["PGS000001"]
    assert parse_research_scores("") == ()


def test_the_shipped_research_scores() -> None:
    """Every `[[research_score]]` of the packaged list parses (`why` included) and none
    repeats a chosen score."""
    research = load_research_scores()
    assert len(research) == 27
    chosen = {s.pgs_id for s in load_scores()}
    assert not {s.pgs_id for s in research} & chosen
    assert all(s.why.strip() and s.evaluations for s in research)
    assert sum(s.condition == "ADHD" for s in research) == 2  # several per condition
    assert {s.condition for s in research} >= {
        "autism spectrum disorder",
        "intelligence quotient",
        "years of education",
        "hand grip strength",
        "overweight and obesity",
    }
