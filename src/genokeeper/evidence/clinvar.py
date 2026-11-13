"""ClinVar germline claims graded into clinical validity, kind and topic.

`CLNSIG` is never split on commas (`Pathogenic,_low_penetrance` is one term). Its part before
the first `|` is the germline classification, which CLNREVSTAT reviews; its terms joined with
`/` span one classification. Each later `|` part is a separate assertion (`risk_factor`,
`drug_response`, ...) with no review status in the VCF. `parse_clnsig` accepts a whitelist and
fails closed: any other value is unmapped, which grades limited with an alert and is shown.
The class and level come from the germline part, except that the legacy `risk_factor` after a
`|`, following a VUS or no germline part, makes a risk allele whose level is capped at
limited; no other term after a `|` does. The level combines the review stars with ClinGen
gene-disease validity (matched to the claim's disease where a curation of it exists; an
expert-panel or practice-guideline claim stands without it) and the gnomAD grpmax filtering
allele frequency; the kind and topic follow the gene's inheritance and the user's dosage,
except on a haploid locus (a male's chrX or chrY outside the PARs), where one copy is
hemizygous whatever the gene's inheritance. A claim that would be a finding where the user's dosage is unknown is
unassessed, never a finding.

A `drug_response` claim (a whole CLNSIG of that term; slice 1b) is never graded by its stars,
which are the submitter's review status (owner's ruling A3s): it is limited, with no ClinGen
gate and no frequency filter, its reason naming the one authority that grades the variant:
the ClinPGx annotation kept at its allele, else the gene's PharmCAT call, else none. Carried,
it is a finding under drugs; at dosage 0 a non-finding.

`sf_actionable` says whether ACMG SF v3.3 would report a claim, from its class, its placement
and the gene's report rule in `acmg_sf_v3.3.tsv`.
"""

from collections.abc import Mapping
from dataclasses import dataclass

from genokeeper.config import Grading
from genokeeper.evidence.clingen import NEGATIVE_VALIDITY, POSITIVE_VALIDITY
from genokeeper.genelists import (
    INHERITANCE,
    REPORT_ALL,
    REPORT_BIALLELIC,
    REPORT_TRUNCATING,
    GeneList,
    carrier_genes,
    gene_mode,
)
from genokeeper.pgx.pharmcat import PHARMCAT_GENES
from genokeeper.store.db import Kind, Topic, Validity

BASE_CLASSES = frozenset(
    {"pathogenic", "benign", "vus", "risk", "drug_response", "other", "conflicting"}
)

# Every germline CLNSIG term, including those tallied in ClinVar's GRCh37 VCF of 2026-10-11
# (evidence-sources research note of 2026-10-18, "ClinVar").
CLNSIG_MAP: dict[str, str] = {
    "Pathogenic": "pathogenic",
    "Likely_pathogenic": "pathogenic",
    "Pathogenic,_low_penetrance": "pathogenic",
    "Likely_pathogenic,_low_penetrance": "pathogenic",
    "Uncertain_significance": "vus",
    "VUS-high": "vus",
    "VUS-mid": "vus",
    "VUS-low": "vus",
    "Likely_benign": "benign",
    "Benign": "benign",
    "Conflicting_classifications_of_pathogenicity": "conflicting",
    "Conflicting_interpretations_of_pathogenicity": "conflicting",  # before 2024
    "Established_risk_allele": "risk",
    "Likely_risk_allele": "risk",
    "Uncertain_risk_allele": "vus",  # findings take established and likely risk alleles only
    "risk_factor": "risk",  # legacy term for a risk allele
    "drug_response": "drug_response",
    "confers_sensitivity": "other",
    "protective": "other",
    "association": "other",
    "association_not_found": "other",
    "Affects": "other",
    "other": "other",
    "not_provided": "other",
    "no_classification_for_the_single_variant": "other",
    "no_classifications_from_unflagged_records": "other",
}

_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,  # before 2024
}
_CONFLICT_STATUS = frozenset(
    {
        "criteria_provided,_conflicting_classifications",
        "criteria_provided,_conflicting_interpretations",
    }
)
# Classes graded only as watch or non-finding carry no level.
_UNGRADED = frozenset({"vus", "benign", "other"})
# The terms a germline classification (the part before the first `|`) may join with `/`.
GERMLINE_TERMS = frozenset(
    term
    for term, cls in CLNSIG_MAP.items()
    if cls in {"pathogenic", "benign", "vus", "conflicting"}
) | {"Established_risk_allele", "Likely_risk_allele"}
# A whole CLNSIG of one such term is a single non-germline assertion (legacy risk_factor,
# drug_response, ...), which CLNREVSTAT then reviews.
_ASSERTION_CLASSES = frozenset({"risk", "drug_response", "other"})
# The only terms allowed after a `|`: separate non-germline assertions. A germline
# classification term there (Pathogenic, Benign, a VUS or a risk-allele term) is unmapped.
LEGACY_RISK_TERM = "risk_factor"
ASSERTION_TERMS = frozenset(
    term for term, cls in CLNSIG_MAP.items() if cls in {"drug_response", "other"}
) | {LEGACY_RISK_TERM}

_DOSAGES = frozenset({"0", "1", "2", "unknown"})
# ClinVar's molecular consequences (MC) that truncate for ACMG SF's TTN rule: nonsense,
# frameshift, splice acceptor and splice donor.
SF_TRUNCATING = frozenset({"SO:0001587", "SO:0001589", "SO:0001574", "SO:0001575"})

CONFLICT_REASON = "conflicting classifications"
BA1_REASON = "population frequency above BA1"
BS1_REASON = "population frequency above BS1"
NO_CRITERIA_REASON = "no assertion criteria"
NO_CURATION_REASON = "no ClinGen gene validity"
NOT_MATCHED_REASON = "ClinGen curation not matched to this disease; gene-level"
RISK_PART_REASON = "risk_factor has no review status"
INHERITANCE_REASON = "inheritance unknown"
GENOTYPE_REASON = "genotype unknown"
HEMIZYGOUS_REASON = "hemizygous (haploid locus)"
NO_CLINPGX_REASON = (
    "ClinVar's stars on a drug-response record are the submitter's review status; no"
    " ClinPGx annotation at levels 1B–2B is located at this variant"
)


@dataclass(frozen=True)
class ValidityGrade:
    """A ClinVar claim's kind, mapped level, the source's own level, reason, topic and alert."""

    kind: Kind
    level: Validity | None
    source_level: str
    reason: str | None
    topic: Topic | None
    alert: str | None


def stars(clnrevstat: str) -> int:
    """ClinVar's review stars for a CLNREVSTAT value; anything unlisted is 0."""
    return _STARS.get(clnrevstat, 0)


@dataclass(frozen=True)
class ParsedClnsig:
    """A CLNSIG value's base class, or None when the value is unmapped; `risk_part` is True
    when the class is a risk term after the germline part, not the germline classification."""

    cls: str | None
    risk_part: bool = False


_UNMAPPED = ParsedClnsig(None)


def _one_class(terms: list[str]) -> str | None:
    """One class for mapped terms joined by `/`: all alike, or P with B (conflicting)."""
    classes = {CLNSIG_MAP[t] for t in terms}
    if len(classes) == 1:
        return classes.pop()
    return "conflicting" if classes == {"pathogenic", "benign"} else None


def parse_clnsig(clnsig: str) -> ParsedClnsig:
    """Classify a CLNSIG value against the whitelist; anything outside it is unmapped.

    - Every term after a `|` must be in ASSERTION_TERMS (legacy `risk_factor`,
      drug_response or an other-class term); a germline term there is unmapped.
    - A whole value of one non-germline term (`risk_factor`, `drug_response`, ...) is that
      assertion.
    - Otherwise the germline part's `/`-terms must all be GERMLINE_TERMS of one class (P with
      B is conflicting). An empty germline part counts as absent when other parts follow.
    - `risk_factor` after a VUS or an absent germline part makes a risk allele
      (`risk_part`); after any other germline class it changes nothing. An absent germline
      part with no `risk_factor` after it is unmapped.
    """
    germline, bar, rest = clnsig.partition("|")
    others = rest.split("|") if bar else []
    if any(t not in ASSERTION_TERMS for t in others):
        return _UNMAPPED
    risk_after = LEGACY_RISK_TERM in others
    if not germline:
        # Absent: only a risk_factor after it makes a claim.
        return ParsedClnsig("risk", risk_part=True) if risk_after else _UNMAPPED
    if not others and CLNSIG_MAP.get(germline) in _ASSERTION_CLASSES:
        return ParsedClnsig(CLNSIG_MAP[germline])
    terms = germline.split("/")
    if any(t not in GERMLINE_TERMS for t in terms):
        return _UNMAPPED
    cls = _one_class(terms)
    if cls == "vus" and risk_after:
        return ParsedClnsig("risk", risk_part=True)
    return ParsedClnsig(cls)


def clnsig_class(clnsig: str) -> str | None:
    """The base class of a CLNSIG value, or None when it is unmapped."""
    return parse_clnsig(clnsig).cls


def _level(
    *,
    cls: str | None,
    clnrevstat: str,
    gene: str,
    gene_validity: str | None,
    faf_grpmax: float | None,
    grading: Grading,
    vid: str | None,
    mode: str | None,
    ba1_exceptions: frozenset[str],
    validity_matched: bool,
    validity_source: str | None = None,
) -> tuple[Validity, list[str]]:
    """The mapped level and its reasons, checked in this order:

    1. conflicting: a ClinVar conflict; BA1/BS1; a Disputed or Refuted curation, when it is
       matched to the claim's disease (`validity_matched`) or the claim has under 3 stars;
    2. strong at 3 or 4 stars, the reason naming a gene-level negative curation;
    3. limited without a positive curation (Definitive, Strong, Moderate);
    4. strong at 2 stars, moderate at 1, else limited.

    A gene-level value graded in place of a disease-matched one adds NOT_MATCHED_REASON.
    BA1 and BS1 apply to a pathogenic claim with a frequency, unless its `vid` is on the
    ClinGen BA1 exception list; BS1 is the gene's (`Grading.bs1_for`), none in an autosomal
    recessive gene outside [bs1_by_gene] (`mode` "AR"). A reason naming the gene's validity
    names its `validity_source` after it, when one is given (GenCC's, slice 3f).
    """
    validity = None if gene_validity is None else gene_validity.strip().lower()
    named = f"gene validity {gene_validity}" + (
        f" ({validity_source})" if validity_source else ""
    )
    n = stars(clnrevstat)
    # The gene-level value stands in for a curation of the claim's own disease: say so.
    unmatched = [] if validity_matched or validity is None else [NOT_MATCHED_REASON]
    negative = [named, *unmatched]
    conflicts: list[str] = []
    if validity in NEGATIVE_VALIDITY and (validity_matched or n < 3):
        conflicts.extend(negative)
    if cls == "conflicting" or clnrevstat in _CONFLICT_STATUS:
        conflicts.append(CONFLICT_REASON)
    if (
        cls == "pathogenic"
        and faf_grpmax is not None
        and (vid is None or vid not in ba1_exceptions)
    ):
        bs1 = grading.bs1_for(gene, autosomal_recessive=mode == "AR")
        if faf_grpmax > grading.ba1:
            conflicts.append(BA1_REASON)
        elif bs1 is not None and faf_grpmax > bs1:
            conflicts.append(BS1_REASON)
    if conflicts:
        return Validity.CONFLICTING, conflicts
    if n >= 3:
        # An expert panel or a practice guideline stands on its own (track design, A3).
        return Validity.STRONG, negative if validity in NEGATIVE_VALIDITY else []
    if validity not in POSITIVE_VALIDITY:
        reason = NO_CURATION_REASON if validity is None else named
        return Validity.LIMITED, [reason, *unmatched]
    if n == 2:
        return Validity.STRONG, unmatched
    if n == 1:
        return Validity.MODERATE, unmatched
    return Validity.LIMITED, [NO_CRITERIA_REASON, *unmatched]


def placement(
    *,
    cls: str | None,
    gene: str,
    lists: dict[str, GeneList],
    dosage: str,
    inheritance: str | None,
    male: bool,
    haploid: bool,
) -> tuple[Kind, Topic | None, list[str]]:
    """Kind, topic and reasons from the class, the gene's inheritance and the dosage.

    A claim that would be a finding at an `unknown` dosage is unassessed instead, keeping
    the topic and reasons it would have had: what is unknown is the user's genotype, so it
    is never a finding. A VUS watch item stays one whatever the dosage.
    """
    kind, topic, reasons = _claim_placement(
        cls=cls,
        gene=gene,
        lists=lists,
        dosage=dosage,
        inheritance=inheritance,
        male=male,
        haploid=haploid,
    )
    if kind is Kind.FINDING and dosage == "unknown":
        return Kind.UNASSESSED, topic, reasons
    return kind, topic, reasons


def _claim_placement(
    *,
    cls: str | None,
    gene: str,
    lists: dict[str, GeneList],
    dosage: str,
    inheritance: str | None,
    male: bool,
    haploid: bool,
) -> tuple[Kind, Topic | None, list[str]]:
    """Kind, topic and reasons as if an `unknown` dosage could be any (`placement` makes
    such a finding unassessed).

    On a haploid locus (`haploid`: a male's chrX or chrY outside the PARs) one copy is the
    whole genotype, so a pathogenic or conflicting claim with a copy is a finding under
    risks whatever the gene's inheritance, and never carrier status.
    """
    sf = lists["acmg_sf"].genes
    tier3 = lists["carrier_tier3"].genes
    if dosage == "0" or cls in {"benign", "other"}:
        return Kind.NON_FINDING, None, []
    if cls == "drug_response":
        return (
            Kind.FINDING,
            Topic.DRUGS,
            [GENOTYPE_REASON] if dosage == "unknown" else [],
        )
    if cls == "vus":
        listed = gene in sf or gene in tier3
        return (Kind.WATCH if listed else Kind.NON_FINDING), None, []
    unknown = [GENOTYPE_REASON] if dosage == "unknown" else []
    if cls is None:
        # An unmapped value is shown (limited, with its alert), never hidden.
        return Kind.FINDING, Topic.RISKS, unknown
    if cls == "risk":
        return Kind.FINDING, Topic.RISKS, unknown
    # A pathogenic or conflicting claim: the lists win over the caller's inheritance.
    if haploid:
        return Kind.FINDING, Topic.RISKS, unknown or [HEMIZYGOUS_REASON]
    mode = gene_mode(gene, lists, inheritance)
    if mode not in INHERITANCE:
        return Kind.FINDING, Topic.RISKS, [INHERITANCE_REASON, *unknown]
    carrier = gene in carrier_genes(lists) and (
        mode == "AR" or (mode == "XL" and not male)
    )
    if dosage == "unknown":
        return Kind.FINDING, (Topic.CARRIER if carrier else Topic.RISKS), unknown
    if mode == "AD" or dosage == "2" or (mode == "XL" and male):
        return Kind.FINDING, Topic.RISKS, []
    if carrier:
        return Kind.FINDING, Topic.CARRIER, []
    if mode == "AR":
        # A single het in a recessive gene outside the carrier lists.
        return Kind.NON_FINDING, None, []
    # An XL het in a female, outside the carrier lists.
    return Kind.FINDING, Topic.RISKS, []


def sf_actionable(
    *,
    clnsig: str,
    gene: str,
    lists: dict[str, GeneList],
    kind: Kind,
    topic: Topic | None,
    vid: str,
    consequence: str,
) -> bool:
    """Whether ACMG SF v3.3 would report a claim placed at `kind` and `topic` (`placement`):
    a P/LP claim in an SF gene, a finding, and one the gene's report rule admits. `all` takes
    any such finding, a woman's carrier copy in GLA or OTC included; `biallelic` a finding
    under risks, which placement gives a recessive gene only at two copies or a haploid
    locus; `truncating` one whose `consequence` (ClinVar's MC) is in SF_TRUNCATING;
    `homozygous:<vid>` a finding under risks at that vid. A pair of single copies is the
    compound-heterozygote promotion's (`sources.promote`)."""
    rule = lists["acmg_sf"].report.get(gene)
    if rule is None or kind is not Kind.FINDING or clnsig_class(clnsig) != "pathogenic":
        return False
    if rule == REPORT_ALL:
        return True
    if topic is not Topic.RISKS:
        return False
    if rule == REPORT_BIALLELIC:
        return True
    if rule == REPORT_TRUNCATING:
        terms = {part.split("|")[0] for part in consequence.split(",")}
        return bool(terms & SF_TRUNCATING)
    return rule.partition(":")[2] == vid


def drug_response_reason(
    gene: str, vid: str | None, clinpgx: Mapping[str, tuple[str, str]]
) -> str:
    """What grades a drug_response claim in `gene` at allele `vid` (A3s): the ClinPGx
    annotation `clinpgx` holds there (`evidence.pharmgkb.by_vid`), else a PharmCAT gene's call."""
    found = clinpgx.get(vid) if vid else None
    if found is not None:
        return f"graded on ClinPGx's annotation {found[0]} (level {found[1]})"
    if gene in PHARMCAT_GENES:
        return f"graded on the {gene} gene call"
    return NO_CLINPGX_REASON


def grade_clinvar(
    *,
    clnsig: str,
    clnrevstat: str,
    gene: str,
    gene_validity: str | None,
    lists: dict[str, GeneList],
    faf_grpmax: float | None,
    dosage: str,
    grading: Grading,
    inheritance: str | None = None,
    male: bool = False,
    haploid: bool = False,
    vid: str | None = None,
    ba1_exceptions: frozenset[str] = frozenset(),
    validity_matched: bool = False,
    validity_source: str | None = None,
    clinpgx: Mapping[str, tuple[str, str]] | None = None,
) -> ValidityGrade:
    """Grade one ClinVar claim at a site where the user has `dosage` copies of the allele.

    `lists` holds `acmg_sf`, `carrier_tier3` and `carrier_tier4`; `inheritance` (AD, AR or
    XL) is used only for a gene on none of them. `haploid` says the site has one copy (a male's chrX or chrY
    outside the PARs, `reliability.on_haploid_contig`): a pathogenic or conflicting claim
    with a copy there is hemizygous, a finding under risks whatever the gene's inheritance.
    `faf_grpmax` is gnomAD's filtering allele frequency at `vid`, the claim's GRCh38
    allele; a `vid` in `ba1_exceptions` (`load_ba1_exceptions`) is exempt from BA1 and BS1.
    `gene_validity` is the ClinGen classification `match_curation` gives the claim, and
    `validity_matched` whether a curation of the claim's own disease gave it;
    `validity_source` who gave it, named in the reasons when given (GenCC's, slice 3f).
    A drug_response claim's reason is `drug_response_reason` of its `vid` and `clinpgx`.
    """
    if dosage not in _DOSAGES:
        raise ValueError(f"dosage {dosage!r}: expected one of {sorted(_DOSAGES)}")
    parsed = parse_clnsig(clnsig)
    cls = parsed.cls
    level: Validity | None
    if cls is None:
        level, level_reasons = Validity.LIMITED, []
    elif cls in _UNGRADED:
        level, level_reasons = None, []
    elif cls == "drug_response":
        reason = drug_response_reason(gene, vid, clinpgx or {})
        level, level_reasons = Validity.LIMITED, [reason]
    elif parsed.risk_part:
        # CLNREVSTAT reviews the germline part, not the risk term: grade with no stars.
        level, level_reasons = _level(
            cls=cls,
            clnrevstat="",
            gene=gene,
            gene_validity=gene_validity,
            faf_grpmax=faf_grpmax,
            grading=grading,
            vid=vid,
            mode=gene_mode(gene, lists, inheritance),
            ba1_exceptions=ba1_exceptions,
            validity_matched=validity_matched,
            validity_source=validity_source,
        )
        others = [r for r in level_reasons if r != NO_CRITERIA_REASON]
        level_reasons = [RISK_PART_REASON, *others]
    else:
        level, level_reasons = _level(
            cls=cls,
            clnrevstat=clnrevstat,
            gene=gene,
            gene_validity=gene_validity,
            faf_grpmax=faf_grpmax,
            grading=grading,
            vid=vid,
            mode=gene_mode(gene, lists, inheritance),
            ba1_exceptions=ba1_exceptions,
            validity_matched=validity_matched,
            validity_source=validity_source,
        )
    kind, topic, kind_reasons = placement(
        cls=cls,
        gene=gene,
        lists=lists,
        dosage=dosage,
        inheritance=inheritance,
        male=male,
        haploid=haploid,
    )
    reasons = level_reasons + kind_reasons
    return ValidityGrade(
        kind=kind,
        level=level,
        source_level=f"{clnsig} ({stars(clnrevstat)}★)",
        reason="; ".join(reasons) or None,
        topic=topic,
        alert=f"unmapped CLNSIG: {clnsig}" if cls is None else None,
    )
