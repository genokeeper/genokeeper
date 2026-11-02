"""A pharmacogenomic gene result graded (track design, slice 1b, "The store" and "Gene-call
reliability"), and CPIC's gene-drug levels.

- **Validity** is the strongest recommendation strength PharmCAT reports for the gene's
  drugs: Strong, Moderate, Optional -> strong, moderate, limited. With none of the three
  (only `No recommendation`, `Unspecified`, or no drug) it is limited, warned
  `NO_STRENGTH` (builder's call).
- **Actionable** when one of the gene's drugs is a CPIC level A or B pair with it (CPIC's
  pair list, `load_cpic_levels`) or carries a DPWG guideline.
- **Reliability** of the gene call: high when two callers agree and no defining position is
  missing; medium for one caller (PharmCAT alone, or one outside caller) with none missing;
  low when a position is missing, when PharmCAT did not take the outside call it was given,
  or when the call is unresolved. There is no `p_wrong`.
- **Unresolved** (`unresolved`): the outside callers disagree (their reason), or PharmCAT
  made no call. Such a gene is a gene-level unassessed statement, never a finding.

CPIC's data comes from its API (`api.cpicpgx.org`) as JSON rows, fetched whole by the
bootstrap (step_cpic) and filed under the fetch date with its sha256: the gene-drug pairs
with their levels, and MT-RNR1's allele names (`load_cpic_alleles`). A file that is not a
non-empty list of rows with the keys read is refused (`ValueError`), so the bootstrap never
files a maintenance page or an error answer, and an empty answer never makes every user's
MT-RNR1 `Reference`. Its anchor comes from the bytes on this disk, which such a response
cannot pass because none is a list of rows with those keys; a republished table with those
keys is taken as it is, never refused.

Nothing here reads a recommendation's text: a `DrugRec` carries none.
"""

import json
from collections.abc import Collection, Mapping
from pathlib import Path
from typing import Any

from genokeeper.pgx.outside import OutsideCall
from genokeeper.pgx.pharmcat import DrugRec, GeneResult
from genokeeper.store.db import Level, Validity

# CPIC's recommendation strengths, as PharmCAT reports them, and their validity.
STRENGTH = {
    "Strong": Validity.STRONG,
    "Moderate": Validity.MODERATE,
    "Optional": Validity.LIMITED,
}
_RANK = {Validity.STRONG: 3, Validity.MODERATE: 2, Validity.LIMITED: 1}
# CPIC levels that make a gene-drug pair actionable.
ACTIONABLE_LEVELS = frozenset({"A", "A/B", "B"})
DPWG = "DPWG"
NO_STRENGTH = "no recommendation strength reported"
# Every unresolved gene call's warning starts with this (query.unassessed_reason reads it).
GENE_UNRESOLVED = "gene call unresolved"
PHARMCAT_NO_CALL = "PharmCAT made no call"
NOT_TAKEN = "PharmCAT did not take the outside call"
MISSING = "defining positions missing"
# CPIC's data a run reads, as the bootstrap files it in `<lab>/cpic/<YYYYMMDD>/` (step_cpic):
# the gene-drug pairs and MT-RNR1's alleles.
CPIC_FILES = {"pairs": "cpic_pairs.json", "MT-RNR1": "MT-RNR1_alleles.json"}


def _rows(path: Path, keys: Collection[str]) -> list[dict[str, Any]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(rows, list)
        or not rows
        or not all(isinstance(r, dict) and set(keys) <= r.keys() for r in rows)
    ):
        raise ValueError(f"{path}: not CPIC's rows of {', '.join(keys)}")
    return rows


def load_cpic_levels(path: Path) -> dict[tuple[str, str], str]:
    """CPIC's level of each (gene, drug) pair, the drug lower-cased, from CPIC's
    `pair_view` rows (`genesymbol`, `drugname`, `cpiclevel`)."""
    return {
        (str(r["genesymbol"]), str(r["drugname"]).lower()): str(r["cpiclevel"])
        for r in _rows(path, ("genesymbol", "drugname", "cpiclevel"))
        if r["genesymbol"] and r["drugname"] and r["cpiclevel"]
    }


def load_cpic_alleles(path: Path) -> list[str]:
    """The allele names of CPIC's `allele` rows for one gene (`name`)."""
    return [str(r["name"]) for r in _rows(path, ("name",)) if r["name"]]


def _pair_actionable(
    gene: str, drug: str, cpic_levels: Mapping[tuple[str, str], str]
) -> bool:
    return cpic_levels.get((gene, drug.lower())) in ACTIONABLE_LEVELS


def cpic_actionable(
    gene: str, condition: str, cpic_levels: Mapping[tuple[str, str], str]
) -> bool:
    """Whether a ClinVar drug_response claim in `gene` is actionable: its `condition` names a
    drug whose pair with the gene is CPIC level A or B."""
    text = condition.lower()
    return any(
        g == gene and drug in text and level in ACTIONABLE_LEVELS
        for (g, drug), level in cpic_levels.items()
    )


def unresolved(result: GeneResult, outside: OutsideCall | None) -> str | None:
    """Why the gene's call is unresolved (`GENE_UNRESOLVED: <why>`), None when it is not."""
    if outside is not None and outside.unresolved_reason is not None:
        return f"{GENE_UNRESOLVED}: {outside.unresolved_reason}"
    if result.diplotype is None:
        return f"{GENE_UNRESOLVED}: {PHARMCAT_NO_CALL}"
    return None


def _strongest(drugs: tuple[DrugRec, ...]) -> tuple[Validity | None, DrugRec | None]:
    best: tuple[Validity | None, DrugRec | None] = (None, None)
    for drug in drugs:
        validity = STRENGTH.get(drug.classification or "")
        if validity is not None and (
            best[0] is None or _RANK[validity] > _RANK[best[0]]
        ):
            best = (validity, drug)
    return best


def source_level(result: GeneResult) -> str:
    """The source's own level: the strongest recommendation's source and strength, e.g.
    `CPIC Strong`, or `NO_STRENGTH`."""
    _, drug = _strongest(result.drugs)
    return NO_STRENGTH if drug is None else f"{drug.source} {drug.classification}"


def grade_gene_result(
    result: GeneResult,
    *,
    cpic_levels: Mapping[tuple[str, str], str],
    outside: OutsideCall | None,
) -> tuple[Validity, bool, Level, tuple[str, ...]]:
    """The gene result's validity, whether it is actionable, its gene call's reliability
    level and the warnings behind them, by the rules above. `outside` is the gene's outside
    call, None when PharmCAT called it from our calls alone."""
    warnings: list[str] = []
    validity, _ = _strongest(result.drugs)
    if validity is None:
        validity = Validity.LIMITED
        warnings.append(NO_STRENGTH)
    actionable = any(
        _pair_actionable(result.gene, d.drug, cpic_levels) or d.source == DPWG
        for d in result.drugs
    )
    why = unresolved(result, outside)
    if outside is not None and why is None and result.call_source != "OUTSIDE":
        warnings.append(NOT_TAKEN)
    if outside is not None and outside.note is not None:
        warnings.append(outside.note)
    if result.missing_positions and result.call_source == "MATCHER":
        warnings.append(f"{MISSING}: {', '.join(result.missing_positions)}")
    if why is not None:
        warnings.append(why)
    if (
        why is not None
        or NOT_TAKEN in warnings
        or any(w.startswith(MISSING) for w in warnings)
    ):
        level = Level.LOW
    elif outside is not None and outside.agreed and len(outside.callers) >= 2:
        level = Level.HIGH
    else:
        level = Level.MEDIUM
    return validity, actionable, level, tuple(warnings)
