"""A ClinVar claim's diseases matched to ClinGen's gene-disease curations by MONDO ID.

ClinVar's `CLNDISDB` separates a record's diseases by `|` and each disease's IDs by `,`, and
writes a MONDO ID `MONDO:MONDO:0011450`; ClinGen writes it `MONDO:0011450`. Only MONDO IDs
match (track design, A3: 0.3% of the disease entries on P/LP records carry only an OMIM ID).
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from genokeeper.evidence import clingen
from genokeeper.evidence.clingen import Curation

_CLINVAR_MONDO = "MONDO:MONDO:"


def parse_clndisdb(value: str) -> list[frozenset[str]]:
    """Each disease's MONDO IDs, written as ClinGen writes them, one set per `|` part (an
    empty set for a disease without one, or for `.`)."""
    return [
        frozenset(
            "MONDO:" + part.strip().removeprefix(_CLINVAR_MONDO)
            for part in disease.split(",")
            if part.strip().startswith(_CLINVAR_MONDO)
        )
        for disease in value.split("|")
    ]


@dataclass(frozen=True)
class CurationMatch:
    """The classification that grades a claim, and whether a curation of one of its own
    diseases gave it (else it is the gene-level value, or None for an uncurated gene)."""

    classification: str | None
    matched: bool


def match_curation(
    gene: str,
    diseases: Sequence[frozenset[str]],
    curations: Mapping[str, list[Curation]],
    gene_level: str | None,
) -> CurationMatch:
    """The curation of `gene` that grades a claim about `diseases` (`parse_clndisdb`).

    Of the gene's curations whose MONDO ID is one of the claim's, the one `clingen.gene_level`
    picks: positive, then negative, then the rest. Diseases with no curation are ignored
    while one matches. With no match the gene-level value `gene_level` applies, unmatched.
    """
    wanted = frozenset().union(*diseases)
    matched = [c for c in curations.get(gene, []) if c.mondo and c.mondo in wanted]
    if not matched:
        return CurationMatch(gene_level, matched=False)
    return CurationMatch(clingen.gene_level(matched), matched=True)
