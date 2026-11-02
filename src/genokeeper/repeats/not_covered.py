"""The STRchive v2.28.0 repeat loci this pipeline does not assess (slice 3e), each written
as an unassessed statement with the reason no validated short-read method covers it.

The 32 loci left out altogether are listed in resources/repeats/SOURCE.md ("Out").
"""

from dataclasses import dataclass

from genokeeper.store.db import Topic

PARALOGS = (
    "not covered: paralogs (NOTCH2NLA, NOTCH2NLB, NOTCH2NLR) take its reads, and GIAB's"
    " HG002 truth has no genotype there"
)
COMMON_EXPANSIONS = (
    "not covered: expansions are common in healthy people and rarely cause disease, so"
    " even a confirmed one is hard to interpret"
)
DODECAMER = (
    "not covered: the GC-rich dodecamer expansion is far longer than the fragment, and no"
    " short-read method is validated for it"
)
INSERTED = (
    "not covered: the disease allele inserts a motif the reference lacks ({motif}), which"
    " no validated short-read method detects"
)
CGG = "not covered: a GC-rich CGG expansion that short reads underestimate"
PROVISIONAL = (
    "not covered: added to STRchive in 2026 as provisional; no validated method yet"
)


@dataclass(frozen=True)
class NotCovered:
    """One locus not assessed: its gene, the condition (STRchive v2.28.0's disease name,
    shortened), STRchive's evidence, the dashboard topic and why it is not covered."""

    gene: str
    condition: str
    evidence: str
    topic: Topic
    reason: str


def _inserted(gene: str, condition: str, evidence: str, motif: str) -> NotCovered:
    return NotCovered(
        gene, condition, evidence, Topic.RISKS, INSERTED.format(motif=motif)
    )


def _fame(gene: str, kind: int, evidence: str) -> NotCovered:
    return _inserted(
        gene, f"familial adult myoclonic epilepsy type {kind}", evidence, "TTTCA"
    )


def _risk(gene: str, condition: str, evidence: str, reason: str) -> NotCovered:
    return NotCovered(gene, condition, evidence, Topic.RISKS, reason)


NOT_COVERED: tuple[NotCovered, ...] = (
    NotCovered(
        "NOTCH2NLC",
        "neuronal intranuclear inclusion disease",
        "Definitive",
        Topic.RISKS,
        PARALOGS,
    ),
    NotCovered(
        "ATXN8OS",
        "spinocerebellar ataxia type 8",
        "Definitive",
        Topic.RISKS,
        COMMON_EXPANSIONS,
    ),
    NotCovered(
        "CSTB",
        "progressive myoclonic epilepsy type 1 (Unverricht-Lundborg disease)",
        "Definitive",
        Topic.CARRIER,
        DODECAMER,
    ),
    _inserted("DAB1", "spinocerebellar ataxia type 37", "Strong", "ATTTC"),
    _inserted("BEAN1", "spinocerebellar ataxia type 31", "Definitive", "TGGAA"),
    _fame("SAMD12", 1, "Definitive"),
    _fame("STARD7", 2, "Strong"),
    _fame("MARCHF6", 3, "Definitive"),
    _fame("YEATS2", 4, "Limited"),
    _fame("TNRC6A", 6, "Limited"),
    _fame("RAPGEF2", 7, "Limited"),
    _fame("RAI1", 8, "Limited"),
    _risk("GIPC1", "oculopharyngodistal myopathy type 2", "Definitive", CGG),
    _risk("LRP12", "oculopharyngodistal myopathy type 1", "Strong", CGG),
    _risk("RILPL1", "oculopharyngodistal myopathy type 4", "Moderate", CGG),
    _risk("ABCD3", "oculopharyngodistal myopathy type 5", "Moderate", CGG),
    _risk(
        "NUTM2B-AS1",
        "oculopharyngeal myopathy with leukoencephalopathy type 1",
        "Moderate",
        CGG,
    ),
    _risk(
        "GOLGA8A",
        "atypical frontotemporal lobar degeneration",
        "Provisional",
        PROVISIONAL,
    ),
    _risk("EP400", "spinocerebellar ataxia", "Provisional", PROVISIONAL),
    _risk(
        "BCLAF3",
        "neurodevelopmental disorder with fragile site FRAXG",
        "Provisional",
        PROVISIONAL,
    ),
    _risk("FAM193B", "oculopharyngodistal myopathy", "Provisional", PROVISIONAL),
    _risk("TBC1D7", "oculopharyngodistal myopathy", "Provisional", PROVISIONAL),
)
