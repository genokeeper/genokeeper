"""The research layers (slice 3f, "Research layers, one predicate"): sources whose statements
are research associations, hidden from a query that names no gene, marker or text, never
trusted in the change feed (so they neither enter nor leave it) and counted, not listed, on
the dashboard's Changes tab. Every such test calls `is_research_layer`, so a new layer is one
name added here."""

from genokeeper.sources.gwas import GWAS_TOOL
from genokeeper.sources.hla import HLA_ASSOCIATION_TOOL
from genokeeper.sources.pharmgkb import PHARMGKB_RESEARCH_TOOL

RESEARCH_LAYERS = frozenset({GWAS_TOOL, PHARMGKB_RESEARCH_TOOL, HLA_ASSOCIATION_TOOL})


def is_research_layer(source: str | None) -> bool:
    """Whether a statement's `source` is a research layer's (a statement without one is
    not)."""
    return source in RESEARCH_LAYERS
