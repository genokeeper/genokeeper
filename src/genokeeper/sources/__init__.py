"""The source seam: what each source of statements hands the run, and how the run writes it.

A source (`Source`) collects, from its own artifacts and the run's shared evidence
(`RunContext`), every row it would write: per variant its `variant`, `call` and `statement`
(`Row`), its extra tables (slice 1b's `gene_call`), and where its input rows
went that no row stands for (`Collected.dropped`). It reads the store never and writes
nothing. The pure steps across sources come next (`promote.promote_compound_hets`), and
then `write_collected` writes them all with one accounting row per tool, counted from the
rows' kinds after those steps. So a new source is one more `Source` in the run's list.
"""

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Protocol

import duckdb

from genokeeper.config import Grading
from genokeeper.evidence.clingen import Curation
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import GeneList
from genokeeper.store.db import Kind, record_accounting
from genokeeper.store.writer import StoreWriter

type Fields = dict[str, object]

# Slice 3f's "information, not calibrated": in the reasons of a statement whose truth is not
# measured here (the HHV-6 band's width, ancestry, the blood-group phenotypes, the traits'
# phenotypes).
INFORMATION_NOTE = "information from the tools' own models; not checked against a reference genome here"


def present(fields: Mapping[str, object]) -> Fields | None:
    """`fields` without what no tool gave, or None when that is all of them."""
    return {k: v for k, v in fields.items() if v not in (None, (), [], {})} or None


@dataclass(frozen=True)
class Row:
    """One statement of a source, with the variant and the call it was graded on (None
    where it has none: a gene-level statement), and what its tools measured, whatever its
    grade (`present`; None: nothing), which only `sources.withheld` reads: never stored."""

    variant: Fields | None
    call: Fields | None
    statement: Fields
    measured: Fields | None = None


@dataclass(frozen=True)
class Collected:
    """Everything one source writes: its rows, its extra tables (table -> rows), and its
    accounting's input rows and drop reasons."""

    tool: str
    rows: list[Row]
    extra: dict[str, list[Fields]]
    input_rows: int
    dropped: dict[str, int]


@dataclass(frozen=True)
class RunContext:
    """The evidence every source grades with, shared by the run."""

    run_id: str
    lists: Mapping[str, GeneList]
    grading: Grading
    curations: Mapping[str, list[Curation]]
    gene_level: Mapping[str, str]
    faf: FafTable
    ba1_exceptions: frozenset[str]
    inheritance: Mapping[str, str]
    male: bool
    haploid: Callable[[str, int], bool]
    source_version: str  # the ClinVar release date
    # CPIC's level of each (gene, drug) pair (`evidence.pgx.load_cpic_levels`).
    cpic_levels: Mapping[tuple[str, str], str] = field(default_factory=dict)
    # vid -> the ClinPGx annotation the run stores at that allele, its id and level
    # (`evidence.pharmgkb.by_vid`), which grades a ClinVar drug_response claim at it.
    clinpgx: Mapping[str, tuple[str, str]] = field(default_factory=dict)
    # The genes this run's PharmCAT report calls, unresolved ones left out
    # (`evidence.pgx.unresolved`): a drug_response claim there is graded on the gene's call.
    gene_calls: frozenset[str] = frozenset()
    # Gene -> who gave its validity when it is not ClinGen (`evidence.gencc.merge`).
    validity_source: Mapping[str, str] = field(default_factory=dict)


class Source(Protocol):
    @property
    def tool(self) -> str: ...

    def collect(self, ctx: RunContext) -> Collected: ...


def write_collected(
    writer: StoreWriter,
    con: duckdb.DuckDBPyConnection,
    run_id: str,
    collected: Sequence[Collected],
) -> None:
    """Write each source's rows and extra tables through `writer`, flushed before it
    returns, and one accounting row per tool: its statements' kinds plus its `dropped`.

    Sources write in order, so where two rows hold a variant or a call at one vid the
    earlier source's is kept (`StoreWriter`)."""
    for c in collected:
        for row in c.rows:
            if row.variant is not None:
                writer.add("variant", row.variant)
            if row.call is not None:
                writer.add("call", row.call)
            writer.add("statement", row.statement)
        for table, rows in c.extra.items():
            for extra in rows:
                writer.add(table, extra)
    writer.flush()
    for c in collected:
        kinds = dict.fromkeys(Kind, 0)
        for row in c.rows:
            kinds[Kind(str(row.statement["kind"]))] += 1
        record_accounting(
            con,
            run_id,
            c.tool,
            input_rows=c.input_rows,
            findings=kinds[Kind.FINDING],
            watch=kinds[Kind.WATCH],
            non_findings=kinds[Kind.NON_FINDING],
            unassessed=kinds[Kind.UNASSESSED],
            dropped=c.dropped,
        )
