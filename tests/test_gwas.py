"""Single-SNP research associations (`evidence.gwas`, `sources.gwas`): the GWAS Catalog's
genome-wide significant single-SNP rows, located on the HGDP+1kGP panel's alleles, and
each a statement on our call there: carried a finding (also when the risk allele is the
reference), not carried a non-finding, unknown unassessed; validity capped at limited, no
absolute risk, the topic from the trait's EFO."""

import json
import zipfile
from compression import zstd
from dataclasses import replace
from pathlib import Path
from typing import Any

import pysam
import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.gwas_mini import EFO, OBO, Assoc, write_gwas

from genokeeper import pipeline
from genokeeper.carrier import Answer, answer, genotype
from genokeeper.config import load_grading
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.gwas import (
    Association,
    load_associations,
    locate,
    panel_alleles,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import NOCALL_REASON, OTHER_ALLELE_WARNING, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.query import Filters, query
from genokeeper.sources import RunContext, write_collected
from genokeeper.sources.gwas import GwasSource
from genokeeper.store.db import (
    Level,
    RowType,
    Topic,
    Validity,
    connect,
    start_run,
    transaction,
)
from genokeeper.store.writer import StoreWriter

# chr1 of the test FASTA repeats ACGT: position p holds "ACGT"[(p - 1) % 4].
SEQUENCE = "ACGT" * 500
# The panel's alleles: (C, A) at 102, a palindromic (G, C) at 103, (A, G) at 105, a REF the
# FASTA does not have at 110, two ALTs at 113 and an indel at 102, which no SNP reads.
# Position 200 is not in the panel.
PANEL = [(102, "C", "A"), (103, "G", "C"), (105, "A", "G"), (110, "T", "C")]
PANEL += [(113, "A", "G"), (113, "A", "C"), (102, "C", "CT")]


def _fasta(tmp_path: Path) -> Path:
    fasta = tmp_path / "ref.fa"
    lines = [SEQUENCE[i : i + 60] for i in range(0, len(SEQUENCE), 60)]
    fasta.write_text(">chr1\n" + "\n".join(lines) + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _pvar(tmp_path: Path) -> Path:
    text = "#CHROM\tPOS\tID\tREF\tALT\n" + "".join(
        f"1\t{pos}\tv{pos}\t{ref}\t{alt}\n" for pos, ref, alt in PANEL
    )
    path = tmp_path / "ancestry.pvar.zst"
    path.write_bytes(zstd.compress(text.encode()))
    return path


def _located(
    tmp_path: Path, rows: list[Assoc]
) -> tuple[Any, list[Any], dict[str, int]]:
    found = load_associations(write_gwas(tmp_path / "gwas.zip", rows))
    panel = panel_alleles(
        _pvar(tmp_path), {(a.chrom, a.pos) for a in found.associations}
    )
    located, dropped = locate(found.associations, panel, fasta=_fasta(tmp_path))
    return found, located, dropped


def test_reader_keeps_single_snp_genome_wide_associations(tmp_path: Path) -> None:
    """A row stands when its p-value is at most 5e-8, it names one SNP and one risk base
    of it, and the catalog maps it to GRCh38; of one SNP, risk allele and trait set the
    strongest study stands. Every other row is counted by why."""
    rows = [
        Assoc(102, "A", p="2E-9", study="GCST000002", pubmed="2"),
        Assoc(102, "A", p="1E-12", study="GCST000003", pubmed="3"),  # the strongest
        Assoc(102, "A", p="5E-8"),  # at the threshold: kept, then a duplicate
        Assoc(102, "A", p="6E-8"),
        Assoc(102, "A", p="NR"),
        Assoc(102, "A", snps="rs102 x rs105"),
        Assoc(102, "A", snps="rs102; rs105"),
        Assoc(102, "?"),
        Assoc(102, "AT"),
        Assoc(102, "A", chrom=""),
        Assoc(102, "A", trait="methadone dose measurement", uris=f"{EFO}EFO_0011047"),
        # Another trait set: another association of the same allele.
        Assoc(
            102,
            "A",
            trait="body height",
            uris=f"{EFO}EFO_0004339",
            gene="PPARG - CAND1",
        ),
    ]
    found = load_associations(write_gwas(tmp_path / "gwas.zip", rows))
    assert (found.input_rows, found.release) == (12, "20261003")
    assert [(a.study, a.pubmed_id, a.p_value) for a in found.associations] == [
        ("GCST000003", "3", "1E-12"),
        ("GCST000001", "17463246", "1E-10"),
    ]
    assert {k: v for k, v in found.dropped.items() if v} == {
        "duplicate_association": 2,
        "p_value_above_5e-8": 2,
        "not_one_snp": 2,
        "risk_allele_unknown": 2,
        "no_grch38_position": 1,
        "dosing_trait": 1,
    }
    (first, _) = found.associations
    assert (first.rsid, first.chrom, first.pos, first.risk_allele) == (
        "rs102",
        "chr1",
        102,
        "A",
    )
    assert (first.trait, first.efo_ids, first.gene) == (
        "type 2 diabetes mellitus",
        ("MONDO_0005148",),
        "PPARG",
    )
    assert found.associations[1].gene is None  # between two genes


def test_one_association_per_trait_set_in_any_order(tmp_path: Path) -> None:
    """Two EFO IDs in either order are one trait set; rows mapped to no EFO are told
    apart by their trait."""
    pair = f"{OBO}MONDO_0005148, {EFO}EFO_0000400"
    rows = [
        Assoc(102, "A", uris=pair),
        Assoc(102, "A", uris=", ".join(reversed(pair.split(", "))), p="1E-9"),
        Assoc(102, "A", uris="", trait="body height"),
        Assoc(102, "A", uris="", trait="body mass index"),
    ]
    found = load_associations(write_gwas(tmp_path / "gwas.zip", rows))
    assert [(a.trait, a.efo_ids) for a in found.associations] == [
        ("type 2 diabetes mellitus", ("MONDO_0005148", "EFO_0000400")),
        ("body height", ()),
        ("body mass index", ()),
    ]
    assert found.dropped["duplicate_association"] == 1


def test_reader_refuses_another_file(tmp_path: Path) -> None:
    path = tmp_path / "gwas.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("other.tsv", "a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="not the GWAS Catalog associations download"):
        load_associations(path)
    path.write_text("not a zip")
    with pytest.raises(ValueError, match="not the GWAS Catalog associations download"):
        load_associations(path)


def test_effect_as_the_catalog_gives_it(tmp_path: Path) -> None:
    """An odds ratio as it stands; a beta, which the catalog gives unsigned with its
    direction for the risk allele, negative for a decrease; none where none is given."""
    rows = [
        Assoc(102, "A"),
        Assoc(105, "G", effect="0.02", ci="[0.01-0.03] unit decrease", uris=f"{EFO}a"),
        Assoc(105, "G", effect="0.5", ci="z score increase", uris=f"{EFO}b"),
        Assoc(105, "G", effect="", ci="", uris=f"{EFO}c"),
    ]
    found = load_associations(write_gwas(tmp_path / "gwas.zip", rows))
    assert [(a.measure, a.value, a.ci) for a in found.associations] == [
        ("OR", 1.14, (1.08, 1.2)),
        ("beta", -0.02, (-0.03, -0.01)),
        ("beta", 0.5, None),
        (None, None, None),
    ]


def test_strand_ambiguous_snp_skipped(tmp_path: Path) -> None:
    """A risk allele whose complement is also among the site's alleles (a palindromic SNP,
    or a site with both) cannot be told from the other strand's, and one that is none of
    the panel's alleles on the forward strand is another strand's or another allele: both
    are skipped, counted."""
    _, located, dropped = _located(
        tmp_path,
        [Assoc(103, "G"), Assoc(113, "C"), Assoc(102, "G"), Assoc(102, "T")],
    )
    assert located == []
    assert {k: v for k, v in dropped.items() if v} == {
        "strand_ambiguous": 2,
        "risk_allele_mismatch": 2,
    }


def test_located_on_the_panel(tmp_path: Path) -> None:
    """The panel gives each SNP's REF (checked against the FASTA) and ALTs: a risk allele
    that is the REF is located at the panel's first ALT, also where it has two; a SNP the
    panel lacks (a REF risk allele among them) or a REF the FASTA does not have is
    counted."""
    _, located, dropped = _located(
        tmp_path,
        [Assoc(102, "A"), Assoc(105, "A"), Assoc(200, "C"), Assoc(113, "A")]
        + [Assoc(110, "C")],
    )
    assert [(a.rsid, allele) for a, allele in located] == [
        ("rs102", ("chr1", 102, "C", "A")),
        ("rs105", ("chr1", 105, "A", "G")),
        ("rs113", ("chr1", 113, "A", "G")),
    ]
    assert {k: v for k, v in dropped.items() if v} == {
        "not_in_panel": 1,
        "reference_differs": 1,
    }


# --- the source -------------------------------------------------------------------------


def _called(gt: str | None, dosage: str, *warnings: str) -> Called:
    row = RowType.NONE if dosage == "unknown" else RowType.PASS
    return Called(gt, dosage, row, 40, 30, None, None, tuple(warnings))


def _reliability(called: Called) -> Reliability:
    level = Level.LOW if called.dosage == "unknown" else Level.HIGH
    return Reliability(p_wrong=None, level=level, flags=(), warnings=called.warnings)


def _association(pos: int, risk: str, **kw: Any) -> Association:
    fields: dict[str, Any] = {
        "rsid": f"rs{pos}",
        "chrom": "chr1",
        "pos": pos,
        "risk_allele": risk,
        "trait": "type 2 diabetes mellitus",
        "gene": "PPARG",
        "efo_ids": ("MONDO_0005148",),
        "measure": "OR",
        "value": 1.14,
        "ci": (1.08, 1.2),
        "p_value": "1E-10",
        "pubmed_id": "17463246",
        "study": "GCST000001",
    }
    return Association(**(fields | kw))


def _store(tmp_path: Path, located: list[Any]) -> tuple[Any, str]:
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {"gwas_catalog": "20261003"})
    ctx = RunContext(
        run_id=run_id,
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=True,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
    )
    rows = [(a, allele, c, _reliability(c)) for a, allele, c in located]
    source = GwasSource(rows, len(rows) + 1, {"not_in_panel": 1}, "20261003")
    with transaction(con):
        write_collected(
            StoreWriter(con, tmp_path / "work"), con, run_id, [source.collect(ctx)]
        )
        pipeline.check_calls(con, run_id)
    return con, run_id


def _statements(con: Any, run_id: str) -> dict[str, dict[str, Any]]:
    cursor = con.execute("SELECT * FROM statement WHERE run_id = ?", [run_id])
    names = [c[0] for c in cursor.description]
    return {row[1]: dict(zip(names, row, strict=True)) for row in cursor.fetchall()}


C_A = ("chr1", 102, "C", "A")


OTHER = OTHER_ALLELE_WARNING
# Carried a finding, not carried a non-finding, unknown unassessed. An ALT's copies are its
# dosage (0 when another allele is called); a reference risk allele's are the ploidy less
# the located ALT's dosage, and are not counted (None: the association is dropped) when
# the GT carries another allele, SNP or indel, or the risk allele is neither of the vid's.
CASES = [
    (1, "A", _called("0/1", "1"), ("finding", 1)),
    (2, "A", _called("0/0", "0"), ("non_finding", 0)),
    (3, "A", _called(None, "unknown", NOCALL_REASON), ("unassessed", None)),
    (4, "C", _called("0/0", "0"), ("finding", 2)),
    (5, "C", _called("1/1", "2"), ("non_finding", 0)),
    (6, "C", _called("2/2", "0", OTHER), None),
    (7, "C", _called("1/2", "1"), None),  # the ALT and a deletion CAA>C
    # One copy of chrX on a male (a haploid call): one reference copy, or none.
    (8, "C", _called("0", "0"), ("finding", 1)),
    (9, "C", _called("1", "1"), ("non_finding", 0)),
    (10, "C", _called("0/2", "0", OTHER), None),
    (11, "A", _called("0/2", "0", OTHER), ("non_finding", 0)),  # Task 34's row
    (12, "C", _called("1/1", "0", OTHER), None),  # a hom deletion CAA>C
    (13, "C", _called("0/1", "0", OTHER), None),  # a het insertion C>CT
    (14, "C", _called("1/1", "0", OTHER), None),  # a hom insertion C>CT
    (15, "C", _called("0/1", "1"), ("finding", 1)),
    (16, "G", _called("0/1", "1"), None),  # neither allele of the vid
    (17, "A", _called("1/1", "0", OTHER), ("non_finding", 0)),  # an indel called there
]
LOCATED = [
    (_association(pos, risk), ("chr1", pos, "C", "A"), called)
    for pos, risk, called, _ in CASES
]


def test_kinds_and_copies_from_the_call(tmp_path: Path) -> None:
    con, run_id = _store(tmp_path, LOCATED)
    found = {
        sid: (s["kind"], json.loads(s["claim"])["risk_allele_copies"])
        for sid, s in _statements(con, run_id).items()
    }
    assert found == {
        f"gwas:rs{pos}-{risk}:MONDO_0005148": expected
        for pos, risk, _, expected in CASES
        if expected is not None
    }
    (accounting,) = con.execute(
        "SELECT input_rows, findings, non_findings, unassessed, dropped FROM accounting"
    ).fetchall()
    assert accounting[:4] == (18, 4, 5, 1)
    assert json.loads(accounting[4]) == {"not_in_panel": 1, "other_allele_at_site": 7}


def test_query_reads_every_topic_of_the_run(tmp_path: Path) -> None:
    """The regression of Task 34: no unassessed claim on a known call reaches a query."""
    con, run_id = _store(tmp_path, LOCATED)
    every = Filters(
        min_validity=Validity.LIMITED, min_reliability=Level.LOW, gene="PPARG"
    )
    kept = [query(con, run_id, replace(every, topic=t)).findings for t in Topic]
    assert sum(map(len, kept)) == 4


def test_validity_capped_limited(tmp_path: Path) -> None:
    """However strong the study, a single-SNP association is limited."""
    located = [(_association(102, "A", p_value="1E-300"), C_A, _called("0/1", "1"))]
    con, run_id = _store(tmp_path, located)
    (statement,) = _statements(con, run_id).values()
    assert (statement["validity"], statement["source_level"]) == (
        "limited",
        "GWAS Catalog p = 1E-300",
    )
    narrowed = Filters(min_validity=Validity.LIMITED, marker="rs102")
    (finding,) = query(con, run_id, narrowed).findings
    assert finding["claim"] == {
        "schema": "genomelab.claim/1",
        "type": "gwas_association",
        "condition": "type 2 diabetes mellitus",
        "trait_efo": "MONDO_0005148",
        "risk_allele": "A",
        "risk_allele_copies": 1,
        "p_value": "1E-300",
        "pubmed_id": "17463246",
        "study_accession": "GCST000001",
    }
    assert finding["sources"] == [
        {"source": "gwas_catalog", "id": "GCST000001", "version": "20261003"}
    ]


def test_only_a_narrowed_query_returns_an_association(tmp_path: Path) -> None:
    """Off the feed: a query naming no gene, marker or text, the dashboard's and explain's
    included, counts an association under `research_association` before any other filter;
    naming its gene, its marker or a text of it returns it."""
    located = [(_association(102, "A"), C_A, _called("0/1", "1"))]
    con, run_id = _store(tmp_path, located)
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    for f in (Filters(), loosest, replace(loosest, topic=Topic.RISKS)):
        result = query(con, run_id, f)
        assert (result.findings, result.hidden) == ([], {"research_association": 1})
    for narrowed in (
        replace(loosest, gene="pparg"),
        replace(loosest, marker="rs102"),
        replace(loosest, text="diabetes"),
    ):
        assert len(query(con, run_id, narrowed).findings) == 1


def test_no_absolute_effect(tmp_path: Path) -> None:
    located = [
        (_association(102, "A"), C_A, _called("0/1", "1")),
        (
            _association(102, "A", efo_ids=("E",), measure=None, value=None, ci=None),
            C_A,
            _called("0/1", "1"),
        ),
    ]
    con, run_id = _store(tmp_path, located)
    effects = {
        sid: json.loads(s["effect"]) for sid, s in _statements(con, run_id).items()
    }
    assert effects == {
        "gwas:rs102-A:MONDO_0005148": {
            "relative": {"measure": "OR", "value": 1.14, "ci": [1.08, 1.2]},
            "absolute": None,
        },
        "gwas:rs102-A:E": {"relative": None, "absolute": None},
    }


def test_topic_from_efo(tmp_path: Path) -> None:
    """The first of its traits the EFO table maps to a topic other than research, else
    research."""
    located = [
        (
            _association(102, "A", efo_ids=("EFO_0004339", "MONDO_0005148")),
            C_A,
            _called("0/1", "1"),
        ),
        (_association(102, "A", efo_ids=("EFO_0004611",)), C_A, _called("0/1", "1")),
        (_association(102, "A", efo_ids=("EFO_0004339",)), C_A, _called("0/1", "1")),
        (_association(102, "A", efo_ids=(), trait="height"), C_A, _called("0/1", "1")),
    ]
    con, run_id = _store(tmp_path, located)
    topics = {sid: s["topic"] for sid, s in _statements(con, run_id).items()}
    assert topics == {
        "gwas:rs102-A:EFO_0004339+MONDO_0005148": "risks",
        "gwas:rs102-A:EFO_0004611": "wellness",
        "gwas:rs102-A:EFO_0004339": "research",
        "gwas:rs102-A:height": "research",  # no EFO: its id names the trait
    }


def test_genotype_answers_with_the_risk_allele(tmp_path: Path) -> None:
    """`genokeeper genotype` answers a research association by its risk allele's copies:
    a hom-ref user carries a reference risk allele twice, and a hom-ALT one not at
    all, whatever the ALT's dosage says; beside another allele it is not answerable."""
    located = [
        (_association(4, "C"), ("chr1", 4, "C", "A"), _called("0/0", "0")),
        (_association(5, "C"), ("chr1", 5, "C", "A"), _called("1/1", "2")),
        (_association(6, "C"), ("chr1", 6, "C", "A"), _called("0/2", "0", OTHER)),
    ]
    con, run_id = _store(tmp_path, located)
    sites = genotype(con, run_id, gene="PPARG").to_json()["sites"]
    assert [(s["id"], s["answer"], s["not_a_finding_because"]) for s in sites] == [
        ("gwas:rs4-C:MONDO_0005148", "carrier", None),
        (
            "gwas:rs5-C:MONDO_0005148",
            "not_a_carrier",
            "the allele is not carried (dosage 0)",
        ),
    ]  # rs6 is dropped
    site = sites[0] | {"genotype": sites[0]["genotype"] | {"gt": "0/2", "dosage": "0"}}
    assert answer(site) == Answer.UNKNOWN


# --- check_calls on every association ------------------------------------------------

VID = "GRCh38:chr1:102:C:A"
CALLS = {"0/0": "0", "0/1": "1", "1/1": "2", "1/2": "1", "0/2": "0", "./.": "unknown"}
CALLS |= {"0": "0", "1": "1"}


def _check_store(tmp_path: Path, kind: str, gt: str, risk: str, copies: Any) -> Any:
    con = connect(tmp_path / "findings.duckdb")
    run_id = pipeline.start_run(con, {}, {})
    con.execute(
        "INSERT INTO call (run_id, vid, gt, dosage, level, row_type, warnings)"
        " VALUES (?, ?, ?, ?, 'high', 'PASS', [])",
        [run_id, VID, gt, CALLS[gt]],
    )
    claim = {
        "type": "gwas_association",
        "risk_allele": risk,
        "risk_allele_copies": copies,
    }
    con.execute(
        "INSERT INTO statement (run_id, statement_id, kind, vid, claim, validity,"
        " source_level, source, source_version, topic)"
        " VALUES (?, 'gwas:x', ?, ?, ?, 'limited', 'l', 'gwas_catalog', 'v', 'research')",
        [run_id, kind, VID, json.dumps(claim)],
    )
    return con, run_id


@pytest.mark.parametrize(
    ("kind", "gt", "risk", "copies"),
    [
        ("finding", "0/0", "C", 2),
        ("finding", "0/1", "C", 1),
        ("non_finding", "1/1", "C", 0),
        ("finding", "0", "C", 1),
        ("non_finding", "1", "C", 0),
        ("finding", "0/1", "A", 1),
        ("non_finding", "0/2", "A", 0),
        ("unassessed", "./.", "A", None),
    ],
)
def test_check_calls_passes_a_risk_allele_its_call_carries(
    tmp_path: Path, kind: str, gt: str, risk: str, copies: Any
) -> None:
    con, run_id = _check_store(tmp_path, kind, gt, risk, copies)
    pipeline.check_calls(con, run_id)


@pytest.mark.parametrize(
    ("kind", "gt", "risk", "copies"),
    [
        ("finding", "1/1", "C", 0),  # no reference copy
        ("finding", "1/2", "C", 1),  # the ALT and another allele: no reference copy
        ("finding", "0/1", "C", 2),  # the claim's count is not the call's
        ("non_finding", "0/1", "C", 0),
        ("unassessed", "0/0", "C", None),
        ("unassessed", "0/2", "C", None),  # another allele settles the copies
        (
            "finding",
            "0/2",
            "C",
            1,
        ),  # another allele: a reference's copies not countable
        ("non_finding", "1/2", "C", 0),
        ("finding", "0/1", "G", 1),  # neither of the vid's alleles
        ("non_finding", "0/0", "c", 0),
        ("finding", "0/0", "A", 0),
        ("non_finding", "0/1", "A", 0),  # an ALT risk allele faces the same rule
        ("finding", "0/1", "A", 2),
        ("unassessed", "0/2", "A", None),
    ],
)
def test_check_calls_refuses_a_kind_its_call_contradicts(
    tmp_path: Path, kind: str, gt: str, risk: str, copies: Any
) -> None:
    con, run_id = _check_store(tmp_path, kind, gt, risk, copies)
    with pytest.raises(pipeline.CallMismatch):
        pipeline.check_calls(con, run_id)
