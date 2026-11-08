"""The research layer of rare variants in the genes of Genomics England PanelApp's
Intellectual disability panel: the panel file, which variants qualify, what each statement
says, and that a `predicted:` finding is not repeated."""

import json
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.panelapp_mini import GENES, GeneRow, write_panel

from genokeeper.annotate.vep import VepRecord
from genokeeper.config import load_grading
from genokeeper.evidence.clinvar import GENOTYPE_REASON
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.predicted import EscapeRegion
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.sources import Collected, Row, RunContext
from genokeeper.sources.panelapp import (
    ALREADY_PREDICTED,
    FAILED_STATEMENT,
    NOT_A_DIAGNOSIS,
    ONE_COPY_NOTE,
    PANEL_VERSION,
    PANELAPP_TOOL,
    GenePanel,
    PanelAppSource,
    PanelGene,
    collect_panel,
    failed_panel,
    load_panel,
    without_repeats,
)
from genokeeper.sources.predicted import PREDICTED_TOOL, PredictedSource
from genokeeper.sources.research import is_research_layer
from genokeeper.store.db import Kind, Level, RowType

GRADING = load_grading()
LISTS = load_gene_lists()
VID = "GRCh38:chr17:12000:C:T"
GREEN = "GENEGREEN"
AR_TEXT = "BIALLELIC, autosomal or pseudoautosomal"
# As `test_predicted`'s: MANE transcript ENST01, its last exon's coding bases escaping NMD.
REGIONS = [
    EscapeRegion("chr17", 15950, 16000, "+"),
    EscapeRegion("chr17", 20000, 20400, "+"),
]
HET = Called("0/1", "1", RowType.PASS, 40, 30, (15, 15), 0.5, ())
RELIABLE = Reliability(1e-4, Level.HIGH, (), ())


def _panel(tmp_path: Path) -> GenePanel:
    return load_panel(write_panel(tmp_path / "panel.json"))


def _record(**kw: Any) -> VepRecord:
    """A MANE stop-gain that passes every LoF gate, in the panel's green gene."""
    fields: dict[str, Any] = {
        "vid": VID,
        "gene": GREEN,
        "gene_id": "ENSG01",
        "transcript": "ENST01.1",
        "consequence": "stop_gained",
        "mane": True,
        "hgvsc": "ENST01.1:c.898C>T",
        "hgvsp": "ENSP01.1:p.Gln300Ter",
        "exon": "5/9",
        "cds_length": 1500,
        "lof": True,
        "revel": None,
        "alphamissense": None,
        "spliceai": None,
        "splicevault_in_frame": None,
        "pext": 0.9,
    }
    return VepRecord(**(fields | kw))


def _missense(**kw: Any) -> VepRecord:
    return _record(
        **{
            "consequence": "missense_variant",
            "lof": False,
            "revel": GRADING.predicted.revel,
            **kw,
        }
    )


def _ctx(
    tmp_path: Path,
    faf: dict[str, float] | None = None,
    inheritance: dict[str, str] | None = None,
    haploid: bool = False,
) -> RunContext:
    return RunContext(
        run_id="run",
        lists=LISTS,
        grading=GRADING,
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf", faf)),
        ba1_exceptions=frozenset(),
        inheritance=inheritance or {},
        male=True,
        haploid=lambda chrom, pos: haploid,
        source_version="20261011",
    )


def _collect(
    tmp_path: Path,
    *records: VepRecord,
    called: Called = HET,
    faf: dict[str, float] | None = None,
    inheritance: dict[str, str] | None = None,
    genes: Sequence[GeneRow] = GENES,
    haploid: bool = False,
) -> Collected:
    predicted = PredictedSource(
        sites=[(r, called, RELIABLE, 0) for r in records],
        input_rows=len(records),
        escapes={r.transcript.split(".")[0]: REGIONS for r in records},
        truncations=[],
        haploinsufficient=frozenset(),
    )
    panel = load_panel(write_panel(tmp_path / "panel.json", genes))
    return collect_panel(panel, predicted, _ctx(tmp_path, faf, inheritance, haploid))


def _one(tmp_path: Path, *records: VepRecord, **kw: Any) -> dict[str, object]:
    (row,) = _collect(tmp_path, *records, **kw).rows
    return row.statement


def _dropped(collected: Collected) -> dict[str, int]:
    return {reason: n for reason, n in collected.dropped.items() if n}


def test_the_panel_keeps_each_genes_confidence_and_inheritance(tmp_path: Path) -> None:
    """Green, amber and red (`0` and `1` are both red); an empty mode of inheritance reads as not stated, and a gene is found by its symbol or, as a
    symbol can have moved since the panel's HGNC snapshot, by its Ensembl id."""
    panel = _panel(tmp_path)
    assert panel.name == "Intellectual disability"
    assert panel.version == PANEL_VERSION
    assert set(panel.genes) == {g[0] for g in GENES}
    assert panel.genes[GREEN] == PanelGene("green", AR_TEXT)
    assert panel.genes["GENEAMBER"] == PanelGene("amber", "Unknown")
    assert panel.genes["GENERED"] == PanelGene("red", "not stated by the panel")
    assert panel.genes["GENENONE"].confidence == "red"
    renamed = _record(gene="NEWNAME", gene_id="ENSG01")
    assert panel.gene(renamed) == panel.genes[GREEN]
    assert panel.gene(_record(gene="ELSEWHERE", gene_id="ENSG99")) is None


def test_a_missing_panel_names_the_bootstrap(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="scripts/bootstrap.sh"):
        load_panel(tmp_path / "panel.json")


@pytest.mark.parametrize(
    ("fields", "match"),
    [
        ({"id": 286}, "panel 286"),
        ({"version": "11.94"}, "version 11.94"),
        ({"genes": []}, "no genes"),
    ],
)
def test_a_wrong_panel_is_refused(
    tmp_path: Path, fields: dict[str, Any], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        load_panel(write_panel(tmp_path / "panel.json", **fields))


def test_a_duplicate_gene_or_unknown_confidence_is_refused(tmp_path: Path) -> None:
    twice = (*GENES, GENES[0])
    with pytest.raises(ValueError, match="duplicate gene"):
        load_panel(write_panel(tmp_path / "panel.json", twice))
    odd = [(GREEN, "7", AR_TEXT, "ENSG01")]
    with pytest.raises(ValueError, match="confidence"):
        load_panel(write_panel(tmp_path / "panel.json", odd))


def test_a_rare_loss_of_function_in_a_panel_gene_is_one_research_statement(
    tmp_path: Path,
) -> None:
    """Limited, under research, never actionable, predicted-not-classified, with the panel,
    its confidence and inheritance in the claim and the reason; the id carries the vid and
    the transcript, as a variant can have records on two."""
    s = _one(tmp_path, _record())
    assert s["statement_id"] == f"panelapp:{VID}:ENST01"
    assert (s["kind"], s["validity"], s["topic"], s["actionable"], s["effect"]) == (
        Kind.FINDING,
        "limited",
        "research",
        False,
        None,
    )
    assert (s["source"], s["source_id"], s["source_version"]) == (
        PANELAPP_TOOL,
        "ENST01.1",
        PANEL_VERSION,
    )
    assert s["source_level"] == "stop_gained on ENST01.1"
    claim = json.loads(str(s["claim"]))
    assert claim["classification"] == "predicted, not classified"
    assert claim["type"] == "panelapp_variant"
    assert claim["panel"] == "Intellectual disability"
    assert claim["panel_version"] == PANEL_VERSION
    assert claim["panel_confidence"] == "green"
    assert claim["panel_inheritance"] == AR_TEXT
    reason = str(s["reason"])
    assert "predicted, not classified" in reason
    assert "Intellectual disability" in reason
    assert (
        "gene rated green on the panel: its link to the condition is well established;"
        " this says nothing about this variant"
    ) in reason
    assert "diagnostic grade" not in reason
    assert AR_TEXT in reason
    assert "discuss with a geneticist" in reason
    assert NOT_A_DIAGNOSIS in reason
    assert ONE_COPY_NOTE in reason


@pytest.mark.parametrize(
    ("gene", "level"),
    [("GENEAMBER", "amber"), ("GENERED", "red")],
)
def test_amber_and_red_genes_name_the_genes_evidence(
    tmp_path: Path, gene: str, level: str
) -> None:
    reason = str(_one(tmp_path, _record(gene=gene, gene_id="ENSG02"))["reason"])
    assert f"gene rated {level} on the panel: the evidence for its link" in reason
    assert "this says nothing about this variant" in reason
    assert "diagnostic grade" not in reason


MONO = "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted"
XLR = "X-LINKED: hemizygous mutation in males, biallelic mutations in females"
BOTH = "BOTH monoallelic and biallelic, autosomal or pseudoautosomal"


@pytest.mark.parametrize("inheritance", [MONO, XLR, BOTH, "Unknown", ""])
def test_the_one_copy_note_is_only_for_a_biallelic_gene(
    tmp_path: Path, inheritance: str
) -> None:
    """Every statement says it is no diagnosis; the one-copy sentence only where it is true."""
    genes = [(GREEN, "3", inheritance, "ENSG01")]
    reason = str(_one(tmp_path, _record(), genes=genes)["reason"])
    assert NOT_A_DIAGNOSIS in reason
    assert ONE_COPY_NOTE not in reason


@pytest.mark.parametrize(
    ("called", "haploid"),
    [
        (Called("1/1", "2", RowType.PASS, 40, 30, (0, 30), 1.0, ()), False),
        (HET, True),
        (Called(None, "unknown", RowType.PASS, 40, 30, None, None, ("x",)), False),
    ],
)
def test_the_one_copy_note_is_only_for_a_heterozygous_call(
    tmp_path: Path, called: Called, haploid: bool
) -> None:
    """A homozygote, a hemizygous male (a haploid site) and an unknown genotype are no het:
    the sentence is dropped, in a gene the panel calls biallelic."""
    reason = str(_one(tmp_path, _record(), called=called, haploid=haploid)["reason"])
    assert NOT_A_DIAGNOSIS in reason
    assert ONE_COPY_NOTE not in reason


def test_a_missense_qualifies_at_its_calibrated_score(tmp_path: Path) -> None:
    """At REVEL's threshold it qualifies; one point under, no; AlphaMissense only where
    REVEL has no score: one calibrated tool per variant."""
    t = GRADING.predicted
    assert _one(tmp_path, _missense())["source_level"] == f"missense, REVEL {t.revel:g}"
    alpha = _missense(revel=None, alphamissense=t.alphamissense)
    assert (
        _one(tmp_path, alpha)["source_level"]
        == f"missense, AlphaMissense {t.alphamissense:g}"
    )
    below = _missense(revel=t.revel - 0.01, alphamissense=0.99)
    collected = _collect(tmp_path, below)
    assert (collected.rows, _dropped(collected)) == ([], {"not_a_panel_rule": 1})


@pytest.mark.parametrize(
    ("record", "reason"),
    [
        (_record(gene="ELSEWHERE", gene_id="ENSG99"), "not_in_panel"),
        (_record(consequence="start_lost", lof=False), "not_a_panel_rule"),
        (
            _record(consequence="intron_variant", lof=False, spliceai=0.9),
            "not_a_panel_rule",
        ),
        (_record(mane=False), "failed_lof_gate"),
        (_record(pext=0.01), "failed_lof_gate"),
    ],
)
def test_a_record_that_does_not_qualify_is_counted(
    tmp_path: Path, record: VepRecord, reason: str
) -> None:
    collected = _collect(tmp_path, record)
    assert (collected.rows, collected.input_rows) == ([], 1)
    assert _dropped(collected) == {reason: 1}


def test_a_common_variant_or_no_copy_is_counted_out(tmp_path: Path) -> None:
    """Rare by the rule the predicted findings use (`grade_predicted`'s BA1 and BS1)."""
    common = _collect(tmp_path, _record(), faf={VID: 0.5})
    assert _dropped(common) == {"population_frequency": 1}
    absent = Called("0/0", "0", RowType.PASS, 40, 30, (30, 0), 0.0, ())
    assert _dropped(_collect(tmp_path, _record(), called=absent)) == {"not_carried": 1}


@pytest.mark.parametrize("mode", ["AR", "AD", None])
@pytest.mark.parametrize(("faf", "kept"), [(0.0099, 1), (0.01, 0), (0.02, 0)])
def test_rare_is_one_limit_at_any_inheritance(
    tmp_path: Path, mode: str | None, faf: float, kept: int
) -> None:
    """Below 0.01 whatever the gene's inheritance: `grade_predicted` lets a recessive gene's
    allele through up to BA1's 0.05, which is not rare."""
    inheritance = {GREEN: mode} if mode else {}
    collected = _collect(tmp_path, _record(), faf={VID: faf}, inheritance=inheritance)
    assert len(collected.rows) == kept
    assert _dropped(collected) == ({} if kept else {"population_frequency": 1})


def test_a_genes_own_limit_still_applies(tmp_path: Path) -> None:
    """Under 0.01 is not enough: a gene's `[bs1_by_gene]` limit (SMAD4's 1e-06) applies too,
    through the frequency rule `grade_predicted` uses."""
    assert GRADING.bs1_by_gene["SMAD4"] == 1e-06
    genes = [("SMAD4", "3", MONO, "ENSG05")]
    record = _record(gene="SMAD4", gene_id="ENSG05")
    collected = _collect(tmp_path, record, faf={VID: 5e-06}, genes=genes)
    assert collected.rows == []
    assert _dropped(collected) == {"population_frequency": 1}


def test_an_unknown_genotype_is_unassessed(tmp_path: Path) -> None:
    unknown = Called(None, "unknown", RowType.PASS, 40, 30, None, None, ("lowq",))
    s = _one(tmp_path, _record(), called=unknown)
    assert s["kind"] == Kind.UNASSESSED
    assert GENOTYPE_REASON in str(s["reason"])


def test_the_splice_loss_says_what_spliceai_scored(tmp_path: Path) -> None:
    """A canonical splice LoF in a gene SpliceAI was not run on says so; one in a listed
    gene, which SpliceAI did score, does not."""
    splice = {
        "consequence": "splice_donor_variant",
        "hgvsp": None,
        "splicevault_in_frame": False,
    }
    outside = str(_one(tmp_path, _record(**splice))["reason"])
    assert "SpliceAI" in outside and "did not score" in outside
    listed = _record(**splice, gene="BRCA1", gene_id="ENSG00000012048", spliceai=0.9)
    assert "SpliceAI" not in str(_one(tmp_path, listed)["reason"])
    assert "SpliceAI" not in str(_one(tmp_path, _record())["reason"])


def test_two_transcripts_of_one_variant_are_two_statements(tmp_path: Path) -> None:
    """A variant can have records on two transcripts: the transcript keeps their ids apart,
    as a shared id would refuse the whole run at the store's primary key."""
    two = _collect(
        tmp_path,
        _record(transcript="ENST01.1"),
        _record(transcript="ENST02.1"),
    )
    ids = [row.statement["statement_id"] for row in two.rows]
    assert ids == [f"panelapp:{VID}:ENST01", f"panelapp:{VID}:ENST02"]


def test_accounting_balances(tmp_path: Path) -> None:
    collected = _collect(
        tmp_path,
        _record(),
        _record(gene="ELSEWHERE", gene_id="ENSG99", vid="GRCh38:chr17:12001:C:T"),
    )
    assert collected.input_rows == 2
    assert len(collected.rows) + sum(collected.dropped.values()) == 2


def _panel_row(vid: str) -> Row:
    return Row(
        {"vid": vid},
        {"vid": vid},
        {"statement_id": f"panelapp:{vid}:ENST01", "vid": vid, "kind": Kind.FINDING},
    )


def _predicted_row(vid: str, kind: Kind) -> Row:
    return Row(
        None,
        None,
        {"statement_id": f"predicted:{vid}:ENST01", "vid": vid, "kind": kind},
    )


def _collected(tool: str, rows: list[Row]) -> Collected:
    return Collected(tool, rows, {}, len(rows), {ALREADY_PREDICTED: 0})


def test_a_predicted_finding_is_not_repeated() -> None:
    """A variant a `predicted:` statement states, on any of its transcripts, keeps that one;
    a predicted statement kept only as a non-finding (one recessive copy) states nothing to
    the user, so the panel's is kept; other tools' rows are left alone."""
    panel = _collected(
        PANELAPP_TOOL,
        [_panel_row("a"), _panel_row("b"), _panel_row("c"), _panel_row("d")],
    )
    predicted = _collected(
        PREDICTED_TOOL,
        [
            _predicted_row("a", Kind.FINDING),
            _predicted_row("b", Kind.UNASSESSED),
            _predicted_row("c", Kind.NON_FINDING),
        ],
    )
    other = _collected("clinvar", [_panel_row("a")])
    kept_panel, kept_predicted, kept_other = without_repeats([panel, predicted, other])
    assert [r.statement["vid"] for r in kept_panel.rows] == ["c", "d"]
    assert kept_panel.dropped == {ALREADY_PREDICTED: 2}
    assert kept_panel.input_rows == 4
    assert kept_predicted == predicted
    assert kept_other == other


def test_a_compound_het_the_promotion_made_is_not_repeated() -> None:
    """The promotion runs before: a one-copy `predicted:` non-finding it makes a finding
    stands, and the panel's statement of that variant goes."""
    promoted = _predicted_row("c", Kind.FINDING)
    panel = _collected(PANELAPP_TOOL, [_panel_row("c")])
    kept, _ = without_repeats([panel, _collected(PREDICTED_TOOL, [promoted])])
    assert (kept.rows, kept.dropped) == ([], {ALREADY_PREDICTED: 1})


def test_a_failed_layer_is_one_unassessed_statement(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    pending = replace(ctx, run_id="pending")
    collected = PanelAppSource(failed_panel(pending, "no panel file")).collect(ctx)
    (row,) = collected.rows
    assert row.statement["run_id"] == "run"
    s = row.statement
    assert (row.variant, row.call) == (None, None)
    assert (s["statement_id"], s["kind"], s["vid"], s["gene"]) == (
        FAILED_STATEMENT,
        Kind.UNASSESSED,
        None,
        None,
    )
    assert (s["validity"], s["topic"], s["actionable"], s["source"]) == (
        "limited",
        "research",
        False,
        PANELAPP_TOOL,
    )
    assert "no panel file" in str(s["reason"])
    assert collected.input_rows == 1


def test_a_panelapp_row_is_a_research_layer() -> None:
    """Hidden, and counted on the Changes tab, as every research layer is."""
    assert is_research_layer(PANELAPP_TOOL)
