"""Predicted findings for variants ClinVar has not classified (track design, slice 1e): the
loss-of-function rule with its six gates and the carrier rule, missense and splice scores,
frequency, ClinVar's own alleles left out, phasing of two carrier LoFs, and the feed order."""

import json
from dataclasses import replace
from pathlib import Path

import duckdb
import pytest
from fixtures.clinvar_mini import write_clinvar
from fixtures.faf_mini import write_faf_table
from fixtures.genome_mini import BRCA1_NMD_ROWS
from fixtures.gvcf_mini import (
    HOMOPOLYMER,
    other_base,
    ref,
    variant,
    write_fasta,
    write_gvcf,
)
from fixtures.store_mini import Row as StoreRow
from fixtures.store_mini import seed

from genokeeper.annotate.vep import VepRecord, vep_input
from genokeeper.config import load_grading
from genokeeper.evidence.clinvar import ValidityGrade
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.predicted import (
    BOTH_COPIES,
    GATE_HOMOPOLYMER,
    GATE_MAPPABILITY,
    GATE_NMD,
    GATE_NO_NMD,
    GATE_NO_PEXT,
    GATE_NO_SPLICEVAULT,
    GATE_NO_STOP,
    GATE_PEXT,
    GATE_SPLICEVAULT,
    GATE_START_LOSS,
    GATE_TRANSCRIPT,
    GENETICIST_NOTE,
    NOT_CARRIED,
    ONE_RECESSIVE_COPY,
    PREDICTED_REASON,
    START_PROXIMAL,
    STOP_CONSEQUENCES,
    TRACK,
    Escape,
    EscapeRegion,
    Site,
    escape,
    grade_predicted,
    homopolymer_runs,
    load_escapes,
    truncations_3prime,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import clinvar_vids
from genokeeper.phasing import IN_CIS_REASON, IN_TRANS_REASON, Phase
from genokeeper.query import Filters, query
from genokeeper.sources import Collected, Row, RunContext, promote
from genokeeper.sources.predicted import DROP_REASONS, PredictedSource
from genokeeper.store.db import Kind, Level, RowType, Topic, Validity, connect

LISTS = load_gene_lists()
GRADING = load_grading()
# The lists without MUTYH's Tier 4 row: an ACMG SF recessive gene outside the carrier lists,
# the case each of ACMG SF's five AR genes outside Tier 3 was before Tier 4 (all five are
# Tier 4 genes now).
OUTSIDE_CARRIER = {
    **LISTS,
    "carrier_tier4": replace(
        LISTS["carrier_tier4"],
        genes={g: m for g, m in LISTS["carrier_tier4"].genes.items() if g != "MUTYH"},
    ),
}
HAPLOINSUFFICIENT = frozenset({"BRCA1"})
# MANE transcript ENST01 on chr17's + strand, a CDS of 1500 bases (its stop included):
# UCSC's escape regions are its last exon's 400 coding bases and the 50 before the last
# junction (rule 1, CDS 1051-1500).
REGIONS = [
    EscapeRegion("chr17", 15950, 16000, "+"),
    EscapeRegion("chr17", 20000, 20400, "+"),
]


def _record(**kw: object) -> VepRecord:
    """A MANE stop-gain at codon 300 of BRCA1's ENST01, expressed, in no escape region."""
    fields: dict[str, object] = {
        "vid": "GRCh38:chr17:12000:C:T",
        "gene": "BRCA1",
        "gene_id": "ENSG01",
        "transcript": "ENST01",
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
    return VepRecord(**(fields | kw))  # type: ignore[arg-type]


def _site(record: VepRecord, regions: object = REGIONS, **kw: object) -> Site:
    stops = record.terms & STOP_CONSEQUENCES
    found = escape(record, regions) if stops else None  # type: ignore[arg-type]
    fields: dict[str, object] = {
        "record": record,
        "dosage": "1",
        "faf": None,
        "homopolymer": 0,
        "flags": (),
        "escape": found,
        "truncations_3prime": 0,
    }
    return Site(**(fields | kw))  # type: ignore[arg-type]


def _grade(site: Site, **kw: object) -> ValidityGrade:
    return grade_predicted(
        site,
        haploinsufficient=HAPLOINSUFFICIENT,
        lists=LISTS,
        grading=GRADING,
        **kw,  # type: ignore[arg-type]
    )


def test_lof_in_haploinsufficient_listed_gene() -> None:
    """A predicted LoF in BRCA1 (ACMG SF, ClinGen dosage score 3) that passes every gate is
    a finding under risks, validity limited, "predicted, not classified"."""
    grade = _grade(_site(_record()))
    assert (grade.kind, grade.level, grade.topic) == (
        Kind.FINDING,
        Validity.LIMITED,
        Topic.RISKS,
    )
    assert grade.reason is not None and grade.reason.startswith(PREDICTED_REASON)
    assert grade.source_level == "stop_gained on ENST01"


def test_lof_in_tolerant_gene_not_finding() -> None:
    """A predicted LoF in a listed gene ClinGen does not rate haploinsufficient (MYH7, whose
    disease is not loss of function), and in an unlisted gene, is no statement."""
    for gene in ("MYH7", "GENEZ"):
        grade = _grade(_site(_record(gene=gene)))
        assert grade.kind is Kind.NON_FINDING, gene
        assert grade.topic is None


def test_homozygous_lof_in_recessive_gene_is_research() -> None:
    """Both copies lost in a recessive gene neither haploinsufficient nor a carrier gene
    (GENEZ, AR by CGD) go under research, the reason naming both copies; a dominant gene's
    (MYH7) stays no statement. MUTYH, AR in ACMG SF, is a Tier 4 carrier gene: both copies
    there are a finding under risks, as in a Tier 3 gene."""
    grade = _grade(_site(_record(gene="GENEZ"), dosage="2"), inheritance="AR")
    assert (grade.kind, grade.level, grade.topic) == (
        Kind.FINDING,
        Validity.LIMITED,
        Topic.RESEARCH,
    )
    assert grade.reason is not None and BOTH_COPIES in grade.reason
    mutyh = _grade(_site(_record(gene="MUTYH"), dosage="2"))
    assert (mutyh.kind, mutyh.topic) == (Kind.FINDING, Topic.RISKS)
    assert _grade(_site(_record(gene="MYH7"), dosage="2")).kind is Kind.NON_FINDING
    # An unknown genotype could be both copies: unassessed (stage-1 amendment 5).
    unknown = _grade(_site(_record(gene="GENEZ"), dosage="unknown"), inheritance="AR")
    assert (unknown.kind, unknown.topic) == (Kind.UNASSESSED, Topic.RESEARCH)


def test_lof_in_carrier_gene_is_carrier_finding() -> None:
    """In a carrier Tier 3 gene a predicted LoF is a carrier finding without the
    haploinsufficiency gate (design decision, 2026-10-21); two copies are under risks."""
    het = _grade(_site(_record(gene="CFTR")))
    assert (het.kind, het.level, het.topic) == (
        Kind.FINDING,
        Validity.LIMITED,
        Topic.CARRIER,
    )
    hom = _grade(_site(_record(gene="CFTR"), dosage="2"))
    assert (hom.kind, hom.topic) == (Kind.FINDING, Topic.RISKS)


def test_tier4_predicted_lof_carrier() -> None:
    """A predicted stop-gain in a carrier Tier 4 gene outside Tier 3 (ABCA4, AR Definitive
    in ClinGen), absent from gnomAD, is a carrier finding as a Tier 3 gene's is; before
    Tier 4 it was one copy in a recessive gene, a non-finding."""
    grade = _grade(_site(_record(gene="ABCA4"), faf=None))
    assert (grade.kind, grade.level, grade.topic) == (
        Kind.FINDING,
        Validity.LIMITED,
        Topic.CARRIER,
    )
    assert grade.reason is not None and grade.reason.startswith(PREDICTED_REASON)
    # SYNE1 is Tier 4 but AD in CGD (`gene_mode`): the gate stays, no statement.
    assert "SYNE1" in LISTS["carrier_tier4"].genes
    dominant = _grade(_site(_record(gene="SYNE1")), inheritance="AD")
    assert (dominant.kind, dominant.topic) == (Kind.NON_FINDING, None)


STOP_IN_LAST_EXON = {"hgvsp": "ENSP01.1:p.Gln450Ter"}  # CDS 1348: 20248
STOP_NEAR_JUNCTION = {"hgvsp": "ENSP01.1:p.Gln355Ter"}  # CDS 1063: 15963
STOP_EARLY = {"hgvsp": "ENSP01.1:p.Gln20Ter"}  # CDS 58
EARLY = f"{GATE_NMD} ({START_PROXIMAL})"
# A single coding exon (UCSC's rule 2), as GJB2's: its region is the whole CDS.
SINGLE_EXON = [EscapeRegion("chr17", 10000, 11500, "+")]


@pytest.mark.parametrize(
    ("record", "site", "gate"),
    [
        (
            {"vid": "GRCh38:chr17:12000:CA:C", "consequence": "frameshift_variant"},
            {"homopolymer": 6},
            GATE_HOMOPOLYMER,
        ),
        ({}, {"flags": ("low_mappability",)}, GATE_MAPPABILITY),
        ({"mane": False}, {}, GATE_TRANSCRIPT),
        (STOP_IN_LAST_EXON, {}, GATE_NMD),
        (STOP_NEAR_JUNCTION, {}, GATE_NMD),
        (STOP_EARLY, {}, GATE_NMD),
        (
            {
                "consequence": "frameshift_variant",
                "hgvsp": "ENSP01.1:p.Lys300AsnfsTer150",  # the stop at codon 449
            },
            {},
            GATE_NMD,
        ),
        (STOP_IN_LAST_EXON, {"truncations_3prime": 1}, GATE_NMD),
        # Never rescued: a stop near the start may re-initiate translation.
        (STOP_EARLY, {"truncations_3prime": 2}, GATE_NMD),
        # Never rescued either where another rule matches too: a single coding exon, VEP's
        # last exon.
        (STOP_EARLY, {"regions": SINGLE_EXON, "truncations_3prime": 2}, EARLY),
        ({"exon": "9/9", **STOP_EARLY}, {"truncations_3prime": 2}, EARLY),
        # VEP's last exon (9/9) though the track does not place the stop there.
        ({"exon": "9/9"}, {}, GATE_NMD),
        # A transcript the track does not list fails, whatever VEP says (fail closed).
        ({"exon": "9/9", **STOP_IN_LAST_EXON}, {"regions": []}, GATE_NO_NMD),
        (
            {
                "consequence": "frameshift_variant",
                "hgvsp": "ENSP01.1:p.Lys300AsnfsTer?",
            },
            {},
            GATE_NO_STOP,
        ),
        (
            {"consequence": "start_lost", "lof": False, "hgvsp": None},
            {},
            GATE_START_LOSS,
        ),
        ({"pext": 0.3}, {}, GATE_PEXT),
        ({"pext": None}, {}, GATE_NO_PEXT),
        (
            {"consequence": "splice_donor_variant", "splicevault_in_frame": True},
            {},
            GATE_SPLICEVAULT,
        ),
        ({"consequence": "splice_donor_variant"}, {}, GATE_NO_SPLICEVAULT),
    ],
)
def test_lof_gates(
    record: dict[str, object], site: dict[str, object], gate: str
) -> None:
    """Each gate of the LoF filter (design decision, 2026-10-23), failed alone, puts the
    variant under research with the gate named, still validity limited."""
    grade = _grade(_site(_record(**record), **site))
    assert (grade.kind, grade.level, grade.topic) == (
        Kind.FINDING,
        Validity.LIMITED,
        Topic.RESEARCH,
    )
    assert grade.reason is not None and f"failed LoF gate {gate}" in grade.reason


def test_escaping_stop_rescued_by_two_downstream_clinvar_truncations() -> None:
    """A stop in the last exon passes gate 3 when ClinVar holds two P/LP truncations 3' of
    it on the transcript (the decision tree's critical region); ones 5' of it, or past the
    CDS end, do not count. A frameshift's new stop is placed from its fsTer."""
    record = _record(**STOP_IN_LAST_EXON)
    found = escape(record, REGIONS)
    assert found is not None
    assert (found.rule, found.stop) == (TRACK, ("chr17", 20248))
    truncations = [
        ("BRCA1", "chr17", 20300),
        ("BRCA1", "chr17", 20350),
        ("BRCA1", "chr17", 20100),  # 5' of the stop
        ("BRCA1", "chr17", 20500),  # past the CDS end
        ("BRCA2", "chr17", 20310),  # another gene
    ]
    count = truncations_3prime(found, REGIONS, "BRCA1", truncations)
    assert count == 2
    grade = _grade(_site(record, truncations_3prime=count))
    assert (grade.kind, grade.topic) == (Kind.FINDING, Topic.RISKS)
    # A GJB2-like single coding exon: a stop past base 100 is rescued the same way.
    gjb2 = _record(gene="GJB2", **STOP_IN_LAST_EXON)
    found = escape(gjb2, SINGLE_EXON)
    assert found is not None and found.rule == TRACK
    downstream = [("GJB2", "chr17", 11400), ("GJB2", "chr17", 11450)]
    count = truncations_3prime(found, SINGLE_EXON, "GJB2", downstream)
    grade = _grade(_site(gjb2, SINGLE_EXON, truncations_3prime=count))
    assert (count, grade.kind, grade.topic) == (2, Kind.FINDING, Topic.CARRIER)
    frameshift = _record(
        consequence="frameshift_variant", hgvsp="ENSP01.1:p.Lys300AsnfsTer12"
    )
    assert escape(frameshift, REGIONS) is None  # the stop at codon 311, CDS 931
    minus = [replace(r, strand="-") for r in REGIONS[:2]]
    stop = escape(_record(**STOP_IN_LAST_EXON), minus)
    assert stop is not None and stop.stop == ("chr17", 20400 - 297)
    # The cross-check: a stop-gain whose CDS length VEP does not give, in a track region.
    inside = escape(_record(vid="GRCh38:chr17:20100:C:T", cds_length=None), REGIONS)
    assert inside == Escape(TRACK, ("chr17", 20100))


def test_rescue_counts_only_the_stop_transcripts_coding_bases() -> None:
    """BR-014: a stop near the last junction counts only ClinVar truncations in its
    transcript's coding bases 3' of it, the track's regions: two in the last intron (on an
    exon only another transcript has) rescue nothing, two in the regions do, on either
    strand."""
    record = _record(**STOP_NEAR_JUNCTION)
    found = escape(record, REGIONS)
    assert found is not None and found.stop == ("chr17", 15963)
    intron = [("BRCA1", "chr17", 16001), ("BRCA1", "chr17", 18000)]
    count = truncations_3prime(found, REGIONS, "BRCA1", intron)
    grade = _grade(_site(record, truncations_3prime=count))
    assert (count, grade.topic) == (0, Topic.RESEARCH)
    exonic = [("BRCA1", "chr17", 16000), ("BRCA1", "chr17", 20001)]
    assert truncations_3prime(found, REGIONS, "BRCA1", exonic) == 2
    minus = [replace(r, strand="-") for r in REGIONS]
    stop = escape(_record(**STOP_IN_LAST_EXON), minus)
    assert stop is not None and stop.stop == ("chr17", 20103)
    downstream = [("BRCA1", "chr17", 20050), ("BRCA1", "chr17", 18000)]
    downstream.append(("BRCA1", "chr17", 15990))
    assert truncations_3prime(stop, minus, "BRCA1", downstream) == 2


@pytest.mark.parametrize(
    ("truncations", "topic"),
    [
        ([("BRCA1", "chr17", 16000), ("BRCA1", "chr17", 20001)], Topic.RISKS),
        ([("BRCA1", "chr17", 16001), ("BRCA1", "chr17", 18000)], Topic.RESEARCH),
    ],
)
def test_collect_rescues_on_the_stop_transcripts_regions(
    tmp_path: Path, truncations: list[tuple[str, str, int]], topic: Topic
) -> None:
    """BR-014 through `collect`: a stop near the last junction is rescued by two ClinVar
    truncations in its own transcript's escape regions, not by two in the last intron,
    though another transcript's regions cover them."""
    record = _record(**STOP_NEAR_JUNCTION)
    source = PredictedSource(
        sites=[(record, CALLED, RELIABLE, 0)],
        input_rows=1,
        escapes={
            "ENST01": REGIONS,
            "ENST99": [EscapeRegion("chr17", 15000, 21000, "+")],
        },
        truncations=truncations,
        haploinsufficient=HAPLOINSUFFICIENT,
    )
    (row,) = source.collect(_context(tmp_path)).rows
    assert row.statement["topic"] == topic


def test_escape_regions_read_in_the_track_layout(tmp_path: Path) -> None:
    """UCSC's nmdEscMane rows as bigBedToBed writes them (BRCA1's, read 2026-10-23 through
    the pinned image): the rule in the mouseover column, one versioned transcript in the
    next; rules 1 and 2 kept per unversioned transcript, 3 and 4 not."""
    bed = tmp_path / "nmd.bed"
    bed.write_text("".join(f"{row}\n" for row in BRCA1_NMD_ROWS))
    assert load_escapes(bed) == {
        "ENST00000357654": [
            EscapeRegion("chr17", 43045677, 43045802, "-"),
            EscapeRegion("chr17", 43047642, 43047692, "-"),
        ]
    }


def test_homopolymer_read_from_the_reference(tmp_path: Path) -> None:
    """An indel of one base in the reference's run of eight A's sits in a homopolymer of
    eight; an SNV and an indel of mixed bases in none."""
    fasta = write_fasta(tmp_path / "ref.fa")
    contig, pos, _ = HOMOPOLYMER
    deletion = f"GRCh38:{contig}:{pos}:CA:C"
    mixed = f"GRCh38:{contig}:{pos}:CAA:CTA"
    snv = f"GRCh38:{contig}:{pos}:C:T"
    assert homopolymer_runs(fasta, [deletion, mixed, snv]) == {
        deletion: 8,
        mixed: 0,
        snv: 0,
    }


def test_predicted_finding_carries_geneticist_note() -> None:
    """Every predicted finding, dominant, carrier or research, says to discuss it with a
    geneticist and that it would need confirmation by a clinical lab."""
    for grade in (
        _grade(_site(_record())),
        _grade(_site(_record(gene="CFTR"))),
        _grade(_site(_record(mane=False))),
    ):
        assert grade.reason is not None
        assert PREDICTED_REASON in grade.reason and GENETICIST_NOTE in grade.reason
        assert "would need confirmation by a clinical lab" in grade.reason


def _missense(**kw: object) -> VepRecord:
    return _record(
        consequence="missense_variant",
        lof=False,
        hgvsp="ENSP01.1:p.Arg300His",
        pext=None,
        **kw,
    )


def test_missense_above_revel_threshold() -> None:
    """A missense in a listed gene at or above REVEL 0.773 is a finding (AlphaMissense
    0.906 only where REVEL has no score); below, a non-finding; unlisted, research; a
    splice change at SpliceAI 0.2 likewise."""
    above = _grade(_site(_missense(revel=0.8)))
    assert (above.kind, above.topic, above.source_level) == (
        Kind.FINDING,
        Topic.RISKS,
        "missense, REVEL 0.8",
    )
    assert _grade(_site(_missense(revel=0.7))).kind is Kind.NON_FINDING
    # REVEL decides where it has a score, whatever AlphaMissense says.
    assert _grade(_site(_missense(revel=0.5, alphamissense=0.99))).kind is (
        Kind.NON_FINDING
    )
    alone = _grade(_site(_missense(alphamissense=0.95)))
    assert (alone.kind, alone.source_level) == (
        Kind.FINDING,
        "missense, AlphaMissense 0.95",
    )
    elsewhere = _grade(_site(_missense(revel=0.9, gene="GENEZ")))
    assert (elsewhere.kind, elsewhere.topic) == (Kind.FINDING, Topic.RESEARCH)
    splice = _grade(_site(_missense(revel=0.1, spliceai=0.2)))
    assert (splice.kind, splice.source_level) == (Kind.FINDING, "splice, SpliceAI 0.2")
    # One copy in an ACMG SF recessive gene outside the carrier lists is no statement; in
    # MUTYH, a Tier 4 gene, it is carrier status.
    outside = grade_predicted(
        _site(_missense(revel=0.9, gene="MUTYH")),
        haploinsufficient=HAPLOINSUFFICIENT,
        lists=OUTSIDE_CARRIER,
        grading=GRADING,
    )
    assert outside.kind is Kind.NON_FINDING
    carrier = _grade(_site(_missense(revel=0.9, gene="MUTYH")))
    assert (carrier.kind, carrier.topic) == (Kind.FINDING, Topic.CARRIER)


def test_research_honours_the_dosage() -> None:
    """An unknown genotype is no finding (stage-1 amendment 5): research it would be is
    unassessed, keeping its topic. No copy is no statement, in a listed gene too, and not
    "one copy in a recessive gene"."""
    for record in (_missense(revel=0.95, gene="GENEZ"), _record(mane=False)):
        unknown = _grade(_site(record, dosage="unknown"))
        assert (unknown.kind, unknown.topic) == (Kind.UNASSESSED, Topic.RESEARCH)
    for record in (
        _missense(revel=0.95, gene="GENEZ"),
        _missense(revel=0.9, gene="CFTR"),
    ):
        none = _grade(_site(record, dosage="0"))
        assert (none.kind, none.reason) == (Kind.NON_FINDING, NOT_CARRIED)
    assert DROP_REASONS[NOT_CARRIED] == "not_carried"


def test_common_variant_excluded() -> None:
    """A prediction must be absent from gnomAD or at or below its gene's BS1 (BRCA1's
    0.0001) and BA1; a recessive gene outside [bs1_by_gene] has only BA1."""
    assert _grade(_site(_record(), faf=0.00005)).kind is Kind.FINDING
    assert _grade(_site(_record(), faf=0.0002)).kind is Kind.NON_FINDING
    assert _grade(_site(_record(gene="CFTR"), faf=0.02)).kind is Kind.FINDING
    assert _grade(_site(_record(gene="CFTR"), faf=0.06)).kind is Kind.NON_FINDING


def test_vep_is_given_every_pass_allele(tmp_path: Path) -> None:
    """VEP is given every PASS allele the sample carries, normalised, whatever ClinVar
    says: a benign record's, a multiallelic record's and one with no germline
    classification (a somatic-only record, no CLNSIG); `clinvar_vids` still names the
    classified ones, which grading filters."""
    fasta = write_fasta(tmp_path / "ref.fa")
    snv = ("chr1", 200)
    benign = ("chr1", 300)
    multi = f"{other_base(ref('chr1', 500))},{other_base(ref('chr1', 500), 1)}"
    calls = [
        variant(*snv, 1, [other_base(ref(*snv))], gt="0/1", gq=40, ad=[10, 10]),
        variant(*benign, 1, [other_base(ref(*benign))], gt="1/1", gq=40, ad=[0, 20]),
        variant(
            "chr1", 400, 1, [other_base(ref("chr1", 400))], gt="0/0", gq=30,
            ad=[20, 2], filt="RefCall",
        ),
        variant("chr1", 500, 1, [multi.split(",")[0]], gt="0/1", gq=40, ad=[10, 10]),
    ]  # fmt: skip
    vcf = write_gvcf(tmp_path / "sample.vcf.gz", calls)
    clinvar = write_clinvar(
        tmp_path / "clinvar.vcf.gz",
        "GRCh38",
        [
            (
                "1",
                benign[1],
                9001,
                ref(*benign),
                other_base(ref(*benign)),
                {
                    "CLNSIG": "Benign",
                    "CLNREVSTAT": "criteria_provided,_single_submitter",
                },
            ),
            (
                "1",
                snv[1],
                9002,
                ref(*snv),
                other_base(ref(*snv)),
                {"GENEINFO": "TP53:1"},
            ),
            # Multiallelic: no claim, but classified.
            ("1", 500, 9003, ref("chr1", 500), multi, {"CLNSIG": "Pathogenic"}),
        ],
    )
    classified = clinvar_vids(clinvar, fasta=fasta)
    out = tmp_path / "work" / "input.vcf"
    out.parent.mkdir()
    given = vep_input(vcf, out, fasta=fasta)
    rows = [line.split("\t") for line in out.read_text().splitlines() if line[0] != "#"]
    snv_vid = f"GRCh38:chr1:200:{ref(*snv)}:{other_base(ref(*snv))}"
    benign_vid = f"GRCh38:chr1:300:{ref(*benign)}:{other_base(ref(*benign))}"
    multi_vid = f"GRCh38:chr1:500:{ref('chr1', 500)}:{multi.split(',')[0]}"
    assert [row[2] for row in rows] == given == [snv_vid, benign_vid, multi_vid]
    assert classified >= {benign_vid, multi_vid} and snv_vid not in classified
    assert not list(out.parent.glob("*.pass.bcf")) and not list(
        out.parent.glob("*.norm*")
    )


def _context(tmp_path: Path) -> RunContext:
    return RunContext(
        run_id="run",
        lists=LISTS,
        grading=GRADING,
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=False,
        haploid=lambda _c, _p: False,
        source_version="20261012",
    )


CALLED = Called("0/1", "1", RowType.PASS, 40, 30, (15, 15), 0.5, ())
RELIABLE = Reliability(1e-4, Level.HIGH, (), ())


def test_two_carrier_lofs_are_phased(tmp_path: Path) -> None:
    """Two predicted LoFs in one carrier gene are a candidate pair for read-backed
    phasing; in trans they become a finding under risks, in cis they stay carrier
    findings. The accounting counts every VEP input variant."""
    first = _record(gene="CFTR", vid="GRCh38:chr7:1000:C:T")
    second = _record(gene="CFTR", vid="GRCh38:chr7:2000:G:A")
    source = PredictedSource(
        sites=[(first, CALLED, RELIABLE, 0), (second, CALLED, RELIABLE, 0)],
        input_rows=5,
        escapes={"ENST01": REGIONS},
        truncations=[],
        haploinsufficient=HAPLOINSUFFICIENT,
    )
    ctx = _context(tmp_path)
    collected = [source.collect(ctx)]
    assert [r.statement["topic"] for r in collected[0].rows] == ["carrier", "carrier"]
    assert collected[0].dropped["no_predicted_effect"] == 3
    pair = (("chr7", 1000, "C", "T"), ("chr7", 2000, "G", "A"))
    assert promote.candidate_pairs(collected, OUTSIDE_CARRIER, {}, ctx.haploid) == [
        pair
    ]
    trans = promote.promote_compound_hets(
        collected, LISTS, {}, ctx.haploid, {pair: Phase.IN_TRANS}
    )
    for row in trans[0].rows:
        assert (row.statement["kind"], row.statement["topic"]) == ("finding", "risks")
        assert IN_TRANS_REASON in str(row.statement["reason"])
    cis = promote.promote_compound_hets(
        collected, LISTS, {}, ctx.haploid, {pair: Phase.IN_CIS}
    )
    assert [r.statement["topic"] for r in cis[0].rows] == ["carrier", "carrier"]


def test_two_records_of_one_variant_are_two_statements(tmp_path: Path) -> None:
    """BR-013: a variant with records on two transcripts (a MANE Select missense, a MANE
    Plus Clinical stop) gives a statement for each, its id naming the transcript; the
    accounting counts each graded record once and each other variant once, so it
    balances, and the two records make no compound heterozygote."""
    stop = _record(gene="CFTR", vid="GRCh38:chr7:1000:C:T", transcript="ENST02.4")
    missense = _missense(revel=0.95, gene="CFTR", vid=stop.vid, transcript="ENST01.7")
    source = PredictedSource(
        sites=[(r, CALLED, RELIABLE, 0) for r in (missense, stop)],
        input_rows=3,
        escapes={"ENST02": REGIONS},
        truncations=[],
        haploinsufficient=HAPLOINSUFFICIENT,
    )
    ctx = _context(tmp_path)
    collected = source.collect(ctx)
    ids = [r.statement["statement_id"] for r in collected.rows]
    assert ids == [f"predicted:{stop.vid}:ENST01", f"predicted:{stop.vid}:ENST02"]
    assert collected.input_rows == len(ids) + sum(collected.dropped.values()) == 4
    assert promote.candidate_pairs([collected], LISTS, {}, ctx.haploid) == []


def test_one_copy_hets_outside_tier3_reach_the_promotion(tmp_path: Path) -> None:
    """One copy in a recessive gene outside the carrier lists (MUTYH without its Tier 4
    row: a missense, a LoF passing every gate) is stored as a non-finding, no longer
    dropped, so two make a candidate pair; a het LoF failing a gate is not
    pathogenic-class and is still dropped."""
    missense = _missense(revel=0.9, gene="MUTYH", vid="GRCh38:chr1:1000:G:A")
    lof = _record(gene="MUTYH", vid="GRCh38:chr1:2000:C:T")
    failing = _record(gene="MUTYH", vid="GRCh38:chr1:3000:C:T", mane=False)
    source = PredictedSource(
        sites=[(r, CALLED, RELIABLE, 0) for r in (missense, lof, failing)],
        input_rows=3,
        escapes={"ENST01": REGIONS},
        truncations=[],
        haploinsufficient=HAPLOINSUFFICIENT,
    )
    ctx = replace(_context(tmp_path), lists=OUTSIDE_CARRIER)
    collected = [source.collect(ctx)]
    assert [r.statement["vid"] for r in collected[0].rows] == [missense.vid, lof.vid]
    for row in collected[0].rows:
        assert row.statement["kind"] == "non_finding"
        assert ONE_RECESSIVE_COPY in str(row.statement["reason"])
    assert collected[0].dropped["not_haploinsufficient"] == 1
    pair = (("chr1", 1000, "G", "A"), ("chr1", 2000, "C", "T"))
    assert promote.candidate_pairs(collected, OUTSIDE_CARRIER, {}, ctx.haploid) == [
        pair
    ]


def test_research_hets_in_cgd_recessive_genes_pair(tmp_path: Path) -> None:
    """A missense het in a recessive gene on neither list (AR by CGD) stays research, yet
    pairs (design decision, 2026-10-28): with a ClinVar P het in GENEZ, with another
    predicted het in GENEY. In trans both become risks; in cis research stays research. A
    research statement that failed a LoF gate (CFTR, beside a carrier finding) pairs not."""
    genez = _missense(revel=0.95, gene="GENEZ", vid="GRCh38:chr2:1000:G:A")
    geney = [
        _missense(revel=0.95, gene="GENEY", vid=f"GRCh38:chr3:{pos}:G:A")
        for pos in (1000, 2000)
    ]
    cftr = [
        _record(gene="CFTR", vid="GRCh38:chr7:1000:C:T"),
        _record(gene="CFTR", vid="GRCh38:chr7:2000:C:T", mane=False),
    ]
    source = PredictedSource(
        sites=[(r, CALLED, RELIABLE, 0) for r in (genez, *geney, *cftr)],
        input_rows=5,
        escapes={"ENST01": REGIONS},
        truncations=[],
        haploinsufficient=HAPLOINSUFFICIENT,
    )
    inheritance = {"GENEZ": "AR", "GENEY": "AR"}
    ctx = replace(_context(tmp_path), inheritance=inheritance)
    predicted = source.collect(ctx)
    topics = [r.statement["topic"] for r in predicted.rows]
    assert topics == ["research", "research", "research", "carrier", "research"]
    vid = "GRCh38:chr2:500:A:G"
    clinvar = Collected(
        "clinvar",
        [
            Row(
                {"run_id": "run", "vid": vid, "rsid": None},
                {"run_id": "run", "vid": vid, "dosage": "1"},
                {
                    "statement_id": "clinvar:1",
                    "kind": "non_finding",
                    "source": "clinvar",
                    "vid": vid,
                    "gene": "GENEZ",
                    "topic": None,
                    "reason": None,
                    "claim": json.dumps({"classification": "Pathogenic"}),
                },
            )
        ],
        {},
        1,
        {},
    )
    collected = [clinvar, predicted]
    z = (("chr2", 500, "A", "G"), ("chr2", 1000, "G", "A"))
    y = (("chr3", 1000, "G", "A"), ("chr3", 2000, "G", "A"))
    assert promote.candidate_pairs(collected, LISTS, inheritance, ctx.haploid) == [z, y]
    phases = {z: Phase.IN_TRANS, y: Phase.IN_CIS}
    out = promote.promote_compound_hets(
        collected, LISTS, inheritance, ctx.haploid, phases
    )
    found = {r.statement["vid"]: r.statement for c in out for r in c.rows}
    for at in (vid, genez.vid):
        assert (found[at]["kind"], found[at]["topic"]) == ("finding", "risks")
        assert str(found[at]["reason"]).endswith(IN_TRANS_REASON)
    for record in geney:
        statement = found[record.vid]
        assert (statement["kind"], statement["topic"]) == ("finding", "research")
        assert str(statement["reason"]).endswith(IN_CIS_REASON)
    assert found[cftr[1].vid]["reason"] == predicted.rows[4].statement["reason"]


def test_predicted_never_outranks_classified(tmp_path: Path) -> None:
    """In the feed a predicted finding comes after every ClinVar-classified one, even one
    of lower validity in a later gene and not actionable."""
    con: duckdb.DuckDBPyConnection = connect(tmp_path / "genomelab.duckdb")
    seed(
        con,
        [
            StoreRow("predicted:1", validity="strong", actionable=True, gene="AAA"),
            StoreRow("clinvar:2", validity="limited", gene="ZZZ", pos=1002),
        ],
    )
    con.execute(
        "UPDATE statement SET source = 'vep' WHERE statement_id = 'predicted:1'"
    )
    found = query(con, None, Filters(min_validity=Validity.LIMITED)).findings
    assert [f["id"] for f in found] == ["clinvar:2", "predicted:1"]
