"""PharmGKB's clinical annotations (the ClinPGx download): levels 1B, 2A and 2B in genes
PharmCAT's report gives no drug, each located at its rsID's allele in the run's ClinVar GRCh38 release and
graded on the user's call there, as findings under drugs when carried. No annotation text,
and so no dosing text, is read.

The downloads here are `fixtures.pharmgkb_mini`'s, built per test; no ClinPGx data is
committed.
"""

import json
import zipfile
from dataclasses import replace
from pathlib import Path

import pysam
import pytest
from fixtures.clinvar_mini import write_clinvar
from fixtures.faf_mini import write_faf_table
from fixtures.grch38_mini import write_reference
from fixtures.pharmgkb_mini import Ann, Layout, write_annotations

from genokeeper.config import load_grading
from genokeeper.evidence.clinvar import NO_CLINPGX_REASON
from genokeeper.evidence.frequency import FafTable
from genokeeper.evidence.pharmgkb import (
    DOWNLOADS,
    RESEARCH_LEVEL,
    ClinicalAnnotations,
    by_vid,
    clinvar_alleles,
    left_by_pharmcat,
    locate,
    read_clinical_annotations,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import NOCALL_REASON, Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim
from genokeeper.pgx.pharmcat import DOSING, PHARMCAT_GENES
from genokeeper.sources import RunContext, write_collected
from genokeeper.sources.panel import grade_claim
from genokeeper.sources.pharmgkb import (
    PHARMGKB_RESEARCH_TOOL,
    PHARMGKB_TOOL,
    RESEARCH_REASON,
    PharmGkbSource,
)
from genokeeper.store.db import (
    Level,
    RowType,
    Validity,
    connect,
    start_run,
    transaction,
)
from genokeeper.store.writer import StoreWriter

TNF = Ann("655384799", "rs1800629", "TNF", "2B", ("AA", "AG", "GG"))
VKORC1 = Ann("655385028", "rs8050894", "VKORC1", "1B", ("CC", "CG", "GG"))
SCN1A = Ann(
    "613978931", "rs3812718", "SCN1A", "2B", category="Dosage", drugs="carbamazepine"
)
LEVEL_1A = Ann("1449309937", "rs1065852", "MTHFR", "1A")
LEVEL_3 = Ann("655384602", "rs951439", "RGS4", "3")
HAPLOTYPES = Ann("1139506787", "CYP2A6*1, CYP2A6*2", "CYP2A6", "1B", ("*1/*1",))
CHRNA5 = Ann(
    "981204902",
    "rs16969968",
    "CHRNA5",
    "1B",
    ("AA", "AG", "GG"),
    category="Toxicity;Dosage",
    drugs="nicotine",
    phenotypes="Tobacco Use Disorder",
)


def _download(
    tmp_path: Path, annotations: list[Ann], layout: Layout = "summary"
) -> Path:
    name = DOWNLOADS[0] if layout == "summary" else DOWNLOADS[1]
    return write_annotations(tmp_path / layout / name, annotations, layout=layout)


def _graded(read: ClinicalAnnotations) -> tuple[object, ...]:
    """Everything a download's reading gives the grading: the annotations kept (id,
    variant, gene, level, categories, drugs, phenotypes, genotypes), how many rows there
    were and where the others went."""
    return (read.annotations, read.input_rows, read.dropped, read.created)


def test_level_1a_ignored(tmp_path: Path) -> None:
    """1A restates a guideline and is never used; levels 3 and 4 are not either. Every
    row is counted: kept, or dropped by its level."""
    path = _download(tmp_path, [LEVEL_1A, LEVEL_3, TNF, CHRNA5])
    read = read_clinical_annotations(path)
    assert [a.annotation_id for a in read.annotations] == [TNF.id, CHRNA5.id]
    assert read.input_rows == 4
    assert {k: v for k, v in read.dropped.items() if v} == {"level_1a": 1, "level_3": 1}
    assert read.created == "20250705"
    (tnf, chrna5) = read.annotations
    assert (tnf.variant, tnf.gene, tnf.level, tnf.genotypes) == (
        "rs1800629",
        "TNF",
        "2B",
        ("AA", "AG", "GG"),
    )
    assert chrna5.drugs == ("nicotine",)
    assert chrna5.phenotypes == ("Tobacco Use Disorder",)


def test_the_live_download_comes_first() -> None:
    assert DOWNLOADS == ("summaryAnnotations.zip", "clinicalAnnotations.zip")


def test_reads_the_live_summary_layout(tmp_path: Path) -> None:
    """The live download (`summary_annotations.tsv` keyed by `Summary Annotation ID`) is
    read as the frozen one was: the same columns, the id kept as the annotation's id."""
    path = _download(tmp_path, [LEVEL_1A, LEVEL_3, TNF, CHRNA5], "summary")
    read = read_clinical_annotations(path)
    assert [a.annotation_id for a in read.annotations] == [TNF.id, CHRNA5.id]
    assert read.input_rows == 4
    assert {k: v for k, v in read.dropped.items() if v} == {"level_1a": 1, "level_3": 1}
    assert read.created == "20250705"
    (tnf, chrna5) = read.annotations
    assert (tnf.variant, tnf.gene, tnf.level, tnf.genotypes) == (
        "rs1800629",
        "TNF",
        "2B",
        ("AA", "AG", "GG"),
    )
    assert (chrna5.categories, chrna5.drugs) == (("Toxicity",), ("nicotine",))


def test_still_reads_the_frozen_clinical_layout(tmp_path: Path) -> None:
    """An older dated folder holds the frozen download (`clinical_annotations.tsv` keyed by
    `Clinical Annotation ID`); it is still read."""
    path = _download(tmp_path, [LEVEL_1A, TNF, CHRNA5], "clinical")
    assert path.name == "clinicalAnnotations.zip"
    read = read_clinical_annotations(path)
    assert [a.annotation_id for a in read.annotations] == [TNF.id, CHRNA5.id]
    assert read.input_rows == 3
    assert read.dropped["level_1a"] == 1


def test_same_grades_for_ids_in_both_layouts(tmp_path: Path) -> None:
    """The same rows in either layout give the same annotations, counts and ids, so a
    statement (`pharmgkb:<id>`) is graded alike across the switch."""
    rows = [LEVEL_1A, LEVEL_3, TNF, VKORC1, SCN1A, HAPLOTYPES, CHRNA5]
    old = read_clinical_annotations(_download(tmp_path, rows, "clinical"))
    new = read_clinical_annotations(_download(tmp_path, rows, "summary"))
    assert old.annotations
    assert _graded(old) == _graded(new)


def test_a_zip_with_neither_layout_names_both(tmp_path: Path) -> None:
    path = tmp_path / "other.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("CREATED_2026-10-05.txt", "Created.\n")
        z.writestr("annotations.tsv", "x\n")
    with pytest.raises(
        ValueError, match="summary_annotations.tsv.*clinical_annotations.tsv"
    ):
        read_clinical_annotations(path)


def test_a_zip_mixing_the_layouts_is_refused(tmp_path: Path) -> None:
    """A table of one layout beside the other's allele table is neither download."""
    path = tmp_path / "mixed.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("CREATED_2026-10-05.txt", "Created.\n")
        z.writestr("summary_annotations.tsv", "x\n")
        z.writestr("clinical_ann_alleles.tsv", "x\n")
    with pytest.raises(ValueError, match="not the ClinPGx"):
        read_clinical_annotations(path)


def test_pharmcat_genes_left_only_where_pharmcat_gives_no_drug(tmp_path: Path) -> None:
    """In a gene PharmCAT calls, an annotation is PharmCAT's to answer when its report gives
    that gene a drug (CYP2C19), and PharmGKB's when it gives none (F5, which PharmCAT 3.4.0
    lists only as an unannotated gene call): dropped, counted, or kept."""
    assert {"VKORC1", "CYP2D6", "HLA-B", "MT-RNR1", "F5"} <= PHARMCAT_GENES
    f5 = Ann(
        "1183689558", "rs6025", "F5", "2B", ("CC", "CT", "TT"), category="Toxicity"
    )
    cyp2c19 = Ann("1043880", "rs4244285", "CYP2C19", "1B")
    read = read_clinical_annotations(_download(tmp_path, [f5, cyp2c19, TNF]))
    assert read.dropped["pharmcat_gene"] == 0
    left, answered = left_by_pharmcat(read.annotations, {"CYP2C19", "CYP2D6"})
    assert [a.gene for a in left] == ["F5", "TNF"]
    assert answered == 1


def test_a_download_without_its_tables_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "clinicalAnnotations.zip"
    path.write_bytes(b"<html>maintenance</html>")
    with pytest.raises(ValueError, match="not the ClinPGx clinical annotations"):
        read_clinical_annotations(path)


@pytest.mark.tools
def test_annotations_located_at_their_clinvar_allele(tmp_path: Path) -> None:
    """An rsID's allele comes from the run's ClinVar GRCh38 release, normalised; it is kept
    only where its REF and ALT are both among the annotation's genotype bases (both on the
    plus strand), and only one allele may match."""
    fasta = write_reference(tmp_path / "ref.fa")
    with pysam.FastaFile(str(fasta)) as f:
        base = {pos: f.fetch("chr1", pos - 1, pos) for pos in (100, 200, 300, 400)}

    def other(b: str, skip: str = "") -> str:
        return next(x for x in "ACGT" if x not in (b, skip))

    flip = str.maketrans("ACGT", "TGCA")
    alt = {pos: other(b) for pos, b in base.items()}
    records = [
        ("1", 100, 1, base[100], alt[100], {"CLNSIG": "drug_response", "RS": "100"}),
        ("1", 200, 2, base[200], alt[200], {"CLNSIG": "Benign", "RS": "200"}),
        ("1", 300, 3, base[300], alt[300], {"CLNSIG": "Benign", "RS": "300"}),
        ("1", 300, 4, base[300], other(base[300], alt[300]), {"RS": "300"}),
        ("1", 400, 5, base[400], f"{alt[400]}A", {"RS": "400"}),
    ]
    clinvar = write_clinvar(tmp_path / "clinvar.vcf.gz", "GRCh38", records)

    def labels(pos: int, *, flipped: bool = False) -> tuple[str, ...]:
        r, a = base[pos], alt[pos]
        if flipped:
            r, a = r.translate(flip), a.translate(flip)
        return (r + r, r + a, a + a)

    annotations = [
        Ann("1", "rs100", "GENEA", "2A", labels(100)),
        Ann("2", "rs200", "GENEB", "2B", labels(200, flipped=True)),
        Ann("3", "rs300", "GENEC", "2B", (base[300] * 2, alt[300] + base[300])),
        Ann("4", "rs400", "GENED", "2B", labels(400)),
        Ann("5", "rs999", "GENEE", "2B"),
        HAPLOTYPES,
    ]
    read = read_clinical_annotations(_download(tmp_path, annotations))
    index = clinvar_alleles(clinvar, {a.variant for a in read.annotations}, fasta=fasta)
    assert index["rs100"] == [("chr1", 100, base[100], alt[100])]
    located, dropped = locate(read.annotations, index)
    assert [(a.annotation_id, allele) for a, allele in located] == [
        ("1", ("chr1", 100, base[100], alt[100])),
        ("3", ("chr1", 300, base[300], alt[300])),
    ]
    assert {k: v for k, v in dropped.items() if v} == {
        "alleles_differ": 1,  # the minus-strand bases of rs200
        "not_a_snv": 1,  # rs400 is an insertion
        "rsid_not_in_clinvar": 1,
        "not_an_rsid": 1,
    }


def _context(tmp_path: Path) -> RunContext:
    return RunContext(
        run_id="run-1",
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
        cpic_levels={("CHRNA5", "nicotine"): "B"},
    )


def _called(dosage: str) -> Called:
    gt = {"0": "0/0", "1": "0/1", "2": "1/1", "unknown": "./."}[dosage]
    warnings = (NOCALL_REASON,) if dosage == "unknown" else ()
    row = RowType.NOCALL if dosage == "unknown" else RowType.PASS
    return Called(gt, dosage, row, 40, 30, None, None, warnings)


RELIABLE = Reliability(p_wrong=1e-6, level=Level.HIGH, flags=(), warnings=())


def test_carried_annotation_is_finding(tmp_path: Path) -> None:
    """Carried (one or two copies): a finding under drugs; not carried: a non-finding;
    unknown: unassessed, never a finding. 1B is strong, 2A and 2B moderate. The claim
    holds no annotation text and no dosing word: a category naming dose is left out."""
    read = read_clinical_annotations(
        _download(tmp_path, [TNF, CHRNA5, SCN1A, Ann("9", "rs9", "GENEX", "2A")])
    )
    tnf, chrna5, scn1a, genex = read.annotations
    rows = [
        (chrna5, ("chr15", 78590583, "G", "A"), _called("2"), RELIABLE),
        (tnf, ("chr6", 31575254, "G", "A"), _called("1"), RELIABLE),
        (scn1a, ("chr2", 166078000, "C", "T"), _called("0"), RELIABLE),
        (genex, ("chr1", 1000, "C", "T"), _called("unknown"), RELIABLE),
    ]
    source = PharmGkbSource(
        rows, input_rows=6, dropped={"level_1a": 2}, version=read.created
    )
    collected = source.collect(_context(tmp_path))
    assert collected.tool == PHARMGKB_TOOL == "pharmgkb"
    by_id = {r.statement["statement_id"]: r for r in collected.rows}
    chr15 = by_id[f"pharmgkb:{chrna5.annotation_id}"]
    s = chr15.statement
    assert (s["kind"], s["topic"], s["validity"]) == ("finding", "drugs", "strong")
    assert (s["source"], s["source_id"], s["source_version"]) == (
        "pharmgkb",
        chrna5.annotation_id,
        "20250705",
    )
    assert s["source_level"] == "ClinPGx level 1B"
    assert s["gene"] == "CHRNA5" and s["vid"] == "GRCh38:chr15:78590583:G:A"
    assert s["actionable"] is True  # CPIC level B for (CHRNA5, nicotine)
    assert chr15.variant == {
        "run_id": "run-1",
        "vid": "GRCh38:chr15:78590583:G:A",
        "rsid": "rs16969968",
    }
    assert chr15.call is not None and chr15.call["dosage"] == "2"
    claim = json.loads(str(s["claim"]))
    assert claim == {
        "schema": "genomelab.claim/1",
        "type": "pharmgkb_clinical_annotation",
        "annotation_id": chrna5.annotation_id,
        "gene": "CHRNA5",
        "category": "Toxicity",
        "condition": "Tobacco Use Disorder",
        "drugs": [{"drug": "nicotine"}],
    }
    tnf_row = by_id[f"pharmgkb:{tnf.annotation_id}"].statement
    assert (tnf_row["kind"], tnf_row["validity"], tnf_row["actionable"]) == (
        "finding",
        "moderate",
        False,
    )
    scn1a_row = by_id[f"pharmgkb:{scn1a.annotation_id}"].statement
    assert scn1a_row["kind"] == "non_finding"
    assert json.loads(str(scn1a_row["claim"]))["category"] is None
    assert by_id[f"pharmgkb:{genex.annotation_id}"].statement["kind"] == "unassessed"
    for row in collected.rows:
        text = json.dumps(row.statement)
        assert not DOSING.search(text), text
        assert "decreased" not in text  # no annotation text
    assert (collected.input_rows, collected.dropped) == (6, {"level_1a": 2})


# --- slice 3f: level 3 as a research layer ------------------------------------------------

RESEARCH = Ann(
    "655384602",
    "rs951439",
    "RGS4",
    "3",
    ("CC", "CT", "TT"),
    category="Efficacy;Dosage",
    drugs="risperidone",
    phenotypes="Schizophrenia",
)


def test_level_3_is_research(tmp_path: Path) -> None:
    """A level-3 row is not the `pharmgkb` tool's statement (still counted `level_3`
    there) but is carried to the research layer: a `pharmgkb_research` finding under
    research, limited, never actionable, with its own reason and source level."""
    read = read_clinical_annotations(_download(tmp_path, [RESEARCH, TNF]))
    assert RESEARCH_LEVEL == "3"
    assert [a.annotation_id for a in read.annotations] == [TNF.id]
    assert [a.annotation_id for a in read.research] == [RESEARCH.id]
    assert read.dropped["level_3"] == 1
    [rgs4] = read.research
    source = PharmGkbSource(
        [(rgs4, ("chr1", 163150000, "C", "T"), _called("1"), RELIABLE)],
        input_rows=1,
        dropped={},
        version=read.created,
        tool=PHARMGKB_RESEARCH_TOOL,
    )
    collected = source.collect(_context(tmp_path))
    assert collected.tool == PHARMGKB_RESEARCH_TOOL == "pharmgkb_research"
    [row] = collected.rows
    s = row.statement
    assert s["statement_id"] == f"pharmgkb_research:{RESEARCH.id}"
    assert (s["kind"], s["topic"], s["validity"], s["actionable"]) == (
        "finding",
        "research",
        "limited",
        False,
    )
    assert (s["reason"], s["source_level"]) == (RESEARCH_REASON, "ClinPGx level 3")
    assert (s["source"], s["source_id"]) == ("pharmgkb_research", RESEARCH.id)
    assert RESEARCH_REASON == "one study or an unconfirmed association"


def test_level_3_dropped_in_pharmcat_gene_and_dosing(tmp_path: Path) -> None:
    """In a gene PharmCAT gives a drug, a level-3 row is PharmCAT's (`pharmcat_gene`), as a
    1B–2B one is; its dosing category is left out as theirs is."""
    cyp2c19 = Ann("1043999", "rs12248560", "CYP2C19", "3", category="Dosage;Toxicity")
    read = read_clinical_annotations(_download(tmp_path, [cyp2c19, RESEARCH]))
    left, answered = left_by_pharmcat(read.research, {"CYP2C19"})
    assert [a.annotation_id for a in left] == [RESEARCH.id]
    assert answered == 1
    by_id = {a.annotation_id: a for a in read.research}
    assert by_id[cyp2c19.id].categories == ("Toxicity",)
    assert by_id[RESEARCH.id].categories == ("Efficacy",)
    for a in read.research:
        assert not any(DOSING.search(x) for x in (*a.categories, *a.drugs))


def test_level_4_still_dropped(tmp_path: Path) -> None:
    level_4 = Ann("655384999", "rs4", "GENE4", "4")
    read = read_clinical_annotations(_download(tmp_path, [level_4, TNF]))
    assert [a.annotation_id for a in read.annotations] == [TNF.id]
    assert read.research == []
    assert {k: v for k, v in read.dropped.items() if v} == {"level_4": 1}


def test_research_accounting_balances(tmp_path: Path) -> None:
    """Each tool's statements plus its drops are its input rows: `pharmgkb` counts the
    level-3 row as `level_3`, `pharmgkb_research` its level-3 rows by where they went."""
    cyp2c19 = Ann("1043999", "rs12248560", "CYP2C19", "3")
    unlocated = Ann("1043998", "rs999", "GENEZ", "3")
    read = read_clinical_annotations(
        _download(tmp_path, [RESEARCH, cyp2c19, unlocated, TNF, LEVEL_1A])
    )
    alleles = {
        "rs951439": [("chr1", 163150000, "C", "T")],
        "rs1800629": [("chr6", 31575254, "G", "A")],
    }
    sources = []
    for annotations, tool, input_rows, drops in (
        (read.annotations, PHARMGKB_TOOL, read.input_rows, read.dropped),
        (read.research, PHARMGKB_RESEARCH_TOOL, len(read.research), {}),
    ):
        left, answered = left_by_pharmcat(annotations, {"CYP2C19"})
        located, dropped = locate(left, alleles)
        dropped["pharmcat_gene"] += answered
        sources.append(
            PharmGkbSource(
                [(a, allele, _called("1"), RELIABLE) for a, allele in located],
                input_rows=input_rows,
                dropped={r: drops.get(r, 0) + dropped[r] for r in dropped},
                version=read.created,
                tool=tool,
            )
        )
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {"pharmgkb": read.created})
    ctx = replace(_context(tmp_path), run_id=run_id)
    with transaction(con):
        write_collected(
            StoreWriter(con, tmp_path / "work"),
            con,
            run_id,
            [s.collect(ctx) for s in sources],
        )
    rows = con.execute(
        "SELECT tool, input_rows, findings, dropped FROM accounting ORDER BY tool"
    ).fetchall()
    tallies = {
        tool: (n, found, {k: v for k, v in json.loads(d).items() if v})
        for tool, n, found, d in rows
    }
    assert tallies == {
        "pharmgkb": (5, 1, {"level_1a": 1, "level_3": 3}),
        "pharmgkb_research": (
            3,
            1,
            {"pharmcat_gene": 1, "rsid_not_in_clinvar": 1},
        ),
    }


def test_aldh2_rs671_reaches_drugs(tmp_path: Path) -> None:
    """Slice 3f, "Traits": ALDH2 rs671 is in ClinPGx at 2B (annotation 1450810451,
    ethanol, Toxicity), outside PharmCAT's genes, so the existing path already gives a
    carrier a statement under drugs; no trait statement repeats it."""
    aldh2 = Ann(
        "1450810451",
        "rs671",
        "ALDH2",
        "2B",
        ("AA", "AG", "GG"),
        category="Toxicity",
        drugs="ethanol",
        phenotypes="Alcohol abuse",
    )
    read = read_clinical_annotations(_download(tmp_path, [aldh2]))
    (kept,) = read.annotations
    left, answered = left_by_pharmcat(read.annotations, set())
    assert (left, answered) == ([kept], 0)
    located, _ = locate(left, {"rs671": [("chr12", 111803962, "G", "A")]})
    source = PharmGkbSource(
        [(a, allele, _called("1"), RELIABLE) for a, allele in located],
        input_rows=1,
        dropped={},
        version=read.created,
    )
    (row,) = source.collect(_context(tmp_path)).rows
    s = row.statement
    assert (s["statement_id"], s["kind"], s["topic"], s["validity"]) == (
        "pharmgkb:1450810451",
        "finding",
        "drugs",
        "moderate",
    )
    assert s["gene"] == "ALDH2" and s["vid"] == "GRCh38:chr12:111803962:G:A"


def test_one_row_per_variant_on_clinpgx(tmp_path: Path) -> None:
    """Owner's ruling A3s: a ClinVar drug_response claim at an rsID with a kept ClinPGx
    annotation is graded on it (limited, naming it), so the variant has one row at its
    ClinPGx level; the strongest kept level answers for an allele, and it does not grade a
    claim at another allele of the same rsID."""
    mthfr = Ann("1000000002", "rs900002", "MTHFR", "2A", ("CC", "CT", "TT"))
    weaker = Ann("1000000003", "rs900002", "MTHFR", "2B", ("CC", "CT", "TT"))
    read = read_clinical_annotations(_download(tmp_path, [weaker, mthfr]))
    allele = ("chr1", 11796321, "G", "A")
    ctx = replace(
        _context(tmp_path), clinpgx=by_vid((a, allele) for a in read.annotations)
    )
    assert ctx.clinpgx == {"GRCh38:chr1:11796321:G:A": ("1000000002", "2A")}
    kept = next(a for a in read.annotations if a.level == "2A")
    rows = (
        PharmGkbSource([(kept, allele, _called("1"), RELIABLE)], 1, {}, read.created)
        .collect(ctx)
        .rows
    )
    claim = Claim(
        variation_id=3520,
        gene="MTHFR",
        clnsig="drug_response",
        clnrevstat="reviewed_by_expert_panel",
        disease_ids="",
        condition="methotrexate response",
        chrom="chr1",
        pos=11796321,
        ref="G",
        alt="A",
        vid="GRCh38:chr1:11796321:G:A",
        rsid="rs900002",
    )
    grade = grade_claim(claim, _called("1"), ctx, None)
    assert [r.statement["validity"] for r in rows] == ["moderate"]
    assert (grade.level, grade.reason) == (
        Validity.LIMITED,
        "graded on ClinPGx's annotation 1000000002 (level 2A)",
    )
    other = replace(claim, alt="T", vid="GRCh38:chr1:11796321:G:T")
    assert grade_claim(other, _called("1"), ctx, None).reason == NO_CLINPGX_REASON
