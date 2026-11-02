"""PharmCAT: our own calls written at every position PharmCAT defines, the arguments it
runs with, the run in its container, and its JSON report read without dosing text.

The positions are a few rows in the shape of PharmCAT 3.4.0's
`pharmcat_positions_3.4.0.vcf.bgz` (a release asset): a SNV, a multi-allelic SNV, a deletion
and G6PD's chrX positions, with PharmCAT's own REF and ALT spellings.

`fixtures/pharmcat_report_mini.json` is PharmCAT's published example report
(`docs/examples/pharmcat.example.report.json` at tag v3.4.0, a synthetic sample) trimmed to
five genes (CYP2C19, CYP2D6, DPYD, HLA-A, HLA-B), seven drug annotations of three sources,
three variants per gene, one citation per drug and one unannotated gene call; one edit: the
call of CYP2C19's second variant (chr10:94762706) is blanked, which is how the report lists
a missing position (PharmCAT's `VariantReport.isMissing`).
"""

import gzip
import json
from pathlib import Path
from typing import Any

import pytest
from fixtures.gvcf_mini import block, other_base, ref, variant, write_fasta, write_gvcf

from genokeeper.genotype.gvcf import (
    NOCALL_REASON,
    OTHER_ALLELE_WARNING,
    Allele,
    Called,
    calls_at,
)
from genokeeper.pgx import pharmcat
from genokeeper.pgx.pharmcat import (
    DOSING,
    DrugRec,
    GeneResult,
    PharmcatFailed,
    parse_report,
    pharmcat_args,
    pharmcat_positions,
    run_pharmcat,
    write_pharmcat_vcf,
)
from genokeeper.store.db import Level, RowType

REPORT = Path(__file__).parent / "fixtures" / "pharmcat_report_mini.json"

POSITIONS = """\
##fileformat=VCFv4.2
##source=PharmCAT allele definitions
##contig=<ID=chr6,assembly=GRCh38.p14,species="Homo sapiens">
##contig=<ID=chr10,assembly=GRCh38.p14,species="Homo sapiens">
##contig=<ID=chrX,assembly=GRCh38.p14,species="Homo sapiens">
##INFO=<ID=PX,Number=.,Type=String,Description="Gene">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT
chr6\t18130687\trs1142345\tT\tC,G\t.\tPASS\tPX=TPMT\tGT\t0/0
chr10\t94761900\trs12248560\tC\tT\t.\tPASS\tPX=CYP2C19\tGT\t0/0
chr10\t94949281\trs9332131\tGA\tG\t.\tPASS\tPX=CYP2C9\tGT\t0/0
chrX\t154532046\t.\tA\tC\t.\tPASS\tPX=G6PD\tGT\t0/0
chrX\t154532055\t.\tCTCT\tC\t.\tPASS\tPX=G6PD\tGT\t0/0
"""
TPMT_C = ("chr6", 18130687, "T", "C")
TPMT_G = ("chr6", 18130687, "T", "G")
CYP2C19_17 = ("chr10", 94761900, "C", "T")
CYP2C9_DEL = ("chr10", 94949281, "GA", "G")
G6PD_SNV = ("chrX", 154532046, "A", "C")
G6PD_DEL = ("chrX", 154532055, "CTCT", "C")


def _called(
    dosage: str,
    gt: str | None,
    *,
    row_type: RowType = RowType.PASS,
    warnings: tuple[str, ...] = (),
) -> Called:
    return Called(
        gt=gt,
        dosage=dosage,
        row_type=row_type,
        gq=40,
        dp=30,
        ad=None,
        vaf=None,
        warnings=warnings,
    )


HET = _called("1", "0/1")
HOM_REF = _called("0", "0/0", row_type=RowType.REF_BLOCK)
NOCALL = _called("unknown", "./.", row_type=RowType.NOCALL, warnings=(NOCALL_REASON,))


@pytest.fixture
def positions(tmp_path: Path) -> Path:
    path = tmp_path / "pharmcat_positions_3.4.0.vcf.bgz"
    path.write_bytes(gzip.compress(POSITIONS.encode()))
    return path


def _rows(vcf: Path) -> list[list[str]]:
    lines = vcf.read_text().splitlines()
    header = [line for line in lines if line.startswith("#CHROM")]
    assert header == ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
    assert lines[0] == "##fileformat=VCFv4.2"
    return [line.split("\t") for line in lines if not line.startswith("#")]


def _write(
    calls: dict[Allele, Called], positions: list[Allele], out: Path, *, male: bool
) -> Path:
    """`write_pharmcat_vcf` with every call graded high."""
    levels = dict.fromkeys(calls, Level.HIGH)
    return write_pharmcat_vcf(calls, positions, out, male=male, levels=levels)


def _genotypes(vcf: Path) -> dict[tuple[str, int], str]:
    return {(row[0], int(row[1])): row[9] for row in _rows(vcf)}


def test_positions_split_by_alt(positions: Path) -> None:
    assert pharmcat_positions(positions) == [
        TPMT_C,
        TPMT_G,
        CYP2C19_17,
        CYP2C9_DEL,
        G6PD_SNV,
        G6PD_DEL,
    ]


def test_every_position_written(positions: Path, tmp_path: Path) -> None:
    """One row per PharmCAT position, in its order and with its own REF and ALTs, whether
    or not a call was made there: a position left out would be PharmCAT's no-call too, but
    written as one it says so."""
    out = _write(
        {CYP2C19_17: HET},
        pharmcat_positions(positions),
        tmp_path / "pharmcat.vcf",
        male=False,
    )
    assert out == tmp_path / "pharmcat.vcf"
    rows = _rows(out)
    assert [(r[0], int(r[1]), r[3], r[4], r[8]) for r in rows] == [
        ("chr6", 18130687, "T", "C,G", "GT"),
        ("chr10", 94761900, "C", "T", "GT"),
        ("chr10", 94949281, "GA", "G", "GT"),
        ("chrX", 154532046, "A", "C", "GT"),
        ("chrX", 154532055, "CTCT", "C", "GT"),
    ]
    assert [r[9] for r in rows] == ["./.", "0/1", "./.", "./.", "./."]


def test_multiallelic_position_from_each_alts_dosage(
    positions: Path, tmp_path: Path
) -> None:
    calls = {TPMT_C: HET, TPMT_G: HET, CYP2C19_17: _called("2", "1/1")}
    calls[CYP2C9_DEL] = _called("1", "0/1")
    out = _write(calls, pharmcat_positions(positions), tmp_path / "p.vcf", male=False)
    genotypes = _genotypes(out)
    assert genotypes[("chr6", 18130687)] == "1/2"
    assert genotypes[("chr10", 94761900)] == "1/1"
    assert genotypes[("chr10", 94949281)] == "0/1"


@pytest.mark.tools
def test_het_at_a_multiallelic_position_from_the_call_layer(tmp_path: Path) -> None:
    """A `0/1` record of one listed ALT at a position PharmCAT lists several ALTs for is that
    het: `calls_at` gives each sibling ALT dosage 0 with `OTHER_ALLELE_WARNING`, which here
    means the listed ALT, not an allele PharmCAT does not list (TPMT rs1142345 `T>C,G` and
    UGT1A1 rs3064744 `CAT>C,CATAT,CATATAT` shapes, on the mini reference). An allele the
    position does not list is still `./.`; `1/2` and `1/1` stay as called."""
    snv, indel, unlisted, both, hom = 2001, 2201, 2401, 2601, 2801
    base = {p: ref("chr1", p) for p in (snv, unlisted, both, hom)}
    c, g, t = (other_base(base[snv], k) for k in range(3))
    cat = ref("chr1", indel, 3)
    positions = [
        ("chr1", snv, base[snv], c),
        ("chr1", snv, base[snv], g),
        *[("chr1", indel, cat, alt) for alt in (cat[0], cat + cat[1:], cat + 2 * cat[1:])],
        *[("chr1", p, base[p], other_base(base[p], k)) for p in (unlisted, both, hom)
          for k in (0, 1)],
    ]  # fmt: skip
    records = [
        block("chr1", 1900, snv - 1),
        variant("chr1", snv, 1, [c, t], gt="0/1", gq=43, ad=[10, 10, 6]),
        block("chr1", snv + 1, indel - 1),
        variant("chr1", indel, 3, [cat + cat[1:], cat[0]], gt="0/1", gq=40, ad=[9, 9, 2]),
        block("chr1", indel + 3, unlisted - 1),
        variant("chr1", unlisted, 1, [other_base(base[unlisted], 2)], gt="0/1", gq=40,
                ad=[10, 10]),
        block("chr1", unlisted + 1, both - 1),
        variant("chr1", both, 1, [other_base(base[both], k) for k in (0, 2, 1)],
                gt="1/3", gq=40, ad=[0, 10, 1, 10]),
        block("chr1", both + 1, hom - 1),
        variant("chr1", hom, 1, [other_base(base[hom], 1)], gt="1/1", gq=40, ad=[0, 20]),
        block("chr1", hom + 1, 3000),
    ]  # fmt: skip
    fasta = write_fasta(tmp_path / "ref.fa")
    gvcf = write_gvcf(tmp_path / "mini.g.vcf.gz", records)
    calls = calls_at(
        gvcf, positions, fasta=fasta, haploid=lambda *_: False, scratch=tmp_path
    )
    assert OTHER_ALLELE_WARNING in calls[positions[1]].warnings  # the sibling G
    levels = dict.fromkeys(positions, Level.HIGH)
    out = write_pharmcat_vcf(
        calls, positions, tmp_path / "p.vcf", male=False, levels=levels
    )
    assert _genotypes(out) == {
        ("chr1", snv): "0/1",
        ("chr1", indel): "0/2",
        ("chr1", unlisted): "./.",
        ("chr1", both): "1/2",
        ("chr1", hom): "2/2",
    }


def test_unknown_is_nocall_never_ref(positions: Path, tmp_path: Path) -> None:
    """`0/0` only where the call says no copy of any of the position's ALTs and nothing
    else was called there; an unknown, a missing call, or another allele PharmCAT does not
    list is `./.` (`--missing-to-ref` would read it as reference)."""
    other = _called("0", "0/2", warnings=(OTHER_ALLELE_WARNING,))
    calls = {
        TPMT_C: HOM_REF,
        TPMT_G: HOM_REF,
        CYP2C19_17: NOCALL,
        CYP2C9_DEL: other,
    }
    out = _write(calls, pharmcat_positions(positions), tmp_path / "p.vcf", male=False)
    genotypes = _genotypes(out)
    assert genotypes[("chr6", 18130687)] == "0/0"
    assert genotypes[("chr10", 94761900)] == "./."
    assert genotypes[("chr10", 94949281)] == "./."
    # One ALT of a multi-allelic position unknown: the whole position is.
    calls[TPMT_G] = NOCALL
    out = _write(calls, pharmcat_positions(positions), tmp_path / "p.vcf", male=False)
    assert _genotypes(out)[("chr6", 18130687)] == "./."


def test_low_reliability_is_nocall(positions: Path, tmp_path: Path) -> None:
    """A call graded low (or not graded) is written as no call: PharmCAT then reports the
    position missing, and the gene call's reliability drops with it; medium and high calls
    are written as called."""
    alleles = pharmcat_positions(positions)
    calls = {
        TPMT_C: HET,
        TPMT_G: HOM_REF,
        CYP2C19_17: HET,
        CYP2C9_DEL: _called("1", "0/1"),
        G6PD_SNV: _called("1", "1"),
        G6PD_DEL: _called("0", "0"),
    }
    levels = {
        TPMT_C: Level.HIGH,
        TPMT_G: Level.LOW,  # one ALT of the position low: the position is
        CYP2C19_17: Level.MEDIUM,
        G6PD_SNV: Level.LOW,
        G6PD_DEL: Level.HIGH,
    }  # CYP2C9_DEL not graded
    out = tmp_path / "p.vcf"
    written = write_pharmcat_vcf(calls, alleles, out, male=True, levels=levels)
    assert _genotypes(written) == {
        ("chr6", 18130687): "./.",
        ("chr10", 94761900): "0/1",
        ("chr10", 94949281): "./.",
        ("chrX", 154532046): ".",
        ("chrX", 154532055): "0",
    }


def test_male_g6pd_haploid(positions: Path, tmp_path: Path) -> None:
    """G6PD sits on chrX outside the PARs: a male's genotype there is one allele. The call
    layer writes his homozygote as one index; a heterozygote there is no call."""
    calls = {
        G6PD_SNV: _called("1", "1"),
        G6PD_DEL: _called("0", "0", row_type=RowType.REF_BLOCK),
        CYP2C19_17: HET,
    }
    alleles = pharmcat_positions(positions)
    male = _genotypes(_write(calls, alleles, tmp_path / "m.vcf", male=True))
    assert male[("chrX", 154532046)] == "1"
    assert male[("chrX", 154532055)] == "0"
    assert male[("chr10", 94761900)] == "0/1"  # autosomes stay diploid
    assert male[("chr6", 18130687)] == "./."
    calls[G6PD_SNV] = NOCALL
    calls[G6PD_DEL] = _called("1", "0/1")
    male = _genotypes(_write(calls, alleles, tmp_path / "m.vcf", male=True))
    assert male[("chrX", 154532046)] == "."
    assert male[("chrX", 154532055)] == "."
    female = {G6PD_SNV: HET, G6PD_DEL: _called("0", "0/0")}
    written = _write(female, alleles, tmp_path / "f.vcf", male=False)
    genotypes = _genotypes(written)
    assert genotypes[("chrX", 154532046)] == "0/1"
    assert genotypes[("chrX", 154532055)] == "0/0"


@pytest.mark.parametrize("outside", [None, "outside.tsv", "-0.tsv"])
@pytest.mark.parametrize("vcf", ["pharmcat.vcf", "-0.vcf", "--missing-to-ref.vcf"])
def test_args_never_missing_to_ref(
    tmp_path: Path, vcf: str, outside: str | None
) -> None:
    """PharmCAT's own `pharmcat` reads the VCF as written: a missing position is a no-call.
    `--missing-to-ref` (`-0`) is the preprocessor's, which is never run, and never passed."""
    args = pharmcat_args(
        tmp_path / vcf,
        None if outside is None else tmp_path / outside,
        tmp_path / "out",
    )
    assert "--missing-to-ref" not in args and "-0" not in args
    assert not any(a.startswith("--missing") for a in args)
    assert args[0] == "pharmcat"
    assert "-reporterJson" in args and "-reporterHtml" in args
    assert ("-po" in args) == (outside is not None)


# --- the run and the report ----------------------------------------------------------------


def _fake_container(calls: list[dict[str, Any]], *, report: bool = True) -> Any:
    """Stands in for `run_container`: records the call and writes PharmCAT's outputs into
    the read-write mount."""

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        calls.append(
            {"name": name, "args": args, "mounts": mounts, "platform": platform}
            | {"label": label}
        )
        (out,) = [host for host, _, read_only in mounts if not read_only]
        if report:
            (out / "pharmcat.report.json").write_text(REPORT.read_text())
            (out / "pharmcat.report.html").write_text("<html></html>")

    return fake


def test_run_offline_in_the_pinned_image(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[dict[str, Any]] = []
    monkeypatch.setattr(pharmcat, "run_container", _fake_container(calls))
    vcf = tmp_path / "in" / "pharmcat.vcf"
    outside = tmp_path / "outside" / "outside.tsv"
    for path in (vcf, outside):
        path.parent.mkdir()
        path.write_text("x")
    out = tmp_path / "out"
    report = run_pharmcat(vcf, outside, out)
    assert report == out / "pharmcat.report.json"
    assert (out / "pharmcat.report.html").is_file()  # kept for the user
    (call,) = calls
    assert (call["name"], call["platform"], call["label"]) == (
        "pharmcat",
        None,
        "dna:pgx",
    )
    assert call["args"] == pharmcat_args(vcf, outside, out)
    assert sorted(call["mounts"]) == sorted(
        [
            (vcf.parent, "/pharmcat-in", True),
            (outside.parent, "/pharmcat-outside", True),
            (out, "/pharmcat-out", False),
        ]
    )
    calls.clear()
    run_pharmcat(vcf, None, out)
    assert [m[1] for m in calls[0]["mounts"]] == ["/pharmcat-in", "/pharmcat-out"]


def test_run_without_a_report_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[dict[str, Any]] = []
    monkeypatch.setattr(pharmcat, "run_container", _fake_container(calls, report=False))
    vcf = tmp_path / "pharmcat.vcf"
    vcf.write_text("x")
    with pytest.raises(PharmcatFailed, match=r"pharmcat\.report\.json"):
        run_pharmcat(vcf, None, tmp_path / "out")


CYP2C19_IMPLICATION = (
    "CYP2C19: Normal clopidogrel active metabolite formation; normal on-treatment"
    " platelet reactivity"
)
DPWG_CLOPIDOGREL = (
    "The guideline does not provide a description of the impact of a normal metabolizer"
    " phenotype on clopidogrel."
)


def test_parse_gene_and_drugs() -> None:
    results = {r.gene: r for r in parse_report(REPORT)}
    assert list(results) == ["CYP2C19", "CYP2D6", "DPYD", "HLA-A", "HLA-B", "F2"]
    # A gene PharmCAT calls without a drug (`unannotatedGeneCalls`): its call, no drugs.
    assert results["F2"] == GeneResult(
        gene="F2",
        diplotype="rs1799963 reference (G)/rs1799963 reference (G)",
        phenotype=None,
        activity_score=None,
        call_source="MATCHER",
        missing_positions=(),
        drugs=(),
    )
    assert results["CYP2C19"] == GeneResult(
        gene="CYP2C19",
        diplotype="*38/*38",
        phenotype="Normal Metabolizer",
        activity_score=None,
        call_source="MATCHER",
        missing_positions=("chr10:94762706",),
        drugs=(
            # Two CPIC annotations (populations) of one strength: one recommendation.
            DrugRec("clopidogrel", "CPIC", "Strong", (CYP2C19_IMPLICATION,)),
            DrugRec("clopidogrel", "DPWG", "No recommendation", (DPWG_CLOPIDOGREL,)),
        ),
    )
    cyp2d6 = results["CYP2D6"]
    assert (cyp2d6.diplotype, cyp2d6.phenotype, cyp2d6.activity_score) == (
        "*1/*3",
        "Intermediate Metabolizer",
        "1.0",
    )
    assert (cyp2d6.call_source, cyp2d6.missing_positions) == ("OUTSIDE", ())
    assert [(d.drug, d.source, d.classification) for d in cyp2d6.drugs] == [
        ("codeine", "CPIC", "Moderate"),
        ("atomoxetine", "DPWG", "Unspecified"),
    ]
    assert results["DPYD"].activity_score == "2.0"
    assert results["HLA-B"].phenotype == (
        "*15:02 positive; *57:01 positive; *58:01 negative"
    )
    assert results["HLA-B"].drugs == (
        DrugRec(
            "abacavir",
            "CPIC",
            "Strong",
            ("HLA-B: Significantly increased risk of abacavir hypersensitivity",),
        ),
    )
    # PharmCAT had no data and no outside call: no diplotype, no phenotype.
    hla_a = results["HLA-A"]
    assert (hla_a.diplotype, hla_a.phenotype, hla_a.call_source) == (None, None, "NONE")


def test_haploid_unknown_is_no_call(tmp_path: Path) -> None:
    """PharmCAT labels a one-allele gene it cannot call `Unknown` (a diploid one
    `Unknown/Unknown`), with the phenotype `No Result`, as its run on a male's G6PD with
    every position missing showed: no diplotype, no phenotype."""
    report = json.loads(REPORT.read_text())
    (diplotype,) = report["genes"]["DPYD"]["sourceDiplotypes"]
    diplotype.update(
        label="Unknown",
        allele2=None,
        phenotypes=["No Result"],
        activityScore=None,
        lookupKey=["No Result"],
    )
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    dpyd = {r.gene: r for r in parse_report(path)}["DPYD"]
    assert (dpyd.diplotype, dpyd.phenotype, dpyd.activity_score) == (None, None, None)


def _dosing_texts(value: Any) -> list[str]:
    """Every `drugRecommendation` string in a report."""
    if isinstance(value, dict):
        found = [v for k, v in value.items() if k == "drugRecommendation"]
        return found + [t for v in value.values() for t in _dosing_texts(v)]
    if isinstance(value, list):
        return [t for v in value for t in _dosing_texts(v)]
    return []


def test_dosing_text_absent_from_results(tmp_path: Path) -> None:
    """No recommendation text reaches a result, and no implication that speaks of a dose;
    the report's `drugRecommendation` and `dosingInformation` are never read: without them
    the results are the same."""
    report = json.loads(REPORT.read_text())
    texts = _dosing_texts(report)
    assert len(texts) >= 5 and any("mg" in t for t in texts)
    results = parse_report(REPORT)
    shown = repr(results)
    assert not [t for t in texts if t in shown]
    implications = [i for r in results for d in r.drugs for i in d.implications]
    assert implications and not [i for i in implications if DOSING.search(i)]
    # The DPWG atomoxetine implication speaks of the dose requirement: dropped.
    (atomoxetine,) = [d for r in results for d in r.drugs if d.drug == "atomoxetine"]
    assert atomoxetine.implications == ()

    def strip(value: Any) -> Any:
        if isinstance(value, dict):
            return {
                k: strip(v)
                for k, v in value.items()
                if k not in ("drugRecommendation", "dosingInformation")
            }
        if isinstance(value, list):
            return [strip(v) for v in value]
        return value

    stripped = tmp_path / "stripped.json"
    stripped.write_text(json.dumps(strip(report)))
    assert parse_report(stripped) == results


# --- slice 3f: gene spans, phase sets and UGT1A1*28 -----------------------------------------


def test_gene_regions_from_px(tmp_path: Path) -> None:
    """One span per gene `INFO/PX` names, from its first to its last position ± the margin;
    a position of interest without PX defines no gene."""
    path = tmp_path / "positions.vcf.bgz"
    rows = [
        "chr10\t94761900\trs12248560\tC\tT\t.\tPASS\tPX=CYP2C19\tGT\t0/0",
        "chr10\t94645745\trs12777823\tG\tA\t.\tPASS\tPOI\tGT\t0/0",
        "chr10\t94852738\trs4244285\tG\tA\t.\tPASS\tPX=CYP2C19\tGT\t0/0",
        "chr6\t18130687\trs1142345\tT\tC,G\t.\tPASS\tPX=TPMT\tGT\t0/0",
    ]
    header = POSITIONS.split("#CHROM")[0] + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
    text = header + "\tINFO\tFORMAT\tPharmCAT\n" + "\n".join(rows) + "\n"
    path.write_bytes(gzip.compress(text.encode()))
    assert pharmcat.gene_regions(path) == [
        ("CYP2C19", "chr10", 94760900, 94853738),
        ("TPMT", "chr6", 18129687, 18131687),
    ]
    assert pharmcat.gene_regions(path, margin=0)[1] == (
        "TPMT",
        "chr6",
        18130687,
        18130687,
    )


UGT1A1_ROW = (
    "chr2\t233760233\trs3064744\tCAT\tC,CATAT,CATATAT\t.\tPASS\tPX=UGT1A1\tGT\t0/0\n"
)


def test_ugt1a1_28_read_by_position(tmp_path: Path) -> None:
    """UGT1A1*28's three ALTs are the positions file's at rs3064744; a file that spells
    them otherwise is refused."""
    assert pharmcat.UGT1A1_28 == (
        ("chr2", 233760233, "CAT", "C"),
        ("chr2", 233760233, "CAT", "CATAT"),
        ("chr2", 233760233, "CAT", "CATATAT"),
    )
    path = tmp_path / "positions.vcf.bgz"
    path.write_bytes(gzip.compress((POSITIONS + UGT1A1_ROW).encode()))
    assert pharmcat_positions(path)[-3:] == list(pharmcat.UGT1A1_28)
    path.write_bytes(
        gzip.compress((POSITIONS + UGT1A1_ROW.replace(",CATATAT", "")).encode())
    )
    with pytest.raises(ValueError, match="rs3064744"):
        pharmcat_positions(path)


def test_vcf_writes_phase_sets(positions: Path, tmp_path: Path) -> None:
    """A heterozygous biallelic row whose allele has a phase set is written `GT:PS` with
    its ALT on the haplotype WhatsHap put it on; every other row keeps `GT` with `.` for
    PS, and the header declares PS. Without phases the VCF is as before."""
    het_x = ("chrX", 154532046, "A", "C")
    calls = {
        CYP2C19_17: HET,
        CYP2C9_DEL: _called("1", "0/1"),
        TPMT_C: _called("2", "1/1"),
        TPMT_G: HOM_REF,
        het_x: HET,
        G6PD_DEL: HOM_REF,
    }
    levels = dict.fromkeys(calls, Level.HIGH)
    phases = {CYP2C19_17: (100, 1), CYP2C9_DEL: (100, 0), TPMT_C: (100, 1)}
    out = write_pharmcat_vcf(
        calls,
        pharmcat_positions(positions),
        tmp_path / "p.vcf",
        male=False,
        levels=levels,
        phases=phases,
    )
    text = out.read_text()
    assert '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">' in text
    rows = _rows(out)
    assert {r[8] for r in rows} == {"GT:PS"}
    assert _genotypes(out) == {
        ("chr6", 18130687): "1/1:.",  # homozygous: not phased
        ("chr10", 94761900): "0|1:100",
        ("chr10", 94949281): "1|0:100",
        ("chrX", 154532046): "0/1:.",  # heterozygous, but no phase set
        ("chrX", 154532055): "0/0:.",
    }
    plain = write_pharmcat_vcf(
        calls,
        pharmcat_positions(positions),
        tmp_path / "q.vcf",
        male=False,
        levels=levels,
    )
    assert "ID=PS" not in plain.read_text()
    assert {r[8] for r in _rows(plain)} == {"GT"}


def test_trusted_low_grade_is_used(tmp_path: Path) -> None:
    """rs3064744 (UGT1A1*28) called heterozygous TA7 at level low: a no-call unless the
    calibration trusts it; an unknown dosage stays a no-call either way."""
    path = tmp_path / "positions.vcf.bgz"
    path.write_bytes(gzip.compress((POSITIONS + UGT1A1_ROW).encode()))
    alleles = pharmcat_positions(path)
    ta5, ta7, ta8 = pharmcat.UGT1A1_28
    other = (OTHER_ALLELE_WARNING,)
    calls = {
        ta5: _called("0", "0/1", warnings=other),
        ta7: _called("1", "0/1"),
        ta8: _called("0", "0/1", warnings=other),
    }
    levels = dict.fromkeys(pharmcat.UGT1A1_28, Level.LOW)

    def written(trusted: frozenset[Allele]) -> str:
        out = write_pharmcat_vcf(
            calls,
            alleles,
            tmp_path / "p.vcf",
            male=False,
            levels=levels,
            trusted=trusted,
        )
        return _genotypes(out)[("chr2", 233760233)]

    assert written(frozenset()) == "./."
    assert written(frozenset(pharmcat.UGT1A1_28)) == "0/2"
    calls[ta7] = NOCALL
    assert written(frozenset()) == "./."
    assert written(frozenset(pharmcat.UGT1A1_28)) == "./."
