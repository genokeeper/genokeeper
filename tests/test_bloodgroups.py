"""Blood groups by RBCeq2 2.4.4 (slice 3f, "Blood groups"; rulings R8, R9, R10).

`fixtures/hg002/hg002.rbceq2_*.tsv.txt` are RBCeq2's three tables on HG002's GIAB v5.0q
truth genotypes at the packaged positions, run without `--HPAs` (Task 5): 51 systems called,
ABO `ABO*O.01.01/ABO*O.01.01`. Their index header `UUID: <uuid>` changes per run."""

import json
import shutil
from pathlib import Path

import pysam
import pytest
from fixtures.run_context import run_context

from genokeeper import bloodgroups
from genokeeper.bloodgroups import (
    OUT_NAME,
    DbItem,
    Rbceq2Positions,
    SystemCall,
    load_positions,
    rbceq2_positions,
    read_rbceq2,
    uncovered,
)
from genokeeper.genotype.gvcf import NOCALL_REASON, Called, calls_at
from genokeeper.pgx.pharmcat import write_pharmcat_vcf
from genokeeper.sources import INFORMATION_NOTE
from genokeeper.sources.bloodgroups import (
    BLOODGROUP_TOOL,
    HPA_REASON,
    MNS_REASON,
    NOT_CLINICAL,
    RH_REASON,
    BloodGroupSource,
)
from genokeeper.store.db import Level, RowType

HG002 = Path(__file__).parent / "fixtures" / "hg002"
# ABO c.261delG: the GRCh38 reference holds the O allele's deletion, so A and B alleles are
# the insertion at chr9:133257521 (db.tsv `133257521_T_TC`); ABO*O.01.01 is `133257521_ref`.
ABO_261 = ("chr9", 133257521, "T", "TC")


def _called(dosage: str, *, haploid: bool = False) -> Called:
    if dosage == "unknown":
        return Called(
            "./.", "unknown", RowType.NOCALL, None, None, None, None, (NOCALL_REASON,)
        )
    gt = {"0": "0/0", "1": "0/1", "2": "1/1"}[dosage]
    if haploid:
        gt = {"0": "0", "1": "1"}[dosage]
    return Called(gt, dosage, RowType.PASS, 40, 30, None, None, ())


def _fasta(tmp_path: Path) -> Path:
    """chr9 and chrX of 60 bases each: `ACGT` repeated."""
    path = tmp_path / "ref.fa"
    path.write_text(">chr9\n" + "ACGT" * 15 + "\n>chrX\n" + "ACGT" * 15 + "\n")
    pysam.faidx(str(path))
    return path


def _hg002_out(tmp_path: Path) -> Path:
    out = tmp_path / "rbceq2"
    out.mkdir()
    for suffix in ("_geno.tsv", "_pheno_alphanumeric.tsv", "_pheno_numeric.tsv"):
        shutil.copy(HG002 / f"hg002.rbceq2{suffix}.txt", out / f"{OUT_NAME}{suffix}")
    return out


def _statements(
    tmp_path: Path, source: BloodGroupSource
) -> dict[str, dict[str, object]]:
    collected = source.collect(run_context(tmp_path))
    assert collected.tool == BLOODGROUP_TOOL == "rbceq2"
    kinds = [r.statement["kind"] for r in collected.rows]
    assert collected.input_rows == len(kinds) + sum(collected.dropped.values())
    out = {}
    for row in collected.rows:
        s = row.statement
        assert (s["topic"], s["actionable"], s["exclusion"]) == ("traits", False, False)
        assert (s["source"], s["validity"]) == ("rbceq2", "limited")
        out[str(s["statement_id"])] = s
    return out


def test_positions_cover_abo_261delG(tmp_path: Path) -> None:
    """The packaged cut holds ABO c.261delG's position (the insertion A and B alleles carry,
    and the O allele's reference-only item), and no platelet, RH or MNS position. Written,
    the reference-only item at a variant's position is left out (one row there), and one
    elsewhere takes the FASTA's base with ALT `.`."""
    items = load_positions()
    abo = {(i.chrom, i.pos, i.ref, i.alt) for i in items if i.system == "ABO"}
    assert ABO_261 in abo
    assert ("chr9", 133257521, None, None) in abo
    names = {i.system for i in items}
    assert not names & {"RHD", "RHCE", "GYPA", "GYPB", "GYP"}
    assert not any(n.startswith("HPA") for n in names)
    assert {"ABO", "KEL", "FY", "JK", "XK", "XG", "KLF"} <= names
    fake = [
        DbItem("ABO", "chr9", 10, "C", "CA"),  # chr9:10 is C in `_fasta`
        DbItem("ABO", "chr9", 10, None, None),
        DbItem("ABO", "chr9", 21, None, None),  # chr9:21 is A
    ]
    written = rbceq2_positions(_fasta(tmp_path), fake)
    assert written.alleles == (("chr9", 10, "C", "CA"), ("chr9", 21, "A", "."))
    assert written.systems == {
        ("chr9", 10): frozenset({"ABO"}),
        ("chr9", 21): frozenset({"ABO"}),
    }
    assert written.unchecked == {}


def test_ref_mismatch_leaves_the_allele_out(tmp_path: Path) -> None:
    """R14 (R9 per allele): an allele whose REF the GRCh38 FASTA does not have is not
    written and is named under its system; the system stays shown, its statement naming
    the allele in its reasons, and the other systems are untouched."""
    fake = [
        DbItem("FY", "chr9", 10, "GT", "G"),  # chr9:10-11 is CG
        DbItem("FY", "chr9", 13, "A", "G"),
        DbItem("JK", "chr9", 14, "C", "T"),
    ]
    written = rbceq2_positions(_fasta(tmp_path), fake)
    assert written.alleles == (("chr9", 13, "A", "G"), ("chr9", 14, "C", "T"))
    assert written.unchecked == {"FY": ("chr9:10:GT>G",)}
    statements = _statements(
        tmp_path,
        BloodGroupSource(
            calls={
                "FY": SystemCall("FY*02/FY*02", "Fy(a-b+)"),
                "JK": SystemCall("JK*01/JK*02", "Jk(a+b+)"),
            },
            not_run=None,
            withheld={},
            uncovered={},
            unchecked=written.unchecked,
        ),
    )
    fy, jk = statements["bloodgroup:FY"], statements["bloodgroup:JK"]
    assert (fy["kind"], jk["kind"]) == ("finding", "finding")
    assert (
        "RBCeq2 allele chr9:10:GT>G not checked: its GRCh38 REF disagrees with the"
        " reference"
    ) in str(fy["reason"])
    assert "RBCeq2 allele" not in str(jk["reason"])


def test_long_ref_named_by_its_length(tmp_path: Path) -> None:
    """A REF longer than ten bases (the long deletions at chr1:159205735 and
    chr4:88131048) is named by its length."""
    fake = [DbItem("FY", "chr9", 1, "ACGTACGTACGA", "A")]  # chr9:12 is T, not A
    written = rbceq2_positions(_fasta(tmp_path), fake)
    assert written.alleles == ()
    assert written.unchecked == {"FY": ("chr9:1:12bp>A",)}


def test_uncovered_position_written_as_nocall(tmp_path: Path) -> None:
    """R10: a position with no genotype is written `./.` by PharmCAT's writer, and its
    system is unassessed from that no-call, even though RBCeq2 (which reads `./.` as the
    reference) names a genotype for it."""
    positions = Rbceq2Positions(
        alleles=(ABO_261, ("chr9", 133255928, "C", "G"), ("chr7", 142952543, "C", "T")),
        systems={
            ("chr9", 133257521): frozenset({"ABO"}),
            ("chr9", 133255928): frozenset({"ABO"}),
            ("chr7", 142952543): frozenset({"KEL"}),
        },
        unchecked={},
    )
    calls = {
        ABO_261: _called("unknown"),
        ("chr9", 133255928, "C", "G"): _called("0"),
        ("chr7", 142952543, "C", "T"): _called("0"),
    }
    levels = dict.fromkeys(calls, Level.HIGH)
    vcf = write_pharmcat_vcf(
        calls,
        positions.alleles,
        tmp_path / "rbceq2_sites.vcf",
        male=True,
        levels=levels,
    )
    rows = [line.split("\t") for line in vcf.read_text().splitlines() if line[0] != "#"]
    assert [(r[1], r[3], r[4], r[9]) for r in rows] == [
        ("133257521", "T", "TC", "./."),
        ("133255928", "C", "G", "0/0"),
        ("142952543", "C", "T", "0/0"),
    ]
    missing = uncovered(vcf, positions)
    assert missing == {"ABO": (133257521,)}
    statements = _statements(
        tmp_path,
        BloodGroupSource(
            calls={
                "ABO": SystemCall("ABO*A1.01/ABO*A1.01", "A1"),
                "KEL": SystemCall("KEL*02/KEL*02", "K-k+"),
            },
            not_run=None,
            withheld={},
            uncovered=missing,
        ),
    )
    abo = statements["bloodgroup:ABO"]
    assert abo["kind"] == "unassessed"
    assert "133257521" in str(abo["reason"])
    assert "A1" not in str(abo["claim"])
    assert statements["bloodgroup:KEL"]["kind"] == "finding"


def test_male_chrx_haploid(tmp_path: Path) -> None:
    """A male's chrX outside the PARs is one allele (XK, chrX:37686082); XG's
    db.tsv position (chrX:2748343) lies in PAR1, so it stays diploid."""
    xk, xg = ("chrX", 37686082, "T", "G"), ("chrX", 2748343, "G", "C")
    calls = {xk: _called("1", haploid=True), xg: _called("1")}
    levels = dict.fromkeys(calls, Level.HIGH)
    vcf = write_pharmcat_vcf(
        calls, [xg, xk], tmp_path / "s.vcf", male=True, levels=levels
    )
    gts = {
        int(f[1]): f[9]
        for f in (line.split("\t") for line in vcf.read_text().splitlines())
        if f[0] == "chrX"
    }
    assert gts == {37686082: "1", 2748343: "0/1"}
    nocall = write_pharmcat_vcf(
        {xk: _called("unknown")},
        [xk],
        tmp_path / "n.vcf",
        male=True,
        levels={xk: Level.HIGH},
    )
    assert nocall.read_text().splitlines()[-1].endswith("\t.")


def test_read_rbceq2_hg002(tmp_path: Path) -> None:
    """Every system of HG002's tables, read by position whatever the UUID: ABO O, and the
    MNS genes RBCeq2 reports on its own (which the source never shows)."""
    calls = read_rbceq2(_hg002_out(tmp_path))
    assert len(calls) == 51
    assert calls["ABO"] == SystemCall("ABO*O.01.01/ABO*O.01.01", "O")
    assert calls["GYP"] == SystemCall("GYP*00/GYP*00", None)
    assert {"GYPA", "GYPB", "KLF", "XG", "XK", "FUT2"} <= set(calls)
    assert not any(name.startswith("HPA") for name in calls)
    assert all(c.genotype is not None for c in calls.values())


def test_read_rbceq2_undetermined_is_no_call(tmp_path: Path) -> None:
    out = tmp_path / "o"
    out.mkdir()
    head = "UUID: 1\tABO\tKEL\n"
    (out / f"{OUT_NAME}_geno.tsv").write_text(
        head + "s\tUndetermined/Undetermined\tKEL*02/KEL*02\n"
    )
    (out / f"{OUT_NAME}_pheno_alphanumeric.tsv").write_text(head + "s\t\tK-k+\n")
    assert read_rbceq2(out) == {
        "ABO": SystemCall(None, None),
        "KEL": SystemCall("KEL*02/KEL*02", "K-k+"),
    }


def test_withheld_system_unassessed(tmp_path: Path) -> None:
    """A system the calibration withholds is unassessed with the table's reason; a system
    RBCeq2 did not call, or a run without RBCeq2, is unassessed too, never a finding."""
    calls = read_rbceq2(_hg002_out(tmp_path))
    statements = _statements(
        tmp_path,
        BloodGroupSource(
            calls=calls, not_run=None, withheld={"KEL": "not covered: x"}, uncovered={}
        ),
    )
    assert (
        statements["bloodgroup:KEL"]["kind"],
        statements["bloodgroup:KEL"]["reason"],
    ) == (
        "unassessed",
        "not covered: x",
    )
    assert statements["bloodgroup:ABO"]["kind"] == "finding"
    missing = {k: v for k, v in calls.items() if k != "JK"}
    assert (
        _statements(
            tmp_path,
            BloodGroupSource(calls=missing, not_run=None, withheld={}, uncovered={}),
        )["bloodgroup:JK"]["kind"]
        == "unassessed"
    )
    for s in _statements(
        tmp_path,
        BloodGroupSource(
            calls=None, not_run="not calibrated", withheld={}, uncovered={}
        ),
    ).values():
        assert s["kind"] == "unassessed"


def test_rh_and_mns_never_called(tmp_path: Path) -> None:
    """RH, MNS and the platelet antigens are each one unassessed statement every run, with
    their reasons, even when RBCeq2's output holds calls of the MNS genes; no GYPA, GYPB,
    GYP or HPA-<n> statement is written."""
    calls = read_rbceq2(_hg002_out(tmp_path))
    statements = _statements(
        tmp_path, BloodGroupSource(calls=calls, not_run=None, withheld={}, uncovered={})
    )
    for sid, reason in (
        ("bloodgroup:RH", RH_REASON),
        ("bloodgroup:MNS", MNS_REASON),
        ("bloodgroup:HPA", HPA_REASON),
    ):
        assert (statements[sid]["kind"], statements[sid]["reason"]) == (
            "unassessed",
            reason,
        )
    assert not {"bloodgroup:GYPA", "bloodgroup:GYPB", "bloodgroup:GYP"} & set(
        statements
    )
    assert not any(sid.startswith("bloodgroup:HPA-") for sid in statements)
    assert set(statements) == {
        *(f"bloodgroup:{s}" for s in bloodgroups.systems()),
        "bloodgroup:RH",
        "bloodgroup:MNS",
        "bloodgroup:HPA",
    }


def test_bloodgroup_findings_are_information(tmp_path: Path) -> None:
    """A called system is a finding under traits, never actionable, its reasons the
    serology rule and `INFORMATION_NOTE` (RBCeq2's typing is not checked here); the claim
    holds the genotype and phenotype."""
    calls = read_rbceq2(_hg002_out(tmp_path))
    statements = _statements(
        tmp_path, BloodGroupSource(calls=calls, not_run=None, withheld={}, uncovered={})
    )
    abo = statements["bloodgroup:ABO"]
    assert abo["kind"] == "finding"
    assert NOT_CLINICAL in str(abo["reason"]) and INFORMATION_NOTE in str(abo["reason"])
    claim = json.loads(str(abo["claim"]))
    assert (claim["system"], claim["genotype"], claim["phenotype"]) == (
        "ABO",
        "ABO*O.01.01/ABO*O.01.01",
        "O",
    )
    findings = [s for s in statements.values() if s["kind"] == "finding"]
    assert len(findings) == len(bloodgroups.systems())
    assert all(INFORMATION_NOTE in str(s["reason"]) for s in findings)


def test_run_rbceq2_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Offline in its image, gated `dna:bloodgroups`, without `--HPAs`, `--RH` or
    `--phased`."""
    seen: dict[str, object] = {}

    def fake(name: str, args: list[str], **kwargs: object) -> None:
        seen.update(name=name, args=args, **kwargs)
        out = tmp_path / "out"
        for suffix in ("_geno.tsv", "_pheno_alphanumeric.tsv"):
            (out / f"{OUT_NAME}{suffix}").write_text("")

    monkeypatch.setattr(bloodgroups, "run_container", fake)
    vcf = tmp_path / "in" / "rbceq2_sites.vcf"
    vcf.parent.mkdir()
    vcf.write_text("")
    assert bloodgroups.run_rbceq2(vcf, tmp_path / "out") == tmp_path / "out"
    assert seen["name"] == "rbceq2" and seen["label"] == "dna:bloodgroups"
    assert seen["args"] == [
        "rbceq2", "--vcf", "/in/rbceq2_sites.vcf", "--out", "/out/rbceq2",
        "--reference_genome", "GRCh38", "--processes", "1",
    ]  # fmt: skip


GVCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##contig=<ID=chr9,length=60>
##contig=<ID=chrX,length=60>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmini
chr9\t1\t.\tA\t<*>\t0\t.\tEND=30\tGT:GQ:MIN_DP:PL\t0/0:50:20:0,51,509
chr9\t39\t.\tG\tA,<*>\t40\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t0/1:45:30:15,15,0:0.5,0:0,990,990,990,990,990
"""


def test_reference_only_positions_through_the_call_layer(tmp_path: Path) -> None:
    """A reference-only position (ALT `.`) goes through `calls_at` and PharmCAT's writer
    as a real gVCF gives it: a reference block's `0/0`, a position no record covers `./.`
    (a no-call its system is unassessed from, R10), and a variant allele at another
    position its genotype."""
    fasta = _fasta(tmp_path)
    plain = tmp_path / "s.g.vcf"
    plain.write_text(GVCF_HEADER)
    gvcf = Path(pysam.tabix_index(str(plain), preset="vcf", force=True))
    in_block = ("chr9", 10, "C", ".")
    no_record = ("chr9", 45, "A", ".")
    variant = ("chr9", 39, "G", "A")
    positions = Rbceq2Positions(
        alleles=(in_block, variant, no_record),
        systems={
            ("chr9", 10): frozenset({"ABO"}),
            ("chr9", 39): frozenset({"ABO"}),
            ("chr9", 45): frozenset({"KEL"}),
        },
        unchecked={},
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    calls = calls_at(
        gvcf,
        positions.alleles,
        fasta=fasta,
        haploid=lambda _chrom, _pos: False,
        scratch=scratch,
    )
    assert (calls[in_block].dosage, calls[in_block].row_type) == (
        "0",
        RowType.REF_BLOCK,
    )
    assert calls[no_record].dosage == "unknown"
    assert calls[variant].dosage == "1"
    vcf = write_pharmcat_vcf(
        calls,
        positions.alleles,
        tmp_path / "rbceq2_sites.vcf",
        male=True,
        levels=dict.fromkeys(calls, Level.HIGH),
    )
    rows = [line.split("\t") for line in vcf.read_text().splitlines() if line[0] != "#"]
    assert [(r[1], r[3], r[4], r[9]) for r in rows] == [
        ("10", "C", ".", "0/0"),
        ("39", "G", "A", "0/1"),
        ("45", "A", ".", "./."),
    ]
    assert uncovered(vcf, positions) == {"KEL": (45,)}
