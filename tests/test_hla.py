"""HLA-A and HLA-B drug alleles where OptiType and T1K agree: the MHC-region and unmapped
reads extracted from the CRAM, each caller run offline in its pinned image and read, and the
two calls compared at each PGx allele.

`fixtures/optitype_na11995_result.tsv` is OptiType 1.5.0's result on its own public test
reads (`test/exome/NA11995_SRR766010_*_fished.fastq` at tag v1.5.0, 1000 Genomes NA11995)
from a run of the pinned image, offline, on 2026-10-20. `fixtures/t1k_genotype.tsv` is
written in T1K 1.0.10's genotype format (gene, allele count, then per allele its equivalent
alleles, abundance and quality), as its run on those reads wrote it, with other values: the
exome-fished reads gave T1K no allele of positive quality.
"""

from pathlib import Path
from typing import Any

import pytest

from genokeeper.containers import IMAGES, PLATFORM, PLATFORMS
from genokeeper.pgx import hla
from genokeeper.pgx.hla import (
    MHC_REGION,
    PGX_HLA,
    extract_hla_reads,
    hla_outside,
    read_optitype,
    read_t1k,
    run_optitype,
    run_t1k,
)
from genokeeper.pgx.outside import outside_lines

FIXTURES = Path(__file__).parent / "fixtures"
HLA_A = ("*01:01", "*01:01")
NEGATIVE_B = "*15:02 negative; *57:01 negative; *58:01 negative"


def test_both_positive() -> None:
    """Both callers carry B*57:01: positive, agreed by two callers; PharmCAT gets the
    diplotype, from which it finds each allele's status itself."""
    opt = {"HLA-A": HLA_A, "HLA-B": ("*08:01", "*57:01")}
    t1k = {"HLA-A": HLA_A, "HLA-B": ("*57:01", "*08:01")}
    a, b = hla_outside(opt, t1k)
    assert (b.gene, b.diplotype, b.unresolved_reason, b.agreed) == (
        "HLA-B",
        "*08:01/*57:01",
        None,
        True,
    )
    assert b.phenotype == "*15:02 negative; *57:01 positive; *58:01 negative"
    assert b.callers == {"OptiType": "*08:01/*57:01", "T1K": "*57:01/*08:01"}
    assert (a.diplotype, a.phenotype) == ("*01:01/*01:01", "*31:01 negative")
    # Two-field calls that differ away from the PGx alleles still agree on them; PharmCAT
    # reads only those alleles' status from the diplotype it is given (OptiType's).
    t1k["HLA-B"] = ("*08:02", "*57:01")
    (_, b) = hla_outside(opt, t1k)
    assert (b.diplotype, b.unresolved_reason) == ("*08:01/*57:01", None)
    assert b.callers == {"OptiType": "*08:01/*57:01", "T1K": "*08:02/*57:01"}


def test_disagreement_unresolved() -> None:
    """An allele the callers do not both call positive, or both negative, is unresolved;
    the gene is then unresolved (one line per gene), the allele detail in its reason."""
    opt = {"HLA-A": HLA_A, "HLA-B": ("*08:01", "*57:01")}
    t1k = {"HLA-A": HLA_A, "HLA-B": ("*08:01", "*57:03")}
    _, b = hla_outside(opt, t1k)
    assert (b.diplotype, b.agreed) == (None, False)
    assert b.unresolved_reason == (
        "OptiType and T1K differ on *57:01 (OptiType positive, T1K negative)"
    )
    assert b.phenotype == "*15:02 negative; *57:01 unresolved; *58:01 negative"
    assert outside_lines([b]) == ""
    # T1K's allele dropped for its quality, or several equivalent alleles: not negative.
    for t1k_b in ((None, "*08:01"), ("*57:01,*57:03", "*08:01")):
        _, b = hla_outside({"HLA-B": ("*08:01", "*44:02")}, {"HLA-B": t1k_b})
        assert b.unresolved_reason is not None
        assert "*57:01 (OptiType negative, T1K unknown)" in b.unresolved_reason
    a, b = hla_outside({"HLA-B": ("*08:01", "*44:02")}, {})
    assert a.unresolved_reason == "OptiType and T1K made no call"
    assert b.unresolved_reason is not None and "T1K made no call" in b.unresolved_reason


def test_truncated_to_two_fields() -> None:
    """Both callers' alleles are compared at two fields, an expression suffix kept; T1K's
    equivalent alleles that share two fields are one; an allele of quality 0 or less is
    dropped (None), and one distinct allele is a homozygote."""
    assert read_t1k(FIXTURES / "t1k_genotype.tsv") == {
        "HLA-A": ("*01:01", "*01:01"),
        "HLA-B": ("*08:01", "*57:01"),
        "HLA-C": (None, "*07:01,*07:04"),
    }
    assert read_optitype(FIXTURES / "optitype_na11995_result.tsv") == {
        "HLA-A": ("*01:01", "*01:01"),
        "HLA-B": ("*08:01", "*57:01"),
        "HLA-C": ("*06:02", "*07:01"),
    }
    assert hla.two_fields("HLA-A*01:11N") == "*01:11N"
    assert hla.two_fields("A*02:01:01:01") == "*02:01"
    assert hla.two_fields("HLA-A*01:01:38L") == "*01:01"


def test_one_line_per_gene() -> None:
    opt = read_optitype(FIXTURES / "optitype_na11995_result.tsv")
    t1k = read_t1k(FIXTURES / "t1k_genotype.tsv")
    calls = hla_outside(opt, t1k)
    assert [c.gene for c in calls] == list(PGX_HLA) == ["HLA-A", "HLA-B"]
    assert outside_lines(calls) == "HLA-A\t*01:01/*01:01\nHLA-B\t*08:01/*57:01\n"


def _fake_container(calls: list[dict[str, Any]], fixture: Path, name: str):
    def run(
        image: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        calls.append(
            {
                "name": image,
                "args": args,
                "mounts": mounts,
                "platform": platform,
                "label": label,
            }
        )
        (host,) = [h for h, _, read_only in mounts if not read_only]
        (host / name).write_text(fixture.read_text())

    return run


def _reads(tmp_path: Path) -> tuple[Path, Path]:
    fq1, fq2 = tmp_path / "reads" / "hla_1.fq", tmp_path / "reads" / "hla_2.fq"
    fq1.parent.mkdir()
    fq1.write_text("x")
    fq2.write_text("x")
    return fq1, fq2


def test_callers_run_offline_in_their_images(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[dict[str, Any]] = []
    fq1, fq2 = _reads(tmp_path)
    fixture = FIXTURES / "optitype_na11995_result.tsv"
    monkeypatch.setattr(
        hla, "run_container", _fake_container(calls, fixture, "optitype_result.tsv")
    )
    out = tmp_path / "optitype"
    assert run_optitype(fq1, fq2, out) == read_optitype(fixture)
    (run,) = calls
    assert (run["name"], run["platform"], run["label"]) == (
        "optitype",
        PLATFORM,
        "dna:pgx",
    )
    assert run["args"][:2] == ["optitype", "run"]
    assert " ".join(run["args"][2:]) == (
        "-i /hla-in/hla_1.fq -i /hla-in/hla_2.fq --dna -o /optitype-out -p optitype"
        " --solver glpk"
    )
    calls.clear()
    fixture = FIXTURES / "t1k_genotype.tsv"
    monkeypatch.setattr(
        hla, "run_container", _fake_container(calls, fixture, "t1k_genotype.tsv")
    )
    index = tmp_path / "t1k-index" / "hla_dna_seq.fa"
    index.parent.mkdir()
    index.write_text(">x\n")
    assert run_t1k(fq1, fq2, tmp_path / "t1k", index=index) == read_t1k(fixture)
    (run,) = calls
    assert (run["name"], run["platform"], run["label"]) == ("t1k", PLATFORM, "dna:pgx")
    assert " ".join(run["args"]) == (
        "run-t1k -1 /hla-in/hla_1.fq -2 /hla-in/hla_2.fq --preset hla-wgs"
        " -f /t1k-index/hla_dna_seq.fa -o t1k --od /t1k-out"
    )
    assert (index.parent, "/t1k-index", True) in run["mounts"]
    for name in ("optitype", "t1k"):
        assert IMAGES[name].startswith(f"quay.io/biocontainers/{name}@sha256:")
        assert PLATFORMS[name] == PLATFORM


def test_mhc_and_unmapped_reads_extracted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One gated pipeline reads the MHC region and the unplaced unmapped reads, keeps
    pairs (no secondary, supplementary or duplicate records) and writes them as FASTQ,
    renamed into place only once complete."""
    commands: list[list[str]] = []

    def run_tool(args: list[str]) -> str:
        commands.append(args)
        script_args = args[args.index("--") + 1 :] if "--" in args else args
        Path(script_args[-2]).write_text("@r/1\n")
        Path(script_args[-1]).write_text("@r/2\n")
        return ""

    monkeypatch.setattr(hla, "run_tool", run_tool)
    monkeypatch.setattr(hla, "gated", lambda label, args: [label, "--", *args])
    cram, fasta = tmp_path / "s.cram", tmp_path / "g.fa"
    out = tmp_path / "hla"
    fq1, fq2 = extract_hla_reads(cram, fasta, out)
    assert (fq1, fq2) == (out / "hla_1.fq", out / "hla_2.fq")
    assert fq1.read_text() == "@r/1\n"
    (command,) = commands
    assert command[:3] == ["dna:pgx", "--", "bash"]
    script = command[4]
    assert "samtools view -u -h -M -F 0xD00" in script and '"*"' in script
    chrom, start, end = MHC_REGION
    assert f"{chrom}:{start}-{end}" in command
    assert str(cram) in command and str(fasta) in command
    assert sorted(p.name for p in out.iterdir()) == ["hla_1.fq", "hla_2.fq"]


HG002 = FIXTURES / "hg002"


def test_read_t1k_all_genes_on_hg002() -> None:
    """Every `HLA-` gene T1K typed on HG002, at two fields; a gene with no allele kept (none
    typed, or each of quality 0: HLA-DRB5) and the genes outside `HLA-` are absent. The
    default still reads only class I, for PharmCAT."""
    path = HG002 / "hg002.t1k_genotype.tsv.txt"
    hom = ("*01:01", "*01:01")
    assert hla.read_t1k(path, genes=None) == {
        "HLA-A": ("*26:01", "*01:01"),
        "HLA-B": ("*38:01", "*35:08"),
        "HLA-C": ("*04:01", "*12:03"),
        "HLA-E": ("*01:01", "*01:03"),
        "HLA-F": hom,
        "HLA-G": ("*01:01", "*01:06"),
        "HLA-H": ("*01:02", "*02:01"),
        "HLA-J": hom,
        "HLA-K": hom,
        "HLA-L": ("*01:01", "*01:03"),
        "HLA-N": hom,
        "HLA-P": hom,
        "HLA-S": ("*01:02", "*01:02"),
        "HLA-T": hom,
        "HLA-U": hom,
        "HLA-V": hom,
        "HLA-W": hom,
        "HLA-DRA": hom,
        "HLA-DRB1": ("*10:01", "*04:02"),
        "HLA-DRB4": ("*01:03", "*01:03"),
        "HLA-DQA1": ("*03:01", "*01:05"),
        "HLA-DQA2": hom,
        "HLA-DQB1": ("*03:02", "*05:01"),
        "HLA-DQB2": ("*01:02", "*01:02"),
        "HLA-DPA1": ("*01:03", "*01:03"),
        "HLA-DPA2": hom,
        "HLA-DPB1": ("*04:01", "*04:01"),
        "HLA-DPB2": ("*03:01", "*01:01"),
        "HLA-DMA": hom,
        "HLA-DMB": ("*01:01", "*01:02"),
        "HLA-DOA": hom,
        "HLA-DOB": hom,
    }
    assert read_t1k(path) == {
        "HLA-A": ("*26:01", "*01:01"),
        "HLA-B": ("*38:01", "*35:08"),
        "HLA-C": ("*04:01", "*12:03"),
    }
    # Equivalents joined and a quality-0 allele None, as for class I.
    assert hla.read_t1k(FIXTURES / "t1k_genotype.tsv", genes=None) == read_t1k(
        FIXTURES / "t1k_genotype.tsv"
    )


def test_read_t1k_skips_one_field_names(tmp_path: Path) -> None:
    """A gene T1K names at one field (`MICB*034`, as its run showed) is skipped, not an
    error, when every gene is read."""
    path = tmp_path / "t1k_genotype.tsv"
    path.write_text(
        "HLA-Z\t1\tHLA-Z*034\t50.0\t60\t.\t0\t-1\n"
        "HLA-DQB1\t2\tHLA-DQB1*02:01:01\t40.0\t30\tHLA-DQB1*05:01:01\t30.0\t20\n"
    )
    assert hla.read_t1k(path, genes=None) == {"HLA-DQB1": ("*02:01", "*05:01")}


def test_status_group_match() -> None:
    """A one-field group is positive when every equivalent name of one allele is in it,
    unknown when only some are or an allele was dropped; a two-field allele gives what
    `_status` gives."""
    assert hla.status(("*27:05", "*07:02"), "*27") == "positive"
    assert hla.status(("*27:05,*07:02", "*08:01"), "*27") == "unknown"
    assert hla.status(("*07:02", "*08:01"), "*27") == "negative"
    assert hla.status((None, "*08:01"), "*27") == "unknown"
    assert hla.status((None, "*27:05"), "*27") == "positive"
    assert hla.status(("*27:05,*27:02", "*08:01"), "*27") == "positive"
    assert hla.status(None, "*27") == "unknown"
    # `*2` is no prefix of `*27:05`: the group ends at the colon.
    assert hla.status(("*27:05", "*07:02"), "*2") == "negative"
    calls: list[hla.HlaCall | None] = [
        None,
        ("*57:01", "*08:01"),
        ("*57:01,*57:03", "*08:01"),
        (None, "*08:01"),
        ("*08:01", "*44:02"),
    ]
    for call in calls:
        assert hla.status(call, "*57:01") == hla._status(call, "*57:01")  # pyright: ignore[reportPrivateUsage]


def test_hla_outside_unchanged() -> None:
    """PharmCAT's HLA lines on HG002's two callers: every PGx allele negative, agreed, and
    OptiType's diplotype given."""
    opt = read_optitype(HG002 / "hg002.optitype_result.tsv.txt")
    t1k = read_t1k(HG002 / "hg002.t1k_genotype.tsv.txt")
    a, b = hla_outside(opt, t1k)
    assert (a.diplotype, a.phenotype, a.agreed) == (
        "*01:01/*26:01",
        "*31:01 negative",
        True,
    )
    assert (b.diplotype, b.phenotype, b.agreed) == ("*35:08/*38:01", NEGATIVE_B, True)
