import gzip
import stat
from pathlib import Path

import pysam
import pytest
from fixtures.grch38_mini import bwa_index, write_fastq_pair, write_reference

from genokeeper import tools
from genokeeper.align import Platform, PlatformRefused, align, check_platform
from genokeeper.tools import ToolFailed


@pytest.fixture(autouse=True)
def no_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")


@pytest.fixture(scope="session")
def indexed_reference(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The mini GRCh38 reference with its bwa index, built once per session; read only."""
    return bwa_index(write_reference(tmp_path_factory.mktemp("grch38") / "mini.fa"))


def test_t7_platform_accepted(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta)
    assert check_platform(fq1) == Platform(flowcell_prefix="E", quality_levels=40)


def test_binned_qualities_refused(tmp_path: Path) -> None:
    # DNBSEQ-T7+ bins base qualities to four levels; the T7 checkpoint never saw that.
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta, qualities="#+5?")
    with pytest.raises(PlatformRefused, match="4 base quality levels"):
        check_platform(fq1)


@pytest.mark.parametrize(
    ("flowcell", "prefix"), [("V300012345", "V"), ("ML1000123", "ML")]
)
def test_g400_flowcell_refused(tmp_path: Path, flowcell: str, prefix: str) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta, names=flowcell)
    with pytest.raises(PlatformRefused, match=f"flowcell prefix {prefix}"):
        check_platform(fq1)


def test_reads_from_two_flowcell_kinds_refused(tmp_path: Path) -> None:
    # A plain (not gzipped) FASTQ whose reads come from a T7 and a G400 flowcell.
    fasta = write_reference(tmp_path / "ref.fa")
    records: list[str] = []
    for kind, flowcell in (("t7", None), ("g400", "V300012345")):
        (tmp_path / kind).mkdir()
        fq1, _ = write_fastq_pair(tmp_path / kind, fasta, pairs=10, names=flowcell)
        records += gzip.decompress(fq1.read_bytes()).decode().splitlines()
    mixed = tmp_path / "mixed.fq"
    mixed.write_text("\n".join(records) + "\n")
    with pytest.raises(PlatformRefused, match="E and V"):
        check_platform(mixed)


def _archived(fastq: Path, tmp_path: Path, original: bool) -> Path:
    """`fastq` as ENA serves a run: each name replaced by the run accession and a counter,
    with the instrument's name kept as the second field when `original`."""
    lines = gzip.decompress(fastq.read_bytes()).decode().splitlines()
    for i in range(0, len(lines), 4):
        kept = f" {lines[i][1:]}" if original else f" {i // 4 + 1} length=150"
        lines[i] = f"@SRR38517577.{i // 4 + 1}{kept}"
    out = tmp_path / "archived.fq"
    out.write_text("\n".join(lines) + "\n")
    return out


def test_archive_accession_reads_the_instruments_name(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta)
    archived = _archived(fq1, tmp_path, original=True)
    assert check_platform(archived) == Platform(flowcell_prefix="E", quality_levels=40)


def test_archive_accession_without_the_instruments_name_refused(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta)
    with pytest.raises(PlatformRefused, match="flowcell prefix SRR"):
        check_platform(_archived(fq1, tmp_path, original=False))


def test_names_false_skips_the_flowcell_keeps_the_qualities(tmp_path: Path) -> None:
    """A FASTQ written from a CRAM of ENA's reads carries the accession alone
    (`@SRR38517577.12`): `names=False` skips the flowcell and keeps the quality test."""
    fasta = write_reference(tmp_path / "ref.fa")
    fq1, _ = write_fastq_pair(tmp_path, fasta)
    archived = _archived(fq1, tmp_path, original=False)
    assert archived.read_text().startswith("@SRR38517577.1 ")
    assert check_platform(archived, names=False) == Platform(
        flowcell_prefix="E", quality_levels=40
    )
    with pytest.raises(PlatformRefused, match="flowcell prefix SRR"):
        check_platform(archived)
    (tmp_path / "binned").mkdir()
    binned, _ = write_fastq_pair(tmp_path / "binned", fasta, qualities="#+5?")
    with pytest.raises(PlatformRefused, match="4 base quality levels"):
        check_platform(
            _archived(binned, tmp_path / "binned", original=False), names=False
        )


@pytest.mark.tools
def test_align_writes_sorted_markdup_cram(
    indexed_reference: Path, tmp_path: Path
) -> None:
    fq1, fq2 = write_fastq_pair(tmp_path, indexed_reference)
    out = tmp_path / "out" / "mini.cram"
    out.parent.mkdir()
    assert align(fq1, fq2, indexed_reference, out, sample="mini", threads=2) == out
    assert out.read_bytes()[:4] == b"CRAM"
    with pysam.AlignmentFile(
        str(out), "rc", reference_filename=str(indexed_reference)
    ) as cram:
        header = cram.header.to_dict()
        assert header["HD"]["SO"] == "coordinate"
        assert any("markdup" in pg.get("CL", "") for pg in header["PG"])
        assert header["RG"] == [{"ID": "mini", "SM": "mini", "PL": "DNBSEQ"}]
        reads = list(cram.fetch("chr1"))  # needs the index
    assert len(reads) == 800
    assert sum(r.is_duplicate for r in reads) == 2  # the repeated pair
    assert sorted(p.name for p in out.parent.iterdir()) == [
        "mini.cram",
        "mini.cram.crai",
    ]


@pytest.mark.tools
def test_failing_stage_raises_and_leaves_no_cram(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "unindexed.fa")  # no bwa index: bwa mem fails
    fq1, fq2 = write_fastq_pair(tmp_path, fasta)
    out = tmp_path / "out" / "mini.cram"
    out.parent.mkdir()
    with pytest.raises(ToolFailed, match="bwa") as failed:
        align(fq1, fq2, fasta, out, sample="mini", threads=1)
    assert failed.value.cmd[0] == "bwa"
    assert list(out.parent.iterdir()) == []


@pytest.mark.tools
def test_one_gate_wraps_the_whole_pipeline(
    indexed_reference: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gate_log = tmp_path / "gate.log"
    gate = tmp_path / "coreslot-gate"
    gate.write_text(
        f'#!/usr/bin/env bash\necho "$1 $2 $3" >> "{gate_log}"\nshift 2\nexec "$@"\n'
    )
    gate.chmod(gate.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setattr(tools, "GATE", gate)
    fq1, fq2 = write_fastq_pair(tmp_path, indexed_reference)
    out = tmp_path / "mini.cram"
    align(fq1, fq2, indexed_reference, out, sample="mini", threads=2)
    # One slot for bwa | fixmate | sort | markdup and the index: gating each stage of a
    # pipe makes the later ones wait for the slot the first holds while it blocks on them.
    assert gate_log.read_text().splitlines() == ["dna:align -- bash"]
