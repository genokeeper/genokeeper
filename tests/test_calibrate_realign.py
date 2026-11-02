"""The calibration samples re-aligned from their CRAMs: the reads written back as FASTQ, counted
against the old CRAM before and after the alignment, and the record that names what the
re-alignment replaced."""

import gzip
import json
import shutil
from dataclasses import replace
from pathlib import Path

import pysam
import pytest
from fixtures.genome_mini import MINI_QC, FakeVerifyBamID, write_lab
from fixtures.grch38_mini import READ, bwa_index, write_fastq_pair, write_reference

from genokeeper import pipeline, qc, tools
from genokeeper.align import align
from genokeeper.calibrate import (
    calibration_dir,
    hg002,
    hg003_dir,
    layout,
    realign,
    sample,
)
from genokeeper.calibrate.hg002 import downsampled_pair, library_dir
from genokeeper.calibrate.realign import (
    cram_counts,
    cram_to_fastq,
    fastq_counts,
    old_record,
)
from genokeeper.calibrate.records import (
    ALIGNED,
    Record,
    on_this_reference,
    read_record,
    write_record,
)
from genokeeper.calibrate.sample import SampleRun, align_sample
from genokeeper.checksum import sha256_file
from genokeeper.config import grch38_fasta, load_grading, unmasked_fasta
from genokeeper.pipeline import IDENTITIES, RunRefused

GRADING = replace(load_grading(), qc=MINI_QC)


@pytest.fixture(autouse=True)
def _no_gate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(tools, "GATE", tmp_path / "no-gate")


def _lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A lab with both FASTAs (text the tests never decode), `GENOKEEPER_HOME` on it; the
    masked FASTA's index and the checkpoint, which `realign` asks for first, taken as
    there (the real check runs in the end-to-end test)."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    monkeypatch.setattr(hg002, "calling_inputs", lambda _lab: None)
    for fasta, seq in (
        (grch38_fasta(lab), "ACGTNNNNAC"),
        (unmasked_fasta(lab), "ACGTAC"),
    ):
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(f">chr1\n{seq}\n")
    return lab


def _legacy(lab: Path, key: str, cram: Path, cram_key: str) -> Path:
    """3a's record of `key` (no reference digest) naming `cram`; its path."""
    path = library_dir(lab, key) / ALIGNED
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {"cram": str(cram), "key": cram_key, "metrics": {"total_bases": 1.0}}
        )
    )
    return path


def _old_cram(lab: Path, work: str, cram_key: str) -> Path:
    cram = calibration_dir(lab) / work / "align" / cram_key / "sample.cram"
    cram.parent.mkdir(parents=True)
    cram.write_bytes(b"old cram")
    return cram


def _reads(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return f.read().splitlines()[0::4]


@pytest.mark.tools
def test_cram_to_fastq_round_trips(tmp_path: Path) -> None:
    """Every primary read of a CRAM `align.align` wrote comes back once, under its name, and
    the FASTQ's reads and bases are the CRAM's."""
    fasta = bwa_index(write_reference(tmp_path / "ref" / "mini.fa"))
    fq1, fq2 = write_fastq_pair(tmp_path, fasta)
    cram = tmp_path / "old" / "sample.cram"
    cram.parent.mkdir()
    align(fq1, fq2, fasta, cram, sample="sample", threads=2)
    out = tmp_path / "pair"
    pair = (out / "x_1.fq.gz", out / "x_2.fq.gz")
    assert cram_to_fastq(cram, fasta, pair, threads=2) == pair
    first, second = (_reads(p) for p in pair)
    assert len(first) + len(second) == 800
    assert first == second  # mates in one order, with no /1 or /2
    with pysam.AlignmentFile(str(cram), "rc", reference_filename=str(fasta)) as f:
        names = sorted({str(r.query_name) for r in f if not r.is_secondary})
    assert sorted(n.removeprefix("@") for n in first) == names
    assert (
        fastq_counts(pair) == cram_counts(cram, fasta, threads=2) == (800, 800 * READ)
    )
    assert sorted(p.name for p in out.iterdir()) == ["x_1.fq.gz", "x_2.fq.gz"]


def test_old_record_reads_the_legacy_shape(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    cram = _old_cram(lab, "hg002-pcrfree", "abc")
    _legacy(lab, "pcrfree", cram, "abc")
    record = old_record(lab, "pcrfree")
    assert (record.cram, record.key, record.reference_sha256) == (
        str(cram),
        "abc",
        None,
    )
    assert not on_this_reference(lab, record)


def test_old_record_refuses_a_cram_elsewhere(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    _old_cram(lab, "hg002-pcr", "abc")  # where `align` writes: there, and not named
    elsewhere = tmp_path / "typed" / "sample.cram"
    elsewhere.parent.mkdir()
    elsewhere.write_bytes(b"")
    _legacy(lab, "pcr", elsewhere, "abc")
    with pytest.raises(RunRefused, match="names the CRAM .*hg002-pcr/align/abc/sample"):
        old_record(lab, "pcr")
    # The right folder, another key's: refused too.
    _legacy(lab, "pcr", _old_cram(lab, "hg002-pcr", "other"), "abc")
    with pytest.raises(RunRefused, match="names the CRAM .*/other/sample.cram, not"):
        old_record(lab, "pcr")
    # No record at all: HG002 has nothing to start from.
    (library_dir(lab, "pcr") / ALIGNED).unlink()
    with pytest.raises(RunRefused, match="no record"):
        old_record(lab, "pcr")


def test_old_record_finds_hg003s_single_cram(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    cram = _old_cram(lab, "hg003", "k1")
    building = hg003_dir(lab) / "align" / ".building-1a2b" / "sample.cram"
    building.parent.mkdir()
    building.write_bytes(b"a killed run's")
    counted: list[Path] = []

    def counts(cram: Path, fasta: Path, *, threads: int) -> tuple[int, int]:
        counted.append(fasta)
        return 10, 1000

    monkeypatch.setattr(realign, "cram_counts", counts)
    record = old_record(lab, "hg003")
    assert record == Record(
        cram=str(cram),
        key="k1",
        metrics={"total_bases": 1000},
        gvcf=None,
        reference_sha256=sha256_file(unmasked_fasta(lab)),
        accepted=True,
        snp_f1=None,
        replaced=None,
    )
    assert counted == [unmasked_fasta(lab)]
    assert read_record(lab, "hg003") == record
    # Two CRAMs and no record: refused naming both, nothing written.
    other = _lab(tmp_path / "two", monkeypatch)
    crams = [_old_cram(other, "hg003", k) for k in ("k1", "k2")]
    with pytest.raises(RunRefused) as e:
        old_record(other, "hg003")
    assert all(str(c) in str(e.value) for c in crams)
    assert read_record(other, "hg003") is None
    shutil.rmtree(hg003_dir(other) / "align")
    with pytest.raises(RunRefused, match="no CRAM"):
        old_record(other, "hg003")


def _never(*_a: object, **_k: object) -> object:
    raise AssertionError("called")


def test_realign_already_current_is_a_no_op(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = _lab(tmp_path, monkeypatch)
    cram = _old_cram(lab, "hg002-pcrfree", "new")
    write_record(
        lab,
        "pcrfree",
        Record(
            cram=str(cram),
            key="new",
            metrics={"total_bases": 1.0},
            gvcf=None,
            reference_sha256=sha256_file(grch38_fasta(lab)),
            accepted=True,
            snp_f1=None,
            replaced=None,
        ),
    )
    for name in ("align_sample", "run_sample"):
        monkeypatch.setattr(hg002, name, _never)
    monkeypatch.setattr(realign, "cram_counts", _never)
    assert hg002.main(["realign", "pcrfree"]) == 0
    captured = capsys.readouterr()
    assert len((captured.out + captured.err).strip().splitlines()) == 1
    assert "this reference" in captured.out


def test_realign_refuses_without_the_masked_index_before_counting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Hours of conversion are not spent on a lab that could not align the pair: the
    masked FASTA's bwa index is asked for before the old CRAM is read."""
    lab = _lab(tmp_path, monkeypatch)
    monkeypatch.setattr(hg002, "calling_inputs", pipeline.calling_inputs)
    _legacy(lab, "pcr", _old_cram(lab, "hg002-pcr", "abc"), "abc")
    monkeypatch.setattr(realign, "cram_counts", _never)
    monkeypatch.setattr(realign, "cram_to_fastq", _never)
    assert hg002.main(["realign", "pcr"]) == 2
    assert "bootstrap" in capsys.readouterr().err


@pytest.mark.tools
def test_realign_count_mismatch_refuses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """A FASTQ holding other counts than the old CRAM is removed, the old record kept, and
    nothing aligned."""
    lab = tmp_path / "lab"
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    monkeypatch.setattr(hg002, "calling_inputs", lambda _lab: None)
    fasta = write_reference(unmasked_fasta(lab))
    (tmp_path / "reads").mkdir()
    fq1, fq2 = write_fastq_pair(tmp_path / "reads", bwa_index(fasta))
    cram = calibration_dir(lab) / "hg002-pcr" / "align" / "abc" / "sample.cram"
    cram.parent.mkdir(parents=True)
    align(fq1, fq2, fasta, cram, sample="sample", threads=2)
    record = _legacy(lab, "pcr", cram, "abc")
    before = record.read_bytes()
    real = cram_counts

    def one_more(cram: Path, fasta: Path, *, threads: int) -> tuple[int, int]:
        reads, bases = real(cram, fasta, threads=threads)
        return reads + 1, bases

    monkeypatch.setattr(realign, "cram_counts", one_more)
    for name in ("align_sample", "run_sample"):
        monkeypatch.setattr(hg002, name, _never)
    assert hg002.main(["realign", "pcr", "--threads", "2"]) == 2
    err = capsys.readouterr().err
    assert "800 reads and 80000 bases, the CRAM 801 and 80000" in err, err
    assert len(err.strip().splitlines()) == 1
    pair = downsampled_pair(lab, "pcr")
    assert list(pair[0].parent.iterdir()) == []
    assert record.read_bytes() == before


def _hg003_identities(lab: Path) -> tuple[Path, Path]:
    """HG003's identities as 2b's run left them: its downsampled pair's sha256 beside the
    FASTA's; the pair itself deleted."""
    ds = layout(lab)["downsampled"]
    pair = (
        ds / "E200029822_L01_UDB-451_1.fq.gz",
        ds / "E200029822_L01_UDB-451_2.fq.gz",
    )
    known = {str(p): {"sha256": "0" * 64} for p in (*pair, grch38_fasta(lab))}
    (hg003_dir(lab) / IDENTITIES).write_text(json.dumps(known))
    return pair


def test_realign_hg003_twice_is_a_no_op(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    lab = _lab(tmp_path, monkeypatch).resolve()
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    old = _old_cram(lab, "hg003", "k1")
    pair = _hg003_identities(lab)
    monkeypatch.setattr(realign, "cram_counts", lambda *_a, **_k: (4, 400))
    monkeypatch.setattr(realign, "fastq_counts", lambda _p: (4, 400))

    def to_fastq(
        cram: Path, fasta: Path, out: tuple[Path, Path], *, threads: int
    ) -> tuple[Path, Path]:
        assert (cram, fasta, out) == (old, unmasked_fasta(lab), pair)
        out[0].parent.mkdir(parents=True, exist_ok=True)
        for p in out:
            p.write_bytes(b"reads")
        return out

    asked: list[dict[str, object]] = []
    new = hg003_dir(lab) / "align" / "k2" / "sample.cram"

    def run(_paths: object, **kwargs: object) -> SampleRun:
        asked.append(kwargs)
        new.parent.mkdir()
        new.write_bytes(b"new")
        gvcf = hg003_dir(lab) / "call" / "c" / "sample.g.vcf.gz"
        gvcf.parent.mkdir(parents=True)
        gvcf.write_bytes(b"")
        return SampleRun(cram=new, key="k2", metrics={"total_bases": 400.0}, gvcf=gvcf)

    monkeypatch.setattr(realign, "cram_to_fastq", to_fastq)
    monkeypatch.setattr(hg002, "run_sample", run)
    monkeypatch.setattr(hg002, "align_sample", _never)
    assert hg002.main(["realign", "hg003"]) == 0, capsys.readouterr().err
    (kwargs,) = asked
    assert (kwargs["name"], kwargs["sex"], kwargs["from_cram"]) == (
        "hg003",
        "male",
        True,
    )
    assert (kwargs["fastq1"], kwargs["fastq2"]) == pair
    record = read_record(lab, "hg003")
    assert record is not None and on_this_reference(lab, record)
    assert (record.cram, record.key, record.accepted) == (str(new), "k2", True)
    assert record.replaced == {"cram": str(old.parent), "pair": [str(p) for p in pair]}
    capsys.readouterr()
    monkeypatch.setattr(realign, "cram_counts", _never)
    monkeypatch.setattr(hg002, "run_sample", _never)
    assert hg002.main(["realign", "hg003"]) == 0
    captured = capsys.readouterr()
    assert len((captured.out + captured.err).strip().splitlines()) == 1


def test_realign_refuses_bases_the_alignment_lost(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """The new CRAM's `total_bases` must equal the old CRAM's bases from `cram_counts`, not
    the old record's stored metrics; otherwise the old record stays."""
    lab = _lab(tmp_path, monkeypatch)
    cram = _old_cram(lab, "hg002-pcr", "abc")
    record = _legacy(lab, "pcr", cram, "abc")  # its stored total_bases is 1.0
    before = record.read_bytes()
    monkeypatch.setattr(realign, "cram_counts", lambda *_a, **_k: (4, 400))
    monkeypatch.setattr(realign, "fastq_counts", lambda _p: (4, 400))

    def to_fastq(_c: Path, _f: Path, out: tuple[Path, Path], **_k: object) -> object:
        out[0].parent.mkdir(parents=True, exist_ok=True)
        for p in out:
            p.write_bytes(b"reads")
        return out

    monkeypatch.setattr(realign, "cram_to_fastq", to_fastq)
    lost = (tmp_path / "new.cram", "k2", {"total_bases": 300.0})
    monkeypatch.setattr(hg002, "align_sample", lambda *_a, **_k: lost)
    monkeypatch.setattr(hg002, "run_sample", _never)
    assert hg002.main(["realign", "pcr"]) == 2
    err = capsys.readouterr().err
    assert "300" in err and "400" in err, err
    assert record.read_bytes() == before


@pytest.mark.tools
def test_realign_replaces_the_old_cram_and_clean_frees_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """HG002 pcr's 3a CRAM of ENA-named reads, re-aligned for real on the mini lab: its
    reads become the pair, the run's align step aligns them (the platform gate reading the
    qualities only), the bases match, and the record names what it replaced, which `clean`
    then deletes."""
    mini = write_lab(tmp_path)
    lab = mini.paths.lab
    monkeypatch.setenv("GENOKEEPER_HOME", str(lab))
    monkeypatch.setattr(sample, "load_grading", lambda: GRADING)
    monkeypatch.setattr(qc, "run_container", FakeVerifyBamID())
    monkeypatch.setattr(pipeline, "infer_sex_from_cram", lambda *_a, **_k: "male")
    for suffix in ("", ".fai"):
        shutil.copy(f"{grch38_fasta(lab)}{suffix}", f"{unmasked_fasta(lab)}{suffix}")
    # ENA's names: the run accession and a counter, the instrument's name gone.
    archived = []
    for fastq in (mini.fastq1, mini.fastq2):
        lines = gzip.decompress(fastq.read_bytes()).decode().splitlines()
        for i in range(0, len(lines), 4):
            lines[i] = f"@SRR38517578.{i // 4 + 1} {i // 4 + 1} length={READ}"
        out = tmp_path / "ena" / fastq.name
        out.parent.mkdir(exist_ok=True)
        out.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
        archived.append(out)
    cram, key, _ = align_sample(
        mini.paths,
        name="hg002-pcr",
        fastq1=archived[0],
        fastq2=archived[1],
        sex="male",
        from_cram=True,
    )
    _legacy(lab, "pcr", cram, key)
    assert hg002.main(["realign", "pcr", "--threads", "2"]) == 0, capsys.readouterr()
    record = read_record(lab, "pcr")
    assert record is not None and on_this_reference(lab, record)
    pair = downsampled_pair(lab, "pcr")
    assert record.cram != str(cram) and Path(record.cram).is_file()
    assert record.metrics["total_bases"] == 800 * READ
    assert record.replaced == {"cram": str(cram.parent), "pair": [str(p) for p in pair]}
    assert [len(_reads(p)) for p in pair] == [400, 400]
    assert _reads(pair[0])[0].startswith("@SRR38517578.")
    capsys.readouterr()
    assert hg002.main(["clean", "pcr"]) == 0
    assert "GB freed" in capsys.readouterr().out
    assert not cram.parent.exists() and not any(p.exists() for p in pair)
    assert Path(record.cram).is_file()
    assert hg002.main(["clean", "pcr"]) == 0
    assert "nothing" in capsys.readouterr().out
