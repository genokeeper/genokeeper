"""The gnomAD filtering allele frequency table: read from a miniature written here, and built
from a miniature joint sites VCF that a fake bucket serves (`fixtures.gnomad_mini`)."""

import contextlib
import json
import threading
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from email.message import Message
from pathlib import Path
from typing import Any

import duckdb
import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.gnomad_mini import (
    CHROMS,
    Bucket,
    Site,
    base,
    md5_b64,
    other,
    snv,
    vcf_bytes,
    write_reference,
)

from genokeeper import tools
from genokeeper.checksum import sha256_file
from genokeeper.cli import main
from genokeeper.config import Grading, load_grading, load_paths, unmasked_fasta
from genokeeper.evidence import frequency
from genokeeper.evidence.frequency import (
    GNOMAD_CHROMS,
    SPAN_PAD,
    FafMissing,
    FafTable,
    StreamMismatch,
    TwoTier,
    _commands,
    _record,
    build_chromosome,
    faf_file,
    faf_floor,
    faf_folder,
    gnomad_url,
    metadata_url,
    object_metadata,
    stale_chromosomes,
    two_tier,
)
from genokeeper.pipeline import RunRefused, run_genome
from genokeeper.variants.key import make_vid

ROWS = {
    "GRCh38:chr1:1000:A:G": 0.02,
    "GRCh38:chrX:5000:C:T": 0.0003,
    "GRCh38:chrY:7000:G:A": 0.06,
}


def _grading(bs1_by_gene: dict[str, float]) -> Grading:
    """A Grading built in code: BA1 0.05, BS1 0.01 and these genes' BS1s; the sections
    without defaults are the shipped ones."""
    shipped = load_grading()
    return Grading(
        ba1=0.05,
        bs1=0.01,
        bs1_by_gene=bs1_by_gene,
        paralogs=shipped.paralogs,
        lpa=shipped.lpa,
    )


def test_gnomad_chromosomes_are_the_24_files() -> None:
    assert GNOMAD_CHROMS == (*(f"chr{n}" for n in range(1, 23)), "chrX", "chrY")


def test_the_table_lives_under_the_lab(tmp_path: Path) -> None:
    assert faf_folder(tmp_path) == tmp_path / "gnomad" / "4.1"


def test_open_refuses_incomplete_manifest(tmp_path: Path) -> None:
    with pytest.raises(FafMissing, match="manifest"):
        FafTable.open(tmp_path / "absent")
    folder = write_faf_table(tmp_path / "t", ROWS, chroms=GNOMAD_CHROMS[:-1])
    with pytest.raises(FafMissing, match="chrY"):
        FafTable.open(folder)
    folder = write_faf_table(tmp_path / "t", ROWS)
    (folder / "faf-chr7.parquet").unlink()
    with pytest.raises(FafMissing, match="faf-chr7.parquet"):
        FafTable.open(folder)


@pytest.mark.parametrize(
    ("key", "value"),
    [("release", "4.0"), ("field", "AF"), ("floor", "0.0001"), ("chromosomes", [])],
)
def test_open_refuses_a_manifest_of_another_shape(
    tmp_path: Path, key: str, value: object
) -> None:
    folder = write_faf_table(tmp_path / "t", ROWS)
    manifest = folder / "manifest.json"
    data = json.loads(manifest.read_text())
    data[key] = value
    manifest.write_text(json.dumps(data))
    with pytest.raises(FafMissing, match=key):
        FafTable.open(folder)
    manifest.write_text("{not json")
    with pytest.raises(FafMissing, match="manifest"):
        FafTable.open(folder)


def test_lookup_returns_only_present_vids(tmp_path: Path) -> None:
    table = FafTable.open(write_faf_table(tmp_path / "t", ROWS))
    asked = ["GRCh38:chr1:1000:A:G", "GRCh38:chrY:7000:G:A", "GRCh38:chr2:1:A:C"]
    assert table.lookup(asked) == {
        "GRCh38:chr1:1000:A:G": 0.02,
        "GRCh38:chrY:7000:G:A": 0.06,
    }
    # Another allele at a listed position is another vid.
    assert table.lookup(["GRCh38:chr1:1000:A:T"]) == {}
    assert table.lookup([]) == {}


def test_the_table_names_its_manifest(tmp_path: Path) -> None:
    folder = write_faf_table(tmp_path / "t", ROWS)
    table = FafTable.open(folder)
    assert table.sha256 == sha256_file(folder / "manifest.json")
    described = table.describe()
    assert "gnomAD v4.1" in described and "fafmax_faf95_max_joint" in described
    assert "3 variants" in described


# --- the builder ------------------------------------------------------------------------

FLOOR = 0.000001


def _sites(faf10: float = 0.02) -> list[Site]:
    return [
        snv("chr22", 10, faf10),  # kept
        snv("chr22", 20, 0.0000005),  # below the floor
        snv("chr22", 30, 0.03, filter="BOTH_FILTERED"),  # filtered in both: not used
        snv("chr22", 40, 0.004, filter="EXOMES_FILTERED"),  # one side filtered: kept
        snv("chr22", 50, None),  # no FAF
        # An un-normalised deletion of one A of the run at 101-106: left-aligns to 100.
        Site("chr22", 105, "AA", "A", 0.001),
        snv("chr22", 150, 0.003),  # kept
    ]


def _bucket(tmp_path: Path, **sites: list[Site]) -> Bucket:
    tmp_path.mkdir(parents=True, exist_ok=True)
    return Bucket({c: vcf_bytes(tmp_path, c, s) for c, s in sites.items()})


def _rows(folder: Path, chrom: str) -> dict[str, tuple[float, str]]:
    with duckdb.connect() as con:
        found = con.execute(
            "SELECT vid, faf95_max, gen_anc FROM read_parquet(?)",
            [str(faf_file(folder, chrom))],
        ).fetchall()
    return {vid: (faf, anc) for vid, faf, anc in found}


def _manifest(folder: Path) -> dict[str, Any]:
    return json.loads((folder / "manifest.json").read_text())


def _assert_the_control(bucket: Bucket) -> None:
    """Every URL is a chromosome's file or its metadata, fetched by a GET with no body; the
    first request for a URL carries no header, and each later one exactly
    `Range: bytes=<the bytes the earlier responses served>-`."""
    allowed = {u for c in GNOMAD_CHROMS for u in (gnomad_url(c), metadata_url(c))}
    served = {id(r): s.delivered for r, s in bucket.log}
    received: dict[str, int] = {}
    for request in bucket.requests:
        url = request.full_url
        assert url in allowed
        assert request.get_method() == "GET"
        assert request.data is None
        if url not in received:
            assert request.header_items() == []
        else:
            assert request.header_items() == [("Range", f"bytes={received[url]}-")]
        received[url] = received.get(url, 0) + served.get(id(request), 0)


def _ranges(bucket: Bucket, chrom: str) -> list[str | None]:
    return [
        r.get_header("Range") for r, _ in bucket.log if r.full_url == gnomad_url(chrom)
    ]


@pytest.fixture
def grch38(tmp_path: Path) -> Path:
    return write_reference(tmp_path / "GRCh38.fa")


@pytest.fixture(autouse=True)
def no_wait(monkeypatch: pytest.MonkeyPatch) -> None:
    """Resume at once: the builder's pause between resumes is for the real bucket."""
    monkeypatch.setattr(frequency, "_WAIT", 0)


def test_requests_are_whole_file_gets(tmp_path: Path, grch38: Path) -> None:
    # The control's anchor: every URL comes from a chromosome name alone, which no input
    # of the user's can move; the first request for it asks for no byte range, and a
    # resume only for the bytes after those already received, which no position moves.
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    size = len(bucket.files["chr22"])
    bucket.faults["chr22"] = [(size // 3, ConnectionResetError)]
    folder = tmp_path / "faf"
    build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert [r.full_url for r in bucket.requests] == [
        metadata_url("chr22"),
        gnomad_url("chr22"),
        gnomad_url("chr22"),
    ]
    assert _ranges(bucket, "chr22") == [None, f"bytes={size // 3}-"]
    _assert_the_control(bucket)
    assert gnomad_url("chr22") == (
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
        "gnomad.joint.v4.1.sites.chr22.vcf.bgz"
    )
    assert metadata_url("chr22") == (
        "https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o/"
        "release%2F4.1%2Fvcf%2Fjoint%2Fgnomad.joint.v4.1.sites.chr22.vcf.bgz"
        "?fields=size,md5Hash"
    )
    with pytest.raises(ValueError, match="chr22_KI270731v1_random"):
        gnomad_url("chr22_KI270731v1_random")
    with pytest.raises(ValueError, match="chrM"):
        build_chromosome("chrM", folder, floor=FLOOR, fasta=grch38, opener=bucket)


def test_object_metadata_is_size_and_base64_md5(tmp_path: Path) -> None:
    bucket = _bucket(tmp_path, chr22=_sites())
    data = bucket.files["chr22"]
    assert object_metadata("chr22", opener=bucket) == (len(data), md5_b64(data))


def test_rows_below_floor_and_both_filtered_left_out(
    tmp_path: Path, grch38: Path
) -> None:
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    folder = tmp_path / "faf"
    rows = build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    kept = _rows(folder, "chr22")
    assert rows == len(kept) == 4
    snv10 = make_vid("chr22", 10, base("chr22", 10), other("chr22", 10))
    assert kept[snv10] == (0.02, "nfe")
    assert make_vid("chr22", 40, base("chr22", 40), other("chr22", 40)) in kept
    assert not [v for v in kept if int(v.split(":")[2]) in (20, 30, 50)]
    data = bucket.files["chr22"]
    assert _manifest(folder) == {
        "release": "4.1",
        "field": "fafmax_faf95_max_joint",
        "floor": FLOOR,
        "chromosomes": {
            "chr22": {
                "url": gnomad_url("chr22"),
                "size": len(data),
                "md5": md5_b64(data),
                "rows": 4,
            }
        },
    }
    assert sorted(p.name for p in folder.iterdir()) == [
        ".manifest.lock",
        "faf-chr22.parquet",
        "manifest.json",
    ]


def test_vid_is_normalised(tmp_path: Path, grch38: Path) -> None:
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    folder = tmp_path / "faf"
    build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    kept = _rows(folder, "chr22")
    assert kept[make_vid("chr22", 100, "CA", "C")] == (0.001, "nfe")
    assert make_vid("chr22", 105, "AA", "A") not in kept


def test_ref_disagreeing_with_the_fasta_is_skipped(
    tmp_path: Path, grch38: Path
) -> None:
    # A gnomAD REF on an N of NCBI's no-alt FASTA makes `bcftools norm` abort the whole
    # stream; the record is left out instead, and the chromosome still builds.
    wrong = other("chr22", 60)
    third = next(b for b in "ACGT" if b not in (wrong, base("chr22", 60)))
    bucket = _bucket(
        tmp_path / "served", chr22=[*_sites(), Site("chr22", 60, wrong, third, 0.02)]
    )
    folder = tmp_path / "faf"
    rows = build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    kept = _rows(folder, "chr22")
    assert rows == len(kept) == 4
    assert not [v for v in kept if int(v.split(":")[2]) == 60]


def test_the_norm_query_pipeline_is_gated_once(
    tmp_path: Path, grch38: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Gating both stages makes `query` wait for the slot `norm` holds while `norm` blocks
    # on the full pipe, until the gate's timeout. The gate is the one set when it runs.
    gate = tmp_path / "coreslot-gate"
    gate.write_text('#!/bin/sh\necho "$1" >> "$0.log"\nshift 2\nexec "$@"\n')
    gate.chmod(0o755)
    monkeypatch.setattr(tools, "GATE", gate)
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    build_chromosome(
        "chr22", tmp_path / "faf", floor=FLOOR, fasta=grch38, opener=bucket
    )
    assert (tmp_path / "coreslot-gate.log").read_text() == "dna:faf\n"
    norm, query = _commands("bcftools", tmp_path / "ref.fa", FLOOR)
    assert norm[:4] == [str(gate), "dna:faf", "--", "bcftools"]
    assert norm[4] == "norm"
    assert query[:2] == ["bcftools", "query"]


def test_interrupted_stream_leaves_no_file(tmp_path: Path, grch38: Path) -> None:
    # Every response breaks at the same byte: the resumes run out.
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    cut = len(bucket.files["chr22"]) // 2
    bucket.cut["chr22"] = cut
    folder = tmp_path / "faf"
    with pytest.raises(ConnectionResetError):
        build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert _ranges(bucket, "chr22") == [None] + [f"bytes={cut}-"] * frequency._RESUMES
    _assert_the_control(bucket)
    # No Parquet file, no temporary file and no manifest entry.
    assert list(folder.iterdir()) == []
    # Another process restarts the chromosome from its first byte, with no range.
    del bucket.cut["chr22"]
    bucket.requests.clear()
    bucket.log.clear()
    rows = build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert rows == 4
    assert _ranges(bucket, "chr22") == [None]
    _assert_the_control(bucket)


def test_stream_resumes_after_a_drop(tmp_path: Path, grch38: Path) -> None:
    whole = _bucket(tmp_path / "served", chr22=_sites())
    build_chromosome(
        "chr22", tmp_path / "whole", floor=FLOOR, fasta=grch38, opener=whole
    )
    bucket = Bucket(dict(whole.files))
    size = len(bucket.files["chr22"])
    # A reset, a stall the read timeout ends, and a connection closed early.
    breaks = [size // 4, size // 2, 3 * size // 4]
    bucket.faults["chr22"] = [
        (breaks[0], ConnectionResetError),
        (breaks[1], TimeoutError),
        (breaks[2], None),
    ]
    folder = tmp_path / "faf"
    rows = build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert rows == 4
    assert _rows(folder, "chr22") == _rows(tmp_path / "whole", "chr22")
    assert _manifest(folder) == _manifest(tmp_path / "whole")
    assert _ranges(bucket, "chr22") == [None, *(f"bytes={b}-" for b in breaks)]
    _assert_the_control(bucket)


@pytest.mark.parametrize(("codes", "resumed"), [((503, 429), True), ((403,), False)])
def test_a_transient_http_error_is_resumed(
    tmp_path: Path, grch38: Path, codes: tuple[int, ...], resumed: bool
) -> None:
    # 5xx and 429 are the bucket's own transient answers, met as a drop; another 4xx stays
    # a refusal. Here they answer the first resume after a reset.
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    size = len(bucket.files["chr22"])
    bucket.faults["chr22"] = [(size // 2, ConnectionResetError)]
    answers = list(codes)

    def opener(request: urllib.request.Request, timeout: float | None = None) -> Any:
        if request.get_header("Range") and answers:
            bucket.requests.append(request)
            code = answers.pop(0)
            raise urllib.error.HTTPError(request.full_url, code, "x", Message(), None)
        return bucket(request, timeout)

    folder = tmp_path / "faf"
    if resumed:
        rows = build_chromosome(
            "chr22", folder, floor=FLOOR, fasta=grch38, opener=opener
        )
        assert rows == 4
        _assert_the_control(bucket)
    else:
        with pytest.raises(StreamMismatch, match="answered 403"):
            build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=opener)
        assert list(folder.iterdir()) == []


@pytest.mark.parametrize("answer", ["ignored", "elsewhere", "other-total"])
def test_resume_with_wrong_content_range_refused(
    tmp_path: Path, grch38: Path, answer: str
) -> None:
    # A 200 for the whole file, or a 206 from another byte or of another object size,
    # would feed bytes twice or skip some: the chromosome fails instead.
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    bucket.faults["chr22"] = [(len(bucket.files["chr22"]) // 2, ConnectionResetError)]
    bucket.answer = answer
    folder = tmp_path / "faf"
    with pytest.raises(StreamMismatch, match="resuming"):
        build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert len(_ranges(bucket, "chr22")) == 2
    assert list(folder.iterdir()) == []


def test_md5_mismatch_refused(tmp_path: Path, grch38: Path) -> None:
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    folder = tmp_path / "faf"
    data = bucket.files["chr22"]
    for expected, what in (
        ((len(data), md5_b64(b"another file")), "md5"),
        ((len(data) + 1, md5_b64(data)), "size"),
    ):
        with pytest.raises(StreamMismatch, match=what):
            build_chromosome(
                "chr22",
                folder,
                floor=FLOOR,
                fasta=grch38,
                opener=bucket,
                expected=expected,
            )
    assert list(folder.iterdir()) == []


def test_changed_md5_rebuilds_chromosome(tmp_path: Path, grch38: Path) -> None:
    bucket = _bucket(tmp_path / "v1", chr21=[snv("chr21", 10, 0.5)], chr22=_sites())
    folder = tmp_path / "faf"
    for chrom in CHROMS:
        build_chromosome(chrom, folder, floor=FLOOR, fasta=grch38, opener=bucket)
    metadata = {c: object_metadata(c, opener=bucket) for c in CHROMS}
    stale = stale_chromosomes(folder, floor=FLOOR, metadata=metadata)
    assert not set(stale) & set(CHROMS)
    # The bucket now serves another chr22, whose row at 10 has another value.
    bucket.files["chr22"] = _bucket(tmp_path / "v2", chr22=_sites(0.07)).files["chr22"]
    metadata = {c: object_metadata(c, opener=bucket) for c in CHROMS}
    stale = stale_chromosomes(folder, floor=FLOOR, metadata=metadata)
    assert set(stale) & set(CHROMS) == {"chr22"}
    build_chromosome(
        "chr22",
        folder,
        floor=FLOOR,
        fasta=grch38,
        opener=bucket,
        expected=metadata["chr22"],
    )
    snv10 = make_vid("chr22", 10, base("chr22", 10), other("chr22", 10))
    assert _rows(folder, "chr22")[snv10] == (0.07, "nfe")
    manifest = _manifest(folder)
    assert manifest["chromosomes"]["chr22"]["md5"] == metadata["chr22"][1]
    # The other chromosome is untouched.
    assert manifest["chromosomes"]["chr21"]["md5"] == metadata["chr21"][1]
    snv21 = make_vid("chr21", 10, base("chr21", 10), other("chr21", 10))
    assert _rows(folder, "chr21") == {snv21: (0.5, "nfe")}


def test_stale_without_metadata_lists_what_the_folder_lacks(
    tmp_path: Path, grch38: Path
) -> None:
    folder = tmp_path / "faf"
    assert stale_chromosomes(folder, floor=FLOOR) == list(GNOMAD_CHROMS)
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    build_chromosome("chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    assert "chr22" not in stale_chromosomes(folder, floor=FLOOR)
    faf_file(folder, "chr22").unlink()
    assert "chr22" in stale_chromosomes(folder, floor=FLOOR)
    # A table built under another floor is stale everywhere.
    write_faf_table(folder, floor=0.001)
    assert stale_chromosomes(folder, floor=FLOOR) == list(GNOMAD_CHROMS)
    assert stale_chromosomes(folder, floor=0.001) == []


def test_another_floor_starts_a_new_manifest(tmp_path: Path, grch38: Path) -> None:
    bucket = _bucket(tmp_path / "served", chr21=[snv("chr21", 10, 0.5)], chr22=_sites())
    folder = tmp_path / "faf"
    build_chromosome("chr21", folder, floor=FLOOR, fasta=grch38, opener=bucket)
    build_chromosome("chr22", folder, floor=0.002, fasta=grch38, opener=bucket)
    manifest = _manifest(folder)
    assert manifest["floor"] == 0.002
    assert list(manifest["chromosomes"]) == ["chr22"]
    assert manifest["chromosomes"]["chr22"]["rows"] == 3  # 10, 40 and 150


def test_two_tier_keeps_low_rows_only_inside_listed_genes(
    tmp_path: Path, grch38: Path
) -> None:
    bucket = _bucket(tmp_path / "served", chr22=_sites())
    folder = tmp_path / "faf"
    tiers = TwoTier(genome_floor=0.01, spans={"GENE1": ("chr22", 1, 120)})
    rows = build_chromosome(
        "chr22", folder, floor=FLOOR, fasta=grch38, opener=bucket, two_tier=tiers
    )
    kept = _rows(folder, "chr22")
    # 10 (0.02) reaches the genome-wide floor; 40 and the deletion lie inside the span;
    # 150 (0.003) lies outside it.
    assert rows == 3
    assert sorted(int(v.split(":")[2]) for v in kept) == [10, 40, 100]
    assert _manifest(folder)["two_tier"] == {
        "genome_floor": 0.01,
        "spans": {"GENE1": ["chr22", 1, 120]},
    }
    assert "chr22" not in stale_chromosomes(folder, floor=FLOOR, two_tier=tiers)
    assert "chr22" in stale_chromosomes(folder, floor=FLOOR)


def test_floor_is_the_lowest_bs1() -> None:
    grading = load_grading()
    assert faf_floor(grading) == min(grading.bs1, *grading.bs1_by_gene.values())
    assert faf_floor(grading) == 0.000001


def test_every_bs1_gene_has_a_span() -> None:
    grading = load_grading()
    tiers = two_tier(grading)
    assert tiers.genome_floor == 0.01
    assert set(tiers.spans) == set(grading.bs1_by_gene)
    # NCBI Gene's GRCh38 span of BRCA1, chr17:43044295-43170327, padded.
    assert tiers.spans["BRCA1"] == (
        "chr17",
        43044295 - SPAN_PAD,
        43170327 + SPAN_PAD,
    )
    with pytest.raises(ValueError, match="NOSUCHGENE"):
        two_tier(_grading({"NOSUCHGENE": 0.001}))


def test_build_faf_check_lists_the_chromosomes_the_table_lacks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path))
    write_faf_table(faf_folder(tmp_path), chroms=GNOMAD_CHROMS[:-2])
    assert main(["build-faf", "--check"]) == 0
    assert capsys.readouterr().out == "chrX\nchrY\n"


def test_build_faf_refuses_without_the_grch38_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path))
    assert main(["build-faf", "--chrom", "chr22"]) == 2
    err = capsys.readouterr().err
    assert str(unmasked_fasta(tmp_path)) in err
    assert "scripts/bootstrap.sh" in err
    assert not faf_folder(tmp_path).exists()


def test_build_faf_takes_only_gnomad_chromosomes(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        main(["build-faf", "--chrom", "chrM"])
    assert "chrM" in capsys.readouterr().err


def test_concurrent_builds_keep_every_manifest_entry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Both writers read the manifest before either writes, unless a lock serialises them.
    read, both_read = frequency._listed, threading.Barrier(2, timeout=1)

    def listed(folder: Path, policy: dict[str, Any]) -> dict[str, Any]:
        found = read(folder, policy)
        with contextlib.suppress(threading.BrokenBarrierError):
            both_read.wait()
        return found

    monkeypatch.setattr(frequency, "_listed", listed)
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda c: _record(tmp_path, c, {}, {}), ["chr1", "chr2"]))
    assert list(_manifest(tmp_path)["chromosomes"]) == ["chr1", "chr2"]


@pytest.mark.parametrize("built", ["floor", "tiers"])
def test_a_table_that_misses_a_bs1_is_uncovered(tmp_path: Path, built: str) -> None:
    # A BS1 lowered below the table's floor, or a gene outside its tiers, would read the
    # alleles between as "no frequency", and BS1 would never fire on them.
    grading = _grading({"BRCA1": 0.0001})
    folder = write_faf_table(
        tmp_path / "faf", floor=0.001 if built == "floor" else FLOOR
    )
    if built == "tiers":
        manifest = json.loads((folder / "manifest.json").read_text())
        manifest["two_tier"] = {"genome_floor": 0.01, "spans": {}}
        (folder / "manifest.json").write_text(json.dumps(manifest))
    assert FafTable.open(folder).uncovered(grading) is not None
    assert FafTable.open(write_faf_table(tmp_path / "ok")).uncovered(grading) is None


def test_a_run_refuses_a_table_that_misses_a_bs1(tmp_path: Path) -> None:
    grading = _grading({"BRCA1": 0.0001})
    folder = write_faf_table(tmp_path / "faf", floor=0.001)
    with pytest.raises(RunRefused, match=r"floor 0\.001 .*genokeeper build-faf"):
        run_genome(
            load_paths(tmp_path),
            fastq1=tmp_path / "r1.fq.gz",
            fastq2=tmp_path / "r2.fq.gz",
            clinvar38=tmp_path / "clinvar.vcf.gz",
            grading=grading,
            faf=FafTable.open(folder),
        )
