"""The GRC false-duplication mask, read from two FASTAs that differ only where the masked
one holds N; the region files every slice adds, none of which may lie in a masked block."""

import json
from pathlib import Path

import pytest

from genokeeper import masking
from genokeeper.masking import (
    WINDOWS,
    load_mask,
    masked_intervals,
    overlaps,
    region_resources,
)
from genokeeper.repeats.expansionhunter import REGION_EXTENSION


def _fasta(path: Path, contigs: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in contigs.items()))
    return path


def _pair(
    tmp_path: Path, unmasked: dict[str, str], masked: dict[str, str]
) -> tuple[Path, Path]:
    return (
        _fasta(tmp_path / "unmasked.fa", unmasked),
        _fasta(tmp_path / "masked.fa", masked),
    )


THREE = {name: "ACGTACGTAC" for name in ("chr1", "chr2", "chr3")}


def test_masked_intervals_found(tmp_path: Path) -> None:
    masked = {**THREE, "chr2": "AC" + "N" * 5 + "TAC"}  # chr2:3-7, 1-based
    assert masked_intervals(*_pair(tmp_path, THREE, masked)) == [("chr2", 2, 7)]


def test_existing_n_is_not_a_mask(tmp_path: Path) -> None:
    both = {**THREE, "chr1": "ACGNACGTAC"}
    assert masked_intervals(*_pair(tmp_path, both, both)) == []


def test_substitution_refused(tmp_path: Path) -> None:
    changed = {**THREE, "chr1": "ACGTCCGTAC"}  # the A at offset 4 is a C
    with pytest.raises(ValueError, match="chr1:4"):
        masked_intervals(*_pair(tmp_path, THREE, changed))


def test_other_contigs_refused(tmp_path: Path) -> None:
    longer = {**THREE, "chr3": "ACGTACGTACG"}
    with pytest.raises(ValueError, match="chr3"):
        masked_intervals(*_pair(tmp_path, THREE, longer))


def test_masked_intervals_ignore_case(tmp_path: Path) -> None:
    unmasked, masked = _pair(tmp_path, {"chr1": "acgtACGT"}, {"chr1": "ACGTNNGT"})
    assert masked_intervals(unmasked, masked) == [("chr1", 4, 6)]


def test_masked_intervals_merge_across_windows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(masking, "WINDOW", 4)
    unmasked, masked = _pair(tmp_path, {"chr1": "ACGTACGTAC"}, {"chr1": "ACGNNNNNAC"})
    assert masked_intervals(unmasked, masked) == [("chr1", 3, 8)]


def test_overlaps_names_each_region() -> None:
    mask = [("chr21", 5_010_000, 5_166_246)]
    regions = {"smn": [WINDOWS["smn"]], "x": [("chr21", 5_100_000, 5_100_001)]}
    assert overlaps(mask, regions) == {"smn": [], "x": mask}


def test_region_resources_reads_both_shapes(tmp_path: Path) -> None:
    root = tmp_path / "resources"
    (root / "repeats").mkdir(parents=True)
    catalog = [{"LocusId": "HTT", "ReferenceRegion": "chr4:3074876-3074933"}]
    (root / "repeats" / "a.GRCh38.json").write_text(json.dumps(catalog))
    (root / "b.GRCh38.tsv").write_text("# a comment\nchrom\tstart\tend\nchr1\t10\t20\n")
    (root / "c.GRCh38.tsv").write_text("gene\tchrom\tpos\nG\tchr2\t5\n")
    (root / "spans.tsv").write_text("chrom\tstart\tend\nchr3\t0\t9\n")
    found = region_resources(root)
    assert found == {
        str(root / "repeats" / "a.GRCh38.json"): [
            ("chr4", 3074876 - REGION_EXTENSION, 3074933 + REGION_EXTENSION)
        ],
        str(root / "b.GRCh38.tsv"): [("chr1", 10, 20)],
        str(root / "c.GRCh38.tsv"): [("chr2", 4, 5)],
    }
    (root / "d.GRCh38.tsv").write_text("gene\tscore\nG\t1\n")
    with pytest.raises(ValueError, match=r"d\.GRCh38\.tsv"):
        region_resources(root)


def test_verify_refuses_another_mask(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    masked = {**THREE, "chr2": "AC" + "N" * 5 + "TAC"}
    unmasked, masked_fa = _pair(tmp_path, THREE, masked)
    out = tmp_path / "mask.bed"
    monkeypatch.setattr(masking, "load_mask", lambda: [("chr2", 2, 8)])
    assert masking.main(["verify", str(unmasked), str(masked_fa), str(out)]) == 2
    assert "differ" in capsys.readouterr().err
    assert not out.exists()
    monkeypatch.setattr(masking, "load_mask", lambda: [("chr2", 2, 7)])
    assert masking.main(["verify", str(unmasked), str(masked_fa), str(out)]) == 0
    assert out.read_text() == "chr2\t2\t7\n"
    assert "1 masked intervals, 5 bp" in capsys.readouterr().out


def test_no_region_resource_touches_the_mask() -> None:
    """The CI guard every later slice runs: no region file under resources/ and no stage-3
    window lies in a masked block, where it would read zero depth."""
    mask = load_mask()
    assert mask  # the committed mask is not empty
    regions = {
        **{name: [window] for name, window in WINDOWS.items()},
        **region_resources(),
    }
    assert {
        "src/genokeeper/resources/repeats/stage3a.GRCh38.json",
        "src/genokeeper/resources/deletions/exons.GRCh38.tsv",
    } <= set(regions)
    assert all(regions.values())  # a scan that read nothing fails
    touched = {name: hit for name, hit in overlaps(mask, regions).items() if hit}
    assert touched == {}
