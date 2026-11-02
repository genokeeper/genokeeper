import pytest

from genokeeper.variants.key import make_vid


def test_vid_adds_chr_and_uppercases() -> None:
    assert make_vid("1", 169549811, "c", "t") == "GRCh38:chr1:169549811:C:T"
    assert make_vid("chr1", 169549811, "C", "T") == "GRCh38:chr1:169549811:C:T"


def test_vid_names_the_mitochondrion_chrm() -> None:
    assert make_vid("MT", 3243, "A", "G") == "GRCh38:chrM:3243:A:G"


@pytest.mark.parametrize(
    ("ref", "alt"), [("A", "C,G"), ("", "C"), ("A", ""), ("A", ".")]
)
def test_vid_refuses_empty_or_multiallelic(ref: str, alt: str) -> None:
    with pytest.raises(ValueError):
        make_vid("chr1", 5, ref, alt)
