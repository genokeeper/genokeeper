"""haplocheck as a contamination QC on mutserve's chrM VCF (slice 3f, "haplocheck").

`fixtures/hg002/hg002.haplocheck.txt` is haplocheck 1.3.3's `--raw` table on HG002's
mutserve VCF (Task 5): 18 quoted columns, not contaminated, level `ND`."""

from collections.abc import Sequence
from pathlib import Path

import pytest

from genokeeper import haplocheck
from genokeeper.containers import PLATFORMS
from genokeeper.haplocheck import Haplocheck, demoted, read_haplocheck, run_haplocheck

HG002 = Path(__file__).parent / "fixtures" / "hg002" / "hg002.haplocheck.txt"


def test_read_haplocheck_hg002() -> None:
    assert read_haplocheck(HG002) == Haplocheck(
        contaminated=False, level=None, major="H5a7", minor="H5a7"
    )


def test_read_haplocheck_contaminated(tmp_path: Path) -> None:
    lines = HG002.read_text().splitlines()
    fields = lines[1].split("\t")
    fields[1:3] = ['"YES"', '"0.052"']
    fields[11] = '"U5a1"'
    path = tmp_path / "haplocheck.raw.txt"
    path.write_text(lines[0] + "\n" + "\t".join(fields) + "\n")
    assert read_haplocheck(path) == Haplocheck(
        contaminated=True, level=0.052, major="H5a7", minor="U5a1"
    )


def test_demoted_reason() -> None:
    """Why every heteroplasmy is unknown: haplocheck's failure, else its flag with the
    level it reports (never "nan"), else nothing."""
    flag = Haplocheck(True, 0.052, "H5a7", "U5a1")
    assert demoted(flag, None) == "haplocheck: possible contamination at level 0.052"
    unleveled = Haplocheck(True, None, "H5a7", "U5a1")
    assert demoted(unleveled, None) == (
        "haplocheck: possible contamination (level not reported)"
    )
    assert demoted(None, "haplocheck failed: exit 1") == "haplocheck failed: exit 1"
    assert demoted(read_haplocheck(HG002), None) is None


def test_read_haplocheck_refuses_other_shapes(tmp_path: Path) -> None:
    lines = HG002.read_text().splitlines()
    path = tmp_path / "haplocheck.raw.txt"
    path.write_text(lines[0] + "\n")
    with pytest.raises(ValueError, match="one sample"):
        read_haplocheck(path)
    path.write_text("\n".join([lines[0], lines[1], lines[1]]) + "\n")
    with pytest.raises(ValueError, match="one sample"):
        read_haplocheck(path)
    path.write_text(lines[0] + "\n" + lines[1].replace('"NO"', '"MAYBE"') + "\n")
    with pytest.raises(ValueError, match="MAYBE"):
        read_haplocheck(path)


def test_run_haplocheck_offline_raw(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """haplocheck runs in its image behind the core gate as `dna:mito`, with `--raw` (the
    only way to its 18 columns) and `--out`, whose `.raw.txt` it reads."""
    seen: list[tuple[str, list[str], Sequence[tuple[Path, str, bool]], str | None]] = []

    def fake(
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), mounts, label))
        assert platform == PLATFORMS["haplocheck"]
        (out,) = [host for host, target, ro in mounts if target == "/out"]
        (out / "haplocheck.raw.txt").write_text(HG002.read_text())

    monkeypatch.setattr(haplocheck, "run_container", fake)
    vcf = tmp_path / "mito" / "mutserve.vcf.gz"
    vcf.parent.mkdir()
    vcf.write_bytes(b"x")
    out = tmp_path / "haplocheck"
    found = run_haplocheck(vcf, out)
    assert found == read_haplocheck(HG002)
    ((name, args, mounts, label),) = seen
    assert (name, label) == ("haplocheck", "dna:mito")
    assert args == [
        "haplocheck",
        "--raw",
        "--out",
        "/out/haplocheck",
        "/in/mutserve.vcf.gz",
    ]
    assert (vcf.parent, "/in", True) in mounts and (out, "/out", False) in mounts
