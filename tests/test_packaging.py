"""The genome-ask skill is one file in the package, linked from `.claude/skills/` and carried
in the wheel as a regular file."""

import subprocess
import zipfile
from pathlib import Path

from genokeeper import config

REPO = Path(__file__).resolve().parents[1]
LINK = REPO / ".claude" / "skills" / "genome-ask" / "SKILL.md"
PACKAGED = "genokeeper/resources/skills/genome-ask/SKILL.md"


def test_skill_symlink_points_into_the_package() -> None:
    assert LINK.is_symlink()
    assert LINK.resolve() == REPO / "src" / PACKAGED


def test_wheel_carries_skill(tmp_path: Path) -> None:
    subprocess.run(
        ["uv", "build", "--wheel", "--out-dir", str(tmp_path)],
        cwd=REPO,
        check=True,
        capture_output=True,
    )
    (wheel,) = tmp_path.glob("*.whl")
    with zipfile.ZipFile(wheel) as zf:
        assert zf.read(PACKAGED) == (REPO / "src" / PACKAGED).read_bytes()


def test_skill_resource_reads_the_packaged_file() -> None:
    assert (
        config.skill_resource().read_bytes() == (REPO / "src" / PACKAGED).read_bytes()
    )
