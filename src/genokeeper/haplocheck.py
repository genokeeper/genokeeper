"""haplocheck as a contamination QC on mutserve's chrM VCF (slice 3f, "haplocheck").

haplocheck 1.3.3 runs offline in its image (`containers.IMAGES["haplocheck"]`, built from
docker/) behind the core gate as `dna:mito`, on `mito.VCF_FILE`. Only `--raw` gives its full
table (major and minor haplogroups and their levels); `--out <name>` writes `<name>`,
`<name>.html` and `<name>.raw.txt`, of which the last is read. Its flag only ever demotes:
a contaminated run's heteroplasmies are unknown (`demoted`, passed as `mito.mito_called`'s
`heteroplasmy_withheld`), so a false flag costs an unassessed answer, never a wrong one; a
haplocheck that failed demotes them the same way.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

from genokeeper.containers import PLATFORMS, run_container

LABEL = "dna:mito"
OUT_NAME = "haplocheck"
RAW_SUFFIX = ".raw.txt"
# "Contamination Level" when haplocheck finds none.
NOT_DETERMINED = "ND"
_STATUS = {"YES": True, "NO": False}


@dataclass(frozen=True)
class Haplocheck:
    """haplocheck's verdict on one sample: contaminated or not, the contamination level
    (None when it gives none), and the major and minor haplogroups it read."""

    contaminated: bool
    level: float | None
    major: str | None
    minor: str | None


def read_haplocheck(path: Path) -> Haplocheck:
    """The one sample of haplocheck's `--raw` table at `path`; `ValueError` for no sample or
    several, or a contamination status other than YES or NO."""
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: need one sample, found {len(rows)}")
    row = rows[0]
    status = row["Contamination Status"]
    if status not in _STATUS:
        raise ValueError(f"{path}: contamination status {status!r}")
    level = row["Contamination Level"]
    return Haplocheck(
        contaminated=_STATUS[status],
        level=None if level in ("", NOT_DETERMINED) else float(level),
        major=row["Major Haplogroup"] or None,
        minor=row["Minor Haplogroup"] or None,
    )


def demoted(flagged: Haplocheck | None, failed: str | None) -> str | None:
    """Why every heteroplasmy is unknown on haplocheck's account: its failure `failed`, or
    its flag with the level it reports; None when it ran and found no contamination."""
    if failed is not None or flagged is None or not flagged.contaminated:
        return failed
    if flagged.level is None:
        return "haplocheck: possible contamination (level not reported)"
    return f"haplocheck: possible contamination at level {flagged.level:.3g}"


def run_haplocheck(vcf: Path, out_dir: Path) -> Haplocheck:
    """haplocheck on mutserve's `vcf`, its files written into `out_dir`; its verdict."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "haplocheck",
        ["haplocheck", "--raw", "--out", f"/out/{OUT_NAME}", f"/in/{vcf.name}"],
        mounts=[(vcf.parent, "/in", True), (out_dir, "/out", False)],
        platform=PLATFORMS["haplocheck"],
        label=LABEL,
    )
    raw = out_dir / f"{OUT_NAME}{RAW_SUFFIX}"
    if not raw.is_file():
        raise RuntimeError(f"haplocheck exited 0 without {raw.name}")
    return read_haplocheck(raw)
