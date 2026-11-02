"""Where the bootstrap puts slice 3b's inputs (`scripts/bootstrap.sh`'s `sv` and
`hg002_sv_truth` steps).

`python -m genokeeper.deletions.files layout` prints the `sv` paths (`SvFiles.required`),
one per line, for the bootstrap.
"""

import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from genokeeper.config import load_paths

GNOMAD_SV_FILE = "gnomad.v4.1.sv.sites.vcf.gz"
# 3a's truth folder (the bootstrap's hg002 step) and the v5.0q structural-variant files.
TRUTH_FOLDER = "HG002-v5.0q-GRCh38"
TRUTH_VCF = "HG002_GRCh38_v5.0q_stvar.vcf.gz"
TRUTH_BED = "HG002_GRCh38_v5.0q_stvar.benchmark.bed"


@dataclass(frozen=True)
class SvFiles:
    """gnomAD v4.1's SV sites VCF, with its `.tbi` beside it."""

    gnomad: Path

    def required(self) -> list[Path]:
        """The VCF and its index."""
        return [self.gnomad, self.gnomad.with_name(f"{self.gnomad.name}.tbi")]


def sv_files(lab: Path) -> SvFiles:
    """`<lab>/sv/`, where the bootstrap's sv step writes gnomAD-SV."""
    return SvFiles(gnomad=lab / "sv" / GNOMAD_SV_FILE)


def truth_files(lab: Path) -> tuple[Path, Path]:
    """HG002's v5.0q structural-variant VCF (its `.tbi` beside it) and benchmark BED, in
    3a's truth folder."""
    folder = lab / "giab" / TRUTH_FOLDER
    return folder / TRUTH_VCF, folder / TRUTH_BED


def main(argv: Sequence[str] | None = None) -> int:
    """`layout` (`SvFiles.required()` under $GENOKEEPER_HOME, one per line), for
    scripts/bootstrap.sh."""
    args = list(sys.argv[1:] if argv is None else argv)
    if args == ["layout"]:
        print("\n".join(map(str, sv_files(load_paths().lab).required())))
        return 0
    print("usage: python -m genokeeper.deletions.files layout", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
