"""Slice 3b's input paths: where the bootstrap puts gnomAD-SV and HG002's v5.0q SV truth."""

import os
import subprocess
import sys
from pathlib import Path

from genokeeper.deletions.files import SvFiles, sv_files, truth_files


def test_layout_lists_required(tmp_path: Path) -> None:
    """`layout` prints the two `sv` paths, one per line, under $GENOKEEPER_HOME."""
    files = sv_files(tmp_path)
    assert files == SvFiles(gnomad=tmp_path / "sv" / "gnomad.v4.1.sv.sites.vcf.gz")
    assert files.required() == [files.gnomad, Path(f"{files.gnomad}.tbi")]
    shown = subprocess.run(
        [sys.executable, "-I", "-m", "genokeeper.deletions.files", "layout"],
        env={**os.environ, "GENOKEEPER_HOME": str(tmp_path)},
        capture_output=True,
        text=True,
        check=True,
    )
    assert shown.stdout.splitlines() == [str(p) for p in files.required()]


def test_truth_files_are_3as_folder(tmp_path: Path) -> None:
    vcf, bed = truth_files(tmp_path)
    folder = tmp_path / "giab" / "HG002-v5.0q-GRCh38"
    assert vcf == folder / "HG002_GRCh38_v5.0q_stvar.vcf.gz"
    assert bed == folder / "HG002_GRCh38_v5.0q_stvar.benchmark.bed"
