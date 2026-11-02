"""Supervised admixture proportions by fastmixture 1.3.0 (slice 3f, "Admixture").

The panel's reference (`typing_files.build_panel_reference`) is its unrelated samples'
LD-pruned SNPs (`pruned.bim`), fastmixture's supervised P file with the superpopulations as
labels 1..K (AMR unlabelled) and `groups.json` naming each label. The user's genotypes at
the pruned SNPs (`pruned_alleles`) are written by `pgsc.write_target_vcf` (reference rows
explicit, an unknown call `./.`: the projection needs every SNP of the P file), turned into a
bed by one plink2 call in pgsc_calc's image (`to_bed`; variant IDs `chrom:pos:ref:alt`, as
the panel's), and projected on the P file (`run_projection`). Both containers run offline,
gated `dna:ancestry`.
"""

import json
from collections.abc import Sequence
from pathlib import Path

from genokeeper.containers import PLATFORMS, run_container
from genokeeper.genotype.gvcf import Allele
from genokeeper.typing_files import LABEL, SEED

PROJECTION = "projection"


def pruned_alleles(bim: Path) -> list[Allele]:
    """The pruned SNPs of the panel's `pruned.bim` in its order, as `(chr<N>, pos, REF,
    ALT)`: plink2 writes ALT as A1 and REF as A2."""
    alleles: list[Allele] = []
    for line in bim.read_text().splitlines():
        chrom, _id, _cm, pos, a1, a2 = line.split("\t")
        alleles.append((f"chr{chrom.removeprefix('chr')}", int(pos), a2, a1))
    return alleles


def read_groups(path: Path) -> tuple[str, ...]:
    """`groups.json`'s labelled groups in label order (1..K)."""
    groups: dict[str, int] = json.loads(path.read_text())["groups"]
    return tuple(sorted(groups, key=groups.__getitem__))


def read_q(path: Path, groups: Sequence[str]) -> dict[str, float]:
    """fastmixture's one-sample Q file, column k the k-th of `groups`; `ValueError` for
    other than one row of `len(groups)` columns."""
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(rows) != 1 or len(rows[0]) != len(groups):
        width = len(rows[0]) if rows else 0
        raise ValueError(
            f"{path}: {len(rows)} rows of {width} columns, not one of {len(groups)}"
        )
    return dict(zip(groups, (float(v) for v in rows[0]), strict=True))


def to_bed(vcf: Path, out: Path) -> Path:
    """`vcf` (the user's calls at `pruned_alleles`) as the PLINK bed `out` with its bim
    and fam, by plink2 in pgsc_calc's image; `out`."""
    out.parent.mkdir(parents=True, exist_ok=True)
    run_container(
        "pgsc_plink2",
        [
            "plink2", "--vcf", f"/in/{vcf.name}", "--set-all-var-ids", "@:#:$r:$a",
            "--make-bed", "--out", f"/out/{out.stem}",
        ],
        mounts=[(vcf.parent, "/in", True), (out.parent, "/out", False)],
        platform=PLATFORMS["pgsc_plink2"],
        label=LABEL,
    )  # fmt: skip
    if not out.is_file():
        raise RuntimeError(f"plink2 wrote no {out.name}")
    return out


def run_projection(bed: Path, p_file: Path, k: int, out: Path) -> dict[str, float]:
    """The user's `bed` projected on the supervised `p_file` (K `k`, seed `SEED`) by
    fastmixture, its Q file written into `out`; the proportions labelled by the
    `groups.json` beside `p_file`. `ValueError` when `k` is not that file's K."""
    groups = read_groups(p_file.parent / "groups.json")
    if len(groups) != k:
        raise ValueError(
            f"K {k}, but {p_file.parent / 'groups.json'} has {len(groups)}"
        )
    out.mkdir(parents=True, exist_ok=True)
    run_container(
        "fastmixture",
        [
            "fastmixture", "--bfile", f"/in/{bed.stem}", "--K", str(k),
            "--seed", str(SEED), "--threads", "1",
            "--projection", f"/ref/{p_file.name}", "--out", f"/out/{PROJECTION}",
        ],
        mounts=[(bed.parent, "/in", True), (p_file.parent, "/ref", True), (out, "/out", False)],
        platform=PLATFORMS["fastmixture"],
        label=LABEL,
    )  # fmt: skip
    q = out / f"{PROJECTION}.K{k}.s{SEED}.Q"
    if not q.is_file():
        raise RuntimeError(f"fastmixture wrote no {q.name}")
    return read_q(q, groups)
