"""The lab's typing inputs (`genokeeper.typing_files`): where the bootstrap puts them, and the
two that are built from the HGDP+1kGP panel (the supervised-admixture reference and the
per-superpopulation allele frequencies for ROH)."""

import io
import json
import tarfile
from collections.abc import Sequence
from compression import zstd  # noqa: F401  (tarfile's "w:zst" needs it)
from pathlib import Path
from typing import Any

import pysam
import pytest

from genokeeper import typing_files
from genokeeper.containers import IMAGES, PLATFORMS
from genokeeper.prs.pgsc import PANEL_NAME, panel_path
from genokeeper.typing_files import (
    HHV6_UNIQUE,
    build_panel_reference,
    build_roh_af,
)
from genokeeper.typing_files import (
    typing_files as layout_of,
)

PREFIX = "GRCh38_HGDP+1kGP_ALL"
POPS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")
HMMIX_HELPERS = (
    "hg38_Outgroup_1000g_HGDP.txt",
    "hg38_ancestral.tar.gz",
    "hg38_strick_callability_mask.bed",
    "hg38_mutationrate.bed",
)


def _psam() -> str:
    """Twelve samples, two per superpopulation, as the panel's psam lays them out."""
    rows = ["#IID\tSEX\tSuperPop\tPopulation\tProject"]
    for i, pop in enumerate(POPS):
        rows += [
            f"{pop}{i}a\t1\t{pop}\tpop{i}\tHGDP",
            f"{pop}{i}b\t2\t{pop}\tpop{i}\t1kG",
        ]
    return "\n".join(rows) + "\n"


def _panel(path: Path) -> Path:
    members = {
        "GRCh37_HGDP+1kGP_ALL.psam": b"not this one",
        f"{PREFIX}.pgen": b"pgen",
        f"{PREFIX}.psam": _psam().encode(),
        f"{PREFIX}.pvar.zst": b"pvar",
        "GRCh38_HGDP+1kGP.king.cutoff.out.id": b"king",
    }
    with tarfile.open(path, mode="w:zst") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


# --- the layout -----------------------------------------------------------------------------


def test_layout_lists_required(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    files = layout_of(tmp_path)
    assert files.hhv6_fasta == tmp_path / "hhv6" / "hhv6.fa"
    assert files.hmmix == tmp_path / "hmmix" / "0.9.2"
    assert files.panel_ref == tmp_path / "ancestry"
    assert files.roh_af == tmp_path / "roh"
    required = files.required()
    assert len(required) == len(set(required))
    fasta = files.hhv6_fasta
    assert required[0] == fasta
    assert [p.name for p in required[1:6]] == [
        f"hhv6.fa.{ext}" for ext in ("amb", "ann", "bwt", "pac", "sa")
    ]
    assert all(p.parent == fasta.parent for p in required[:6])
    for name in HMMIX_HELPERS:
        assert files.hmmix / name in required
    archaic = [p for p in required if p.parent == files.hmmix / "archaic"]
    assert sorted(p.name for p in archaic) == sorted(
        f"individuals_highcov.{c}.bcf{ext}"
        for c in range(1, 23)
        for ext in ("", ".csi")
    )
    # K is the labelled groups: the six superpopulations less AMR, which is label 0.
    assert files.panel_ref / "reference.K5.s1.P" in required
    for name in ("pruned.bed", "pruned.bim", "pruned.fam", "labels.txt", "groups.json"):
        assert files.panel_ref / name in required
    for pop in POPS:
        assert files.roh_af / f"af.{pop}.tsv.gz" in required
        assert files.roh_af / f"af.{pop}.tsv.gz.tbi" in required
    assert len(required) == 6 + 4 + 44 + 1 + 5 + 2 * len(POPS)
    assert not any(" " in str(p) for p in required)  # one per line for the bootstrap
    monkeypatch.setenv("GENOKEEPER_HOME", str(tmp_path))
    assert typing_files.main(["layout"]) == 0
    assert capsys.readouterr().out.splitlines() == [str(p) for p in required]


def test_hhv6_unique_regions() -> None:
    assert HHV6_UNIQUE == {"NC_001664.4": (8090, 151288), "NC_000898.1": (8794, 153321)}


# --- ROH allele frequencies ----------------------------------------------------------------

AFREQ = (
    "#CHROM\tPOS\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
    "{c}\t30\trs3\tT\tC\t0.125\t12\n"
    "{c}\t10\trs1\tA\tG\t0.25\t12\n"
    "{c}\t20\trs2\tC\tT\tnan\t0\n"
    "{c}\t40\trs4\tG\tA\t1\t12\n"
)


@pytest.mark.parametrize("contig", ["1", "chr1"])
def test_roh_af_format(tmp_path: Path, contig: str) -> None:
    """bcftools' AF file: `chr1<TAB>pos<TAB>REF,ALT<TAB>alt frequency`, by position,
    bgzipped and indexed, whether the panel names its contigs `1` or `chr1`; a site with no
    frequency is left out."""
    afreq = tmp_path / "freq.EUR.afreq"
    afreq.write_text(AFREQ.format(c=contig))
    out = typing_files.write_af_file(afreq, tmp_path / "af.EUR.tsv.gz")
    assert out == tmp_path / "af.EUR.tsv.gz"
    with pysam.BGZFile(str(out)) as f:  # type: ignore[no-untyped-call]
        lines = f.read().decode().splitlines()
    assert lines == [
        "chr1\t10\tA,G\t0.25",
        "chr1\t30\tT,C\t0.125",
        "chr1\t40\tG,A\t1",
    ]
    with pysam.TabixFile(str(out)) as tbx:
        assert list(tbx.contigs) == ["chr1"]
        assert [r for r in tbx.fetch("chr1", 25, 35)] == ["chr1\t30\tT,C\t0.125"]


def test_roh_af_one_file_per_category(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    panel = _panel(tmp_path / PANEL_NAME)
    seen: list[dict[str, Any]] = []

    def fake(
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        (_, _, read_only), (out, _, _) = mounts
        seen.append(
            {"name": name, "args": list(args), "platform": platform, "label": label}
        )
        assert read_only
        prefix = args[args.index("--out") + 1].rsplit("/", 1)[1]
        for pop in POPS:
            (out / f"{prefix}.{pop}.afreq").write_text(AFREQ.format(c="1"))

    monkeypatch.setattr(typing_files, "run_container", fake)
    out_dir = tmp_path / "roh"
    written = build_roh_af(panel, out_dir)
    assert written == [out_dir / f"af.{pop}.tsv.gz" for pop in POPS]
    assert all(Path(f"{p}.tbi").exists() for p in written)
    (call,) = seen
    assert (call["name"], call["platform"], call["label"]) == (
        "pgsc_plink2",
        PLATFORMS["pgsc_plink2"],
        "dna:ancestry",
    )
    args = " ".join(call["args"])
    assert "--loop-cats SuperPop" in args
    assert "--freq cols=chrom,pos,ref,alt,altfreq,nobs" in args
    assert "--autosome" in args and "--max-alleles 2" in args
    # Only the AF files remain beside each other: no extracted panel, no plink2 output.
    assert sorted(p.name for p in out_dir.iterdir()) == sorted(
        [f"af.{pop}.tsv.gz{ext}" for pop in POPS for ext in ("", ".tbi")]
    )


# --- the supervised-admixture reference -----------------------------------------------------


def test_panel_reference_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """plink2 prunes and writes a bed in pgsc_calc's image, then fastmixture runs
    `--supervised` in its own: both offline, gated `dna:ancestry`; AMR is label 0 and K is
    the number of labelled groups; the label file follows the bed's `.fam` order."""
    panel = _panel(tmp_path / PANEL_NAME)
    calls: list[dict[str, Any]] = []
    fam_order = [f"{pop}{i}{s}" for i, pop in enumerate(POPS) for s in "ab"][::-1]

    def fake(
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        calls.append(
            {
                "name": name,
                "args": list(args),
                "mounts": mounts,
                "platform": platform,
                "label": label,
            }
        )
        out = next(host for host, target, ro in mounts if target == "/out" and not ro)
        if name == "pgsc_plink2" and "--indep-pairwise" in args:
            (out / "prune.prune.in").write_text("rs1\n")
        elif name == "pgsc_plink2":
            assert "--extract" in args and "--make-bed" in args
            (out / "pruned.bed").write_bytes(b"bed")
            (out / "pruned.bim").write_text("1\trs1\t0\t10\tA\tG\n")
            (out / "pruned.fam").write_text(
                "".join(f"0\t{iid}\t0\t0\t0\t-9\n" for iid in fam_order)
            )
        else:
            k, seed = args[args.index("--K") + 1], args[args.index("--seed") + 1]
            (out / f"reference.K{k}.s{seed}.P").write_text("0.5\n")
            (out / f"reference.K{k}.s{seed}.Q").write_text("1\n")

    monkeypatch.setattr(typing_files, "run_container", fake)
    out_dir = tmp_path / "ancestry"
    p_file = build_panel_reference(panel, out_dir, threads=3)
    assert p_file == out_dir / "reference.K5.s1.P"
    prune, bed, fast = calls
    assert [c["name"] for c in calls] == ["pgsc_plink2", "pgsc_plink2", "fastmixture"]
    assert all(c["label"] == "dna:ancestry" for c in calls)
    assert [c["platform"] for c in calls] == [
        PLATFORMS["pgsc_plink2"],
        PLATFORMS["pgsc_plink2"],
        PLATFORMS["fastmixture"],
    ]
    assert {"pgsc_plink2", "fastmixture"} <= IMAGES.keys()
    pruning = " ".join(prune["args"])
    assert "--indep-pairwise 1000 50 0.1 --maf 0.01 --autosome" in pruning
    assert f"--pfile /panel/{PREFIX} vzs" in pruning
    assert "--remove /panel/GRCh38_HGDP+1kGP.king.cutoff.out.id" in pruning
    assert "--extract /out/prune.prune.in" in " ".join(bed["args"])
    assert "--out /out/pruned" in " ".join(bed["args"])
    args = fast["args"]
    assert args[0] == "fastmixture"
    for option, value in (
        ("--bfile", "/out/pruned"),
        ("--supervised", "/out/labels.txt"),
        ("--seed", "1"),
        ("--K", "5"),
        ("--threads", "3"),
        ("--out", "/out/reference"),
    ):
        assert args[args.index(option) + 1] == value
    # Labels follow the .fam order; AMR is 0 and the other five groups are 1..5.
    groups = json.loads((out_dir / "groups.json").read_text())
    assert groups == {
        "K": 5,
        "groups": {"AFR": 1, "CSA": 2, "EAS": 3, "EUR": 4, "MID": 5},
        "unlabelled": ["AMR"],
    }
    by_pop = {"AMR": 0, **groups["groups"]}
    labels = (out_dir / "labels.txt").read_text().splitlines()
    assert labels == [str(by_pop[iid.rstrip("0123456789ab")]) for iid in fam_order]
    assert sorted(p.name for p in out_dir.iterdir()) == sorted(
        [
            "pruned.bed",
            "pruned.bim",
            "pruned.fam",
            "labels.txt",
            "groups.json",
            "reference.K5.s1.P",
            "reference.K5.s1.Q",
        ]
    )


def test_panel_reference_refuses_other_groups(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A panel whose superpopulations are not the six this module names is refused before
    anything runs."""
    panel = tmp_path / PANEL_NAME
    members = {
        f"{PREFIX}.pgen": b"pgen",
        f"{PREFIX}.psam": b"#IID\tSEX\tSuperPop\nx\t1\tAFR\ny\t1\tOCE\n",
        f"{PREFIX}.pvar.zst": b"pvar",
        "GRCh38_HGDP+1kGP.king.cutoff.out.id": b"king",
    }
    with tarfile.open(panel, mode="w:zst") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))

    def boom(*args: object, **kwargs: object) -> None:
        raise AssertionError("ran a container")

    monkeypatch.setattr(typing_files, "run_container", boom)
    with pytest.raises(ValueError, match="superpopulations"):
        build_panel_reference(panel, tmp_path / "ancestry", threads=1)


def test_cli_builds_from_the_labs_panel(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[tuple[str, Path, Path]] = []
    monkeypatch.setattr(
        typing_files,
        "build_panel_reference",
        lambda panel, out, *, threads: seen.append(("panel", panel, out)) or out,
    )
    monkeypatch.setattr(
        typing_files,
        "build_roh_af",
        lambda panel, out: seen.append(("af", panel, out)) or [],
    )
    assert typing_files.main(["build-panel", str(tmp_path)]) == 0
    assert typing_files.main(["build-af", str(tmp_path)]) == 0
    assert seen == [
        ("panel", panel_path(tmp_path), tmp_path / "ancestry"),
        ("af", panel_path(tmp_path), tmp_path / "roh"),
    ]
    assert typing_files.main(["nonsense"]) == 2
