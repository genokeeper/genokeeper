"""Where the bootstrap puts slice 3f's inputs (`scripts/bootstrap.sh`'s `typing` step), and
the two it builds from the HGDP+1kGP reference panel.

`python -m genokeeper.typing_files layout` prints the paths (`TypingFiles.required`), one per
line, for the bootstrap; `build-panel LAB` and `build-af LAB` build the panel-derived inputs
(both CPU-heavy: each streams the 16 GB panel, so they run behind the core gate).

- `hhv6/hhv6.fa`: HHV-6A (`NC_001664.4`) and HHV-6B (`NC_000898.1`) concatenated, with its
  bwa index beside it;
- `hmmix/0.9.2/`: the hg38 outgroup, ancestral, mask and mutation-rate files (Zenodo
  11212339) and `archaic/`, the 22 archaic BCFs with their `.csi` (Zenodo 13368126);
- `ancestry/`: the panel pruned to a PLINK bed, its label file and fastmixture's
  supervised P file, the reference an ancestry projection reads;
- `roh/`: one `bcftools roh --AF-file` file per superpopulation.
"""

import json
import os
import sys
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.config import load_paths
from genokeeper.containers import PLATFORMS, run_container
from genokeeper.prs.pgsc import extract_panel_members, panel_path

# The core gate's label for the panel runs.
LABEL = "dna:ancestry"
HMMIX_VERSION = "0.9.2"
HMMIX_HELPERS = (
    "hg38_Outgroup_1000g_HGDP.txt",
    "hg38_ancestral.tar.gz",
    "hg38_strick_callability_mask.bed",
    "hg38_mutationrate.bed",
)
# Zenodo 13368126: one BCF with its index per autosome.
ARCHAIC = tuple(
    f"individuals_highcov.{chrom}.bcf{ext}"
    for chrom in range(1, 23)
    for ext in ("", ".csi")
)
BWA_INDEX_EXTENSIONS = ("amb", "ann", "bwt", "pac", "sa")
# Contig -> (first, last) 1-based inclusive: each genome's unique region, between its
# direct-repeat ends (RefSeq `repeat_region` DR_L and DR_R).
HHV6_UNIQUE: dict[str, tuple[int, int]] = {
    "NC_001664.4": (8090, 151288),
    "NC_000898.1": (8794, 153321),
}
# The panel's psam column and its labels (read 2026-10-30): the five that are labelled in
# the supervised run, and the one that is not. AMR is admixed, so it is label 0.
SUPERPOP_COLUMN = "SuperPop"
UNLABELLED = "AMR"
SUPERPOPS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")
SEED = 1
# plink2's pruning of the reference: window 1000 variants, step 50, r2 0.1.
PRUNE = ("--indep-pairwise", "1000", "50", "0.1", "--maf", "0.01", "--autosome")


def _labelled(superpops: Iterable[str]) -> list[str]:
    return sorted(set(superpops) - {UNLABELLED})


@dataclass(frozen=True)
class TypingFiles:
    """The folders of the typing inputs under the lab."""

    hhv6_fasta: Path
    hmmix: Path
    panel_ref: Path
    roh_af: Path

    def reference_p(self) -> Path:
        """fastmixture's P file for the supervised run (`K` labelled groups, seed `SEED`)."""
        return self.panel_ref / f"reference.K{len(_labelled(SUPERPOPS))}.s{SEED}.P"

    def af_file(self, pop: str) -> Path:
        return self.roh_af / f"af.{pop}.tsv.gz"

    def required(self) -> list[Path]:
        """The FASTA and its five bwa index files, hmmix's files, the panel reference
        (bed, labels, groups and P file) and each superpopulation's AF file with its index."""
        fasta = self.hhv6_fasta
        return [
            fasta,
            *(fasta.with_name(f"{fasta.name}.{ext}") for ext in BWA_INDEX_EXTENSIONS),
            *(self.hmmix / name for name in HMMIX_HELPERS),
            *(self.hmmix / "archaic" / name for name in ARCHAIC),
            *(self.panel_ref / f"pruned.{ext}" for ext in ("bed", "bim", "fam")),
            self.panel_ref / "labels.txt",
            self.panel_ref / "groups.json",
            self.reference_p(),
            *(
                path
                for pop in SUPERPOPS
                for path in (
                    self.af_file(pop),
                    self.af_file(pop).with_name(f"{self.af_file(pop).name}.tbi"),
                )
            ),
        ]


def typing_files(lab: Path) -> TypingFiles:
    return TypingFiles(
        hhv6_fasta=lab / "hhv6" / "hhv6.fa",
        hmmix=lab / "hmmix" / HMMIX_VERSION,
        panel_ref=lab / "ancestry",
        roh_af=lab / "roh",
    )


def read_superpops(psam: Path) -> dict[str, str]:
    """The psam's sample id -> `SUPERPOP_COLUMN` label. `ValueError` when the panel does not
    have exactly the superpopulations `SUPERPOPS` names (it is pinned by its md5, so this
    refuses another panel, not a drift)."""
    lines = psam.read_text().splitlines()
    header = lines[0].lstrip("#").split("\t")
    if "IID" not in header or SUPERPOP_COLUMN not in header:
        raise ValueError(f"{psam}: no IID and {SUPERPOP_COLUMN} columns")
    iid, pop = header.index("IID"), header.index(SUPERPOP_COLUMN)
    labels = {}
    for line in lines[1:]:
        fields = line.split("\t")
        labels[fields[iid]] = fields[pop]
    if set(labels.values()) != set(SUPERPOPS):
        raise ValueError(
            f"{psam}: superpopulations {sorted(set(labels.values()))},"
            f" not {list(SUPERPOPS)}"
        )
    return labels


def build_panel_reference(panel: Path, out: Path, *, threads: int) -> Path:
    """The supervised-admixture reference in `out`: the panel's unrelated samples' autosomal
    biallelic SNPs, pruned (`PRUNE`) and written as `pruned.{bed,bim,fam}` by plink2 in
    pgsc_calc's image, `labels.txt` (the `.fam` order; `UNLABELLED` is 0, the other groups
    1..K by name) and `groups.json`, then fastmixture's `--supervised` run in its image,
    which writes the P file returned. Both containers run offline, gated `LABEL`."""
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out, prefix=".panel-") as tmp:
        work = Path(tmp)
        found = extract_panel_members(panel, work)
        pops = read_superpops(work / found[".psam"])
        groups = {pop: i for i, pop in enumerate(_labelled(pops.values()), start=1)}
        prefix = found[".pgen"].removesuffix(".pgen")
        select = [
            "plink2", "--pfile", f"/panel/{prefix}", "vzs",
            "--remove", f"/panel/{found['.king.cutoff.out.id']}",
            "--max-alleles", "2", "--snps-only", "just-acgt",
        ]  # fmt: skip
        mounts = [(work, "/panel", True), (out, "/out", False)]

        def plink2(*args: str) -> None:
            run_container(
                "pgsc_plink2",
                [*select, *args],
                mounts=mounts,
                platform=PLATFORMS["pgsc_plink2"],
                label=LABEL,
            )

        plink2(*PRUNE, "--out", "/out/prune")
        plink2("--extract", "/out/prune.prune.in", "--make-bed", "--out", "/out/pruned")
    labels = []
    for line in (out / "pruned.fam").read_text().splitlines():
        iid = line.split()[1]
        if iid not in pops:
            raise ValueError(f"{iid} is in the bed and not in the panel's psam")
        labels.append(groups.get(pops[iid], 0))
    (out / "labels.txt").write_text("".join(f"{label}\n" for label in labels))
    k = len(groups)
    (out / "groups.json").write_text(
        json.dumps({"K": k, "groups": groups, "unlabelled": [UNLABELLED]}, indent=2)
        + "\n"
    )
    run_container(
        "fastmixture",
        [
            "fastmixture", "--bfile", "/out/pruned", "--K", str(k),
            "--seed", str(SEED), "--threads", str(threads),
            "--supervised", "/out/labels.txt", "--out", "/out/reference",
        ],
        mounts=[(out, "/out", False)],
        platform=PLATFORMS["fastmixture"],
        label=LABEL,
    )  # fmt: skip
    for scratch in (*out.glob("prune.*"), *out.glob("pruned.log")):
        scratch.unlink()
    p_file = out / f"reference.K{k}.s{SEED}.P"
    if not p_file.is_file():
        raise ValueError(f"fastmixture wrote no {p_file.name}")
    return p_file


def _afreq_rows(afreq: Path) -> Iterator[tuple[int, int, str, str, str]]:
    """(contig number, position, REF, ALT, ALT frequency as plink2 wrote it) per row with a
    frequency; the contig `1` or `chr1`."""
    lines = afreq.read_text().splitlines()
    header = lines[0].lstrip("#").split("\t")
    chrom, pos, ref, alt, freq = (
        header.index(c) for c in ("CHROM", "POS", "REF", "ALT", "ALT_FREQS")
    )
    for line in lines[1:]:
        fields = line.split("\t")
        if fields[freq].lower() in ("nan", "", "."):
            continue
        number = fields[chrom].removeprefix("chr")
        yield int(number), int(fields[pos]), fields[ref], fields[alt], fields[freq]


def write_af_file(afreq: Path, dest: Path) -> Path:
    """`bcftools roh --AF-file`'s file from plink2's `.afreq`: `chr<N>`, position,
    `REF,ALT` and the ALT frequency, by position; bgzipped and tabix-indexed at `dest`
    (`.tsv.gz`)."""
    rows = sorted(_afreq_rows(afreq))
    plain = dest.with_suffix("")
    plain.write_text("".join(f"chr{c}\t{p}\t{r},{a}\t{f}\n" for c, p, r, a, f in rows))
    # Compresses `plain` to `dest` and indexes it (.tbi), removing `plain`.
    pysam.tabix_index(str(plain), seq_col=0, start_col=1, end_col=1, force=True)  # type: ignore[no-untyped-call]
    return dest


def build_roh_af(panel: Path, out: Path) -> list[Path]:
    """One AF file per superpopulation in `out` (`af.<POP>.tsv.gz`, with its `.tbi`): plink2
    `--freq --loop-cats SuperPop` over the panel's unrelated samples' autosomal biallelic
    SNPs in pgsc_calc's image (offline, gated `LABEL`), written by `write_af_file`."""
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out, prefix=".panel-") as tmp:
        work = Path(tmp)
        found = extract_panel_members(panel, work)
        prefix = found[".pgen"].removesuffix(".pgen")
        run_container(
            "pgsc_plink2",
            [
                "plink2", "--pfile", f"/panel/{prefix}", "vzs",
                "--remove", f"/panel/{found['.king.cutoff.out.id']}",
                "--autosome", "--max-alleles", "2", "--snps-only", "just-acgt",
                "--freq", "cols=chrom,pos,ref,alt,altfreq,nobs",
                "--loop-cats", SUPERPOP_COLUMN, "--out", "/out/.freq",
            ],
            mounts=[(work, "/panel", True), (out, "/out", False)],
            platform=PLATFORMS["pgsc_plink2"],
            label=LABEL,
        )  # fmt: skip
        written = []
        for pop in SUPERPOPS:
            afreq = out / f".freq.{pop}.afreq"
            if not afreq.is_file():
                raise ValueError(f"plink2 wrote no frequencies for {pop}")
            written.append(write_af_file(afreq, out / f"af.{pop}.tsv.gz"))
            afreq.unlink()
    for scratch in out.glob(".freq.*"):
        scratch.unlink()
    return written


def main(argv: Sequence[str] | None = None) -> int:
    """`layout` (`TypingFiles.required()` under $GENOKEEPER_HOME, one per line), `build-panel
    LAB` and `build-af LAB`, for scripts/bootstrap.sh."""
    args = list(sys.argv[1:] if argv is None else argv)
    if args == ["layout"]:
        print("\n".join(map(str, typing_files(load_paths().lab).required())))
        return 0
    if len(args) == 2 and args[0] in ("build-panel", "build-af"):
        lab = Path(args[1])
        files = typing_files(lab)
        if args[0] == "build-panel":
            build_panel_reference(
                panel_path(lab), files.panel_ref, threads=min(8, os.cpu_count() or 1)
            )
        else:
            build_roh_af(panel_path(lab), files.roh_af)
        return 0
    print(
        "usage: python -m genokeeper.typing_files layout | build-panel LAB | build-af LAB",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
