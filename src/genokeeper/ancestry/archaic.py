"""The archaic share by hmmix 0.9.2 (slice 3f, "Admixture"; ruling R5).

`run_hmmix` runs, in one hmmix container (offline, gated `dna:ancestry`), the steps Task 1
recorded: the user's PASS calls split one VCF per autosome (hmmix pairs each with its
ancestral FASTA by the varying part of the names), `create_ingroup` against the 1000G+HGDP
outgroup on the strict callability mask (refused when the observations come back nearly
empty, as they do silently without vcftools), `train`, and `decode -admixpop` with the four
high-coverage archaic genomes, which writes `decoded.diploid.txt`.

A decoded `Archaic` segment is Neanderthal when the most derived SNPs it shares with any of
the three Neanderthals (Altai, Vindija, Chagyrskaya) outnumber those it shares with
Denisova, Denisovan when Denisova's are more, and unclassified on a tie (`Segment.origin`;
hmmix's README reads its example segment this way, and the published anchor kept fragments
that share more with a Neanderthal than with Denisova). `archaic_share` takes the
Neanderthal segments' raw span over `ARCHAIC_GENOME_BP` (R15: the anchor and the statement
are the Neanderthal share; the Denisovan span is kept beside it), the published
anchor's fixed denominator (Coll Macià et al. 2026, Mol Biol Evol 43:msag134: 2 ×
2,550,757,832 bp, autosomes, both haplotypes), not the strict mask's sum: the unphased
segments' span stands for both haplotypes', as archaic segments rarely overlap between them.
"""

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from genokeeper.containers import PLATFORMS, run_container
from genokeeper.typing_files import LABEL

ARCHAIC_GENOME_BP = 2 * 2_550_757_832
DECODED = "decoded.diploid.txt"
NEANDERTHALS = ("AltaiNeandertal", "Vindija33.19", "Chagyrskaya-Phalanx")
DENISOVA = "Denisova"
# Fewer observation rows than this means create_ingroup failed silently (Task 5's HG002
# truth gave 133,532).
MIN_OBS_ROWS = 1000

_SCRIPT = """\
set -euo pipefail
vcf="$1"; work="$2"; min_obs="$3"
cd "$work"
mkdir -p anc vcf
chroms=$(seq 1 22)
tar -xzf /h/hg38_ancestral.tar.gz -C anc \
  $(for c in $chroms; do printf 'hg38_ancestral/homo_sapiens_ancestor_%s.fa ' "$c"; done)
ind=$(bcftools query -l "$vcf" | head -n 1)
for c in $chroms; do
  bcftools view -f PASS -r "chr$c" "$vcf" -Oz -o "vcf/owner.chr$c.vcf.gz"
  bcftools index -f "vcf/owner.chr$c.vcf.gz"
done
join() { local IFS=,; echo "$*"; }
vcfs=$(join $(for c in $chroms; do echo "vcf/owner.chr$c.vcf.gz"; done))
ancs=$(join $(for c in $chroms; do echo "anc/hg38_ancestral/homo_sapiens_ancestor_$c.fa"; done))
arch=$(join $(for c in $chroms; do echo "/h/archaic/individuals_highcov.$c.bcf"; done))
hmmix create_ingroup -ind="$ind" -vcf="$vcfs" -weights=/h/hg38_strick_callability_mask.bed \
  -outgroup=/h/hg38_Outgroup_1000g_HGDP.txt -ancestral="$ancs" -out=obs
obs="obs.$ind.txt"
[ "$(wc -l < "$obs")" -gt "$min_obs" ] || { echo "hmmix: $obs nearly empty" >&2; exit 3; }
hmmix train -obs="$obs" -weights=/h/hg38_strick_callability_mask.bed -mutrates=/h/hg38_mutationrate.bed \
  -out=trained.json
hmmix decode -obs="$obs" -weights=/h/hg38_strick_callability_mask.bed -mutrates=/h/hg38_mutationrate.bed \
  -param=trained.json -admixpop="$arch" -out=decoded
rm -rf anc vcf
"""


@dataclass(frozen=True)
class Segment:
    """One decoded segment: its span (0-based, end exclusive), state (`Human` or
    `Archaic`), and the most derived SNPs it shares with one Neanderthal and with
    Denisova."""

    chrom: str
    start: int
    end: int
    state: str
    neanderthal: int
    denisova: int

    @property
    def origin(self) -> Literal["neanderthal", "denisovan"] | None:
        """Which archaic group an `Archaic` segment shares more derived SNPs with; None for
        a tie or a `Human` segment."""
        if self.state != "Archaic" or self.neanderthal == self.denisova:
            return None
        return "neanderthal" if self.neanderthal > self.denisova else "denisovan"


@dataclass(frozen=True)
class ArchaicShare:
    """The Neanderthal and Denisovan segments' summed spans, and the Neanderthal span's
    share of the denominator."""

    neanderthal_bp: int
    denisovan_bp: int
    share: float


def run_hmmix(vcf: Path, files: Path, out: Path) -> Path:
    """hmmix on the user's indexed `vcf` with its helper files `files`
    (`typing_files.TypingFiles.hmmix`), its work and `DECODED` in `out`; the decoded file.
    `RuntimeError` when it exits 0 without it, and what `run_container` raises."""
    out.mkdir(parents=True, exist_ok=True)
    run_container(
        "hmmix",
        ["bash", "-c", _SCRIPT, "hmmix", f"/in/{vcf.name}", "/w", str(MIN_OBS_ROWS)],
        mounts=[(vcf.parent, "/in", True), (files, "/h", True), (out, "/w", False)],
        platform=PLATFORMS["hmmix"],
        label=LABEL,
    )
    decoded = out / DECODED
    if not decoded.is_file():
        raise RuntimeError(f"hmmix exited 0 without {DECODED}")
    return decoded


def read_segments(path: Path) -> list[Segment]:
    """Each segment of hmmix's `decode -admixpop` output."""
    with path.open(newline="", encoding="utf-8") as f:
        return [
            Segment(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                state=row["state"],
                neanderthal=max(int(row[name]) for name in NEANDERTHALS),
                denisova=int(row[DENISOVA]),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]


def archaic_share(
    segments: Sequence[Segment], genome_bp: int = ARCHAIC_GENOME_BP
) -> ArchaicShare:
    """The Neanderthal and Denisovan segments' raw spans, and the Neanderthal span over
    `genome_bp`; unclassified segments count for neither."""
    spans = {"neanderthal": 0, "denisovan": 0}
    for s in segments:
        if s.origin is not None:
            spans[s.origin] += s.end - s.start
    return ArchaicShare(
        neanderthal_bp=spans["neanderthal"],
        denisovan_bp=spans["denisovan"],
        share=spans["neanderthal"] / genome_bp,
    )
