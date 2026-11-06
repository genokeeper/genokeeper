"""chrM with mutserve: each variant's heteroplasmy level, and a chrM claim's call from them.

DeepVariant's chrM calls are not used (track design, 2a "Mitochondria"): a diploid caller
reads a 20% heteroplasmy as noise or a het. mutserve 2.0.3 (the mtDNA-Server 2 caller) runs
in its pinned container (`containers.IMAGES["mutserve"]`, offline, behind the core gate as
`dna:mito`) on the chrM reads alone, which samtools extracts from the CRAM natively, against
the FASTA's chrM (rCRS in the GRCh38 no-alt set). It writes one row per variant from level
0.01 (`MUTSERVE_LEVEL`), with its level and coverage, and with `--write-raw` the coverage at
every position. Its own filters stand: mapping and base quality at least 20 (its defaults,
passed explicitly). Indels are not called (mutserve's `--deletions` and `--insertions` are
beta), so a chrM indel claim stays unknown.

A chrM claim's call (`mito_called`): mutserve's variant with the claim's allele, on at least
[reliability] min_dp reads (else unknown), is one copy of the mitochondrial genome, carried at
its level (`Called.vaf`), a heteroplasmy below [mito] min_homoplasmy; another base called
there is dosage 0 with the other-allele warning; no variant where at least [reliability]
min_dp reads cover the position is the reference, as a reference block (no candidate at a
covered position); fewer reads is unknown. Grading
(`grade_mito`): mutserve gives no error estimate and stage 2b calibrates no chrM, so a call
is never high: medium for a homoplasmic call or the reference, low for a heteroplasmy or an
unknown; a known call says that mtDNA comes from the mother.

Slice 3f ("chrM floor", "haplocheck"): given the run's heteroplasmy floor, the claim's allele
below it is unknown with its level kept, and one at or above it, still a heteroplasmy, says
it depends on the tissue. A heteroplasmy is unknown whatever the floor when haplocheck flags
contamination (a contaminant's variant would pass as one) or fails, or the calibration
withholds heteroplasmies. Each such unknown's first warning starts `genotype unknown`, as `query` reads it, and the next keeps the level.
Without these keywords a call is graded exactly as before 3f.
"""

import csv
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from genokeeper.containers import PLATFORM, run_container
from genokeeper.genotype.gvcf import NO_COVERAGE_REASON, OTHER_ALLELE_WARNING, Called
from genokeeper.genotype.reliability import OUTSIDE_CALIBRATED, Reliability, StrataIndex
from genokeeper.store.db import Level, RowType
from genokeeper.tools import run_tool

CHRM = "chrM"
# mutserve's lowest reported level (its --level default), its mapping and base quality.
MUTSERVE_LEVEL = 0.01
MUTSERVE_MAPQ = 20
MUTSERVE_BASEQ = 20
CALLS_FILE = "mutserve.txt"
RAW_FILE = "mutserve_raw.txt"
# mutserve's VCF, for haplocheck; with it, mutserve names its text table from the part of
# this name before the first dot, which is CALLS_FILE's.
VCF_FILE = "mutserve.vcf.gz"
HETEROPLASMY_WARNING = "heteroplasmy"
INDEL_REASON = "genotype unknown: mutserve calls no indels on chrM"
MUTSERVE_WARNING = "mutserve level; no error estimate, not calibrated"
TISSUE_NOTE = "tissue-dependent; confirm in blood or urine"
MATERNAL_NOTE = "mtDNA is inherited from the mother; a father passes on none"


@dataclass(frozen=True)
class MitoCall:
    """One mutserve variant: position on chrM, REF, the variant base, its level (the share
    of reads carrying it) and the reads covering the position."""

    pos: int
    ref: str
    alt: str
    level: float
    dp: int


def read_calls(out_dir: Path) -> list[MitoCall]:
    """mutserve's variants from `out_dir/mutserve.txt`."""
    with (out_dir / CALLS_FILE).open(newline="") as f:
        return [
            MitoCall(
                pos=int(row["Pos"]),
                ref=row["Ref"].upper(),
                alt=row["Variant"].upper(),
                level=float(row["VariantLevel"]),
                dp=int(row["Coverage"]),
            )
            for row in csv.DictReader(f, delimiter="\t")
        ]


def chrm_coverage(out_dir: Path) -> dict[int, int]:
    """Reads covering each chrM position, from mutserve's raw table (`--write-raw`)."""
    with (out_dir / RAW_FILE).open(newline="") as f:
        return {
            int(row["POS"]): int(row["COV-TOTAL"])
            for row in csv.DictReader(f, delimiter="\t")
        }


def call_chrm(
    cram: Path, fasta: Path, out_dir: Path, *, threads: int = 2
) -> list[MitoCall]:
    """mutserve on the chrM reads of `cram` (indexed, on `fasta`), its two tables and its
    VCF (`VCF_FILE`, from the same run) written into `out_dir` (created if missing);
    returns its variants. The extracted reads are removed afterwards."""
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".mito-", dir=out_dir))
    try:
        bam, chrm = work / "chrM.bam", work / "chrM.fa"
        run_tool(
            [
                "samtools",
                "view",
                "-b",
                "-T",
                str(fasta),
                "-o",
                str(bam),
                str(cram),
                CHRM,
            ]
        )
        run_tool(["samtools", "index", str(bam)])
        chrm.write_text(run_tool(["samtools", "faidx", str(fasta), CHRM]))
        run_tool(["samtools", "faidx", str(chrm)])
        run_container(
            "mutserve",
            [
                "mutserve",
                "call",
                "--reference",
                "/work/chrM.fa",
                "--output",
                f"/work/{VCF_FILE}",
                "--level",
                str(MUTSERVE_LEVEL),
                "--mapQ",
                str(MUTSERVE_MAPQ),
                "--baseQ",
                str(MUTSERVE_BASEQ),
                "--threads",
                str(threads),
                "--write-raw",
                "--no-ansi",
                "/work/chrM.bam",
            ],
            mounts=[(work, "/work", False)],
            platform=PLATFORM,
            label="dna:mito",
        )
        for name in (CALLS_FILE, RAW_FILE, VCF_FILE):
            if not (work / name).is_file():
                raise RuntimeError(f"mutserve exited 0 without {name}")
            (work / name).replace(out_dir / name)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return read_calls(out_dir)


def mito_called(
    allele: tuple[str, int, str, str],
    calls: Sequence[MitoCall],
    coverage: Mapping[int, int],
    *,
    min_dp: int,
    min_homoplasmy: float,
    floor: float | None = None,
    heteroplasmy_withheld: str | None = None,
) -> Called:
    """The call at a chrM claim's `(chrom, pos, ref, alt)` from mutserve's variants and
    coverage (see the module docstring). A heteroplasmy (a level below `min_homoplasmy`) is
    unknown with `heteroplasmy_withheld`'s reason when given (the calibration's, or
    haplocheck's flag or failure), else when below `floor`; one at or above `floor` is
    warned `TISSUE_NOTE`. Each keyword left out leaves the call as before slice 3f."""
    _, pos, ref, alt = allele
    if len(ref) != 1 or len(alt) != 1:
        return _unknown(INDEL_REASON)
    at = [c for c in calls if c.pos == pos]
    own = next((c for c in at if c.alt == alt.upper()), None)
    if own is not None:
        if own.dp < min_dp:  # as in the haplogroup's VCF: too few reads are no evidence
            return _unknown(NO_COVERAGE_REASON, dp=own.dp)
        level = own.level
        if level < min_homoplasmy:
            why = heteroplasmy_withheld
            if why is None and floor is not None and level < floor:
                why = f"heteroplasmy at level {level:.3g}, below the floor {floor:g}"
            if why is not None:
                return Called(
                    gt=None, dosage="unknown", row_type=RowType.PASS, gq=None,
                    dp=own.dp, ad=None, vaf=level,
                    warnings=(
                        f"genotype unknown: {why}",
                        f"{HETEROPLASMY_WARNING} at level {level:.3g}",
                    ),
                )  # fmt: skip
        warnings = (
            ()
            if level >= min_homoplasmy
            else (f"{HETEROPLASMY_WARNING} at level {level:.3g}",)
            + ((TISSUE_NOTE,) if floor is not None else ())
        )
        return Called(
            gt="1",
            dosage="1",
            row_type=RowType.PASS,
            gq=None,
            dp=own.dp,
            ad=None,
            vaf=own.level,
            warnings=warnings,
        )
    if at:
        return Called(
            gt="0",
            dosage="0",
            row_type=RowType.PASS,
            gq=None,
            dp=at[0].dp,
            ad=None,
            vaf=None,
            warnings=(OTHER_ALLELE_WARNING,),
        )
    depth = coverage.get(pos, 0)
    if depth < min_dp:
        return _unknown(NO_COVERAGE_REASON, dp=depth)
    return Called(
        gt="0",
        dosage="0",
        row_type=RowType.REF_BLOCK,
        gq=None,
        dp=depth,
        ad=None,
        vaf=None,
        warnings=(),
    )


def _unknown(reason: str, *, dp: int | None = None) -> Called:
    return Called(
        gt=None, dosage="unknown", row_type=RowType.NONE, gq=None, dp=dp, ad=None,
        vaf=None, warnings=(reason,),
    )  # fmt: skip


def grade_mito(
    called: Called,
    chrom: str,
    pos: int,
    *,
    strata: StrataIndex,
    min_homoplasmy: float,
) -> Reliability:
    """A chrM call's reliability: never high (no error estimate, no calibration: a known
    genotype is warned `OUTSIDE_CALIBRATED`, as HG003's truth covers chr1-22 only, and
    `MATERNAL_NOTE`); medium for a homoplasmic call, another base called, or the reference;
    low for a heteroplasmy or an unknown genotype."""
    flags = strata.flags_at(chrom, pos)
    warnings = list(called.warnings)
    if called.dosage == "unknown":
        return Reliability(
            p_wrong=None, level=Level.LOW, flags=flags, warnings=tuple(warnings)
        )
    warnings[:0] = [MUTSERVE_WARNING, OUTSIDE_CALIBRATED, MATERNAL_NOTE]
    heteroplasmy = called.vaf is not None and called.vaf < min_homoplasmy
    level = Level.LOW if heteroplasmy else Level.MEDIUM
    return Reliability(p_wrong=None, level=level, flags=flags, warnings=tuple(warnings))
