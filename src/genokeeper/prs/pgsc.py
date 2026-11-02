"""pgsc_calc v2.3.0 on our own calls (track design, slice 1d), and where the bootstrap files
what it needs (`scripts/bootstrap.sh`, step_pgs).

The lab keeps, under `<lab>/pgs/`: the chosen scores' harmonised GRCh38 scoring files
(`scoring_file`), the HGDP+1kGP reference panel (`panel_path`) and the ancestry positions
built from it (`ancestry_pvar`), pgsc_calc at the commit tag v2.3.0 points to
(`pgsc_calc_dir`), Nextflow (`nextflow_path`) with its nf-schema plugin under its own
NXF_HOME (`nxf_home`), and the Java it runs on (`java_home`). None of it holds genotypes.

The ancestry positions are the panel's GRCh38 variants pgsc_calc's ancestry step can use:
its FILTER_VARIANTS keeps the reference's unrelated samples' autosomal biallelic ACGT SNPs
with a frequency of 5% or more (`--maf 0.05`, its `maf_ref` default) before it intersects
them with the target. The panel's variant file lists every variant with no frequency, so
plink2 in pgsc_calc's own image computes them once, at bootstrap (`ancestry_positions`),
and our calls are written at these positions, not at every panel variant.

pgsc_calc's input is a VCF of our calls at every score and ancestry position
(`target_positions`, `write_target_vcf`): reference and no-call rows explicit, each
position's genotype by the rule PharmCAT's input follows (`genotype.sites`), autosomes
only. `run_pgsc_calc` runs it with `--run_ancestry` on the panel, `--scorefile` on the
local scoring files (never `--pgs_id`, which downloads), `-profile docker,arm` and our
config (`nextflow_config`), behind the core gate as `LABEL`.

Nextflow is the orchestrator and runs on this machine (the track's ruling R16): it reads
only paths and starts every process in a container, offline (`NXF_OFFLINE`, its plugin
prefetched). Its launch dir, work dir and results lie in the run's step directory under
`reports/work/`, which the deny rules cover; its NXF_HOME is the lab's, which holds no
genotype. pgsc_calc's nextflow.config (at the pinned commit) would fetch nf-core's
institutional configs from `custom_config_base`, but only while NXF_OFFLINE is unset: that
fetch carries no genotype, and its anchor is the environment our command sets, which the
config cannot change because Nextflow reads the environment before the config. Two
controls, each with its anchor:

- `--network none` and the pinned images, against a pipeline process's own network call
  carrying genotype-derived data: our config appends `--network none` to the arm profile's
  `docker.runOptions` and pins each container label to a digest written in
  `containers.IMAGES`; its anchor comes from the `docker run` command line Nextflow builds
  from that config, which the process inside cannot change because the network namespace
  is set before it starts.
- The containment check (`check_containment`), against a configuration that did not take
  (a profile or a pgsc_calc release overriding our options, a process with no container):
  after the run, every task's `.command.run` docker line must carry `--network none` and a
  pinned digest, or the results are refused. Its anchor is Nextflow's own record of what it
  ran, written before each task's container starts, which the pipeline's processes cannot
  rewrite after they ran. A process can write in its task directory while it runs, but one
  started with `--network none` has no network to hide; the check is not a defence against
  an image hostile to us, which the pinned digests are.
"""

import csv
import gzip
import io
import itertools
import math
import re
import shutil
import sys
import tarfile
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from compression import zstd
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.containers import IMAGES, PLATFORMS, docker_memory, run_container
from genokeeper.genotype.gvcf import Allele, Called
from genokeeper.genotype.sites import site_genotype
from genokeeper.phasing import Phase
from genokeeper.store.db import Level
from genokeeper.tools import gated, run_tool
from genokeeper.variants.key import ucsc_contig

PGSC_CALC_VERSION = "2.3.0"
# The commit the tag v2.3.0 points to (github.com/PGScatalog/pgsc_calc, read 2026-10-21).
PGSC_CALC_COMMIT = "72ee54f2119d2fe65ac671327342dfc412bd0a3b"
# The Nextflow release pgsc_calc v2.3.0's CI tests (.github/workflows/standard-test.yml),
# and the nf-schema plugin its nextflow.config names.
NEXTFLOW_VERSION = "25.10.0"
NF_SCHEMA_VERSION = "2.0.0"
PANEL_NAME = "pgsc_HGDP+1kGP_v1.tar.zst"
ANCESTRY_NAME = "GRCh38_HGDP+1kGP_ancestry"
# The panel's GRCh38 variant file, every variant, kept beside them (GWAS location, R22).
PANEL_PVAR_NAME = "GRCh38_HGDP+1kGP_panel.pvar.zst"
# The core gate's label for pgsc_calc's run and the panel's plink2 filter.
LABEL = "dna:prs"
# The panel's GRCh38 members plink2 reads: its genotypes and the related samples to drop.
_PANEL_MEMBERS = (".pgen", ".psam", ".pvar.zst", ".king.cutoff.out.id")
# A score our target matches below this share of its variants is unassessed (track design,
# slice 1d): applied once, to the match rate pgsc_calc reports (`sources.prs`). pgsc_calc
# itself runs with `--min_overlap 0` (its schema's minimum), so it scores every score and
# never fails a run on overlap.
MIN_OVERLAP = 0.75
# pgsc_calc's ancestry step needs a 32 GB Docker VM (track design, slice 1d, a user step).
# Below it a run warns, naming the setting; no adversary, so no refusal: an out-of-memory
# task fails in Nextflow.
MIN_DOCKER_MEMORY = 32 * 2**30
# The target's sampleset (pgsc_calc allows letters and digits only) and its one sample. The
# sampleset names pgsc_calc's output files, cached under reports/work/pgsc/, so it keeps its
# name from before the rename to genokeeper.
SAMPLESET = "genomelab"
SAMPLE = "SAMPLE"
# pgsc_calc v2.3.0's container labels (conf/modules.config) -> the image each runs.
PGSC_LABELS = {
    "pgscatalog_utils": "pgsc_utils",
    "plink2": "pgsc_plink2",
    "zstd": "pgsc_zstd",
    "report": "pgsc_report",
    "pyyaml": "pgsc_pyyaml",
    "fraposa": "pgsc_fraposa",
}
# The arm profile's `docker.runOptions` (nextflow.config at the pinned commit), and ours.
ARM_RUN_OPTIONS = "-u $(id -u):$(id -g) --platform=linux/arm64"
RUN_OPTIONS = f"{ARM_RUN_OPTIONS} --network none"
# APOE's two SNPs on GRCh38 (dbSNP): rs429358 T>C (C is ε4's, Arg112) and rs7412 C>T (T is
# ε2's, Cys158). Copies of each ALT -> the diplotype, where no phase is needed; a double
# heterozygote is ε2/ε4 in trans and ε1/ε3 in cis.
APOE_E4: Allele = ("chr19", 44908684, "T", "C")
APOE_E2: Allele = ("chr19", 44908822, "C", "T")
_APOE = {
    (0, 0): "ε3/ε3",
    (1, 0): "ε3/ε4",
    (2, 0): "ε4/ε4",
    (0, 1): "ε2/ε3",
    (0, 2): "ε2/ε2",
    (2, 1): "ε1/ε4",
    (1, 2): "ε1/ε2",
    (2, 2): "ε1/ε1",
}
_APOE_PHASED = {Phase.IN_TRANS: "ε2/ε4", Phase.IN_CIS: "ε1/ε3"}
_AUTOSOMES = frozenset(str(n) for n in range(1, 23))
_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
# A Nextflow task's work directory: two hex digits, then thirty.
_TASK = re.compile(r"[0-9a-f]{2}/[0-9a-f]{30}")


class ContainmentRefused(RuntimeError):
    """A task pgsc_calc ran cannot be shown to have run offline in a pinned image: its
    results are refused."""


class PgscFailed(RuntimeError):
    """pgsc_calc's results are missing or not in the shape this module reads."""


def pgs_dir(lab: Path) -> Path:
    return lab / "pgs"


def scoring_file(lab: Path, pgs_id: str) -> Path:
    """A score's harmonised GRCh38 scoring file, as the PGS Catalog names it."""
    return pgs_dir(lab) / "scores" / f"{pgs_id}_hmPOS_GRCh38.txt.gz"


def panel_path(lab: Path) -> Path:
    return pgs_dir(lab) / PANEL_NAME


def ancestry_pvar(lab: Path) -> Path:
    return pgs_dir(lab) / "ancestry" / f"{ANCESTRY_NAME}.pvar.zst"


def panel_pvar(lab: Path) -> Path:
    return pgs_dir(lab) / "ancestry" / PANEL_PVAR_NAME


def pgsc_calc_dir(lab: Path) -> Path:
    return pgs_dir(lab) / "pgsc_calc" / PGSC_CALC_COMMIT


def nextflow_path(lab: Path) -> Path:
    return pgs_dir(lab) / "nextflow" / f"nextflow-{NEXTFLOW_VERSION}-dist"


def nxf_home(lab: Path) -> Path:
    return pgs_dir(lab) / "nxf-home"


def java_home(lab: Path) -> Path:
    """A link to Homebrew's openjdk@21: Nextflow 25.10 runs on Java 17 to 25."""
    return pgs_dir(lab) / "java"


def extract_panel_members(panel: Path, work: Path) -> dict[str, str]:
    """The panel archive's GRCh38 genotype members (`_PANEL_MEMBERS`) extracted into `work`:
    member suffix -> file name. `ValueError` for a panel without exactly one of each."""
    found: dict[str, str] = {}
    with tarfile.open(panel, mode="r|zst") as tar:
        for member in tar:
            name = Path(member.name).name
            kind = next((k for k in _PANEL_MEMBERS if name.endswith(k)), None)
            if kind is None or not name.startswith("GRCh38_") or not member.isfile():
                continue
            if kind in found:
                raise ValueError(f"{panel}: two GRCh38 {kind} members")
            found[kind] = name
            source = tar.extractfile(member)
            if source is None:
                raise ValueError(f"{panel}: unreadable member {member.name}")
            with source, (work / name).open("wb") as target:
                shutil.copyfileobj(source, target)
    if set(found) != set(_PANEL_MEMBERS):
        raise ValueError(f"{panel}: GRCh38 members {sorted(found.values())}")
    return found


def ancestry_positions(panel: Path, out_dir: Path) -> Path:
    """The panel's ancestry positions (`ANCESTRY_NAME`.pvar.zst in `out_dir`): its GRCh38
    genotypes extracted from the archive, then filtered by plink2 offline in pgsc_calc's
    pinned image; its GRCh38 variant file, every variant, is kept beside them
    (`PANEL_PVAR_NAME`). `ValueError` for a panel without exactly one of each GRCh38 member."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=".panel-") as tmp:
        work = Path(tmp)
        found = extract_panel_members(panel, work)
        prefix = found[".pgen"].removesuffix(".pgen")
        run_container(
            "pgsc_plink2",
            [
                "plink2", "--pfile", f"/panel/{prefix}", "vzs",
                "--remove", f"/panel/{found['.king.cutoff.out.id']}",
                "--autosome", "--max-alleles", "2", "--snps-only", "just-acgt",
                "--maf", "0.05", "--make-just-pvar", "zs",
                "--out", f"/out/{ANCESTRY_NAME}",
            ],
            mounts=[(work, "/panel", True), (out_dir, "/out", False)],
            platform=PLATFORMS["pgsc_plink2"],
            label=LABEL,
        )  # fmt: skip
        (work / found[".pvar.zst"]).replace(out_dir / PANEL_PVAR_NAME)
    (out_dir / f"{ANCESTRY_NAME}.log").unlink(missing_ok=True)
    return out_dir / f"{ANCESTRY_NAME}.pvar.zst"


def _scoring_rows(path: Path) -> Iterator[tuple[str, int, str, list[str]]]:
    """A harmonised scoring file's variants: (chromosome as the file names it, GRCh38
    position, effect allele, the other alleles), rows without a GRCh38 position left out.
    The other allele is `other_allele`, or the harmonisation's inferred ones."""
    names = ("hm_chr", "hm_pos", "effect_allele", "other_allele", "hm_inferOtherAllele")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        header = next(lines, "").rstrip("\n").split("\t")
        at = [header.index(n) if n in header else None for n in names]
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            chrom, pos, effect, other, inferred = (
                fields[i] if i is not None and i < len(fields) else "" for i in at
            )
            if pos:
                others = [other] if other else [a for a in inferred.split("/") if a]
                yield chrom, int(pos), effect, others


def read_pvar(path: Path) -> Iterator[Allele]:
    """Each `(chrom, pos, ref, alt)` of a plink2 variant file (`.pvar.zst`, as the panel's
    and `ancestry_positions`' are), one per ALT, in the file's order and the GRCh38 FASTA's
    contig names (`chr1`)."""
    with zstd.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                chrom, pos, _id, ref, alts = line.rstrip("\n").split("\t")[:5]
                for alt in alts.split(","):
                    yield ucsc_contig(chrom), int(pos), ref, alt


def target_positions(
    scoring_files: Sequence[Path], ref_pvar: Path, *, fasta: Path
) -> list[Allele]:
    """Every autosomal `(chrom, pos, ref, alt)` pgsc_calc's input needs, in the GRCh38
    FASTA's contig names, sorted: each score's effect allele against each of its other
    alleles, and each of the panel's ancestry positions (`ref_pvar`), the allele the FASTA
    has at the position as REF; a score's SNV whose complement the FASTA has, at a position
    no pair matches as written, is on the reverse strand: written complemented, pgsc_calc's
    flip matching sees it. A pair the FASTA has neither of, or an allele of other letters
    than ACGT, is no position: pgsc_calc counts it unmatched."""
    wanted: dict[str, set[tuple[int, str, str, bool]]] = defaultdict(set)
    for path in scoring_files:
        for chrom, pos, effect, others in _scoring_rows(path):
            for other in others:
                wanted[chrom.removeprefix("chr")].add((pos, effect, other, True))
    for chrom, pos, ref, alt in read_pvar(ref_pvar):
        wanted[chrom.removeprefix("chr")].add((pos, ref, alt, False))
    found: list[Allele] = []
    with pysam.FastaFile(str(fasta)) as reference:
        for chrom in sorted(_AUTOSOMES & set(wanted), key=int):
            contig = ucsc_contig(chrom)
            if contig not in reference.references:
                continue
            sequence = reference.fetch(contig).upper()
            flips: list[Allele] = []
            for pos, a, b, score in wanted[chrom]:
                if not set(a + b) <= _BASES or not a or not b:
                    continue
                pairs = [(a, b), (b, a)]
                if score and len(a) == len(b) == 1:  # a palindrome flips onto itself
                    x, y = a.translate(_COMPLEMENT), b.translate(_COMPLEMENT)
                    pairs += [(x, y), (y, x)]
                for i, (ref, alt) in enumerate(pairs):
                    if sequence[pos - 1 : pos - 1 + len(ref)] == ref:
                        (found if i < 2 else flips).append((contig, pos, ref, alt))
                        break
            # An inferred other allele that matches only flipped is no strand flip.
            direct = {p for c, p, _, _ in found if c == contig}
            found += [f for f in flips if f[1] not in direct]
    return sorted(set(found), key=lambda a: (int(a[0][3:]), a[1], a[2], a[3]))


def write_target_vcf(
    calls: Iterable[tuple[Allele, Called | None, Level | None]], out: Path
) -> Path:
    """Write pgsc_calc's input to `out` (gzipped, byte for byte the same for the same
    calls): one row per position of `calls` (each allele with our call there and its grade,
    in `target_positions`' order, an allele's position's ALTs next to each other), its
    genotype by `genotype.sites.site_genotype` on an autosome.

    The calls stream through: a genome's millions of positions are never held at once.
    `out` holds the user's genotypes: the caller puts it under the run's work directory."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with (
        gzip.GzipFile(out, "wb", compresslevel=1, mtime=0) as raw,
        io.TextIOWrapper(raw, encoding="utf-8") as handle,
    ):
        handle.write(
            "##fileformat=VCFv4.2\n"
            # `genomelab`, the name before the rename: this file's sha256 is in a step key.
            f"##source=genomelab (our calls for pgsc_calc {PGSC_CALC_VERSION})\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n"
        )
        for (chrom, pos, ref), site in itertools.groupby(calls, key=lambda c: c[0][:3]):
            alleles, called, levels = zip(*site, strict=True)
            gt = site_genotype(called, levels, haploid=False)
            alts = ",".join(a[3] for a in alleles)
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alts}\t.\tPASS\t.\tGT\t{gt}\n")
    return out


def nextflow_config(digests: Mapping[str, str]) -> str:
    """Our Nextflow config: the arm profile's `docker.runOptions` with `--network none`
    appended, and each of pgsc_calc's container labels (`PGSC_LABELS`) pinned to its image
    by digest (`digests`, label -> `repo@sha256:<digest>`). `ValueError` for a label
    missing, one too many, or an image without a digest."""
    unpinned = sorted(label for label, ref in digests.items() if "@sha256:" not in ref)
    if unpinned or set(digests) != set(PGSC_LABELS):
        raise ValueError(f"labels {sorted(digests)}, unpinned {unpinned}")
    blocks = "".join(
        f"    withLabel: {label} {{\n        container = '{ref}'\n    }}\n"
        for label, ref in digests.items()
    )
    return f"docker.runOptions = '{RUN_OPTIONS}'\nprocess {{\n{blocks}}}\n"


def run_pgsc_calc(
    target_vcf: Path, scores: Sequence[Path], ref: Path, out_dir: Path, *, lab: Path
) -> Path:
    """Run pgsc_calc on `target_vcf` (`*.vcf.gz`) with the scoring files `scores` and the
    reference panel `ref`, in `out_dir` (its launch dir, with its work dir and results
    inside), behind the core gate as `LABEL`; the results' directory, once
    `check_containment` passed. Warns on stderr when Docker's VM has less than
    `MIN_DOCKER_MEMORY`. Raises `ContainmentRefused`, and what `run_tool` raises.

    Multi-allelic matches keep pgsc_calc's default (dropped), the path its authors run
    with the HGDP+1kGP panel: with `--keep_multiallelic true` the panel's scoring lost
    matched variants at its multi-allelic sites and SCORE_AGGREGATE failed
    (`Missing variants`); a position our target writes with two ALTs goes unmatched and
    shows in the match rate."""
    memory = docker_memory()
    if memory < MIN_DOCKER_MEMORY:
        print(
            f"genokeeper run: Docker gives its containers {memory / 2**30:.1f} GiB of"
            f" memory, below the {MIN_DOCKER_MEMORY // 2**30} GiB pgsc_calc's ancestry"
            " step needs: raise it in Docker Desktop (Settings, Resources)",
            file=sys.stderr,
        )
    scored = out_dir / "scores"
    scored.mkdir(parents=True)
    (out_dir / "tmp").mkdir()  # NXF_TEMP: Java's tmpdir, which Nextflow does not create
    for path in scores:
        (scored / path.name).symlink_to(path)
    sheet = out_dir / "samplesheet.csv"
    prefix = str(target_vcf).removesuffix(".vcf.gz")
    sheet.write_text(f"sampleset,path_prefix,chrom,format\n{SAMPLESET},{prefix},,vcf\n")
    config = out_dir / "genokeeper.config"
    config.write_text(
        nextflow_config({label: IMAGES[name] for label, name in PGSC_LABELS.items()})
    )
    command = [
        str(nextflow_path(lab)), "run", str(pgsc_calc_dir(lab) / "main.nf"),
        "-profile", "docker,arm", "-c", str(config), "-w", str(out_dir / "work"),
        "-ansi-log", "false",
        "--input", str(sheet), "--target_build", "GRCh38",
        "--scorefile", f"{scored}/*_hmPOS_GRCh38.txt.gz",
        "--run_ancestry", str(ref), "--min_overlap", "0",
        "--outdir", str(out_dir / "results"),
    ]  # fmt: skip
    env = {
        "NXF_HOME": str(nxf_home(lab)),
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
        "NXF_JAVA_HOME": str(java_home(lab)),
        "NXF_TEMP": str(out_dir / "tmp"),
    }
    run_tool(gated(LABEL, command), env=env, cwd=out_dir)
    check_containment(out_dir / "work")
    return out_dir / "results"


def check_containment(work_dir: Path) -> None:
    """Raise `ContainmentRefused` unless every task in Nextflow's `work_dir` has a
    `.command.run` whose every `docker run` line sets the network to `none` (and to nothing
    else) and names its image by a digest among pgsc_calc's pinned ones. A work dir with no
    task, a task without a launcher, and a launcher with no docker line are refused too:
    they verify nothing."""
    pinned = {IMAGES[name] for name in PGSC_LABELS.values()}
    tasks = sorted(
        d
        for d in work_dir.glob("*/*")
        if d.is_dir() and _TASK.fullmatch(f"{d.parent.name}/{d.name}")
    )
    if not tasks:
        raise ContainmentRefused(
            f"pgsc_calc's work dir {work_dir} has no task to check"
        )
    for task in tasks:
        launcher = task / ".command.run"
        if not launcher.is_file():
            raise ContainmentRefused(f"pgsc_calc's task {task} has no .command.run")
        lines = [
            line.split()
            for line in launcher.read_text(encoding="utf-8").splitlines()
            if "docker run" in line
        ]
        if not lines:
            raise ContainmentRefused(f"{launcher} has no docker run line")
        for tokens in lines:
            networks = [
                value
                for flag, value in itertools.pairwise(tokens)
                if flag in ("--network", "--net")
            ] + [
                t.split("=", 1)[1]
                for t in tokens
                if t.startswith(("--network=", "--net="))
            ]
            if not networks or set(networks) != {"none"}:
                raise ContainmentRefused(
                    f"{launcher} starts a container without --network none ({networks})"
                )
            images = [t for t in tokens if "@sha256:" in t] or ["an image by tag"]
            if any(image not in pinned for image in images):
                raise ContainmentRefused(
                    f"{launcher} starts {images[0]}, not one of the pinned digests"
                )


@dataclass(frozen=True)
class ScoreResult:
    """One score's result for the user: the percentile and Z against the population
    pgsc_calc finds the user closest to (None when it scored nothing: the overlap was too
    low), that population, whether it assigned it with low confidence, and the share of the
    score's variants our target matched; `z_norm2`, the `Z_norm2` column pgsc_calc's
    default normalisation writes (None when the table has no such column)."""

    pgs_id: str
    percentile: float | None
    z: float | None
    population: str
    low_confidence: bool
    match_rate: float
    z_norm2: float | None = None


def _one(results: Path, pattern: str) -> Path:
    found = sorted(results.glob(pattern))
    if len(found) != 1:
        raise PgscFailed(f"{results}: {len(found)} files {pattern}, expected one")
    return found[0]


def _pgs_id(accession: str) -> str:
    """The PGS ID an accession starts with (pgsc_calc may suffix the file's name)."""
    match = re.match(r"PGS\d{6}", accession)
    if match is None:
        raise PgscFailed(f"pgsc_calc reported a score {accession!r}")
    return match.group()


def _number(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    return None if math.isnan(value) else value


def _ours(path: Path) -> list[dict[str, str]]:
    """The target sampleset's rows of a gzipped tab-separated pgsc_calc table."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if r["sampleset"] == SAMPLESET
        ]


# pandas' spellings of the low-confidence flag (pgscatalog-utils 2.0.4,
# lib/_ancestry/tools.py): a bool, or NaN, written empty, when no p-value threshold was
# given, as pgsc_calc v2.3.0's ancestry_analysis.nf gives none: its assignment stands
# unflagged.
_LOW_CONFIDENCE = {"True": True, "False": False, "": False}


def _low_confidence(value: str) -> bool:
    if value not in _LOW_CONFIDENCE:
        raise PgscFailed(f"pgsc_calc's low-confidence flag {value!r}")
    return _LOW_CONFIDENCE[value]


def parse_results(results: Path) -> dict[str, ScoreResult]:
    """Each score's result from pgsc_calc's `results` directory: the match rate from the
    sampleset's match summary (matched variants over all the score's), the percentile and Z
    from its scores, the closest population from its population similarity. `PgscFailed`
    for a missing table or a sample not there exactly once."""
    counts: dict[str, list[int]] = {}
    with _one(results, f"{SAMPLESET}/match/*_summary.csv").open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            count = counts.setdefault(_pgs_id(row["accession"]), [0, 0])
            count[1] += int(row["count"])
            count[0] += int(row["count"]) if row["match_status"] == "matched" else 0
    rows = _ours(_one(results, f"{SAMPLESET}/score/*_pgs.txt.gz"))
    scored = {
        _pgs_id(row["PGS"]): (
            _number(row["percentile_MostSimilarPop"]),
            _number(row["Z_MostSimilarPop"]),
        )
        for row in rows
    }
    z_norm2 = {
        _pgs_id(row["PGS"]): _number(row["Z_norm2"]) for row in rows if "Z_norm2" in row
    }
    similar = _ours(_one(results, f"{SAMPLESET}/score/*_popsimilarity.txt.gz"))
    if len(similar) != 1:
        raise PgscFailed(
            f"{len(similar)} rows of {SAMPLESET} in its population similarity"
        )
    (owner,) = similar
    return {
        pgs_id: ScoreResult(
            pgs_id,
            *scored.get(pgs_id, (None, None)),
            population=owner["MostSimilarPop"],
            low_confidence=_low_confidence(owner["MostSimilarPop_LowConfidence"]),
            match_rate=matched / total,
            z_norm2=z_norm2.get(pgs_id),
        )
        for pgs_id, (matched, total) in counts.items()
    }


@dataclass(frozen=True)
class PopSimilarity:
    """pgsc_calc's random-forest similarity of the user to the panel's superpopulations
    (`*_popsimilarity.txt.gz`): the closest one, whether it was assigned with low
    confidence, and each superpopulation's probability (`RF_P_<POP>`). A similarity, not
    ancestry proportions."""

    population: str
    low_confidence: bool
    probabilities: dict[str, float]


def parse_similarity(results: Path) -> PopSimilarity:
    """The sampleset's row of pgsc_calc's population similarity in `results`; `PgscFailed`
    for a missing table, a sample not there exactly once, or a probability not a number."""
    similar = _ours(_one(results, f"{SAMPLESET}/score/*_popsimilarity.txt.gz"))
    if len(similar) != 1:
        raise PgscFailed(
            f"{len(similar)} rows of {SAMPLESET} in its population similarity"
        )
    (owner,) = similar
    probabilities: dict[str, float] = {}
    for column, value in owner.items():
        if column.startswith("RF_P_"):
            number = _number(value)
            if number is None:
                raise PgscFailed(f"pgsc_calc's {column} {value!r}")
            probabilities[column.removeprefix("RF_P_")] = number
    return PopSimilarity(
        population=owner["MostSimilarPop"],
        low_confidence=_low_confidence(owner["MostSimilarPop_LowConfidence"]),
        probabilities=probabilities,
    )


def apoe_genotype(
    dosages: Mapping[Allele, str], phase: Phase | None = None
) -> str | None:
    """APOE's diplotype from the dosages (`0`, `1`, `2`, `unknown`) at `APOE_E4` and
    `APOE_E2`, and for a double heterozygote the read-backed `phase` of the two; None when
    a dosage is unknown, or a double heterozygote is not phased."""
    e4, e2 = dosages[APOE_E4], dosages[APOE_E2]
    if "unknown" in (e4, e2):
        return None
    if (e4, e2) == ("1", "1"):
        return _APOE_PHASED.get(phase) if phase is not None else None
    return _APOE[(int(e4), int(e2))]


def main(argv: Sequence[str] | None = None) -> int:
    """`python -m genokeeper.prs.pgsc ancestry-positions PANEL OUT_DIR` (the bootstrap's
    step_pgs)."""
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3 or args[0] != "ancestry-positions":
        print(
            "usage: python -m genokeeper.prs.pgsc ancestry-positions PANEL OUT_DIR",
            file=sys.stderr,
        )
        return 2
    ancestry_positions(Path(args[1]), Path(args[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
