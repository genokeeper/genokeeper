"""Where things live on this machine, and the grading thresholds from config/grading.toml."""

import tomllib
from dataclasses import asdict, dataclass, field
from importlib import resources
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any

from genokeeper.env import getenv
from genokeeper.repeats.loci import RepeatLocus, Rule
from genokeeper.store.db import Topic, Validity

# src/genokeeper/config.py -> the repository root.
_REPO = Path(__file__).resolve().parents[2]
GRADING_FILE = _REPO / "config" / "grading.toml"


@dataclass(frozen=True)
class Paths:
    """The repository, the user's raw data and reports inside it, and the tool lab outside it."""

    repo: Path
    raw: Path
    reports: Path
    lab: Path


DEMO_MARKER = "demo.json"


def demo_lab(lab: Path | None) -> bool:
    """Whether `lab` is a demo folder's own `lab/` (`load_paths`), which has no update check."""
    return (
        lab is not None and lab.name == "lab" and (lab.parent / DEMO_MARKER).is_file()
    )


def skill_resource() -> Traversable:
    """The packaged genome-ask SKILL.md, which the wheel carries."""
    return (
        resources.files("genokeeper")
        / "resources"
        / "skills"
        / "genome-ask"
        / "SKILL.md"
    )


def load_paths(root: Path | None = None) -> Paths:
    """Paths under `root`; by default this repository when it is a checkout, else the
    working folder (the installed tool's). `$GENOKEEPER_HOME` (`$GENOMELAB_HOME` before
    the rename) moves the lab, except in a demo folder (it holds `DEMO_MARKER`), whose
    lab is its own `lab/`."""
    if root is not None:
        repo = root
    elif (_REPO / "pyproject.toml").is_file() and (
        _REPO / "src" / "genokeeper"
    ).is_dir():
        repo = _REPO
    else:
        repo = Path.cwd()
    home = getenv("HOME")
    if (repo / DEMO_MARKER).is_file():
        lab = repo / "lab"
    else:
        lab = Path(home) if home else Path.home() / "genome-lab"
    return Paths(
        repo=repo, raw=repo / "data" / "raw", reports=repo / "reports", lab=lab
    )


MASKED_REFERENCE = (
    "GRCh38 no-alt analysis set, GRC false duplications masked"
    " (GIAB maskedGRC_exclusions_v2)"
)


def grch38_fasta(lab: Path) -> Path:
    """The lab's GRCh38 FASTA, `MASKED_REFERENCE` (`chr` names), with its `.fai`, bwa index
    and rtg SDF beside it: `<lab>/reference/GRCh38_no_alt_maskedGRC.fa`. Every alignment,
    call and calibration uses it."""
    return lab / "reference" / "GRCh38_no_alt_maskedGRC.fa"


def unmasked_fasta(lab: Path) -> Path:
    """NCBI's GRCh38 no-alt analysis set as published, with its `.fai`:
    `<lab>/reference/GRCh38_no_alt.fa`. The bootstrap (`scripts/bootstrap.sh`,
    step_grch38_fasta) fetches it here; CRAMs written before the masked reference decode
    only against it, and the gnomAD FAF builder normalises with it."""
    return lab / "reference" / "GRCh38_no_alt.fa"


@dataclass(frozen=True, kw_only=True)
class QcThresholds:
    """The stage 2 QC gate ([qc]): a sample outside any of these is refused before calling.
    The defaults are config/grading.toml's values, whose comments give their sources."""

    min_mean_depth: float = 15.0
    depth_mapq: int = 20
    min_fraction_20x: float = 0.80
    min_mapped_fraction: float = 0.98
    max_duplicate_fraction: float = 0.10
    min_median_insert: float = 250.0
    max_freemix: float = 0.03
    min_titv: float = 1.9
    max_titv: float = 2.3
    max_het_hom: float = 3.0


@dataclass(frozen=True, kw_only=True)
class SexThresholds:
    """The sex-chromosome depths ([sex]): reads at MAPQ >= `mapq` count. The X copy count
    (`DosageThresholds`) decides the sex used for calling."""

    mapq: int = 30


# The lowest chrM heteroplasmy floor the config may set (slice 3f, "chrM floor").
MIN_HETEROPLASMY_FLOOR = 0.02


@dataclass(frozen=True, kw_only=True)
class MitoThresholds:
    """chrM from mutserve ([mito]): a level at or above `min_homoplasmy` is homoplasmic,
    the whole of the user's mitochondria; below it, a heteroplasmy. A chrM claim's allele
    below `heteroplasmy_floor` is unknown (slice 3f), which the calibration may raise but
    never lower."""

    min_homoplasmy: float = 0.9
    heteroplasmy_floor: float = 0.10


@dataclass(frozen=True, kw_only=True)
class Hhv6Thresholds:
    """Inherited integrated HHV-6 ([hhv6], slice 3f): the viral depth over a genome's
    unique region divided by the haploid human depth is about 1 for one integrated copy per
    cell; a ratio from `min_ratio` to `max_ratio` is a finding. The calibration's highest
    background ratio times `background_factor` is the ceiling below which a ratio is a
    non-finding."""

    min_ratio: float = 0.5
    max_ratio: float = 1.5
    background_factor: float = 3.0


@dataclass(frozen=True, kw_only=True)
class RohThresholds:
    """Runs of homozygosity ([roh], slice 3f): runs of `min_run` bp or more are summed and
    counted; a single run of `isodisomy_min` bp or more may suggest isodisomy."""

    min_run: int = 5_000_000
    isodisomy_min: int = 10_000_000


@dataclass(frozen=True, kw_only=True)
class PredictedThresholds:
    """Slice 1e's thresholds for variants ClinVar has not classified ([predicted]): a
    missense at or above `revel` (AlphaMissense's `alphamissense` where REVEL has none), a
    SpliceAI delta at or above `spliceai`; a loss of function's pext at least `min_pext` of
    its gene's largest, and no indel in a homopolymer of `homopolymer` bases or more."""

    revel: float = 0.773
    alphamissense: float = 0.906
    spliceai: float = 0.2
    min_pext: float = 0.5
    homopolymer: int = 6


@dataclass(frozen=True, kw_only=True)
class RepeatThresholds:
    """The repeat expansions ([repeats], [repeat_loci]): the reads an allele needs to carry
    a finding, and each locus's STRchive thresholds and grading rule (slices 3a and 3e).
    The built-in `loci` is empty: the shipped file holds them."""

    min_allele_reads: int = 5
    loci: dict[str, RepeatLocus] = field(default_factory=dict)


# One or two shared reads never exclude an expansion: the floor of [exclusion] min_spanning.
# Its anchor comes from this constant, which an edited grading.toml cannot move because the
# floor is code, not read from the file.
MIN_SPANNING_FLOOR = 3


@dataclass(frozen=True, kw_only=True)
class ExclusionThresholds:
    """Slice 3e's exclusion rule ([exclusion]): the spanning reads every allele needs before
    a non-standard repeat locus is called normal."""

    min_spanning: int = 3


@dataclass(frozen=True, kw_only=True)
class SmnThresholds:
    """SMN1 from two callers ([smn]): Parascopy's call counts only at this quality or above."""

    min_parascopy_qual: float = 20.0


@dataclass(frozen=True, kw_only=True)
class YleafThresholds:
    """The Y haplogroup from Yleaf ([yleaf]): below this QC score no haplogroup is stated."""

    min_qc_score: float = 0.95


@dataclass(frozen=True, kw_only=True)
class DeletionThresholds:
    """Slice 3b's deletion screen ([deletions]): the size floor and the long band's start,
    the exon margin, what is common in gnomAD-SV (AF above `max_gnomad_af` among at least
    `min_gnomad_an` alleles, at `min_reciprocal_overlap`), DELLY's depth ratio below which
    depth supports a deletion, the length from which the depth decides, and the chance of
    being false below which it is "likely"."""

    min_length: int = 50
    long_from: int = 1000
    exon_margin: int = 20
    max_gnomad_af: float = 0.01
    min_gnomad_an: int = 2000
    min_reciprocal_overlap: float = 0.5
    max_depth_ratio: float = 0.7
    depth_decides_from: int = 1_000_000
    likely_below: float = 0.1


@dataclass(frozen=True, kw_only=True)
class IdentityThresholds:
    """Slice 3c's sample identity ([identity]): the run's genotypes against the delivered
    23andMe-format file at the SNPs both hold. Below `min_concordance` over at least
    `min_compared` SNPs the run refuses; with fewer compared, identity is not checked."""

    min_concordance: float = 0.98
    min_compared: int = 1000


@dataclass(frozen=True, kw_only=True)
class CoverageThresholds:
    """Slice 3c's per-gene callability ([coverage]): the share of a gene's MANE CDS
    (+/- `cds_margin` bp) covered at MAPQ >= `mapq` by at least `min_depth_per_copy` reads
    per expected copy; a gene below `min_callable` is not reliably assessed."""

    mapq: int = 20
    min_depth_per_copy: int = 10
    cds_margin: int = 2
    min_callable: float = 0.95


@dataclass(frozen=True, kw_only=True)
class DosageThresholds:
    """Slice 3c's sex-chromosome dosage ([dosage]): a depth ratio is n copies within
    `tolerance` of n/2; the chrX allele balance counts sites whose ALT fraction is at least
    `het_low`, balanced at most `het_high`, over at least `min_x_sites`; `x_het_min` of
    them balanced in one X suggests an extra X; a Y ratio at most `y_absent_max` is no Y.
    `calibrated_on` empty: the bounds are not yet set from the calibration report (the
    statement's `source_version`; the `foundation` table's staleness decides whether the
    statement grades, and FMR1 and AFF2 read the dosage only with a value and a current
    table)."""

    # Half the smallest mosaic share the screen reports (30%): a design value, not
    # calibrated.
    tolerance: float = 0.15
    het_low: float = 0.2
    het_high: float = 0.8
    x_het_min: float = 0.33
    min_x_sites: int = 1000
    y_absent_max: float = 0.06
    calibrated_on: str = "2026-11-05 HG002 pcrfree/pcr, HG003, HG004"


# The genes clonal haematopoiesis mutates most (the closing design, "Allele-fraction guard").
CLONAL_GENES = (
    "TP53", "CHEK2", "DNMT3A", "TET2", "ASXL1", "PPM1D", "JAK2", "SF3B1", "SRSF2",
)  # fmt: skip


@dataclass(frozen=True, kw_only=True)
class AfGuardThresholds:
    """Slice 3c's allele-fraction guard ([af_guard]): a heterozygous call whose VAF is below
    `min_het_vaf` is warned, and in one of `clonal_genes` its level is capped at low.
    `calibrated_on` empty: the bound is the design's, not yet set from the calibration
    report."""

    min_het_vaf: float = 0.25
    clonal_genes: tuple[str, ...] = CLONAL_GENES
    calibrated_on: str = "2026-11-05 HG002 pcrfree"


@dataclass(frozen=True, kw_only=True)
class InsertionThresholds(DeletionThresholds):
    """Slice 3d's insertion screen ([insertions]): `[deletions]`'s values, and the distance
    within which a gnomAD-SV insertion is the same one (Truvari's `--refdist` default)."""

    max_distance: int = 500


@dataclass(frozen=True, kw_only=True)
class CnvThresholds:
    """Slice 3d's depth copy-number screen ([cnv]): the size floor, where the `5kb+` band
    starts, where the `large` band starts, and CNVpytor's bin size. The exon margin, gnomAD-SV
    rules, reciprocal overlap and the "likely" threshold are `[deletions]`'s."""

    min_length: int = 1000
    band_from: int = 5000
    large_from: int = 500_000
    bin_size: int = 500


@dataclass(frozen=True, kw_only=True)
class ParalogThresholds:
    """Slice 3d's paralog loci ([paralogs]): Parascopy's call counts at this quality or
    above; STRC's one copy needs CATSPER2's unique depth over the run's mean within
    `unique_loss_ratio`; Gauchian's carrier call needs `gba_min_reads` reads supporting its
    variant (no default: the shipped file holds it with its source)."""

    min_qual: float = 20.0
    unique_loss_ratio: tuple[float, float] = (0.3, 0.7)
    gba_min_reads: int


@dataclass(frozen=True, kw_only=True)
class LpaThresholds:
    """Slice 3d's Lp(a) KIV-2 count ([lpa]): the calibration withholds LPA when KILDA is
    off HG002's assembly-derived count by more than `max_kiv2_error` of that count (a
    fraction, not copies; no default: the shipped file holds it with its source)."""

    max_kiv2_error: float


@dataclass(frozen=True, kw_only=True)
class Grading:
    """Display-level cut-offs, GRCh38 PAR intervals, frequency thresholds, GRCh38 stratum
    BED names, and the stage 2 QC, sex and chrM thresholds. `ba1`, `bs1`, `paralogs` and
    `lpa` have no defaults: a Grading built in code names them."""

    p_high: float = 1e-3
    p_low: float = 1e-2
    min_dp: int = 10
    max_dp_factor: float = 2.0
    par: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    ba1: float
    bs1: float
    bs1_by_gene: dict[str, float] = field(default_factory=dict)
    strata: dict[str, str] = field(default_factory=dict)
    qc: QcThresholds = field(default_factory=QcThresholds)
    sex: SexThresholds = field(default_factory=SexThresholds)
    mito: MitoThresholds = field(default_factory=MitoThresholds)
    hhv6: Hhv6Thresholds = field(default_factory=Hhv6Thresholds)
    roh: RohThresholds = field(default_factory=RohThresholds)
    predicted: PredictedThresholds = field(default_factory=PredictedThresholds)
    repeats: RepeatThresholds = field(default_factory=RepeatThresholds)
    exclusion: ExclusionThresholds = field(default_factory=ExclusionThresholds)
    smn: SmnThresholds = field(default_factory=SmnThresholds)
    yleaf: YleafThresholds = field(default_factory=YleafThresholds)
    deletions: DeletionThresholds = field(default_factory=DeletionThresholds)
    identity: IdentityThresholds = field(default_factory=IdentityThresholds)
    coverage: CoverageThresholds = field(default_factory=CoverageThresholds)
    dosage: DosageThresholds = field(default_factory=DosageThresholds)
    af_guard: AfGuardThresholds = field(default_factory=AfGuardThresholds)
    insertions: InsertionThresholds = field(default_factory=InsertionThresholds)
    cnv: CnvThresholds = field(default_factory=CnvThresholds)
    paralogs: ParalogThresholds
    lpa: LpaThresholds

    def bs1_for(self, gene: str, *, autosomal_recessive: bool = False) -> float | None:
        """The gene's BS1 from [bs1_by_gene]; else None for an autosomal recessive gene
        (carriers are healthy, so a pathogenic allele may be common: only BA1 applies);
        else the fallback `bs1` (X-linked and unknown inheritance included)."""
        if gene in self.bs1_by_gene:
            return self.bs1_by_gene[gene]
        return None if autosomal_recessive else self.bs1


def _number(where: str, value: object) -> float:
    # TOML booleans are Python bools, which are ints: refuse them explicitly.
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{where}: expected a number, got {value!r}")
    return float(value)


def _integer(where: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{where}: expected an integer, got {value!r}")
    return value


def _pair(where: str, value: object) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise TypeError(f"{where}: expected [low, high], got {value!r}")
    low, high = value
    return (_number(where, low), _number(where, high))


def _string(where: str, value: object) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{where}: expected a string, got {value!r}")
    return value


def _strings(where: str, value: object) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise TypeError(f"{where}: expected a list of strings, got {value!r}")
    return tuple(value)


_RELIABILITY = {
    "p_high": _number,
    "p_low": _number,
    "min_dp": _integer,
    "max_dp_factor": _number,
}
_FREQUENCY = {"ba1": _number, "bs1": _number}
_QC = {
    "min_mean_depth": _number,
    "depth_mapq": _integer,
    "min_fraction_20x": _number,
    "min_mapped_fraction": _number,
    "max_duplicate_fraction": _number,
    "min_median_insert": _number,
    "max_freemix": _number,
    "min_titv": _number,
    "max_titv": _number,
    "max_het_hom": _number,
}
_SEX = {"mapq": _integer}
_MITO = {"min_homoplasmy": _number, "heteroplasmy_floor": _number}
_HHV6 = {"min_ratio": _number, "max_ratio": _number, "background_factor": _number}
_ROH = {"min_run": _integer, "isodisomy_min": _integer}
_PREDICTED = {
    "revel": _number,
    "alphamissense": _number,
    "spliceai": _number,
    "min_pext": _number,
    "homopolymer": _integer,
}
_REPEATS = {"min_allele_reads": _integer}
_EXCLUSION = {"min_spanning": _integer}
_SMN = {"min_parascopy_qual": _number}
_YLEAF = {"min_qc_score": _number}
_DELETIONS = {
    "min_length": _integer,
    "long_from": _integer,
    "exon_margin": _integer,
    "max_gnomad_af": _number,
    "min_gnomad_an": _integer,
    "min_reciprocal_overlap": _number,
    "max_depth_ratio": _number,
    "depth_decides_from": _integer,
    "likely_below": _number,
}
_IDENTITY = {"min_concordance": _number, "min_compared": _integer}
_COVERAGE = {
    "mapq": _integer,
    "min_depth_per_copy": _integer,
    "cds_margin": _integer,
    "min_callable": _number,
}
_DOSAGE = {
    "tolerance": _number,
    "het_low": _number,
    "het_high": _number,
    "x_het_min": _number,
    "min_x_sites": _integer,
    "y_absent_max": _number,
    "calibrated_on": _string,
}
_AF_GUARD = {
    "min_het_vaf": _number,
    "clonal_genes": _strings,
    "calibrated_on": _string,
}
_INSERTIONS = {"max_distance": _integer}
_CNV = {
    "min_length": _integer,
    "band_from": _integer,
    "large_from": _integer,
    "bin_size": _integer,
}
_PARALOGS = {
    "min_qual": _number,
    "unique_loss_ratio": _pair,
    "gba_min_reads": _integer,
}
_LPA = {"max_kiv2_error": _number}
_LOCUS_REQUIRED = (
    "locus_id",
    "variant_id",
    "chrom",
    "motif",
    "normal_max",
    "pathogenic_min",
    "evidence",
    "topic",
)
_LOCUS_OPTIONAL = (
    "intermediate",
    "reduced_penetrance",
    "note",
    "rule",
    "validity",
    "recessive_at",
)
_FRACTIONS = (
    "min_fraction_20x",
    "min_mapped_fraction",
    "max_duplicate_fraction",
    "max_freemix",
)
_SECTIONS: dict[str, frozenset[str] | None] = {
    "reliability": frozenset(_RELIABILITY),
    "frequency": frozenset(_FREQUENCY),
    "bs1_by_gene": None,  # any gene symbol
    "par": None,  # any contig name
    "strata": None,  # any flag name
    "qc": frozenset(_QC),
    "sex": frozenset(_SEX),
    "mito": frozenset(_MITO),
    "hhv6": frozenset(_HHV6),
    "roh": frozenset(_ROH),
    "predicted": frozenset(_PREDICTED),
    "repeats": frozenset(_REPEATS),
    "repeat_loci": None,  # any gene symbol
    "exclusion": frozenset(_EXCLUSION),
    "smn": frozenset(_SMN),
    "yleaf": frozenset(_YLEAF),
    "deletions": frozenset(_DELETIONS),
    "identity": frozenset(_IDENTITY),
    "coverage": frozenset(_COVERAGE),
    "dosage": frozenset(_DOSAGE),
    "af_guard": frozenset(_AF_GUARD),
    "insertions": frozenset(_INSERTIONS),
    "cnv": frozenset(_CNV),
    "paralogs": frozenset(_PARALOGS),
    "lpa": frozenset(_LPA),
}


def load_grading(path: Path | None = None) -> Grading:
    """Read the grading config (default: config/grading.toml in this repository).

    Unknown sections or keys, and values of the wrong type, are errors.
    """
    path = path if path is not None else GRADING_FILE
    with path.open("rb") as f:
        data = tomllib.load(f)
    sections = _sections(path, data)
    deletions = DeletionThresholds(
        **_typed(path, "deletions", sections["deletions"], _DELETIONS)
    )
    reliability = _typed(path, "reliability", sections["reliability"], _RELIABILITY)
    frequency = _typed(path, "frequency", sections["frequency"], _FREQUENCY)
    grading = Grading(
        **reliability,
        par={
            contig: _intervals(f"{path}: [par] {contig}", intervals)
            for contig, intervals in sections["par"].items()
        },
        **frequency,
        bs1_by_gene={
            gene: _number(f"{path}: [bs1_by_gene] {gene}", value)
            for gene, value in sections.get("bs1_by_gene", {}).items()
        },
        strata={
            flag: _string(f"{path}: [strata] {flag}", bed)
            for flag, bed in sections["strata"].items()
        },
        qc=QcThresholds(**_typed(path, "qc", sections["qc"], _QC)),
        sex=SexThresholds(**_typed(path, "sex", sections["sex"], _SEX)),
        mito=MitoThresholds(**_typed(path, "mito", sections["mito"], _MITO)),
        hhv6=Hhv6Thresholds(**_typed(path, "hhv6", sections["hhv6"], _HHV6)),
        roh=RohThresholds(**_typed(path, "roh", sections["roh"], _ROH)),
        predicted=PredictedThresholds(
            **_typed(path, "predicted", sections["predicted"], _PREDICTED)
        ),
        repeats=RepeatThresholds(
            **_typed(path, "repeats", sections["repeats"], _REPEATS),
            loci={
                gene: _repeat_locus(path, gene, table)
                for gene, table in sections["repeat_loci"].items()
            },
        ),
        exclusion=ExclusionThresholds(
            **_typed(path, "exclusion", sections["exclusion"], _EXCLUSION)
        ),
        smn=SmnThresholds(**_typed(path, "smn", sections["smn"], _SMN)),
        yleaf=YleafThresholds(**_typed(path, "yleaf", sections["yleaf"], _YLEAF)),
        deletions=deletions,
        insertions=InsertionThresholds(
            **asdict(deletions),
            **_typed(path, "insertions", sections["insertions"], _INSERTIONS),
        ),
        identity=IdentityThresholds(
            **_typed(path, "identity", sections["identity"], _IDENTITY)
        ),
        coverage=CoverageThresholds(
            **_typed(path, "coverage", sections["coverage"], _COVERAGE)
        ),
        dosage=DosageThresholds(**_typed(path, "dosage", sections["dosage"], _DOSAGE)),
        af_guard=AfGuardThresholds(
            **_typed(path, "af_guard", sections["af_guard"], _AF_GUARD)
        ),
        cnv=CnvThresholds(**_typed(path, "cnv", sections["cnv"], _CNV)),
        paralogs=ParalogThresholds(
            **_typed(path, "paralogs", sections["paralogs"], _PARALOGS)
        ),
        lpa=LpaThresholds(**_typed(path, "lpa", sections["lpa"], _LPA)),
    )
    if not 0 < grading.p_high <= grading.p_low < 1:
        raise ValueError(f"{path}: need 0 < p_high <= p_low < 1")
    if not 0 < grading.bs1 <= grading.ba1 <= 1:
        raise ValueError(f"{path}: need 0 < bs1 <= ba1 <= 1")
    for gene, value in grading.bs1_by_gene.items():
        if not 0 < value <= grading.ba1:
            raise ValueError(f"{path}: [bs1_by_gene] {gene}: need 0 < BS1 <= ba1")
    _check_gates(path, grading.qc)
    if not 0 < grading.mito.min_homoplasmy <= 1:
        raise ValueError(f"{path}: [mito] min_homoplasmy: need a fraction in (0, 1]")
    # Its anchor comes from MIN_HETEROPLASMY_FLOOR, which an edited grading.toml cannot move
    # because the bound is code: a lower floor would grade mutserve's noise as heteroplasmy.
    if (
        not MIN_HETEROPLASMY_FLOOR
        <= grading.mito.heteroplasmy_floor
        < (grading.mito.min_homoplasmy)
    ):
        raise ValueError(
            f"{path}: [mito] heteroplasmy_floor: need at least"
            f" {MIN_HETEROPLASMY_FLOOR:g} and below min_homoplasmy"
        )
    h = grading.hhv6
    if not 0 < h.min_ratio <= 1 <= h.max_ratio or h.background_factor < 1:
        raise ValueError(
            f"{path}: [hhv6] need 0 < min_ratio <= 1 <= max_ratio (one copy is 1)"
            " and background_factor >= 1"
        )
    if not 0 < grading.roh.min_run <= grading.roh.isodisomy_min:
        raise ValueError(f"{path}: [roh] need 0 < min_run <= isodisomy_min")
    p = grading.predicted
    if not all(0 < v <= 1 for v in (p.revel, p.alphamissense, p.spliceai, p.min_pext)):
        raise ValueError(
            f"{path}: [predicted] scores and min_pext: need a fraction in (0, 1]"
        )
    # A run of one base is every single-base indel. Its anchor comes from this bound, which
    # an edited grading.toml cannot move because the bound is code, not read from the file.
    if p.homopolymer < 2:
        raise ValueError(f"{path}: [predicted] homopolymer: need at least 2 bases")
    if grading.repeats.min_allele_reads < 1:
        raise ValueError(f"{path}: [repeats] min_allele_reads: need at least 1")
    if grading.exclusion.min_spanning < MIN_SPANNING_FLOOR:
        raise ValueError(
            f"{path}: [exclusion] min_spanning: need at least {MIN_SPANNING_FLOOR}"
        )
    if grading.smn.min_parascopy_qual < 0:
        raise ValueError(
            f"{path}: [smn] min_parascopy_qual: need a non-negative number"
        )
    d = grading.deletions
    fractions = (
        d.max_gnomad_af,
        d.min_reciprocal_overlap,
        d.max_depth_ratio,
        d.likely_below,
    )
    if (
        not 0 < d.min_length < d.long_from <= d.depth_decides_from
        or min(d.exon_margin, d.min_gnomad_an) < 0
    ):
        raise ValueError(
            f"{path}: [deletions] need 0 < min_length < long_from <= depth_decides_from,"
            " exon_margin >= 0"
            " and min_gnomad_an >= 0"
        )
    # Exons at most twice the margin apart merge into one target: a kept deletion must not
    # fit wholly in the intron between them, or `screen.flank_only` reads it as an exon hit.
    if 2 * d.exon_margin >= d.min_length:
        raise ValueError(
            f"{path}: [deletions] need 2 x exon_margin < min_length, so no deletion the"
            " screen keeps fits in an intron inside a target"
        )
    if not all(0 < v <= 1 for v in fractions):
        raise ValueError(
            f"{path}: [deletions] rates and ratios: need a fraction in (0, 1]"
        )
    if not 0 < grading.identity.min_concordance <= 1:
        raise ValueError(
            f"{path}: [identity] min_concordance: need a fraction in (0, 1]"
        )
    if grading.identity.min_compared < 1:
        raise ValueError(f"{path}: [identity] min_compared: need at least 1")
    c = grading.coverage
    if min(c.mapq, c.cds_margin) < 0 or c.min_depth_per_copy < 1:
        raise ValueError(
            f"{path}: [coverage] need mapq >= 0, cds_margin >= 0 and"
            " min_depth_per_copy >= 1"
        )
    if not 0 < c.min_callable <= 1:
        raise ValueError(f"{path}: [coverage] min_callable: need a fraction in (0, 1]")
    _check_dosage(path, grading.dosage)
    a = grading.af_guard
    # A bound of one half or more would warn every balanced heterozygote.
    if not 0 < a.min_het_vaf < 0.5:
        raise ValueError(f"{path}: [af_guard] min_het_vaf: need 0 < min_het_vaf < 0.5")
    if not a.clonal_genes or len(set(a.clonal_genes)) != len(a.clonal_genes):
        raise ValueError(f"{path}: [af_guard] clonal_genes: need distinct gene names")
    # A negative distance would match no gnomAD-SV insertion, so every common one would
    # show. Its anchor comes from this bound, which an edited grading.toml cannot move
    # because the bound is code, not read from the file.
    if grading.insertions.max_distance < 0:
        raise ValueError(f"{path}: [insertions] max_distance: need at least 0 bp")
    c = grading.cnv
    if not 0 < c.min_length < c.band_from < c.large_from or c.bin_size < 1:
        raise ValueError(
            f"{path}: [cnv] need 0 < min_length < band_from < large_from and bin_size >= 1"
        )
    _check_paralogs(path, grading.paralogs, grading.lpa)
    return grading


def _check_paralogs(path: Path, p: ParalogThresholds, lpa: LpaThresholds) -> None:
    # Each bound's anchor comes from this function, which an edited grading.toml cannot
    # move because the bounds are code, not read from the file. A floor of 0 reads would
    # let a carrier call no read supports through; a band outside (0, 1], or low not below
    # high, would let a depth that a one-copy loss cannot give confirm it; an error of 0
    # withholds LPA on any difference, one above 1 passes a count off by more than the
    # truth itself.
    if p.min_qual < 0:
        raise ValueError(f"{path}: [paralogs] min_qual: need a non-negative number")
    low, high = p.unique_loss_ratio
    if not 0 < low < high <= 1:
        raise ValueError(
            f"{path}: [paralogs] unique_loss_ratio: need 0 < low < high <= 1"
        )
    if p.gba_min_reads < 1:
        raise ValueError(f"{path}: [paralogs] gba_min_reads: need at least 1 read")
    if not 0 < lpa.max_kiv2_error <= 1:
        raise ValueError(
            f"{path}: [lpa] max_kiv2_error: need a fraction of the count in (0, 1]"
        )


def _repeat_locus(path: Path, gene: str, table: object) -> RepeatLocus:
    where = f"{path}: [repeat_loci] {gene}"
    if not isinstance(table, dict):
        raise TypeError(f"{where}: expected an inline table")
    absent = [key for key in _LOCUS_REQUIRED if key not in table]
    if absent:
        raise ValueError(f"{where} lacks {absent}")
    extra = set(table) - set(_LOCUS_REQUIRED) - set(_LOCUS_OPTIONAL)
    if extra:
        raise ValueError(f"{where}: unknown keys {sorted(extra)}")
    locus = RepeatLocus(
        gene=gene,
        locus_id=_string(f"{where} locus_id", table["locus_id"]),
        variant_id=_string(f"{where} variant_id", table["variant_id"]),
        chrom=_string(f"{where} chrom", table["chrom"]),
        motif=_string(f"{where} motif", table["motif"]),
        normal_max=_integer(f"{where} normal_max", table["normal_max"]),
        intermediate=_span(f"{where} intermediate", table.get("intermediate")),
        reduced_penetrance=_span(
            f"{where} reduced_penetrance", table.get("reduced_penetrance")
        ),
        pathogenic_min=_integer(f"{where} pathogenic_min", table["pathogenic_min"]),
        evidence=_string(f"{where} evidence", table["evidence"]),
        topic=Topic(_string(f"{where} topic", table["topic"])),
        note=None if "note" not in table else _string(f"{where} note", table["note"]),
        rule=_choice(f"{where} rule", Rule, table.get("rule", "standard")),
        validity=None
        if "validity" not in table
        else _choice(f"{where} validity", Validity, table["validity"]),
        recessive_at=None
        if "recessive_at" not in table
        else _integer(f"{where} recessive_at", table["recessive_at"]),
    )
    if not 0 <= locus.normal_max < locus.pathogenic_min:
        raise ValueError(f"{where}: need 0 <= normal_max < pathogenic_min")
    _check_rule(where, locus)
    return locus


def _choice[E: (Rule, Validity)](where: str, kind: type[E], value: object) -> E:
    text = _string(where, value)
    try:
        return kind(text)
    except ValueError:
        raise ValueError(
            f"{where}: {text!r} is not one of {[str(m) for m in kind]}"
        ) from None


def _check_rule(where: str, locus: RepeatLocus) -> None:
    # The X-dosage rule reads one or two X copies, and the two-allele rules a diploid call:
    # a rule on the wrong chromosome would grade a haploid call as diploid or the reverse.
    # Its anchor comes from these chromosome checks, which an edited grading.toml cannot
    # move because they are code, not read from the file.
    if locus.rule is Rule.X_LINKED and locus.chrom != "chrX":
        raise ValueError(f"{where}: an x_linked locus must be on chrX")
    if locus.rule in (Rule.EXCLUSION, Rule.DISTINCT) and locus.chrom in (
        "chrX",
        "chrY",
    ):
        raise ValueError(f"{where}: an {locus.rule} locus must not be on chrX or chrY")
    # Both alleles at `recessive_at` must lie between the normal and the disease ranges.
    recessive = locus.recessive_at
    if (
        recessive is not None
        and not locus.normal_max < recessive < locus.pathogenic_min
    ):
        raise ValueError(
            f"{where}: recessive_at must lie above normal_max and below pathogenic_min"
        )


def _span(where: str, value: object) -> tuple[int, int] | None:
    if value is None:
        return None
    (start, end), *_ = _intervals(where, [value])
    return (start, end)


def _typed(
    path: Path, section: str, table: dict[str, Any], types: dict[str, Any]
) -> dict[str, Any]:
    # Every key must be present: a gate threshold never falls back silently. Its anchor
    # comes from `types`, this module's key table, which an edited grading.toml cannot move
    # because the table is code, not read from the file.
    absent = sorted(set(types) - set(table))
    if absent:
        raise ValueError(f"{path}: [{section}] lacks {absent}")
    return {
        key: types[key](f"{path}: [{section}] {key}", value)
        for key, value in table.items()
    }


def _check_gates(path: Path, qc: QcThresholds) -> None:
    for key in _FRACTIONS:
        if not 0 <= getattr(qc, key) <= 1:
            raise ValueError(f"{path}: [qc] {key}: need a fraction between 0 and 1")
    if not 0 < qc.min_titv < qc.max_titv:
        raise ValueError(f"{path}: [qc] need 0 < min_titv < max_titv")


def _check_dosage(path: Path, d: DosageThresholds) -> None:
    # A tolerance of a quarter or more would let one ratio be two copy counts.
    if not 0 < d.tolerance < 0.25:
        raise ValueError(f"{path}: [dosage] need 0 < tolerance < 0.25")
    if not 0 < d.het_low < d.het_high < 1 or not 0 <= d.x_het_min <= 1:
        raise ValueError(
            f"{path}: [dosage] need 0 < het_low < het_high < 1 and x_het_min a fraction"
        )
    if d.min_x_sites < 1:
        raise ValueError(f"{path}: [dosage] min_x_sites: need at least 1")
    if not 0 < d.y_absent_max < 0.5 - d.tolerance:
        raise ValueError(f"{path}: [dosage] need 0 < y_absent_max < 0.5 - tolerance")


def _sections(path: Path, data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    # Its anchor comes from this set of required sections, which an edited grading.toml
    # cannot move because the set is code, not read from the file.
    missing = {
        "reliability",
        "frequency",
        "par",
        "strata",
        "qc",
        "sex",
        "mito",
        "hhv6",
        "roh",
        "predicted",
        "repeats",
        "repeat_loci",
        "exclusion",
        "smn",
        "yleaf",
        "deletions",
        "identity",
        "coverage",
        "dosage",
        "af_guard",
        "insertions",
        "cnv",
        "paralogs",
        "lpa",
    } - set(data)
    if missing:
        raise ValueError(f"{path}: missing sections {sorted(missing)}")
    unknown = set(data) - set(_SECTIONS)
    if unknown:
        raise ValueError(f"{path}: unknown sections {sorted(unknown)}")
    sections: dict[str, dict[str, Any]] = {}
    for section, table in data.items():
        if not isinstance(table, dict):
            raise TypeError(f"{path}: [{section}] must be a table")
        keys = _SECTIONS[section]
        extra = set() if keys is None else set(table) - keys
        if extra:
            raise ValueError(f"{path}: unknown keys in [{section}]: {sorted(extra)}")
        sections[section] = table
    return sections


def _intervals(where: str, value: object) -> list[tuple[int, int]]:
    if not isinstance(value, list):
        raise TypeError(f"{where}: expected a list of [start, end] pairs")
    intervals: list[tuple[int, int]] = []
    for pair in value:
        if not isinstance(pair, list) or len(pair) != 2:
            raise TypeError(f"{where}: expected [start, end], got {pair!r}")
        start, end = _integer(where, pair[0]), _integer(where, pair[1])
        if start > end:
            raise ValueError(f"{where}: interval {pair!r} starts after it ends")
        intervals.append((start, end))
    return intervals
