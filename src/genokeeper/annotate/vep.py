"""Ensembl VEP and SpliceAI on the run's variants (track design, slice 1e).

`annotate` runs VEP 116 (`containers.IMAGES["vep"]`, Ensembl's own image, its plugins baked
in) offline on the lab's indexed GRCh38 cache, with MANE Select and MANE Plus Clinical
flags, HGVS from the lab's FASTA, exon numbers and total lengths, and three plugins on files
the bootstrap fetched whole: REVEL, AlphaMissense and SpliceVault. It keeps one record per
variant (`VepRecord`): its MANE Select row where that has a consequence of interest, else
the row of its most severe one, a MANE transcript first among equals. gnomAD's base-level
pext (`PextTable`) is then read for each loss-of-function record. SpliceAI (`containers.IMAGES["spliceai"]`, built from
docker/spliceai.Dockerfile) runs locally on the records of `genes_for_spliceai` only: its
precomputed scores need an Illumina account.

VEP runs behind the core gate as `dna:annotate`, each of SpliceAI's shards (`run_spliceai`)
as `dna:spliceai`. The containment, against a tool's own network call carrying the user's
variants: its anchor comes from `run_container`'s `--network none`, which neither VEP nor
SpliceAI can change because the network namespace is set before their process starts; VEP's
`--offline` only keeps it from trying.

`python -m genokeeper.annotate.vep layout` prints where the bootstrap puts each file
(`VepFiles.required`); `nmd-bed BB OUT` converts UCSC's NMD Escape track for MANE from bigBed to BED
with UCSC's bigBedToBed in its pinned image, offline.
"""

import gzip
import re
import shutil
import sys
from collections.abc import Collection, Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

import pysam

from genokeeper.config import load_paths
from genokeeper.containers import PLATFORMS, run_container
from genokeeper.tools import run_tool
from genokeeper.variants.key import make_vid

VEP_VERSION = "116.3"
CACHE_VERSION = "116"
LABEL = "dna:annotate"
# Where the bootstrap puts slice 1e's files, under <lab>/vep/116/.
CACHE = "cache"
CACHE_INFO = "cache/homo_sapiens/116_GRCh38/info.txt"
REVEL_FILE = "revel_v1.3_grch38.tsv.gz"
ALPHAMISSENSE_FILE = "AlphaMissense_hg38.tsv.gz"
SPLICEVAULT_FILE = "SpliceVault_data_GRCh38.tsv.gz"
PEXT_FILE = "pext_base_level.bed.gz"
PEXT_MAX_FILE = "pext_gene_max.tsv"
NMD_FILE = "nmdEscMane.bed"

VEP_INPUT = "input.vcf"
VEP_OUTPUT = "vep.tsv.gz"
SPLICEAI_INPUT = "spliceai_input.vcf"
SPLICEAI_OUTPUT = "spliceai.vcf"
# Where `run_spliceai` puts its shards' inputs and outputs while they run.
SPLICEAI_SHARDS = "spliceai-shards"
SPLICEAI_LABEL = "dna:spliceai"
# The cores one of several SpliceAI processes keeps busy; `run_spliceai` runs one per this
# many of the run's threads. One process scores one allele at a time, five Keras `predict`
# calls on one sequence each, and holds about 5 cores whatever it is given. On 1,757
# of HG003's sites (2026-10-29) one process took 1224 s at about 520%, two 741 s,
# three 649 s, and four 597 s at about 250% each.
SPLICEAI_CORES = 3

# The consequences grading reads (Sequence Ontology terms as VEP writes them).
LOF_CONSEQUENCES = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
    }
)
START_LOST = "start_lost"
MISSENSE = "missense_variant"
# The columns VEP writes (its --fields).
FIELDS = (
    "Uploaded_variation",
    "Gene",
    "SYMBOL",
    "Feature",
    "Feature_type",
    "Consequence",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "EXON",
    "CDS_position",
    "HGVSc",
    "HGVSp",
    "REVEL",
    "am_pathogenicity",
    "SpliceVault_top_events",
)
# A pext value is read at the variant's base, else at the nearest base of its gene within
# this many bases: a ±1/2 splice variant sits in the intron, the base it spoils beside it.
PEXT_WINDOW = 3


@dataclass(frozen=True)
class VepFiles:
    """slice 1e's files in the lab (`vep_files`), each tabix-indexed where VEP or the
    pext lookup reads it by region."""

    folder: Path
    revel: Path
    alphamissense: Path
    splicevault: Path
    pext: Path
    pext_max: Path
    nmd: Path

    def required(self) -> list[Path]:
        """Every file a run reads: the cache's info.txt for the cache, each indexed file
        with its .tbi."""
        indexed = [self.revel, self.alphamissense, self.splicevault, self.pext]
        return [
            self.folder / CACHE_INFO,
            *(p for f in indexed for p in (f, f.with_name(f"{f.name}.tbi"))),
            self.pext_max,
            self.nmd,
        ]


def vep_files(lab: Path) -> VepFiles:
    """`<lab>/vep/116/`, where the bootstrap's vep step writes slice 1e's files."""
    folder = lab / "vep" / CACHE_VERSION
    return VepFiles(
        folder=folder,
        revel=folder / REVEL_FILE,
        alphamissense=folder / ALPHAMISSENSE_FILE,
        splicevault=folder / SPLICEVAULT_FILE,
        pext=folder / PEXT_FILE,
        pext_max=folder / PEXT_MAX_FILE,
        nmd=folder / NMD_FILE,
    )


@dataclass(frozen=True)
class VepRecord:
    """One variant's annotation on its chosen transcript (`parse_vep`).

    `consequence` is VEP's `&`-joined terms; `mane` says the transcript is MANE Select or
    MANE Plus Clinical; `exon` is VEP's `n/total` and `cds_length` the transcript's CDS
    (stop codon included). `revel` and `alphamissense` score a missense allele, `spliceai` is
    SpliceAI's largest delta score on the variant's gene, `splicevault_in_frame` whether
    any of SpliceVault's top four events at a lost splice site keeps the frame (None: no
    SpliceVault data at the site), and `pext` the base's mean pext over its gene's
    largest (None: gnomAD has no pext there)."""

    vid: str
    gene: str
    gene_id: str
    transcript: str
    consequence: str
    mane: bool
    hgvsc: str | None
    hgvsp: str | None
    exon: str | None
    cds_length: int | None
    lof: bool
    revel: float | None
    alphamissense: float | None
    spliceai: float | None
    splicevault_in_frame: bool | None
    pext: float | None

    @property
    def terms(self) -> frozenset[str]:
        return frozenset(self.consequence.split("&"))


def vep_command(*, fasta: Path, files: VepFiles, threads: int) -> list[str]:
    """VEP's command in its container: /work holds the input and output, /vep the lab's
    slice 1e folder (`files`), /ref the FASTA's folder."""
    return [
        "vep",
        "--input_file", f"/work/{VEP_INPUT}", "--format", "vcf",
        "--output_file", f"/work/{VEP_OUTPUT}", "--tab", "--compress_output", "gzip",
        "--fields", ",".join(FIELDS),
        "--offline", "--cache", "--dir_cache", f"/vep/{CACHE}",
        "--cache_version", CACHE_VERSION, "--species", "homo_sapiens",
        "--assembly", "GRCh38", "--fasta", f"/ref/{fasta.name}",
        "--mane", "--hgvs", "--numbers", "--total_length", "--symbol",
        "--dir_plugins", "/plugins",
        "--plugin", f"REVEL,file=/vep/{files.revel.name}",
        "--plugin", f"AlphaMissense,file=/vep/{files.alphamissense.name}",
        "--plugin", f"SpliceVault,file=/vep/{files.splicevault.name}",
        "--fork", str(max(1, threads)), "--no_stats", "--force_overwrite",
    ]  # fmt: skip


def _value(text: str) -> str | None:
    return None if text in ("", "-") else text


def _number(text: str) -> float | None:
    value = _value(text)
    return None if value is None else float(value)


def _total(text: str) -> int | None:
    value = _value(text)
    tail = None if value is None or "/" not in value else value.split("/")[1]
    return int(tail) if tail and tail.isdigit() else None


def _in_frame(text: str) -> bool | None:
    """Whether a SpliceVault top event keeps the frame: each event is
    `TopN:type:impact:percent:frame`, frame `inFrame` or `Frameshift`."""
    value = _value(text)
    if value is None:
        return None
    events = [e for e in re.split(r"[,|&]", value) if e]
    return any(e.rsplit(":", 1)[-1].lower() == "inframe" for e in events)


def _rank(terms: frozenset[str]) -> int:
    if terms & LOF_CONSEQUENCES:
        return 3
    if START_LOST in terms:
        return 2
    return 1 if MISSENSE in terms else 0


def _rows(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as f:
        header: list[str] | None = None
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [fields[0].lstrip("#"), *fields[1:]]
                continue
            if header is None:
                raise ValueError(f"{path}: a row before VEP's column header")
            yield dict(zip(header, fields, strict=True))


def parse_vep(path: Path, *, keep: Collection[str]) -> list[VepRecord]:
    """One record per variant of VEP's tab output at `path` (in order of the variants'
    first row): its MANE Select row when that row's consequence is one the grading reads (a
    loss of function, start-loss or missense; design decision, 2026-10-28: a minor isoform's
    LoF never hides the MANE missense, and a LoF only on a non-MANE isoform still reaches
    gate 2); else its row of the highest rank (a loss-of-function consequence, then
    start-loss, then missense, then anything else), among equals a gene in `keep` first,
    then MANE Select, then MANE Plus Clinical. A variant with no row of rank above 0 is kept only in a gene of `keep`, for
    SpliceAI."""
    best: dict[str, tuple[tuple[bool, int, bool, bool, bool], dict[str, str]]] = {}
    for row in _rows(path):
        if row["Feature_type"] != "Transcript" or _value(row["SYMBOL"]) is None:
            continue
        rank = _rank(frozenset(row["Consequence"].split("&")))
        listed = row["SYMBOL"] in keep
        if rank == 0 and not listed:
            continue
        mane_select = _value(row["MANE_SELECT"]) is not None
        key = (
            mane_select and rank > 0,
            rank,
            listed,
            mane_select,
            _value(row["MANE_PLUS_CLINICAL"]) is not None,
        )
        vid = row["Uploaded_variation"]
        if vid not in best or key > best[vid][0]:
            best[vid] = (key, row)
    return [_record(row) for _, row in best.values()]


def _record(row: dict[str, str]) -> VepRecord:
    terms = frozenset(row["Consequence"].split("&"))
    return VepRecord(
        vid=row["Uploaded_variation"],
        gene=row["SYMBOL"],
        gene_id=row["Gene"],
        transcript=row["Feature"],
        consequence=row["Consequence"],
        mane=_value(row["MANE_SELECT"]) is not None
        or _value(row["MANE_PLUS_CLINICAL"]) is not None,
        hgvsc=_value(row["HGVSc"]),
        hgvsp=_value(row["HGVSp"]),
        exon=_value(row["EXON"]),
        cds_length=_total(row["CDS_position"]),
        lof=bool(terms & LOF_CONSEQUENCES),
        revel=_number(row["REVEL"]),
        alphamissense=_number(row["am_pathogenicity"]),
        spliceai=None,
        splicevault_in_frame=_in_frame(row["SpliceVault_top_events"]),
        pext=None,
    )


def locus(vid: str) -> tuple[str, int, str, str]:
    """A vid's (chrom, pos, ref, alt)."""
    _, chrom, pos, ref, alt = vid.split(":")
    return chrom, int(pos), ref, alt


class PextTable:
    """gnomAD v4.1's base-level pext (GTEx v10) as the bootstrap derived it: a tabix-indexed
    BED of `chrom, start, end, gene_id, exp_prop_mean`, and each gene's largest mean."""

    def __init__(self, bed: Path, gene_max: Path) -> None:
        self.bed = bed
        self.gene_max: dict[str, float] = {}
        for line in gene_max.read_text().splitlines():
            gene, value = line.split("\t")
            self.gene_max[gene] = float(value)

    def relative(self, chrom: str, pos: int, gene_id: str) -> float | None:
        """The mean pext at `pos` (else at the nearest base of the gene within
        `PEXT_WINDOW`) over the gene's largest; None where gnomAD has none."""
        top = self.gene_max.get(gene_id)
        if not top:
            return None
        with pysam.TabixFile(str(self.bed)) as tabix:
            if chrom not in tabix.contigs:
                return None
            found: list[tuple[int, float]] = []
            for line in tabix.fetch(
                chrom, max(0, pos - 1 - PEXT_WINDOW), pos + PEXT_WINDOW
            ):
                _, _, end, gene, value = line.split("\t")
                if gene == gene_id:
                    found.append((abs(int(end) - pos), float(value)))
        return None if not found else min(found)[1] / top


def fasta_contigs(fasta: Path) -> dict[str, int]:
    """The FASTA's contigs and their lengths, in its order."""
    with pysam.FastaFile(str(fasta)) as f:
        return dict(zip(f.references, f.lengths, strict=True))


def write_sites(path: Path, vids: Sequence[str], fasta: Path) -> None:
    """A sites-only VCF of `vids` (sorted by the FASTA's contig order), each record's ID
    its vid."""
    contigs = fasta_contigs(fasta)
    order = {name: i for i, name in enumerate(contigs)}
    loci = sorted(
        (locus(v) for v in vids), key=lambda a: (order[a[0]], a[1], a[2], a[3])
    )
    header = "##fileformat=VCFv4.2\n" + "".join(
        f"##contig=<ID={name},length={length}>\n" for name, length in contigs.items()
    )
    path.write_text(
        header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"{c}\t{p}\t{make_vid(c, p, r, a)}\t{r}\t{a}\t.\t.\t.\n"
            for c, p, r, a in loci
        )
    )


def vep_input(vcf: Path, out: Path, *, fasta: Path) -> list[str]:
    """VEP's input at `out`, a sites-only VCF with each vid as its ID: every PASS allele
    the sample carries in DeepVariant's VCF `vcf`, normalised and split against `fasta`
    (`bcftools norm -f fasta -m -both`), whatever ClinVar says of it (the run filters the
    classified ones at grading, so a new ClinVar release does not re-key VEP). Returns the
    vids it wrote. The intermediate files sit beside `out`, which holds the user's
    variants: the caller puts it under the run's work directory."""
    passed = out.with_name(f"{out.name}.pass.bcf")
    normalised = out.with_name(f"{out.name}.norm.vcf")
    run_tool(
        [
            "bcftools",
            "view",
            "--no-version",
            "-f",
            "PASS",
            "-Ob",
            "-o",
            str(passed),
            str(vcf),
        ]
    )
    run_tool(
        [
            "bcftools", "norm", "--no-version", "-f", str(fasta), "-m", "-both",
            "-Ov", "-o", str(normalised), str(passed),
        ]
    )  # fmt: skip
    vids: list[str] = []
    with pysam.VariantFile(str(normalised)) as found:
        for rec in found:
            alt = (rec.alts or ("",))[0]
            carried = 1 in (rec.samples[0]["GT"] or ())
            if rec.ref is None or not carried or not alt.isalpha():
                continue
            vids.append(make_vid(rec.chrom, rec.pos, rec.ref, alt))
    passed.unlink()
    normalised.unlink()
    given = list(dict.fromkeys(vids))
    write_sites(out, given, fasta)
    return given


def parse_spliceai(path: Path) -> dict[tuple[str, str], float]:
    """SpliceAI's largest delta score per (vid, gene) from its output VCF: its INFO
    `SpliceAI=ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`, one entry per
    gene, comma-separated."""
    scores: dict[tuple[str, str], float] = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        chrom, pos, _, ref, alt, _, _, info = line.split("\t")[:8]
        vid = make_vid(chrom, int(pos), ref, alt)
        for item in info.split(";"):
            if not item.startswith("SpliceAI="):
                continue
            for entry in item.removeprefix("SpliceAI=").split(","):
                fields = entry.split("|")
                if len(fields) < 6 or fields[0] != alt:
                    continue
                deltas = [float(x) for x in fields[2:6] if x not in ("", ".")]
                if deltas:
                    key = (vid, fields[1])
                    scores[key] = max(scores.get(key, 0.0), *deltas)
    return scores


def run_spliceai(out_dir: Path, fasta: Path, *, threads: int) -> None:
    """SpliceAI on `out_dir/SPLICEAI_INPUT`, written to `out_dir/SPLICEAI_OUTPUT`: the
    input's records split in order into `threads // SPLICEAI_CORES` shards (at least
    one, none empty) under `out_dir/SPLICEAI_SHARDS/`, one container each behind the gate as
    `SPLICEAI_LABEL`, all started at once; their outputs joined in that order, the first
    one's header once, and the shards removed. SpliceAI scores each record alone and writes
    the records in input order under the input's header and one INFO line, so the output is
    the bytes one process writes (four shards on HG003) and the shard count stays out
    of the step's key. A shard that fails fails the step once the running shards finish."""
    lines = (out_dir / SPLICEAI_INPUT).read_text().splitlines(keepends=True)
    header = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    count = max(1, min(len(body), threads // SPLICEAI_CORES))
    bounds = [len(body) * i // count for i in range(count + 1)]
    folder = out_dir / SPLICEAI_SHARDS
    folder.mkdir(exist_ok=True)
    for i in range(count):
        (folder / f"{i}.vcf").write_text(
            "".join(header + body[bounds[i] : bounds[i + 1]])
        )

    def shard(i: int) -> list[str]:
        run_container(
            "spliceai",
            [
                "spliceai", "-I", f"/work/{i}.vcf", "-O", f"/work/{i}.out.vcf",
                "-R", f"/ref/{fasta.name}", "-A", "grch38",
            ],
            mounts=[(folder, "/work", False), (fasta.parent, "/ref", True)],
            label=SPLICEAI_LABEL,
        )  # fmt: skip
        return (folder / f"{i}.out.vcf").read_text().splitlines(keepends=True)

    pool = ThreadPoolExecutor(max_workers=count)
    try:
        outputs = list(pool.map(shard, range(count)))
    finally:
        pool.shutdown(cancel_futures=True)
    joined = [line for line in outputs[0] if line.startswith("#")] + [
        line for output in outputs for line in output if not line.startswith("#")
    ]
    (out_dir / SPLICEAI_OUTPUT).write_text("".join(joined))
    shutil.rmtree(folder)


def annotate(
    out_dir: Path,
    *,
    genes_for_spliceai: Collection[str],
    files: VepFiles,
    fasta: Path,
    threads: int,
) -> list[VepRecord]:
    """VEP's records (`parse_vep`) for the sites-only VCF `out_dir/VEP_INPUT`
    (`vep_input`), VEP writing its output beside it; each loss-of-function record's
    `pext`; and SpliceAI's score on the records of `genes_for_spliceai` (SpliceAI is not
    run when there are none). `out_dir` holds the user's variants: the caller puts it
    under the run's work directory."""
    mounts = [
        (out_dir, "/work", False),
        (files.folder, "/vep", True),
        (fasta.parent, "/ref", True),
    ]
    run_container(
        "vep",
        vep_command(fasta=fasta, files=files, threads=threads),
        mounts=mounts,
        label=LABEL,
    )
    records = parse_vep(out_dir / VEP_OUTPUT, keep=genes_for_spliceai)
    pext = PextTable(files.pext, files.pext_max)
    records = [
        replace(r, pext=pext.relative(*locus(r.vid)[:2], r.gene_id)) if r.lof else r
        for r in records
    ]
    listed = [r.vid for r in records if r.gene in genes_for_spliceai]
    if not listed:
        return records
    write_sites(out_dir / SPLICEAI_INPUT, listed, fasta)
    run_spliceai(out_dir, fasta, threads=threads)
    scores = parse_spliceai(out_dir / SPLICEAI_OUTPUT)
    return [
        replace(r, spliceai=scores.get((r.vid, r.gene)))
        if r.gene in genes_for_spliceai
        else r
        for r in records
    ]


def nmd_bed(bigbed: Path, out: Path) -> None:
    """UCSC's NMD Escape track (`nmdEscMane.bb`) as BED at `out`, by bigBedToBed in its
    pinned image, offline."""
    out.parent.mkdir(parents=True, exist_ok=True)
    run_container(
        "ucsc_bigbedtobed",
        ["bigBedToBed", f"/in/{bigbed.name}", f"/out/{out.name}"],
        mounts=[(bigbed.parent, "/in", True), (out.parent, "/out", False)],
        platform=PLATFORMS["ucsc_bigbedtobed"],
    )


def main(argv: Sequence[str] | None = None) -> int:
    """`layout` (`VepFiles.required()` under $GENOKEEPER_HOME, one per line, in order) or
    `nmd-bed BB OUT`, for scripts/bootstrap.sh."""
    args = list(sys.argv[1:] if argv is None else argv)
    if args == ["layout"]:
        print("\n".join(map(str, vep_files(load_paths().lab).required())))
        return 0
    if len(args) == 3 and args[0] == "nmd-bed":
        nmd_bed(Path(args[1]), Path(args[2]))
        return 0
    print(
        "usage: python -m genokeeper.annotate.vep {layout | nmd-bed BB OUT}",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
