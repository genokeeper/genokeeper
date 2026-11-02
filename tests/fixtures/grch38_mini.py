"""A small GRCh38-shaped reference (chr1, chrX, chrY of random sequence) and reads for it,
written into a test's tmp_path: paired FASTQ named as a DNBSEQ-T7 run names its reads, and
coordinate-sorted CRAMs of chosen depth per contig. Nothing here is committed as a genomic
file; the bwa index is built in the test session.
"""

import gzip
import random
import subprocess
from collections.abc import Mapping
from pathlib import Path

import pysam

CONTIGS = {"chr1": 30_000, "chrX": 12_000, "chrY": 12_000}
READ = 100
INSERT = 300
FLOWCELL = "E100000001"  # invented; the E prefix is the original DNBSEQ-T7's
# Phred+33 qualities over many levels, as an unbinned T7 run writes them.
QUALITIES = "".join(chr(33 + q) for q in range(2, 42))

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def write_reference(path: Path, *, seed: int = 38) -> Path:
    """Random sequence for each of CONTIGS, with its .fai."""
    rng = random.Random(seed)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for name, length in CONTIGS.items():
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            f.write(f">{name}\n")
            for i in range(0, length, 60):
                f.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def bwa_index(fasta: Path) -> Path:
    """`bwa index` the FASTA in place (<fasta>.bwt and the rest); needs bwa on PATH."""
    subprocess.run(["bwa", "index", str(fasta)], check=True, capture_output=True)
    return fasta


def read_name(i: int) -> str:
    """A DNBSEQ-T7 read name: flowcell, lane, column, row, read number."""
    return f"{FLOWCELL}L1C{1 + i // 1000:03d}R{1 + i % 1000:03d}{i:07d}"


def write_fastq_pair(
    folder: Path,
    fasta: Path,
    *,
    pairs: int = 400,
    names: str | None = None,
    qualities: str = QUALITIES,
    seed: int = 7,
) -> tuple[Path, Path]:
    """Gzipped paired FASTQ of `pairs` fragments from chr1 (read 2 the reverse complement
    of the fragment's end). `names` replaces the flowcell in every read name; quality
    characters are drawn from `qualities`. The last pair repeats the first: one duplicate."""
    rng = random.Random(seed)
    with pysam.FastaFile(str(fasta)) as ref:
        chrom = ref.fetch("chr1")
    fq1, fq2 = folder / "reads_1.fq.gz", folder / "reads_2.fq.gz"
    # Distinct starts: with a fixed insert size, two pairs at one start are duplicates.
    starts = rng.sample(range(len(chrom) - INSERT), pairs - 1)
    starts.append(starts[0])
    with gzip.open(fq1, "wt") as r1, gzip.open(fq2, "wt") as r2:
        for i, start in enumerate(starts):
            name = read_name(i)
            if names is not None:
                name = name.replace(FLOWCELL, names)
            fragment = chrom[start : start + INSERT]
            for out, seq, mate in (
                (r1, fragment[:READ], 1),
                (r2, _revcomp(fragment[-READ:]), 2),
            ):
                qual = "".join(rng.choice(qualities) for _ in range(READ))
                out.write(f"@{name}/{mate}\n{seq}\n+\n{qual}\n")
    return fq1, fq2


def write_cram(
    path: Path,
    fasta: Path,
    depths: Mapping[str, float],
    *,
    mapq: int = 60,
    seed: int = 11,
) -> Path:
    """A coordinate-sorted, indexed CRAM of proper pairs (insert INSERT) tiling each contig
    of `depths` at about that depth; needs samtools on PATH."""
    rng = random.Random(seed)
    unsorted = path.with_name(f"{path.name}.unsorted.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS.items()],
        "RG": [{"ID": "mini", "SM": "mini", "PL": "DNBSEQ"}],
    }
    with (
        pysam.FastaFile(str(fasta)) as ref,
        pysam.AlignmentFile(str(unsorted), "wb", header=header) as out,
    ):
        n = 0
        for contig, depth in depths.items():
            length = CONTIGS[contig]
            seq = ref.fetch(contig)
            for _ in range(round(depth * length / (2 * READ))):
                start = rng.randrange(0, length - INSERT)
                name = read_name(n)
                n += 1
                for flag, pos, mate_pos, tlen in (
                    (99, start, start + INSERT - READ, INSERT),
                    (147, start + INSERT - READ, start, -INSERT),
                ):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = name
                    read.flag = flag
                    read.reference_name = contig
                    read.reference_start = pos
                    read.mapping_quality = mapq
                    read.cigarstring = f"{READ}M"
                    read.next_reference_name = contig
                    read.next_reference_start = mate_pos
                    read.template_length = tlen
                    read.query_sequence = seq[pos : pos + READ]
                    read.query_qualities = pysam.qualitystring_to_array("I" * READ)
                    read.set_tag("RG", "mini")
                    out.write(read)
    subprocess.run(
        [
            "samtools",
            "sort",
            "-O",
            "cram",
            "--reference",
            str(fasta),
            "--write-index",
            "-o",
            str(path),
            str(unsorted),
        ],
        check=True,
        capture_output=True,
    )
    unsorted.unlink()
    return path


def write_regions(path: Path, rows: list[tuple[str, int, int, str]]) -> Path:
    """A named BED (chrom, start, end, name), as genokeeper.regions derives them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in rows))
    return path
