"""Slice 3d's spike-ins (`genokeeper.calibrate.spikein`): which events are spiked, and a
spiked CRAM written from a small BAM the test writes with pysam, through native samtools."""

import itertools
import random
import shutil
import zlib
from collections import Counter
from pathlib import Path

import pysam
import pytest

from genokeeper import tools
from genokeeper.calibrate.spikein import pick_spikes, spike_names, write_spiked
from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import Deletion
from genokeeper.regions import Interval

SAMTOOLS = shutil.which("samtools")


def _genes(n: int) -> dict[str, list[Target]]:
    """`n` genes 300 kb apart on chr1, each of six 200 bp exons 3 kb apart."""
    genes: dict[str, list[Target]] = {}
    for g in range(n):
        start = 1_000_000 + g * 300_000
        genes[f"G{g}"] = [
            Target("chr1", start + i * 3_200, start + i * 3_200 + 200, (i + 1,))
            for i in range(6)
        ]
    return genes


def _overlaps(d: Deletion, intervals: list[Interval]) -> bool:
    return any(c == d.chrom and s < d.end and d.pos < e for c, s, e in intervals)


def test_pick_spikes_rules() -> None:
    targets = _genes(300)
    mask = [
        ("chr1", 1_000_000 + g * 300_000, 1_000_000 + g * 300_000 + 20_000)
        for g in range(0, 300, 7)
    ]
    owned = [
        ("chr1", 1_000_000 + g * 300_000 + 5_000, 1_000_000 + g * 300_000 + 6_000)
        for g in range(3, 300, 5)
    ]
    losses, gains = pick_spikes(targets, mask, owned)
    assert (losses, gains) == pick_spikes(targets, mask, owned, n_per=20, seed=3)
    assert pick_spikes(targets, mask, owned, seed=4) != (losses, gains)
    for found, svtype in ((losses, "DEL"), (gains, "DUP")):
        assert {d.svtype for d in found} == {svtype}
        exons = Counter(
            len(d.callers["spikein"].split(":")[1].split(",")) for d in found
        )
        assert exons == {1: 20, 2: 20, 3: 20}
    spikes = sorted(losses + gains, key=lambda d: d.pos)
    for d in spikes:
        assert d.length >= 1000
        assert not _overlaps(d, mask) and not _overlaps(d, owned)
        gene, numbers = d.callers["spikein"].split(":")
        listed = [int(n) for n in numbers.split(",")]
        first, last = listed[0], listed[-1]
        exons = targets[gene][first - 1], targets[gene][last - 1]
        # From the middle of the intron before the first exon to the middle after the last.
        assert d.pos == (targets[gene][first - 2].end + exons[0].start) // 2
        assert d.end == (exons[1].end + targets[gene][last].start) // 2
    for a, b in itertools.pairwise(spikes):
        assert b.pos - a.end >= 100_000


def test_pick_spikes_caps_the_flanks() -> None:
    """An intron of 100 kb: the event reaches at most 20 kb past its exons."""
    targets = {
        "G": [
            Target("chr1", 1_000_000, 1_000_200, (1,)),
            Target("chr1", 1_100_000, 1_100_200, (2,)),
            Target("chr1", 1_200_000, 1_200_200, (3,)),
        ]
    }
    losses, gains = pick_spikes(targets, [], [], n_per=1)
    (d,) = [*losses, *gains]
    assert (d.pos, d.end) == (1_080_000, 1_120_200)


# --- the spiked CRAM ----------------------------------------------------------------------

LENGTH = 60_000
LOSS = Deletion("chr1", 10_000, 14_000, "N", "0/1", {"spikein": "G:2"}, None)
GAIN = Deletion(
    "chr1", 30_000, 34_000, "N", "0/1", {"spikein": "G:2"}, None, svtype="DUP"
)


def _even(name: str) -> bool:
    return zlib.crc32(name.encode()) % 2 == 0


def _write_bam(path: Path, sequence: str) -> list[tuple[str, int, int]]:
    """Pairs of 100 bp reads: forty inside the loss, forty inside the gain, ten with one
    mate inside the loss and the other 30 kb away, forty elsewhere. Returns (name, first
    mate's start, second mate's start)."""
    rng = random.Random(5)
    pairs: list[tuple[str, int, int]] = []
    for i in range(40):
        start = rng.randrange(10_000, 13_700)
        pairs.append((f"loss{i}", start, start + 150))
        start = rng.randrange(30_000, 33_700)
        pairs.append((f"gain{i}", start, start + 150))
        start = rng.randrange(40_000, 55_000)
        pairs.append((f"other{i}", start, start + 150))
    for i in range(10):
        start = rng.randrange(11_000, 12_000)
        pairs.append((f"far{i}", start, start + 30_000))
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": LENGTH}],
    }
    reads: list[pysam.AlignedSegment] = []
    with pysam.AlignmentFile(str(path) + ".unsorted.bam", "wb", header=header) as out:
        for name, one, two in pairs:
            for first, (start, mate) in enumerate(((one, two), (two, one))):
                read = pysam.AlignedSegment(out.header)
                read.query_name = name
                read.reference_id = 0
                read.reference_start = start
                read.next_reference_id = 0
                read.next_reference_start = mate
                read.mapping_quality = 60
                read.cigarstring = "100M"
                read.query_sequence = sequence[start : start + 100]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.flag = 1 | 2 | (64 if first == 0 else 128) | (16 if first else 32)
                read.template_length = (two + 100 - one) * (1 if first == 0 else -1)
                reads.append(read)
                out.write(read)
    pysam.sort("-o", str(path), str(path) + ".unsorted.bam")
    pysam.index(str(path))
    return pairs


def _reads(path: Path, fasta: Path) -> Counter[tuple[str, int, int, str]]:
    with pysam.AlignmentFile(str(path), reference_filename=str(fasta)) as f:
        return Counter(
            (r.query_name or "", r.flag, r.reference_start, r.query_sequence or "")
            for r in f.fetch(until_eof=True)
        )


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools is not on PATH")
def test_write_spiked_depth(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The pairs over a loss with an even crc32 are gone, mates included wherever they
    lie; those over a gain appear twice, the copies `<name>:spike`; every other read once.
    The whole unit runs behind one gate, `dna:spikein` (a fake gate records its label)."""
    gate = tmp_path / "gate"
    gate.write_text(
        f'#!/bin/sh\necho "$1" >> {tmp_path / "gated.txt"}\nshift 2\nexec "$@"\n'
    )
    gate.chmod(0o755)
    monkeypatch.setattr(tools, "GATE", gate)
    sequence = "".join(random.Random(9).choice("ACGT") for _ in range(LENGTH))
    fasta = tmp_path / "ref.fa"
    fasta.write_text(
        ">chr1\n" + "\n".join(sequence[i : i + 60] for i in range(0, LENGTH, 60)) + "\n"
    )
    pysam.faidx(str(fasta))
    bam = tmp_path / "in.bam"
    pairs = _write_bam(bam, sequence)
    drop, duplicate = spike_names(bam, fasta, [LOSS], [GAIN])
    assert drop == {
        n for n, _, _ in pairs if n.startswith(("loss", "far")) and _even(n)
    }
    assert duplicate == {n for n, _, _ in pairs if n.startswith("gain") and _even(n)}
    assert drop and duplicate and any(n.startswith("far") for n in drop)
    out = write_spiked(
        bam, fasta, [LOSS], [GAIN], tmp_path / "spiked" / "out.cram", threads=2
    )
    assert out == tmp_path / "spiked" / "out.cram"
    assert Path(f"{out}.crai").is_file()
    assert (tmp_path / "gated.txt").read_text() == "dna:spikein\n"
    assert sorted(p.name for p in out.parent.iterdir()) == ["out.cram", "out.cram.crai"]
    before, after = _reads(bam, fasta), _reads(out, fasta)
    for (name, flag, start, seq), count in before.items():
        assert count == 1
        if name in drop:
            assert after[(name, flag, start, seq)] == 0
        elif name in duplicate:
            assert after[(name, flag, start, seq)] == 1
            assert after[(f"{name}:spike", flag, start, seq)] == 1
        else:
            assert after[(name, flag, start, seq)] == 1
    assert sum(after.values()) == len(before) - 2 * len(drop) + 2 * len(duplicate)
