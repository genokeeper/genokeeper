"""Slice 3d's spike-ins (the plan's Decisions, "Spike-ins"): copy-number losses and gains
written into a copy of an HG002 library's CRAM, so the depth callers' recall per band can be
measured where HG002 carries no such event.

`pick_spikes` chooses the events: per direction `n_per` each of 1, 2 and 3 consecutive
exons of a screen gene, none touching the mask (3b's and 3c's masked intervals, as the
caller joins them) or a region a locus-specific caller owns, all `SPACING` apart, each from
the middle of the intron before its first exon to the middle of the intron after its last,
capped at `FLANK_CAP` a side and at least `MIN_SPIKE` long. `write_spiked` spikes by read
name: every name with a read overlapping an event and `zlib.crc32(name) % 2 == 0`
(`spike_names`) is dropped for a loss and written twice for a gain, both mates wherever the
mate lies (C22), the copies named `<name>:spike`. It runs as one unit behind the core gate
(`dna:spikein`): this module's command line, which selects the named pairs with native
samtools, rewrites them with pysam, merges them back and indexes the result.
"""

import argparse
import json
import random
import sys
import tempfile
import zlib
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from pathlib import Path

import pysam

from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import Deletion, _Index
from genokeeper.regions import Interval
from genokeeper.tools import gated, run_tool

SPIKEIN_LABEL = "dna:spikein"
EXON_COUNTS = (1, 2, 3)
FLANK_CAP = 20_000
MIN_SPIKE = 1_000
SPACING = 100_000
N_PER = 20
SEED = 3
SPIKE_CALLER = "spikein"
COPY_SUFFIX = ":spike"


def _windows(
    gene: str, ts: Sequence[Target], k: int, blocked: _Index[None]
) -> list[tuple[str, int, int, str]]:
    """Each run of `k` consecutive targets with a target on either side, as its event span
    and its `<gene>:<exon numbers>` label, unless it is short or touches `blocked`."""
    found: list[tuple[str, int, int, str]] = []
    for i in range(1, len(ts) - k):
        before, first, last, after = ts[i - 1], ts[i], ts[i + k - 1], ts[i + k]
        start = max((before.end + first.start) // 2, first.start - FLANK_CAP)
        end = min((last.end + after.start) // 2, last.end + FLANK_CAP)
        if end - start < MIN_SPIKE or blocked.hits(first.chrom, start, end):
            continue
        numbers = sorted({n for t in ts[i : i + k] for n in t.exons})
        found.append((first.chrom, start, end, f"{gene}:{','.join(map(str, numbers))}"))
    return found


def pick_spikes(
    targets: Mapping[str, Sequence[Target]],
    mask: Sequence[Interval],
    owned: Sequence[Interval],
    *,
    n_per: int = N_PER,
    seed: int = SEED,
) -> tuple[list[Deletion], list[Deletion]]:
    """The losses and the gains to spike (see the module docstring), drawn by `seed`; a
    shorter list when the genes hold too few runs `SPACING` apart. Each event's caller is
    `spikein`, naming its gene and exon numbers."""
    blocked = _Index((c, s, e, None) for c, s, e in [*mask, *owned])
    candidates: dict[int, list[tuple[str, int, int, str]]] = {}
    for k in EXON_COUNTS:
        candidates[k] = [
            w
            for gene in sorted(targets)
            for w in _windows(
                gene, sorted(targets[gene], key=lambda t: t.start), k, blocked
            )
        ]
    rng = random.Random(seed)
    for k in EXON_COUNTS:
        rng.shuffle(candidates[k])
    chosen: list[Interval] = []
    losses: list[Deletion] = []
    gains: list[Deletion] = []
    for k in EXON_COUNTS:
        for svtype, out in (("DEL", losses), ("DUP", gains)):
            picked = 0
            for chrom, start, end, label in candidates[k]:
                if picked == n_per:
                    break
                if any(
                    c == chrom and max(start - e, s - end) < SPACING
                    for c, s, e in chosen
                ):
                    continue
                chosen.append((chrom, start, end))
                called = {SPIKE_CALLER: label}
                out.append(
                    Deletion(chrom, start, end, "N", "0/1", called, None, svtype=svtype)
                )
                picked += 1
    return losses, gains


def _spiked(name: str) -> bool:
    return zlib.crc32(name.encode()) % 2 == 0


def spike_names(
    cram: Path, fasta: Path, losses: Sequence[Deletion], gains: Sequence[Deletion]
) -> tuple[set[str], set[str]]:
    """The read names to drop (a read over a loss) and to write twice (over a gain), each
    with `crc32 % 2 == 0`; a name over both is dropped."""
    with pysam.AlignmentFile(str(cram), reference_filename=str(fasta)) as f:

        def names(events: Sequence[Deletion]) -> set[str]:
            return {
                r.query_name
                for d in events
                for r in f.fetch(d.chrom, d.pos, d.end)
                if r.query_name is not None and _spiked(r.query_name)
            }

        drop = names(losses)
        return drop, names(gains) - drop


def _rewrite(
    selected: Path, fasta: Path, drop: set[str], duplicate: set[str], out: Path
) -> None:
    """`selected`'s reads into the BAM `out`: a dropped name's skipped, any other's kept,
    and a duplicated one's written again as `<name>:spike`."""
    with (
        pysam.AlignmentFile(str(selected), reference_filename=str(fasta)) as f,
        pysam.AlignmentFile(str(out), "wb", template=f) as written,
    ):
        for read in f.fetch(until_eof=True):
            name = read.query_name or ""
            if name in drop:
                continue
            written.write(read)
            if name not in duplicate:
                continue
            copy = pysam.AlignedSegment.from_dict(read.to_dict(), f.header)
            copy.query_name = f"{name}{COPY_SUFFIX}"
            written.write(copy)


def _spike(
    cram: Path,
    fasta: Path,
    losses: Sequence[Deletion],
    gains: Sequence[Deletion],
    out: Path,
    work: Path,
    threads: int = 1,
) -> None:
    """The unit `write_spiked` gates: select the named pairs, rewrite, merge, index, each
    samtools run on `threads`."""
    drop, duplicate = spike_names(cram, fasta, losses, gains)
    names, selected = work / "names.txt", work / "selected.cram"
    rest, rewritten = work / "rest.cram", work / "spiked.bam"
    names.write_text("".join(f"{n}\n" for n in sorted(drop | duplicate)))
    run_tool(
        [
            "samtools", "view", "-@", str(threads), "-N", str(names), "-U", str(rest),
            "-O", "cram",
            "-T", str(fasta), "-o", str(selected), str(cram),
        ]
    )  # fmt: skip
    _rewrite(selected, fasta, drop, duplicate, rewritten)
    try:
        run_tool(
            [
                "samtools", "merge", "-@", str(threads), "-f", "-O", "cram",
                "--reference", str(fasta),
                "-o", str(out), str(rest), str(rewritten),
            ]
        )  # fmt: skip
        run_tool(["samtools", "index", "-@", str(threads), str(out)])
    except BaseException:
        out.unlink(missing_ok=True)
        Path(f"{out}.crai").unlink(missing_ok=True)
        raise


def write_spiked(
    cram: Path,
    fasta: Path,
    losses: Sequence[Deletion],
    gains: Sequence[Deletion],
    out: Path,
    *,
    threads: int = 1,
) -> Path:
    """`cram` with `losses` and `gains` spiked in (the module docstring), as the indexed
    CRAM `out`, one unit behind the core gate, its samtools runs on `threads`; its
    temporaries, in a hidden folder beside `out`, are removed. Returns `out`."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent, prefix=".spikein-") as folder:
        spikes = Path(folder) / "spikes.json"
        spikes.write_text(
            json.dumps(
                {
                    "losses": [asdict(d) for d in losses],
                    "gains": [asdict(d) for d in gains],
                }
            )
        )
        command = [sys.executable, "-m", "genokeeper.calibrate.spikein"]
        run_tool(
            gated(
                SPIKEIN_LABEL,
                [*command, str(cram), str(fasta), str(spikes), str(out), folder]
                + ["--threads", str(threads)],
            )
        )
    return out


def main(argv: Sequence[str] | None = None) -> int:
    """`write_spiked`'s gated unit: CRAM, FASTA, the spikes' JSON, the output and a work
    folder."""
    parser = argparse.ArgumentParser(prog="python -m genokeeper.calibrate.spikein")
    for name in ("cram", "fasta", "spikes", "out", "work"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    spikes = json.loads(args.spikes.read_text())
    losses = [Deletion(**d) for d in spikes["losses"]]
    gains = [Deletion(**d) for d in spikes["gains"]]
    _spike(args.cram, args.fasta, losses, gains, args.out, args.work, args.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
