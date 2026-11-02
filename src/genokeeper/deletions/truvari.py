"""Truvari 5.5.0 between the two callers and against HG002's truth (slice 3b).

`truvari bench` runs offline in its image (`containers.IMAGES`, built from docker/) behind the
core gate as `dna:truvari`, each input's folder mounted read-only and one output directory
read-write, with no reference (no `-f`), so no sequence is compared. Truvari writes a folder
of `tp-base`, `tp-comp`, `fn` and `fp` VCFs; the readers here take the calls back from it.
`write_calls` writes our own VCF of symbolic `<DEL>` records that carry their class in
`INFO/CLASS`, so nothing is joined back after Truvari (final re-review R4). The same runs
match another pair of callers over another type (slice 3d: `names`, `svtype`).
"""

from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

from genokeeper.containers import run_container
from genokeeper.deletions.callers import (
    caller_min_length,
    read_deletions,
    read_records,
    sort_index,
)
from genokeeper.deletions.screen import Deletion

LABEL = "dna:truvari"
# The two callers' union, matched: Dysgu's deletions as base and DELLY's as comparison, with
# no reference, so no sequence is compared (`--pctseq 0`); `--pick single` pairs each call
# once. The plan's "Decisions": `truvari collapse` would merge without saying which caller
# found what (Truvari wiki, v5.5.0, 2026-09-19).
BETWEEN_CALLERS: tuple[str, ...] = (
    "--pctseq", "0", "--pctsize", "0.5", "--pctovl", "0.5", "--passonly",
    "--pick", "single", "--sizemax", "-1",
)  # fmt: skip
# Against the truth: GIAB's own v5.0q command (its README, updated 2026-07-17; `--pick ac
# --passonly -r 2000 -C 5000` with `--includebed` the benchmark BED), plus `--pctseq 0
# --sizemin 35 --sizemax -1`, without `-f` and `--refine`: every comparison record is ours
# and symbolic, which refine skips. `--sizemin 35` (0.7 x the 50 bp floor) lets a true
# deletion a little under 50 bp match a call of 50 bp, so that call is not counted false.
AGAINST_TRUTH: tuple[str, ...] = (
    "--pctseq", "0", "--pick", "ac", "--passonly", "-r", "2000", "-C", "5000",
    "--sizemin", "35", "--sizemax", "-1",
)  # fmt: skip
# The symbolic ALTs `write_calls` declares, by `Deletion.svtype`.
_ALTS = {"DEL": "Deletion", "DUP": "Duplication", "INS": "Insertion"}
# Characters that would break an INFO value.
_NOT_IN_CLASS = frozenset(" \t;=,")


def _part(folder: Path, name: str) -> Path:
    return folder / f"{name}.vcf.gz"


def match_callers(
    base: Path,
    comp: Path,
    out_dir: Path,
    *,
    min_length: int,
    names: tuple[str, str] = ("dysgu", "delly"),
) -> Path:
    """`truvari bench` of the first of `names`' calls `base` (3b: Dysgu's deletions) against
    the second's `comp` (DELLY's), calls from `min_length` bp, a partner from the second
    caller's own floor (`--sizefilt`, `callers.caller_min_length`); returns the folder it
    writes, `out_dir/both`."""
    run_container(
        "truvari",
        [
            "truvari", "bench", "-b", f"/base/{base.name}", "-c", f"/comp/{comp.name}",
            "-o", "/out/both", *BETWEEN_CALLERS,
            "--sizemin", str(min_length),
            "--sizefilt", str(caller_min_length(names[1], min_length)),
        ],
        mounts=[
            (base.parent, "/base", True),
            (comp.parent, "/comp", True),
            (out_dir, "/out", False),
        ],
        label=LABEL,
    )  # fmt: skip
    return out_dir / "both"


def read_union(
    folder: Path, names: tuple[str, str] = ("dysgu", "delly"), svtype: str = "DEL"
) -> list[Deletion]:
    """The callers' union of `svtype` from `match_callers`'s folder: `tp-base`'s records of
    the first of `names` (3b: Dysgu), each with its partner's id, QUAL and `depth_ratio`
    (the `tp-comp` record of the same `MatchId`), then `fn`'s (the first only) and `fp`'s
    (the second only, with their own)."""
    base, comp = names
    partners = {
        info["MatchId"]: d
        for d, info in read_records(
            _part(folder, "tp-comp"), caller=comp, svtype=svtype
        )
    }
    both: list[Deletion] = []
    for d, info in read_records(_part(folder, "tp-base"), caller=base, svtype=svtype):
        partner = partners.get(info["MatchId"])
        if partner is None:
            raise ValueError(
                f"{folder}: no tp-comp record of MatchId {info['MatchId']}"
            )
        both.append(
            replace(
                d,
                callers={**d.callers, **partner.callers},
                depth_ratio=partner.depth_ratio,
            )
        )
    return [
        *both,
        *read_deletions(_part(folder, "fn"), caller=base, svtype=svtype),
        *read_deletions(_part(folder, "fp"), caller=comp, svtype=svtype),
    ]


def write_calls(calls: Sequence[tuple[Deletion, str]], fasta: Path, out: Path) -> Path:
    """Our own VCF `out` of `calls`, each a symbolic `<SVTYPE>` record (FILTER PASS, its GT)
    with `INFO/END`, `SVTYPE`, `SVLEN` (a deletion's negative; an insertion's `END` its POS)
    and `CLASS=<its key>`, the contigs from the FASTA's `.fai`; sorted, bgzipped and indexed
    natively."""
    for _, key in calls:
        if not key or _NOT_IN_CLASS & set(key):
            raise ValueError(f"CLASS {key!r} would break the INFO field")
    contigs = [
        line.split("\t")[:2]
        for line in Path(f"{fasta}.fai").read_text().splitlines()
        if line
    ]
    lines = [
        "##fileformat=VCFv4.2",
        *(f"##contig=<ID={name},length={length}>" for name, length in contigs),
        *(
            f'##ALT=<ID={svtype},Description="{description}">'
            for svtype, description in _ALTS.items()
            if any(d.svtype == svtype for d, _ in calls)
        ),
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##INFO=<ID=CLASS,Number=1,Type=String,Description="Call class">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    ]
    lines += [
        f"{d.chrom}\t{d.pos}\t{n}\t{d.ref}\t<{d.svtype}>\t.\tPASS\t"
        f"END={_end(d)};SVTYPE={d.svtype};SVLEN={_svlen(d)};CLASS={key}\tGT\t{d.gt}"
        for n, (d, key) in enumerate(calls)
    ]
    plain = out.with_name(f"{out.name}.unsorted.vcf")
    plain.write_text("\n".join(lines) + "\n")
    return sort_index(plain, out)


def _end(d: Deletion) -> int:
    return d.pos if d.svtype == "INS" else d.end


def _svlen(d: Deletion) -> int:
    return -d.length if d.svtype == "DEL" else d.length


def against_truth(calls: Path, truth: Path, bench: Path, out_dir: Path) -> Path:
    """`truvari bench` of the truth (base) against `calls` (comparison, `write_calls`'s),
    inside the benchmark BED `bench`; returns the folder it writes, `out_dir/truth`."""
    run_container(
        "truvari",
        [
            "truvari", "bench", "-b", f"/truth/{truth.name}", "-c", f"/calls/{calls.name}",
            "--includebed", f"/bench/{bench.name}", "-o", "/out/truth", *AGAINST_TRUTH,
        ],
        mounts=[
            (truth.parent, "/truth", True),
            (calls.parent, "/calls", True),
            (bench.parent, "/bench", True),
            (out_dir, "/out", False),
        ],
        label=LABEL,
    )  # fmt: skip
    return out_dir / "truth"


def judged(
    folder: Path, *, svtype: str = "DEL"
) -> tuple[list[tuple[str, int, bool]], list[tuple[int, bool]]]:
    """From `against_truth`'s folder: each comparison call's `CLASS`, length and whether the
    truth confirms it (`tp-comp` true, `fp` false), and each truth event's length and
    whether a call found it (`tp-base` true, `fn` false); the records of `svtype`."""
    calls = [
        (info["CLASS"], d.length, confirmed)
        for part, confirmed in (("tp-comp", True), ("fp", False))
        for d, info in read_records(_part(folder, part), caller="calls", svtype=svtype)
    ]
    truth = [
        (d.length, found)
        for part, found in (("tp-base", True), ("fn", False))
        for d in read_deletions(_part(folder, part), caller="truth", svtype=svtype)
    ]
    return calls, truth
