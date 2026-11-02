"""Inherited chromosomally integrated HHV-6 (slice 3f, "HHV-6 band"): the viral depth over
each genome's unique region divided by the haploid human depth, classified against the
calibration's background and the band around one copy, and the `hhv6:ici` statement.

`fixtures/hg002/hg002.hhv6.coverage.txt` is `samtools coverage -r <contig>:<unique region>`
on HG002's unmapped pairs aligned to both genomes (Task 5): no read in either region."""

import json
import math
from pathlib import Path

import pytest
from fixtures.run_context import run_context

from genokeeper import hhv6
from genokeeper.config import Hhv6Thresholds
from genokeeper.hhv6 import Hhv6Call, classify, read_coverage
from genokeeper.sources import INFORMATION_NOTE
from genokeeper.sources.hhv6 import (
    BETWEEN,
    HHV6_TOOL,
    ICI,
    ICI_CONFIRM,
    Hhv6Source,
)
from genokeeper.store.db import Kind

HG002 = Path(__file__).parent / "fixtures" / "hg002" / "hg002.hhv6.coverage.txt"
A, B = "NC_001664.4", "NC_000898.1"
GRADING = Hhv6Thresholds()


def test_read_coverage_hg002() -> None:
    depths = read_coverage(HG002)
    assert depths == {A: 0.0, B: 0.0}
    # HG002 at its own depth (about 36x, so 18x haploid): ratio 0, below any background.
    call = classify(depths, 18.0, grading=GRADING, background=0.03)
    assert (call.ratio, call.kind, call.species) == (0.0, Kind.NON_FINDING, None)
    # A zero background (no viral reads in any calibration library): ratio 0 is still a
    # non-finding; any viral depth above it is unassessed (needs ddPCR).
    zero = classify(depths, 18.0, grading=GRADING, background=0.0)
    assert zero.kind == Kind.NON_FINDING
    some = classify({A: 0.18, B: 0.0}, 18.0, grading=GRADING, background=0.0)
    assert some.kind == Kind.UNASSESSED


@pytest.mark.parametrize(
    ("ratio", "kind"),
    [
        (0.04, Kind.NON_FINDING),
        (0.3, Kind.UNASSESSED),
        (0.5, Kind.FINDING),
        (1.0, Kind.FINDING),
        (1.5, Kind.FINDING),
        (1.6, Kind.UNASSESSED),
    ],
)
def test_classify_band_edges(ratio: float, kind: Kind) -> None:
    haploid = 20.0
    call = classify(
        {A: ratio * haploid / 4, B: ratio * haploid},
        haploid,
        grading=GRADING,
        background=0.05,
    )
    assert call.ratio == pytest.approx(ratio)
    assert call.kind == kind
    assert call.species == "HHV-6B"  # the higher contig
    swapped = classify(
        {A: ratio * haploid, B: 0.0}, haploid, grading=GRADING, background=0.05
    )
    assert (swapped.species, swapped.kind) == ("HHV-6A", kind)


@pytest.mark.parametrize("haploid", [math.nan, 0.0])
def test_nan_depth_unassessed(haploid: float, tmp_path: Path) -> None:
    call = classify({A: 10.0, B: 0.0}, haploid, grading=GRADING, background=0.05)
    assert call.kind == Kind.UNASSESSED
    (row,) = (
        Hhv6Source(call=call, not_run=None, withheld=None)
        .collect(run_context(tmp_path))
        .rows
    )
    assert row.statement["kind"] == "unassessed"
    assert row.statement["reason"] == hhv6.NO_DEPTH


def _statement(tmp_path: Path, source: Hhv6Source) -> dict[str, object]:
    collected = source.collect(run_context(tmp_path))
    assert collected.tool == HHV6_TOOL
    (row,) = collected.rows
    assert row.statement["statement_id"] == "hhv6:ici"
    assert (row.statement["topic"], row.statement["actionable"]) == ("wellness", False)
    assert row.statement["exclusion"] is False
    return row.statement


def test_hhv6_finding_is_information(tmp_path: Path) -> None:
    """A ratio in the band is a finding of information: the species named, what it means
    for blood tests and children, ddPCR as the confirming test, and `INFORMATION_NOTE`
    (the band's width is arithmetic, not measured)."""
    call = Hhv6Call(species="HHV-6A", ratio=0.95, kind=Kind.FINDING)
    statement = _statement(tmp_path, Hhv6Source(call=call, not_run=None, withheld=None))
    reason = str(statement["reason"])
    assert statement["kind"] == "finding"
    assert statement["validity"] == "moderate"
    assert ICI.format(species="HHV-6A") in reason
    assert ICI_CONFIRM in reason
    assert INFORMATION_NOTE in reason
    claim = json.loads(str(statement["claim"]))
    assert (claim["species"], claim["confirm"]) == ("HHV-6A", "ddPCR")
    between = Hhv6Call(species="HHV-6B", ratio=0.3, kind=Kind.UNASSESSED)
    statement = _statement(
        tmp_path, Hhv6Source(call=between, not_run=None, withheld=None)
    )
    assert (statement["kind"], statement["reason"]) == ("unassessed", BETWEEN)
    below = Hhv6Call(species=None, ratio=0.0, kind=Kind.NON_FINDING)
    statement = _statement(
        tmp_path, Hhv6Source(call=below, not_run=None, withheld=None)
    )
    assert statement["kind"] == "non_finding"
    for not_run, withheld in (("not calibrated", None), (None, "not covered: x")):
        statement = _statement(
            tmp_path, Hhv6Source(call=call, not_run=not_run, withheld=withheld)
        )
        assert (statement["kind"], statement["reason"]) == (
            "unassessed",
            not_run or withheld,
        )


def test_extract_and_align_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Pairs with either read unmapped come out of one gated samtools pipeline; bwa mem
    aligns them natively, gated `dna:hhv6`, into a sorted indexed BAM; each file is
    renamed into place only once its pipeline has finished."""
    ran: list[list[str]] = []

    def fake(args: list[str], **_: object) -> str:
        ran.append(args)
        script = args[2]
        if "samtools fastq" in script:
            for part in args[-2:]:
                Path(part).write_text("@r\nA\n+\nI\n")
        else:
            Path(args[-1]).write_text("bam")
            Path(args[-1] + ".bai").write_text("bai")
        return ""

    monkeypatch.setattr(hhv6, "run_tool", fake)
    fq1, fq2 = hhv6.extract_unmapped(
        tmp_path / "s.cram", tmp_path / "ref.fa", tmp_path / "out"
    )
    assert (fq1.name, fq2.name) == ("unmapped_1.fq", "unmapped_2.fq")
    assert fq1.is_file() and fq2.is_file()
    assert "flag.unmap || flag.munmap" in ran[0][2]
    assert "-F 0xD00" in ran[0][2]
    bam = hhv6.align_viral(
        fq1, fq2, tmp_path / "hhv6.fa", tmp_path / "out" / "hhv6.bam", threads=4
    )
    assert bam == tmp_path / "out" / "hhv6.bam"
    assert bam.is_file() and bam.with_name("hhv6.bam.bai").is_file()
    assert "bwa mem -t 4" in ran[1][2] and "samtools sort" in ran[1][2]
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        "hhv6.bam",
        "hhv6.bam.bai",
        "unmapped_1.fq",
        "unmapped_2.fq",
    ]
