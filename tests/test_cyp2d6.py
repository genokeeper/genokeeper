"""CYP2D6 from BCyrius and Aldy: each caller run offline in its image and read, with the
CPIC activity score each states, and the two agreed on the score.

`fixtures/aldy_cyp2d6_na10860.aldy` is Aldy 4.8.3's output on its own test sample
(`aldy/tests/resources/NA10860_hg38.bam` in its sdist, the public GeT-RM sample NA10860)
from a run of the built image, offline, on 2026-10-20: the header, each solution's line and
its first two rows. `fixtures/bcyrius_nocall.tsv` is BCyrius 2.2.0's output on that BAM,
which holds only the CYP2D6 region, so its depth normalisation finds nothing: a no-call.
"""

from pathlib import Path
from typing import Any

import pytest

from genokeeper.pgx import cyp2d6
from genokeeper.pgx.cyp2d6 import agree, read_aldy, read_bcyrius, run_aldy, run_bcyrius
from genokeeper.pgx.outside import CallerCall, cpic_score, outside_lines

FIXTURES = Path(__file__).parent / "fixtures"


def _call(
    caller: str,
    diplotype: str | None,
    score: str | None = None,
    *,
    abstained: bool = False,
) -> CallerCall:
    return CallerCall(
        caller=caller,
        diplotype=diplotype,
        phenotype=None,
        activity_score=None if score is None else cpic_score(score),
        abstained=abstained,
    )


def test_agree_on_activity_not_star_string() -> None:
    """`*1/*4` and `*2/*4` are different star alleles of one score: the callers agree on
    the score each states, which is what PharmCAT is given, without either diplotype."""
    a, b = _call("BCyrius", "*1/*4", "1.0"), _call("Aldy", "*2/*4", "1")
    agreed = agree(a, b, "CYP2D6")
    assert (agreed.gene, agreed.agreed, agreed.unresolved_reason) == (
        "CYP2D6",
        True,
        None,
    )
    assert (agreed.diplotype, agreed.activity_score) == (None, "1.0")
    assert agreed.callers == {
        "BCyrius": "*1/*4 (activity score 1.0)",
        "Aldy": "*2/*4 (activity score 1.0)",
    }
    assert outside_lines([agreed]) == "CYP2D6\t\t\t1.0\n"
    # One diplotype however each caller writes it (order, sub-alleles, Aldy's markers):
    # it goes to PharmCAT with the score.
    same = agree(
        _call("BCyrius", "*4/*1", "1.0"),
        _call("Aldy", "*1.018/*4.021X", "1.0"),
        "CYP2D6",
    )
    assert (same.diplotype, same.activity_score) == ("*1/*4", "1.0")
    assert outside_lines([same]) == "CYP2D6\t*1/*4\t\t1.0\n"
    # Three and four copies of a function allele are CPIC's one bin, ≥3.0.
    bin_ = agree(
        _call("BCyrius", "*1x3/*4", "3.0"), _call("Aldy", "*1x4/*4", "4.0"), "CYP2D6"
    )
    assert (bin_.agreed, bin_.diplotype, bin_.activity_score) == (True, None, "≥3.0")
    assert outside_lines([bin_]) == "CYP2D6\t\t\t≥3.0\n"


def test_disagree_unresolved_shows_both() -> None:
    """Different scores, or a score one caller lacks, leave the gene unresolved with both
    calls shown; no outside line is written for it."""
    a, b = _call("BCyrius", "*1/*4", "1.0"), _call("Aldy", "*1/*41", "1.5")
    split = agree(a, b, "CYP2D6")
    assert (split.agreed, split.diplotype, split.activity_score) == (False, None, None)
    assert split.unresolved_reason == (
        "BCyrius and Aldy differ on the activity score (1.0 against 1.5)"
    )
    assert split.callers == {
        "BCyrius": "*1/*4 (activity score 1.0)",
        "Aldy": "*1/*41 (activity score 1.5)",
    }
    assert outside_lines([split]) == ""
    nocall = agree(_call("BCyrius", None), b, "CYP2D6")
    assert nocall.unresolved_reason == "BCyrius made no call"
    assert nocall.callers == {
        "BCyrius": "no call",
        "Aldy": "*1/*41 (activity score 1.5)",
    }
    both_unscored = agree(
        _call("BCyrius", "*1/*22", "n/a"), _call("Aldy", "*1/*41", "n/a"), "CYP2D6"
    )
    assert both_unscored.unresolved_reason == (
        "no activity score for BCyrius's *1/*22; no activity score for Aldy's *1/*41"
    )
    assert outside_lines([both_unscored]) == ""


def test_one_unscored_caller_leaves_the_others_score_unagreed() -> None:
    """A caller that names a diplotype but no score (Aldy's `cpic=indeterminate`) leaves the
    other caller's score to stand alone, not agreed, so the gene call is at most medium; the
    unscored call stays shown."""
    lone = agree(
        _call("BCyrius", "*6/*68+*68+*68+*4", "0.0"),
        _call("Aldy", "*68/*68+*4C", "n/a", abstained=True),
        "CYP2D6",
    )
    assert (lone.agreed, lone.diplotype, lone.activity_score) == (False, None, "0.0")
    assert lone.unresolved_reason is None
    assert lone.note == (
        "one caller's score: Aldy stated no activity score for *68/*68+*4C,"
        " so BCyrius's 0.0 was not checked"
    )
    assert lone.callers["Aldy"] == "*68/*68+*4C (no activity score)"
    assert outside_lines([lone]) == "CYP2D6\t\t\t0.0\n"


def test_several_scores_are_not_an_abstention() -> None:
    """A caller whose solutions state different scores has no one score but did not
    abstain: the gene stays unresolved, as for two different scores."""
    several = agree(
        _call("BCyrius", "*1/*1", "2.0"), _call("Aldy", "*1/*4 or *1/*1"), "CYP2D6"
    )
    assert (several.agreed, several.activity_score, several.note) == (False, None, None)
    assert several.unresolved_reason == "no activity score for Aldy's *1/*4 or *1/*1"


def test_readers_mark_an_abstention(tmp_path: Path) -> None:
    """Aldy's solutions all without `cpic_score`, or BCyrius's `n/a`, abstain; Aldy's
    solutions with different scores do not."""
    head = "#Sample\tGene\tSolutionID\tMajor\tMinor\n"
    lone = tmp_path / "lone.aldy"
    lone.write_text(
        head + "#Solution 1: *(4.004), *(68.001); gene=CYP2D6; cpic=indeterminate\n"
        "s\tCYP2D6\t1\t*4/*68\t\n",
        encoding="utf-8",
    )
    split = tmp_path / "split.aldy"
    split.write_text(
        head + "#Solution 1: *(1), *(4); cpic=IM; cpic_score=1.0\n"
        "#Solution 2: *(1), *(1); cpic=NM; cpic_score=2.0\n"
        "s\tCYP2D6\t1\t*1/*4\t\ns\tCYP2D6\t2\t*1/*1\t\n",
        encoding="utf-8",
    )
    assert read_aldy(lone).abstained
    assert not read_aldy(split).abstained
    tsv = tmp_path / "b.tsv"
    tsv.write_text(
        "Sample\tGenotype\tFilter\tActivity score\tPredicted phenotype\n"
        "s\t*1/*22\tPASS\tn/a\tIndeterminate\n",
        encoding="utf-8",
    )
    assert read_bcyrius(tsv).abstained


def test_three_copies_agreed_are_binned() -> None:
    """Three copies written either way are one diplotype, which PharmCAT is given in
    CPIC's `x≥3` bin, with the bin's score."""
    agreed = agree(
        _call("BCyrius", "*1x3/*4", "3.0"), _call("Aldy", "*1+*1+*1/*4", "3"), "CYP2D6"
    )
    assert (agreed.diplotype, agreed.activity_score) == ("*1x3/*4", "≥3.0")
    assert outside_lines([agreed]) == "CYP2D6\t*1x≥3/*4\t\t≥3.0\n"


def _fake_container(calls: list[dict[str, Any]], output: Path | None):
    """A `run_container` that records its call and leaves `output` in the mounted output
    directory, under the name the tool writes."""

    def run(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        calls.append(
            {
                "name": name,
                "args": args,
                "mounts": mounts,
                "platform": platform,
                "label": label,
            }
        )
        (host,) = [h for h, _, read_only in mounts if not read_only]
        if output is not None:
            written = Path(args[args.index("-o") + 1]).name if "-o" in args else None
            target = written or f"{args[args.index('--id') + 1]}.tsv"
            (host / target).write_text(output.read_text())

    return run


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    cram = tmp_path / "work" / "sample.cram"
    fasta = tmp_path / "ref" / "genome.fa"
    for path in (cram, fasta):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x")
    return cram, fasta


def test_aldy_command_has_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Aldy reads the CRAM with the reference it was written against, as WGS on hg38,
    offline in its image; its solutions are read back."""
    calls: list[dict[str, Any]] = []
    fixture = FIXTURES / "aldy_cyp2d6_na10860.aldy"
    monkeypatch.setattr(cyp2d6, "run_container", _fake_container(calls, fixture))
    cram, fasta = _inputs(tmp_path)
    out = tmp_path / "out"
    call = run_aldy(cram, fasta, out, "CYP2D6")
    (run,) = calls
    assert (run["name"], run["platform"], run["label"]) == ("aldy", None, "dna:pgx")
    assert run["args"] == [
        "aldy",
        "genotype",
        "-p",
        "wgs",
        "-g",
        "CYP2D6",
        "--genome",
        "hg38",
        "-r",
        "/aldy-ref/genome.fa",
        "-o",
        "/aldy-out/CYP2D6.aldy",
        "/aldy-in/sample.cram",
    ]
    assert sorted(run["mounts"]) == sorted(
        [
            (cram.parent, "/aldy-in", True),
            (fasta.parent, "/aldy-ref", True),
            (out, "/aldy-out", False),
        ]
    )
    assert call == read_aldy(out / "CYP2D6.aldy")
    assert call.caller == "Aldy"
    assert call.diplotype == "*1/*4+*4.021X or *1/*4.021+*4NX"
    assert (call.phenotype, call.activity_score) == ("intermediate", "1.0")


def test_bcyrius_run_and_read(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[dict[str, Any]] = []
    monkeypatch.setattr(
        cyp2d6, "run_container", _fake_container(calls, FIXTURES / "bcyrius_nocall.tsv")
    )
    cram, fasta = _inputs(tmp_path)
    out = tmp_path / "out"
    call = run_bcyrius(cram, fasta, out)
    (run,) = calls
    assert (run["name"], run["platform"], run["label"]) == ("bcyrius", None, "dna:pgx")
    assert run["args"] == [
        "bcyrius",
        "--input",
        "/bcyrius-in/sample.cram",
        "--reference",
        "/bcyrius-ref/genome.fa",
        "--genome",
        "chr38",
        "--outDir",
        "/bcyrius-out",
        "--id",
        "bcyrius",
    ]
    # BCyrius's no-call: genotype `None`.
    assert (call.caller, call.diplotype, call.activity_score) == ("BCyrius", None, None)
    called = tmp_path / "called.tsv"
    lines = (FIXTURES / "bcyrius_nocall.tsv").read_text().splitlines()
    # Sample, Genotype, Filter, Call method, Preferred diplotype, Other candidates,
    # Activity score, Predicted phenotype.
    row = "S\t*2x2/*4\tPASS\tRead_backed_phase\t\t\t2.0\tNormal Metabolizer"
    called.write_text(f"{lines[0]}\n{row}\n")
    read = read_bcyrius(called)
    assert (read.diplotype, read.phenotype, read.activity_score) == (
        "*2x2/*4",
        "Normal Metabolizer",
        "2.0",
    )


def test_a_filtered_bcyrius_call_is_no_agreeing_caller(tmp_path: Path) -> None:
    """BCyrius flags a call it doubts (`LowQ_high_CN`): Aldy's equal score alone does not
    make two agreeing callers, and the gene stays unresolved with the filter named."""
    tsv = tmp_path / "b.tsv"
    tsv.write_text(
        "Sample\tGenotype\tFilter\tActivity score\tPredicted phenotype\n"
        "s\t*1/*1x3\tLowQ_high_CN\t3.0\tUltrarapid Metabolizer\n",
        encoding="utf-8",
    )
    called = agree(read_bcyrius(tsv), _call("Aldy", "*1/*1x3", "3"), "CYP2D6")
    assert (called.agreed, called.activity_score) == (False, None)
    assert called.unresolved_reason == "BCyrius filtered its *1/*1x3 as LowQ_high_CN"
    assert outside_lines([called]) == ""
