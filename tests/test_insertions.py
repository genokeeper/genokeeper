"""Insertions of 50 bp or more (`genokeeper.insertions`): INSurVeyor's and Scramble's HG002
outputs (tests/fixtures/hg002, Task 6) read, the Dysgu+INSurVeyor union read back from
Truvari's folder, Scramble's mobile elements joined to it, and the shown insertions stated
through 3b's screen and source. Containers and native tools are stubs that record their
arguments and write the output.

The counts are the tests' own reading of the fixtures (`_count`), never the code under
test."""

import gzip
import json
import shutil
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

import pytest
from fixtures.run_context import run_context

from genokeeper import insertions
from genokeeper.calibrate.deletions import class_of
from genokeeper.config import load_grading
from genokeeper.deletions import callers, truvari
from genokeeper.deletions.callers import read_records
from genokeeper.deletions.exons import Target
from genokeeper.deletions.screen import (
    CallClass,
    Deletion,
    Screen,
    Shown,
    classify,
    key,
    screen_all,
)
from genokeeper.genelists import DPYD_NOTE, screen_overrides
from genokeeper.insertions import (
    INSERTION_CALLERS,
    INSERTION_WORDING,
    INSERTIONS_TOOL,
    INSURVEYOR_LABEL,
    MEI_SCOPE,
    SCRAMBLE_LABEL,
    add_scramble,
    insertion_union,
    mei_family,
    read_scramble,
    run_insurveyor,
    run_scramble,
)
from genokeeper.sources import Collected
from genokeeper.sources.deletions import CONSEQUENCE_WARNING, DeletionSource
from genokeeper.store.db import Kind

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
INSURVEYOR = FIXTURES / "hg002.insurveyor.vcf.txt"
SCRAMBLE = FIXTURES / "hg002.scramble.meis.txt"
PARTS = ("tp-base", "tp-comp", "fn", "fp")
CONFIRM = "needs confirmation by long reads or PCR across the insertion"
# BRCA2's exon 2 (MANE v1.5) widened by 3b's 20 bp margin, and DPYD's exon 22.
BRCA2 = {"BRCA2": [Target("chr13", 32316401, 32316547, (2,))]}
DPYD = {"DPYD": [Target("chr1", 97082309, 97082490, (22,))]}
INS_CALLERS = {"dysgu": "75357 QUAL 7", "insurveyor": "T_INS_3831 QUAL ."}


def _count(path: Path, *, pass_only: bool = False) -> int:
    rows = [line.split("\t") for line in path.open() if not line.startswith("#")]
    return sum(1 for r in rows if not pass_only or r[6] == "PASS")


def _insertion(
    pos: int,
    *,
    chrom: str = "chr13",
    svlen: int | None = 300,
    found: Mapping[str, str] = INS_CALLERS,
) -> Deletion:
    return Deletion(
        chrom, pos, pos + 1, "A", "0/1", found, None, svtype="INS", svlen=svlen
    )


def _screen(targets: Mapping[str, list[Target]] = BRCA2) -> Screen:
    return Screen(
        targets=targets,
        mask={},
        repeats={},
        population=[],
        thresholds=load_grading().insertions,
        callers=INSERTION_CALLERS,
    )


def _collect(
    tmp_path: Path,
    shown: list[Shown],
    *,
    rate: float = 0.5,
    overrides: Mapping[str, tuple[str, str]] | None = None,
) -> Collected:
    source = DeletionSource(
        shown=shown,
        called=len(shown),
        dropped={},
        not_run=None,
        p_wrong=lambda c: (rate, rate / 2, min(rate * 2, 1.0)),
        haploinsufficient=frozenset({"BRCA2"}),
        recall={},
        wording=INSERTION_WORDING,
        placement_overrides=overrides or {},
        screen_claim={"not_covered": [MEI_SCOPE]},
    )
    return source.collect(run_context(tmp_path))


def _gene_rows(collected: Collected) -> list[dict[str, object]]:
    return [
        r.statement
        for r in collected.rows
        if r.statement["statement_id"] != f"{INSERTIONS_TOOL}:screen"
    ]


def test_read_insurveyor_hg002() -> None:
    found = read_records(INSURVEYOR, caller="insurveyor", svtype="INS")
    assert len(found) == _count(INSURVEYOR, pass_only=True) == 92
    for d, _ in found:
        assert d.svtype == "INS"
        assert d.svlen is not None and d.svlen >= 50
        assert d.end == d.pos + 1
        assert set(d.callers) == {"insurveyor"}


def _bcftools(seen: list[list[str]]) -> object:
    """Stands for native bcftools: records each call and writes its `-o` file."""

    def tool(args: list[str], **_: object) -> str:
        seen.append(list(args))
        if "-o" in args:
            Path(args[args.index("-o") + 1]).write_bytes(b"")
        return ""

    return tool


def test_insertion_union_hg002(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dysgu = tmp_path / "in" / "dysgu.ins.vcf.gz"
    raw = tmp_path / "in" / "insurveyor.vcf.gz"
    out = tmp_path / "out"
    for path in (dysgu, raw):
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(b"")
    out.mkdir()
    tools: list[list[str]] = []
    benches: list[tuple[list[str], list[tuple[Path, str, bool]]]] = []

    def bench(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, label) == ("truvari", "dna:truvari")
        benches.append((list(args), list(mounts)))
        folder = out / "both"
        folder.mkdir()
        for part in PARTS:
            text = (FIXTURES / f"hg002.ins.{part}.vcf.txt").read_bytes()
            (folder / f"{part}.vcf.gz").write_bytes(gzip.compress(text))

    monkeypatch.setattr(callers, "run_tool", _bcftools(tools))
    monkeypatch.setattr(truvari, "run_container", bench)
    union = insertion_union(dysgu, raw, out, min_length=50)
    expected = sum(
        _count(FIXTURES / f"hg002.ins.{p}.vcf.txt") for p in ("tp-base", "fn", "fp")
    )
    assert len(union) == expected == 34 + 137 + 49
    assert all(set(d.callers) <= set(INSERTION_CALLERS) for d in union)
    assert all(d.svtype == "INS" and d.end == d.pos + 1 for d in union)
    keys = [frozenset(d.callers) for d in union]
    assert keys.count(frozenset(INSERTION_CALLERS)) == 34
    # INSurVeyor's raw calls kept to its PASS insertions of 50 bp or more, then matched.
    view = tools[0]
    assert view[:2] == ["bcftools", "view"]
    assert view[-1] == str(raw)
    expression = view[view.index("-i") + 1]
    assert 'INFO/SVTYPE="INS"' in expression and ">=50" in expression
    ((args, mounts),) = benches
    assert args[args.index("-b") + 1] == f"/base/{dysgu.name}"
    kept = Path(view[view.index("-o") + 1]).parent
    assert (kept, "/comp", True) in mounts
    assert args[args.index("--sizefilt") + 1] == "50"


def test_insertion_statement_wording(tmp_path: Path) -> None:
    d = _insertion(32316450)
    ((shown,), dropped) = screen_all([d], _screen())
    assert sum(dropped.values()) == 0
    collected = _collect(tmp_path, [shown])
    (statement,) = _gene_rows(collected)
    assert statement["statement_id"] == "insertions:BRCA2:chr13:32316450-32316451"
    assert statement["vid"] == "GRCh38:chr13:32316450:A:<INS>"
    assert statement["kind"] == str(Kind.FINDING)
    assert statement["validity"] == "limited"
    assert statement["source"] == INSERTIONS_TOOL == "insertions"
    assert statement["source_level"] == (
        "Dysgu 1.9.0 and INSurVeyor 1.1.3, matched by Truvari 5.5.0;"
        " Scramble 1.0.2 for mobile elements"
    )
    assert str(statement["reason"]).startswith(f"possible insertion, {CONFIRM}")
    claim = json.loads(str(statement["claim"]))
    assert (claim["type"], claim["length"]) == ("insertion", 300)
    likely = _gene_rows(_collect(tmp_path, [shown], rate=0.01))[0]
    assert str(likely["reason"]).startswith(f"likely present, {CONFIRM}")


def test_insertion_screen_statement_names_mei_scope(tmp_path: Path) -> None:
    (screen,) = [
        r.statement
        for r in _collect(tmp_path, []).rows
        if r.statement["statement_id"] == "insertions:screen"
    ]
    claim = json.loads(str(screen["claim"]))
    assert claim["type"] == "insertion_screen"
    assert claim["not_covered"] == [MEI_SCOPE]
    assert MEI_SCOPE == (
        "mobile elements labelled by Scramble on the target exons only; MELT and xTea"
        " are not open to a private user"
    )
    assert str(screen["reason"]).startswith("screened: 0 insertions called, 0 shown;")


def test_read_scramble_hg002() -> None:
    meis = read_scramble(SCRAMBLE)
    rows = [line for line in SCRAMBLE.read_text().splitlines()[1:] if line.strip()]
    assert len(meis) == len(rows) == 3
    for d in meis:
        assert (d.svtype, d.svlen, d.lower_bound) == ("INS", None, True)
        assert d.end == d.pos + 1
        assert mei_family(d) in {"Alu", "L1", "SVA"}
    # By hand: the first row, chr22:50623930, alu, Plus, 5 clipped reads.
    first = meis[0]
    assert (first.chrom, first.pos) == ("chr22", 50623930)
    assert first.callers == {"scramble": "Alu Plus 5"}


def test_read_scramble_header_only(tmp_path: Path) -> None:
    path = tmp_path / "meis.txt"
    path.write_text(SCRAMBLE.read_text().splitlines()[0] + "\n")
    assert read_scramble(path) == []
    other = SCRAMBLE.read_text().replace("\talu\t", "\tHERV\t", 1)
    path.write_text(other)
    with pytest.raises(ValueError, match="HERV"):
        read_scramble(path)


def test_add_scramble(tmp_path: Path) -> None:
    s = _screen()
    call = _insertion(32316450)
    before = key(classify(call, s))
    near = replace(
        _insertion(32316650, svlen=None, found={"scramble": "Alu Plus 6"}),
        lower_bound=True,
    )
    far = replace(near, pos=32318450, end=32318451)
    joined = add_scramble([call], [near, far], max_distance=500)
    assert len(joined) == 2
    labelled, alone = joined
    assert labelled.callers == {**INS_CALLERS, "scramble": "Alu Plus 6"}
    assert (labelled.pos, labelled.svlen) == (call.pos, call.svlen)
    assert key(classify(labelled, s)) == before
    assert mei_family(labelled) == "Alu"
    assert alone == far
    c = classify(alone, s)
    assert (c.band, c.evidence) == ("mei", "scramble")
    # Through the screen and the source: the labelled call reads its family; the MEI alone
    # is not cut short (its length is unknown) and lies off every target here.
    shown, dropped = screen_all(joined, s)
    assert [x.deletion for x in shown] == [labelled]
    assert dropped["no listed exon"] == 1
    (statement,) = _gene_rows(_collect(tmp_path, shown))
    assert str(statement["reason"]).startswith(f"possible Alu insertion, {CONFIRM}")
    inside = replace(far, pos=32316460, end=32316461)
    ((only,), _) = screen_all([inside], s)
    assert (only.call_class.band, only.call_class.evidence) == ("mei", "scramble")
    (mei,) = _gene_rows(_collect(tmp_path, [only]))
    assert str(mei["reason"]).startswith(f"possible Alu insertion, {CONFIRM}")
    assert json.loads(str(mei["claim"]))["length"] is None
    call_row = _collect(tmp_path, [only]).rows[1].call
    assert call_row is not None
    assert call_row["warnings"] == [
        "found by one caller (Scramble)",
        CONSEQUENCE_WARNING,
    ]


def test_add_scramble_nearest_labels() -> None:
    """Two MEIs near one call: the nearer labels it; an MEI on another chromosome is its
    own event."""
    call = _insertion(10_000)
    mei = _insertion(10_300, svlen=None, found={"scramble": "L1 Minus 4"})
    nearer = _insertion(9_900, svlen=None, found={"scramble": "SVA Plus 9"})
    other = replace(mei, chrom="chr14")
    labelled, alone = add_scramble([call], [mei, nearer, other], max_distance=500)
    assert labelled.callers["scramble"] == "SVA Plus 9"
    assert alone == other


def _cram_and_fasta(tmp_path: Path) -> tuple[Path, Path, Path]:
    cram, fasta, out = (
        tmp_path / "in" / "x.cram",
        tmp_path / "ref" / "g.fa",
        tmp_path / "out",
    )
    for path in (cram, fasta):
        path.parent.mkdir()
        path.write_bytes(b"")
    out.mkdir()
    return cram, fasta, out


type Seen = list[tuple[str, list[str], list[tuple[Path, str, bool]], str | None]]


def test_run_scramble_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path)
    regions = tmp_path / "exons.bed"
    regions.write_text("chr13\t32315007\t32316167\n")
    seen: Seen = []
    tools: list[list[str]] = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), list(mounts), label))
        shutil.copy(SCRAMBLE, out / "scramble_MEIs.txt")

    monkeypatch.setattr(insertions, "run_container", fake)
    monkeypatch.setattr(insertions, "run_tool", _bcftools(tools))
    meis = run_scramble(cram, fasta, out, regions=regions, threads=3)
    assert meis == out / "meis.txt" and meis.read_text() == SCRAMBLE.read_text()
    view, index = tools
    assert view[:2] == ["samtools", "view"]
    assert "-M" in view and "-b" in view
    assert view[view.index("-L") + 1] == str(regions)
    assert view[view.index("-T") + 1] == str(fasta)
    assert view[view.index("-o") + 1] == str(out / "exons.bam")
    assert view[view.index("-@") + 1] == "3"
    assert view[-1] == str(cram)
    assert index == ["samtools", "index", str(out / "exons.bam")]
    ((name, args, mounts, label),) = seen
    assert (name, label) == ("scramble", SCRAMBLE_LABEL)
    assert SCRAMBLE_LABEL == "dna:scramble"
    assert args[:2] == ["bash", "-c"]
    steps = [s.strip() for s in args[2].split("&&")]
    assert "cluster_identifier /out/exons.bam > /out/clusters.txt" in steps
    assert "cd /opt/scramble/cluster_analysis/bin" in steps
    (rscript,) = [s for s in steps if s.startswith("Rscript")]
    assert "--eval-meis" in rscript and "--eval-dels" not in rscript
    assert "--ref /ref/g.fa" in rscript
    assert "--cluster-file /out/clusters.txt" in rscript
    assert "--out-name /out/scramble" in rscript
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }


def test_run_insurveyor_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cram, fasta, out = _cram_and_fasta(tmp_path)
    seen: Seen = []
    tools: list[list[str]] = []

    def fake(
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        seen.append((name, list(args), list(mounts), label))
        (out / "work").mkdir()
        (out / "work" / "out.pass.vcf.gz").write_bytes(b"")

    monkeypatch.setattr(insertions, "run_container", fake)
    monkeypatch.setattr(callers, "run_tool", _bcftools(tools))
    result = run_insurveyor(cram, fasta, out, threads=4)
    assert result == out / "insurveyor.vcf.gz"
    ((name, args, mounts, label),) = seen
    assert (name, label) == ("insurveyor", INSURVEYOR_LABEL)
    assert INSURVEYOR_LABEL == "dna:insurveyor"
    assert args == [
        "insurveyor.py", "--threads", "4", "/cram/x.cram", "/out/work", "/ref/g.fa",
    ]  # fmt: skip
    assert (out, "/out", False) in mounts
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/out", False),
    }
    sort, index = tools
    assert sort[:2] == ["bcftools", "sort"]
    assert sort[-1] == str(out / "work" / "out.pass.vcf.gz")
    assert sort[sort.index("-o") + 1] == str(result)
    assert index[-1] == str(result)


def test_dpyd_insertion_drugs(tmp_path: Path) -> None:
    d = _insertion(97082400, chrom="chr1")
    ((shown,), _) = screen_all([d], _screen(DPYD))
    (gated,) = _gene_rows(_collect(tmp_path, [shown]))
    assert gated["kind"] == str(Kind.NON_FINDING)
    (statement,) = _gene_rows(_collect(tmp_path, [shown], overrides=screen_overrides()))
    assert (statement["kind"], statement["topic"]) == ("finding", "drugs")
    assert DPYD_NOTE in str(statement["reason"])
    assert str(statement["reason"]).startswith(f"possible insertion, {CONFIRM}")


def test_insertion_class_keys_roundtrip() -> None:
    s = _screen()
    for band in ("short", "long", "mei"):
        for evidence in ("both", *INSERTION_CALLERS, "scramble"):
            for repeat in (True, False):
                for masked in (True, False):
                    c = CallClass(band, evidence, repeat, masked)
                    assert class_of(key(c)) == c
    # Every class `classify` gives an insertion is one of those.
    short = classify(_insertion(32316450, svlen=60), s)
    long = classify(_insertion(32316450, svlen=6000), s)
    assert (short.band, long.band) == ("short", "long")
    assert short.junction is False
