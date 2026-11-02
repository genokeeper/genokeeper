"""Slice 3d's three calibration tables built by `genokeeper calibrate` (`calibrate.build`'s
`build_cnv`, `build_insertions`, `build_paralogs`) on a scratch lab of four calibration
libraries (HG002's two, HG003, HG004), every container faked with what the tools wrote on
HG002 (`fixtures.genome_mini`), and 3b's judging loop moved into `judge_library`."""

import gzip
import hashlib
import json
import platform
import sys
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pytest
from fixtures.clingen_mini import write_dosage
from fixtures.genome_mini import (
    DOSAGE_DATE,
    LARGE_REGION,
    FakeDeletions,
    FakeSlice3d,
    lab_paths,
    write_clingen_regions,
    write_paralog_models,
    write_parascopy,
)
from test_calibrate_deletions import _gnomad, _hg002_lab, _strata

from genokeeper import insertions, paralogs, smn, tools
from genokeeper.calibrate import build, calibration_dir, tables
from genokeeper.calibrate.build import (
    ROWS,
    build_cnv,
    build_deletions,
    build_insertions,
    build_paralogs,
)
from genokeeper.calibrate.classes import open_class_table
from genokeeper.calibrate.records import (
    RECORDED,
    Record,
    library_dir,
    read_aligned,
    read_called,
    write_record,
)
from genokeeper.calibrate.slice3d import (
    INSERTION_BANDS,
    TRIO,
    cnv_versions,
    insertion_versions,
    open_paralog_table,
    paralog_versions,
)
from genokeeper.calibrate.tables import REPORT, TABLES, calibrate, folder
from genokeeper.checksum import sha256_file
from genokeeper.cnv import callers as cnv_callers
from genokeeper.config import grch38_fasta, load_grading
from genokeeper.deletions import callers as deletion_callers
from genokeeper.deletions import truvari as deletion_truvari
from genokeeper.deletions.exons import load_exons, targets
from genokeeper.deletions.files import sv_files, truth_files
from genokeeper.deletions.screen import Deletion
from genokeeper.evidence.clingen import DOSAGE_FILE
from genokeeper.masking import masked_bed
from genokeeper.pgx import cyp2d6
from genokeeper.pipeline import RunRefused

SLICE3D_TABLES = ("cnv", "insertions", "paralogs")
# The scratch lab's GRCh38 contigs: chr20 holds every HG002 fixture call, chrX and chrY tell
# a male's shards from a female's.
FAI = {"chr20": 64_444_167, "chrX": 156_040_895, "chrY": 57_227_415}
# v5.0q's benchmark covers chr20's first half: a gain past it is judged by the trio.
BENCH_END = 30_000_000
# A container call: the library whose CRAM it read, its arguments, its name, its `/out`.
type Call = tuple[str, list[str], str, str]


def _truth(path: Path, extra: Sequence[str] = ()) -> None:
    """v5.0q's SV truth: two deletions of 1 kb or more, a short one, two insertions, and
    `extra`."""
    records = [
        "chr20\t5000000\td1\tN\t<DEL>\t.\tPASS\tEND=5002000;SVTYPE=DEL\tGT\t0/1",
        "chr20\t6000000\td2\tN\t<DEL>\t.\tPASS\tEND=6003000;SVTYPE=DEL\tGT\t0/1",
        "chr20\t6500000\td3\tN\t<DEL>\t.\tPASS\tEND=6500200;SVTYPE=DEL\tGT\t0/1",
        "chr20\t7000000\ti1\tN\t<INS>\t.\tPASS\tEND=7000000;SVTYPE=INS;SVLEN=300\tGT\t0/1",
        "chr20\t7500000\ti2\tN\t<INS>\t.\tPASS\tEND=7500000;SVTYPE=INS;SVLEN=90\tGT\t0/1",
        *extra,
    ]
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress((header + "\n".join(records) + "\n").encode()))


class Recorded:
    """Every container call, with the calibration library whose CRAM it read (its `/cram`
    mount's library folder; empty for a call that reads none), before `run` takes it."""

    def __init__(self, lab: Path, run: Any) -> None:
        self.run = run
        self.calls: list[Call] = []
        self._folders = {library_dir(lab, key): key for key in RECORDED}

    def library(self, mounts: Sequence[tuple[Path, str, bool]]) -> str:
        for source, target, _ in mounts:
            if target == "/cram":
                for parent in (source, *source.parents):
                    if parent in self._folders:
                        return self._folders[parent]
        return ""

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        out = next((s for s, t, _ in mounts if t == "/out"), Path())
        self.calls.append((self.library(mounts), list(args), name, out.name))
        self.run(name, args, mounts=mounts, platform=platform, label=label)


class Lab:
    """The scratch lab, its fakes, and what the spike-in writer was asked."""

    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        paths = lab_paths(root)
        self.lab = lab = paths.lab
        _strata(lab)
        _gnomad(sv_files(lab).gnomad, [("away", 100, 3000, "DEL", "PASS", 5000, 0.5)])
        fasta = self.fasta = grch38_fasta(lab)
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(">chr20\nN\n")
        Path(f"{fasta}.fai").write_text(
            "".join(f"{c}\t{n}\t7\t60\t61\n" for c, n in FAI.items())
        )
        masked_bed(lab).write_text("")
        truth, bench = truth_files(lab)
        _truth(truth)
        bench.write_text(f"chr20\t0\t{BENCH_END}\n")
        write_parascopy(lab)
        write_paralog_models(lab)
        write_clingen_regions(lab)
        write_dosage(
            lab / "clingen-dosage" / DOSAGE_DATE / DOSAGE_FILE, [("BRCA1", "3")]
        )
        for key in TRIO:
            self.record(key)
        monkeypatch.setattr(tools, "GATE", root / "no-gate")
        self.deletions = FakeDeletions()
        self.slice3d = FakeSlice3d()
        self.recorded = Recorded(lab, self.slice3d.over(self.deletions))
        for module in (deletion_callers, deletion_truvari, cnv_callers, smn, cyp2d6):
            monkeypatch.setattr(module, "run_container", self.recorded)
        for module in (insertions, paralogs):
            monkeypatch.setattr(module, "run_container", self.recorded)
            monkeypatch.setattr(module, "run_tool", self.slice3d.tool)
        self.spiked: list[tuple[str, int, int, int]] = []

        def write_spiked(
            cram: Path,
            fasta: Path,
            losses: Sequence[Deletion],
            gains: Sequence[Deletion],
            out: Path,
            *,
            threads: int = 1,
        ) -> Path:
            self.spiked.append((cram.stem, len(losses), len(gains), threads))
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(b"")
            return out

        monkeypatch.setattr(build, "write_spiked", write_spiked)

    def record(self, key: str) -> None:
        folder = library_dir(self.lab, key)
        folder.mkdir(parents=True, exist_ok=True)
        cram = folder / f"{key}.cram"
        cram.write_bytes(b"")
        record = Record(
            cram=str(cram),
            key=key * 8,
            metrics={"mean_depth": 30.0},
            gvcf=None,
            reference_sha256=sha256_file(self.fasta),
            accepted=True,
            snp_f1=None,
            replaced=None,
        )
        write_record(self.lab, key, record)

    def fasta_sha(self) -> str:
        return sha256_file(self.fasta)


@pytest.fixture
def lab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Lab:
    return Lab(tmp_path, monkeypatch)


def _built(lab: Path) -> set[str]:
    return {name for name in TABLES if folder(lab, name).is_dir()}


def test_slice3d_tables_registered() -> None:
    """3d's three tables are rows of 3c's registry with no table built before them, and 3b's
    deletion table never waits on one (C6)."""
    for name, builder in zip(
        SLICE3D_TABLES, (build_cnv, build_insertions, build_paralogs), strict=True
    ):
        assert TABLES[name].needs == ()
        assert ROWS[name][1] is builder
    assert not set(TABLES["deletions"].needs) & set(SLICE3D_TABLES)


def test_build_only_paralogs(lab: Lab, capsys: pytest.CaptureFixture[str]) -> None:
    """`genokeeper calibrate --only paralogs` builds that table alone, opened under the run's
    own versions, its report printed; no other table's step runs."""
    assert calibrate(lab.lab, only="paralogs", report=True, threads=2) == 0
    assert _built(lab.lab) == {"paralogs"}
    versions = paralog_versions(lab.lab, lab.fasta_sha(), load_grading())
    table, reason = open_paralog_table(lab.lab, versions)
    assert reason is None and table is not None
    assert set(table.measured) == set(TRIO)
    out = capsys.readouterr().out
    assert "paralogs table" in out
    # Only the tools a verdict reads: Parascopy on every library, Gauchian and KILDA on
    # HG002's two (`paralog_measure`), Aldy nowhere.
    ran = {(key, name) for key, _, name, _ in lab.recorded.calls}
    assert {name for _, name in ran} == {"parascopy", "gauchian", "kilda"}
    assert {k for k, n in ran if n in ("gauchian", "kilda")} == {"pcr", "pcrfree"}
    # Each library's CRAM read (Parascopy mounts its own elsewhere: empty).
    assert {key for key, _ in ran} - {""} == set(TRIO)


@pytest.mark.parametrize(
    ("builder", "missing"),
    [(build_paralogs, "hg004"), (build_cnv, "hg003"), (build_insertions, "pcr")],
)
def test_build_refuses_missing_sample(lab: Lab, builder: Any, missing: str) -> None:
    """A library not aligned refuses the table with 3c's own reason, naming its align
    command, before any step runs; nothing is written."""
    (library_dir(lab.lab, missing) / "aligned.json").unlink()
    with pytest.raises(RunRefused) as e:
        builder(lab.lab, 2)
    # 3c's own refusal: its wording and command (`align` or `realign`) are 3c's.
    with pytest.raises(RunRefused) as own:
        read_aligned(lab.lab, missing)
    assert str(e.value) == str(own.value)
    assert lab.recorded.calls == []
    assert _built(lab.lab) == set()


def test_step_failure_refuses_that_table_only(
    lab: Lab, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Gauchian failing on the first library refuses the paralogs table, naming the library
    and the tool; `genokeeper calibrate` writes the other two."""
    monkeypatch.setattr(
        tables, "TABLES", {name: TABLES[name] for name in SLICE3D_TABLES}
    )
    lab.slice3d.fail.add("gauchian")
    assert calibrate(lab.lab, only=None, report=False, threads=2) == 2
    assert _built(lab.lab) == {"cnv", "insertions"}
    err = capsys.readouterr().err
    assert "genokeeper calibrate: paralogs: HG002 pcrfree: gauchian failed" in err


# The table's shares come from `math.lgamma`, whose last digits differ by libm: the same table
# hashes to another digest on Linux x86_64 and another again on Linux arm64.
@pytest.mark.skipif(
    (sys.platform, platform.machine()) != ("darwin", "arm64"),
    reason="the digest pins floats as macOS arm64's libm computes them",
)
def test_build_deletions_uses_judge_library(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """3b's table, its judging moved into `judge_library` (once per library), is the table
    3b's own loop measured on its fixture library (its digest taken at bb93de1)."""
    paths = _hg002_lab(tmp_path, monkeypatch, (2000, 200, 40))
    judged: list[Path] = []
    judge = build.judge_library

    def spy(work: Path, *args: Any, **kw: Any) -> Any:
        judged.append(work)
        return judge(work, *args, **kw)

    monkeypatch.setattr(build, "judge_library", spy)
    table = build_deletions(paths, threads=2)
    assert [p.parent.name for p in judged] == ["pcrfree", "pcr"]
    document = json.dumps(
        {
            "bands": {k: asdict(v) for k, v in table.bands.items()},
            "lk": {k: asdict(v) for k, v in table.likelihoods.items()},
            "measured": table.measured,
        },
        sort_keys=True,
    )
    assert hashlib.sha256(document.encode()).hexdigest() == (
        "18c847d75386b6820713fc3d747e5fed2dbef4dc05e0cbea2b794671e15eb080"
    )
    assert {band: terms.calls for band, terms in table.bands.items()} == {
        "long": 68,
        "short": 186,
    }


def test_hg004_runs_female(lab: Lab) -> None:
    """HG004 (female) runs its depth caller without `--haploid-chr` and every caller without
    a chrY shard; HG003 (male) with both."""
    build_cnv(lab.lab, 2)

    def shards(library: str, *, depth: bool) -> set[str]:
        return {
            out
            for key, args, name, out in lab.recorded.calls
            if key == library
            and name == "delly"
            and args[:1] != ["cp"]
            and (args[:2] == ["delly", "cnv"]) == depth
        }

    def haploid(library: str) -> list[bool]:
        return [
            "--haploid-chr" in args
            for key, args, name, _ in lab.recorded.calls
            if key == library and args[:2] == ["delly", "cnv"]
        ]

    assert shards("hg004", depth=True) == {"chr20", "chrX"}
    assert shards("hg004", depth=False) == {"chr20", "chrX"}
    assert not any(haploid("hg004"))
    assert shards("hg003", depth=True) == {"chr20", "chrX", "chrY"}
    assert haploid("hg003") and all(haploid("hg003"))


class _Large(FakeSlice3d):
    """`FakeSlice3d` with a 600 kb loss both depth callers make over the ClinGen region on
    the libraries in `on` only."""

    def __init__(self, on: set[str], recorded: Recorded) -> None:
        super().__init__()
        self.on, self.recorded = on, recorded

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        chrom, start, _ = LARGE_REGION
        end = start + 600_000
        here = self.recorded.library(mounts) in self.on
        self.dellycnv = (
            [
                (
                    f"{chrom}\t{start + 1}\tBIG\tN\t<CNV>\t60\tPASS"
                    f"\tIMPRECISE;SVTYPE=CNV;END={end}\tGT:CN:FT\t./.:1:PASS"
                )
            ]
            if here
            else []
        )
        self.cnvpytor = (
            [f"deletion\t{chrom}:{start + 1}-{end}\t600000\t0.5\t0.0"] if here else []
        )
        super().__call__(name, args, mounts=mounts, platform=platform, label=label)


def test_large_withheld_recorded(lab: Lab, monkeypatch: pytest.MonkeyPatch) -> None:
    """A large loss on HG002 (one library) and on HG003 counts once per genome: two of the
    three genomes call one, so `loss_large` is withheld in the `cnv` table's
    `extra["large"]`, and the report says so."""
    large = _Large({"pcrfree", "hg003"}, lab.recorded)
    lab.recorded.run = large.over(lab.deletions)
    for module in (insertions, paralogs):
        monkeypatch.setattr(module, "run_tool", large.tool)
    build_cnv(lab.lab, 2)
    versions = cnv_versions(lab.lab, lab.fasta_sha(), load_grading())
    table, reason = open_class_table(lab.lab, "cnv", versions)
    assert reason is None and table is not None
    pooled = table.extra["large"]["calls"]
    assert {sample: counts["loss_large"] for sample, counts in pooled.items()} == {
        "HG002": 1,
        "HG003": 1,
        "HG004": 0,
    }
    reason = "not covered: large loss calls on 2 of 3 calibration genomes"
    assert table.extra["large"]["withheld"] == {"loss_large": reason}
    assert table.extra["cnvpytor"]["decision"] in ("kept", "dropped")
    report = (folder(lab.lab, "cnv") / REPORT).read_text()
    assert f"loss_large: {reason}" in report
    assert "large calls on HG003: loss_1-5kb" not in report
    assert "CNVpytor " in report


def test_pcr_library_read_as_cram_only(lab: Lab) -> None:
    """HG002 PCR's record has no calls (3c never calls `pcr`): every 3d table is built from
    its CRAM alone, each under the run's own versions; the spike-ins are written on
    HG002's two libraries with the builder's threads."""
    with pytest.raises(RunRefused, match="not called"):
        read_called(lab.lab, "pcr")
    for builder in (build_cnv, build_insertions, build_paralogs):
        builder(lab.lab, 2)
    assert _built(lab.lab) == set(SLICE3D_TABLES)
    grading, sha = load_grading(), lab.fasta_sha()
    cnv, why = open_class_table(lab.lab, "cnv", cnv_versions(lab.lab, sha, grading))
    assert why is None and cnv is not None and set(cnv.measured) == {"pcrfree", "pcr"}
    # Losses judged against the truth's deletions (the fake confirms each), the one gain
    # past the benchmark by the trio: both parents call it too, so false (C18); the
    # spike-ins are the recall's truth.
    for library in ("pcrfree", "pcr"):
        loss = cnv.measured[library]["loss_1-5kb"]
        assert loss["calls"] > 0 and loss["false"] == 0
        gain = cnv.measured[library]["gain_1-5kb"]
        assert (gain["calls"], gain["false"]) == (1, 1)
        spiked = sum(cnv.measured[library][band]["truth"] for band in cnv.bands)
        assert spiked == 120
    ins, why = open_class_table(
        lab.lab, "insertions", insertion_versions(lab.lab, sha, grading)
    )
    assert why is None and ins is not None
    assert set(ins.bands) == set(INSERTION_BANDS) and set(ins.measured) == {
        "pcrfree",
        "pcr",
    }
    short = ins.measured["pcrfree"]["short"]
    assert short["calls"] > 0 and (short["truth"], short["found"]) == (2, 2)
    assert [(cram, threads) for cram, _, _, threads in lab.spiked] == [
        ("pcrfree", 2),
        ("pcr", 2),
    ]
    assert not (
        calibration_dir(lab.lab) / "hg002" / "pcr" / "spikein" / "spiked.cram"
    ).exists()


def test_scramble_meis_judged_by_position(lab: Lab) -> None:
    """Scramble's lone mobile elements (band `mei`) inside the benchmark are judged by
    position: one within `max_distance` of a truth insertion is true and finds it, one near
    none is false; one outside the benchmark is not counted."""
    starts = sorted(
        x.start
        for ts in targets(load_exons(), 20).values()
        for x in ts
        if x.chrom == "chrX"
    )
    near, alone, outside = starts[0] + 10, starts[50] + 10, starts[100] + 10
    truth, bench = truth_files(lab.lab)
    info = f"END={near + 200};SVTYPE=INS;SVLEN=300"
    _truth(truth, [f"chrX\t{near + 200}\ti3\tN\t<INS>\t.\tPASS\t{info}\tGT\t0/1"])
    bench.write_text(f"chr20\t0\t{BENCH_END}\nchrX\t0\t{starts[99]}\n")
    lab.slice3d.meis += [
        f"chrX:{at}\talu\tPlus\t7\t130.0\t100\t100\tGCC\tleft\t2\t64"
        + "\tNone Found" * 4
        + "\tNA"
        for at in (near, alone, outside)
    ]
    build_insertions(lab.lab, 2)
    versions = insertion_versions(lab.lab, lab.fasta_sha(), load_grading())
    table, reason = open_class_table(lab.lab, "insertions", versions)
    assert reason is None and table is not None
    for library in ("pcrfree", "pcr"):
        mei = table.measured[library]["mei"]
        assert (mei["calls"], mei["false"], mei["truth"], mei["found"]) == (2, 1, 1, 1)
        assert mei["classes"] == ["mei:scramble:norepeat:unmasked:nojunction"]
    assert table.bands["mei"].calls == 4
