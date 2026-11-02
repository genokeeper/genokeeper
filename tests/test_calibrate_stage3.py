"""The stage-3 calibration table (`genokeeper.calibrate.stage3`): what HG002's two libraries
withhold or cap, HG002's repeat counts from GIAB's v5.0q phased truth, and the table written
whole and refused when missing or made with other versions. The tools' results are their real
output on HG002 (tests/fixtures/hg002)."""

import json
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pysam
import pytest

from genokeeper.calibrate import build, repeat_rules, stage3, tables
from genokeeper.calibrate.records import RECORDED, Record, library_dir, write_record
from genokeeper.calibrate.repeat_rules import RFC1_WITHHELD
from genokeeper.calibrate.stage3 import (
    FALSE_CALL,
    FORMAT,
    RULES,
    SMN_KEY,
    SMN_WITHHELD,
    TABLE,
    YLEAF_KEY,
    YLEAF_WITHHELD,
    Stage3Table,
    hg002_truth,
    measure,
    open_stage3,
    stage3_dir,
    stage3_versions,
    truth_repeats,
    write_table,
)
from genokeeper.calibrate.tables import NOT_CALIBRATED
from genokeeper.checksum import sha256_file
from genokeeper.config import Paths, grch38_fasta, load_grading
from genokeeper.pipeline import RunRefused, Stage3Results
from genokeeper.repeats.expansionhunter import read_calls, region
from genokeeper.repeats.loci import RepeatCall
from genokeeper.repeats.rfc1 import RFC1_LOCUS, read_rfc1
from genokeeper.smn import ParascopyCall, read_parascopy, read_smn_caller
from genokeeper.store.db import Level, Topic, Validity
from genokeeper.yhaplogroup import read_yleaf

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
LOCI = load_grading().repeats.loci
YLEAF_OUT = FIXTURES / "hg002.yleaf.out.txt"
EH_JSON = FIXTURES / "hg002.eh3e.json"


def _hg002() -> Stage3Results:
    calls = read_calls(FIXTURES / "hg002.eh3e.vcf.txt", EH_JSON)
    return Stage3Results(
        repeats=calls,
        rfc1=read_rfc1(EH_JSON, calls),
        images={},
        caller=read_smn_caller(FIXTURES / "hg002.smncaller.json"),
        parascopy=read_parascopy(FIXTURES / "hg002.parascopy.paralog.bed.txt"),
        yleaf=read_yleaf(FIXTURES / "hg002.yleaf.hg_prediction.txt"),
        yleaf_out=YLEAF_OUT,
        flank=1.0,
        failed={},
    )


def _sized(results: Stage3Results, gene: str, *counts: int) -> Stage3Results:
    """`results` with `gene`'s alleles at `counts`, each interval the count alone."""
    assert results.repeats is not None
    variant = LOCI[gene].variant_id
    call = results.repeats[variant]
    alleles = tuple(
        replace(a, repeats=n, ci=(n, n))
        for a, n in zip(call.alleles, counts, strict=True)
    )
    calls = dict(results.repeats) | {variant: replace(call, alleles=alleles)}
    return replace(results, repeats=calls)


def _truth(results: Stage3Results) -> dict[str, tuple[int, ...] | None]:
    """Each called locus's truth as the results call it: no size error anywhere (the
    fixture is slice 3e's run, on every locus)."""
    assert results.repeats is not None
    calls: dict[str, RepeatCall] = results.repeats
    return {
        gene: tuple(sorted(a.repeats for a in calls[locus.variant_id].alleles))
        for gene, locus in LOCI.items()
        if locus.variant_id in calls
    }


def test_measure_withholds_what_the_run_would_not_grade_normal() -> None:
    hg002 = _hg002()
    alone, _, _ = measure({"pcrfree": hg002}, _truth(hg002), LOCI)
    # ATXN2 23/27 (CI 26-66) and CACNA1A 12/12 (CI 12-27): normal estimates whose
    # intervals reach past the normal range, which the run grades as watches.
    assert alone["ATXN2"] == alone["CACNA1A"] == FALSE_CALL
    kept = {"HTT", "DMPK", "C9orf72", "ATXN1", "ATXN3", "ATXN7", "ATN1", "FXN", "AR"}
    assert not kept & set(alone)
    results = {"pcrfree": hg002, "pcr": _sized(hg002, "HTT", 17, 31)}
    withheld, _, measured = measure(results, _truth(hg002), LOCI)
    assert withheld["HTT"] == FALSE_CALL
    assert measured["pcr"]["repeats"]["HTT"]["alleles"] == [17, 31]


def test_measure_caps_size_error() -> None:
    hg002 = _hg002()
    truth = _truth(hg002) | {"HTT": (17, 24)}
    off = _sized(hg002, "HTT", 20, 24)
    _, caps, measured = measure({"pcrfree": hg002, "pcr": off}, truth, LOCI)
    assert caps["HTT"] == Level.LOW
    assert measured["pcr"]["repeats"]["HTT"]["size_error"] == 3
    within = _sized(hg002, "HTT", 19, 24)
    _, kept, _ = measure({"pcrfree": hg002, "pcr": within}, truth, LOCI)
    assert "HTT" not in kept
    # Outside the benchmark regions nothing is compared.
    _, unknown, _ = measure({"pcr": off}, truth | {"HTT": None}, LOCI)
    assert "HTT" not in unknown


def test_measure_caps_a_one_allele_truth_on_its_one_allele() -> None:
    hg002 = _hg002()
    truth = _truth(hg002) | {"AR": (27,)}
    _, caps, _ = measure({"pcrfree": hg002}, truth, LOCI)
    assert "AR" not in caps
    _, capped, _ = measure({"pcrfree": hg002}, truth | {"AR": (23,)}, LOCI)
    assert capped["AR"] == Level.LOW


def test_measure_smn_requires_two_two() -> None:
    hg002 = _hg002()
    assert hg002.caller is not None
    withheld, _, _ = measure({"pcrfree": hg002, "pcr": hg002}, _truth(hg002), LOCI)
    assert SMN_KEY not in withheld
    for other in (
        replace(hg002, parascopy=ParascopyCall(1, 2, 156.0, True)),
        replace(hg002, caller=replace(hg002.caller, smn2=3)),
        replace(hg002, parascopy=None),
    ):
        found, _, _ = measure({"pcrfree": hg002, "pcr": other}, _truth(hg002), LOCI)
        assert SMN_KEY in found


def test_measure_yleaf_requires_l816(tmp_path: Path) -> None:
    hg002 = _hg002()
    withheld, _, measured = measure({"pcrfree": hg002}, _truth(hg002), LOCI)
    assert YLEAF_KEY not in withheld
    assert measured["pcrfree"]["yleaf"]["l816"] is True
    ancestral = tmp_path / "sample.out"
    ancestral.write_text(
        "".join(
            line.replace("\tT\tD\t", "\tA\tA\t") if "\tL816\t" in line else line
            for line in YLEAF_OUT.read_text().splitlines(keepends=True)
        )
    )
    for other in (replace(hg002, yleaf_out=ancestral), replace(hg002, yleaf_out=None)):
        found, _, _ = measure({"pcrfree": hg002, "pcr": other}, _truth(hg002), LOCI)
        assert YLEAF_KEY in found


def test_measure_yleaf_requires_one_haplogroup_on_both_libraries() -> None:
    hg002 = _hg002()
    assert hg002.yleaf is not None
    same, _, _ = measure({"pcrfree": hg002, "pcr": hg002}, _truth(hg002), LOCI)
    assert YLEAF_KEY not in same
    # Both derive L816, but Yleaf predicts another branch below it on one library.
    other = replace(hg002, yleaf=replace(hg002.yleaf, name="J-L816"))
    found, _, _ = measure({"pcrfree": hg002, "pcr": other}, _truth(hg002), LOCI)
    assert YLEAF_KEY in found


def _truth_vcf(tmp_path: Path, records: list[str]) -> Path:
    """A VCF shaped as GIAB's HG002 v5.0q small-variant truth (`GT:AD`, phased, a male's
    chrX haploid), bgzipped and indexed."""
    plain = tmp_path / "truth.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        "##contig=<ID=chr4,length=190214555>\n"
        "##contig=<ID=chrX,length=156040895>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + "".join(f"{r}\n" for r in records)
    )
    return Path(pysam.tabix_index(str(plain), preset="vcf", force=True))


def test_truth_repeats_phased(tmp_path: Path) -> None:
    htt = region(LOCI["HTT"])
    assert htt == ("chr4", 3074877, 3074933)  # 19 CAG in GRCh38
    vcf = _truth_vcf(
        tmp_path,
        [
            "chr4\t3074876\t.\tCCAG\tC\t30\t.\t.\tGT:AD\t1|0:0,1",
            "chr4\t3080000\t.\tCAG\tC\t30\t.\t.\tGT:AD\t0|1:0,1",
            "chrX\t67545316\t.\tT\tTGCAGCAGCAGCA\t30\t.\t.\tGT:AD\t1:0,1",
            # HG002's FXN insertions as GIAB writes them: left-normalised into the
            # poly-A flank, anchored two bases before the region. The first is three GAA
            # units; the second also lengthens the poly-A, so it gives no GAA truth.
            "chr9\t69037284\t.\tA\tAAAGAAGAAG,AAAAAAGAAGAAG\t30\t.\t.\tGT:AD\t2|1:0,1,1",
            "chr9\t69100284\t.\tA\tAAAGAAGAAG\t30\t.\t.\tGT:AD\t1|1:0,1",
        ],
    )
    bench = tmp_path / "bench.bed"
    bench.write_text(
        "chr4\t3059966\t3227137\nchr9\t69000000\t69200000\nchrX\t67415899\t67762860\n"
    )
    fxn = region(LOCI["FXN"])
    assert fxn == ("chr9", 69037287, 69037304)  # 6 GAA in GRCh38
    assert truth_repeats(vcf, bench, region=fxn, motif="GAA") is None
    # The same three-unit insertion alone, two bases before a region of the same shape.
    shifted = ("chr9", fxn[1] + 63000, fxn[2] + 63000)
    assert truth_repeats(vcf, bench, region=shifted, motif="GAA") == (9, 9)
    assert truth_repeats(vcf, bench, region=htt, motif="CAG") == (18, 19)
    ar = region(LOCI["AR"])
    assert ar == ("chrX", 67545317, 67545385)  # 23 GCA in GRCh38: 22 CAG and the CAA
    assert truth_repeats(vcf, bench, region=ar, motif="CAG", copies=1) == (27,)
    # A haploid record where two haplotypes are expected says nothing.
    assert truth_repeats(vcf, bench, region=ar, motif="CAG") is None
    outside = tmp_path / "outside.bed"
    outside.write_text("chr4\t3074900\t3227137\n")
    assert truth_repeats(vcf, outside, region=htt, motif="CAG") is None


def _table(
    versions: dict[str, str], measured: dict[str, Any] | None = None
) -> Stage3Table:
    return Stage3Table(
        versions=versions,
        withheld={"AR": FALSE_CALL},
        caps={"HTT": Level.LOW},
        flank={"pcrfree": 1.1, "pcr": 0.6},
        measured={} if measured is None else measured,
    )


VERSIONS = {
    "expansionhunter": "genomelab/expansionhunter:5.0.0-a",
    "yleaf": "genomelab/yleaf:4.1.6-b",
    "catalog": "c" * 64,
}


def test_open_stage3_missing_and_stale(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    table, reason = open_stage3(lab, VERSIONS)
    assert table is None and reason == NOT_CALIBRATED
    assert NOT_CALIBRATED.endswith("run `genokeeper calibrate`")
    _put(lab, _table(VERSIONS))
    table, reason = open_stage3(lab, VERSIONS)
    assert reason is None and table == _table(VERSIONS)
    other = VERSIONS | {"expansionhunter": "genomelab/expansionhunter:5.0.0-z"}
    table, reason = open_stage3(lab, other)
    assert table is None and reason is not None
    assert "expansionhunter" in reason and "yleaf" not in reason
    assert reason.endswith("run `genokeeper calibrate`")


def _put(lab: Path, table: Stage3Table) -> None:
    tables.write_table(lab, "stage3", lambda s: write_table(s, table), table.versions)


def test_write_table_atomic(tmp_path: Path) -> None:
    folder = stage3_dir(tmp_path / "lab")
    broken = _table(VERSIONS, {"pcrfree": object()})
    with pytest.raises(TypeError):
        _put(tmp_path / "lab", broken)
    assert not folder.exists() and list(folder.parent.iterdir()) == []
    _put(tmp_path / "lab", _table(VERSIONS))
    with pytest.raises(TypeError):
        _put(tmp_path / "lab", broken)
    assert open_stage3(tmp_path / "lab", VERSIONS) == (_table(VERSIONS), None)
    assert [p.name for p in folder.parent.iterdir()] == [folder.name]


def test_versions_follow_the_grading(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(stage3, "parascopy_md5s", lambda _lab: "md5s")
    grading = load_grading()
    base = stage3_versions(tmp_path, "f" * 64, grading)
    assert json.loads(base["exclusion"]) == {"min_spanning": 3}
    assert json.loads(base["repeats"])["loci"]["CNBP"]["normal_max"] == 26
    spanning = replace(grading, exclusion=replace(grading.exclusion, min_spanning=4))
    moved = stage3_versions(tmp_path, "f" * 64, spanning)
    assert {k for k in base if base[k] != moved[k]} == {"exclusion"}
    loci = dict(grading.repeats.loci)
    loci["CNBP"] = replace(loci["CNBP"], normal_max=27)
    row = replace(grading, repeats=replace(grading.repeats, loci=loci))
    moved = stage3_versions(tmp_path, "f" * 64, row)
    assert {k for k in base if base[k] != moved[k]} == {"repeats"}
    # A wording field no measurement reads leaves the table current.
    for name, value in (
        ("note", "another note"),
        ("topic", Topic.RESEARCH),
        ("evidence", "Strong"),
        ("validity", Validity.LIMITED),
    ):
        loci = dict(grading.repeats.loci)
        loci["CNBP"] = replace(loci["CNBP"], **{name: value})
        worded = replace(grading, repeats=replace(grading.repeats, loci=loci))
        assert stage3_versions(tmp_path, "f" * 64, worded) == base, name
    # The trio checks' tolerance is in the versions too.
    monkeypatch.setattr(repeat_rules, "MENDEL_TOLERANCE", 3)
    moved = stage3_versions(tmp_path, "f" * 64, grading)
    assert {k for k in base if base[k] != moved[k]} == {"rules"}


def test_stage3_table_of_an_earlier_format_not_used(tmp_path: Path) -> None:
    lab = tmp_path / "lab"
    table = _table(VERSIONS)
    _put(lab, table)
    document = json.loads((stage3_dir(lab) / TABLE).read_text())
    assert document["format"] == FORMAT == 3
    assert "homozygous_ok" not in document
    assert open_stage3(lab, VERSIONS) == (table, None)
    # A format-2 table (with the false-normal trials' `homozygous_ok`) is not used.
    document["format"] = 2
    (stage3_dir(lab) / TABLE).write_text(json.dumps(document))
    found, reason = open_stage3(lab, VERSIONS)
    assert found is None and reason is not None and "format 2" in reason


# `build_stage3` on faked steps: the trio's records, the truth files and a FASTA on disk.


def _build_lab(tmp_path: Path, keys: tuple[str, ...]) -> Paths:
    """A lab with GIAB's HG002 truth files (empty), a FASTA, and a current record of each
    of `keys`, its CRAM on disk."""
    lab = tmp_path / "lab"
    fasta = grch38_fasta(lab)
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\nACGT\n")
    vcf, bench = hg002_truth(lab)
    vcf.parent.mkdir(parents=True)
    for path in (vcf, Path(f"{vcf}.tbi"), bench):
        path.write_bytes(b"")
    for key in keys:
        folder = library_dir(lab, key)
        folder.mkdir(parents=True, exist_ok=True)
        cram = folder / f"{key}.cram"
        cram.write_bytes(b"")
        record = Record(
            cram=str(cram),
            key=f"align-{key}",
            metrics={"mean_depth": 30.0},
            gvcf=None,
            reference_sha256=sha256_file(fasta),
            accepted=True,
            snp_f1=None,
            replaced=None,
        )
        write_record(lab, key, record)
    return Paths(repo=tmp_path, raw=tmp_path / "raw", reports=tmp_path / "r", lab=lab)


class _Steps:
    """Stands for `stage3_steps` and `repeat_steps` in `build`: records each call's sample
    folder, step and sex, and gives HG002's slice-3e output (no SMN or Yleaf result)."""

    def __init__(self, lab: Path) -> None:
        self.lab = lab
        self.calls: list[tuple[str, str, object]] = []
        hg002 = _hg002()
        self.results = replace(
            hg002, caller=None, parascopy=None, yleaf=None, yleaf_out=None
        )

    def _key(self, work: Path) -> str:
        (key,) = [k for k in RECORDED if library_dir(self.lab, k) == work]
        return key

    def stage3(self, work: Path, **kwargs: Any) -> tuple[Any, dict[str, Any]]:
        self.calls.append((self._key(work), "stage3_steps", kwargs["sex"]))
        assert "lab" in kwargs
        return self.results, {"repeats": SimpleNamespace(key=f"r-{work.name}")}

    def repeats(self, work: Path, **kwargs: Any) -> tuple[Any, dict[str, Any]]:
        self.calls.append((self._key(work), "repeat_steps", kwargs["sex"]))
        assert "lab" not in kwargs
        found = SimpleNamespace(
            repeats=self.results.repeats,
            rfc1=self.results.rfc1,
            images={},
            flank=1.0,
            failed={},
        )
        return found, {"repeats": SimpleNamespace(key=f"r-{work.name}")}


def _fake_build(monkeypatch: pytest.MonkeyPatch, paths: Paths) -> _Steps:
    steps = _Steps(paths.lab)
    monkeypatch.setattr(build, "stage3_versions", lambda *_a: dict(VERSIONS))
    monkeypatch.setattr(build, "stage3_steps", steps.stage3)
    monkeypatch.setattr(build, "repeat_steps", steps.repeats)
    monkeypatch.setattr(build, "truth_repeats", lambda *_a, **_k: None)
    return steps


def test_build_stage3_refuses_unaligned_sample(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = _build_lab(tmp_path, ("pcrfree", "pcr", "hg003"))
    steps = _fake_build(monkeypatch, paths)
    with pytest.raises(RunRefused, match=r"HG004 is not aligned; run .*align hg004"):
        build.build_stage3(paths, threads=1)
    assert steps.calls == []  # refused before any step
    assert not stage3_dir(paths.lab).exists()


def test_build_stage3_runs_each_sample_under_its_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every sample's steps run inside `run_lock(library_dir(lab, key))`."""
    paths = _build_lab(tmp_path, RECORDED)
    steps = _fake_build(monkeypatch, paths)
    held: list[Path] = []

    @contextmanager
    def lock(work: Path) -> Iterator[None]:
        held.append(work)
        try:
            yield
        finally:
            held.remove(work)

    def locked(run: Any) -> Any:
        def step(work: Path, **kwargs: Any) -> Any:
            assert held == [work]
            return run(work, **kwargs)

        return step

    monkeypatch.setattr(build, "run_lock", lock)
    monkeypatch.setattr(build, "stage3_steps", locked(steps.stage3))
    monkeypatch.setattr(build, "repeat_steps", locked(steps.repeats))
    build.build_stage3(paths, threads=1)
    assert [key for key, _, _ in steps.calls] == list(RECORDED)


def test_build_stage3_refuses_a_parents_failed_repeat_step(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = _build_lab(tmp_path, RECORDED)
    steps = _fake_build(monkeypatch, paths)

    def repeats(work: Path, **kwargs: Any) -> Any:
        found, kept = steps.repeats(work, **kwargs)
        if work == library_dir(paths.lab, "hg003"):
            found.failed = {"expansionhunter": "a truncated CRAM"}
        return found, kept

    monkeypatch.setattr(build, "repeat_steps", repeats)
    with pytest.raises(
        RunRefused,
        match="HG003: expansionhunter failed: a truncated CRAM; no stage 3 table",
    ):
        build.build_stage3(paths, threads=1)
    assert not stage3_dir(paths.lab).exists()


def test_build_stage3_parents_run_repeats_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """HG002's libraries get the three stage-3 steps, HG003 and HG004 the repeat step alone
    with their own sex; one table holds 3a's and the rules' withheld keys merged, the lower
    cap of both and the rules' measurements under `RULES`."""
    paths = _build_lab(tmp_path, RECORDED)
    steps = _fake_build(monkeypatch, paths)
    measure_3a = build.measure

    def measure(*args: Any) -> Any:
        withheld, caps, measured = measure_3a(*args)
        return withheld, caps | {"HTT": Level.MEDIUM}, measured

    measure_rules = build.measure_repeat_rules

    def rules(*args: Any, **kwargs: Any) -> Any:
        withheld, caps, measured = measure_rules(*args, **kwargs)
        return withheld, caps | {"HTT": Level.LOW, "TBP": Level.LOW}, measured

    monkeypatch.setattr(build, "measure", measure)
    monkeypatch.setattr(build, "measure_repeat_rules", rules)
    table = build.build_stage3(paths, threads=1)
    assert steps.calls == [
        ("pcrfree", "stage3_steps", "male"),
        ("pcr", "stage3_steps", "male"),
        ("hg003", "repeat_steps", "male"),
        ("hg004", "repeat_steps", "female"),
    ]
    # 3a's SMN and Y haplogroup reasons beside the rules' (no SMN or Yleaf result here).
    assert table.withheld[SMN_KEY] == SMN_WITHHELD
    assert table.withheld[YLEAF_KEY] == YLEAF_WITHHELD
    # HG002's fixture does not see its benign ~115-unit AAAAG allele: RFC1 withheld.
    assert table.withheld[RFC1_LOCUS] == RFC1_WITHHELD
    assert table.caps["HTT"] == Level.LOW and table.caps["TBP"] == Level.LOW
    assert set(table.measured) == {"pcrfree", "pcr", RULES}
    assert set(table.measured[RULES]) == set(RECORDED)
    assert open_stage3(paths.lab, VERSIONS) == (table, None)
