"""Slice 3d's calibration, its pure part (`genokeeper.calibrate.slice3d`): the gain truth's
k-mer rule on a hand-made FASTA, trio judging, spike-in recall, the CNVpytor decision, the
paralog loci withheld against HG002's truth, the large bands pooled per genome, and which
inputs each table's versions name."""

import gzip
import random
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from genokeeper.calibrate import slice3d, spikein
from genokeeper.calibrate.classes import ClassTable, Geometry, measure, summary
from genokeeper.calibrate.deletions import deletion_versions
from genokeeper.calibrate.records import RECORDED
from genokeeper.calibrate.slice3d import (
    CNV_BANDS,
    CNV_FEATURES,
    INSERTION_BANDS,
    INSERTION_FEATURES,
    TRIO,
    ParalogTable,
    Truth,
    classify_gain_truth,
    cnv_versions,
    cnvpytor_kept,
    insertion_versions,
    large_withheld,
    load_kiv2_truth,
    load_truth,
    open_paralog_table,
    paralog_measure,
    paralog_versions,
    pool_large,
    read_truth,
    report,
    spike_recall,
    trio_judged,
    write_paralog_table,
)
from genokeeper.cnv.screen import LargeEvent
from genokeeper.config import CnvThresholds, Grading, LpaThresholds, load_grading
from genokeeper.containers import IMAGES
from genokeeper.deletions.files import sv_files
from genokeeper.deletions.screen import CallClass, Deletion
from genokeeper.evidence.clingen import RECURRENT_FILE, REGIONS_FILE
from genokeeper.masking import masked_bed
from genokeeper.paralogs import GauchianCall
from genokeeper.regions import giab38_dir

SAMTOOLS = shutil.which("samtools")
G = Geometry((("chr1", 100_000, 100_250),), 1_000_000)


def _event(
    chrom: str, pos: int, end: int, svtype: str = "DEL", caller: str = "truth"
) -> Deletion:
    return Deletion(chrom, pos, end, "N", "0/1", {caller: "x"}, None, svtype=svtype)


def _insertion(pos: int, length: int) -> Deletion:
    return Deletion(
        "chr1",
        pos,
        pos + 1,
        "N",
        "0/1",
        {"truth": "x"},
        None,
        svtype="INS",
        svlen=length,
    )


def test_constants() -> None:
    assert tuple(TRIO) == RECORDED  # 3c's library table
    assert TRIO == {
        "pcrfree": "HG002",
        "pcr": "HG002",
        "hg003": "HG003",
        "hg004": "HG004",
    }
    assert CNV_FEATURES == {
        "evidence": ("both", "dellycnv", "cnvpytor"),
        "junction": ("true", "false"),
        "masked": ("true", "false"),
    }
    assert INSERTION_FEATURES == {
        "evidence": ("both", "dysgu", "insurveyor", "scramble"),
        "repeat": ("true", "false"),
        "masked": ("true", "false"),
    }
    assert INSERTION_BANDS == ("short", "long", "mei")
    assert CNV_BANDS == ("loss_1-5kb", "loss_5kb+", "gain_1-5kb", "gain_5kb+")


# --- gain truth ---------------------------------------------------------------------------


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools is not on PATH")
def test_gain_truth_kmer_rule(tmp_path: Path) -> None:
    """A copy of the next 1,200 bp is a DUP after POS, a copy of the previous 1,200 bp one
    before it, a novel sequence none; 79% of the k-mers contained is none, 80% one."""
    rng = random.Random(7)
    sequence = "".join(rng.choice("ACGT") for _ in range(30_000))
    fasta = tmp_path / "ref.fa"
    lines = [sequence[i : i + 60] for i in range(0, len(sequence), 60)]
    fasta.write_text(">chr1\n" + "\n".join(lines) + "\n")
    novel = "".join(rng.choice("ACGT") for _ in range(1200))
    # VCF POS p: the anchor base is sequence[p - 1]; the insertion follows it.
    after, before, new, partial, enough = 5_000, 12_000, 20_000, 25_000, 8_000
    found = classify_gain_truth(
        [
            (_insertion(after, 1200), sequence[after : after + 1200]),
            (_insertion(before, 1200), sequence[before - 1200 : before]),
            (_insertion(new, 1200), novel),
            # 954 copied bases: 924 of 1,170 k-mers (79%) are the reference's.
            (
                _insertion(partial, 1200),
                sequence[partial : partial + 954] + novel[:246],
            ),
            # 966 copied bases: 936 of 1,170 (80%).
            (_insertion(enough, 1200), sequence[enough : enough + 966] + novel[:234]),
        ],
        fasta,
    )
    assert [(d.pos, d.end, d.svtype) for d in found] == [
        (after, after + 1200, "DUP"),
        (before - 1200, before, "DUP"),
        (enough, enough + 1200, "DUP"),
    ]


def test_gain_truth_skips_short_and_bounded(tmp_path: Path) -> None:
    """An insertion under 1 kb or of a lower-bound length is never a gain truth (no
    reference read)."""
    short = _insertion(100, 999)
    bounded = replace(_insertion(100, 1500), lower_bound=True)
    assert (
        classify_gain_truth([(short, "A" * 999), (bounded, "A" * 1500)], tmp_path) == []
    )


# --- trio, spike-ins, CNVpytor ------------------------------------------------------------


def test_trio_judged() -> None:
    loss = _event("chr1", 1_000, 5_000)
    gain = _event("chr1", 1_000, 5_000, "DUP")
    other = _event("chr1", 50_000, 60_000)
    common = _event("chr1", 90_000, 95_000)
    hg003 = [_event("chr1", 1_100, 5_000, "DUP"), common]
    hg004 = [_event("chr1", 1_000, 4_900), _event("chr1", 50_000, 60_000), common]
    hg003 += [_event("chr1", 50_000, 59_000)]
    judged = trio_judged(
        [loss, gain, other, common],
        [hg003, hg004],
        min_overlap=0.5,
        common=lambda d: d.pos == 90_000,
    )
    assert judged == [
        (loss, True),  # in HG004 only (HG003's is a gain)
        (gain, True),  # in HG003 only
        (other, False),  # in both parents
        (common, False),  # common
    ]
    assert trio_judged([other], [[], []], min_overlap=0.5, common=lambda _: False) == [
        (other, False)
    ]


def test_spike_recall_by_band() -> None:
    t = CnvThresholds()
    one = _event("chr1", 100_000, 102_000)
    two = _event("chr1", 300_000, 302_000)
    three = _event("chr1", 500_000, 502_000, "DUP")
    found = [
        _event("chr1", 100_100, 102_000),
        _event("chr1", 300_000, 302_000, "DUP"),  # the other direction
        _event("chr1", 500_000, 509_000, "DUP"),  # overlap under 0.5
    ]
    assert spike_recall(found, [one, two, three], t, min_overlap=0.5) == [
        ("loss_1-5kb", True),
        ("loss_1-5kb", False),
        ("gain_1-5kb", False),
    ]


def _cnv_calls(evidence: str, n: int, false: int) -> list[tuple[CallClass, int, bool]]:
    c = CallClass("loss_1-5kb", evidence, False, False, False)
    return [(c, 2000, False)] * false + [(c, 2000, True)] * (n - false)


def test_cnvpytor_kept_rule() -> None:
    """Both callers' false share below DELLY's (every call DELLY made) keeps CNVpytor;
    equal or above drops it; the four counts are recorded. CNVpytor-only calls count in
    neither."""
    pytor = _cnv_calls("cnvpytor", 50, 40)
    kept, record = cnvpytor_kept(
        {
            "pcrfree": _cnv_calls("both", 100, 2) + pytor,
            "pcr": _cnv_calls("dellycnv", 100, 28),
        }
    )
    assert kept
    assert record == {
        "decision": "kept",
        "both": {"calls": 100, "false": 2},
        "delly": {"calls": 200, "false": 30},
    }
    for both_false, decision in ((15, "dropped"), (16, "dropped"), (14, "kept")):
        kept, record = cnvpytor_kept(
            {
                "pcrfree": _cnv_calls("both", 100, both_false),
                "pcr": _cnv_calls("dellycnv", 100, 30 - both_false),
            }
        )
        assert record["decision"] == decision and kept == (decision == "kept")
        assert record["delly"] == {"calls": 200, "false": 30}


# --- paralogs -----------------------------------------------------------------------------

WILD = GauchianCall(4, False, False, (), (), None)
CARRIER = GauchianCall(4, True, False, ("L483P",), (), 12)
LPA = LpaThresholds(max_kiv2_error=0.20)
LOCI = {"STRC": 2, "RHD": 2, "AMY1": 6, "C4": 4, "CYP2A6": 2, "GBA1": 4}


def _consistent() -> dict[str, dict[str, int | None]]:
    return {library: dict(LOCI) for library in TRIO}


def _measure(
    per_sample: dict[str, dict[str, int | None]] | None = None,
    truth: dict[tuple[str, str], Truth] | None = None,
    gauchian: dict[str, GauchianCall | None] | None = None,
    kiv2: dict[str, float | None] | None = None,
    kiv2_truth: float | None = 38.0,
) -> dict[str, str]:
    withheld, _ = paralog_measure(
        _consistent() if per_sample is None else per_sample,
        load_truth() if truth is None else truth,
        {"pcrfree": WILD, "pcr": WILD} if gauchian is None else gauchian,
        {"pcrfree": 38.0, "pcr": 41.0} if kiv2 is None else kiv2,
        kiv2_truth,
        LPA,
    )
    return withheld


def test_paralog_measure_withholds() -> None:
    assert _measure() == {}
    _, measured = paralog_measure(
        _consistent(), load_truth(), {"pcrfree": WILD, "pcr": WILD},
        {"pcrfree": 38.0, "pcr": 41.0}, 38.0, LPA,
    )  # fmt: skip
    assert measured["pcr"]["copies"]["AMY1"] == 6 and measured["pcr"]["kiv2"] == 41.0
    assert measured["pcrfree"]["gba1_wild_type"] is True

    disagree = _consistent()
    disagree["pcr"]["RHD"] = 1
    assert _measure(disagree) == {"RHD": "not covered: HG002's libraries disagree"}

    off = _consistent()
    off["pcrfree"]["C4"] = off["pcr"]["C4"] = 3
    assert _measure(off) == {"C4": "not covered: pcrfree gave 3, truth 4"}

    exceeds = _consistent()
    for library, cn in (("pcrfree", 5), ("pcr", 5), ("hg003", 2), ("hg004", 2)):
        exceeds[library]["AMY1"] = cn
    five = load_truth() | {("AMY1", "HG002"): Truth(5, None)}
    assert _measure(exceeds, five) == {
        "AMY1": "not covered: HG002 exceeds its parents' copies"
    }

    parent = _consistent()
    parent["hg003"]["STRC"] = 1
    assert _measure(parent) == {"STRC": "not covered: hg003 gave 1, truth 2"}

    missing = _consistent()
    missing["hg004"]["CYP2A6"] = None
    assert _measure(missing) == {"CYP2A6": "not covered: no call in calibration"}
    del missing["hg004"]
    assert set(_measure(missing).values()) == {"not covered: no call in calibration"}

    assert _measure(gauchian={"pcrfree": WILD, "pcr": CARRIER}) == {
        "GBA1": "not covered: GBA1 not wild type on HG002 pcr"
    }
    assert _measure(gauchian={"pcrfree": WILD, "pcr": None}) == {
        "GBA1": "not covered: no call in calibration"
    }

    assert _measure(kiv2={"pcrfree": 38.0, "pcr": 46.0}) == {
        "LPA": "not covered: KILDA off HG002's count by 8.0 (21%)"
    }
    assert _measure(kiv2={"pcrfree": 38.0, "pcr": None}) == {
        "LPA": "not covered: no call in calibration"
    }

    no_rhd = {k: v for k, v in load_truth().items() if k[0] != "RHD"}
    assert _measure(truth=no_rhd) == {
        "RHD": "not calibrated: no HG002 truth: no row in the truth resource"
    }


TRUTH_HEADER = (
    "# a comment\n# another\n"
    "locus\tsample\thaplotype_1\thaplotype_2\ttotal\tsource\tmethod\treason\n"
)


def test_load_truth_empty_total() -> None:
    text = TRUTH_HEADER + (
        "STRC\tHG002\t1\t1\t2\tGIAB v5.0q\tlook-up\t\n"
        "KIV2\tHG002\t\t\t\tHG002 Q100 v1.1\tlook-up\tassembly broken inside LPA\n"
    )
    truth = read_truth(text)
    assert truth[("STRC", "HG002")] == Truth(2, None)
    assert truth[("KIV2", "HG002")] == Truth(None, "assembly broken inside LPA")
    withheld = _measure(truth=load_truth() | truth, kiv2_truth=None)
    assert withheld == {
        "LPA": "not calibrated: no HG002 truth: assembly broken inside LPA"
    }
    with pytest.raises(ValueError, match="KIV2"):
        read_truth(TRUTH_HEADER + "KIV2\tHG002\t\t\t\tQ100\tlook-up\t\n")
    packaged = load_truth()
    assert packaged[("STRC", "HG004")] == Truth(2, None)
    assert load_kiv2_truth() == packaged[("KIV2", "HG002")] == Truth(38, None)


def test_paralog_table_roundtrip(tmp_path: Path) -> None:
    table = ParalogTable(
        versions={"gauchian": "g"}, withheld={"RHD": "why"}, measured={"pcr": {"x": 1}}
    )
    write_paralog_table(tmp_path / "lab", table, report="r")
    assert open_paralog_table(tmp_path / "lab", {"gauchian": "g"}) == (table, None)
    found, reason = open_paralog_table(tmp_path / "lab", {"gauchian": "h"})
    assert (
        found is None
        and reason == "stale: gauchian changed; run `genokeeper calibrate`"
    )


# --- large events -------------------------------------------------------------------------


def _large(chrom: str, pos: int, end: int, svtype: str = "DEL") -> LargeEvent:
    band = "loss_large" if svtype == "DEL" else "gain_large"
    event = _event(chrom, pos, end, svtype, "dellycnv")
    return LargeEvent(event, CallClass(band, "both", False, False, False), ("22q11.2",))


def test_pool_large() -> None:
    pcrfree = _large("chr22", 18_000_000, 18_600_000)
    pcr = _large("chr22", 18_000_000, 18_540_000)  # reciprocal overlap 0.9
    pooled = pool_large(
        {"pcrfree": [pcrfree], "pcr": [pcr], "hg003": [], "hg004": []},
        min_overlap=0.5,
    )
    assert pooled["HG002"] == {"loss_large": 1, "gain_large": 0}
    assert pooled["HG003"] == pooled["HG004"] == {"loss_large": 0, "gain_large": 0}
    assert large_withheld(pooled) == {}
    pooled = pool_large(
        {
            "pcrfree": [pcrfree],
            "pcr": [pcr, _large("chr22", 18_000_000, 18_600_000, "DUP")],
            "hg003": [_large("chr1", 1_000_000, 1_700_000)],
            "hg004": [],
        },
        min_overlap=0.5,
    )
    assert pooled["HG002"] == {"loss_large": 1, "gain_large": 1}
    assert large_withheld(pooled) == {
        "loss_large": "not covered: large loss calls on 2 of 3 calibration genomes"
    }


# --- versions -----------------------------------------------------------------------------


def _lab(tmp_path: Path) -> Path:
    """The files the four versions functions hash, with any content."""
    lab = tmp_path / "lab"
    giab = giab38_dir(lab)
    giab.mkdir(parents=True)
    for name in load_grading().strata.values():
        (giab / name).write_bytes(gzip.compress(b"chr20\t100\t200\n"))
    gnomad = sv_files(lab).gnomad
    gnomad.parent.mkdir(parents=True)
    gnomad.write_bytes(b"gnomad")
    bed = masked_bed(lab)
    bed.parent.mkdir(parents=True)
    bed.write_text("chr21\t5010000\t5166246\n")
    table = lab / "parascopy" / "homology_table"
    table.mkdir(parents=True)
    (table / "GRCh38.md5").write_text("abc GRCh38.bed.gz\n")
    (lab / "parascopy" / ".unpacked.models").write_text("def\n")
    return lab


def _all(lab: Path, grading: Grading) -> dict[str, dict[str, str]]:
    sha = "f" * 64
    return {
        "cnv": cnv_versions(lab, sha, grading),
        "insertions": insertion_versions(lab, sha, grading),
        "paralogs": paralog_versions(lab, sha, grading),
        "deletions": deletion_versions(lab, sha, grading),
    }


def test_versions_name_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    lab, grading = _lab(tmp_path), load_grading()
    before = _all(lab, grading)

    def changed(after: dict[str, dict[str, str]]) -> set[str]:
        return {name for name in before if after[name] != before[name]}

    assert before["cnv"]["cnv"] == '{"band_from": 5000, "bin_size": 500,' + (
        ' "large_from": 500000, "min_length": 1000}'
    )
    monkeypatch.setitem(IMAGES, "gauchian", "genomelab/gauchian:other")
    assert changed(_all(lab, grading)) == {"paralogs"}
    monkeypatch.undo()
    cnv = replace(grading, cnv=replace(grading.cnv, min_length=2000))
    assert changed(_all(lab, cnv)) == {"cnv"}
    masked_bed(lab).write_text("chr21\t5010000\t5166247\n")
    assert changed(_all(lab, grading)) == {"cnv"}
    masked_bed(lab).write_text("chr21\t5010000\t5166246\n")
    # `[par]` decides a male's loss or gain, and the spike-in rule the depth recall.
    par = replace(grading, par={**grading.par, "chrX": [(10_001, 2_781_479)]})
    assert changed(_all(lab, par)) == {"cnv"}
    monkeypatch.setattr(spikein, "N_PER", 19)
    assert changed(_all(lab, grading)) == {"cnv"}
    monkeypatch.undo()
    monkeypatch.setattr(slice3d, "SCRAMBLE_MARGIN", 600)
    assert changed(_all(lab, grading)) == {"insertions"}
    monkeypatch.undo()
    clingen = lab / "clingen-dosage" / "20261101"
    clingen.mkdir(parents=True)
    (clingen / REGIONS_FILE).write_text("a new region list\n")
    (clingen / RECURRENT_FILE).write_text("a new recurrent list\n")
    assert changed(_all(lab, grading)) == set()
    masked_bed(lab).unlink()
    with pytest.raises(FileNotFoundError):
        cnv_versions(lab, "f" * 64, grading)


def test_slice3d_imports_no_pipeline() -> None:
    """`pipeline` may import the module: it imports nothing that imports `pipeline`."""
    code = (
        "import sys, genokeeper.calibrate.slice3d, genokeeper.calibrate.spikein;"
        " sys.exit('genokeeper.pipeline' in sys.modules)"
    )
    assert subprocess.run([sys.executable, "-c", code], check=False).returncode == 0


def test_report() -> None:
    """The report names the CNVpytor decision with its counts, the large calls per genome
    and the band withheld, each class's `p_wrong`, and each withheld paralog locus."""
    c = CallClass("loss_1-5kb", "both", False, False, False)
    calls = [(c, 2000, False)] * 2 + [(c, 2000, True)] * 98
    truth: list[tuple[str | None, bool]] = [("loss_1-5kb", True)] * 9
    geometry = {b: G for b in CNV_BANDS}
    bands, likelihoods = measure(
        {"pcrfree": calls}, {"pcrfree": truth}, geometry, dict.fromkeys(CNV_BANDS, 0.1),
        CNV_FEATURES, band=lambda x: x.band,
    )  # fmt: skip
    _, record = cnvpytor_kept({"pcrfree": calls})
    pooled = {
        "HG002": {"loss_large": 1, "gain_large": 0},
        "HG003": {"loss_large": 1, "gain_large": 0},
    }
    measured = summary(
        calls, truth, CNV_FEATURES, CNV_BANDS, band=lambda x: x.band,
        key=lambda x: "loss_1-5kb:both:norepeat:unmasked:nojunction",
    )  # fmt: skip
    cnv = ClassTable(
        {"delly": "d"}, bands, likelihoods, {"pcrfree": measured}, CNV_FEATURES,
        {"cnvpytor": record, "large": {"calls": pooled, "withheld": large_withheld(pooled)}},
    )  # fmt: skip
    paralogs = ParalogTable({"gauchian": "g"}, {"RHD": "not covered: why"}, {})
    text = report(cnv, None, paralogs)
    assert "CNVpytor dropped: calls both depth callers made 2 false of 100" in text
    assert "large calls on HG003: loss_large 1, gain_large 0" in text
    assert (
        "loss_large: not covered: large loss calls on 2 of 2 calibration genomes"
        in text
    )
    assert "loss_1-5kb:both:norepeat:unmasked:nojunction: p_wrong " in text
    assert "  RHD: not covered: why" in text
    assert "insertions table" not in text
