"""ExpansionHunter and REViewer (`genokeeper.repeats.expansionhunter`): the catalog cut of the
eleven loci, ExpansionHunter's calls as read from its HG002 output (tests/fixtures/hg002),
and the two container runs, whose container is replaced by a stub that checks the arguments
and copies the fixture into the output directory."""

import json
import os
import re
import shutil
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any

import pytest

from genokeeper.config import load_grading
from genokeeper.repeats import expansionhunter
from genokeeper.repeats.expansionhunter import (
    ADDED,
    CATALOG,
    EXTRA,
    REPLACED,
    cut_catalog,
    flank_ratio,
    read_calls,
    region,
    review,
    run_expansionhunter,
)
from genokeeper.repeats.not_covered import NOT_COVERED
from genokeeper.repeats.rfc1 import RFC1_LOCUS
from genokeeper.sources.repeats import CFTR_T, CFTR_TG

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
VCF = FIXTURES / "hg002.eh3e.vcf.txt"
JSON = FIXTURES / "hg002.eh3e.json"
SVG = FIXTURES / "hg002.reviewer.HTT.svg"
LOCI = load_grading().repeats.loci
ELEVEN = (
    "HTT", "AR", "DMPK", "C9orf72", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A",
    "ATN1", "FXN",
)  # fmt: skip
SOURCE_MD = files("genokeeper").joinpath("resources", "repeats", "SOURCE.md")


def _catalog() -> list[dict[str, Any]]:
    return json.loads(CATALOG.read_text())


def _variant_ids(entry: dict[str, Any]) -> set[str]:
    """The VARIDs ExpansionHunter 5.0.0 writes for a catalog entry: its `VariantId`s, or
    its LocusId for a one-repeat entry without them."""
    return set(entry.get("VariantId") or [entry["LocusId"]])


def test_catalog_has_the_loci() -> None:
    catalog = _catalog()
    assert len(catalog) == 26
    assert [e["LocusId"] for e in catalog].count("RFC1") == 1
    assert {e["LocusId"] for e in catalog} == {l.locus_id for l in LOCI.values()} | set(
        EXTRA
    )
    for locus in LOCI.values():
        (entry,) = (e for e in catalog if e["LocusId"] == locus.locus_id)
        assert locus.variant_id in _variant_ids(entry)
        regions = entry["ReferenceRegion"]
        first = regions if isinstance(regions, str) else regions[0]
        assert first.split(":")[0] == locus.chrom
    assert all("OfftargetRegions" not in e for e in catalog)
    (pabpn1,) = (e for e in catalog if e["LocusId"] == "PABPN1")
    assert pabpn1["LocusStructure"] == "(GCN)*"
    # RFC1 and CFTR close the file, in that order.
    assert [e["LocusId"] for e in catalog[-2:]] == ["RFC1", "CFTR"]


# The five entries the Illumina catalog lacks or counts short, from the slice's values: the
# four trimmed to whole units from STRchive v2.28.0's ExpansionHunter asset, and CFTR's
# TG/T tract read off GRCh38 and ClinVar.
ADDED_ENTRIES = [
    {
        "LocusId": "PABPN1",
        "LocusStructure": "(GCN)*",
        "ReferenceRegion": "chr14:23321472-23321502",
        "VariantType": "Repeat",
    },
    {
        "LocusId": "ZFHX3",
        "LocusStructure": "(CCG)*",
        "ReferenceRegion": "chr16:72787695-72787755",
        "VariantType": "Repeat",
    },
    {
        "LocusId": "THAP11",
        "LocusStructure": "(CAG)*",
        "ReferenceRegion": "chr16:67842863-67842950",
        "VariantType": "Repeat",
    },
    {
        "LocusId": "FGF14",
        "LocusStructure": "(AAG)*",
        "ReferenceRegion": "chr13:102161574-102161724",
        "VariantType": "Repeat",
    },
    {
        "LocusId": "CFTR",
        "LocusStructure": "(TG)*(T)*",
        "ReferenceRegion": ["chr7:117548606-117548628", "chr7:117548628-117548635"],
        "VariantId": ["CFTR_TG", "CFTR_T"],
        "VariantType": ["Repeat", "Repeat"],
    },
]


def _units(structure: str) -> list[str]:
    return re.findall(r"\(([A-Z]+)\)\*", structure)


def test_added_regions() -> None:
    added = json.loads(ADDED.read_text())
    assert added == ADDED_ENTRIES
    for entry in added:
        found = entry["ReferenceRegion"]
        regions = [found] if isinstance(found, str) else found
        for span, unit in zip(regions, _units(entry["LocusStructure"]), strict=True):
            start, end = (int(x) for x in span.split(":")[1].split("-"))
            assert (end - start) % len(unit) == 0, (entry["LocusId"], span)


def test_added_regions_start_and_end_on_a_unit() -> None:
    """Each added region reads whole units of its motif in the reference (N any base)."""
    path = os.environ.get("GRCH38_FASTA")
    if path is None:
        pytest.skip("GRCH38_FASTA names no local GRCh38 FASTA")
    import pysam

    with pysam.FastaFile(path) as fasta:
        for entry in json.loads(ADDED.read_text()):
            found = entry["ReferenceRegion"]
            regions = [found] if isinstance(found, str) else found
            units = _units(entry["LocusStructure"])
            for span, unit in zip(regions, units, strict=True):
                chrom, bounds = span.split(":")
                start, end = (int(x) for x in bounds.split("-"))
                seq = fasta.fetch(chrom, start, end).upper()
                pattern = unit.replace("N", "[ACGT]")
                head, tail = seq[: len(unit)], seq[-len(unit) :]
                assert re.fullmatch(pattern, head), (entry["LocusId"], span, head)
                assert re.fullmatch(pattern, tail), (entry["LocusId"], span, tail)


def test_rfc1_entry() -> None:
    (entry,) = (e for e in _catalog() if e["LocusId"] == "RFC1")
    assert entry == {
        "LocusId": "RFC1",
        "LocusStructure": "(AARRG)*",
        "ReferenceRegion": "chr4:39348424-39348479",
        "VariantType": "Repeat",
        "RFC1MotifAnalysis": True,
    }


def test_strchive_thresholds_match_config() -> None:
    path = os.environ.get("STRCHIVE_LOCI")
    if path is None:
        pytest.skip("STRCHIVE_LOCI names no local STRchive-loci.json")
    strchive = {e["gene"]: e for e in json.loads(Path(path).read_text())}
    for gene, locus in LOCI.items():
        entry = strchive[gene]
        assert entry["benign_max"] == locus.normal_max, gene
        # STRchive's intermediate range is the locus's reduced-penetrance one where the
        # config has no intermediate range (AR); STRchive's (None, None) is none.
        span = locus.intermediate or locus.reduced_penetrance
        found = (entry["intermediate_min"], entry["intermediate_max"])
        if gene == "ATXN10":
            # STRchive's 33-850 is cut below its own pathogenic minimum of 800.
            assert (found[0], locus.intermediate) == (33, (33, 799))
        else:
            assert (None if found == (None, None) else found) == span, gene
        if gene != "HTT":
            assert entry["pathogenic_min"] == locus.pathogenic_min, gene


RFC1_SOURCE = {
    "LocusId": "RFC1",
    "LocusStructure": "(AARRG)*",
    "ReferenceRegion": "4:39348424-39348479",
    "VariantType": "Repeat",
}
CFTR_ADDED = ADDED_ENTRIES[-1]


def _files(
    tmp_path: Path, source: list[dict[str, Any]], added: list[dict[str, Any]]
) -> tuple[Path, Path]:
    paths = tmp_path / "variant_catalog.json", tmp_path / "added.json"
    for path, entries in zip(paths, (source, added), strict=True):
        path.write_text(json.dumps(entries))
    return paths


def test_cut_keeps_the_loci_unchanged(tmp_path: Path) -> None:
    other = {
        "LocusId": "OTHER_X",
        "LocusStructure": "(A)*",
        "ReferenceRegion": "chr1:1-2",
    }
    wanted = {
        "LocusId": "HTT",
        "LocusStructure": "(CAG)*CAACAG(CCG)*",
        "ReferenceRegion": ["4:1-2", "4:8-9"],
        "VariantId": ["HTT", "HTT_CCG"],
        "VariantType": ["Repeat", "Repeat"],
        "OfftargetRegions": ["chr1:5-6"],
    }
    source, added = _files(tmp_path, [other, wanted, other, RFC1_SOURCE], [CFTR_ADDED])
    cut = cut_catalog(source, [LOCI["HTT"]], added=added)
    # Only the off-target regions go, and the contigs take the reference's `chr` names;
    # RFC1 gains its motif analysis and CFTR comes last.
    assert cut == [
        {k: v for k, v in wanted.items() if k != "OfftargetRegions"}
        | {"ReferenceRegion": ["chr4:1-2", "chr4:8-9"]},
        RFC1_SOURCE
        | {"ReferenceRegion": "chr4:39348424-39348479", "RFC1MotifAnalysis": True},
        CFTR_ADDED,
    ]


def test_cut_replaces_pabpn1(tmp_path: Path) -> None:
    assert frozenset({"PABPN1"}) == REPLACED
    illumina_pabpn1 = {
        "LocusId": "PABPN1",
        "LocusStructure": "(GCG)*",
        "ReferenceRegion": "14:23321472-23321490",
        "VariantType": "Repeat",
    }
    htt = {"LocusId": "HTT", "LocusStructure": "(CAG)*", "ReferenceRegion": "4:1-2"}
    fgf14, pabpn1 = ADDED_ENTRIES[3], ADDED_ENTRIES[0]
    source, added = _files(
        tmp_path,
        [illumina_pabpn1, htt, RFC1_SOURCE],
        [fgf14, CFTR_ADDED, pabpn1],
    )
    loci = [LOCI["HTT"], LOCI["PABPN1"], LOCI["FGF14"]]
    cut = cut_catalog(source, loci, added=added)
    # PABPN1 in Illumina's place, then the added loci in their file's order.
    assert [e["LocusId"] for e in cut] == ["PABPN1", "HTT", "FGF14", "RFC1", "CFTR"]
    assert cut[0] == pabpn1


@pytest.mark.parametrize(
    ("source", "added", "loci", "named"),
    [
        # A LocusId in both files outside REPLACED.
        (
            [{"LocusId": "HTT", "ReferenceRegion": "4:1-2"}, RFC1_SOURCE],
            [{"LocusId": "HTT", "ReferenceRegion": "chr4:1-2"}, CFTR_ADDED],
            ["HTT"],
            "HTT",
        ),
        # A locus in neither file.
        ([RFC1_SOURCE], [CFTR_ADDED], ["HTT"], "HTT"),
        # A REPLACED name missing from the added file.
        (
            [{"LocusId": "PABPN1", "ReferenceRegion": "14:1-2"}, RFC1_SOURCE],
            [CFTR_ADDED],
            ["PABPN1"],
            "PABPN1",
        ),
        # An EXTRA entry missing.
        (
            [{"LocusId": "HTT", "ReferenceRegion": "4:1-2"}],
            [CFTR_ADDED],
            ["HTT"],
            "RFC1",
        ),
        (
            [{"LocusId": "HTT", "ReferenceRegion": "4:1-2"}, RFC1_SOURCE],
            [],
            ["HTT"],
            "CFTR",
        ),
    ],
    ids=["duplicate", "absent", "replaced-missing", "rfc1-missing", "cftr-missing"],
)
def test_cut_refuses_other_duplicates(
    tmp_path: Path,
    source: list[dict[str, Any]],
    added: list[dict[str, Any]],
    loci: list[str],
    named: str,
) -> None:
    paths = _files(tmp_path, source, added)
    with pytest.raises(ValueError, match=named):
        cut_catalog(paths[0], [LOCI[g] for g in loci], added=paths[1])


def test_regions_are_illuminas() -> None:
    """HTT's CAG region is GRCh38's whole 19-unit tract and AR's the 23-unit GCA run (22 CAG
    and the CAA), as Illumina's ExpansionHunter v5.0.0 catalog has them; STRchive's stopped
    two CAG short at HTT and reached past AR's tract."""
    assert region(LOCI["HTT"]) == ("chr4", 3074877, 3074933)
    assert region(LOCI["AR"]) == ("chrX", 67545317, 67545385)
    assert region(LOCI["C9orf72"]) == ("chr9", 27573529, 27573546)


def test_config_motif_is_the_catalog_motif() -> None:
    """A locus whose catalog structure is one repeat has the config's motif, in its phase."""
    for locus in LOCI.values():
        (entry,) = (e for e in _catalog() if e["LocusId"] == locus.locus_id)
        structure = entry["LocusStructure"]
        if structure.count("(") == 1 and structure.endswith(")*"):
            assert structure == f"({locus.motif})*", locus.gene


def test_catalog_resource_is_the_cut_of_illumina() -> None:
    path = os.environ.get("ILLUMINA_CATALOG")
    if path is None:
        pytest.skip("ILLUMINA_CATALOG names no local copy of Illumina's v5.0.0 catalog")
    with as_file(ADDED) as added:
        assert cut_catalog(Path(path), LOCI.values(), added=added) == _catalog()


def test_out_list() -> None:
    """SOURCE.md's "Out" table: the 32 STRchive genes left out, each with a reason, none
    a `[repeat_loci]` row and RFC1 not among them."""
    out = _out_table()
    assert len(out) == 32
    assert all(reason.strip() for reason in out.values())
    assert not set(out) & set(LOCI)
    assert "RFC1" not in out
    # Nor is a locus written as not covered.
    assert not set(out) & {entry.gene for entry in NOT_COVERED}


def test_every_strchive_locus_has_one_state() -> None:
    """Every STRchive v2.28.0 gene is in exactly one of `[repeat_loci]`, RFC1 (graded by
    its motif call), the not-covered statements and SOURCE.md's Out table."""
    path = os.environ.get("STRCHIVE_LOCI")
    if path is None:
        pytest.skip("STRCHIVE_LOCI names no local STRchive-loci.json")
    genes = {e["gene"] for e in json.loads(Path(path).read_text())}
    states = [
        set(LOCI),
        {"RFC1"},
        {entry.gene for entry in NOT_COVERED},
        set(_out_table()),
    ]
    assert sum(len(s) for s in states) == len(genes) == 79
    assert set().union(*states) == genes


def _out_table() -> dict[str, str]:
    text = SOURCE_MD.read_text()
    section = text[text.index("## Out") :]
    rows = [
        line.split("|")[1:-1]
        for line in section.splitlines()
        if line.startswith("| ") and not line.startswith("| Gene")
    ]
    return {gene.strip(): reason.strip() for gene, reason in rows}


def test_read_calls_hg002() -> None:
    calls = read_calls(VCF, JSON)
    # Slice 3e's run on every locus; the eleven's calls are slice 3a's run's, unchanged.
    for gene in ELEVEN:
        assert LOCI[gene].variant_id in calls, gene
    htt = calls[LOCI["HTT"].variant_id]
    assert [a.repeats for a in htt.alleles] == [17, 24]
    assert all(15 <= a.repeats <= 30 for a in htt.alleles)
    assert htt.alleles[1].ci == (24, 24)
    assert (
        htt.alleles[0].spanning,
        htt.alleles[0].flanking,
        htt.alleles[0].inrepeat,
    ) == (4, 4, 0)
    assert (htt.coverage, htt.read_length, htt.low_depth) == (19.864865, 150, False)
    ar = calls[LOCI["AR"].variant_id]
    assert len(ar.alleles) == 1
    assert (ar.alleles[0].repeats, ar.alleles[0].spanning, ar.alleles[0].flanking) == (
        27,
        3,
        6,
    )
    assert calls[LOCI["ATXN2"].variant_id].alleles[1].ci == (26, 64)
    assert len(calls[LOCI["C9orf72"].variant_id].alleles) == 2


def test_read_calls_hg002_new_loci() -> None:
    calls = read_calls(VCF, JSON)
    assert {locus.variant_id for locus in LOCI.values()} <= set(calls)
    fmr1 = calls[LOCI["FMR1"].variant_id]
    # GIAB's HG002 v5.0q has about 31-32 CGG.
    assert len(fmr1.alleles) == 1 and 29 <= fmr1.alleles[0].repeats <= 34
    fgf14 = [a.repeats for a in calls[LOCI["FGF14"].variant_id].alleles]
    assert len(fgf14) == 2 and fgf14[0] != fgf14[1] and max(fgf14) <= 20
    for variant in (CFTR_TG, CFTR_T):
        assert calls[variant].alleles
        assert all(a.spanning > 0 for a in calls[variant].alleles)
    shown = [locus.variant_id for locus in LOCI.values()] + [
        RFC1_LOCUS,
        CFTR_TG,
        CFTR_T,
    ]
    assert all(calls[v].longest_flanking > 0 for v in shown)


def test_read_calls_low_depth_and_missing_count(tmp_path: Path) -> None:
    lines = VCF.read_text().splitlines()
    header = [l for l in lines if l.startswith("#")]
    htt = next(l for l in lines if "VARID=HTT;" in l)
    atxn1 = next(l for l in lines if "VARID=ATXN1;" in l)
    low = htt.replace("\tPASS\t", "\tLowDepth\t", 1)
    missing = atxn1.replace("30/31:30-30/31-31", "./.:./.")
    vcf = tmp_path / "eh.vcf"
    vcf.write_text("\n".join([*header, low, missing]) + "\n")
    calls = read_calls(vcf, JSON)
    assert set(calls) == {LOCI["HTT"].variant_id}
    assert calls[LOCI["HTT"].variant_id].low_depth


def test_read_calls_longest_flanking(tmp_path: Path) -> None:
    data = json.loads(JSON.read_text())
    variants = data["LocusResults"]["AR"]["Variants"]["AR"]
    variants["CountsOfFlankingReads"] = "(1, 1), (6, 1), (25, 1)"
    htt = data["LocusResults"]["HTT"]["Variants"]
    for variant in htt.values():
        variant["CountsOfFlankingReads"] = "()"
    found = tmp_path / "eh.json"
    found.write_text(json.dumps(data))
    calls = read_calls(VCF, found)
    assert calls["AR"].longest_flanking == 25
    assert calls["HTT"].longest_flanking == 0
    # A unit count listed with a count of zero reads holds no read.
    variants["CountsOfFlankingReads"] = "(1, 1), (30, 0)"
    found.write_text(json.dumps(data))
    assert read_calls(VCF, found)["AR"].longest_flanking == 1


def test_read_calls_inrepeat(tmp_path: Path) -> None:
    """A variant's in-repeat reads are the total of the JSON's `CountsOfInrepeatReads`,
    which counts reads `ADIR` leaves out (it counts only reads consistent with an allele at
    least a read long)."""
    assert all(c.inrepeat == 0 for c in read_calls(VCF, JSON).values())
    data = json.loads(JSON.read_text())
    data["LocusResults"]["CNBP"]["Variants"]["CNBP"]["CountsOfInrepeatReads"] = (
        "(37, 2), (40, 1)"
    )
    found = tmp_path / "eh.json"
    found.write_text(json.dumps(data))
    calls = read_calls(VCF, found)
    assert calls["CNBP"].inrepeat == 3
    assert all(a.inrepeat == 0 for a in calls["CNBP"].alleles)


class _Fake:
    """Stands for `run_container` and `run_tool`: records each call and writes the output
    ExpansionHunter and REViewer wrote on HG002 into the mounted output directory."""

    def __init__(self) -> None:
        self.containers: list[
            tuple[str, list[str], list[tuple[Path, str, bool]], str | None]
        ] = []
        self.tools: list[list[str]] = []

    def container(
        self,
        name: str,
        args: list[str],
        *,
        mounts: list[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        self.containers.append((name, list(args), list(mounts), label))
        (out,) = (host for host, target, ro in mounts if target == "/out" and not ro)
        if name == "expansionhunter":
            shutil.copy(VCF, out / "eh.vcf")
            shutil.copy(JSON, out / "eh.json")
            (out / "eh_realigned.bam").write_bytes(b"")
        else:
            locus = args[args.index("--locus") + 1]
            for locus_id in locus.split(","):
                shutil.copy(SVG, out / f"rev.{locus_id}.svg")

    def tool(self, args: list[str], **_: Any) -> str:
        self.tools.append(list(args))
        if args[:2] == ["samtools", "sort"]:
            Path(args[args.index("-o") + 1]).write_bytes(b"")
        return ""


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> _Fake:
    stub = _Fake()
    monkeypatch.setattr(expansionhunter, "run_container", stub.container)
    monkeypatch.setattr(expansionhunter, "run_tool", stub.tool)
    return stub


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
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


def test_run_expansionhunter_args(tmp_path: Path, fake: _Fake) -> None:
    cram, fasta, out = _inputs(tmp_path)
    vcf = run_expansionhunter(cram, fasta, out, sex="male", threads=3)
    assert vcf == out / "eh.vcf" and vcf.is_file()
    ((name, args, mounts, label),) = fake.containers
    assert (name, label) == ("expansionhunter", "dna:repeats")
    assert args[0] == "ExpansionHunter"
    assert args[args.index("--sex") + 1] == "male"
    assert args[args.index("--analysis-mode") + 1] == "seeking"
    assert "-m" not in args
    assert args[args.index("--threads") + 1] == "3"
    assert args[args.index("--reads") + 1] == "/cram/x.cram"
    assert args[args.index("--reference") + 1] == "/ref/g.fa"
    assert args[args.index("--output-prefix") + 1] == "/out/eh"
    catalog = next(host for host, target, _ in mounts if target == "/catalog.json")
    assert args[args.index("--variant-catalog") + 1] == "/catalog.json"
    assert catalog.read_text() == CATALOG.read_text()
    assert {(t, ro) for _, t, ro in mounts} == {
        ("/cram", True),
        ("/ref", True),
        ("/catalog.json", True),
        ("/out", False),
    }
    assert [t[:2] for t in fake.tools] == [["samtools", "sort"], ["samtools", "index"]]
    assert fake.tools[0][-3:] == [
        "-o",
        str(out / "eh_realigned.sorted.bam"),
        str(out / "eh_realigned.bam"),
    ]
    assert fake.tools[1][-1] == str(out / "eh_realigned.sorted.bam")


def test_review_args(tmp_path: Path, fake: _Fake) -> None:
    _, fasta, out = _inputs(tmp_path)
    vcf = out / "eh.vcf"
    shutil.copy(VCF, vcf)
    images = review(out, vcf=vcf, fasta=fasta, locus_ids=["HTT", "AR"])
    ((name, args, mounts, label),) = fake.containers
    assert (name, label) == ("reviewer", "dna:repeats")
    assert args[args.index("--locus") + 1] == "HTT,AR"
    assert args[args.index("--region-extension-length") + 1] == "1000"
    assert args[args.index("--reads") + 1] == "/out/eh_realigned.sorted.bam"
    assert args[args.index("--vcf") + 1] == "/out/eh.vcf"
    assert args[args.index("--output-prefix") + 1] == "/out/rev"
    assert (fasta.parent, "/ref", True) in mounts
    assert (out, "/out", False) in mounts
    assert images == {
        "HTT": out / "rev.HTT.svg",
        "AR": out / "rev.AR.svg",
    }


def test_flank_ratio(
    tmp_path: Path, fake: _Fake, monkeypatch: pytest.MonkeyPatch
) -> None:
    cram, fasta, _ = _inputs(tmp_path)
    beds: list[str] = []

    def bedcov(args: list[str], **_: Any) -> str:
        fake.tools.append(list(args))
        bed = Path(args[-2]).read_text()
        beds.append(bed)
        return "chr9\t27572983\t27573483\t10000\nchr9\t27573546\t27574046\t14000\n"

    monkeypatch.setattr(expansionhunter, "run_tool", bedcov)
    ratio = flank_ratio(
        cram, fasta, region=("chr9", 27573484, 27573546), mean_depth=30.0
    )
    assert ratio == pytest.approx(24.0 / 30.0)
    assert beds == ["chr9\t27572983\t27573483\nchr9\t27573546\t27574046\n"]
    assert fake.tools[0][:4] == ["samtools", "bedcov", "--reference", str(fasta)]
    assert fake.tools[0][-1] == str(cram)
