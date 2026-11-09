"""Ensembl VEP and SpliceAI (`genokeeper.annotate.vep`), their containers faked: VEP's
command, its tab output read into records per variant and gene, and SpliceAI on the listed
genes' variants only."""

import gzip
from collections.abc import Sequence
from pathlib import Path

import pysam
import pytest
from fixtures.gvcf_mini import write_fasta

from genokeeper import pipeline
from genokeeper.annotate import vep
from genokeeper.annotate.vep import (
    VEP_INPUT,
    VEP_OUTPUT,
    VepFiles,
    annotate,
    locus,
    nmd_bed,
    parse_vep,
    vep_files,
)
from genokeeper.containers import IMAGES, PLATFORM
from genokeeper.genelists import LISTED, load_gene_lists

FIELDS = (
    "Uploaded_variation\tGene\tSYMBOL\tFeature\tFeature_type\tConsequence\tMANE_SELECT"
    "\tMANE_PLUS_CLINICAL\tEXON\tCDS_position\tHGVSc\tHGVSp\tREVEL"
    "\tam_pathogenicity\tSpliceVault_top_events"
)
STOP = "GRCh38:chr1:200:C:T"
MISSENSE_TERM = "missense_variant"
MISSENSE = "GRCh38:chr1:300:G:A"
SPLICE = "GRCh38:chr1:400:G:C"
INTRONIC = "GRCh38:chr1:500:A:G"
ELSEWHERE = "GRCh38:chr1:600:T:C"
# VEP 116's tab output with --fields, --mane, --hgvs, --numbers and --total_length: `##`
# lines, the column header, then one row per variant and transcript (`-` for no value).
_ = "-"
ROWS = [
    "\t".join(row)
    for row in (
        # The stop: synonymous on a non-MANE transcript, stop_gained on MANE Select.
        (STOP, "ENSG01", "BRCA1", "ENST02", "Transcript", "synonymous_variant", _, _)
        + ("9/22", "1230/5500", "ENST02.1:c.1230C>T")
        + ("ENSP02.1:p.Gln410=", _, _, _),
        (STOP, "ENSG01", "BRCA1", "ENST01", "Transcript", "stop_gained", "NM_007294.4")
        + (_, "10/23", "1234/5592", "ENST01.9:c.1234C>T")
        + ("ENSP01.3:p.Gln412Ter", _, _, _),
        (MISSENSE, "ENSG01", "BRCA1", "ENST01", "Transcript", "missense_variant")
        + ("NM_007294.4", _, "11/23", "2000/5592", "ENST01.9:c.2000G>A")
        + ("ENSP01.3:p.Arg667His", "0.812", "0.95", _),
        # A splice donor on a non-MANE transcript outranks a MANE intron row: a MANE Select
        # row wins only with a consequence of interest (gate 2 then names the transcript).
        (SPLICE, "ENSG01", "BRCA1", "ENST01", "Transcript", "intron_variant")
        + ("NM_007294.4", _, _, _, _, _, _, _, _),
        (SPLICE, "ENSG01", "BRCA1", "ENST03", "Transcript", "splice_donor_variant")
        + (_, _, _, _, "ENST03.1:c.100+1G>C", _, _, _)
        + ("Top1:CD:-35:31%:Frameshift,Top2:ES:3:2%:inFrame",),
        (INTRONIC, "ENSG01", "BRCA1", "ENST01", "Transcript", "intron_variant")
        + ("NM_007294.4", _, _, _, _, _, _, _, _),
        (ELSEWHERE, "ENSG09", "GENEZ", "ENST09", "Transcript", "intron_variant")
        + (_, _, _, _, _, _, _, _, _),
        (ELSEWHERE, "ENSG09", "GENEZ", "ENSR09", "RegulatoryFeature")
        + ("regulatory_region_variant", _, _, _, _, _, _, _, _, _),
    )
]


def _write_vep_output(path: Path, rows: Sequence[str] = ROWS) -> Path:
    text = (
        "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## Column descriptions:\n"
        f"#{FIELDS}\n" + "".join(f"{r}\n" for r in rows)
    )
    path.write_bytes(gzip.compress(text.encode()))
    return path


def _files(lab: Path) -> VepFiles:
    """`vep_files(lab)` with a pext table: BRCA1's stop base at 0.4 of its gene's 0.8,
    and a base one beside the splice donor at the gene's largest."""
    files = vep_files(lab)
    files.folder.mkdir(parents=True)
    plain = files.pext.with_suffix("")
    plain.write_text("chr1\t199\t200\tENSG01\t0.4\nchr1\t401\t402\tENSG01\t0.8\n")
    pysam.tabix_index(str(plain), preset="bed", force=True)
    files.pext_max.write_text("ENSG01\t0.8\n")
    _revel_file(files.revel, [("300", "G", "A", "R", "H", "0.81", "ENST09")])
    return files


REVEL_HEADER = (
    "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\ttranscripts\n"
)


def _revel_file(path: Path, rows: Sequence[tuple[str, ...]]) -> None:
    """REVEL's file as the bootstrap indexes it, its rows on chr `1` at their GRCh38 base."""
    plain = path.with_suffix("")
    body = "".join(f"1\t{r[0]}\t" + "\t".join(r) + "\n" for r in rows)
    plain.write_text(REVEL_HEADER + body)
    pysam.tabix_index(str(plain), seq_col=0, start_col=2, end_col=2, force=True)


class FakeContainers:
    """Stands in for `run_container` in `genokeeper.annotate.vep`: VEP writes `ROWS`, BRCA1
    named `gene`, SpliceAI scores each site it is handed 0.3 on `spliceai_gene`."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, list[str], str | None]] = []
        self.spliceai_sites: list[str] = []
        self.gene = self.spliceai_gene = "BRCA1"

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        self.calls.append((name, list(args), label))
        work = next(host for host, target, _ in mounts if target == "/work")
        if name == "vep":
            named = [row.replace("\tBRCA1\t", f"\t{self.gene}\t") for row in ROWS]
            _write_vep_output(work / VEP_OUTPUT, named)
            return
        assert name == "spliceai"
        given, written = (args[args.index(flag) + 1] for flag in ("-I", "-O"))
        lines = (work / Path(given).name).read_text().splitlines()
        out: list[str] = []
        for line in lines:
            if line.startswith("#"):
                out.append(line)
                continue
            fields = line.split("\t")
            self.spliceai_sites.append(fields[2])
            info = (
                f"SpliceAI={fields[4]}|{self.spliceai_gene}|0.00|0.30|0.10|0.00|1|2|3|4"
            )
            out.append("\t".join([*fields[:7], info]))
        (work / Path(written).name).write_text("\n".join(out) + "\n")


def _sites(path: Path, vids: Sequence[str]) -> Path:
    rows = []
    for vid in vids:
        chrom, pos, ref, alt = locus(vid)
        rows.append(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(rows)
    )
    return path


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> FakeContainers:
    found = FakeContainers()
    monkeypatch.setattr(vep, "run_container", found)
    return found


def _annotate(
    tmp_path: Path, genes: Sequence[str], *, threads: int = 4
) -> list[vep.VepRecord]:
    fasta = write_fasta(tmp_path / "ref" / "GRCh38_no_alt.fa")
    out = tmp_path / "work"
    out.mkdir()
    _sites(out / VEP_INPUT, [STOP, MISSENSE, SPLICE, INTRONIC, ELSEWHERE])
    return annotate(
        out,
        genes_for_spliceai=genes,
        files=_files(tmp_path / "lab"),
        fasta=fasta,
        threads=threads,
    )


def test_vep_args_offline_mane(fake: FakeContainers, tmp_path: Path) -> None:
    """VEP runs offline on the lab's cache with MANE Select and MANE Plus Clinical, HGVS
    from the lab's FASTA, and the REVEL, AlphaMissense and SpliceVault plugins on the lab's
    files, behind the gate as dna:annotate; its image is pinned by digest."""
    _annotate(tmp_path, [])
    (name, args, label), *_ = fake.calls
    assert (name, label) == ("vep", "dna:annotate")
    assert "@sha256:" in IMAGES["vep"]
    for flag in ("--offline", "--cache", "--mane", "--hgvs", "--tab", "--total_length"):
        assert flag in args
    assert "--database" not in args
    assert args[args.index("--fasta") + 1] == "/ref/GRCh38_no_alt.fa"
    assert args[args.index("--dir_cache") + 1] == "/vep/cache"
    plugins = [args[i + 1] for i, a in enumerate(args) if a == "--plugin"]
    assert plugins == [
        "REVEL,file=/vep/revel_v1.3_grch38.tsv.gz",
        "AlphaMissense,file=/vep/AlphaMissense_hg38.tsv.gz",
        "SpliceVault,file=/vep/SpliceVault_data_GRCh38.tsv.gz",
    ]
    assert args[args.index("--fork") + 1] == "4"


def test_vep_output_parsed(tmp_path: Path) -> None:
    """One gene's record per variant: its most severe consequence of interest, MANE first
    among equals; a variant of no interest only in a kept gene."""
    path = _write_vep_output(tmp_path / VEP_OUTPUT)
    records = {r.vid: r for r in parse_vep(path, keep={"BRCA1"})}
    assert list(records) == [STOP, MISSENSE, SPLICE, INTRONIC]
    stop = records[STOP]
    assert (stop.transcript, stop.consequence, stop.mane, stop.lof) == (
        "ENST01",
        "stop_gained",
        True,
        True,
    )
    assert (stop.exon, stop.cds_length) == ("10/23", 5592)
    assert stop.hgvsp == "ENSP01.3:p.Gln412Ter"
    assert (records[MISSENSE].revel, records[MISSENSE].alphamissense) == (0.812, 0.95)
    splice = records[SPLICE]
    assert (splice.transcript, splice.mane, splice.lof) == ("ENST03", False, True)
    assert splice.splicevault_in_frame is True
    assert records[INTRONIC].splicevault_in_frame is None
    assert records[INTRONIC].cds_length is None
    # Without BRCA1 kept, its intronic variant goes.
    assert INTRONIC not in {r.vid for r in parse_vep(path, keep=set())}


def _choice_row(
    vid: str, gene: str, transcript: str, consequence: str, **kw: str
) -> str:
    """One VEP row with only the columns the row choice reads set."""
    mane, plus, revel = kw.get("mane", _), kw.get("plus", _), kw.get("revel", _)
    return "\t".join(
        (
            vid,
            kw.get("gene_id", f"ENSG-{gene}"),
            gene,
            transcript,
            "Transcript",
            consequence,
            mane,
            plus,
        )
        + (_, _, _, _, revel, _, _)
    )


# Each variant's rows; for each, the (gene, transcript) of every record it should give.
MANE_MISSENSE = "GRCh38:chr1:700:C:T"
TWO_GENES_STOP = "GRCh38:chr1:800:C:T"
PLUS_STOP = "GRCh38:chr1:900:C:T"
PLUS_MISSENSE = "GRCh38:chr1:1100:C:T"
SELECT_AND_PLUS = "GRCh38:chr1:1200:C:T"
OVERLAP_PLUS_STOP = "GRCh38:chr1:1300:C:T"
OVERLAP_MANE_STOP = "GRCh38:chr1:1400:C:T"
SHARED_EXON = "GRCh38:chr1:1500:C:T"
SHARED_SYMBOL = "GRCh38:chr1:1600:C:T"
SELECT, PLUS = {"mane": "NM_007294.4"}, {"plus": "NM_000001.1"}
CHOICE_ROWS = [
    _choice_row(MANE_MISSENSE, "BRCA1", "ENST05", "splice_donor_variant"),
    _choice_row(MANE_MISSENSE, "BRCA1", "ENST01", MISSENSE_TERM, **SELECT),
    _choice_row(TWO_GENES_STOP, "GENEZ", "ENST09", "stop_gained"),
    _choice_row(TWO_GENES_STOP, "BRCA1", "ENST06", "stop_gained"),
    _choice_row(PLUS_STOP, "BRCA1", "ENST07", "stop_gained"),
    _choice_row(PLUS_STOP, "BRCA1", "ENST08", "stop_gained", **PLUS),
    _choice_row(PLUS_MISSENSE, "BRCA1", "ENST01", "intron_variant", **SELECT),
    _choice_row(PLUS_MISSENSE, "BRCA1", "ENST02", MISSENSE_TERM, **PLUS),
    _choice_row(PLUS_MISSENSE, "BRCA1", "ENST03", "stop_gained"),
    _choice_row(SELECT_AND_PLUS, "BRCA1", "ENST01", MISSENSE_TERM, **SELECT),
    _choice_row(SELECT_AND_PLUS, "BRCA1", "ENST02", "stop_gained", **PLUS),
    _choice_row(OVERLAP_PLUS_STOP, "GENEZ", "ENST09", MISSENSE_TERM, mane="NM_9.1"),
    _choice_row(OVERLAP_PLUS_STOP, "BRCA1", "ENST02", "stop_gained", **PLUS),
    _choice_row(OVERLAP_MANE_STOP, "GENEZ", "ENST09", "stop_gained", mane="NM_9.1"),
    _choice_row(OVERLAP_MANE_STOP, "BRCA1", "ENST01", MISSENSE_TERM, **SELECT),
    _choice_row(SHARED_EXON, "BRCA1", "ENST02", MISSENSE_TERM, **PLUS),
    _choice_row(SHARED_EXON, "BRCA1", "ENST01", MISSENSE_TERM, **SELECT),
    _choice_row(SHARED_SYMBOL, "BRCA1", "ENST01", MISSENSE_TERM, **SELECT),
    _choice_row(SHARED_SYMBOL, "BRCA1", "ENST31", "stop_gained", gene_id="ENSG-B"),
]


@pytest.mark.parametrize(
    ("vid", "expected"),
    [
        # MANE Select first (design decision, 2026-10-28): its missense is graded, not a
        # minor isoform's splice donor.
        (MANE_MISSENSE, {("BRCA1", "ENST01")}),
        # No row of one gene hides another gene's (BR-013), listed or not.
        (TWO_GENES_STOP, {("GENEZ", "ENST09"), ("BRCA1", "ENST06")}),
        # Equal rank, no MANE Select row: MANE Plus Clinical.
        (PLUS_STOP, {("BRCA1", "ENST08")}),
        # BR-012: a minor isoform's stop no longer hides the MANE Plus Clinical missense.
        (PLUS_MISSENSE, {("BRCA1", "ENST02")}),
        # BR-013: the MANE Select missense and the MANE Plus Clinical stop are both graded.
        (SELECT_AND_PLUS, {("BRCA1", "ENST01"), ("BRCA1", "ENST02")}),
        (OVERLAP_PLUS_STOP, {("GENEZ", "ENST09"), ("BRCA1", "ENST02")}),
        (OVERLAP_MANE_STOP, {("GENEZ", "ENST09"), ("BRCA1", "ENST01")}),
        # A MANE Plus Clinical row of no higher rank adds nothing: one record, on Select.
        (SHARED_EXON, {("BRCA1", "ENST01")}),
        # Two Ensembl genes under one symbol are two genes.
        (SHARED_SYMBOL, {("BRCA1", "ENST01"), ("BRCA1", "ENST31")}),
    ],
)
def test_vep_records_per_gene(
    tmp_path: Path, vid: str, expected: set[tuple[str, str]]
) -> None:
    """Each gene of a variant gives its own records: its MANE Select row with a
    consequence of interest, and a MANE Plus Clinical row of a higher rank; a gene with
    neither gives its row of the highest rank, MANE Plus Clinical first among equals."""
    path = _write_vep_output(tmp_path / VEP_OUTPUT, CHOICE_ROWS)
    records = parse_vep(path, keep={"BRCA1"})
    assert {(r.gene, r.transcript) for r in records if r.vid == vid} == expected


def test_vep_parse_holds_only_rows_it_can_choose() -> None:
    """VEP writes a row per variant and transcript, millions on a genome: a rank-0 row of a
    gene outside `keep` is never held, so memory follows the variants with a graded or
    listed row, and one gene's many transcripts collapse as they stream by into the rows
    `parse_vep` can choose (its MANE Select, its MANE Plus Clinical and its best row)."""
    lines = [
        _choice_row(f"GRCh38:chr1:{pos}:C:T", "GENEZ", "ENST09", "intron_variant")
        for pos in range(2000, 3000)
    ]
    lines += [
        _choice_row(MANE_MISSENSE, "BRCA1", f"ENST{i:02}", "intron_variant")
        for i in range(50)
    ]
    header = FIELDS.split("\t")
    held = vep._candidates(
        (dict(zip(header, line.split("\t"), strict=True)) for line in lines),
        keep={"BRCA1"},
    )
    assert list(held) == [(MANE_MISSENSE, "ENSG-BRCA1")]


def test_spliceai_only_listed_genes(fake: FakeContainers, tmp_path: Path) -> None:
    """SpliceAI scores the records of the listed genes only, in its own image behind the
    gate; pext is read for each loss-of-function record (a splice site's at the nearest
    base of its gene); with no listed record SpliceAI does not run."""
    records = {r.vid: r for r in _annotate(tmp_path, ["BRCA1"])}
    names = [(name, label) for name, _, label in fake.calls]
    assert names == [("vep", "dna:annotate"), ("spliceai", "dna:spliceai")]
    assert fake.spliceai_sites == [STOP, MISSENSE, SPLICE, INTRONIC]
    assert {vid: r.spliceai for vid, r in records.items()} == {
        STOP: 0.3,
        MISSENSE: 0.3,
        SPLICE: 0.3,
        INTRONIC: 0.3,
    }
    assert (records[STOP].pext, records[SPLICE].pext) == (0.5, 1.0)
    assert records[MISSENSE].pext is None
    fake.calls.clear()
    unlisted = _annotate(tmp_path / "unlisted", [])
    assert [name for name, _, _ in fake.calls] == ["vep"]
    assert all(r.spliceai is None for r in unlisted)


@pytest.mark.parametrize(("named", "score"), [("MUT", 0.3), ("GENEZ", None)])
def test_spliceai_score_of_a_renamed_gene(
    fake: FakeContainers, tmp_path: Path, named: str, score: float | None
) -> None:
    """BR-023: MMUT, which SpliceAI's GENCODE v24 annotation still names MUT, takes
    SpliceAI's MUT entry; another gene's entry at the variant is never taken."""
    fake.gene, fake.spliceai_gene = "MMUT", named
    records = _annotate(tmp_path, ["MMUT"])
    assert {r.gene for r in records} == {"MMUT"}
    assert {r.spliceai for r in records} == {score}


def test_spliceai_names_are_listed_genes_and_old_names() -> None:
    """Each renamed entry maps a gene SpliceAI is run on to a name no listed gene bears, so
    no listed gene's own entry is read for another; that the old name is its gene's own
    rests on the derivation `SPLICEAI_NAMES` states."""
    lists = load_gene_lists()
    listed = {gene for name in LISTED for gene in lists[name].genes}
    assert set(vep.SPLICEAI_NAMES) <= listed
    assert not set(vep.SPLICEAI_NAMES.values()) & listed


@pytest.mark.parametrize("threads", [vep.SPLICEAI_CORES * 3, 64])
def test_spliceai_shards_join_to_one_process_output(
    fake: FakeContainers, tmp_path: Path, threads: int
) -> None:
    """With the threads for several SpliceAI processes, the listed records go in order to
    one container each (never more containers than records), and their outputs join into
    the bytes one process writes; the shards are removed."""
    _annotate(tmp_path / "one", ["BRCA1"], threads=1)
    one = (tmp_path / "one" / "work" / vep.SPLICEAI_OUTPUT).read_bytes()
    fake.calls.clear()
    fake.spliceai_sites.clear()
    _annotate(tmp_path / "many", ["BRCA1"], threads=threads)
    work = tmp_path / "many" / "work"
    spliceai = [args for name, args, _ in fake.calls if name == "spliceai"]
    # Four listed records: three shards, or four of one record each.
    assert len(spliceai) == min(threads // vep.SPLICEAI_CORES, 4)
    # The shards run at once, so the order they are handed in is the output's to show.
    assert sorted(fake.spliceai_sites) == sorted([STOP, MISSENSE, SPLICE, INTRONIC])
    assert (work / vep.SPLICEAI_OUTPUT).read_bytes() == one
    assert not (work / vep.SPLICEAI_SHARDS).exists()


def test_nmd_track_converted_offline_in_its_image(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The bootstrap's NMD Escape track goes to BED by bigBedToBed in its pinned amd64
    image, the bigBed mounted read-only."""
    calls: list[tuple[str, list[str], object, object]] = []
    monkeypatch.setattr(
        vep,
        "run_container",
        lambda name, args, *, mounts, platform=None: calls.append(
            (name, list(args), mounts, platform)
        ),
    )
    bigbed = tmp_path / "in" / "nmdEscMane.bb"
    nmd_bed(bigbed, tmp_path / "out" / "nmd.bed")
    assert calls == [
        (
            "ucsc_bigbedtobed",
            ["bigBedToBed", "/in/nmdEscMane.bb", "/out/nmd.bed"],
            [(bigbed.parent, "/in", True), (tmp_path / "out", "/out", False)],
            PLATFORM,
        )
    ]
    assert "@sha256:" in IMAGES["ucsc_bigbedtobed"]


def test_vep_step_key_holds_the_spliceai_names(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The vep step's cache key holds `SPLICEAI_NAMES`: an edit to the map runs VEP and
    SpliceAI again without a step-version bump, and the same map reuses the cached step."""
    monkeypatch.setattr(pipeline, "vep_input", lambda *_, **__: [])
    monkeypatch.setattr(pipeline, "annotate", lambda *_, **__: [])

    def step() -> pipeline._Step:
        return pipeline.vep_step(
            tmp_path / "work",
            vcf=tmp_path / "calls.vcf.gz",
            fasta=tmp_path / "ref.fa",
            files=vep_files(tmp_path),
            listed=["MMUT"],
            threads=1,
            parts=["call", "images"],
        )

    first = step()
    again = step()
    monkeypatch.setitem(vep.SPLICEAI_NAMES, "PLP1", "PLP")
    renamed = step()
    assert (again.key, again.reused) == (first.key, True)
    assert renamed.key != first.key


def _missense_record(
    vid: str, hgvsp: str = "ENSP9.1:p.Ala5Val", **kw: object
) -> vep.VepRecord:
    fields: dict[str, object] = {
        "vid": vid,
        "gene": "GENEZ",
        "gene_id": "ENSG-Z",
        "transcript": "ENST01",
        "consequence": MISSENSE_TERM,
        "mane": True,
        "hgvsc": None,
        "hgvsp": hgvsp,
        "exon": None,
        "cds_length": None,
        "lof": False,
        "revel": None,
        "alphamissense": None,
        "spliceai": None,
        "splicevault_in_frame": None,
        "pext": None,
    }
    return vep.VepRecord(**(fields | kw))  # type: ignore[arg-type]


REVEL_ROWS = [
    ("100", "C", "T", "A", "V", "0.85", "ENST08"),
    # One allele and protein change, scored differently on two transcripts: no single score.
    ("200", "C", "T", "A", "V", "0.90", "ENST08"),
    ("200", "C", "T", "A", "V", "0.30", "ENST09"),
    ("300", "C", "T", "A", "V", "0.50", "ENST08"),
    ("300", "C", "T", "A", "V", "0.50", "ENST09"),
    ("400", "C", "T", "A", "V", "0.85", "ENST08"),
    ("500", "C", "G", "A", "V", "0.85", "ENST08"),
]


@pytest.mark.parametrize(
    ("vid", "hgvsp", "expected"),
    [
        ("GRCh38:chr1:100:C:T", "ENSP9.1:p.Ala5Val", 0.85),
        ("GRCh38:chr1:200:C:T", "ENSP9.1:p.Ala5Val", None),
        ("GRCh38:chr1:300:C:T", "ENSP9.1:p.Ala5Val", 0.50),
        # Another protein change at the same allele is not this one's score.
        ("GRCh38:chr1:400:C:T", "ENSP9.1:p.Ala5Thr", None),
        # Another allele at the base.
        ("GRCh38:chr1:500:C:T", "ENSP9.1:p.Ala5Val", None),
        # A contig REVEL does not hold.
        ("GRCh38:chrUn_KI270742v1:100:C:T", "ENSP9.1:p.Ala5Val", None),
    ],
)
def test_revel_by_allele(
    tmp_path: Path, vid: str, hgvsp: str, expected: float | None
) -> None:
    """BR-042: a missense VEP gave no REVEL takes REVEL's score for its allele and protein
    change on another transcript, and only when every such row agrees."""
    path = tmp_path / "revel.tsv.gz"
    _revel_file(path, REVEL_ROWS)
    record = vep.revel_by_allele([_missense_record(vid, hgvsp)], path)[0]
    assert (record.revel, record.revel_other_transcript) == (
        expected,
        expected is not None,
    )


@pytest.mark.parametrize(
    "kw",
    [{"revel": 0.2}, {"consequence": "stop_gained", "hgvsp": "ENSP9.1:p.Ala5Ter"}],
)
def test_revel_by_allele_leaves_the_others(
    tmp_path: Path, kw: dict[str, object]
) -> None:
    """A score VEP matched on the transcript stays; only a missense is looked up."""
    path = tmp_path / "revel.tsv.gz"
    _revel_file(path, REVEL_ROWS)
    record = _missense_record("GRCh38:chr1:100:C:T", **kw)
    assert vep.revel_by_allele([record], path) == [record]


def test_annotate_fills_revel_by_allele(
    fake: FakeContainers, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`annotate` gives a missense VEP left without REVEL its score by allele."""
    monkeypatch.setitem(
        globals(), "ROWS", [r.replace("\t0.812\t", "\t-\t") for r in ROWS]
    )
    found = {r.vid: r for r in _annotate(tmp_path, [])}[MISSENSE]
    assert (found.revel, found.revel_other_transcript) == (0.81, True)
