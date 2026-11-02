"""Read-backed phasing of candidate compound heterozygotes: the pairs the promotion would make,
phased by WhatsHap at those sites and the run's heterozygous SNVs between them, and the
promotion's three outcomes.

The sites handed to WhatsHap are checked with a stub that writes its output in the shape
WhatsHap 2.8 wrote on synthetic reads (run on 2026-10-20): a pair in trans as `0|1` and `1|0`
under one PS, a pair in cis as `0|1` twice, and a site no read joins to another left `0/1`
without PS. The phases themselves are checked with the pinned image on synthetic reads
(a stub phases a pair whether or not any read could join it)."""

import json
import random
from collections.abc import Sequence
from pathlib import Path

import pysam
import pytest

from genokeeper import phasing
from genokeeper.annotate.vep import fasta_contigs
from genokeeper.containers import PLATFORM
from genokeeper.evidence.predicted import (
    FAILED_GATE,
    GATE_TRANSCRIPT,
    GENETICIST_NOTE,
    ONE_RECESSIVE_COPY,
    PREDICTED_REASON,
)
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import Called
from genokeeper.phasing import (
    IN_CIS_REASON,
    IN_TRANS_REASON,
    Phase,
    phase_candidates,
    phase_pairs,
)
from genokeeper.sources import Collected, Row
from genokeeper.sources.promote import (
    COMPOUND_HET_REASON,
    candidate_pairs,
    promote_compound_hets,
)
from genokeeper.store.db import Kind, RowType, Topic

type Allele = tuple[str, int, str, str]

AR_GENE = "GENEAR"  # recessive by CGD, on neither list
# How one het in a carrier gene (Tier 3 or 4: CFTR, MUTYH, HFE) is graded.
CARRIER = {"kind": Kind.FINDING, "topic": Topic.CARRIER}
A: Allele = ("chr1", 1000, "A", "C")
B: Allele = ("chr1", 1150, "C", "G")
C: Allele = ("chr1", 3000, "T", "A")


def _vid(allele: Allele) -> str:
    return "GRCh38:" + ":".join(map(str, allele))


def _row(
    n: int,
    allele: Allele,
    *,
    gene: str = AR_GENE,
    dosage: str = "1",
    cls: str = "Pathogenic",
    kind: Kind = Kind.NON_FINDING,
    topic: Topic | None = None,
    source: str = "clinvar",
    reason: str | None = None,
) -> Row:
    vid = _vid(allele)
    return Row(
        variant={"run_id": "r", "vid": vid, "rsid": None},
        call={"run_id": "r", "vid": vid, "dosage": dosage},
        statement={
            "statement_id": f"{source}:{n}",
            "kind": str(kind),
            "source": source,
            "vid": vid,
            "gene": gene,
            "topic": None if topic is None else str(topic),
            "reason": reason,
            "claim": json.dumps({"classification": cls}),
            "actionable": False,
        },
    )


def _promote(rows: list[Row], phases: dict[tuple[Allele, Allele], Phase] | None):
    collected = [Collected("clinvar", rows, {}, len(rows), {})]
    (out,) = promote_compound_hets(
        collected,
        load_gene_lists(),
        {AR_GENE: "AR", "GENEAD": "AD"},
        lambda chrom, pos: chrom == "chrX",
        phases=phases,
    )
    return {r.statement["statement_id"]: r.statement for r in out.rows}


def test_in_trans_pair_promoted() -> None:
    found = _promote([_row(1, A), _row(2, B)], {(A, B): Phase.IN_TRANS})
    for statement in found.values():
        assert statement["kind"] == Kind.FINDING
        assert statement["topic"] == Topic.RISKS
        assert statement["reason"] == IN_TRANS_REASON


def test_in_cis_pair_not_a_finding() -> None:
    found = _promote([_row(1, A), _row(2, B)], {(A, B): Phase.IN_CIS})
    for statement in found.values():
        assert statement["kind"] == Kind.NON_FINDING
        assert statement["topic"] is None
        assert statement["reason"] == IN_CIS_REASON


def test_unphased_keeps_reason() -> None:
    for phases in ({(A, B): Phase.UNPHASED}, {}, None):
        found = _promote([_row(1, A), _row(2, B)], phases)
        assert {s["kind"] for s in found.values()} == {Kind.FINDING}
        assert {s["reason"] for s in found.values()} == {COMPOUND_HET_REASON}
    # Three hets: in trans with one partner is enough; a site in cis with every partner
    # it was phased against, and unphased with none, is no finding.
    found = _promote(
        [_row(1, A), _row(2, B), _row(3, C)],
        {(A, B): Phase.IN_TRANS, (A, C): Phase.IN_CIS, (B, C): Phase.UNPHASED},
    )
    assert found["clinvar:1"]["reason"] == IN_TRANS_REASON
    assert found["clinvar:2"]["reason"] == IN_TRANS_REASON
    assert found["clinvar:3"]["reason"] == COMPOUND_HET_REASON
    found = _promote(
        [_row(1, A), _row(2, B), _row(3, C)],
        {(A, B): Phase.IN_CIS, (A, C): Phase.IN_CIS, (B, C): Phase.IN_TRANS},
    )
    assert found["clinvar:1"]["kind"] == Kind.NON_FINDING
    assert found["clinvar:1"]["reason"] == IN_CIS_REASON


def test_sf_recessive_pair_actionable_unless_in_cis() -> None:
    """ACMG SF v3.3 reports two P/LP variants in a recessive SF gene "even when phase is
    undetermined" (Table 1 footnote a), never two on one haplotype; a pair off the list, a
    predicted het and an HFE pair (p.C282Y homozygotes only) are not reported."""
    mutyh = [_row(1, A, gene="MUTYH", **CARRIER), _row(2, B, gene="MUTYH", **CARRIER)]
    for phases, actionable in (
        ({(A, B): Phase.IN_TRANS}, True),
        ({(A, B): Phase.UNPHASED}, True),
        ({(A, B): Phase.IN_CIS}, False),
    ):
        found = _promote(mutyh, phases)
        assert {s["actionable"] for s in found.values()} == {actionable}, phases
    for gene, single in ((AR_GENE, {}), ("HFE", CARRIER)):
        pair = [_row(1, A, gene=gene, **single), _row(2, B, gene=gene, **single)]
        found = _promote(pair, {(A, B): Phase.IN_TRANS})
        assert {s["kind"] for s in found.values()} == {Kind.FINDING}
        assert {s["actionable"] for s in found.values()} == {False}, gene
    found = _promote(
        [_row(1, A, gene="MUTYH", **CARRIER), _predicted(2, B, "MUTYH", Topic.CARRIER)],
        {(A, B): Phase.IN_TRANS},
    )
    assert (found["clinvar:1"]["actionable"], found["vep:2"]["actionable"]) == (
        True,
        False,
    )


def _predicted(
    n: int, allele: Allele, gene: str, topic: Topic | None = None, why: str = ""
) -> Row:
    """A predicted het as `sources.predicted` stores it: a one-copy non-finding, or a
    finding under `topic` for `why`."""
    why = why or f"missense; {ONE_RECESSIVE_COPY}"
    reason = f"{PREDICTED_REASON}; {why}; {GENETICIST_NOTE}"
    kind = Kind.NON_FINDING if topic is None else Kind.FINDING
    return _row(
        n, allele, gene=gene, cls="predicted, not classified", kind=kind, topic=topic,
        source="vep", reason=reason,
    )  # fmt: skip


def test_tier3_and_predicted_pairs_promoted() -> None:
    """Any two heterozygous pathogenic-class variants in one recessive gene are a pair
    (design decision, 2026-10-28): two classified carrier findings in a Tier 3 gene (CFTR),
    a classified and a predicted het (MUTYH), two predicted hets (HFE), each a carrier
    finding in a Tier 4 gene. A
    research statement that failed a LoF gate is not pathogenic-class and joins none. In cis each keeps
    its single-het grading, a Tier 3 carrier finding included; predicted rows keep their
    own reason before the phase."""
    cftr: tuple[Allele, Allele] = (("chr7", 100, "A", "G"), ("chr7", 200, "C", "T"))
    mutyh: tuple[Allele, Allele] = (("chr1", 100, "A", "G"), ("chr1", 200, "C", "T"))
    hfe: tuple[Allele, Allele] = (("chr6", 100, "A", "G"), ("chr6", 200, "C", "T"))
    carrier = Kind.FINDING, Topic.CARRIER
    failed = f"{FAILED_GATE} {GATE_TRANSCRIPT}"
    rows = [
        _row(1, cftr[0], gene="CFTR", kind=carrier[0], topic=carrier[1]),
        _row(2, cftr[1], gene="CFTR", kind=carrier[0], topic=carrier[1]),
        _row(3, mutyh[0], gene="MUTYH", **CARRIER),
        _predicted(4, mutyh[1], "MUTYH", Topic.CARRIER, "missense"),
        _predicted(5, hfe[0], "HFE", Topic.CARRIER, "missense"),
        _predicted(6, hfe[1], "HFE", Topic.CARRIER, "missense"),
        _predicted(7, ("chr6", 300, "G", "A"), "HFE", Topic.RESEARCH, failed),
    ]
    collected = [
        Collected("clinvar", rows[:3], {}, 3, {}),
        Collected("vep", rows[3:], {}, 4, {}),
    ]
    lists, haploid = load_gene_lists(), lambda chrom, pos: False
    pairs = candidate_pairs(collected, lists, {}, haploid)
    assert pairs == [mutyh, hfe, cftr]
    phases = {cftr: Phase.IN_CIS, mutyh: Phase.IN_TRANS}
    out = promote_compound_hets(collected, lists, {}, haploid, phases)
    found = {r.statement["statement_id"]: r.statement for c in out for r in c.rows}
    for n in (1, 2):
        cis = found[f"clinvar:{n}"]
        assert (cis["kind"], cis["topic"], cis["reason"]) == (
            Kind.FINDING,
            Topic.CARRIER,
            IN_CIS_REASON,
        )
    assert found["clinvar:3"]["reason"] == IN_TRANS_REASON
    for n, why in (
        (4, IN_TRANS_REASON),
        (5, COMPOUND_HET_REASON),
        (6, COMPOUND_HET_REASON),
    ):
        row = found[f"vep:{n}"]
        assert (row["kind"], row["topic"]) == (Kind.FINDING, Topic.RISKS)
        assert str(row["reason"]).startswith(PREDICTED_REASON)
        assert str(row["reason"]).endswith(f"{GENETICIST_NOTE}; {why}")
    assert found["vep:7"]["topic"] == Topic.RESEARCH


# WhatsHap 2.8's output for three sites joined by reads and one alone (FORMAT GT:PS).
WHATSHAP_OUT = {
    A: "0|1:1000",
    B: "1|0:1000",
    C: "0|1:3000",
    ("chr1", 3150, "C", "G"): "0|1:3000",
    ("chr1", 5000, "T", "A"): "0/1",
    # The run's heterozygous SNVs between a pair's sites (the chain).
    ("chr1", 1100, "G", "T"): "1|0:1000",
    ("chr1", 2500, "G", "T"): "0/1",
    ("chr1", 4000, "G", "T"): "0/1",
}


class FakeWhatsHap:
    """Stands in for `run_container` in `genokeeper.phasing`: keeps the sites it was asked
    to phase and writes WHATSHAP_OUT's genotype at each."""

    def __init__(self) -> None:
        self.sites: list[Allele] = []
        self.args: list[str] = []
        self.mounts: dict[str, tuple[Path, bool]] = {}

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, platform, label) == ("whatshap", PLATFORM, "dna:phase")
        self.args = list(args)
        self.mounts = {target: (host, ro) for host, target, ro in mounts}
        work = self.mounts["/out"][0]
        source = work / Path(args[-2]).name
        out = work / Path(args[args.index("-o") + 1]).name
        with pysam.VariantFile(str(source)) as vcf:
            header = str(vcf.header).replace(
                "#CHROM",
                '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n#CHROM',
            )
            lines = []
            for rec in vcf:
                site = (rec.chrom, rec.pos, rec.ref or "", (rec.alts or ("",))[0])
                self.sites.append(site)
                gt = WHATSHAP_OUT[site]
                fmt = "GT:PS" if ":" in gt else "GT"
                lines.append(
                    f"{rec.chrom}\t{rec.pos}\t.\t{site[2]}\t{site[3]}\t.\tPASS\t.\t{fmt}\t{gt}"
                )
        out.write_text(header + "\n".join(lines) + "\n")


@pytest.mark.tools
def test_candidate_sites_and_their_chain_phased(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The candidates: pairs of het pathogenic sites in one recessive gene, nothing else
    # (Tier 3 and predicted pairs: `test_tier3_and_predicted_pairs_promoted`).
    lone, dominant = ("chr1", 5000, "T", "A"), ("chr1", 7000, "G", "C")
    hom, benign = ("chr1", 7200, "A", "G"), ("chr1", 7300, "A", "T")
    x1, x2 = ("chrX", 100, "A", "G"), ("chrX", 200, "A", "G")
    rows = [
        _row(1, A),
        _row(2, B),
        _row(3, C),
        _row(4, ("chr1", 3150, "C", "G")),
        _row(5, dominant, gene="GENEAD"),
        _row(6, ("chr1", 7050, "G", "C"), gene="GENEAD"),
        _row(9, hom, gene="GENEAR2", dosage="2"),
        _row(10, benign, gene="GENEAR2", cls="Benign"),
        _row(11, x1, gene="GENEAR3"),
        _row(12, x2, gene="GENEAR3"),
        _row(13, lone, gene="GENEAR4"),
    ]
    collected = [Collected("clinvar", rows, {}, len(rows), {})]
    pairs = candidate_pairs(
        collected,
        load_gene_lists(),
        {
            AR_GENE: "AR",
            "GENEAD": "AD",
            "GENEAR2": "AR",
            "GENEAR3": "AR",
            "GENEAR4": "AR",
        },
        lambda chrom, pos: chrom == "chrX",
    )
    d = ("chr1", 3150, "C", "G")
    assert pairs == [(A, B), (A, C), (A, d), (B, C), (B, d), (C, d)]

    fake = FakeWhatsHap()
    monkeypatch.setattr(phasing, "run_container", fake)
    fasta = tmp_path / "ref" / "ref.fa"
    fasta.parent.mkdir()
    fasta.write_text(">chr1\n" + "A" * 6000 + "\n")
    pysam.faidx(str(fasta))
    cram = tmp_path / "reads" / "sample.cram"
    cram.parent.mkdir()
    cram.write_text("placeholder\n")
    asked = [*pairs, (C, ("chr1", 5000, "T", "A"))]
    chain = [
        ("chr1", 1100, "G", "T"),
        ("chr1", 2500, "G", "T"),
        ("chr1", 4000, "G", "T"),
    ]
    calls = _calls(
        tmp_path / "sample.vcf",
        fasta,
        [
            (("chr1", 1000, "A", "T"), "PASS", "0/1"),  # at a candidate's position
            (chain[0], "PASS", "0/1"),
            (("chr1", 2000, "G", "T"), "PASS", "1/1"),
            (("chr1", 2100, "G", "T"), "RefCall", "0/1"),
            (("chr1", 2200, "GA", "G"), "PASS", "0/1"),  # an indel
            (("chr1", 2300, "G", "T,C"), "PASS", "1/2"),
            (chain[1], "PASS", "0|1"),
            (chain[2], "PASS", "0/1"),
            (("chr1", 5800, "G", "T"), "PASS", "0/1"),  # past every pair
        ],
    )

    phases = phase_pairs(cram, fasta, calls, asked)

    # The candidate sites and the heterozygous PASS SNVs between a pair's sites, each
    # once, in order.
    assert fake.sites == [
        A,
        chain[0],
        B,
        chain[1],
        C,
        d,
        chain[2],
        ("chr1", 5000, "T", "A"),
    ]
    assert phases[(A, B)] == Phase.IN_TRANS
    assert phases[(C, d)] == Phase.IN_CIS
    assert phases[(A, C)] == Phase.UNPHASED  # two phase sets
    assert phases[(C, ("chr1", 5000, "T", "A"))] == Phase.UNPHASED
    assert set(phases) == set(asked)
    assert fake.mounts["/ref"] == (fasta.parent, True)
    assert fake.mounts["/reads"] == (cram.parent, True)
    assert "--ignore-read-groups" in fake.args
    assert fake.args[fake.args.index("--reference") + 1] == "/ref/ref.fa"
    assert fake.args[-1] == "/reads/sample.cram"
    # The work directory is removed.
    assert sorted(p.name for p in cram.parent.iterdir()) == ["sample.cram"]


# Real WhatsHap on synthetic reads. A random chr19 of 14 kb: heterozygous SNVs every
# 300 bp from 1,000 to 10,000, then two more at 12,000 and 12,300 past a 2 kb gap no fragment
# spans; site i's ALT is on haplotype i % 2. 150 bp read pairs with inserts of 350-450 bp,
# one every 15 bp, alternating haplotypes, all at MAPQ 60.
SITES = [1000 + 300 * i for i in range(31)] + [12000, 12300]


def _synthetic(tmp_path: Path) -> tuple[Path, Path, list[Allele]]:
    """The reference, the CRAM of the reads and the sites' alleles, in SITES' order."""
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(14000)]
    haps = [seq.copy(), seq.copy()]
    alleles: list[Allele] = []
    for i, pos in enumerate(SITES):
        ref = seq[pos - 1]
        alt = "ACGT"[("ACGT".index(ref) + 1) % 4]
        haps[i % 2][pos - 1] = alt
        alleles.append(("chr19", pos, ref, alt))
    fasta = tmp_path / "ref" / "ref.fa"
    fasta.parent.mkdir()
    fasta.write_text(">chr19\n" + "".join(seq) + "\n")
    pysam.faidx(str(fasta))
    cram = tmp_path / "reads" / "sample.cram"
    cram.parent.mkdir()
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr19", "LN": 14000}],
    }
    reads: list[pysam.AlignedSegment] = []
    for n, start in enumerate(range(0, 14000 - 450, 15)):
        hap, insert = haps[n % 2], rng.randint(350, 450)
        mate = start + insert - 150
        for first, at in ((True, start), (False, mate)):
            read = pysam.AlignedSegment()
            read.query_name = f"f{n}"
            read.flag = 99 if first else 147
            read.reference_id = read.next_reference_id = 0
            read.reference_start = at
            read.next_reference_start = mate if first else start
            read.template_length = insert if first else -insert
            read.mapping_quality, read.cigarstring = 60, "150M"
            read.query_sequence = "".join(hap[at : at + 150])
            read.query_qualities = pysam.qualitystring_to_array("I" * 150)
            reads.append(read)
    with pysam.AlignmentFile(
        str(cram), "wc", header=header, reference_filename=str(fasta)
    ) as out:
        for read in sorted(reads, key=lambda r: r.reference_start):
            out.write(read)
    pysam.index(str(cram))
    return fasta, cram, alleles


def _calls(path: Path, fasta: Path, rows: Sequence[tuple[Allele, str, str]]) -> Path:
    """`rows` (allele, FILTER, GT) as the run's bgzipped, indexed calls."""
    contigs = "".join(
        f"##contig=<ID={c},length={n}>\n" for c, n in fasta_contigs(fasta).items()
    )
    path.write_text(
        '##fileformat=VCFv4.2\n##FILTER=<ID=RefCall,Description="x">\n'
        + contigs
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        + "".join(
            f"{c}\t{p}\t.\t{r}\t{a}\t.\t{f}\t.\tGT\t{gt}\n"
            for (c, p, r, a), f, gt in rows
        )
    )
    return Path(pysam.tabix_index(str(path), preset="vcf", force=True))


@pytest.mark.tools
def test_whatshap_phases_pairs_farther_apart_than_a_fragment(tmp_path: Path) -> None:
    fasta, cram, site = _synthetic(tmp_path)
    calls = _calls(tmp_path / "sample.vcf", fasta, [(s, "PASS", "0/1") for s in site])
    far_cis, far_trans = (site[2], site[30]), (site[2], site[29])
    near, gap = (site[4], site[5]), (site[30], site[31])

    phases = phase_pairs(cram, fasta, calls, [far_cis, far_trans, near, gap])

    assert phases == {
        far_cis: Phase.IN_CIS,  # 8.4 kb apart, one haplotype
        far_trans: Phase.IN_TRANS,
        near: Phase.IN_TRANS,  # one fragment spans both
        gap: Phase.UNPHASED,  # one haplotype, but no fragment crosses the gap
    }


# --- phase sets for PharmCAT (slice 3f) -----------------------------------------------------

HG002_PHASED = (
    Path(__file__).parent / "fixtures" / "hg002" / "hg002.pharmcat_phased.vcf.txt"
)


class FixtureWhatsHap:
    """Stands in for `run_container` in `genokeeper.phasing`: keeps the sites VCF it was
    handed and writes HG002's WhatsHap 2.8 output (`HG002_PHASED`) as its result."""

    def __init__(self) -> None:
        self.sites: list[Allele] = []

    def __call__(
        self,
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        assert (name, platform, label) == ("whatshap", PLATFORM, "dna:phase")
        work = {target: host for host, target, _ in mounts}["/out"]
        with pysam.VariantFile(str(work / Path(args[-2]).name)) as vcf:
            self.sites = [
                (r.chrom, r.pos, r.ref or "", (r.alts or ("",))[0]) for r in vcf
            ]
        out = work / Path(args[args.index("-o") + 1]).name
        out.write_text(HG002_PHASED.read_text())


def _hg002_sites() -> list[Allele]:
    with pysam.VariantFile(str(HG002_PHASED)) as vcf:
        return [(r.chrom, r.pos, r.ref or "", (r.alts or ("",))[0]) for r in vcf]


def test_phase_sites_reads_hg002_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """HG002's heterozygous SNVs in PharmCAT's gene spans, phased by WhatsHap on its
    PCR-free CRAM: `phase_sites` reads each phased site's set and haplotype; the sites VCF
    it hands WhatsHap holds the run's heterozygous PASS SNVs inside the regions as the
    chain, and none outside them."""
    fake = FixtureWhatsHap()
    monkeypatch.setattr(phasing, "run_container", fake)
    sites = _hg002_sites()
    lengths = {"chr1": 248956422, "chr2": 242193529, "chr6": 170805979}
    others = [f"chr{n}" for n in (*range(3, 6), *range(7, 23), "X", "Y")]
    fasta = tmp_path / "ref" / "ref.fa"
    fasta.parent.mkdir()
    fasta.write_text("")
    fasta.with_name("ref.fa.fai").write_text(
        "".join(
            f"{c}\t{lengths.get(c, 200_000_000)}\t0\t60\t61\n"
            for c in ["chr1", "chr2", *others[:3], "chr6", *others[3:]]
        )
    )
    cram = tmp_path / "reads" / "sample.cram"
    cram.parent.mkdir()
    cram.write_text("placeholder\n")
    inside = ("chr1", 97100000, "G", "T")  # in DPYD's span, not one of the sites
    calls = _calls(
        tmp_path / "calls.vcf",
        fasta,
        [
            (inside, "PASS", "0/1"),
            (("chr1", 97100100, "G", "T"), "PASS", "1/1"),
            (("chr1", 97100200, "G", "T"), "RefCall", "0/1"),
            (("chr1", 97100300, "GA", "G"), "PASS", "0/1"),
            (("chr1", 120000000, "G", "T"), "PASS", "0/1"),  # outside every region
            (("chr2", 233760000, "A", "C"), "PASS", "0/1"),  # UGT1A1's span
        ],
    )
    regions = [("chr1", 97077000, 97922000), ("chr2", 233759000, 233773000)]

    phased = phasing.phase_sites(cram, fasta, calls, sites, regions)

    assert set(phased) <= set(sites)
    assert len(phased) == 624
    sets: dict[int, int] = {}
    for ps, haplotype in phased.values():
        assert haplotype in (0, 1)
        sets[ps] = sets.get(ps, 0) + 1
    assert max(sets.values()) >= 2
    assert phased[("chr1", 97119090, "A", "G")] == (97119090, 0)  # `1|0`
    assert phased[("chr1", 97125880, "T", "A")] == (97125880, 1)  # `0|1`
    chain = set(fake.sites) - set(sites)
    assert chain == {inside, ("chr2", 233760000, "A", "C")}
    assert set(sites) <= set(fake.sites)
    # The work directory is removed.
    assert sorted(p.name for p in cram.parent.iterdir()) == ["sample.cram"]


def test_phase_candidates_are_het_snvs() -> None:
    """The run phases only heterozygous biallelic SNVs off chrM and haploid loci: the kind
    slice 3f's calibration measured; an indel stays unphased."""

    def call(dosage: str) -> Called:
        return Called("0/1", dosage, RowType.PASS, 40, 30, None, None, ())

    snv, indel = ("chr1", 10, "A", "G"), ("chr1", 20, "AT", "A")
    hom, x, m = ("chr1", 30, "C", "T"), ("chrX", 5, "A", "G"), ("chrM", 9, "A", "G")
    graded = {snv: call("1"), indel: call("1"), hom: call("2")}
    graded |= {x: call("1"), m: call("1")}
    assert phase_candidates(graded, lambda chrom, _pos: chrom == "chrX") == [snv]
