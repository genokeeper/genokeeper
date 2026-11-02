"""Read-backed phasing of candidate compound heterozygotes with WhatsHap.

A candidate pair is two heterozygous pathogenic-class variants in one recessive gene
(`sources.promote.candidate_pairs`, track design 2a "Phasing"). WhatsHap 2.8 phases those
sites and the run's heterozygous PASS SNVs between each pair's two sites: a VCF of them, each
`0/1`, and the CRAM's reads, in its pinned container (`containers.IMAGES["whatshap"]`,
offline, behind the core gate as `dna:phase`). A read pair spans a few hundred bases, so two
sites farther apart are joined only through a chain of heterozygous sites between them that
overlapping fragments link. Two sites WhatsHap puts in one phase set (PS) are in cis
when their ALTs sit on the same haplotype (`0|1` and `0|1`) and in trans otherwise (`0|1` and
`1|0`); a site it leaves unphased, or two sites in different phase sets (no chain of reads
joins them), leave the pair unphased. Only the pair's own two sites are read for its phase.
Two candidate alleles at one position are not handed to it and stay unphased.

PharmCAT's heterozygous positions are phased the same way (`phase_sites`, slice 3f "Phase
sets"): the run's heterozygous PASS SNVs inside each PharmCAT gene's span
(`pharmcat.gene_regions`) are the chain, and each phased position's phase set and haplotype
go into PharmCAT's VCF as `GT:PS` when the calibration trusts phase.

A pair in cis keeps its statements at their single-het grading, unpromoted, so the phase is
a control: its anchor comes from the reads in the CRAM that span the pair's sites and the
chain between them, which neither ClinVar's records nor the promotion can move, because
WhatsHap reads each site's allele from those reads alone; a pair no chain of reads joins is
never called in cis.
"""

import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from enum import StrEnum
from pathlib import Path

import pysam

from genokeeper.annotate.vep import fasta_contigs
from genokeeper.containers import PLATFORM, run_container
from genokeeper.genotype.gvcf import Called

type Allele = tuple[str, int, str, str]
type Pair = tuple[Allele, Allele]

IN_TRANS_REASON = "in trans (read-backed)"
IN_CIS_REASON = "in cis (read-backed)"
SAMPLE = "sample"


class Phase(StrEnum):
    """Where WhatsHap puts a pair's two ALTs."""

    IN_TRANS = "in_trans"
    IN_CIS = "in_cis"
    UNPHASED = "unphased"


def _write_sites(path: Path, sites: Sequence[Allele], contigs: dict[str, int]) -> None:
    header = (
        "##fileformat=VCFv4.2\n"
        + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in contigs.items())
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n"
    )
    rows = [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1\n" for c, p, r, a in sites]
    path.write_text(header + "".join(rows))


def _haplotypes(path: Path) -> dict[Allele, tuple[int, int]]:
    """Each phased site's (phase set, the haplotype its ALT is on)."""
    found: dict[Allele, tuple[int, int]] = {}
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            sample = rec.samples[0]
            gt, ps = tuple(sample.get("GT") or ()), sample.get("PS")
            if not sample.phased or ps is None or sorted(gt) != [0, 1]:
                continue
            site = (rec.chrom, rec.pos, (rec.ref or "").upper(), (rec.alts or ("",))[0])
            found[site] = (int(str(ps)), gt.index(1))
    return found


def phase_candidates(
    graded: Mapping[Allele, Called], haploid: Callable[[str, int], bool]
) -> list[Allele]:
    """The run's slice 3f sites to phase, sorted: those called heterozygous that are
    biallelic SNVs (the kind the calibration's phase measure compared) off chrM and the
    `haploid` loci; an indel stays unphased."""
    return sorted(
        a
        for a, called in graded.items()
        if called.dosage == "1"
        and len(a[2]) == len(a[3]) == 1
        and a[0] != "chrM"
        and not haploid(a[0], a[1])
    )


def het_snvs(calls: Path, regions: Sequence[tuple[str, int, int]]) -> set[Allele]:
    """The heterozygous PASS biallelic SNVs in `calls` (bgzipped, indexed) inside each
    region (chrom, start, end; 1-based, inclusive)."""
    found: set[Allele] = set()
    with pysam.VariantFile(str(calls)) as vcf:
        for chrom, start, end in regions:
            if chrom not in vcf.header.contigs:
                continue
            for rec in vcf.fetch(chrom, start - 1, end):
                ref, alts = (rec.ref or "").upper(), rec.alts or ()
                gt = tuple(rec.samples[0].get("GT") or ())
                if (
                    list(rec.filter.keys()) == ["PASS"]
                    and len(alts) == 1
                    and len(ref) == len(alts[0]) == 1
                    and gt in ((0, 1), (1, 0))
                ):
                    found.add((rec.chrom, rec.pos, ref, alts[0].upper()))
    return found


def phase_sites(
    cram: Path,
    fasta: Path,
    calls: Path,
    sites: Sequence[Allele],
    regions: Sequence[tuple[str, int, int]],
) -> dict[Allele, tuple[int, int]]:
    """Each of `sites` WhatsHap phases on `cram` (indexed, on `fasta`): its phase set and the
    haplotype its ALT is on. The chain that links them is every heterozygous PASS
    biallelic SNV of the run's `calls` inside `regions` (chrom, start, end; 1-based,
    inclusive). Two of `sites` at one position, and a chain SNV at a site's position, are
    not handed to WhatsHap."""
    contigs = fasta_contigs(fasta)
    order = {name: i for i, name in enumerate(contigs)}
    positions: dict[tuple[str, int], int] = {}
    for chrom, pos, _, _ in set(sites):
        positions[(chrom, pos)] = positions.get((chrom, pos), 0) + 1
    chain = {s for s in het_snvs(calls, regions) if (s[0], s[1]) not in positions}
    phasing = sorted(
        (
            a
            for a in set(sites) | chain
            if positions.get((a[0], a[1]), 1) == 1 and a[0] in order
        ),
        key=lambda a: (order[a[0]], a[1], a[2], a[3]),
    )
    if not phasing:
        return {}
    work = Path(tempfile.mkdtemp(prefix=".phase-", dir=cram.parent))
    try:
        _write_sites(work / "sites.vcf", phasing, contigs)
        run_container(
            "whatshap",
            [
                "whatshap",
                "phase",
                "--ignore-read-groups",
                "--reference",
                f"/ref/{fasta.name}",
                "-o",
                "/out/phased.vcf",
                "/out/sites.vcf",
                f"/reads/{cram.name}",
            ],
            mounts=[
                (fasta.parent, "/ref", True),
                (cram.parent, "/reads", True),
                (work, "/out", False),
            ],
            platform=PLATFORM,
            label="dna:phase",
        )
        phased = _haplotypes(work / "phased.vcf")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    wanted = set(sites)
    return {site: phase for site, phase in phased.items() if site in wanted}


def phase_pairs(
    cram: Path, fasta: Path, calls: Path, pairs: Sequence[Pair]
) -> dict[Pair, Phase]:
    """Each pair's phase from WhatsHap on `cram` (indexed, on `fasta`) at the pairs' sites
    and the heterozygous PASS SNVs of the run's `calls` between them (`phase_sites`, each
    same-chromosome pair's open interval as a region), keyed by the pair as given."""
    alleles = sorted({a for pair in pairs for a in pair})
    regions = [
        (chrom, min(x, y) + 1, max(x, y) - 1)
        for (chrom, x, _, _), (other, y, _, _) in pairs
        if chrom == other
    ]
    phased = phase_sites(cram, fasta, calls, alleles, regions)
    out: dict[Pair, Phase] = {}
    for a, b in pairs:
        pa, pb = phased.get(a), phased.get(b)
        if pa is None or pb is None or pa[0] != pb[0]:
            out[(a, b)] = Phase.UNPHASED
        else:
            out[(a, b)] = Phase.IN_CIS if pa[1] == pb[1] else Phase.IN_TRANS
    return out
