"""Predicted findings for variants ClinVar has not classified (track design, slice 1e).

`grade_predicted` grades one VEP record (`annotate.vep.VepRecord`) at the user's call:

- its rule (`predicted_rule`): a loss-of-function consequence (stop-gain, frameshift,
  canonical splice ±1/2), a start-loss, a missense at or above [predicted]'s REVEL threshold
  (AlphaMissense's where REVEL has no score: one calibrated tool per variant), or a SpliceAI
  delta score at or above [predicted] spliceai. Anything else is a non-finding;
- absent from gnomAD, or at or below its gene's BS1 and BA1 (`Grading.bs1_for`, as ClinVar
  claims are graded), else a non-finding;
- a loss of function (or start-loss) counts only in a gene ClinGen rates haploinsufficient
  (dosage score 3) or in a recessive carrier gene, Tier 3 or Tier 4 (`carrier_genes`; design
  decision, 2026-10-21: recessive genes never get score 3); else a non-finding, except on both copies in a recessive gene, which
  goes under research (`lof_gene_gate`, which slice 3b's deletions share). It must then
  pass every gate of the LoF filter (design decision, 2026-10-23; `failed_gate`), or it goes
  under research with the gate named;
- in a listed gene (ACMG SF or a carrier list) it is placed as a pathogenic ClinVar claim
  would be (`clinvar.placement`): under risks in a dominant gene, carrier status for one copy
  in a recessive carrier gene; elsewhere it goes under research.

The reasons' "Tier 3" names the carrier lists, Tier 4 included. One copy in a recessive
gene outside them is kept as a statement for `promote`'s
compound heterozygotes: a non-finding in an ACMG SF gene (a missense or splice change), or
for a LoF that passes every gate in a gene neither haploinsufficient nor Tier 3; research,
as above, for a missense or splice change in a gene on neither list. No copy (dosage 0) is
a non-finding; a statement at an unknown genotype is unassessed, research included
(stage-1 amendment 5).

Every grade that is a statement is validity limited, with `PREDICTED_REASON` and
`GENETICIST_NOTE` in its reasons, and never actionable: a prediction is not a
classification, and it never outranks a ClinVar-classified finding (`query`).

The gates' anchors: each reads one published resource (the lab's GRCh38 FASTA, GIAB's
low-mappability stratification, VEP's MANE flags, UCSC's NMD Escape track for MANE, gnomAD's
pext, SpliceVault), fetched whole and pinned by the bootstrap, which a variant's own
annotation cannot move because each is read at the variant's position, not from it.
"""

import re
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from genokeeper.annotate.vep import START_LOST, VepRecord, locus
from genokeeper.config import Grading, PredictedThresholds
from genokeeper.evidence.clinvar import (
    BA1_REASON,
    BS1_REASON,
    GENOTYPE_REASON,
    ValidityGrade,
    placement,
)
from genokeeper.genelists import GeneList, carrier_genes, gene_mode
from genokeeper.store.db import Kind, Topic, Validity

PREDICTED_REASON = (
    "predicted, not classified; would need confirmation by a clinical lab"
)
GENETICIST_NOTE = "discuss with a geneticist"
LOF = "loss of function"
START_LOSS = "start-loss"
MISSENSE = "missense"
SPLICE = "splice"
# Why a record is no statement: the reason a non-finding carries.
BELOW_THRESHOLDS = (
    "no loss of function, and missense and splice scores below thresholds"
)
NOT_HAPLOINSUFFICIENT = (
    "loss of function in a gene neither haploinsufficient nor Tier 3"
)
ONE_RECESSIVE_COPY = "one copy in a recessive gene outside carrier Tier 3"
BOTH_COPIES = "both copies, in a recessive gene neither haploinsufficient nor Tier 3"
NOT_CARRIED = "no copy carried"
# The LoF filter's gates, in order; a failing one is named in the reason after FAILED_GATE.
FAILED_GATE = "failed LoF gate"
GATE_HOMOPOLYMER = "1 technical: an indel in a homopolymer"
GATE_MAPPABILITY = "1 technical: in GIAB's low-mappability regions"
GATE_TRANSCRIPT = "2 transcript: not on MANE Select or MANE Plus Clinical"
GATE_NO_STOP = "3 NMD: no new stop found on the transcript"
GATE_NMD = "3 NMD: the new stop escapes nonsense-mediated decay"
GATE_START_LOSS = "4 start-loss: never a finding"
GATE_PEXT = "5 expression: pext below half the gene's largest"
GATE_NO_PEXT = "5 expression: no pext at the base"
GATE_SPLICEVAULT = "6 splice: SpliceVault predicts an in-frame event"
GATE_NO_SPLICEVAULT = "6 splice: no SpliceVault prediction at the site"
LOW_MAPPABILITY = "low_mappability"  # the [strata] flag of GIAB's low-mappability BED
STOP_CONSEQUENCES = frozenset({"stop_gained", "frameshift_variant"})
SPLICE_CONSEQUENCES = frozenset({"splice_donor_variant", "splice_acceptor_variant"})
# Where a new stop escapes NMD (gate 3): within 50 nt of the last exon junction or in the
# last exon by UCSC's track (its rule 1, the CDS from 50 nt before the last junction on, and
# rule 2, a single coding exon); in the last exon by VEP's own EXON; within the first 100
# coding bases by the stop's CDS position. The track's rules 3 and 4 (a stop in an exon over
# 400 nt, not among the spec's) are not read. A MANE transcript the track does not list has
# no NMD data, which fails the gate: the gate never passes on missing data.
TRACK = "in the last exon or within 50 nt of its junction (UCSC's NMD Escape track)"
LAST_EXON = "in the last exon (VEP)"
START_PROXIMAL = "within the first 100 coding bases"
NO_NMD_DATA = "no NMD Escape data for the transcript"
GATE_NO_NMD = f"3 NMD: {NO_NMD_DATA}"
START_BASES = 100
_TRACK_RULES = frozenset({"1", "2"})
# ClinVar's molecular consequences that truncate: nonsense and frameshift.
TRUNCATING = frozenset({"SO:0001587", "SO:0001589"})
# Rescue (gate 3): an escaping stop passes when ClinVar holds this many P/LP truncating
# variants 3' of it (ClinGen SVI's PVS1 tree, "critical region").
RESCUE_TRUNCATIONS = 2
_TER = re.compile(r"p\.[A-Za-z]{3}(\d+)(?:[A-Za-z]{3})?(?:fsTer(\d+)|Ter)$")


@dataclass(frozen=True)
class EscapeRegion:
    """One rule 1 or 2 region of UCSC's NMD Escape track for one MANE transcript: BED
    coordinates (0-based start, end) and strand."""

    chrom: str
    start: int
    end: int
    strand: str


def load_escapes(bed: Path) -> dict[str, list[EscapeRegion]]:
    """The track's rule 1 and 2 regions per transcript (unversioned ENST), from the BED
    bigBedToBed wrote: chrom, start, end, gene, score, strand, thickStart, thickEnd,
    color, mouseover (`<b>Rule N</b> - ...`), transcript (one, versioned), RefSeq."""
    regions: dict[str, list[EscapeRegion]] = {}
    for line in bed.read_text().splitlines():
        fields = line.split("\t")
        if fields[9].removeprefix("<b>Rule ")[0] in _TRACK_RULES:
            region = EscapeRegion(fields[0], int(fields[1]), int(fields[2]), fields[5])
            regions.setdefault(fields[10].split(".")[0], []).append(region)
    return regions


def new_stop(record: VepRecord) -> int | None:
    """The CDS position of the first base of the stop codon a stop-gain or frameshift
    creates: a stop-gain's codon from HGVSp (`p.Gln412Ter`), a frameshift's from its
    `fsTerN` (the first changed codon plus N - 1); None without one (`fsTer?`)."""
    found = _TER.search(record.hgvsp or "")
    if found is None:
        return None
    codon = int(found[1]) + (int(found[2]) - 1 if found[2] else 0)
    return 3 * codon - 2


def _genomic(regions: Sequence[EscapeRegion], offset: int) -> tuple[str, int] | None:
    """The 1-based position `offset` bases into `regions` read 5' to 3'."""
    if not regions:
        return None
    forward = regions[0].strand != "-"
    for region in sorted(regions, key=lambda r: r.start, reverse=not forward):
        size = region.end - region.start
        if offset < size:
            pos = region.start + 1 + offset if forward else region.end - offset
            return region.chrom, pos
        offset -= size
    return None


@dataclass(frozen=True)
class Escape:
    """Where a new stop escapes nonsense-mediated decay (`TRACK`, `LAST_EXON`,
    `START_PROXIMAL` or `NO_NMD_DATA`) and the stop's genomic position (None where neither
    the track nor VEP places it)."""

    rule: str
    stop: tuple[str, int] | None


def escape(record: VepRecord, regions: Sequence[EscapeRegion]) -> Escape | None:
    """How the record's new stop (`new_stop`) escapes NMD on its MANE transcript, or None
    if it does not (see `TRACK`). The track's regions are the CDS's last bases, so the
    stop's CDS position (VEP's, of a CDS that `cds_length` gives with its stop codon, as
    the track counts it) says whether it lies among them, and the regions place it in the
    genome; a stop-gain's own position falling in a region is an escape too (the
    cross-check of the two)."""
    if not regions:
        return Escape(NO_NMD_DATA, None)
    stop, length = new_stop(record), record.cds_length
    if stop is not None and stop <= START_BASES:
        # First, whatever else matches: a stop near the start is never rescued.
        return Escape(START_PROXIMAL, None)
    tail = sum(r.end - r.start for r in regions)
    if stop is not None and length is not None and stop > length - tail:
        return Escape(TRACK, _genomic(regions, stop - (length - tail) - 1))
    chrom, pos, _, _ = locus(record.vid)
    at = (chrom, pos) if "stop_gained" in record.terms else None
    if at and any(r.chrom == chrom and r.start < pos <= r.end for r in regions):
        return Escape(TRACK, at)
    exon, _, exons = (record.exon or "").partition("/")
    if exon and exon == exons:
        return Escape(LAST_EXON, at)
    return None


def truncations_3prime(
    found: Escape,
    regions: Sequence[EscapeRegion],
    gene: str,
    truncations: Iterable[tuple[str, str, int]],
) -> int:
    """How many of ClinVar's P/LP truncating variants (gene, chrom, pos) in `gene` lie 3'
    of the escaping stop in the coding bases of its transcript: in one of the track's
    `regions` for it (its CDS's last bases, the stop's own region among them), so one in
    the last intron, on an exon only another transcript has, is not counted (BR-014)."""
    if found.stop is None or not regions:
        return 0
    chrom, stop = found.stop
    forward = regions[0].strand != "-"
    return sum(
        1
        for g, c, pos in truncations
        if g == gene
        and c == chrom
        and (pos > stop if forward else pos < stop)
        and any(r.chrom == c and r.start < pos <= r.end for r in regions)
    )


def homopolymer_runs(fasta: Path, vids: Iterable[str]) -> dict[str, int]:
    """For each indel vid (normalised, so left-aligned with one anchor base), the length
    of the reference run of its one inserted or deleted base starting right after the
    anchor; 0 for an SNV or an indel of mixed bases."""
    runs: dict[str, int] = {}
    with pysam.FastaFile(str(fasta)) as ref:
        for vid in vids:
            chrom, pos, a, b = locus(vid)
            changed = (a if len(a) > len(b) else b)[1:]
            if len(a) == len(b) or len(set(changed)) != 1:
                runs[vid] = 0
                continue
            seq = ref.fetch(chrom, pos, pos + 100).upper()
            runs[vid] = len(seq) - len(seq.lstrip(changed[0]))
    return runs


@dataclass(frozen=True)
class Site:
    """A record with what its grade reads beside it: the user's dosage, gnomAD's grpmax
    filtering allele frequency, the homopolymer run (`homopolymer_runs`), the call's
    context flags, where its new stop escapes NMD (`escape`) and how many ClinVar
    truncations lie 3' of that stop (`truncations_3prime`)."""

    record: VepRecord
    dosage: str
    faf: float | None
    homopolymer: int
    flags: Collection[str]
    escape: Escape | None
    truncations_3prime: int


def predicted_rule(record: VepRecord, t: PredictedThresholds) -> tuple[str, str] | None:
    """The record's rule and the source level stating it, or None below every one."""
    if record.lof:
        return LOF, f"{record.consequence} on {record.transcript}"
    if START_LOST in record.terms:
        return START_LOSS, f"{record.consequence} on {record.transcript}"
    if "missense_variant" in record.terms:
        if record.revel is not None:
            if record.revel >= t.revel:
                return MISSENSE, f"missense, REVEL {record.revel:g}"
        elif (
            record.alphamissense is not None and record.alphamissense >= t.alphamissense
        ):
            return MISSENSE, f"missense, AlphaMissense {record.alphamissense:g}"
    if record.spliceai is not None and record.spliceai >= t.spliceai:
        return SPLICE, f"splice, SpliceAI {record.spliceai:g}"
    return None


def failed_gate(site: Site, t: PredictedThresholds) -> str | None:
    """The first gate of the LoF filter `site` fails, or None when it passes every one."""
    r = site.record
    if site.homopolymer >= t.homopolymer:
        return GATE_HOMOPOLYMER
    if LOW_MAPPABILITY in site.flags:
        return GATE_MAPPABILITY
    if not r.mane:
        return GATE_TRANSCRIPT
    if r.terms & STOP_CONSEQUENCES:
        if new_stop(r) is None:
            return GATE_NO_STOP
        found = site.escape
        if found is not None and found.rule == NO_NMD_DATA:
            return GATE_NO_NMD
        # The rescue is for a stop near the 3' end; one near the start may re-initiate.
        rescued = site.truncations_3prime >= RESCUE_TRUNCATIONS
        if found is not None and (found.rule == START_PROXIMAL or not rescued):
            return f"{GATE_NMD} ({found.rule})"
    if START_LOST in r.terms:
        return GATE_START_LOSS
    if r.pext is None:
        return GATE_NO_PEXT
    if r.pext < t.min_pext:
        return GATE_PEXT
    if r.terms & SPLICE_CONSEQUENCES:
        if r.splicevault_in_frame is None:
            return GATE_NO_SPLICEVAULT
        if r.splicevault_in_frame:
            return GATE_SPLICEVAULT
    return None


def _non_finding(reason: str) -> ValidityGrade:
    return ValidityGrade(Kind.NON_FINDING, None, "", reason, None, None)


def _statement(
    kind: Kind, topic: Topic | None, source_level: str, reasons: Sequence[str]
) -> ValidityGrade:
    reason = "; ".join([PREDICTED_REASON, *reasons, GENETICIST_NOTE])
    return ValidityGrade(kind, Validity.LIMITED, source_level, reason, topic, None)


def _research(dosage: str, source_level: str, reasons: Sequence[str]) -> ValidityGrade:
    """Under research; unassessed where the genotype is unknown (stage-1 amendment 5)."""
    if dosage == "unknown":
        return _statement(
            Kind.UNASSESSED, Topic.RESEARCH, source_level, [*reasons, GENOTYPE_REASON]
        )
    return _statement(Kind.FINDING, Topic.RESEARCH, source_level, reasons)


def _one_copy(source_level: str, rule: str) -> ValidityGrade:
    """One copy in a recessive gene outside Tier 3: a non-finding stored as a statement,
    so a second het in the gene can make it a compound heterozygote (`promote`)."""
    return _statement(Kind.NON_FINDING, None, source_level, [rule, ONE_RECESSIVE_COPY])


def lof_gene_gate(
    gene: str,
    *,
    mode: str | None,
    dosage: str,
    filtered: bool,
    carrier: Collection[str],
    haploinsufficient: Collection[str],
) -> str | None:
    """Whether a loss of function in `gene` (inheritance `mode`, the user's `dosage`) can
    count: None in a `carrier` gene (`carrier_genes`) whose `mode` is not AD, or one ClinGen
    rates haploinsufficient. A carrier gene whose `mode` is AD (`genelists.gene_mode`:
    SYNE1 and nine other Tier 4 genes are AD in CGD) keeps the gate, as one copy there is
    no carrier status. Else why
    not: `BOTH_COPIES` (research) on both copies, or an unknown genotype, in a recessive
    gene, as recessive genes never get dosage score 3; `ONE_RECESSIVE_COPY` (a non-finding
    kept as a statement) on one copy there when it passed the LoF filter (`filtered`);
    `NOT_HAPLOINSUFFICIENT` (a non-finding) anywhere else. 1e's predicted LoF and 3b's
    deletions both go through here."""
    if (gene in carrier and mode != "AD") or gene in haploinsufficient:
        return None
    if mode == "AR" and dosage in ("2", "unknown"):
        return BOTH_COPIES
    if mode == "AR" and dosage == "1" and filtered:
        return ONE_RECESSIVE_COPY
    return NOT_HAPLOINSUFFICIENT


def grade_predicted(
    site: Site,
    *,
    haploinsufficient: frozenset[str],
    lists: Mapping[str, GeneList],
    grading: Grading,
    inheritance: str | None = None,
    male: bool = False,
    haploid: bool = False,
) -> ValidityGrade:
    """Grade one record (see the module docstring); a non-finding's reason says why it is
    no statement. `inheritance` is the gene's outside the lists (the CGD's), `haploid`
    whether the site has one copy."""
    r = site.record
    t = grading.predicted
    found = predicted_rule(r, t)
    if found is None:
        return _non_finding(BELOW_THRESHOLDS)
    rule, source_level = found
    if site.dosage == "0":
        return _non_finding(NOT_CARRIED)
    sf, carrier = lists["acmg_sf"].genes, carrier_genes(lists)
    mode = gene_mode(r.gene, lists, inheritance)
    bs1 = grading.bs1_for(r.gene, autosomal_recessive=mode == "AR")
    if site.faf is not None and site.faf > grading.ba1:
        return _non_finding(BA1_REASON)
    if site.faf is not None and bs1 is not None and site.faf > bs1:
        return _non_finding(BS1_REASON)
    if rule in (LOF, START_LOSS):
        gate = failed_gate(site, t)
        failed = [] if gate is None else [f"{FAILED_GATE} {gate}"]
        blocked = lof_gene_gate(
            r.gene,
            mode=mode,
            dosage=site.dosage,
            filtered=gate is None,
            carrier=carrier,
            haploinsufficient=haploinsufficient,
        )
        if blocked == BOTH_COPIES:
            # Unassessed where the genotype is unknown: it could be both copies.
            reasons = [rule, BOTH_COPIES, *failed]
            return _research(site.dosage, source_level, reasons)
        if blocked == ONE_RECESSIVE_COPY:
            return _one_copy(source_level, rule)
        if blocked is not None:
            return _non_finding(blocked)
        if gate is not None:
            return _research(site.dosage, source_level, failed)
    if r.gene not in sf and r.gene not in carrier:
        return _research(site.dosage, source_level, [rule])
    kind, topic, reasons = placement(
        cls="pathogenic",
        gene=r.gene,
        lists=dict(lists),
        dosage=site.dosage,
        inheritance=inheritance,
        male=male,
        haploid=haploid,
    )
    if kind is Kind.NON_FINDING:
        return _one_copy(source_level, rule)
    return _statement(kind, topic, source_level, [rule, *reasons])
