"""Statements, calls and variants for the compound-heterozygote promotion, one case each.

`promotion_rows(run_id)` gives `(tool, variant, call, statement)` per statement, in the order
a run writes them: the first source's, then the second's. A later row at a vid an earlier one
holds repeats a variant or call that the store drops (the first at a vid wins). Nothing here
is a genome: every vid, gene and position is made up, apart from the gene-list genes. A
statement's locus is its vid's: the promotion reads no other.
"""

from genokeeper.sources.panel import clinvar_claim

type Fields = dict[str, object]

# Gene -> CGD inheritance. BRCA1 is AD in ACMG SF, which wins over the AR given here.
INHERITANCE = {
    "GENEAR": "AR",
    "GENEVUS": "AR",
    "GENESAME": "AR",
    "GENEFW": "AR",
    "GENEFW2": "AR",
    "GENENOCHR": "AR",
    "BRCA1": "AR",
}


# A single het's grading in a carrier gene (Tier 3 or 4).
CARRIER = {"kind": "finding", "topic": "carrier"}


def haploid(chrom: str, _pos: int) -> bool:
    """A male's X outside the PARs."""
    return chrom == "chrX"


def _vid(chrom: str, pos: int) -> str:
    return f"GRCh38:{chrom}:{pos}:A:G"


def _variant(run_id: str, chrom: str, pos: int, *, at: str | None = None) -> Fields:
    return {"run_id": run_id, "vid": at or _vid(chrom, pos), "rsid": None}


def _call(run_id: str, vid: str, dosage: str) -> Fields:
    return {
        "run_id": run_id,
        "vid": vid,
        "gt": {"0": "0/0", "1": "0/1", "2": "1/1"}.get(dosage, "./."),
        "dosage": dosage,
        "p_wrong": 1e-5,
        "level": "high",
        "gq": 50,
        "dp": 30,
        "ad": [15, 15],
        "vaf": 0.5,
        "row_type": "PASS",
        "flags": [],
        "warnings": [],
    }


def _statement(
    run_id: str,
    n: int,
    vid: str | None,
    gene: str,
    clnsig: str = "Pathogenic",
    *,
    kind: str = "non_finding",
    reason: str | None = None,
    topic: str | None = None,
) -> Fields:
    return {
        "run_id": run_id,
        "statement_id": f"clinvar:{n}",
        "kind": kind,
        "vid": vid,
        "gene": gene,
        "claim": clinvar_claim(
            variation_id=n,
            clnsig=clnsig,
            clnrevstat="criteria_provided,_single_submitter",
            condition="not provided",
        ),
        "validity": "strong",
        "source_level": "criteria_provided,_single_submitter",
        "reason": reason,
        "effect": None,
        "actionable": False,
        "topic": topic or ("risks" if kind != "non_finding" else None),
        "source": "clinvar",
        "source_id": str(n),
        "source_version": "20261012",
        "points": None,
    }


def _site(
    run_id: str,
    tool: str,
    n: int,
    gene: str,
    chrom: str,
    pos: int,
    dosage: str,
    clnsig: str = "Pathogenic",
    *,
    kind: str = "non_finding",
    reason: str | None = None,
    at: str | None = None,
    topic: str | None = None,
) -> tuple[str, Fields | None, Fields | None, Fields]:
    """One statement with its variant and call, at `at` (default: its own vid)."""
    vid = at or _vid(chrom, pos)
    return (
        tool,
        _variant(run_id, chrom, pos, at=vid),
        _call(run_id, vid, dosage),
        _statement(run_id, n, vid, gene, clnsig, kind=kind, reason=reason, topic=topic),
    )


def promotion_rows(
    run_id: str,
) -> list[tuple[str, Fields | None, Fields | None, Fields]]:
    """Every case, keyed by statement number in the comment beside it."""
    panel, oc = "first", "second"
    no_call_vid = _vid("chr1", 500)
    return [
        # 1, 2: two MUTYH (AR, ACMG SF, carrier Tier 4) hets, each a carrier finding: both
        # promoted, a reason kept before ours.
        _site(run_id, panel, 1, "MUTYH", "chr1", 100, "1", **CARRIER),
        _site(
            run_id,
            panel,
            2,
            "MUTYH",
            "chr1",
            200,
            "1",
            "Likely_pathogenic",
            reason="a reason of its own",
            **CARRIER,
        ),
        # 3: unassessed, 4: two copies gone (dosage 0): neither counts.
        _site(run_id, panel, 3, "MUTYH", "chr1", 300, "unknown", kind="unassessed"),
        _site(run_id, panel, 4, "MUTYH", "chr1", 400, "0"),
        # 5: one HFE het alone.
        _site(run_id, panel, 5, "HFE", "chr6", 100, "1"),
        # 6, 7: CFTR is in carrier Tier 3.
        _site(run_id, panel, 6, "CFTR", "chr7", 100, "1"),
        _site(run_id, panel, 7, "CFTR", "chr7", 200, "1"),
        # 8, 9: BRCA1 is AD in ACMG SF, whatever CGD says.
        _site(run_id, panel, 8, "BRCA1", "chr17", 100, "1"),
        _site(run_id, panel, 9, "BRCA1", "chr17", 200, "1"),
        # 10, 11: CASQ2 (AR) on a haploid X.
        _site(run_id, panel, 10, "CASQ2", "chrX", 100, "1"),
        _site(run_id, panel, 11, "CASQ2", "chrX", 200, "1"),
        # 12, 15 VUS (14 its pathogenic partner), 16 and 17 at one vid.
        _site(run_id, panel, 12, "GENEAR", "chr2", 100, "1"),
        _site(run_id, panel, 14, "GENEVUS", "chr3", 100, "1"),
        _site(run_id, panel, 15, "GENEVUS", "chr3", 200, "1", "Uncertain_significance"),
        _site(run_id, panel, 16, "GENESAME", "chr4", 100, "1"),
        _site(run_id, panel, 17, "GENESAME", "chr4", 101, "1", at=_vid("chr4", 100)),
        # 18: a gene-level statement (no vid); 19: a statement whose vid has no call.
        (panel, None, None, _statement(run_id, 18, None, "MUTYH")),
        (
            panel,
            _variant(run_id, "chr1", 500),
            None,
            _statement(run_id, 19, no_call_vid, "MUTYH"),
        ),
        # 20, 21 hets; 22 (the second source) at 21's vid with its own zero-copy call,
        # which the store drops: 22 is graded on 21's het call and is promoted with them.
        _site(run_id, panel, 20, "GENEFW", "chr5", 100, "1"),
        _site(run_id, panel, 21, "GENEFW", "chr5", 200, "1"),
        # 23 at a zero-copy call; 24 a het; 25 (the second source) at 23's vid with a
        # het call the store drops: only 24 counts, alone.
        _site(run_id, panel, 23, "GENEFW2", "chr8", 100, "0"),
        _site(run_id, panel, 24, "GENEFW2", "chr8", 200, "1"),
        # 26 has no variant row (its vid's locus counts); 27 its partner.
        (
            panel,
            None,
            _call(run_id, _vid("chr9", 100), "1"),
            _statement(run_id, 26, _vid("chr9", 100), "GENENOCHR"),
        ),
        _site(run_id, panel, 27, "GENENOCHR", "chr9", 200, "1"),
        # The second source's: 13 pairs with the first's 12 across sources.
        _site(run_id, oc, 13, "GENEAR", "chr2", 200, "1"),
        _site(run_id, oc, 22, "GENEFW", "chrX", 300, "0", at=_vid("chr5", 200)),
        _site(run_id, oc, 25, "GENEFW2", "chr8", 300, "1", at=_vid("chr8", 100)),
    ]
