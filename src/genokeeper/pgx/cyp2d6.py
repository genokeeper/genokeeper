"""CYP2D6 from BCyrius and Aldy, agreed on the CPIC activity score (track design, slice 1b,
"Outside calls").

PharmCAT cannot call CYP2D6 from a VCF (copy number, hybrids), so two callers read the CRAM:
BCyrius 2.2.0 and Aldy 4.8.3 (`aldy genotype -p wgs --genome hg38 -r <fasta>`), each offline
in its image (`containers.IMAGES`), built here from pinned sources (`containers.BUILT`), behind
the core gate as `outside.LABEL`. Each returns a `CallerCall`.

Their star alleles need not match: they are compared on the CPIC activity score each caller
states for its call (BCyrius's `Activity score` column, Aldy's `cpic_score`), read by
`outside.cpic_score`, where a score of 3 or more is CPIC's `≥3.0` bin, and agreed only when
both scores are known and equal (`agree`). The agreed score goes to PharmCAT, with the
diplotype only where both callers name the same one. When one caller names a diplotype but
states no score (Aldy's `cpic=indeterminate` for a star allele carrying variants beyond its
definition) and the other states one, that score goes alone, not agreed, so the gene call is
at most medium (design decision, 2026-10-23); any other difference leaves the gene unresolved
with both calls shown, as does a call BCyrius filters (a Filter other than PASS). Each score
rests on the CPIC tables its caller bundles; tables that have drifted apart leave the gene
unresolved.
"""

import csv
import re
from pathlib import Path

from genokeeper.containers import run_container
from genokeeper.pgx.outside import (
    LABEL,
    CallerCall,
    CallerFailed,
    OutsideCall,
    cpic_score,
)

BCYRIUS_VERSION = "2.2.0"
ALDY_VERSION = "4.8.3"
BCYRIUS, ALDY = "BCyrius", "Aldy"
# Alternatives within one call: Aldy's solutions joined by ` or `, BCyrius's by `;`.
_ALTERNATIVES = re.compile(r"\s+or\s+|;")
# One allele: its star number, a sub-allele number, Aldy's markers, and a copy count.
_ALLELE = re.compile(
    r"^\*(?P<star>\d+)(?:\.\d+)?(?:[A-Z]+\d*)?(?:x(?P<copies>\d+|≥3|N))?$"
)
_NO_CALL = {"", "None", "-", "n/a"}


def alternatives(diplotype: str) -> list[str]:
    """The diplotypes one call offers (Aldy's solutions, BCyrius's `;` alternatives)."""
    return [d.strip() for d in _ALTERNATIVES.split(diplotype)]


def canonical_diplotype(diplotype: str) -> str | None:
    """One diplotype as PharmCAT names it: each copy's star allele, a run of one allele as
    `x<N>`, the haplotypes in star order; None for alternatives, unassigned copies, a bin
    or an unknown count."""
    haplotypes = diplotype.split("/")
    if len(alternatives(diplotype)) > 1 or len(haplotypes) != 2 or "_" in diplotype:
        return None
    named: list[tuple[list[int], str]] = []
    for haplotype in haplotypes:
        stars: list[int] = []
        for part in haplotype.split("+"):
            match = _ALLELE.match(part.strip())
            copies = None if match is None else match["copies"]
            if match is None or (copies is not None and not copies.isdigit()):
                return None
            stars += [int(match["star"])] * (int(copies) if copies else 1)
        runs: list[str] = []
        for i, star in enumerate(stars):
            if i == 0 or star != stars[i - 1]:
                runs.append(f"*{star}")
            else:
                count = int(runs[-1].partition("x")[2] or 1) + 1
                runs[-1] = f"*{star}x{count}"
        named.append((stars, "+".join(runs)))
    return "/".join(name for _, name in sorted(named))


def shown_call(call: CallerCall) -> str:
    if call.diplotype is None:
        return "no call"
    if call.activity_score is None:
        return f"{call.diplotype} (no activity score)"
    return f"{call.diplotype} (activity score {call.activity_score})"


def agree(a: CallerCall, b: CallerCall, gene: str) -> OutsideCall:
    """The two callers' calls of `gene` as one outside call: agreed only when both stated
    scores are known and equal; the one stated score, not agreed, when the other caller
    named a diplotype without a score; otherwise unresolved, with why and both calls."""
    calls = (a, b)
    callers = {c.caller: shown_call(c) for c in calls}
    silent = [c.caller for c in calls if c.diplotype is None]
    filtered = [
        f"{c.caller} filtered its {c.diplotype} as {c.filtered}"
        for c in calls
        if c.diplotype is not None and c.filtered is not None
    ]
    unscored = [
        f"no activity score for {c.caller}'s {c.diplotype}"
        for c in calls
        if c.diplotype is not None and c.activity_score is None
    ]
    scored = [c for c in calls if c.activity_score is not None]
    abstained = [c for c in calls if c.abstained]
    reason = None
    if silent:
        reason = f"{' and '.join(silent)} made no call"
    elif filtered:  # the caller doubts its own call: nothing to weigh against the other
        reason = "; ".join(filtered)
    elif len(scored) == 1 and len(abstained) == 1:
        (only,), (other,) = scored, abstained
        note = (
            f"one caller's score: {other.caller} stated no activity score for"
            f" {other.diplotype}, so {only.caller}'s {only.activity_score} was not checked"
        )
        return OutsideCall(
            gene, None, None, only.activity_score, callers, False, None, note
        )
    elif unscored:
        reason = "; ".join(unscored)
    elif a.activity_score != b.activity_score:
        shown = f"{a.activity_score} against {b.activity_score}"
        reason = f"{a.caller} and {b.caller} differ on the activity score ({shown})"
    if reason is not None:
        return OutsideCall(gene, None, None, None, callers, False, reason)
    names = [canonical_diplotype(c.diplotype or "") for c in calls]
    diplotype = names[0] if names[0] is not None and names[0] == names[1] else None
    return OutsideCall(gene, diplotype, None, a.activity_score, callers, True, None)


# --- BCyrius --------------------------------------------------------------------------------

BCYRIUS_ID = "bcyrius"


def read_bcyrius(tsv: Path) -> CallerCall:
    """BCyrius's call from its TSV (one sample): genotype `None` is its no-call; its
    activity score and phenotype as it states them (`-`, `n/a` or several: none), and its
    Filter other than `PASS` (`LowQ_high_CN`, `Not_assigned_to_haplotypes`)."""
    lines = tsv.read_text(encoding="utf-8").splitlines()
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if len(rows) != 1:
        raise CallerFailed(f"{tsv}: {len(rows)} BCyrius results, not one")
    (row,) = rows
    genotype = row["Genotype"].strip()
    score = row["Activity score"].strip()
    phenotype = row["Predicted phenotype"].strip()
    filtered = row["Filter"].strip()
    return CallerCall(
        caller=BCYRIUS,
        diplotype=None if genotype in _NO_CALL else genotype,
        phenotype=None if phenotype in _NO_CALL or ";" in phenotype else phenotype,
        activity_score=cpic_score(score),
        abstained=genotype not in _NO_CALL and score in _NO_CALL,
        filtered=None if filtered == "PASS" else filtered,
    )


def run_bcyrius(cram: Path, fasta: Path, out_dir: Path) -> CallerCall:
    """BCyrius on `cram` (indexed, written against `fasta`, chr-named contigs), offline;
    its TSV and JSON stay in `out_dir`. Raises `CallerFailed` without a TSV, and what
    `run_container` raises."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run_container(
        "bcyrius",
        [
            "bcyrius",
            "--input",
            f"/bcyrius-in/{cram.name}",
            "--reference",
            f"/bcyrius-ref/{fasta.name}",
            "--genome",
            "chr38",
            "--outDir",
            "/bcyrius-out",
            "--id",
            BCYRIUS_ID,
        ],
        mounts=[
            (cram.parent, "/bcyrius-in", True),
            (fasta.parent, "/bcyrius-ref", True),
            (out_dir, "/bcyrius-out", False),
        ],
        label=LABEL,
    )
    tsv = out_dir / f"{BCYRIUS_ID}.tsv"
    if not tsv.is_file():
        raise CallerFailed(f"BCyrius wrote no {tsv.name} in {out_dir}")
    return read_bcyrius(tsv)


# --- Aldy -----------------------------------------------------------------------------------


def read_aldy(path: Path) -> CallerCall:
    """Aldy's call from its output: each solution's major diplotype (several joined by
    ` or `; none is its no-call), its CPIC phenotype and score where every solution states
    the same."""
    solutions: list[str] = []
    majors: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#Solution "):
            solutions.append(line)
        elif line and not line.startswith("#"):
            fields = line.split("\t")
            majors.setdefault(fields[2], fields[3])
    stated = [
        dict(p.split("=", 1) for p in s.split("; ")[1:] if "=" in p) for s in solutions
    ]
    phenotypes = {s.get("cpic") for s in stated}
    scores = {s.get("cpic_score") for s in stated}
    unscored = scores == {None}
    score = scores.pop() if len(scores) == 1 else None
    phenotype = phenotypes.pop() if len(phenotypes) == 1 else None
    diplotype = " or ".join(dict.fromkeys(majors.values())) or None
    return CallerCall(
        caller=ALDY,
        diplotype=diplotype,
        phenotype=None if phenotype == "indeterminate" else phenotype,
        activity_score=None if score is None else cpic_score(score),
        abstained=diplotype is not None and unscored,
    )


def run_aldy(
    cram: Path, fasta: Path, out_dir: Path, gene: str, *, label: str = LABEL
) -> CallerCall:
    """Aldy on `gene` in `cram` (indexed, written against `fasta`), as WGS on hg38,
    offline behind `label` (slice 3d's CYP2A6 runs as `dna:paralogs`); its output stays in
    `out_dir` as `<gene>.aldy`. Raises `CallerFailed` without it, and what `run_container`
    raises. Aldy's licence (IURTC) permits non-commercial research only and excludes
    diagnostic use (README)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"{gene}.aldy"
    run_container(
        "aldy",
        [
            "aldy",
            "genotype",
            "-p",
            "wgs",
            "-g",
            gene,
            "--genome",
            "hg38",
            "-r",
            f"/aldy-ref/{fasta.name}",
            "-o",
            f"/aldy-out/{name}",
            f"/aldy-in/{cram.name}",
        ],
        mounts=[
            (cram.parent, "/aldy-in", True),
            (fasta.parent, "/aldy-ref", True),
            (out_dir, "/aldy-out", False),
        ],
        label=label,
    )
    if not (out_dir / name).is_file():
        raise CallerFailed(f"Aldy wrote no {name} in {out_dir}")
    return read_aldy(out_dir / name)
