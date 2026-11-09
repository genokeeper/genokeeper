"""StellarPGx, CYP2D6's third caller (roadmap, "Next"): where BCyrius and Aldy do not agree on
the CPIC activity score, its score settles the gene when it equals one of theirs (`vote`).

StellarPGx is a Nextflow pipeline, run offline from the commit `STELLARPGX_COMMIT` (bootstrap)
in its upstream image, pinned by digest, under emulation (no linux/arm64 build). Its tasks get
`--network none`, and `check_containment` reads every task's launcher afterwards: the same
threat and anchor as pgsc_calc's. The image's graphtyper (htslib 1.10) cannot decode our CRAMs
("Failure to decode slice" on all four GIAB libraries), so it reads a BAM of the regions its
tasks use, cut with our samtools (`_REGIONS`). A third caller that fails changes nothing: the gene stays
what the two callers made of it.

Only StellarPGx's activity score is voted, as for the other two. Its score comes from the
function table it bundles, which is not compared with CPIC's here.
"""

import shutil
from dataclasses import replace
from pathlib import Path

from genokeeper.containers import IMAGES
from genokeeper.pgx.cyp2d6 import canonical_diplotype, shown_call
from genokeeper.pgx.outside import (
    LABEL,
    CallerCall,
    CallerFailed,
    OutsideCall,
    cpic_score,
)
from genokeeper.prs.pgsc import (
    ContainmentRefused,
    check_containment,
    java_home,
    nextflow_path,
    nxf_home,
)
from genokeeper.tools import gated, run_tool

STELLARPGX = "StellarPGx"
STELLARPGX_VERSION = "1.2.8"
STELLARPGX_COMMIT = "1952a4494fd6e9f809be87baa0a56c379e1d757b"
ALLELES = "stellarpgx.alleles"
# GRCh38: CYP2D6, its pseudogene and spacer (main.nf's region_a1, res_hg38/test3.bed), and
# test3.bed's two depth controls.
_REGIONS = (
    "chr22:42120000-42150000",
    "chr12:47841537-47905022",
    "chr7:55019017-55211628",
)
_HEADERS = (
    "Result:",
    "Activity score:",
    "Metaboliser status:",
    "Likely background alleles:",
)


def stellarpgx_dir(lab: Path) -> Path:
    return lab / "stellarpgx" / STELLARPGX_COMMIT


def read_stellarpgx(path: Path) -> CallerCall:
    """StellarPGx's call from its `.alleles` report: the one diplotype under `Result:` (a
    possible novel allele, `No_call` or an empty result is a no-call) and the score under
    `Activity score:` (`Indeterminate`: none)."""
    blocks: dict[str, list[str]] = {}
    current: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() in _HEADERS:
            current = blocks.setdefault(line.strip(), [])
        elif line.strip():
            current.append(line.strip())
    if "Result:" not in blocks or "Activity score:" not in blocks:
        raise CallerFailed(f"{path}: not a StellarPGx report")
    result = blocks["Result:"]
    diplotype = result[0] if len(result) == 1 and result[0].startswith("*") else None
    score = " ".join(blocks["Activity score:"])
    activity = None if diplotype is None else cpic_score(score)
    return CallerCall(
        caller=STELLARPGX,
        diplotype=diplotype,
        phenotype=None,
        activity_score=activity,
        abstained=diplotype is not None and activity is None,
    )


def vote(
    pair: OutsideCall, first: CallerCall, second: CallerCall, third: CallerCall | None
) -> OutsideCall:
    """The two callers' `pair` with the `third` caller's call (None: it failed) shown, and
    agreed when they did not agree and its score equals the one of theirs that stands: a caller
    with no score or that filtered its own call has no vote."""
    if third is None:
        return replace(pair, callers={**pair.callers, STELLARPGX: "failed"})
    callers = {**pair.callers, third.caller: shown_call(third)}
    score = third.activity_score
    if pair.agreed:
        dissent = None
        if score is not None and score != pair.activity_score:
            dissent = f"{third.caller} states {score} against the agreed {pair.activity_score}"
        return replace(pair, callers=callers, note=dissent)
    backing = [
        c
        for c in (first, second)
        if score is not None and c.activity_score == score and c.filtered is None
    ]
    if len(backing) != 1:
        if score is not None and pair.activity_score is not None:
            reason = (
                f"{third.caller} states {score} against the one score"
                f" {pair.activity_score} the other two left"
            )
            return replace(
                pair,
                activity_score=None,
                callers=callers,
                unresolved_reason=reason,
                note=None,
            )
        return replace(pair, callers=callers)
    (one,) = backing
    other = second if one is first else first
    names = (
        canonical_diplotype(one.diplotype or ""),
        canonical_diplotype(third.diplotype or ""),
    )
    note = (
        f"two of three callers: {one.caller} and {third.caller} agree on the activity"
        f" score {score}; {other.caller} called {shown_call(other)}"
    )
    return OutsideCall(
        pair.gene,
        names[0] if names[0] is not None and names[0] == names[1] else None,
        None,
        score,
        callers,
        True,
        None,
        note,
    )


def run_stellarpgx(cram: Path, fasta: Path, out_dir: Path, lab: Path) -> CallerCall:
    """StellarPGx on CYP2D6 in `cram` (indexed, written against `fasta`, hg38), its report kept as `ALLELES` in `out_dir`. Raises `CallerFailed` when it is not
    installed or writes no report (its last task ignores its own errors), and what Nextflow's
    run or the containment check raises."""
    pipeline = stellarpgx_dir(lab)
    if not (pipeline / "main.nf").is_file():
        raise CallerFailed(
            f"StellarPGx is not installed at {pipeline}: run the bootstrap"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    config = out_dir / "genokeeper.config"
    config.write_text(
        "singularity.enabled = false\n"
        "docker.enabled = true\n"
        "docker.runOptions = '-u $(id -u):$(id -g) --network none --platform linux/amd64'\n"
        f"process.container = '{IMAGES['stellarpgx']}'\n",
        encoding="utf-8",
    )
    slice_dir = out_dir / "slice"
    slice_dir.mkdir(exist_ok=True)
    bam = slice_dir / "sample.bam"
    run_tool(
        [
            "samtools",
            "view",
            "-b",
            "-T",
            str(fasta),
            "-o",
            str(bam),
            str(cram),
            *_REGIONS,
        ]
    )
    run_tool(["samtools", "index", str(bam)])
    results = out_dir / "results"
    command = [
        str(nextflow_path(lab)), "run", str(pipeline / "main.nf"),
        "-profile", "standard", "-c", str(config), "-w", str(out_dir / "work"),
        "-ansi-log", "false",
        "--gene", "cyp2d6", "--build", "hg38",
        "--ref_file", str(fasta),
        "--in_bam", f"{slice_dir}/sample*{{bam,bai}}",
        "--out_dir", str(results),
        "--db_init", str(pipeline / "database"),
        "--res_init", str(pipeline / "resources"),
        "--caller_init", str(pipeline / "scripts"),
    ]  # fmt: skip
    env = {
        "NXF_HOME": str(nxf_home(lab)),
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
        "NXF_JAVA_HOME": str(java_home(lab)),
        "NXF_TEMP": str(out_dir),
    }
    run_tool(gated(LABEL, command), env=env, cwd=out_dir)
    try:
        check_containment(out_dir / "work", pinned={IMAGES["stellarpgx"]})
    except ContainmentRefused as refused:  # its report is not used; the run goes on
        raise CallerFailed(str(refused)) from refused
    reports = list((results / "cyp2d6" / "alleles").glob("*.alleles"))
    if len(reports) != 1:
        raise CallerFailed(
            f"StellarPGx wrote {len(reports)} reports in {results}, not one"
        )
    (out_dir / ALLELES).write_bytes(reports[0].read_bytes())
    shutil.rmtree(slice_dir)
    return read_stellarpgx(out_dir / ALLELES)
