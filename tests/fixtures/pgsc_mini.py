"""pgsc_calc's own files as a run leaves them, written into a test's tmp_path: a task's
launcher (`.command.run`, as Nextflow 25.10 writes it for the docker executor) and the
results of a sampleset (its match summary, scores and population similarity)."""

import gzip
from pathlib import Path

from genokeeper.containers import IMAGES
from genokeeper.prs import pgsc

PINNED = {label: IMAGES[name] for label, name in pgsc.PGSC_LABELS.items()}
# A task's launcher as Nextflow 25.10.0 wrote it for pgsc_calc v2.3.0's docker executor in
# the smoke run of 2026-10-21 (its header, and nxf_launch's docker line verbatim but for
# the paths).
COMMAND_RUN = """#!/bin/bash
### ---
### name: 'PGSCATALOG_PGSCCALC:PGSCCALC:APPLY_SCORE:PLINK2_SCORE (genomelab chromosome ALL)'
### container: '{image}'
### ...
set -e
nxf_launch() {{
    docker run -i --cpu-shares 2048 --memory 16384m -e "NXF_TASK_WORKDIR" -e "NXF_DEBUG=${{NXF_DEBUG:=0}}" -v {wd}:{wd} -w "$NXF_TASK_WORKDIR" {options} --name $NXF_BOXID {image} /bin/bash -c "eval $(nxf_container_env); /bin/bash -euo pipefail {wd}/.command.run nxf_trace"
}}
on_term() {{
    docker stop $NXF_BOXID
}}
"""
CONTAINED = "-u $(id -u):$(id -g) --platform=linux/arm64 --network none"


def write_task(
    work: Path, image: str = PINNED["plink2"], options: str = CONTAINED, n: int = 0
) -> Path:
    task = work / f"{n:02x}" / f"{n:030x}"
    task.mkdir(parents=True)
    (task / ".command.run").write_text(
        COMMAND_RUN.format(wd=task, options=options, image=image)
    )
    return task


def write_results(
    results: Path, rates: dict[str, float], population: str = "EUR", low: str = "False"
) -> None:
    """pgsc_calc's outputs (`--outdir results`) for sampleset `genomelab` run with
    `--min_overlap 0`: the match summary (one matched and one unmatched row per score),
    a percentile for every score, and the population similarity with the low-confidence
    flag `low`."""
    match = results / pgsc.SAMPLESET / "match"
    score = results / pgsc.SAMPLESET / "score"
    match.mkdir(parents=True)
    score.mkdir(parents=True)
    head = (
        "dataset,accession,score_pass,match_status,ambiguous,is_multiallelic,"
        "duplicate_best_match,duplicate_ID,match_flipped,match_IDs,count,percent\n"
    )
    rows = []
    for pgs_id, rate in rates.items():
        matched, total = round(rate * 1000), 1000
        for status, count in (("matched", matched), ("unmatched", total - matched)):
            rows.append(
                f"genomelab,{pgs_id},True,{status},False,False,False,False,False,"
                f"True,{count},{100 * count / total}\n"
            )
    (match / "genomelab_summary.csv").write_text(head + "".join(rows))
    pgs = ["sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tpercentile_MostSimilarPop"]
    for pgs_id in rates:
        pgs.append(f"reference\tHG1\tHG1\t{pgs_id}\t0.1\t0.0\t50.0")
        pgs.append(f"genomelab\tSAMPLE\tSAMPLE\t{pgs_id}\t0.3\t1.28\t90.0")
    (score / "genomelab_pgs.txt.gz").write_bytes(gzip.compress("\n".join(pgs).encode()))
    similarity = (
        "sampleset\tIID\tPC1\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\n"
        "reference\tHG1\t0.1\tEUR\tFalse\tTrue\n"
        f"genomelab\tSAMPLE\t0.2\t{population}\t{low}\tFalse\n"
    )
    (score / "genomelab_popsimilarity.txt.gz").write_bytes(
        gzip.compress(similarity.encode())
    )
