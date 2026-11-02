"""The mitochondrial haplogroup (`genokeeper.haplogroup`, `sources.mito`): the chrM VCF
Haplogrep3 classifies, its run offline in its image with the tree baked in, its output,
and the statement a run stores. Haplogrep3's container is replaced by a stub that writes
what Haplogrep3 3.3.2 wrote, in the image docker/haplogrep3.Dockerfile builds, on a
synthetic chrM VCF of J1c's mutations (run on 2026-10-21)."""

import json
import re
from pathlib import Path
from typing import Any

import pysam
import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.genome_mini import HAPLOGREP_OUTPUT as OUTPUT
from fixtures.genome_mini import FakeHaplogrep

from genokeeper import haplogroup
from genokeeper.config import load_grading
from genokeeper.containers import DOCKER_DIR
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import load_gene_lists
from genokeeper.haplogroup import (
    RCRS_LENGTH,
    TREE,
    Haplogroup,
    read_haplogroup,
    run_haplogrep,
    uncovered,
    write_chrm_vcf,
)
from genokeeper.mito import MitoCall
from genokeeper.query import Filters, query
from genokeeper.sources import RunContext, write_collected
from genokeeper.sources.mito import MitoSource, haplogroup_validity
from genokeeper.store.db import Validity, connect, start_run, transaction
from genokeeper.store.writer import StoreWriter

MIN_DP = 10


def _covered(depth: int = 2000, **holes: int) -> dict[int, int]:
    found = dict.fromkeys(range(1, RCRS_LENGTH + 1), depth)
    found.update({int(pos.removeprefix("p")): dp for pos, dp in holes.items()})
    return found


def test_heteroplasmic_left_out(tmp_path: Path) -> None:
    """Only homoplasmic SNVs at positions min_dp reads cover reach Haplogrep3: a
    heteroplasmy would read as present, and a call where few reads cover it is no
    evidence."""
    calls = [
        MitoCall(73, "A", "G", 1.0, 2000),
        MitoCall(263, "A", "G", 0.95, 2000),
        MitoCall(5000, "A", "G", 0.2, 2000),  # a heteroplasmy
        MitoCall(7028, "C", "T", 0.99, 5),  # too few reads
    ]
    vcf = write_chrm_vcf(
        calls,
        _covered(p7028=5),
        tmp_path / "chrM.vcf",
        min_homoplasmy=0.9,
        min_dp=MIN_DP,
    )
    with pysam.VariantFile(str(vcf)) as f:
        rows = [(r.chrom, r.pos, r.ref, r.alts, r.samples[0]["GT"]) for r in f]
    assert rows == [
        ("chrM", 73, "A", ("G",), (1,)),
        ("chrM", 263, "A", ("G",), (1,)),
    ]


def test_uncovered_positions_listed() -> None:
    assert uncovered(_covered(p3107=0, p16000=9, p16001=10), min_dp=MIN_DP) == [
        3107,
        16000,
    ]
    # A position mutserve's raw table lacks is no coverage at all.
    assert uncovered({1: 100}, min_dp=MIN_DP) == list(range(2, RCRS_LENGTH + 1))


@pytest.mark.parametrize(
    ("quality", "validity"),
    [
        (1.0, Validity.STRONG),
        (0.9, Validity.STRONG),
        (0.8999, Validity.MODERATE),
        (0.8, Validity.MODERATE),
        (0.7999, Validity.LIMITED),
    ],
)
def test_quality_to_validity(quality: float, validity: Validity) -> None:
    assert haplogroup_validity(quality) is validity


def test_output_parsed(tmp_path: Path) -> None:
    path = tmp_path / "haplogroups.txt"
    path.write_text(OUTPUT)
    assert read_haplogroup(path) == Haplogroup("J1c", 1.0)
    for broken in (OUTPUT.splitlines()[0] + "\n", OUTPUT + OUTPUT.splitlines()[1]):
        path.write_text(broken)
        with pytest.raises(ValueError, match="one sample"):
            read_haplogroup(path)


def test_run_haplogrep_with_the_baked_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeHaplogrep()
    monkeypatch.setattr(haplogroup, "run_container", fake)
    vcf = write_chrm_vcf(
        [MitoCall(73, "A", "G", 1.0, 2000)],
        _covered(),
        tmp_path / "in" / "chrM.vcf",
        min_homoplasmy=0.9,
        min_dp=MIN_DP,
    )
    out = tmp_path / "out"
    out.mkdir()
    assert run_haplogrep(vcf, out, tree=TREE) == Haplogroup("J1c", 1.0)
    assert fake.args[:4] == [
        "java",
        "-jar",
        "/opt/haplogrep3/haplogrep3.jar",
        "classify",
    ]
    assert fake.args[fake.args.index("--tree") + 1] == TREE
    assert "chrM\t73\t.\tA\tG" in fake.vcf


def test_the_baked_tree_is_the_one_recorded() -> None:
    """`TREE`, which `run.tools` records, is the tree the image bakes in and its
    configuration names: an image built with another tree is another image
    (`containers.BUILT`), and this test fails until `TREE` names it."""
    tree_id, version = TREE.split("@")
    dockerfile = (DOCKER_DIR / "haplogrep3.Dockerfile").read_text()
    assert re.search(
        rf"/releases/download/{re.escape(version)}/{re.escape(tree_id)}-"
        rf"{re.escape(version)}\.zip",
        dockerfile,
    )
    assert f"haplogrep3-{haplogroup.HAPLOGREP_VERSION}-linux.tar.gz" in dockerfile
    config = (DOCKER_DIR / "haplogrep3.yaml").read_text()
    assert re.findall(r"^  - (\S+)$", config, re.MULTILINE) == [
        f"/opt/haplogrep3/trees/{tree_id}/{version}/tree.yaml"
    ]


# --- the statement -----------------------------------------------------------------------


def _store(tmp_path: Path, source: MitoSource) -> tuple[Any, str]:
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {"haplogrep3": haplogroup.HAPLOGREP_VERSION})
    ctx = RunContext(
        run_id=run_id,
        lists=load_gene_lists(),
        grading=load_grading(),
        curations={},
        gene_level={},
        faf=FafTable.open(write_faf_table(tmp_path / "faf")),
        ba1_exceptions=frozenset(),
        inheritance={},
        male=True,
        haploid=lambda chrom, pos: False,
        source_version="20261011",
    )
    with transaction(con):
        write_collected(
            StoreWriter(con, tmp_path / "work"), con, run_id, [source.collect(ctx)]
        )
    return con, run_id


def test_haplogroup_stated_under_ancestry(tmp_path: Path) -> None:
    """A gene-level finding (neither vid nor gene) under ancestry: the haplogroup, its
    quality and tree, and the positions Haplogrep3 read as rCRS for want of reads."""
    con, run_id = _store(tmp_path, MitoSource(Haplogroup("J1c", 0.85), [3107], MIN_DP))
    (finding,) = query(con, run_id, Filters()).to_json()["findings"]
    assert (finding["id"], finding["kind"], finding["topic"]) == (
        "haplogrep3:haplogroup",
        "finding",
        "ancestry",
    )
    assert finding["gene"] is None and finding["marker"]["grch38"] is None
    assert finding["claim"] == {
        "schema": "genomelab.claim/1",
        "type": "mito_haplogroup",
        "haplogroup": "J1c",
        "quality": 0.85,
        "tree": TREE,
        "uncovered_positions": [3107],
    }
    assert finding["validity"] == {
        "level": "moderate",
        "source_level": "Haplogrep3 quality 0.85",
        "reason": "1 chrM position below 10 reads read as rCRS",
    }
    assert finding["sources"] == [
        {"source": "haplogrep3", "id": TREE, "version": haplogroup.HAPLOGREP_VERSION}
    ]


def test_uncovered_chrm_unassessed(tmp_path: Path) -> None:
    """Above `MAX_UNCOVERED` of chrM below min_dp no haplogroup is stated: an unassessed
    claim with neither vid nor gene, its reason the coverage (schema 6), counted so."""
    holes = list(range(1, 501))
    con, run_id = _store(tmp_path, MitoSource(None, holes, MIN_DP))
    result = query(con, run_id, Filters(min_validity=Validity.LIMITED))
    reason = "500 of 16569 chrM positions below 10 reads; no haplogroup stated"
    assert result.findings == []
    assert result.unassessed["by_gene"][""]["reasons"] == {reason: 1}
    row = con.execute(
        "SELECT kind, vid, gene, claim FROM statement WHERE run_id = ?", [run_id]
    ).fetchone()
    assert row is not None and row[:3] == ("unassessed", None, None)
    claim = json.loads(row[3])
    assert (claim["haplogroup"], claim["uncovered_positions"]) == (None, None)
    (accounting,) = con.execute(
        "SELECT input_rows, findings, unassessed FROM accounting WHERE run_id = ?",
        [run_id],
    ).fetchall()
    assert accounting == (1, 0, 1)
