"""pgsc_calc on our calls (`genokeeper.prs.pgsc`): the panel's ancestry positions, the
target VCF, the config, the run and its containment check, and its results."""

import gzip
import io
import json
import shutil
import tarfile
from collections.abc import Sequence
from compression import zstd
from dataclasses import replace
from pathlib import Path
from typing import Any

import duckdb
import pysam
import pytest
from fixtures.faf_mini import write_faf_table
from fixtures.pgsc_mini import CONTAINED, PINNED, write_results, write_task

from genokeeper import tools
from genokeeper.config import load_grading
from genokeeper.containers import IMAGES, PLATFORMS
from genokeeper.evidence.frequency import FafTable
from genokeeper.genelists import load_gene_lists
from genokeeper.genotype.gvcf import OTHER_ALLELE_WARNING as OTHER
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import Reliability
from genokeeper.panel.build import Claim
from genokeeper.phasing import Phase
from genokeeper.prs import pgsc
from genokeeper.prs.pgsc import APOE_E2, APOE_E4, ScoreResult
from genokeeper.prs.scores import load_not_covered, load_scores
from genokeeper.query import Filters, query
from genokeeper.sources import RunContext, write_collected
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.sources.prs import DECILE_REASON, WEAK_REASON, PrsSource, decile
from genokeeper.store.db import (
    Level,
    RowType,
    Validity,
    connect,
    start_run,
    transaction,
)
from genokeeper.store.writer import StoreWriter
from genokeeper.variants.key import make_vid


def _tar_zst(path: Path, members: dict[str, bytes]) -> Path:
    with tarfile.open(path, mode="w:zst") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


def test_ancestry_positions_filtered_from_the_panels_grch38_members(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only the panel's GRCh38 genotypes are extracted, and plink2 in pgsc_calc's own image
    keeps the unrelated samples' autosomal biallelic SNPs at 5% or more, the reference
    filter pgsc_calc's ancestry step applies (its FILTER_VARIANTS defaults)."""
    prefix = "GRCh38_HGDP+1kGP_ALL"
    panel = _tar_zst(
        tmp_path / pgsc.PANEL_NAME,
        {
            "meta.txt": b"v0.1",
            "GRCh37_HGDP+1kGP_ALL.pgen": b"37",
            f"{prefix}.pgen": b"pgen",
            f"{prefix}.psam": b"psam",
            f"{prefix}.pvar.zst": b"pvar",
            "GRCh38_HGDP+1kGP.king.cutoff.out.id": b"king",
        },
    )
    seen: dict[str, Any] = {}

    def fake(
        name: str,
        args: Sequence[str],
        *,
        mounts: Sequence[tuple[Path, str, bool]],
        platform: str | None = None,
        label: str | None = None,
    ) -> None:
        (work, _, read_only), (out, _, _) = mounts
        seen.update(
            name=name,
            args=list(args),
            files=sorted(p.name for p in work.iterdir()),
            read_only=read_only,
            platform=platform,
            label=label,
        )
        (out / "GRCh38_HGDP+1kGP_ancestry.pvar.zst").write_bytes(b"kept")

    monkeypatch.setattr(pgsc, "run_container", fake)
    out = pgsc.ancestry_positions(panel, tmp_path / "ancestry")
    assert out.read_bytes() == b"kept"
    assert (out.parent / pgsc.PANEL_PVAR_NAME).read_bytes() == b"pvar"  # every variant
    assert seen["files"] == sorted(
        [f"{prefix}.pgen", f"{prefix}.psam", f"{prefix}.pvar.zst"]
        + ["GRCh38_HGDP+1kGP.king.cutoff.out.id"]
    )
    assert (seen["name"], seen["read_only"], seen["label"]) == (
        "pgsc_plink2",
        True,
        "dna:prs",
    )
    assert seen["platform"] == PLATFORMS["pgsc_plink2"] == "linux/arm64"
    args = " ".join(seen["args"])
    for option in (
        f"--pfile /panel/{prefix} vzs",
        "--remove /panel/GRCh38_HGDP+1kGP.king.cutoff.out.id",
        "--autosome",
        "--max-alleles 2",
        "--snps-only just-acgt",
        "--maf 0.05",
        "--make-just-pvar zs",
    ):
        assert option in args
    assert "pgsc_plink2" in IMAGES
    # Nothing else of the extraction is left beside the positions.
    assert sorted(p.name for p in out.parent.iterdir()) == sorted(
        [out.name, pgsc.PANEL_PVAR_NAME]
    )


# --- the target VCF --------------------------------------------------------------------------


def _called(
    dosage: str, gt: str, row_type: RowType = RowType.PASS, *warnings: str
) -> Called:
    return Called(
        gt=gt,
        dosage=dosage,
        row_type=row_type,
        gq=40,
        dp=30,
        ad=None,
        vaf=None,
        warnings=warnings,
    )


def test_target_vcf_has_ref_rows(tmp_path: Path) -> None:
    """Every target position is a row: a reference call is an explicit `0/0`, an unknown
    or low-graded call `./.` (the rule PharmCAT's input follows), and the ALTs of one
    position share its row."""
    sites = [
        (("chr1", 100, "A", "G"), _called("1", "0/1"), Level.HIGH),
        (("chr1", 200, "C", "T"), _called("0", "0/0", RowType.REF_BLOCK), Level.HIGH),
        (("chr1", 300, "G", "A"), _called("unknown", "./.", RowType.NOCALL), None),
        (("chr1", 400, "T", "C"), _called("1", "0/1"), Level.LOW),
        (("chr2", 500, "A", "C"), _called("1", "0/1"), Level.MEDIUM),
        (
            ("chr2", 500, "A", "G"),
            _called("0", "0/1", RowType.PASS, OTHER),
            Level.MEDIUM,
        ),
    ]
    out = pgsc.write_target_vcf(sites, tmp_path / "target.vcf.gz")
    with gzip.open(out, "rt") as handle:
        lines = handle.read().splitlines()
    header = [line for line in lines if line.startswith("#")]
    assert header[0] == "##fileformat=VCFv4.2"
    assert header[-1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    rows = [line.split("\t") for line in lines if not line.startswith("#")]
    assert [(r[0], r[1], r[3], r[4], r[8], r[9]) for r in rows] == [
        ("chr1", "100", "A", "G", "GT", "0/1"),
        ("chr1", "200", "C", "T", "GT", "0/0"),
        ("chr1", "300", "G", "A", "GT", "./."),
        ("chr1", "400", "T", "C", "GT", "./."),
        ("chr2", "500", "A", "C,G", "GT", "0/1"),
    ]


def test_target_positions_take_their_reference_from_the_fasta(tmp_path: Path) -> None:
    """A scoring file names an effect and another allele, the panel a REF and an ALT: the
    one the GRCh38 FASTA has at the position is the REF. Autosomes only; an allele pair the
    FASTA has neither of is no position."""
    sequence = "ACGTACGTAC" * 10
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr1\n{sequence}\n>chrX\n{sequence}\n")
    pysam.faidx(str(fasta))

    def base(pos: int) -> str:
        return sequence[pos - 1]

    def other(pos: int, k: int = 0) -> str:
        return [b for b in "ACGT" if b != base(pos)][k]

    def flip(b: str) -> str:
        return {"A": "T", "C": "G", "G": "C", "T": "A"}[b]

    columns = "rsID\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos"
    rows = [
        f"rs1\t{base(10)}\t{other(10)}\t0.1\t1\t10",  # effect is the REF
        f"rs2\t{other(20)}\t{base(20)}\t0.1\t1\t20",  # other is the REF
        f"rs3\t{other(30)}\t{other(30, 2)}\t0.1\t1\t30",  # neither, A/T: dropped
        # Reported on the reverse strand: its complement is the REF; pgsc_calc flips it.
        f"rs6\t{flip(other(30))}\t{flip(base(30))}\t0.1\t1\t30",
        f"rs4\t{base(40)}\t{other(40)}\t0.1\tX\t40",  # not an autosome
        f"rs5\t{base(50)}\t{other(50)}\t0.1\t\t",  # not harmonised
    ]
    scoring = tmp_path / "PGS000001_hmPOS_GRCh38.txt.gz"
    scoring.write_bytes(
        gzip.compress(f"#pgs_id=PGS000001\n{columns}\n{chr(10).join(rows)}\n".encode())
    )
    inferred = tmp_path / "PGS000002_hmPOS_GRCh38.txt.gz"
    inferred.write_bytes(
        gzip.compress(
            "effect_allele\thm_chr\thm_pos\thm_inferOtherAllele\n"
            f"{other(60)}\t1\t60\t{other(60, 1)}/{base(60)}\n".encode()
        )
    )
    pvar = tmp_path / "ancestry.pvar.zst"
    pvar.write_bytes(
        zstd.compress(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n"
            f"1\t20\tv1\t{base(20)}\t{other(20)}\n"  # also a score's: one position
            f"1\t70\tv2\t{base(70)}\t{other(70, 2)}\n"
            f"1\t80\tv3\t{other(80)}\t{base(80)}\n"  # the FASTA's base as ALT
            f"1\t90\tv4\t{other(90)}\t{other(90, 1)}\n".encode()  # neither
        )
    )
    assert next(pgsc.read_pvar(pvar)) == ("chr1", 20, base(20), other(20))
    assert pgsc.target_positions([scoring, inferred], pvar, fasta=fasta) == [
        ("chr1", 10, base(10), other(10)),
        ("chr1", 20, base(20), other(20)),
        ("chr1", 30, base(30), other(30)),
        ("chr1", 60, base(60), other(60)),
        ("chr1", 70, base(70), other(70, 2)),
        ("chr1", 80, base(80), other(80)),
    ]


# --- pgsc_calc's run and its containment -----------------------------------------------------


def test_config_keeps_arm_options_and_adds_network_none() -> None:
    """Our config is the arm profile's `docker.runOptions` plus `--network none`, and pins
    each of pgsc_calc's container labels (conf/modules.config) to its digest."""
    text = pgsc.nextflow_config(PINNED)
    assert (
        "docker.runOptions = '-u $(id -u):$(id -g) --platform=linux/arm64 --network none'"
        in text.splitlines()
    )
    assert set(pgsc.PGSC_LABELS) == {
        "pgscatalog_utils",
        "plink2",
        "zstd",
        "report",
        "pyyaml",
        "fraposa",
    }
    for label, ref in PINNED.items():
        assert f"    withLabel: {label} {{\n        container = '{ref}'\n    }}" in text
    with pytest.raises(ValueError, match="plink2"):
        pgsc.nextflow_config(PINNED | {"plink2": "ghcr.io/pgscatalog/plink2:2.00a5.10"})


def test_containment_refuses_unpinned_task(tmp_path: Path) -> None:
    """Every task Nextflow ran must have launched its container with `--network none` and a
    pinned digest, or the run's results are refused; a work dir with no task, a task with
    no launcher or no docker line verifies nothing and is refused too."""
    work = tmp_path / "work"
    write_task(work, PINNED["plink2"], n=1)
    write_task(work, PINNED["pgscatalog_utils"], n=2)
    (work / "ancestry" / "ref_extracted").mkdir(parents=True)  # a storeDir, not a task
    pgsc.check_containment(work)
    other = "ghcr.io/pgscatalog/plink2@sha256:" + "0" * 64
    for options, image, why in (
        ("-u $(id -u):$(id -g) --platform=linux/arm64", PINNED["plink2"], "network"),
        (f"{CONTAINED} --network host", PINNED["plink2"], "network"),
        (CONTAINED, "ghcr.io/pgscatalog/plink2:2.00a5.10", "pinned"),
        (CONTAINED, other, "pinned"),
    ):
        bad = tmp_path / why / options.replace(" ", "")[-12:] / image[-8:]
        write_task(bad, image, options)
        with pytest.raises(pgsc.ContainmentRefused, match=why):
            pgsc.check_containment(bad)
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(pgsc.ContainmentRefused, match="no task"):
        pgsc.check_containment(empty)
    no_launcher = tmp_path / "no-launcher"
    (write_task(no_launcher) / ".command.run").unlink()
    with pytest.raises(pgsc.ContainmentRefused, match=r"\.command\.run"):
        pgsc.check_containment(no_launcher)
    local = tmp_path / "local"
    (write_task(local) / ".command.run").write_text(
        "#!/bin/bash\n/bin/bash -ue .command.sh\n"
    )
    with pytest.raises(pgsc.ContainmentRefused, match="docker"):
        pgsc.check_containment(local)


def test_parse_results_reads_rates_percentiles_and_the_closest_population(
    tmp_path: Path,
) -> None:
    write_results(tmp_path, {"PGS000018": 0.9, "PGS000903": 0.6})
    results = pgsc.parse_results(tmp_path)
    assert results == {
        "PGS000018": pgsc.ScoreResult("PGS000018", 90.0, 1.28, "EUR", False, 0.9),
        "PGS000903": pgsc.ScoreResult("PGS000903", 90.0, 1.28, "EUR", False, 0.6),
    }


@pytest.mark.parametrize(
    ("flag", "low"), [("True", True), ("False", False), ("", False)]
)
def test_low_confidence_read_as_pandas_writes_it(
    tmp_path: Path, flag: str, low: bool
) -> None:
    """pgscatalog-utils writes the flag with pandas: a bool, or empty when pgsc_calc gave
    no p-value threshold (v2.3.0 gives none); any other spelling is refused, never read as
    confident."""
    write_results(tmp_path / flag, {"PGS000018": 0.9}, low=flag)
    assert pgsc.parse_results(tmp_path / flag)["PGS000018"].low_confidence is low
    write_results(tmp_path / "other", {"PGS000018": 0.9}, low="true")
    with pytest.raises(pgsc.PgscFailed, match="low-confidence flag 'true'"):
        pgsc.parse_results(tmp_path / "other")


def test_run_pgsc_calc_offline_on_local_files_inside_its_step(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Nextflow runs from the step's own directory, its work dir and outdir inside it, on
    the local scoring files (never `--pgs_id`), offline, with our config and the arm
    profile, scoring every score whatever its overlap, with pgsc_calc's default for
    multi-allelic matches and NXF_TEMP created; the results are kept only after the
    containment check passed. Docker with less than 32 GiB of memory warns, naming the
    user's setting, and the run goes on."""
    lab = tmp_path / "lab"
    scores = [pgsc.scoring_file(lab, i) for i in ("PGS000018", "PGS000903")]
    for path in scores:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#pgs_id\n")
    target = tmp_path / "step" / "input" / "target.vcf.gz"
    target.parent.mkdir(parents=True)
    target.write_bytes(b"vcf")
    out = tmp_path / "step"
    seen: dict[str, Any] = {}

    def nextflow(
        args: list[str], *, env: dict[str, str] | None = None, cwd: Path | None = None
    ) -> str:
        seen.update(args=args, env=env, cwd=cwd)
        assert env is not None and Path(env["NXF_TEMP"]).is_dir()
        write_task(Path(args[args.index("-w") + 1]), seen["image"])
        rates = {"PGS000018": 0.9, "PGS000903": 0.9}
        write_results(Path(args[args.index("--outdir") + 1]), rates)
        return ""

    monkeypatch.setattr(pgsc, "run_tool", nextflow)
    gate = tmp_path / "coreslot-gate"
    gate.touch()
    monkeypatch.setattr(tools, "GATE", gate)
    monkeypatch.setattr(pgsc, "docker_memory", lambda: int(31.9 * 2**30))
    seen["image"] = PINNED["plink2"]
    results = pgsc.run_pgsc_calc(target, scores, pgsc.panel_path(lab), out, lab=lab)
    warned = capsys.readouterr().err
    assert "31.9 GiB" in warned and "Docker Desktop (Settings, Resources)" in warned
    assert set(pgsc.parse_results(results)) == {"PGS000018", "PGS000903"}
    args = seen["args"]
    assert args[args.index("--") + 1 :][:3] == [
        str(pgsc.nextflow_path(lab)),
        "run",
        str(pgsc.pgsc_calc_dir(lab) / "main.nf"),
    ]
    assert args[: args.index("--")] == [
        str(gate),
        pgsc.LABEL,
    ]  # one gate, the whole run
    joined = " ".join(args)
    assert "-profile docker,arm" in joined and "--target_build GRCh38" in joined
    assert f"--run_ancestry {pgsc.panel_path(lab)}" in joined
    assert "--min_overlap 0" in joined and "--keep_multiallelic" not in joined
    assert "--pgs_id" not in joined and "--efo_id" not in joined
    assert seen["env"]["NXF_OFFLINE"] == "true"
    assert seen["env"]["NXF_HOME"] == str(pgsc.nxf_home(lab))
    assert seen["env"]["NXF_JAVA_HOME"] == str(pgsc.java_home(lab))
    # Everything the run writes is inside its step's directory.
    assert seen["cwd"] == out and Path(seen["env"]["NXF_TEMP"]).is_relative_to(out)
    for option in ("-w", "--outdir", "-c", "--input", "--scorefile"):
        assert Path(args[args.index(option) + 1]).is_relative_to(out), option
    # A task launched without the pinned image: the results are refused.
    seen["image"] = "ghcr.io/pgscatalog/plink2:2.00a5.10"
    other = tmp_path / "other"
    other.mkdir()
    with pytest.raises(pgsc.ContainmentRefused):
        pgsc.run_pgsc_calc(target, scores, pgsc.panel_path(lab), other, lab=lab)


# --- the source: a finding per score, APOE ---------------------------------------------------

CAD = next(s for s in load_scores() if s.pgs_id == "PGS000018")  # EUR, AFR, AMR, MID
SINGLE = "criteria_provided,_single_submitter"


def _context(tmp_path: Path, run_id: str) -> RunContext:
    return RunContext(
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


def _store(tmp_path: Path, source: PrsSource) -> tuple[duckdb.DuckDBPyConnection, str]:
    con = connect(tmp_path / "genomelab.duckdb")
    run_id = start_run(con, {}, {"pgsc_calc": pgsc.PGSC_CALC_VERSION})
    collected = source.collect(_context(tmp_path, run_id))
    with transaction(con):
        write_collected(StoreWriter(con, tmp_path / "work"), con, run_id, [collected])
    return con, run_id


def _statement(con: duckdb.DuckDBPyConnection, statement_id: str) -> dict[str, Any]:
    cursor = con.execute(
        "SELECT * FROM statement WHERE statement_id = ?", [statement_id]
    )
    row = cursor.fetchone()
    assert row is not None, statement_id
    return dict(zip([c[0] for c in cursor.description], row, strict=True))


def _result(rate: float, population: str = "EUR", low: bool = False) -> ScoreResult:
    # pgsc_calc runs with --min_overlap 0: a percentile whatever the rate.
    return ScoreResult(CAD.pgs_id, 90.0, 1.28, population, low, rate)


def test_low_overlap_unassessed(tmp_path: Path) -> None:
    """A score whose variants our target matches below 0.75 is unassessed with its rate:
    neither vid nor gene, its reason the overlap (schema 6), counted under that reason."""
    source = PrsSource(scores=[CAD], results={CAD.pgs_id: _result(0.6)}, apoe=())
    con, run_id = _store(tmp_path, source)
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    assert (row["kind"], row["vid"], row["gene"]) == ("unassessed", None, None)
    assert row["reason"] == "score overlap 0.60, below 0.75"
    assert json.loads(row["effect"])["score"] == {
        "percentile": None,
        "z": None,
        "population": "EUR",
        "match_rate": 0.6,
        "z_norm2": None,
        "decile": None,
    }
    result = query(con, run_id, Filters(min_validity=Validity.LIMITED))
    assert result.unassessed["by_gene"][""]["reasons"] == {
        "score overlap 0.60, below 0.75": 1
    }
    # The APOE result needs both ClinVar claims: none here, one input row dropped.
    (accounting,) = con.execute(
        "SELECT input_rows, unassessed, dropped FROM accounting WHERE tool = 'pgsc_calc'"
    ).fetchall()
    assert accounting[:2] == (2, 1)
    assert json.loads(accounting[2]) == {"apoe_no_clinvar_claim": 1}


def test_no_percentile_at_a_passing_rate_says_so(tmp_path: Path) -> None:
    result = ScoreResult(CAD.pgs_id, None, None, "EUR", False, 0.82)
    source = PrsSource(scores=[CAD], results={CAD.pgs_id: result}, apoe=())
    con, _ = _store(tmp_path, source)
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    assert (row["kind"], row["reason"]) == (
        "unassessed",
        "pgsc_calc gave no percentile",
    )


@pytest.mark.parametrize(
    ("population", "low", "validity", "reason"),
    [
        ("EUR", False, "moderate", None),
        ("MID", False, "moderate", None),
        ("EAS", False, "limited", "no PGS Catalog evaluation in EAS"),
        ("EUR", True, "limited", "closest population EUR assigned with low confidence"),
    ],
)
def test_validity_by_ancestry_match(
    tmp_path: Path, population: str, low: bool, validity: str, reason: str | None
) -> None:
    """Moderate with an evaluation in the population pgsc_calc finds the user closest to,
    limited otherwise; the percentile and Z are against that population."""
    result = _result(0.9, population, low)
    source = PrsSource(scores=[CAD], results={CAD.pgs_id: result}, apoe=())
    con, _ = _store(tmp_path, source)
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    assert (row["kind"], row["topic"], row["validity"], row["reason"]) == (
        "finding",
        "risks",
        validity,
        "; ".join(r for r in (reason, DECILE_REASON) if r),
    )
    assert json.loads(row["effect"]) == {
        "relative": None,
        "absolute": None,
        "score": {
            "percentile": 90.0,
            "z": 1.28,
            "population": population,
            "match_rate": 0.9,
            "z_norm2": None,
            "decile": 10,
        },
    }
    claim = json.loads(row["claim"])
    assert (claim["type"], claim["pgs_id"], claim["condition"]) == (
        "pgs",
        "PGS000018",
        "coronary artery disease",
    )
    assert claim["evaluations"] == [list(e) for e in CAD.evaluations]
    assert claim["weak"] is False


def test_weak_score_carries_the_reason(tmp_path: Path) -> None:
    """A weak score keeps its validity and kind; only its reason (after the decile's) and
    claim say weak."""
    weak = replace(CAD, weak=True)
    source = PrsSource(scores=[weak], results={CAD.pgs_id: _result(0.9)}, apoe=())
    con, _ = _store(tmp_path, source)
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    assert (row["kind"], row["validity"], row["reason"]) == (
        "finding",
        "moderate",
        f"{DECILE_REASON}; {WEAK_REASON}",
    )
    assert json.loads(row["claim"])["weak"] is True


def test_weak_score_below_overlap_keeps_both_reasons(tmp_path: Path) -> None:
    weak = replace(CAD, weak=True)
    source = PrsSource(scores=[weak], results={CAD.pgs_id: _result(0.6)}, apoe=())
    con, run_id = _store(tmp_path, source)
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    reason = "score overlap 0.60, below 0.75; " + WEAK_REASON
    assert (row["kind"], row["reason"]) == ("unassessed", reason)
    result = query(con, run_id, Filters(min_validity=Validity.LIMITED))
    assert result.findings == []
    assert result.unassessed["by_gene"][""]["reasons"] == {reason: 1}


MIGRAINE = next(n for n in load_not_covered() if n.condition == "migraine")


def test_not_covered_statement(tmp_path: Path) -> None:
    """A condition no score qualifies for is one unassessed statement saying so and why,
    with neither vid, gene, percentile nor effect; the accounting counts it."""
    source = PrsSource(scores=[], results={}, apoe=(), not_covered=[MIGRAINE])
    con, _ = _store(tmp_path, source)
    row = _statement(con, "pgsc_calc:not-covered:migraine")
    assert {k: row[k] for k in ("kind", "vid", "gene", "validity", "topic")} == {
        "kind": "unassessed",
        "vid": None,
        "gene": None,
        "validity": "limited",
        "topic": "risks",
    }
    assert row["reason"] == (
        "not covered: no score meets the criteria (evaluated only in its own cohort)"
    )
    assert (row["source_level"], row["effect"], row["actionable"]) == (
        "no qualifying PGS Catalog score",
        None,
        False,
    )
    assert (row["source"], row["source_id"], row["source_version"]) == (
        "pgsc_calc",
        "migraine",
        pgsc.PGSC_CALC_VERSION,
    )
    assert json.loads(row["claim"]) == {
        "schema": CLAIM_SCHEMA,
        "type": "pgs_not_covered",
        "condition": "migraine",
    }
    (accounting,) = con.execute(
        "SELECT input_rows, findings, unassessed FROM accounting"
        " WHERE tool = 'pgsc_calc'"
    ).fetchall()
    assert accounting == (2, 0, 1)  # the condition and APOE, which is dropped


def test_not_covered_found_by_condition(tmp_path: Path) -> None:
    source = PrsSource(
        scores=[CAD],
        results={CAD.pgs_id: _result(0.9)},
        apoe=(),
        not_covered=load_not_covered(),
    )
    con, run_id = _store(tmp_path, source)
    result = query(
        con,
        run_id,
        Filters(min_validity=Validity.LIMITED, text="migraine"),
        list_unassessed=True,
    )
    assert result.findings == []
    assert result.unassessed_sites is not None
    assert [(s["id"], s["validity"]["reason"]) for s in result.unassessed_sites] == [
        (
            "pgsc_calc:not-covered:migraine",
            "not covered: no score meets the criteria (evaluated only in its own cohort)",
        )
    ]
    assert result.unassessed["by_gene"][""]["reasons"] == {
        "not covered: no score meets the criteria (evaluated only in its own cohort)": 1
    }


def _apoe_claim(
    allele: tuple[str, int, str, str], variation_id: int, stars: str
) -> Claim:
    chrom, pos, ref, alt = allele
    return Claim(
        variation_id=variation_id,
        gene="APOE",
        clnsig="drug_response",
        clnrevstat=stars,
        disease_ids="",
        condition="APOE allele",
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        vid=make_vid(chrom, pos, ref, alt),
        rsid=None,
    )


def _graded(
    claim: Claim, dosage: str, level: Level
) -> tuple[Claim, Called, Reliability]:
    gt = {"0": "0/0", "1": "0/1", "2": "1/1"}.get(dosage, "./.")
    return claim, _called(dosage, gt), Reliability(None, level, (), (), None)


def test_apoe_e3e4(tmp_path: Path) -> None:
    """APOE from rs429358 (C, ε4) and rs7412 (T, ε2): one copy of rs429358's C alone is
    ε3/ε4, a gene-level finding under risks with its gene call at the weaker call's level;
    its validity is the weaker of the two SNPs' ClinVar claims, each SNP's strongest."""
    assert pgsc.apoe_genotype({APOE_E4: "1", APOE_E2: "0"}) == "ε3/ε4"
    assert pgsc.apoe_genotype({APOE_E4: "0", APOE_E2: "0"}) == "ε3/ε3"
    assert pgsc.apoe_genotype({APOE_E4: "2", APOE_E2: "0"}) == "ε4/ε4"
    assert pgsc.apoe_genotype({APOE_E4: "0", APOE_E2: "1"}) == "ε2/ε3"
    assert pgsc.apoe_genotype({APOE_E4: "unknown", APOE_E2: "0"}) is None
    # A double heterozygote is ε2/ε4 in trans, ε1/ε3 in cis, unknown unphased.
    double = {APOE_E4: "1", APOE_E2: "1"}
    assert pgsc.apoe_genotype(double, Phase.IN_TRANS) == "ε2/ε4"
    assert pgsc.apoe_genotype(double, Phase.IN_CIS) == "ε1/ε3"
    assert pgsc.apoe_genotype(double, Phase.UNPHASED) is None
    assert pgsc.apoe_genotype(double) is None
    strong = "criteria_provided,_multiple_submitters,_no_conflicts"  # 2★ drug_response
    apoe = (
        _graded(_apoe_claim(APOE_E4, 17864, strong), "1", Level.HIGH),
        _graded(_apoe_claim(APOE_E4, 17865, ""), "1", Level.HIGH),  # 0★: limited
        _graded(_apoe_claim(APOE_E2, 17848, SINGLE), "0", Level.MEDIUM),  # 1★
    )
    con, _ = _store(tmp_path, PrsSource(scores=[], results={}, apoe=apoe))
    row = _statement(con, "apoe:APOE")
    assert (row["kind"], row["gene"], row["vid"], row["topic"]) == (
        "finding",
        "APOE",
        None,
        "risks",
    )
    assert row["validity"] == "moderate"  # rs429358's strongest (strong), rs7412's
    assert json.loads(row["claim"])["diplotype"] == "ε3/ε4"
    gene_call = con.execute(
        "SELECT diplotype, level FROM gene_call WHERE gene = 'APOE'"
    ).fetchone()
    assert gene_call == ("ε3/ε4", "medium")
    # Unphased, the double heterozygote is unassessed and says why; no gene call.
    con.close()
    shutil.rmtree(tmp_path / "work", ignore_errors=True)
    (tmp_path / "genomelab.duckdb").unlink()
    unphased = (apoe[0], _graded(apoe[2][0], "1", Level.MEDIUM))
    con, _ = _store(tmp_path, PrsSource(scores=[], results={}, apoe=unphased))
    row = _statement(con, "apoe:APOE")
    assert (row["kind"], row["reason"]) == (
        "unassessed",
        "ε2/ε4 or ε1/ε3: phase unknown",
    )
    assert con.execute("SELECT count(*) FROM gene_call").fetchone() == (0,)


# --- slice 3f: the population similarity, and Z_norm2 ---------------------------------------

HG002 = Path(__file__).parent / "fixtures" / "hg002"


def _hg002_results(results: Path) -> Path:
    """HG002's pgsc_calc tables (Task 5's fixtures, sampleset `genomelab`) where
    `parse_results` and `parse_similarity` look, with a match summary for its one score."""
    score = results / pgsc.SAMPLESET / "score"
    match = results / pgsc.SAMPLESET / "match"
    score.mkdir(parents=True)
    match.mkdir(parents=True)
    for name in ("pgs", "popsimilarity"):
        text = (HG002 / f"hg002.{name}.txt").read_bytes()
        (score / f"genomelab_{name}.txt.gz").write_bytes(gzip.compress(text))
    (match / "genomelab_summary.csv").write_text(
        "dataset,accession,match_status,count\n"
        "genomelab,PGS000903_hmPOS_GRCh38,matched,1506\n"
        "genomelab,PGS000903_hmPOS_GRCh38,unmatched,299\n"
    )
    return results


def test_parse_similarity_hg002(tmp_path: Path) -> None:
    """The sampleset's row of the population similarity: the closest population, its
    low-confidence flag, and every `RF_P_<POP>` as a probability."""
    similarity = pgsc.parse_similarity(_hg002_results(tmp_path))
    assert similarity == pgsc.PopSimilarity(
        population="EUR",
        low_confidence=False,
        probabilities={
            "AFR": 0.0,
            "AMR": 0.12,
            "CSA": 0.01,
            "EAS": 0.01,
            "EUR": 0.65,
            "MID": 0.21,
        },
    )


def test_z_norm2_parsed_hg002(tmp_path: Path) -> None:
    """pgsc_calc's default normalisation writes `Z_norm2`, read by its column name."""
    (result,) = pgsc.parse_results(_hg002_results(tmp_path)).values()
    assert (result.pgs_id, result.percentile) == ("PGS000903", 31.18440779610195)
    assert result.z_norm2 == pytest.approx(-0.6061704546590493)
    assert result.match_rate == pytest.approx(1506 / 1805)


def test_z_norm2_absent_is_none(tmp_path: Path) -> None:
    write_results(tmp_path, {"PGS000018": 0.9})
    assert pgsc.parse_results(tmp_path)["PGS000018"].z_norm2 is None


@pytest.mark.parametrize(
    ("percentile", "expected"),
    [(0.0, 1), (9.99, 1), (10.0, 2), (50.0, 6), (99.9, 10), (100.0, 10)],
)
def test_decile_bounds(percentile: float, expected: int) -> None:
    """Review Focus 6: a percentile of 100.0 or 0.0 is decile 10 or 1, never 11 or 0."""
    assert decile(percentile) == expected


def test_scored_pgs_shows_decile_and_z_norm2(tmp_path: Path) -> None:
    """A scored statement's effect carries the decile and `Z_norm2` beside the percentile
    (kept for the change feed), and its reason says why deciles are shown."""
    result = ScoreResult(CAD.pgs_id, 91.5, 1.37, "EUR", False, 0.9, z_norm2=1.2)
    con, _ = _store(tmp_path, PrsSource([CAD], {CAD.pgs_id: result}, apoe=()))
    row = _statement(con, f"pgsc_calc:{CAD.pgs_id}")
    score = json.loads(row["effect"])["score"]
    assert (score["percentile"], score["decile"], score["z_norm2"]) == (91.5, 10, 1.2)
    assert row["reason"] == DECILE_REASON
    assert DECILE_REASON == (
        "shown in deciles: one genome's score is less precise than its percentile"
        " (Ding et al., Nat Genet 2022)"
    )
