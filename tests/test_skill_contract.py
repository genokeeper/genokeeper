"""The genome-ask skill documents every key of the query's JSON contract.

A key counts as documented only as a bullet head inside the skill's "The JSON contract"
section, nested under the bullet of its parent key: "- `key` — ..." at the top and
"  - `key`: ..." below it. A word in backticks anywhere else documents nothing.
"""

import argparse
import ast
import json
import re
from dataclasses import replace
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

import pytest
from fixtures import snapshots_mini as sm
from fixtures.contract import CHANGES_CONTRACT, CONTRACT, GENOTYPE_CONTRACT, key_paths
from fixtures.run_context import run_context
from fixtures.store_mini import Row, seed, seed_coverage

from genokeeper import cli
from genokeeper.ancestry.archaic import ArchaicShare
from genokeeper.bloodgroups import SystemCall
from genokeeper.calibrate.tables import NOT_CALIBRATED as DOSAGE_NOT_CALIBRATED
from genokeeper.carrier import genotype
from genokeeper.changes import report
from genokeeper.config import load_grading
from genokeeper.genotype.gvcf import Called
from genokeeper.genotype.reliability import AF_CLONAL, AF_LOW, Reliability
from genokeeper.hhv6 import Hhv6Call
from genokeeper.history import Dates
from genokeeper.mito import MitoCall
from genokeeper.pgx.hla import read_optitype, read_t1k
from genokeeper.pgx.hla_typing import TAG_RSIDS
from genokeeper.pgx.pharmcat import DOSING
from genokeeper.prs.pgsc import PopSimilarity
from genokeeper.query import Filters, query, show
from genokeeper.refresh_log import (
    RefreshRecord,
    SourceOutcome,
    SourceStatus,
    append,
    publications_path,
)
from genokeeper.repeats.loci import Allele, RepeatCall, Rfc1Call
from genokeeper.repeats.not_covered import NOT_COVERED
from genokeeper.roh import RohSummary, Run
from genokeeper.sex import Dosage
from genokeeper.smn import ParascopyCall, SmnCaller
from genokeeper.sources.ancestry import AncestrySource
from genokeeper.sources.bloodgroups import SECRETOR_TRAIT, BloodGroupSource
from genokeeper.sources.dosage import STATEMENT_ID as DOSAGE_STATEMENT
from genokeeper.sources.dosage import DosageSource
from genokeeper.sources.hhv6 import Hhv6Source
from genokeeper.sources.hla import HlaSource
from genokeeper.sources.mitomap import MitomapSource
from genokeeper.sources.repeats import RepeatSource
from genokeeper.sources.roh import RohSource
from genokeeper.sources.smn import SmnSource
from genokeeper.sources.traits import TRAIT_ALLELES, TraitSource, load_traits
from genokeeper.sources.yleaf import YleafSource
from genokeeper.store.db import Kind, Level, RowType, Validity, connect
from genokeeper.yhaplogroup import YHaplogroup

SKILL = Path(__file__).resolve().parents[1] / ".claude/skills/genome-ask/SKILL.md"

_SECTION = "## The JSON contract"
_BULLET = re.compile(r"^(?P<indent> *)- (?P<head>.*)$")
_KEY_HEAD = re.compile(r"`(?P<key>\w+)`(?: —|:)")


def _section(text: str) -> tuple[int, list[str]]:
    """The line index where the contract section's body starts, and its lines."""
    lines = text.splitlines()
    start = lines.index(_SECTION) + 1
    end = next(
        (i for i in range(start, len(lines)) if lines[i].startswith("## ")), len(lines)
    )
    return start, lines[start:end]


def bullets(text: str) -> dict[str, int]:
    """Each key path the contract section documents as a nested bullet head -> its line."""
    offset, lines = _section(text)
    paths: dict[str, int] = {}
    stack: list[tuple[int, str | None]] = []
    for i, line in enumerate(lines, start=offset):
        bullet = _BULLET.match(line)
        if bullet is None:
            continue
        indent = len(bullet["indent"])
        while stack and stack[-1][0] >= indent:
            stack.pop()
        head = _KEY_HEAD.match(bullet["head"])
        stack.append((indent, head["key"] if head else None))
        keys = [key for _, key in stack if key is not None]
        if len(keys) == len(stack):
            paths[".".join(keys)] = i
    return paths


def test_contract_list_matches_the_query_and_show_output(tmp_path: Path) -> None:
    con = connect(tmp_path / "genomelab.duckdb")
    effect = {
        "relative": {"measure": "OR", "value": 2.0, "ci": [1.5, 2.5]},
        "absolute": {"risk": 0.2, "baseline": 0.1, "population": "EUR"},
    }
    seed(
        con,
        [
            Row("clinvar:1", effect=effect),
            Row(
                "clinvar:2",
                kind="unassessed",
                dosage="unknown",
                level="low",
                pos=1002,
                warnings=["genotype unknown: no-call"],
            ),
            Row(
                "pharmcat:CYP2C19",
                gene="CYP2C19",
                level=None,
                topic="drugs",
                gene_call={
                    "diplotype": "*1/*2",
                    "phenotype": "Intermediate Metabolizer",
                    "activity_score": None,
                    "callers": {"PharmCAT": "*1/*2"},
                    "agreed": None,
                    "level": "medium",
                    "warnings": [],
                },
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "pharmcat_gene",
                    "gene": "CYP2C19",
                    "diplotype": "*1/*2",
                    "phenotype": "Intermediate Metabolizer",
                    "activity_score": None,
                    "call_source": "MATCHER",
                    "drugs": [
                        {
                            "drug": "clopidogrel",
                            "source": "CPIC",
                            "classification": "Strong",
                        }
                    ],
                },
            ),
            Row(
                "pharmgkb:655384799",
                gene="TNF",
                topic="drugs",
                validity="moderate",  # level 1B
                pos=1004,
                rsid="rs1800629",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "pharmgkb_clinical_annotation",
                    "annotation_id": "655384799",
                    "gene": "TNF",
                    "category": "Efficacy",
                    "condition": "Psoriasis",
                    "drugs": [{"drug": "etanercept"}],
                },
            ),
            Row(
                "pgsc_calc:PGS000018",
                gene=None,
                level=None,
                validity="moderate",
                effect={
                    "relative": None,
                    "absolute": None,
                    "score": {
                        "percentile": 90.0,
                        "z": 1.28,
                        "population": "EUR",
                        "match_rate": 0.93,
                        "z_norm2": 1.1,
                        "decile": 10,
                    },
                },
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "pgs",
                    "pgs_id": "PGS000018",
                    "condition": "coronary artery disease",
                    "trait_efo": "MONDO_0005010",
                    "publication": "Inouye M et al., J Am Coll Cardiol 2018",
                    "development_ancestry": "GWAS MAE 50.9%, EUR 37%",
                    "evaluations": [["PPM000034", "EUR"]],
                    "licence": "PGS obtained from the Catalog should be cited",
                    "weak": False,
                },
            ),
            Row(
                "haplogrep3:haplogroup",
                gene=None,
                level=None,
                topic="ancestry",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "mito_haplogroup",
                    "haplogroup": "J1c",
                    "quality": 0.93,
                    "tree": "phylotree-rcrs@17.3",
                    "uncovered_positions": [3107],
                },
            ),
            Row(
                "expansionhunter:HTT",
                gene="HTT",
                level=None,
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "repeat_expansion",
                    "gene": "HTT",
                    "locus_id": "HTT",
                    "motif": "CAG",
                    "alleles": [
                        {
                            "repeats": 42,
                            "ci": [41, 43],
                            "spanning": 6,
                            "flanking": 3,
                            "inrepeat": 0,
                        }
                    ],
                    "class": "pathogenic",
                    "thresholds": {
                        "normal_max": 26,
                        "intermediate": [27, 35],
                        "reduced_penetrance": [36, 39],
                        "pathogenic_min": 40,
                    },
                    "image": "rev.HTT.svg",
                    "confirm": None,
                    "rfc1": None,
                },
            ),
            Row(
                "expansionhunter:CFTR",
                gene="CFTR",
                level=None,
                topic="carrier",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "repeat_expansion",
                    "gene": "CFTR",
                    "locus_id": "CFTR",
                    "motif": "TG/T",
                    "alleles": {
                        "tg": [
                            {
                                "repeats": 12,
                                "ci": [12, 12],
                                "spanning": 6,
                                "flanking": 0,
                                "inrepeat": 0,
                            }
                        ],
                        "t": [
                            {
                                "repeats": 5,
                                "ci": [5, 5],
                                "spanning": 6,
                                "flanking": 0,
                                "inrepeat": 0,
                            }
                        ],
                    },
                    "class": None,
                    "thresholds": None,
                    "image": None,
                    "confirm": None,
                    "rfc1": None,
                },
            ),
            Row(
                "smn:SMN1",
                gene="SMN1",
                level=None,
                topic="carrier",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "smn_copy_number",
                    "gene": "SMN1",
                    "smn1": 1,
                    "smn2": 2,
                    "g27134tg": 0,
                    "callers": {"SMNCopyNumberCaller": "SMN1 1, SMN2 2"},
                },
            ),
            Row(
                "yleaf:haplogroup",
                gene=None,
                level=None,
                topic="ancestry",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "y_haplogroup",
                    "haplogroup": "J-C134621",
                    "marker": "C134621",
                    "qc_score": 1.0,
                    "valid_markers": 170393,
                    "tree": "Yleaf 4.1.6 bundled YFull",
                },
            ),
            Row(
                "gwas:rs1801282-C:MONDO_0005148",
                gene="PPARG",
                pos=1006,
                rsid="rs1801282",
                dosage="0",
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "gwas_association",
                    "condition": "type 2 diabetes mellitus",
                    "trait_efo": "MONDO_0005148",
                    "risk_allele": "C",
                    "risk_allele_copies": 2,
                    "p_value": "2E-12",
                    "pubmed_id": "17463246",
                    "study_accession": "GCST000001",
                },
            ),
        ],
    )
    seed(
        con,
        [
            Row(
                "predicted:GRCh38:chr17:1008:C:T:ENST00000357654",
                validity="limited",
                pos=1008,
                claim={
                    "schema": "genomelab.claim/1",
                    "type": "predicted_variant",
                    "classification": "predicted, not classified",
                    "consequence": "stop_gained",
                    "transcript": "ENST00000357654",
                    "hgvsc": "ENST00000357654.9:c.1234C>T",
                    "hgvsp": "ENSP00000350283.3:p.Gln412Ter",
                },
            )
        ],
    )
    con.execute(
        "UPDATE statement SET source = 'vep' WHERE statement_id LIKE 'predicted:%'"
    )
    seed(con, [_dosage_row(tmp_path)])
    # One finding per slice 3f claim type, as its source writes the claim.
    claims = _slice3f_claims(tmp_path / "3f")
    seed(
        con,
        [
            Row(f"slice3f:{t}", gene=None, level=None, claim=claim)
            for t, claim in claims.items()
        ],
    )
    seed_coverage(con, "BRCA1", 0.99)
    day = date(2026, 10, 1)
    dates = {"clinvar:1": Dates(day, day, day, on_or_before=True)}
    run, finding = show(con, None, "clinvar:1", dates=dates)
    shown = {"run_id": run, "finding": finding}
    listed = query(con, None, Filters(), list_unassessed=True).to_json()
    assert listed["unassessed_sites"] and listed["unassessed"]["by_gene"]
    _, predicted = show(con, None, "predicted:GRCh38:chr17:1008:C:T:ENST00000357654")
    gene = query(con, None, Filters(gene="BRCA1")).to_json()
    found = key_paths(listed) | key_paths(shown) | key_paths(predicted)
    assert found | key_paths(gene) == CONTRACT


def test_every_json_key_is_documented() -> None:
    documented = set(bullets(SKILL.read_text()))
    assert sorted((CONTRACT | CHANGES_CONTRACT) - documented) == []


def test_no_bullet_documents_a_key_the_output_lacks() -> None:
    # A key the output dropped (schema 3: `marker.hg19`, `reliability.mapping`) must leave
    # the skill too, or it would read a field that is never there.
    documented = set(bullets(SKILL.read_text()))
    assert sorted(documented - CONTRACT - GENOTYPE_CONTRACT - CHANGES_CONTRACT) == []


def test_skill_documents_genotype_contract(tmp_path: Path) -> None:
    con = connect(tmp_path / "genomelab.duckdb")
    seed(
        con,
        [
            Row("clinvar:1", gene="SMN1", kind="non_finding", dosage="0"),
            Row("clinvar:2", gene="SMN1", kind="unassessed", call=False, pos=1002),
        ],
    )
    document = genotype(con, None, gene="SMN1").to_json()
    assert document["sites"] and document["gene"]["dominant_mechanism_not_assessed"]
    paths = key_paths(document)
    assert paths - CONTRACT == GENOTYPE_CONTRACT
    text = SKILL.read_text()
    assert sorted(paths - set(bullets(text))) == []
    assert "not assessed in slice 1a" not in text


@pytest.mark.parametrize(
    "path", sorted(CONTRACT | GENOTYPE_CONTRACT | CHANGES_CONTRACT)
)
def test_dropping_a_bullet_fails_the_check(path: str) -> None:
    """One bullet documents each path, and nothing else in the skill stands in for it."""
    text = SKILL.read_text()
    lines = text.splitlines()
    i = bullets(text)[path]
    without = "\n".join(lines[:i] + lines[i + 1 :])
    assert path not in set(bullets(without))


def _contract(body: str) -> str:
    return f"# skill\n\n- `p_wrong`: outside\n\n{_SECTION}\n\n{body}\n## Next\n- `gene` —\n"


def test_backticks_that_are_not_bullet_heads_document_nothing() -> None:
    body = "Keys: `findings`, `p_wrong`.\n\n- the `gene` symbol\n- `topic` is a topic\n"
    assert bullets(_contract(body)) == {}


def test_a_nested_key_is_documented_only_under_its_parent() -> None:
    body = (
        "- `validity` — the claim\n  - `level`: its level\n- `reliability` — the call\n"
    )
    assert set(bullets(_contract(body))) == {
        "validity",
        "validity.level",
        "reliability",
    }


def test_a_key_under_a_plain_bullet_is_not_rooted() -> None:
    body = "- Each finding:\n  - `gene` — its gene\n"
    assert bullets(_contract(body)) == {}


def test_p_wrong_is_documented_as_measured_on_hg003(tmp_path: Path) -> None:
    """`p_wrong` is the error rate measured on HG003 (DNBSEQ-T7) with its interval, never
    the caller's own estimate; null outside the calibrated regions."""
    con = connect(tmp_path / "genomelab.duckdb")
    seed(con, [Row("clinvar:1")])
    con.execute("UPDATE call SET p_wrong_lo = 5e-5, p_wrong_hi = 2e-4")
    (finding,) = query(con, None, Filters()).to_json()["findings"]
    assert finding["reliability"]["p_wrong_ci"] == [5e-5, 2e-4]
    lines = SKILL.read_text().splitlines()
    documented = bullets(SKILL.read_text())
    assert "measured on HG003 (DNBSEQ-T7)" in lines[documented["reliability.p_wrong"]]
    text = " ".join(SKILL.read_text().split())
    assert "caller's own estimate" not in text
    assert "`outside calibrated regions`" in text


def test_frequency_is_claimed_only_for_a_run_that_records_it(tmp_path: Path) -> None:
    con = connect(tmp_path / "genomelab.duckdb")
    seed(con, [Row("clinvar:1")])  # a run graded before the FAF table: no `gnomad_faf`
    assert "gnomad_faf" not in query(con, None, Filters()).to_json()["versions"]
    text = " ".join(SKILL.read_text().split())
    assert "Without a `gnomad_faf` key" in text and "frequency was not checked" in text
    assert "**Population frequency is assessed**" not in text


def test_skill_covers_the_traits() -> None:
    """The trait claim names exactly `load_traits()`'s ids, and the Traits bullet quotes
    each of slice 3g's six caveats (`Trait.note`) verbatim."""
    text = SKILL.read_text()
    _, lines = _section(text)
    at = next(i for i, line in enumerate(lines) if "(`trait`," in line)
    end = next(i for i in range(at, len(lines)) if lines[i].endswith(":"))
    header = " ".join(lines[at : end + 1])
    traits = {t.name: t for t in load_traits()}
    assert set(re.findall(r"`(trait:\w+)`", header)) == {f"trait:{n}" for n in traits}
    section = _section_named(text, "## Slice 3f: typing, mitochondria and what remains")
    bullet = section.split("- **Traits**")[1].split("- **")[0]
    for name in ("FUT2", "ABCC11", "MC1R", "HERC2", "AMPD1", "CCR5"):
        assert f"`{traits[name].note}`" in bullet, name


def test_skill_answers_secretor_status_from_one_statement() -> None:
    """BR-026: the routing table reads both FUT2 statements by id (the blood group has no
    gene, so `--gene FUT2` never returns it), secretor status is answered from
    `bloodgroup:FUT2` when it is a finding with a phenotype, otherwise from `trait:FUT2`,
    and the Blood groups bullet quotes FUT2's stored note."""
    text = SKILL.read_text()
    row = next(line for line in text.splitlines() if "secretor" in line.lower())
    assert row.startswith("| ")
    assert "`show bloodgroup:FUT2 --json`" in row and "`show trait:FUT2 --json`" in row
    section = _section_named(text, "## Slice 3f: typing, mitochondria and what remains")
    blood = section.split("- **Blood groups**")[1].split("- **")[0]
    assert f"`{SECRETOR_TRAIT}`" in blood
    traits = section.split("- **Traits**")[1].split("- **")[0]
    assert (
        "answer secretor status from `bloodgroup:FUT2` when it is a finding with a"
        " phenotype, otherwise from `trait:FUT2`, and say which one answered; when both"
        " answer and differ, give the blood group's and name the trait's"
    ) in traits


def _section_named(text: str, heading: str) -> str:
    lines = text.splitlines()
    start = lines.index(heading)
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")),
        len(lines),
    )
    return " ".join(" ".join(lines[start:end]).split())


def test_skill_answers_in_the_users_language() -> None:
    """The skill is not Russian-only: it answers in the language it is asked in, speaks of
    the user, and keeps the store's English text as quoted data."""
    text = SKILL.read_text()
    flat = " ".join(text.split())
    assert "**Answer in the language the user asks in**" in flat
    assert "Answers in the language the user asks in" in flat  # the description
    assert "Answer in Russian" not in flat and "Answers in Russian" not in flat
    # `owner` survives only as the key name and the reason that quotes it.
    assert re.findall(r"\bowner\b", text) == ["owner"] * len(
        re.findall(r"`owner`|<owner>", text)
    )


def test_skill_forbids_dosing() -> None:
    """The skill answers drug questions from the store (PharmCAT's gene results, PharmGKB's
    annotations, ClinVar's drug responses) and never with a dose: it says the store holds
    no dosing text, forbids stating a dose, a dose change or a drug choice, and sends the
    user to a prescriber. Drugs are no longer "not covered"."""
    text = SKILL.read_text()
    section = _section_named(text, "## Drugs")
    for needed in (
        "`--topic drugs`",
        "`pharmcat:<gene>`",
        "`pharmgkb:<annotation_id>`",
        "`drug_response`",
        "The store holds no dosing text",
        "Never state a dose, a dose change, a dosing schedule or a drug choice",
        "prescriber",
        "PharmCAT's HTML report",
    ):
        assert needed in section, needed
    flat = " ".join(text.split())
    assert '"covered in slice 1b"' not in flat
    assert "drug-response classifications are not findings until then" not in flat
    # Every line of the skill that speaks of a dose forbids one (`dosage` is the genotype's
    # copy count, not a dose).
    for line in text.splitlines():
        if DOSING.search(line.replace("dosage", "")):
            assert re.search(r"never|no dosing|holds no", line, re.IGNORECASE), line


def test_skill_labels_predicted() -> None:
    """A predicted finding (slice 1e) is labelled "predicted, not classified", would need a
    clinical lab's confirmation, is never a diagnosis, carries the geneticist note, names a
    failed gate, and comes after the ClinVar findings."""
    section = _section_named(SKILL.read_text(), "## Predicted findings")
    for needed in (
        "`predicted:<vid>:<transcript>`",
        "source `vep`",
        '"predicted, not classified"',
        "confirmation by a clinical lab",
        "Never present it as a diagnosis",
        '"discuss with a geneticist"',
        "`--min-validity limited`",
        "failed LoF gate",
        "never outranks a ClinVar-classified one",
    ):
        assert needed in section, needed


def test_skill_names_the_annotation_limits() -> None:
    """BR-027, BR-028: the predicted findings name the two genes the annotation cannot grade
    in full: SLC37A4, whose only MANE transcript is on a GRC fix patch (its LoF always fails
    gate 2 past gate 1, no reason alone to call it doubtful), and PCDH15, whose last coding
    exon SpliceAI's annotation lacks; Coverage names what SLC37A4 lacks."""
    text = SKILL.read_text()
    section = _section_named(text, "## Predicted findings")
    flat = " ".join(section.split())
    for needed in (
        "SLC37A4 is the exception to the gate-failure rule above, for a loss of function",
        "never paired with a second variant in the gene",
        "A start-loss there fails gate 2 too, but is never a finding anywhere",
        "GRCh38:chr11:119027725:AC:A",
        "it restores MANE's sequence",
        "never graded as a missense",
        "TNNT2's intron 1",
        "GRC fix patch",
        "passes gate 1 always fails LoF gate 2",
        (
            "that gate alone is no reason to call it doubtful: say the prediction cannot be"
            " judged on the primary assembly and a clinical lab can assess it"
        ),
        "PCDH15's last coding exon",
        "no SpliceAI score",
        "still graded as a loss of function",
    ):
        assert needed in flat, needed
    coverage = " ".join(_section_named(text, "## Coverage").split())
    assert "SLC37A4's coverage, deletions, insertions and copy number" in coverage


def test_skill_covers_stage3a() -> None:
    """SMN1, the eleven repeat loci and the Y haplogroup are covered, every non-normal
    result is a screening result that needs confirmation, FMR1 and AFF2 follow slice 3e's
    rule, and the Y haplogroup is no longer listed as not covered."""
    text = SKILL.read_text()
    section = _section_named(text, "## Coverage")
    for needed in (
        "`smn:SMN1`",
        "`expansionhunter:<gene>`",
        "`yleaf:haplogroup`",
        "screening result, needs confirmation",
        "FMR1 and AFF2 are assessed from reads spanning every X allele",
        "stage 3e adds loci beyond 3a's eleven",
    ):
        assert needed in section, needed
    assert "FMR1 is not covered" not in text
    assert "the Y-chromosome haplogroup" not in text


def test_skill_reads_a_repeat_class_from_its_interval() -> None:
    """A repeat claim's `class` is the class its interval's upper bound reaches (BR-009),
    so a watch item's may lie above its estimate's, even `pathogenic` on a borderline one:
    it is worded as the class the interval may reach, never as the allele's."""
    text = " ".join(SKILL.read_text().split())
    assert (
        "`class`: the class the longest allele's interval reaches at its upper bound"
        in text
    )
    assert (
        "word it as the class the interval may reach, never as the allele's class"
        in text
    )
    assert "`class`: the longest allele's class" not in text


def test_skill_reads_a_stage3_finding_as_a_possible_carrier() -> None:
    """A stage-3 gene-level finding answers `possible_carrier` at any validity (BR-009)."""
    text = " ".join(SKILL.read_text().split())
    assert "for a stage-3 gene-level finding at any validity" in text
    assert (
        "`possible_carrier` for a finding, SMN1's one copy or an expanded repeat"
        in text
    )
    assert "still answers `carrier`" not in text


def test_skill_covers_polygenic_breadth() -> None:
    """Every score is found through `--text pgsc_calc` (slice 3g), a weak one's weakness is
    said before its percentile, and a condition no score covers is answered by its stored
    reason, never with a percentile."""
    text = SKILL.read_text()
    section = " ".join(_section_named(text, "## Polygenic scores and APOE").split())
    for needed in (
        "`--text pgsc_calc`",
        "`src/genokeeper/resources/pgs/scores.toml`",
        '"weak predictive power"',
        "before the percentile",
        "`pgsc_calc:not-covered:<slug>`",
        "`--list-unassessed --text <condition> --min-validity limited`",
        '"not covered: no score meets the criteria"',
        "never with a percentile",
    ):
        assert needed in section, needed
    reasons = _bullet_text(text, "unassessed.by_gene.reasons")
    for quoted in ("`weak predictive power: …`", "`not covered: no score meets the"):
        assert quoted in reasons, quoted
    assert "`pgs_not_covered`" in text


def test_skill_covers_stage3e() -> None:
    """Slice 3e's repeats: RFC1 and CFTR's poly-T/TG, the exclusion-only normals as flagged
    non-findings, an abnormal there a "possible" finding at validity limited that the skill
    lowers the filters for, X-linked loci by the run's X dosage, CFTR's 5T never a carrier,
    TCF4 an eye examination, and the loci short reads cannot assess "not covered"."""
    section = _section_named(
        SKILL.read_text(), "## Stage 3e screening: repeats beyond the eleven loci"
    )
    for needed in (
        "`expansionhunter:RFC1`",
        "`expansionhunter:CFTR`",
        "flagged `exclusion`",
        '"possible"',
        "validity limited",
        "X dosage",
        'never "a cystic fibrosis carrier"',
        "eye examination",
        "not covered: <reason>",
        "`genotype --gene`",
    ):
        assert needed in section, needed


def test_skill_covers_deletions() -> None:
    """Deletions of 50 bp or more in the listed genes (slice 3b) are covered: each a
    finding at validity limited, "likely present" or "possible", always needing
    confirmation, with its chance of being false and its warnings; hidden by the default
    view, so the skill lowers the filters itself; "not screened" told from "none found" by
    the run's versions, never "no deletion"."""
    text = SKILL.read_text()
    section = _section_named(text, "## Stage 3b screening: deletions")
    for needed in (
        "`deletions:<gene>:<chrom>:<pos>-<end>`",
        "`deletions:screen`",
        "validity `limited`",
        '"likely present"',
        '"possible"',
        "needs confirmation by MLPA, a microarray, qPCR or long reads",
        "`--min-validity limited --min-reliability low`",
        "measured on HG002",
        "`found by one caller (DELLY)`",
        "`read depth not measured`",
        "`read depth does not show it, which a deletion this short often cannot`",
        "`no calibration call of this class: its chance of being false is extrapolated`",
        "`over a target outside the benchmark the chance of being false was measured in`",
        "`touches only the bases beside an exon, not the exon itself`",
        "`claim.plus_clinical_exons`",
        "`an end in a tandem repeat: may be a repeat-length change`",
        "`an end in a segmental duplication or low-mappability region`",
        "the screen found none; a miss is possible",
        "deletions were not screened in this run",
        'Never "no deletion"',
        "`needs_confirmation`",
        "`delly`",
        "`dysgu`",
    ):
        assert needed in section, needed
    flat = " ".join(text.split())
    row = (
        '| "Do I have a deletion in BRCA1?" |'
        " `--gene BRCA1 --min-validity limited --min-reliability low`"
    )
    assert row in flat
    coverage = _section_named(text, "## Coverage")
    assert "deletions of 50 bp or more" in coverage
    assert "nor are SLC37A4's coverage" in coverage


def test_skill_covers_slice3d() -> None:
    """Slice 3d's screens are covered: each non-normal result needs confirmation, a large
    event is shown, not calibrated, a "possible carrier" answer names the confirmation
    test from its reason, a screen that found none is never "no copy-number change", and
    the question table routes the information loci."""
    text = SKILL.read_text()
    section = _section_named(
        text, "## Stage 3d screening: copy number, insertions, paralogs"
    )
    for needed in (
        "`cnv:<gene>:<chrom>:<pos>-<end>`",
        "`insertions:<gene>:<chrom>:<pos>-<end>`",
        "`cnv:large:<chrom>:<pos>-<end>`",
        "`cnv:named:GJB6-D13S1830`",
        "`paralogs:<gene>`",
        "`lpa:LPA`",
        "Scramble",
        "GBA1",
        "Tier 4",
        '"shown, not calibrated: needs a chromosomal microarray"',
        '"not measured"',
        '"screening result, needs confirmation by <test>"',
        "`--min-validity limited --min-reliability low`",
        '"not covered: <reason>"',
        "`delly cnv`",
        '"the screen found none; a miss is possible"',
        'never "no copy-number change"',
        "`possible_carrier`",
        "the test its `validity.reason` names",
        "`gene.copy_number`",
        "`gene.insertions`",
        "`gene.possible_carriers`",
    ):
        assert needed in section, needed
    flat = " ".join(text.split())
    for row in (
        '| "How many AMY1 copies do I have?" | `--gene AMY1 --min-validity limited',
        '| "What is my RhD?" | `--gene RHD --min-validity limited',
    ):
        assert row in flat, row
    coverage = _section_named(text, "## Coverage")
    assert "Duplications are not covered" not in coverage
    assert "inversions and translocations" in coverage


def test_changes_contract_matches_the_output(tmp_path: Path) -> None:
    root, lab = tmp_path / "snapshots", tmp_path / "lab"
    base, run = sm.run_id(1), sm.run_id(2)
    for rid, validity in ((base, "limited"), (run, "strong")):
        sm.write(
            root,
            rid,
            sm.started(rid),
            {"clinvar": "20261011"},
            {"clinvar38": sm.file(rid)},
            [
                sm.statement("clinvar:1", validity=validity, rsid="rs1"),
                sm.statement(
                    "gwas:rs2-T:E",
                    vid="GRCh38:chr3:5:C:T",
                    gene="PPARG",
                    source="gwas_catalog",
                    topic="research",
                    reason=rid,
                ),
            ],
            [sm.call(), sm.call("GRCh38:chr3:5:C:T")],
        )
    sm.refreshed(tmp_path, lab, base, {"clinvar:1": {"litvar": ["1"], "clinvar": []}})
    started = datetime(2026, 10, 2, 11, tzinfo=UTC)
    run_list = publications_path(tmp_path, started)
    run_list.parent.mkdir(parents=True)
    run_list.write_text(
        json.dumps({"clinvar:1": {"litvar": ["1", "2"], "clinvar": []}})
    )
    append(
        lab,
        RefreshRecord(
            started=started,
            ended=started,
            outcome="partial",
            sources=(
                SourceOutcome("cpic", "20261001", None, SourceStatus.FAILED, "OSError"),
            ),
            base_run=base,
            pre_run=None,
            run=run,
        ),
    )
    (lab / "refresh" / "check.json").write_text(
        json.dumps(
            {
                "ran_at": started.isoformat(),
                "notified": {},
                "sources": [
                    {
                        "name": "clinvar",
                        "kind": "refreshed",
                        "upstream": "20261025",
                        "local": "20261011",
                        "status": "newer",
                        "error": None,
                    },
                    {
                        "name": "gwas",
                        "kind": "refreshed",
                        "upstream": None,
                        "local": None,
                        "status": "failed",
                        "error": "URLError",
                    },
                ],
            }
        )
    )
    document = report(
        tmp_path,
        lab,
        base=None,
        run=None,
        f=Filters(),
        every=True,
        gene=None,
        marker=None,
        text=None,
        now=started,
    )
    assert document["changes"] and document["research"] and document["failed"]
    assert document["check"]["newer"] and document["check"]["failed"] == ["gwas"]
    assert key_paths(document) == CHANGES_CONTRACT


def test_skill_answers_what_changed() -> None:
    """The question "what changed?" is answered from `changes --json`, trusted tier
    first; failed sources and their releases come first; a finding's three dates are
    named; ClinVar's papers are "attached by a lab"; the text is data."""
    text = SKILL.read_text()
    headings = [line for line in text.splitlines() if line.startswith("## ")]
    at = headings.index("## What changed")
    assert headings[at - 1 : at + 2] == [
        "## Answer shape",
        "## What changed",
        "## Coverage",
    ]
    section = _section_named(text, "## What changed")
    for needed in (
        "`uv run genokeeper changes --json`",
        "`--all`",
        "`failed`",
        "`before`",
        "`first_found`",
        "`last_changed`",
        "`assessment_changed`",
        "`on_or_before`",
        "attached by a lab",
        "data",
        "no earlier snapshot",
        "`check.stale`",
        "`check.days`",
        "`check.failed`",
        "has not run for",
    ):
        assert needed in section, needed
    privacy = " ".join(text.split())
    assert "`genokeeper changes`" in privacy[: privacy.index("## Commands")]
    # The PMIDs name the user's variants as a set: never looked up unasked, nor in bulk.
    rule = privacy[privacy.index("**Privacy.**") : privacy.index("## Commands")]
    assert "PMID" in rule and "one at a time" in rule
    # "On or before" can be applied from the JSON alone: the earliest of the three dates.
    assert "earliest of the three dates" in " ".join(section.split())
    # Papers are in the second tier, as the first tier is entered and left only.
    flat = " ".join(section.split())
    second = flat[flat.index("3. **Everything else**") : flat.index("4. **Papers**")]
    assert "`publications`" in second


def _stage3_claims(tmp_path: Path) -> dict[str, dict[str, Any]]:
    """What the stage-3 sources write for a claim, collected from synthetic tool output."""
    ctx = run_context(tmp_path)
    loci = load_grading().repeats.loci
    htt = loci["HTT"]
    allele = Allele(42, (41, 43), 6, 3, 0)
    calls = {
        variant: RepeatCall(variant, alleles, 30.0, False, 150)
        for variant, alleles in (
            (htt.variant_id, (allele,)),
            ("CFTR_TG", (Allele(11, (11, 11), 6, 0, 0), Allele(12, (12, 12), 6, 0, 0))),
            ("CFTR_T", (Allele(5, (5, 5), 6, 0, 0), Allele(7, (7, 7), 6, 0, 0))),
        )
    }
    rfc1 = Rfc1Call("carrier", "2 expanded pathogenic alleles", (allele,), 30.0, False)
    repeats = RepeatSource(
        calls=calls,
        not_run=None,
        loci=loci,
        withheld={},
        caps={},
        images={},
        pcr_like=False,
        insert=400.0,
        male=True,
        min_allele_reads=5,
        min_spanning=3,
        rfc1=rfc1,
        rfc1_not_run=None,
    )
    smn = SmnSource(
        SmnCaller(1, 2, 0, "high"), ParascopyCall(1, 2, 40.0, True), None, None, 20.0
    )
    yleaf = YleafSource(
        YHaplogroup("J-C134621", "C134621", 1.0, 170393), None, None, 1.0, True
    )
    rows = [
        *repeats.collect(ctx).rows,
        *smn.collect(ctx).rows,
        *yleaf.collect(ctx).rows,
    ]
    return {
        str(r.statement["statement_id"]): json.loads(str(r.statement["claim"]))
        for r in rows
    }


HG002 = Path(__file__).parent / "fixtures" / "hg002"


def _slice3f_claims(tmp_path: Path) -> dict[str, dict[str, Any]]:
    """One claim per slice 3f claim type, a result's where the type has one, collected by
    its source from HG002's HLA calls and synthetic tool output."""
    ctx = run_context(tmp_path)
    tags = {rsid: "0" for rsid in TAG_RSIDS} | {"rs2395182": "2", "rs7454108": "1"}
    long = Run("chr11", 48_025_953, 64_604_776, 16_578_824, 72.2)
    called = Called("0/1", "1", RowType.PASS, 40, 30, None, None, ())
    graded = (
        called,
        Reliability(p_wrong=1e-6, level=Level.HIGH, flags=(), warnings=()),
    )
    sources = [
        HlaSource(
            read_optitype(HG002 / "hg002.optitype_result.tsv.txt"),
            read_t1k(HG002 / "hg002.t1k_genotype.tsv.txt", genes=None),
            tags,
            None,
            {},
        ),
        MitomapSource(
            calls=[MitoCall(pos=3243, ref="A", alt="G", level=0.30, dp=3000)],
            coverage={pos: 3000 for pos in range(1, 16570)},
            clinvar_vids=frozenset(),
            floor=0.10,
            min_homoplasmy=0.9,
            min_dp=10,
        ),
        Hhv6Source(Hhv6Call("HHV-6B", 1.0, Kind.FINDING), None, None),
        RohSource(RohSummary(long.length, 1, long, (long,)), None, None),
        BloodGroupSource(
            calls={"FY": SystemCall("FY*02/FY*02", "Fy(a-b+)")},
            not_run=None,
            withheld={},
            uncovered={},
        ),
        AncestrySource(
            PopSimilarity("EUR", False, {"EUR": 0.9, "MID": 0.1}),
            {"EUR": 0.9, "MID": 0.1},
            ArchaicShare(57_768_000, 1_704_000, 0.011324),
            not_run={},
            withheld={},
        ),
        TraitSource(graded={allele: graded for allele in TRAIT_ALLELES}),
    ]
    claims: dict[str, dict[str, Any]] = {}
    for source in sources:
        for row in source.collect(ctx).rows:
            claim = json.loads(str(row.statement["claim"]))
            if claim["type"] not in claims or row.statement["kind"] == Kind.FINDING:
                claims[claim["type"]] = claim
    return claims


SLICE3F_TYPES = (
    *("hla_allele", "hla_association", "hla_coeliac", "hla_table", "hla_kir"),
    *("mitomap", "hhv6_ici", "roh_summary", "blood_group", "ancestry_similarity"),
    *("ancestry_admixture", "ancestry_neanderthal", "trait"),
)


def _claim_doc(text: str, claim_type: str) -> set[str]:
    """The claim keys the contract documents for `claim_type`: the keys its header names
    after "with its", and the bullets under the header, plus the `schema` and `type`."""
    _, lines = _section(text)
    at = next(i for i, line in enumerate(lines) if f"(`{claim_type}`," in line)
    end = next(i for i in range(at, len(lines)) if lines[i].endswith(":"))
    header = " ".join(" ".join(lines[at : end + 1]).split())
    keys = {"schema", "type"} | set(re.findall(r"`(\w+)`", header.split("with its")[1]))
    for line in lines[end + 1 :]:
        if not line.strip() or line.startswith("- "):
            break
        bullet = re.match(r"^  - `(\w+)`:", line)
        if bullet:
            keys.add(bullet[1])
    return keys


@pytest.mark.parametrize(
    ("statement", "claim_type"),
    [
        ("expansionhunter:HTT", "repeat_expansion"),
        ("expansionhunter:RFC1", "repeat_expansion"),
        ("expansionhunter:CFTR", "repeat_expansion"),
        ("smn:SMN1", "smn_copy_number"),
        ("yleaf:haplogroup", "y_haplogroup"),
    ],
)
def test_stage3_claims_carry_the_documented_keys(
    tmp_path: Path, statement: str, claim_type: str
) -> None:
    """What a stage-3 source collects has exactly the keys the skill documents for its
    claim type, so renaming one in the source (`g27134tg`, `smn1`) fails here."""
    claim = _stage3_claims(tmp_path)[statement]
    assert claim["type"] == claim_type
    assert set(claim) == _claim_doc(SKILL.read_text(), claim_type)


def test_slice3f_claim_types_are_the_sources(tmp_path: Path) -> None:
    assert sorted(_slice3f_claims(tmp_path)) == sorted(SLICE3F_TYPES)


@pytest.mark.parametrize("claim_type", SLICE3F_TYPES)
def test_slice3f_claims_carry_the_documented_keys(
    tmp_path: Path, claim_type: str
) -> None:
    """What a slice 3f source collects has exactly the keys the skill documents for its
    claim type."""
    claim = _slice3f_claims(tmp_path)[claim_type]
    assert set(claim) == _claim_doc(SKILL.read_text(), claim_type)


def test_stage3_statements_are_found_by_their_condition(tmp_path: Path) -> None:
    """`--text <disease>` finds a repeat locus, FMR1, RFC1, CFTR, a not-covered locus and
    SMN1 by the condition they name."""
    claims = _stage3_claims(tmp_path)
    not_covered = {entry.gene for entry in NOT_COVERED}
    genes = {*load_grading().repeats.loci, "RFC1", "CFTR", *not_covered}
    assert {c["gene"] for i, c in claims.items() if i.startswith("expansion")} == genes
    con = connect(tmp_path / "genomelab.duckdb")
    seed(
        con,
        [
            Row(
                i,
                kind="unassessed" if c["gene"] in not_covered else "finding",
                gene=c["gene"],
                level=None,
                topic="carrier",
                claim=c,
                # A not-covered locus has no gene call, as in a real run.
                gene_call=None
                if c["gene"] in not_covered
                else {"callers": {}, "level": "medium", "warnings": []},
                reason="not covered" if c["gene"] in not_covered else None,
            )
            for i, c in claims.items()
            if "gene" in c
        ],
    )
    loosest = Filters(min_reliability=Level.LOW, min_validity=Validity.CONFLICTING)
    for text, gene in (
        ("huntington", "HTT"),
        ("spinal muscular atrophy", "SMN1"),
        ("fragile", "FMR1"),
        ("canvas", "RFC1"),
        ("cftr-related", "CFTR"),
        ("neuronal intranuclear", "NOTCH2NLC"),
    ):
        found = query(con, None, replace(loosest, text=text), list_unassessed=True)
        listed = found.to_json()
        assert gene in {
            f["gene"] for f in listed["findings"] + listed["unassessed_sites"]
        }


def _bullet_text(text: str, path: str) -> str:
    """The documented bullet for `path` with its continuation lines, whitespace joined."""
    lines = text.splitlines()
    i = bullets(text)[path]
    end = i + 1
    while end < len(lines) and lines[end].strip() and not _BULLET.match(lines[end]):
        end += 1
    return " ".join(" ".join(lines[i:end]).split())


def test_skill_lists_exactly_the_sources_the_code_emits() -> None:
    emitted: set[str] = set()
    package = Path(cli.__file__).parent
    # Slice 3d's insertions are 3b's source in their own words, named in `insertions.py`.
    for path in [*(package / "sources").glob("*.py"), package / "insertions.py"]:
        for node in ast.walk(ast.parse(path.read_text())):
            if (
                isinstance(node, ast.Assign)
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id.endswith("_TOOL")
                and isinstance(node.value, ast.Constant)
            ):
                emitted.add(str(node.value.value))
    listed = _bullet_text(SKILL.read_text(), "sources.source")
    assert set(re.findall(r"`(\w+)`", listed.split(":", 1)[1])) == emitted


def test_stored_reasons_the_skill_quotes_exist_in_the_code() -> None:
    """A reason or warning in backticks (more than one word) is a string some module
    writes; text around an ellipsis is checked piece by piece."""
    text = SKILL.read_text()
    code = " ".join(p.read_text() for p in Path(cli.__file__).parent.rglob("*.py"))
    missing: list[str] = []
    for path in ("unassessed.by_gene.reasons", "reliability.warnings"):
        quoted = re.findall(r"`([^`]+)`", _bullet_text(text, path).split(":", 1)[1])
        for span in (q for q in quoted if " " in q):
            for piece in re.split(r"…", span):
                piece = piece.strip(" ,")
                if piece and piece not in code:
                    missing.append(piece)
    assert missing == []


def test_commands_block_names_real_subcommands_and_flags() -> None:
    block = SKILL.read_text().split("## Commands")[1].split("```bash")[1]
    block = block.split("```")[0].replace("\\\n", " ")
    parser = cli.build_parser()
    action = next(
        a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
    )  # pyright: ignore[reportPrivateUsage]
    choices: dict[str, argparse.ArgumentParser] = action.choices
    for command in block.strip().splitlines():
        words = command.split()
        assert words[:3] == ["uv", "run", "genokeeper"], command
        sub = choices[words[3]]
        flags = {f for a in sub._actions for f in a.option_strings}  # pyright: ignore[reportPrivateUsage]
        assert set(re.findall(r"--[\w-]+", command)) <= flags, command


def test_privacy_rule_names_every_channel() -> None:
    """Privacy audit C1, I5 and I6: the skill's privacy rule names every place a genotype could
    persist or leave through, and states the deny rules' bound."""
    text = " ".join(SKILL.read_text().split())
    rule = text[text.index("**Privacy.**") : text.index("## Commands")]
    for channel in (
        "auto-memory",
        "commit messages",
        "PR bodies",
        "issues",
        "Artifacts",
        "side-session notes",
        "cloud",
        "routine",
        "Telegram",
        "web search",
        "`grep -r`",
    ):
        assert channel in rule, channel


# --- slice 3c ---------------------------------------------------------------------------------


def _dosage_claim(tmp_path: Path) -> dict[str, Any]:
    """The claim `DosageSource` writes for a measured, calibrated dosage."""
    t = replace(load_grading().dosage, calibrated_on="2026-11-01 HG002, HG003, HG004")
    source = DosageSource(Dosage(0.5, 0.5, 0.05, 2400), None, t)
    (row,) = source.collect(run_context(tmp_path)).rows
    return json.loads(str(row.statement["claim"]))


def _dosage_row(tmp_path: Path) -> Row:
    """The dosage statement as an unassessed claim lists it (its claim's keys)."""
    return Row(
        DOSAGE_STATEMENT,
        kind="unassessed",
        gene=None,
        level=None,
        reason=DOSAGE_NOT_CALIBRATED,
        claim=_dosage_claim(tmp_path),
    )


def test_dosage_claim_carries_the_documented_keys(tmp_path: Path) -> None:
    claim = _dosage_claim(tmp_path)
    assert claim["type"] == "sex_chromosome_dosage" and claim["label"] == "46,XY"
    assert set(claim) == _claim_doc(SKILL.read_text(), "sex_chromosome_dosage")


def test_skill_covers_stage3c() -> None:
    """Identity as a FASTQ mix-up check read from `versions.identity`, the coverage note
    worded as `query` writes it, the dosage statement read with `show`, mosaic loss of Y
    as research, both allele-fraction warnings, and the reference as a change's cause."""
    text = SKILL.read_text()
    section = _section_named(
        text, "## Stage 3c: identity, coverage, sex chromosomes, allele fraction"
    )
    flat = " ".join(section.split())
    for needed in (
        "`versions.identity`",
        "FASTQ mix-up",
        "not a lab",
        "`not reliably assessed: 82% callable`",
        "`assessed by <owner>`",
        "`not measured`",
        f"`uv run genokeeper show {DOSAGE_STATEMENT} --json`",
        "screening result, needs a karyotype",
        "possible mosaic loss of Y",
        "research",
        f"`{AF_LOW}`",
        f"`{AF_CLONAL}`",
        "`reference`",
        "not new knowledge",
    ):
        assert needed in flat, needed
    causes = _bullet_text(text, "changes.causes")
    assert "`reference`" in causes
    headings = [line for line in text.splitlines() if line.startswith("## ")]
    at = headings.index("## Stage 3b screening: deletions")
    assert headings[at + 1].startswith("## Stage 3c")


def test_coverage_note_is_the_one_query_writes(tmp_path: Path) -> None:
    """The skill quotes the note `query --gene` writes, word for word."""
    con = connect(tmp_path / "genomelab.duckdb")
    seed(con, [Row("clinvar:1", gene="PMS2")])
    seed_coverage(con, "PMS2", 0.82)
    note = query(con, None, Filters(gene="PMS2")).to_json()["coverage"]["note"]
    assert note == "not reliably assessed: 82% callable"
    assert f"`{note}`" in " ".join(SKILL.read_text().split())
