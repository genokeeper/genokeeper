"""Ancestry similarity, admixture and the archaic share (slice 3f, "Admixture"; rulings R5,
R6).

`fixtures/hg002/hg002.fastmixture.Q.txt` is fastmixture 1.3.0's projection of HG002's GIAB
v5.0q truth genotypes on the panel's supervised P file (K = 5, labels 1..5 = AFR, CSA, EAS,
EUR, MID). `hg002.hmmix.decoded.txt` is hmmix 0.9.2's `decode -admixpop` on the same truth
(Task 5)."""

import json
from pathlib import Path

import pytest
from fixtures.run_context import run_context

from genokeeper.ancestry import admixture, archaic
from genokeeper.ancestry.archaic import (
    ARCHAIC_GENOME_BP,
    ArchaicShare,
    Segment,
    archaic_share,
    read_segments,
)
from genokeeper.prs.pgsc import PopSimilarity
from genokeeper.sources import INFORMATION_NOTE
from genokeeper.sources.ancestry import (
    ADMIXTURE,
    ANCESTRY_TOOL,
    ARCHAIC,
    SIMILARITY,
    AncestrySource,
)

HG002 = Path(__file__).parent / "fixtures" / "hg002"
GROUPS = ("AFR", "CSA", "EAS", "EUR", "MID")
# Coll Macià et al. 2026 (Mol Biol Evol 43:msag134): NA20810, hmmix on hg38, strict mask.
PUBLISHED_SHARE = 0.0118


def _groups_json(folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    groups = {pop: i for i, pop in enumerate(GROUPS, start=1)}
    path = folder / "groups.json"
    path.write_text(json.dumps({"K": 5, "groups": groups, "unlabelled": ["AMR"]}))
    return path


def test_read_q_hg002(tmp_path: Path) -> None:
    """One row of K proportions, summing to 1, labelled in `groups.json`'s label order."""
    groups = admixture.read_groups(_groups_json(tmp_path))
    assert groups == GROUPS
    q = admixture.read_q(HG002 / "hg002.fastmixture.Q.txt", groups)
    assert tuple(q) == GROUPS
    assert abs(sum(q.values()) - 1) < 1e-6
    assert max(q, key=q.__getitem__) == "EUR"
    assert q["EUR"] == pytest.approx(0.850888)
    with pytest.raises(ValueError, match="5 columns"):
        admixture.read_q(HG002 / "hg002.fastmixture.Q.txt", GROUPS[:4])


def test_pruned_alleles_from_bim(tmp_path: Path) -> None:
    """The pruned SNPs as `write_target_vcf` takes them: `chr<N>`, REF the bim's A2."""
    bim = tmp_path / "pruned.bim"
    bim.write_text("1\t1:13273:G:C\t0\t13273\tC\tG\n22\t22:5:A:T\t0\t5\tT\tA\n")
    assert admixture.pruned_alleles(bim) == [
        ("chr1", 13273, "G", "C"),
        ("chr22", 5, "A", "T"),
    ]


def test_projection_commands(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """One plink2 `--make-bed` in its image, then fastmixture's projection on the panel's
    P file, both offline and gated `dna:ancestry`; the Q file is read against the labels
    of the P file's `groups.json`."""
    seen: list[tuple[str, list[str], dict[str, object]]] = []
    ref = tmp_path / "ref"
    _groups_json(ref)
    p_file = ref / "reference.K5.s1.P"
    p_file.write_text("")

    def fake(name: str, args: list[str], **kwargs: object) -> None:
        seen.append((name, args, kwargs))
        if name == "fastmixture":
            q = tmp_path / "proj" / "projection.K5.s1.Q"
            q.write_text((HG002 / "hg002.fastmixture.Q.txt").read_text())
        else:
            for ext in ("bed", "bim", "fam"):
                (tmp_path / "bed" / f"owner.{ext}").write_text("")

    monkeypatch.setattr(admixture, "run_container", fake)
    vcf = tmp_path / "in" / "pruned.vcf.gz"
    vcf.parent.mkdir()
    vcf.write_text("")
    bed = admixture.to_bed(vcf, tmp_path / "bed" / "owner.bed")
    assert bed == tmp_path / "bed" / "owner.bed"
    q = admixture.run_projection(bed, p_file, 5, tmp_path / "proj")
    assert q["EUR"] == pytest.approx(0.850888)
    (plink, plink_args, plink_kw), (fm, fm_args, fm_kw) = seen
    assert (plink, fm) == ("pgsc_plink2", "fastmixture")
    assert plink_args == [
        "plink2", "--vcf", "/in/pruned.vcf.gz", "--set-all-var-ids", "@:#:$r:$a",
        "--make-bed", "--out", "/out/owner",
    ]  # fmt: skip
    assert fm_args == [
        "fastmixture", "--bfile", "/in/owner", "--K", "5", "--seed", "1",
        "--threads", "1", "--projection", "/ref/reference.K5.s1.P",
        "--out", "/out/projection",
    ]  # fmt: skip
    assert plink_kw["label"] == fm_kw["label"] == "dna:ancestry"
    with pytest.raises(ValueError, match="K 4"):
        admixture.run_projection(bed, p_file, 4, tmp_path / "proj")


def test_archaic_share_arithmetic() -> None:
    """R15: the share is the Neanderthal segments' span alone, so 10 Mb Neanderthal and 1 Mb
    Denisovan over a 1,000 Mb denominator is 0.010 exactly, the Denisovan span kept as a
    field; a Human segment and an unclassified one count for neither."""

    def seg(start: int, length: int, state: str, n: int, d: int) -> Segment:
        return Segment("chr1", start, start + length, state, neanderthal=n, denisova=d)

    segments = [
        seg(0, 10_000_000, "Archaic", 25, 4),
        seg(20_000_000, 1_000_000, "Archaic", 2, 9),
        seg(30_000_000, 5_000_000, "Archaic", 0, 0),
        seg(40_000_000, 50_000_000, "Human", 9, 0),
    ]
    assert [s.origin for s in segments] == [
        "neanderthal",
        "denisovan",
        None,
        None,
    ]
    share = archaic_share(segments, genome_bp=1_000_000_000)
    assert share == ArchaicShare(
        neanderthal_bp=10_000_000, denisovan_bp=1_000_000, share=0.010
    )
    assert ARCHAIC_GENOME_BP == 2 * 2_550_757_832


def test_archaic_share_hg002() -> None:
    """HG002's Neanderthal share over `ARCHAIC_GENOME_BP` is 0.011324, within ±30% of the
    published 0.0118, so a denominator doubled or halved (or the strict mask's sum)
    fails."""
    segments = read_segments(HG002 / "hg002.hmmix.decoded.txt")
    assert sum(s.state == "Archaic" for s in segments) == 1695
    share = archaic_share(segments)
    assert (share.neanderthal_bp, share.denisovan_bp) == (57_768_000, 1_704_000)
    assert round(share.share, 6) == 0.011324
    assert abs(share.share - PUBLISHED_SHARE) <= 0.3 * PUBLISHED_SHARE
    for wrong in (ARCHAIC_GENOME_BP * 2, ARCHAIC_GENOME_BP // 2):
        off = archaic_share(segments, genome_bp=wrong).share
        assert abs(off - PUBLISHED_SHARE) > 0.3 * PUBLISHED_SHARE


def test_run_hmmix_one_container(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """create_ingroup, train and decode with the archaic BCFs in one hmmix container,
    offline and gated `dna:ancestry`, on the strict mask and mutation rates; the decoded
    file returned."""
    seen: dict[str, object] = {}

    def fake(name: str, args: list[str], **kwargs: object) -> None:
        seen.update(name=name, args=args, **kwargs)
        (tmp_path / "out" / "decoded.diploid.txt").write_text("chrom\n")

    monkeypatch.setattr(archaic, "run_container", fake)
    vcf = tmp_path / "calls" / "owner.vcf.gz"
    vcf.parent.mkdir()
    vcf.write_text("")
    decoded = archaic.run_hmmix(vcf, tmp_path / "hmmix", tmp_path / "out")
    assert decoded == tmp_path / "out" / "decoded.diploid.txt"
    assert (seen["name"], seen["label"]) == ("hmmix", "dna:ancestry")
    script = " ".join(str(a) for a in seen["args"])  # type: ignore[union-attr]
    for step in (
        "hmmix create_ingroup",
        "hmmix train",
        "hmmix decode",
        "-admixpop=",
        "-weights=/h/hg38_strick_callability_mask.bed",
        "-mutrates=/h/hg38_mutationrate.bed",
        "-outgroup=/h/hg38_Outgroup_1000g_HGDP.txt",
        "/in/owner.vcf.gz",
    ):
        assert step in script, step


def _statements(tmp_path: Path, source: AncestrySource) -> dict[str, dict[str, object]]:
    collected = source.collect(run_context(tmp_path))
    assert collected.tool == ANCESTRY_TOOL == "ancestry"
    assert collected.input_rows == len(collected.rows) == 3
    out: dict[str, dict[str, object]] = {}
    for row in collected.rows:
        s = row.statement
        assert (s["topic"], s["validity"], s["actionable"], s["exclusion"]) == (
            "ancestry",
            "limited",
            False,
            False,
        )
        out[str(s["statement_id"])] = s
    assert set(out) == {
        "ancestry:similarity",
        "ancestry:admixture",
        "ancestry:neanderthal",
    }
    return out


SIM = PopSimilarity(
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
ADMIX = dict(
    zip(GROUPS, (0.026740, 0.096949, 0.002714, 0.850888, 0.022709), strict=True)
)
SHARE = ArchaicShare(neanderthal_bp=57_768_000, denisovan_bp=1_704_000, share=0.011324)


def _full(**kwargs: object) -> AncestrySource:
    base: dict[str, object] = {
        "similarity": SIM,
        "admixture": ADMIX,
        "archaic": SHARE,
        "not_run": {},
        "withheld": {},
    }
    return AncestrySource(**{**base, **kwargs})  # type: ignore[arg-type]


def test_similarity_never_called_proportion(tmp_path: Path) -> None:
    """`ancestry:similarity` is pgsc_calc's random-forest similarity: its claim names
    probabilities, and neither its claim nor its reason calls it a proportion (the reason
    only says it is not one)."""
    s = _statements(tmp_path, _full())["ancestry:similarity"]
    claim = str(s["claim"])
    reason = str(s["reason"])
    assert "proportion" not in claim.lower()
    assert SIMILARITY in reason
    assert "proportion" not in reason.replace("not ancestry proportions", "")
    assert json.loads(claim)["probabilities"]["EUR"] == 0.65


def test_ancestry_findings_are_information(tmp_path: Path) -> None:
    """All three are findings of information under ancestry at validity limited, each
    with its wording and `INFORMATION_NOTE`, its source level the tool and its version."""
    statements = _statements(tmp_path, _full())
    for sid, wording, tool in (
        ("ancestry:similarity", SIMILARITY, "pgsc_calc 2.3.0"),
        ("ancestry:admixture", ADMIXTURE, "fastmixture 1.3.0"),
        ("ancestry:neanderthal", ARCHAIC, "hmmix 0.9.2"),
    ):
        s = statements[sid]
        assert s["kind"] == "finding", sid
        assert wording in str(s["reason"]) and INFORMATION_NOTE in str(s["reason"])
        assert s["source_level"] == tool
    assert ARCHAIC == "every non-African genome carries some; about 1–2% by this method"
    admix = json.loads(str(statements["ancestry:admixture"]["claim"]))
    assert admix["proportions"]["EUR"] == pytest.approx(0.850888)
    archaic_claim = json.loads(str(statements["ancestry:neanderthal"]["claim"]))
    assert (archaic_claim["neanderthal_bp"], archaic_claim["share"]) == (
        57_768_000,
        0.011324,
    )


def test_not_run_unassessed(tmp_path: Path) -> None:
    """A step not run, a result missing, or a step withheld leaves its statement
    unassessed with the reason; the others stand."""
    statements = _statements(
        tmp_path,
        _full(
            admixture=None,
            not_run={"neanderthal": "hmmix failed"},
            withheld={"similarity": "not covered: x"},
        ),
    )
    assert (
        statements["ancestry:neanderthal"]["kind"],
        statements["ancestry:neanderthal"]["reason"],
    ) == ("unassessed", "hmmix failed")
    assert (
        statements["ancestry:similarity"]["kind"],
        statements["ancestry:similarity"]["reason"],
    ) == ("unassessed", "not covered: x")
    assert statements["ancestry:admixture"]["kind"] == "unassessed"
    none = _statements(
        tmp_path, AncestrySource(None, None, None, not_run={}, withheld={})
    )
    assert {s["kind"] for s in none.values()} == {"unassessed"}
