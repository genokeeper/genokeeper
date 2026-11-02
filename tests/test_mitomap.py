"""MITOMAP's Confirmed variants where ClinVar is silent (`sources.mitomap.MitomapSource`;
slice 3f, "chrM"): one statement per Confirmed variant the user carries at any level that no
ClinVar P/LP claim holds, graded on mutserve's call with the run's floor."""

import json
from pathlib import Path
from typing import Any

from fixtures.run_context import run_context

from genokeeper.mito import MATERNAL_NOTE, TISSUE_NOTE, MitoCall
from genokeeper.sources.mitomap import (
    MITOMAP_TOOL,
    MitomapSource,
    load_mitomap,
)
from genokeeper.variants.key import make_vid

COVERAGE = {pos: 3000 for pos in range(1, 16570)}
MELAS = MitoCall(pos=3243, ref="A", alt="G", level=0.30, dp=3000)
LHON = MitoCall(pos=11778, ref="G", alt="A", level=0.99, dp=3000)


def collect(tmp_path: Path, calls: list[MitoCall], **keywords: Any):
    source = MitomapSource(
        calls=calls,
        coverage=COVERAGE,
        clinvar_vids=frozenset({make_vid("chrM", 11778, "G", "A")}),
        floor=0.10,
        min_homoplasmy=0.9,
        min_dp=10,
        **keywords,
    )
    return source.collect(run_context(tmp_path))


def test_resource_is_the_snv_cut() -> None:
    rows = load_mitomap()
    assert len(rows) == 116
    assert all(len(r.ref) == len(r.alt) == 1 for r in rows)
    melas = next(r for r in rows if (r.pos, r.alt) == (3243, "G"))
    assert (melas.locus, melas.status) == ("MT-TL1", "Cfrm [P]")
    assert melas.diseases.startswith("MELAS")


def test_mitomap_only_where_clinvar_is_silent(tmp_path: Path) -> None:
    # MELAS at 30%, LHON (ClinVar's own claim) homoplasmic, and a common variant MITOMAP
    # does not list.
    other = MitoCall(pos=263, ref="A", alt="G", level=1.0, dp=3000)
    found = collect(tmp_path, [MELAS, LHON, other])
    assert found.tool == MITOMAP_TOOL == "mitomap"
    assert found.input_rows == 116
    assert found.dropped == {"not carried": 114, "ClinVar claim": 1}
    (row,) = found.rows
    s = row.statement
    assert (s["statement_id"], s["kind"], s["topic"], s["validity"]) == (
        "mitomap:3243A>G",
        "finding",
        "risks",
        "moderate",
    )
    assert (s["source"], s["source_level"], s["gene"]) == (
        "mitomap",
        "MITOMAP Confirmed",
        "MT-TL1",
    )
    assert s["vid"] == make_vid("chrM", 3243, "A", "G")
    assert (s["actionable"], s["exclusion"]) == (False, False)
    claim = json.loads(str(s["claim"]))
    assert (claim["type"], claim["allele"], claim["status"]) == (
        "mitomap",
        "m.3243A>G",
        "Cfrm [P]",
    )
    assert row.variant is not None and row.variant["vid"] == s["vid"]
    assert row.call is not None
    assert (row.call["dosage"], row.call["vaf"]) == ("1", 0.30)
    assert TISSUE_NOTE in row.call["warnings"]  # pyright: ignore[reportOperatorIssue]
    assert MATERNAL_NOTE in row.call["warnings"]  # pyright: ignore[reportOperatorIssue]
    # m.583G>A is Cfrm [VUS*] (ClinGen VUS): a finding at validity limited.
    (vus,) = collect(tmp_path, [MitoCall(583, "G", "A", 1.0, 3000)]).rows
    assert vus.statement["validity"] == "limited"


def test_mitomap_below_floor_unassessed(tmp_path: Path) -> None:
    """Carried at any level: below the floor it is a statement, unassessed with its level."""
    low = MitoCall(pos=3243, ref="A", alt="G", level=0.05, dp=3000)
    (row,) = collect(tmp_path, [low]).rows
    assert row.statement["kind"] == "unassessed"
    assert row.call is not None
    assert (row.call["dosage"], row.call["vaf"]) == ("unknown", 0.05)
    assert any("below the floor 0.1" in w for w in row.call["warnings"])  # pyright: ignore[reportGeneralTypeIssues, reportOptionalIterable]
    # Another base at a listed position is not the listed variant: not carried.
    other_base = MitoCall(pos=3243, ref="A", alt="C", level=0.95, dp=3000)
    assert collect(tmp_path, [other_base]).rows == []


def test_mitomap_passes_withheld_on(tmp_path: Path) -> None:
    reason = "not covered: false heteroplasmy at 0.25 in calibration"
    melas = MitoCall(pos=3243, ref="A", alt="G", level=0.40, dp=3000)
    (row,) = collect(tmp_path, [melas], heteroplasmy_withheld=reason).rows
    assert row.statement["kind"] == "unassessed"
    assert row.call is not None
    assert any(reason in w for w in row.call["warnings"])  # pyright: ignore[reportGeneralTypeIssues, reportOptionalIterable]
    homoplasmic = MitoCall(pos=3243, ref="A", alt="G", level=0.95, dp=3000)
    (row,) = collect(tmp_path, [homoplasmic], heteroplasmy_withheld=reason).rows
    assert row.statement["kind"] == "finding"
