"""MITOMAP's Confirmed variants where ClinVar is silent, as a source (`MitomapSource`; slice
3f, "chrM").

`resources/mito/mitomap_confirmed.tsv` is MITOMAP's Confirmed list (CC BY 3.0), cut to the
single-base substitutions mutserve calls. Each listed variant the user carries at any level
(mutserve reported that allele at its position) and that no ClinVar P/LP claim holds (its
vid not in `clinvar_vids`) is one statement `mitomap:<pos><ref>><alt>` under risks, validity
moderate (limited for a row ClinGen classifies VUS, `[VUS*]`), on the call
`mito.mito_called` gives with the run's floor and why heteroplasmies are withheld (the
calibration's or haplocheck's reason): a known copy a finding, an unknown one unassessed.
The rest are
dropped as `not carried` or `ClinVar claim`. Nothing is actionable or flagged `exclusion`.
"""

import csv
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

from genokeeper.genotype.reliability import StrataIndex
from genokeeper.mito import CHRM, MitoCall, grade_mito, mito_called
from genokeeper.sources import Collected, Fields, Row, RunContext, present
from genokeeper.sources.panel import CLAIM_SCHEMA, call_row
from genokeeper.sources.repeats import CONFIRM
from genokeeper.store.db import Kind, Topic, Validity
from genokeeper.variants.key import make_vid

MITOMAP_TOOL = "mitomap"
MITOMAP_FILE = (
    Path(__file__).resolve().parents[1] / "resources" / "mito" / "mitomap_confirmed.tsv"
)
MITOMAP_VERSION = "2026-10-30"
SOURCE_LEVEL = "MITOMAP Confirmed"
MITOMAP_REASON = f"{CONFIRM} by mtDNA testing in a clinical laboratory"
NOT_CARRIED = "not carried"
CLINVAR_CLAIM = "ClinVar claim"
_SNV = re.compile(r"m\.(\d+)([ACGT])>([ACGT])")


@dataclass(frozen=True)
class MitomapVariant:
    """One Confirmed single-base substitution: its position, bases and MITOMAP's columns."""

    pos: int
    ref: str
    alt: str
    allele: str
    locus: str
    diseases: str
    status: str


@cache
def load_mitomap() -> tuple[MitomapVariant, ...]:
    """The packaged Confirmed list; a row that is no single-base substitution is refused."""
    with MITOMAP_FILE.open(newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    found: list[MitomapVariant] = []
    for row in csv.DictReader(lines, delimiter="\t"):
        m = _SNV.fullmatch(row["allele"])
        if m is None:
            raise ValueError(f"{MITOMAP_FILE}: not a substitution: {row['allele']}")
        found.append(
            MitomapVariant(
                pos=int(m[1]),
                ref=m[2],
                alt=m[3],
                allele=row["allele"],
                locus=row["locus"],
                diseases=row["diseases"],
                status=row["status"],
            )
        )
    return tuple(found)


@dataclass(frozen=True)
class MitomapSource:
    """mutserve's variants and coverage, the chrM vids ClinVar's P/LP claims hold, and what
    `mito_called` grades with (the run's floor, `min_homoplasmy`, `min_dp`, why every
    heteroplasmy is withheld: the calibration's or haplocheck's reason); `strata` flags the
    calls."""

    calls: Sequence[MitoCall]
    coverage: Mapping[int, int]
    clinvar_vids: frozenset[str]
    floor: float
    min_homoplasmy: float
    min_dp: int
    heteroplasmy_withheld: str | None = None
    strata: StrataIndex = field(default_factory=lambda: StrataIndex({}))
    tool: str = MITOMAP_TOOL

    def collect(self, ctx: RunContext) -> Collected:
        """One statement per carried Confirmed variant ClinVar does not claim."""
        carried = {(c.pos, c.alt) for c in self.calls}
        rows: list[Row] = []
        dropped = {NOT_CARRIED: 0, CLINVAR_CLAIM: 0}
        variants = load_mitomap()
        for v in variants:
            if (v.pos, v.alt) not in carried:
                dropped[NOT_CARRIED] += 1
                continue
            vid = make_vid(CHRM, v.pos, v.ref, v.alt)
            if vid in self.clinvar_vids:
                dropped[CLINVAR_CLAIM] += 1
                continue
            called = mito_called(
                (CHRM, v.pos, v.ref, v.alt),
                self.calls,
                self.coverage,
                min_dp=self.min_dp,
                min_homoplasmy=self.min_homoplasmy,
                floor=self.floor,
                heteroplasmy_withheld=self.heteroplasmy_withheld,
            )
            reliability = grade_mito(
                called,
                CHRM,
                v.pos,
                strata=self.strata,
                min_homoplasmy=self.min_homoplasmy,
            )
            kind = Kind.UNASSESSED if called.dosage == "unknown" else Kind.FINDING
            claim = {
                "schema": CLAIM_SCHEMA,
                "type": "mitomap",
                "allele": v.allele,
                "locus": v.locus,
                "condition": v.diseases,
                "status": v.status,
                "level": called.vaf,
            }
            statement: Fields = {
                "run_id": ctx.run_id,
                "statement_id": f"{MITOMAP_TOOL}:{v.pos}{v.ref}>{v.alt}",
                "kind": str(kind),
                "vid": vid,
                "gene": v.locus,
                "claim": json.dumps(claim, ensure_ascii=False),
                "validity": str(
                    Validity.LIMITED if "[VUS*]" in v.status else Validity.MODERATE
                ),
                "source_level": SOURCE_LEVEL,
                "reason": MITOMAP_REASON,
                "effect": None,
                "actionable": False,
                "topic": str(Topic.RISKS),
                "source": MITOMAP_TOOL,
                "source_id": v.allele,
                "source_version": MITOMAP_VERSION,
                "points": None,
                "exclusion": False,
            }
            variant: Fields = {"run_id": ctx.run_id, "vid": vid, "rsid": None}
            graded = call_row(ctx.run_id, vid, called, reliability)
            level = present({"level": called.vaf})
            rows.append(Row(variant, graded, statement, level))
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=len(variants),
            dropped={k: n for k, n in dropped.items() if n},
        )
