"""Blood groups as a source (`BloodGroupSource`; slice 3f, "Blood groups"). Tool `rbceq2`,
one statement `bloodgroup:<system>` per system of the packaged positions
(`bloodgroups.systems`) under traits, validity limited, never actionable, never flagged
`exclusion`.

- A system RBCeq2 called is a finding of information: its genotype and phenotype in the
  claim, `NOT_CLINICAL` and `INFORMATION_NOTE` in the reasons (only our genotypes'
  agreement with the truth is measured, not RBCeq2's typing).
- `not_run` (no table, a failed step), `withheld[system]` (the calibration's reason) and `uncovered[system]` (R10: a no-call among its written positions, which RBCeq2
  would read as the reference) make a system unassessed, as does no call from RBCeq2.
- `unchecked[system]` (R14: alleles left out of RBCeq2's input because their GRCh38 REF
  disagrees with the reference, `bloodgroups.allele_name`s) are named in the system's
  reasons (`UNCHECKED`); the system stays shown.
- `bloodgroup:RH`, `bloodgroup:MNS` and `bloodgroup:HPA` are unassessed every run
  (`RH_REASON`, `MNS_REASON`, `HPA_REASON`); RBCeq2's own calls of the MNS genes (`GYPA`,
  `GYPB`, `GYP`) are never shown, and a system it calls with no packaged position (one
  defined only by a large deletion) is dropped.
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass, field

from genokeeper.bloodgroups import MNS_GENES, SOURCE_VERSION, SystemCall, systems
from genokeeper.sources import INFORMATION_NOTE, Collected, Fields, Row, RunContext
from genokeeper.sources.panel import CLAIM_SCHEMA
from genokeeper.store.db import Kind, Topic, Validity

BLOODGROUP_TOOL = "rbceq2"
NOT_CLINICAL = "not for clinical use: a transfusion is always cross-matched by serology"
RH_REASON = (
    "not covered: RhD and RHCE need copy number that short reads do not give"
    " (RBCeq2 2.4.4); RhD zygosity is stated under RHD"
)
MNS_REASON = (
    "not covered: MNS depends on GYPA/GYPB hybrid genes that the truth genomes do not"
    " carry, so it cannot be checked here"
)
HPA_REASON = (
    "not covered: RBCeq2 2.4.4's platelet-antigen positions are GRCh37 coordinates in"
    " its GRCh38 table"
)
UNCOVERED = (
    "genotype unknown at {positions}: RBCeq2 would read a no-call as the reference"
)
UNCHECKED = (
    "RBCeq2 allele {name} not checked: its GRCh38 REF disagrees with the reference"
)
NO_CALL = "RBCeq2 made no call"
NOT_RUN = "not run"
SOURCE_LEVEL = f"RBCeq2 {SOURCE_VERSION}"
# Drop reasons of RBCeq2's calls that no statement stands for.
MNS_GENE = "mns_gene"
NO_POSITION = "no_packaged_position"
_FIXED = {"RH": RH_REASON, "MNS": MNS_REASON, "HPA": HPA_REASON}


@dataclass(frozen=True)
class BloodGroupSource:
    """RBCeq2's calls (`bloodgroups.read_rbceq2`; None when it has none), why it was not
    run, the systems withheld with their reasons, the systems with no-calls among their
    written positions (`bloodgroups.uncovered`), and each system's alleles left out of the
    input (`Rbceq2Positions.unchecked`)."""

    calls: Mapping[str, SystemCall] | None
    not_run: str | None
    withheld: Mapping[str, str]
    uncovered: Mapping[str, tuple[int, ...]]
    unchecked: Mapping[str, tuple[str, ...]] = field(
        default_factory=dict[str, tuple[str, ...]]
    )
    tool: str = BLOODGROUP_TOOL

    def _kind(self, system: str) -> tuple[Kind, str, SystemCall | None]:
        if self.not_run is not None:
            return Kind.UNASSESSED, self.not_run, None
        if system in self.withheld:
            return Kind.UNASSESSED, self.withheld[system], None
        if system in self.uncovered:
            where = ", ".join(str(p) for p in self.uncovered[system])
            return Kind.UNASSESSED, UNCOVERED.format(positions=where), None
        if self.calls is None:
            return Kind.UNASSESSED, NOT_RUN, None
        call = self.calls.get(system)
        if call is None or call.genotype is None:
            return Kind.UNASSESSED, NO_CALL, None
        return Kind.FINDING, f"{NOT_CLINICAL}; {INFORMATION_NOTE}", call

    def _noted(self, system: str) -> tuple[Kind, str, SystemCall | None]:
        kind, reason, call = self._kind(system)
        notes = [UNCHECKED.format(name=n) for n in self.unchecked.get(system, ())]
        return kind, "; ".join([reason, *notes]), call

    def collect(self, ctx: RunContext) -> Collected:
        """One statement per packaged system, and RH, MNS and the platelet antigens."""
        shown = systems()
        rows = [self._row(ctx, s, *self._noted(s)) for s in shown]
        rows += [self._row(ctx, s, Kind.UNASSESSED, r, None) for s, r in _FIXED.items()]
        dropped = {MNS_GENE: 0, NO_POSITION: 0}
        for system in self.calls or {}:
            if system in MNS_GENES:
                dropped[MNS_GENE] += 1
            elif system not in shown:
                dropped[NO_POSITION] += 1
        dropped = {k: n for k, n in dropped.items() if n}
        return Collected(
            tool=self.tool,
            rows=rows,
            extra={},
            input_rows=len(rows) + sum(dropped.values()),
            dropped=dropped,
        )

    def _row(
        self,
        ctx: RunContext,
        system: str,
        kind: Kind,
        reason: str,
        call: SystemCall | None,
    ) -> Row:
        claim = {
            "schema": CLAIM_SCHEMA,
            "type": "blood_group",
            "system": system,
            "genotype": None if call is None else call.genotype,
            "phenotype": None if call is None else call.phenotype,
        }
        statement: Fields = {
            "run_id": ctx.run_id,
            "statement_id": f"bloodgroup:{system}",
            "kind": str(kind),
            "vid": None,
            "gene": None,
            "claim": json.dumps(claim, ensure_ascii=False),
            "validity": str(Validity.LIMITED),
            "source_level": SOURCE_LEVEL,
            "reason": reason,
            "effect": None,
            "actionable": False,
            "topic": str(Topic.TRAITS),
            "source": BLOODGROUP_TOOL,
            "source_id": system,
            "source_version": SOURCE_VERSION,
            "points": None,
            "exclusion": False,
        }
        return Row(None, None, statement)
