"""RFC1's motif analysis (slice 3e): ExpansionHunter's call vocabulary at commit 9e93322f4a.

The labels come from `label()` in ehunter/locus/RFC1Status.hh ("normal", "potential
carrier", "carrier", "affected"), the descriptions from `getRFC1Status` in
ehunter/locus/RFC1MotifAnalysis.cpp (lines 631-675), which ehunter/io/JsonWriter.cpp writes
at `LocusResults.RFC1.Variants.RFC1.RFC1MotifAnalysis.{Call, Description}`; `read_rfc1`
reads them.
"""

import json
from collections.abc import Mapping
from pathlib import Path

from genokeeper.repeats.loci import RepeatCall, Rfc1Call

RFC1_LOCUS = "RFC1"
RFC1_MOTIF = "AARRG"
# Every (label, description) pair `getRFC1Status` returns, in source order.
RFC1_PAIRS: tuple[tuple[str, str], ...] = (
    ("normal", "no expanded allele"),
    ("normal", "expanded allele may exist but not observed"),
    ("carrier", "1 expanded pathogenic allele, 1 short reference allele"),
    ("normal", "1 expanded benign allele, 1 short reference allele"),
    ("normal", "2 expanded alleles may exist but not observed"),
    ("affected", "2 expanded pathogenic alleles"),
    ("carrier", "1 expanded pathogenic allele, 1 expanded benign allele"),
    ("potential carrier", "2 expanded alleles with >30% pathogenic kmers"),
    ("normal", "2 expanded alleles (possibly reference)"),
)


def read_rfc1(json_path: Path, calls: Mapping[str, RepeatCall]) -> Rfc1Call | None:
    """RFC1's motif call from ExpansionHunter's JSON `json_path`
    (`LocusResults.RFC1.Variants.RFC1.RFC1MotifAnalysis.{Call, Description}`), with the
    alleles, coverage and low depth of `calls["RFC1"]` (`read_calls`); None when there is
    no RFC1 call (a no-call, not a failure). `ValueError` naming the file when the analysis
    is missing: the binary is then not the one built."""
    call = calls.get(RFC1_LOCUS)
    if call is None:
        return None
    try:
        results = json.loads(json_path.read_text())["LocusResults"]
        found = results[RFC1_LOCUS]["Variants"][RFC1_LOCUS]["RFC1MotifAnalysis"]
        label, description = found["Call"], found["Description"]
    except KeyError as e:
        raise ValueError(f"{json_path}: no RFC1 motif analysis (key {e})") from e
    return Rfc1Call(label, description, call.alleles, call.coverage, call.low_depth)
