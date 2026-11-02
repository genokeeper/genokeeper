"""RFC1's motif analysis as read from ExpansionHunter's HG002 output (`read_rfc1`,
tests/fixtures/hg002)."""

import json
from pathlib import Path

import pytest

from genokeeper.repeats.expansionhunter import read_calls
from genokeeper.repeats.rfc1 import RFC1_LOCUS, RFC1_PAIRS, read_rfc1

FIXTURES = Path(__file__).parent / "fixtures" / "hg002"
VCF = FIXTURES / "hg002.eh3e.vcf.txt"
JSON = FIXTURES / "hg002.eh3e.json"


def test_read_rfc1_hg002() -> None:
    calls = read_calls(VCF, JSON)
    found = read_rfc1(JSON, calls)
    assert found is not None
    assert (found.call, found.description) == ("normal", "no expanded allele")
    assert (found.call, found.description) in RFC1_PAIRS
    assert [a.repeats for a in found.alleles] == [12, 12]
    assert found.alleles == calls[RFC1_LOCUS].alleles
    assert (found.coverage, found.low_depth) == (
        calls[RFC1_LOCUS].coverage,
        calls[RFC1_LOCUS].low_depth,
    )


def test_read_rfc1_without_analysis_refuses(tmp_path: Path) -> None:
    data = json.loads(JSON.read_text())
    del data["LocusResults"][RFC1_LOCUS]["Variants"][RFC1_LOCUS]["RFC1MotifAnalysis"]
    found = tmp_path / "eh.json"
    found.write_text(json.dumps(data))
    calls = read_calls(VCF, found)
    with pytest.raises(ValueError, match=str(found)):
        read_rfc1(found, calls)


def test_read_rfc1_without_a_call_is_none() -> None:
    """A missing RFC1 call (`REPCN` `.`) is a no-call, not a failure: the run states it
    unassessed "no call at this locus" and the calibration withholds RFC1 alone."""
    calls = read_calls(VCF, JSON)
    del calls[RFC1_LOCUS]
    assert read_rfc1(JSON, calls) is None
