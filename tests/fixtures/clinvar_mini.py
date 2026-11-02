"""A miniature ClinVar GRCh38 VCF, written into a test's tmp_path from the records a test
gives. VCF files are never committed: `write_clinvar` builds them at test time.
"""

from pathlib import Path

import pysam

type Record = tuple[str, int, int, str, str, dict[str, str]]

_HEADER = """\
##fileformat=VCFv4.1
##fileDate=2026-10-12
##source=ClinVar
##reference={assembly}
{contigs}\
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier submitted for germline classifications, e.g. MedGen:NNNNNN">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Aggregate germline classification">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant: symbol:id|symbol:id">
##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

_CONTIGS = ("1", "2", "3", "17", "X", "MT", "NT_187361.1")


def _line(record: Record) -> str:
    chrom, pos, vid, ref, alt, info = record
    fields = ";".join(f"{key}={value}" for key, value in info.items())
    return f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{fields}\n"


def write_clinvar(path: Path, assembly: str, records: list[Record]) -> Path:
    """Write `records` as a bgzipped, tabix-indexed ClinVar-style VCF at `path` (.vcf.gz)."""
    contigs = "".join(f"##contig=<ID={name}>\n" for name in _CONTIGS)
    plain = path.with_suffix("")
    order = {name: i for i, name in enumerate(_CONTIGS)}
    rows = sorted(records, key=lambda r: (order[r[0]], r[1]))
    plain.write_text(
        _HEADER.format(assembly=assembly, contigs=contigs)
        + "".join(_line(r) for r in rows)
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    plain.unlink()
    return path
