"""The variant id (`vid`): `GRCh38:<chr>:<pos>:<REF>:<ALT>`, one biallelic allele."""


def ucsc_contig(chrom: str) -> str:
    """`1` -> `chr1`, `MT`/`chrMT` -> `chrM`; names already in UCSC style pass through."""
    name = chrom.removeprefix("chr")
    if name == "MT":
        name = "M"
    return f"chr{name}"


def make_vid(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Return the GRCh38 variant id; refuse an empty, missing or multi-allelic allele.

    The caller passes a record that is already normalised (ClinVar's VCF is left-aligned and
    biallelic); this function does not left-align.
    """
    for allele in (ref, alt):
        if not allele or allele == "." or "," in allele:
            raise ValueError(f"not one biallelic allele: {chrom}:{pos} {ref!r}>{alt!r}")
    return f"GRCh38:{ucsc_contig(chrom)}:{pos}:{ref.upper()}:{alt.upper()}"
