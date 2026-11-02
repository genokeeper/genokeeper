"""A synthetic ClinGen gene-disease validity export, shaped like the one ClinGen serves at
https://search.clinicalgenome.org/kb/gene-validity/download (read 2026-10-18): a preamble, a
`+++` rule, the header, another rule, then one quoted row per curation."""

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path

HEADER = [
    "GENE SYMBOL",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "MOI",
    "SOP",
    "CLASSIFICATION",
    "ONLINE REPORT",
    "CLASSIFICATION DATE",
    "GCEP",
]


def write_dosage(
    path: Path,
    scores: Sequence[tuple[str, str]],
    *,
    date: str = "23 Oct,2026",
    triplo: Mapping[str, str] | None = None,
    locations: Mapping[str, str] | None = None,
) -> Path:
    """ClinGen's gene curation list (`ClinGen_gene_curation_list_GRCh38.tsv`, shaped as
    ftp.clinicalgenome.org served it on 2026-10-23) with `(gene, haploinsufficiency score)`
    rows, each gene's triplosensitivity score from `triplo` (default 0) and its location
    from `locations` (default `chr1:1-2`)."""
    header = (
        "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score"
        "\tHaploinsufficiency Description\tTriplosensitivity Score\n"
    )
    ts, where = triplo or {}, locations or {}
    rows = "".join(
        f"{gene}\t{i}\t1p1\t{where.get(gene, 'chr1:1-2')}\t{score}\tdescribed"
        f"\t{ts.get(gene, '0')}\n"
        for i, (gene, score) in enumerate(scores, start=1)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"#ClinGen Gene Curation Results\n#{date}\n#Genomic Locations are reported on"
        f" GRCh38 (hg38): GCF_000001405.36\n{header}{rows}"
    )
    return path


def write_gene_validity(
    path: Path,
    curations: Sequence[tuple[str, str] | tuple[str, str, str]],
    *,
    created: str = "2026-10-18",
    extra_column: bool = False,
    mondo_column: bool = True,
) -> Path:
    """Write `(gene, classification)` or `(gene, classification, MONDO ID)` curations as
    ClinGen's CSV export at `path`. A curation without a MONDO ID gets one of its own
    (`MONDO:<n>`, n its 1-based row), which no test's ClinVar disease names;
    `mondo_column=False` leaves the `DISEASE ID (MONDO)` column out."""
    header = ["EXTRA", *HEADER] if extra_column else list(HEADER)
    if not mondo_column:
        header.remove("DISEASE ID (MONDO)")
    width = len(header)
    rule = ["+" * 9] * width
    pad = [""] * (width - 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        out = csv.writer(f, quoting=csv.QUOTE_ALL, lineterminator="\n")
        out.writerow(["CLINGEN GENE DISEASE VALIDITY CURATIONS", *pad])
        out.writerow([f"FILE CREATED: {created}", *pad])
        out.writerow(
            ["WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity", *pad]
        )
        out.writerow(rule)
        out.writerow(header)
        out.writerow(rule)
        for i, (gene, classification, *mondo) in enumerate(curations):
            row = {
                "GENE SYMBOL": gene,
                "GENE ID (HGNC)": f"HGNC:{i + 1}",
                "DISEASE LABEL": f"{gene}-related disease",
                "DISEASE ID (MONDO)": mondo[0] if mondo else f"MONDO:{i + 1:07d}",
                "MOI": "AD",
                "SOP": "SOP10",
                "CLASSIFICATION": classification,
                "ONLINE REPORT": "https://search.clinicalgenome.org/kb/gene-validity/x",
                "CLASSIFICATION DATE": "2024-01-01T00:00:00.000Z",
                "GCEP": "A Gene Curation Expert Panel",
                "EXTRA": "x",
            }
            out.writerow([row[column] for column in header])
    return path


def write_regions(
    path: Path,
    regions: Sequence[tuple[str, str, str, str]],
    *,
    date: str = "28 Apr,2026",
) -> Path:
    """ClinGen's region curation list (`ClinGen_region_curation_list_GRCh38.tsv`, shaped as
    ftp.clinicalgenome.org served it on 2026-10-30) with `(name, location, haploinsufficiency
    score, triplosensitivity score)` rows."""
    header = (
        "#ISCA ID\tISCA Region Name\tcytoBand\tGenomic Location"
        "\tHaploinsufficiency Score\tHaploinsufficiency Description"
        "\tTriplosensitivity Score\tTriplosensitivity Description\n"
    )
    rows = "".join(
        f"ISCA-{46700 + i}\t{name}\t22q11.21\t{location}\t{hi}\tdescribed\t{ts}"
        "\tdescribed\n"
        for i, (name, location, hi, ts) in enumerate(regions, start=1)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"#ClinGen Region Curation Results\n#{date}\n#Genomic Locations are reported on"
        f" GRCh38 (hg38): GCF_000001405.36\n{header}{rows}"
    )
    return path


def write_recurrent(path: Path, rows: Sequence[tuple[str, int, int, str, str]]) -> Path:
    """ClinGen's recurrent CNV BED (`ClinGen_recurrent_CNV_V2.1-hg38.bed`, shaped as served
    on 2026-10-30) with `(chrom, start, end, name, itemRgb)` rows."""
    lines = [
        'track name="ClinGen_recurrent_CNV_V2.0-hg38" itemRGB="On"',
        "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb",
        *(
            f"{c}\t{s}\t{e}\t{name}\t0\t.\t{s}\t{e}\t{rgb}"
            for c, s, e, name, rgb in rows
        ),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path
