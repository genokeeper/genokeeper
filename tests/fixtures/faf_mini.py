"""A miniature gnomAD FAF table, shaped as `genokeeper.evidence.frequency.FafTable` reads it:
one `faf-<chrom>.parquet` per gnomAD chromosome (`vid`, `faf95_max`, `gen_anc`) and a
`manifest.json`. Nothing here is committed as a data file: the table is written at test time.
"""

import json
from collections.abc import Mapping
from pathlib import Path

import duckdb

from genokeeper.evidence.frequency import FAF_FIELD, GNOMAD_CHROMS, GNOMAD_RELEASE

FLOOR = 0.000001


def write_faf_table(
    folder: Path,
    rows: Mapping[str, float] | None = None,
    *,
    floor: float = FLOOR,
    chroms: tuple[str, ...] = GNOMAD_CHROMS,
) -> Path:
    """Write the table of `rows` (vid -> faf95_max, ancestry `nfe`) into `folder`, with a
    manifest listing `chroms` (default: every gnomAD chromosome). Rewrites a table there."""
    folder.mkdir(parents=True, exist_ok=True)
    rows = rows or {}
    entries: dict[str, dict[str, object]] = {}
    with duckdb.connect() as con:
        con.execute("CREATE TABLE faf (vid VARCHAR, faf95_max DOUBLE, gen_anc VARCHAR)")
        if rows:
            con.executemany(
                "INSERT INTO faf VALUES (?, ?, 'nfe')",
                [[v, f] for v, f in rows.items()],
            )
        for chrom in chroms:
            out = folder / f"faf-{chrom}.parquet"
            con.execute(
                f"COPY (SELECT * FROM faf WHERE split_part(vid, ':', 2) = '{chrom}')"
                f" TO '{out}' (FORMAT parquet)"
            )
            entries[chrom] = {
                "url": "https://storage.googleapis.com/gcp-public-data--gnomad/release/"
                f"4.1/vcf/joint/gnomad.joint.v4.1.sites.{chrom}.vcf.bgz",
                "size": 0,
                "md5": "",
                "rows": sum(1 for v in rows if v.split(":")[1] == chrom),
            }
    manifest = {
        "release": GNOMAD_RELEASE,
        "field": FAF_FIELD,
        "floor": floor,
        "chromosomes": entries,
    }
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return folder
