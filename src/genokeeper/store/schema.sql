-- The findings store (spec: "The store"). Every statement is idempotent: connect() runs this
-- file on each open. Enumerated VARCHAR columns are CHECKed against the values of the StrEnum
-- classes in genokeeper/store/db.py; tests/test_store.py fails when the two drift apart.

-- The allele as one run read it: vid = normalised GRCh38:chr:pos:ref:alt, with the rsID of
-- that run's own ClinVar release (another release may name the same allele otherwise).
CREATE TABLE IF NOT EXISTS variant (
    run_id VARCHAR,
    vid VARCHAR,
    rsid VARCHAR,
    PRIMARY KEY (run_id, vid)
);

-- The user's genotype at a variant in one run. p_wrong is the genotype error rate measured on
-- HG003 (DNBSEQ-T7) for the call's cell, [p_wrong_lo, p_wrong_hi] its Jeffreys 95% interval;
-- all three are null where nothing was measured (outside the calibrated regions, chrM).
CREATE TABLE IF NOT EXISTS call (
    run_id VARCHAR,
    vid VARCHAR,
    gt VARCHAR,
    dosage VARCHAR,
    p_wrong DOUBLE,
    p_wrong_lo DOUBLE,
    p_wrong_hi DOUBLE,
    level VARCHAR NOT NULL CHECK (level IN ('high', 'medium', 'low')),
    gq INT,
    dp INT,
    ad INT[],
    vaf DOUBLE,
    row_type VARCHAR NOT NULL
        CHECK (row_type IN ('PASS', 'RefCall', 'ref_block', 'NoCall', 'none')),
    flags VARCHAR[],
    warnings VARCHAR[],
    PRIMARY KEY (run_id, vid)
);

-- A claim about a variant or a gene, with its validity and its source.
CREATE TABLE IF NOT EXISTS statement (
    run_id VARCHAR,
    statement_id VARCHAR,
    kind VARCHAR NOT NULL CHECK (kind IN ('finding', 'watch', 'non_finding', 'unassessed')),
    vid VARCHAR,
    gene VARCHAR,
    claim JSON,
    validity VARCHAR CHECK (validity IN ('strong', 'moderate', 'limited', 'conflicting')),
    source_level VARCHAR,
    reason VARCHAR,
    effect JSON,
    actionable BOOLEAN,
    topic VARCHAR CHECK (topic IN ('risks', 'carrier', 'drugs', 'traits', 'wellness',
        'ancestry', 'research')),
    source VARCHAR,
    source_id VARCHAR,
    source_version VARCHAR,
    points DOUBLE,
    -- An important negative (the closing design's "Where each kind of result goes"): a
    -- non-finding flagged here is shown as an exclusion, with the confirming test its claim
    -- names (`confirm`; the store writer refuses one without it). Slices 3e and 3f set it and
    -- 3f's dashboard block reads it; 3c sets it nowhere. Absent from a schema-6 snapshot.
    exclusion BOOLEAN NOT NULL DEFAULT false,
    PRIMARY KEY (run_id, statement_id),
    CHECK (kind NOT IN ('finding', 'unassessed') OR (validity IS NOT NULL
        AND source_level IS NOT NULL AND source IS NOT NULL AND source_version IS NOT NULL
        AND topic IS NOT NULL)),
    -- A claim waits on a genotype at its variant, or on its gene's call (an unresolved
    -- pharmacogenomic gene, slice 1b): an unassessed claim names one or the other, or, with
    -- neither (a polygenic score below its overlap, slice 1d), carries its reason.
    CHECK (kind <> 'unassessed' OR vid IS NOT NULL OR gene IS NOT NULL
        OR reason IS NOT NULL)
);

-- A pharmacogenomic gene's call in one run (slice 1b): the diplotype, phenotype and activity
-- score PharmCAT reported (from our calls or an outside caller), every caller's own call
-- (caller -> call), whether two callers agreed (null with one), and the call's reliability:
-- high when two callers agree and no defining position is missing, medium for one caller
-- with none missing, low otherwise. No p_wrong: nothing measured it. No dosing text.
CREATE TABLE IF NOT EXISTS gene_call (
    run_id VARCHAR,
    gene VARCHAR,
    diplotype VARCHAR,
    phenotype VARCHAR,
    activity_score DOUBLE,
    callers JSON,
    agreed BOOLEAN,
    level VARCHAR CHECK (level IN ('high', 'medium', 'low')),
    warnings VARCHAR[],
    PRIMARY KEY (run_id, gene)
);

-- A gene's coverage in one run (slice 3c): over its MANE v1.5 CDS +/- [coverage] cds_margin
-- (`bases`), the share at MAPQ >= [coverage] mapq covered by min_depth_per_copy reads per
-- expected copy (`callable`), the mean depth, the shares inside GIAB's difficult strata and
-- outside the calibrated benchmark, the stage-3 source that owns the gene's answer (null
-- for none), and whether `callable` reaches [coverage] min_callable.
CREATE TABLE IF NOT EXISTS coverage (
    run_id VARCHAR,
    gene VARCHAR,
    bases INTEGER,
    callable DOUBLE,
    mean_depth DOUBLE,
    difficult DOUBLE,
    outside_benchmark DOUBLE,
    owner VARCHAR,
    reliable BOOLEAN,
    PRIMARY KEY (run_id, gene)
);

-- One analysis run: its inputs with checksums, tool and database versions.
CREATE TABLE IF NOT EXISTS run (
    run_id VARCHAR PRIMARY KEY,
    started_at TIMESTAMP,
    inputs JSON,
    tools JSON
);

-- Per-tool row accounting: every input row is a finding, a watch item, a non-finding, an
-- unassessed claim or dropped with a reason.
CREATE TABLE IF NOT EXISTS accounting (
    run_id VARCHAR,
    tool VARCHAR,
    input_rows BIGINT,
    findings BIGINT,
    watch BIGINT,
    non_findings BIGINT,
    unassessed BIGINT,
    dropped JSON,
    PRIMARY KEY (run_id, tool)
);

-- The schema version that wrote the store, one row (db.SCHEMA_VERSION; connect() writes it
-- into a new store and refuses one another version wrote).
CREATE TABLE IF NOT EXISTS meta (
    schema_version INTEGER NOT NULL
);

-- Findings with the user's call; gene-level findings (no vid) keep null call columns.
CREATE OR REPLACE VIEW finding AS
SELECT s.*, c.* EXCLUDE (run_id, vid)
FROM statement AS s
LEFT JOIN call AS c ON c.run_id = s.run_id AND c.vid = s.vid
WHERE s.kind = 'finding';

-- Claims that would be findings were the user's genotype at the site known (its dosage is
-- unknown, or it has no call row), with the call: never a finding, never a result.
CREATE OR REPLACE VIEW unassessed AS
SELECT s.*, c.* EXCLUDE (run_id, vid)
FROM statement AS s
LEFT JOIN call AS c ON c.run_id = s.run_id AND c.vid = s.vid
WHERE s.kind = 'unassessed';
