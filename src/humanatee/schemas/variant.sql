CREATE TABLE IF NOT EXISTS Variant (
    variant_id TEXT,
    source TEXT,
    chrom TEXT,
    start INTEGER,
    end INTEGER,
    PRIMARY KEY (variant_id, source)
);
