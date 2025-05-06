CREATE TABLE IF NOT EXISTS Consequence (
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    allele_index INTEGER,
    gene TEXT,
    ensembl_id TEXT NOT NULL,
    csq TEXT NOT NULL,
    ranked_consequence TEXT,
    ranked_impact TEXT,
    PRIMARY KEY(variant_id, source, ensembl_id),
    FOREIGN KEY(variant_id, source) REFERENCES Variant(variant_id, source) ON DELETE CASCADE
    );
