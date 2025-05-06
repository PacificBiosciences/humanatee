CREATE TABLE IF NOT EXISTS Allele (
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    allele_index INTEGER NOT NULL,
    sequence TEXT,
    PRIMARY KEY(variant_id, source, allele_index),
    FOREIGN KEY(variant_id, source) REFERENCES Variant(variant_id, source) ON DELETE CASCADE
    );


CREATE TABLE IF NOT EXISTS SampleAllele (
    sample_id TEXT,
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    genotype_index INTEGER NOT NULL,
    allele_index INTEGER NOT NULL,
    PRIMARY KEY(sample_id, variant_id, source, genotype_index),
    FOREIGN KEY(variant_id, source) REFERENCES Variant(variant_id, source) ON DELETE CASCADE
    );

CREATE TABLE IF NOT EXISTS Gene (
    gene TEXT PRIMARY KEY NOT NULL
    );
