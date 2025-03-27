CREATE TABLE IF NOT EXISTS Variant (
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    priority_rank INTEGER,
    flags TEXT,
    filters TEXT,
    chrom TEXT,
    start INTEGER,
    end INTEGER,
    PRIMARY KEY (variant_id, source)
    );

CREATE TABLE IF NOT EXISTS Allele (
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    allele_index INTEGER NOT NULL,
    sequence TEXT,
    PRIMARY KEY(variant_id, source, allele_index),
    FOREIGN KEY(variant_id, source) REFERENCES Variant(variant_id, source) ON DELETE CASCADE
    );

CREATE TABLE IF NOT EXISTS Consequence (
    variant_id TEXT NOT NULL,
    source TEXT NOT NULL,
    allele_index INTEGER NOT NULL,
    gene TEXT,
    ensembl_id TEXT NOT NULL,
    csq TEXT NOT NULL,
    ranked_consequence TEXT,
    ranked_impact TEXT,
    PRIMARY KEY(variant_id, source, allele_index, ensembl_id),
    FOREIGN KEY(variant_id, source) REFERENCES Variant(variant_id, source) ON DELETE CASCADE,
    FOREIGN KEY(gene) REFERENCES Gene(gene)
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
