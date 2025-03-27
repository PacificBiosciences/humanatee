CREATE TABLE IF NOT EXISTS Locus (
    trid TEXT PRIMARY KEY NOT NULL,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    motifs TEXT NOT NULL,
    struc TEXT NOT NULL,
    lower INTEGER NOT NULL,
    upper INTEGER NOT NULL,
    n_alleles INTEGER NOT NULL,
    recessive_lower INTEGER NOT NULL,
    recessive_upper INTEGER NOT NULL,
    recessive_n_alleles INTEGER NOT NULL,
    allele_length_genotypes TEXT NOT NULL,
    UNIQUE(chrom, start, end, motifs, struc));

CREATE TABLE IF NOT EXISTS Sample (
    sampleid TEXT NOT NULL,
    cohort TEXT NOT NULL,
    UNIQUE(sampleid, cohort));
