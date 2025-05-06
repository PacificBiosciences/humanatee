# TRGT Annotation Tutorial <!-- omit in toc -->

- [Common use cases](#common-use-cases)
  - [Identifying pathogenic expansions using defined thresholds](#identifying-pathogenic-expansions-using-defined-thresholds)
  - [Identifying population outliers](#identifying-population-outliers)
  - [Additional annotation options](#additional-annotation-options)
  - [Example command](#example-command)

## Common use cases

### Identifying pathogenic expansions using defined thresholds

The motif count of known pathogenic repeat loci can be compared to defined pathogenic thresholds to identify pathogenic expansions. To do this, you will need to run [TRGT](https://github.com/PacificBiosciences/trgt/tree/main) using a repeats catalog that has a corresponding annotations file. We have provided a [pathogenic repeats catalog](../resources/pathogenic_repeats.hg38.bed) and a [corresponding annotations file](../resources/pathogenic_repeats.hg38.annotations.tsv) for this purpose. However, you can also create your own repeats catalog and annotations file with your own defined loci and thresholds.

```text
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.phased.vcf \
  --sample-id HG002 \
  --pathogenic-tsv pathogenic_repeats.hg38.annotations.tsv \
  --prefix HG002.GRCh38.trgt.humanatee
```

### Identifying population outliers

For repeats without defined pathogenic motif count thresholds, you can compare the repeat allele lengths (field `AL`) to a reference population to identify outliers. Outliers are defined here as alleles with a length that is larger than the 99th percentile of population. This requires first building a database of repeat allele lengths from a reference population using the same repeats catalog that you use for the sample being annotated.

#### Build population database <!-- omit in toc -->

```text
# merge individual VCFs into a single multi-sample VCF
trgt merge \
  --vcf *.trgt.sorted.vcf.gz \
  --genome human_GRCh38_no_alt_analysis_set.fasta \
  --output merged.trgt.vcf.gz

# index the merged VCF
bcftools index \
  --threads 3 \
  --tbi merged.trgt.vcf.gz

# build the population database
humanatee trgtpop \
  --vcf merged.trgt.vcf.gz \
  --prefix <cohort_id>.trgt

# output = <cohort_id>.trgt.db
```

#### Annotate sample with population database <!-- omit in toc -->

```text
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --annotation-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee
```

### Additional annotation options

#### Gene consequence <!-- omit in toc -->

Gene consequence annotations can be used to annotate and prioritize variants based on their predicted impact on gene function. This requires first building a database of gene consequence annotations using the same repeats catalog that you use for the sample being annotated. If you already have a population database, consequence annotations should be added to it by specifying the path to the existing population database with the `--db` flag.

```text
# download gene annotations from gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz

# build consequence database
humanatee refanno \
  --vcf merged.trgt.vcf.gz \
  --source trgt \
  --prefix <cohort_id>.trgt.refanno \
  --annotation gencode.v47.annotation.gff3.gz csq 1000

# if you already have a population database, use the --db flag to add annotations to it
humanatee refanno \
  --db <cohort_id>.trgt.db
  --vcf merged.trgt.vcf.gz \
  --source trgt \
  --prefix <cohort_id>.trgt.refanno \
  --annotation gencode.v47.annotation.gff3.gz csq 1000

# output = <cohort_id>.trgt.refanno.db
```

#### Gene lookups <!-- omit in toc -->

If your `--annotation-db` file has gene consequence annotations, you can also use the `--gene-lookup` option to annotate the TRGT output with gene information. Multiple gene lookups can be specified by repeating the `--gene-lookup` option.

Example:

```text
--gene-lookup gene_lookups/pli.tsv
```

A gene lookup file should be a headerless 2-column tab-separated file with gene name in the first column and the lookup value in the second column. For example, the `pli.tsv` file might look like this:

```text
SPTLC2  0.78
C1R     0.779
TMEM135 0.779
SRPR    0.778
RAP1B   0.91
ADAMTS9 1e-06
RAB18   0.91
B3GALT2 0.896
PIK3R4  0.0245
```

And would be run like this with an `--annotation-db` file that contains gene consequence annotations:

```text
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --annotation-db population.trgt.refanno.db \
  --gene-lookup gene_lookups/pli.tsv \
  --prefix HG002.GRCh38.trgt.humanatee
```

Additionally, the `--phenotype` option can be used to supply a gene lookup file of phenotype scores. This information will be used to prioritize variants in the output such that variants in genes with a higher phenotype score will be ranked higher.

```text
--phenotype HG002_phrank.tsv \
```

#### Coverage dropouts <!-- omit in toc -->

Coverage dropouts can be annotated by providing a 4-column BED file with the following columns: chromosome, start, end, TRGT info (e.g. `ID=chr1_144502_144606_TA;MOTIFS=TA;STRUC=(TA)n`), and dropout status. Dropout status should be `HaplotypeDropout` if the region is a haplotype dropout or `FullDropout` if the region is a coverage dropout.

```text
--dropouts HG002.GRCh38.trgt.dropouts.txt
```

#### De novo expansions <!-- omit in toc -->

De novo expansions can be annotated by providing the output from [trgt-denovo](https://github.com/PacificBiosciences/trgt-denovo). Both duo and trio modes are supported. If the `--denovo` option is used, it's necessary to indicate which parents were use to run either duo or trio mode.

```text
# trio mode
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --annotation-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee \
  --denovo HG002.GRCh38.trgt.denovo.tsv \
  --mother \
  --father

# duo mode with only the mother
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --annotation-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee \
  --denovo HG002.GRCh38.trgt.denovo.tsv \
  --mother
```

### Example command

To use all available annotation options and produce all possible outputs, the command might look like this:

```text
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --pathogenic-tsv pathogenic_repeats.hg38.annotations.tsv \
  --annotation-db population.trgt.refanno.db \
  --gene-lookup gene_lookups/pli.tsv \
  --gene-lookup gene_lookups/clinvar.tsv \
  --phenotype HG002_phrank.tsv \
  --dropouts HG002.GRCh38.trgt.dropouts.txt \
  --denovo HG002.GRCh38.trgt.denovo.tsv \
  --mother \
  --father \
  --prefix HG002.GRCh38.trgt.humanatee \
  --filter \
  --tsv \
  --xlsx \
  --plot
```

For a description of the output files, see [Interpreting TRGT Annotation](trgt_interpretation.md).
