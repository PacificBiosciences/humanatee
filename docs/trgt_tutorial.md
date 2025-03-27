# TRGT Annotation Tutorial <!-- omit in toc -->

- [Common use cases](#common-use-cases)
  - [Identifying pathogenic expansions using defined thresholds](#identifying-pathogenic-expansions-using-defined-thresholds)
  - [Identifying population outliers](#identifying-population-outliers)
  - [Additional annotation options](#additional-annotation-options)
  - [Example command](#example-command)

## Common use cases

### Identifying pathogenic expansions using defined thresholds

The motif count of known pathogenic repeat loci can be compared to defined pathogenic thresholds to identify pathogenic expansions. To do this, you will need to run [TRGT](https://github.com/PacificBiosciences/trgt/tree/main) using a repeats catalog that has a corresponding annotations file. We have provided a [pathogenic repeats catalog](../resources/pathogenic_repeats.hg38.bed) and a [corresponding annotations file](../resources/pathogenic_repeats.hg38.annotations.tsv) for this purpose. However, you can also create your own repeats catalog and annotations file with your own defined loci and thresholds.

```bash
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.phased.vcf \
  --sample-id HG002 \
  --pathogenic-tsv pathogenic_repeats.hg38.annotations.tsv \
  --prefix HG002.GRCh38.trgt.humanatee
```

### Identifying population outliers

For repeats without defined pathogenic motif count thresholds, you can compare the repeat allele lengths (field `AL`) to a reference population to identify outliers. Outliers are defined here as alleles with a length that is larger than the 99th percentile of population. This requires first building a database of repeat allele lengths from a reference population using the same repeats catalog that you use for the sample being annotated.

#### Build population database <!-- omit in toc -->

```bash
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
  --cohort-id <cohort_id> \
  --prefix <cohort_id>.trgt

# output = <cohort_id>.trgt.db
```

#### Annotate sample with population database <!-- omit in toc -->

```bash
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vcf \
  --sample-id HG002 \
  --population-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee
```

### Additional annotation options

#### Gene consequence <!-- omit in toc -->

Gene consequence annotations can be used to annotate and prioritize variants based on their predicted impact on gene function. To do this, first run [VEP](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html) with the following command. Additional fields can be added to the `--fields` option as needed, but the fields listed below are required for Humanatee to function properly.

```bash
vep --cache \
  --input_file HG002.GRCh38.trgt.sorted.vcf \
  --vcf \
  --output_file HG002.GRCh38.trgt.sorted.vep.vcf \
  --fields "ALLELE_NUM,SYMBOL,Gene,BIOTYPE,Consequence,IMPACT" \
  --symbol \
  --allele_number \
  --no_intergenic \
  --biotype \
  --allow_non_variant \
  --dont_skip \
  --pick_allele_gene \
  --fork 4 \
  --no_stats
```

The resulting VCF can then be annotated with Humanatee using the `--vcf` option as usual.

#### Gene lookups <!-- omit in toc -->

If the VCF has been annotated with gene consequence annotations using VEP, you can also use the `--gene-lookup` option to annotate the TRGT output with gene information. Multiple gene lookups can be specified by repeating the `--gene-lookup` option.

Example:

```bash
--gene-lookup gene_lookups/pli.tsv
```

A gene lookup file should be a headerless 2-column tab-separated file with gene name in the first column and the lookup value in the second column. For example, the `pli.tsv` file might look like this:

```bash
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

And would be run like this on a VEP annotated VCF::

```bash
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vep.vcf \
  --sample-id HG002 \
  --population-db population.trgt.db \
  --gene-lookup gene_lookups/pli.tsv \
  --prefix HG002.GRCh38.trgt.humanatee
```

Additionally, the `--phenotype` option can be used to supply a gene lookup file of phenotype scores. This information will be used to prioritize variants in the output such that variants in genes with a higher phenotype score will be ranked higher.

```bash
--phenotype HG002_phrank.tsv \
```

#### Coverage dropouts <!-- omit in toc -->

Coverage dropouts can be annotated by providing a 4-column BED file with the following columns: chromosome, start, end, TRGT info (e.g. `ID=chr1_144502_144606_TA;MOTIFS=TA;STRUC=(TA)n`), and dropout status. Dropout status should be `HaplotypeDropout` if the region is a haplotype dropout or `FullDropout` if the region is a coverage dropout.

```bash
--dropouts HG002.GRCh38.trgt.dropouts.txt
```

#### De novo expansions <!-- omit in toc -->

De novo expansions can be annotated by providing the output from [trgt-denovo](https://github.com/PacificBiosciences/trgt-denovo). Both duo and trio modes are supported. If the `--denovo` option is used, it's necessary to indicate which parents were use to run either duo or trio mode.

```bash
# trio mode
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vep.vcf \
  --sample-id HG002 \
  --population-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee \
  --denovo HG002.GRCh38.trgt.denovo.tsv \
  --mother \
  --father

# duo mode with only the mother
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vep.vcf \
  --sample-id HG002 \
  --population-db <cohort_id>.trgt.db \
  --prefix HG002.GRCh38.trgt.humanatee \
  --denovo HG002.GRCh38.trgt.denovo.tsv \
  --mother
```

### Example command

To use all available annotation options and produce all possible outputs, the command might look like this:

```bash
humanatee trgt \
  --vcf HG002.GRCh38.trgt.sorted.vep.vcf \
  --sample-id HG002 \
  --pathogenic-tsv pathogenic_repeats.hg38.annotations.tsv \
  --population-db population.trgt.db \
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
