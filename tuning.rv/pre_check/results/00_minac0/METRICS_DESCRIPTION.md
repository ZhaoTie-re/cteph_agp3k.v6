# Metrics Output Documentation

Generated on 2026-05-13 19:05:20

## 1. Sample Metrics (sample_metrics.txt)
**Dimensions:** 3,060 samples x 9 columns

| Column Name | Description |
| :--- | :--- |
| **SampleID** | Unique identifier for the sample (from input VCF/Info file). |
| **Group** | Sample Group (Case/Control) determined by info file outcome column. |
| **TargetDP** | Target Depth (e.g., sequencing depth target from Info file). |
| **MeanDP** | Mean Depth (Actual average coverage from Info file). |
| **SNumHomRef** | Count of sites where the sample is Homozygous Reference. |
| **SNumHet** | Count of sites where the sample is Heterozygous. |
| **SNumHomAlt** | Count of sites where the sample is Homozygous Alternative (relative to REF genome). |
| **SMissCount** | Count of sites where the sample has missing genotype calls. |
| **SMinAC** | Sum of Minor Allele Counts (calculated as N_het + 2 * N_hom_minor using force-accepted major allele). |

## 2. Variant Metrics (variant_metrics.txt)
**Dimensions:** 19,720,391 variants x 12 columns

| Column Name | Description |
| :--- | :--- |
| **#CHROM** | Chromosome identifier. |
| **POS** | Genomic position. |
| **VariantID** | Variant identifier (format: CHROM:POS:REF:ALT). |
| **REF** | Reference allele sequence. |
| **ALT** | Alternative allele sequence. |
| **RefAC** | Reference Allele Count. Total number of REF alleles observed. |
| **AltAC** | Alternative Allele Count. Total number of ALT alleles observed. |
| **MinAC** | Minor Allele Count. The count of the less frequent allele (min(RefAC, AltAC)). |
| **VNumHomRef** | Number of samples that are Homozygous Reference (0/0). |
| **VNumHet** | Number of samples that are Heterozygous (0/1). |
| **VNumHomAlt** | Number of samples that are Homozygous Alternative (1/1). |
| **VMissCount** | Number of samples where the genotype call is missing. |
