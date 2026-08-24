# Gene-based results across the three cohorts

**Figure file:** `cohort_compare_genes.png`

## The question this figure answers

Does a gene-based hit survive relaxing the ancestry filter, and are the scans themselves calibrated?

## Panels

**(a) Scan calibration**

Gene-level genomic-control inflation factor for each impact stratum and test method, per cohort. Computed over the gene P-values, not variant ones, so it measures whether the GENE-level null is calibrated.

**(b) Every tiered gene, in every cohort**

The union of genes reaching a tier in any cohort, each shown in all three. A dot plot of -log10 P rather than an effect forest: the gene-level test reports a P, and the burden odds ratio belongs beside its allele counts, which is where the per-gene figures put it.

**(c) Robustness to the ancestry filter**

How many of the three cohorts call each tiered gene a hit.

## Interpretation

Under a nested design the useful signal in a cross-cohort comparison is instability, not agreement. A gene that appears only when the ancestry filter is relaxed is a candidate artefact; one that holds across all three has survived the only internal check available, which is weaker than replication and should be described as such.

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | narrow_mainland, intermediate_mainland, full_mainland |
| stratum shown | moderate_high |
| method shown | skato |
| genes reaching a tier | 30 |
| genes shown | 18 |
| significant somewhere | 13 |
| gene-level lambda, min | 0.9194 |
| gene-level lambda, max | 1.236 |
| genes called in 1 of 3 cohorts | 16 |
| genes called in 2 of 3 cohorts | 8 |
| genes called in 3 of 3 cohorts | 6 |

## How to read it

1. Read (c) as a sensitivity check, never as replication. The cohorts are nested and share every case, so a gene appearing in all three has not been confirmed by independent data — it has been shown not to depend on the ancestry cut.
2. A gene called in the widest cohort ONLY is the one to distrust: what full adds over narrow is exactly the ancestry-outlier samples.
3. Check (a) before (b). A stratum whose lambda is far from 1 has a mis-calibrated gene-level null, and its tiers are correspondingly optimistic.

## What this figure does *not* establish

- It cannot establish that any gene is real. No independent cohort exists for this phenotype, so nothing here is external validation.
- Tiers come from Benjamini-Hochberg within each scan, so a tier in the HIGH stratum (few genes tested) and one in MODERATE+HIGH are not the same evidence.

## Symbols

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
