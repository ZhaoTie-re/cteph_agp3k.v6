# Calibration across cohorts and variant sets

**Figure file:** `calibration.png`

## The question this figure answers

Is the null calibrated for both statistics in every cohort × variant-set cell of the design, and does calibration change as controls are added or the impact filter is loosened?

## Panels

**(a) LOW+MODERATE+HIGH, narrow**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

**(b) LOW+MODERATE+HIGH, intermediate**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

**(c) LOW+MODERATE+HIGH, full**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

**(d) MODERATE+HIGH, narrow**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

**(e) MODERATE+HIGH, intermediate**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

**(f) MODERATE+HIGH, full**

Observed against expected −log10 P for CMC burden and SKAT-O over the genes each returned in this cell; expected quantiles are (i − 0.5)/n. The grey band is the pointwise 95 % interval of the uniform null, Beta(i, n − i + 1). The corner text gives each statistic's genomic-control inflation factor, read from scan_qc.tsv. No threshold is drawn on a QQ; each cell's own Bonferroni threshold and gene count are in the table below.

## Interpretation

λ_GC spans 1.15–1.25 over the 6 cells. Read across a row: the cohorts are nested and differ almost entirely in controls, so a trend as the sample set widens is a statement about the controls added, not a replication. Read down a column: the variant sets are nested, and the wider set's extra genes carry only LOW-impact variants. Within a cell the two curves share genes, samples and denominator, so their difference is the statistic. λ is reported, not targeted: the statistic is a count over a handful of minor alleles, and its median-based λ is inflated by the two-variant genes (table column lambda_gc_nvar2 against lambda_gc_nvar_ge3); the residual 1.1–1.2 over larger sets is consistent with the platform confound and cannot be separated from it here (METHODS §4.4).

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | 3 |
| variant sets | 2 |
| cells drawn | 6 |
| lambda_GC range | 1.15–1.25 |

## Full statistics

**Denominator and calibration per cell**

| cohort | stratum | n_genes_mapped | bonferroni_threshold | lambda_gc_cmc | lambda_gc_nvar2_cmc | lambda_gc_nvar_ge3_cmc | lambda_gc_skato | lambda_gc_nvar2_skato | lambda_gc_nvar_ge3_skato |
|---|---|---|---|---|---|---|---|---|---|
| narrow_mainland | low_moderate_high | 14,557 | 3.435e-06 | 1.164 | 1.367 | 1.126 | 1.237 | 1.275 | 1.227 |
| intermediate_mainland | low_moderate_high | 15,009 | 3.331e-06 | 1.185 | 1.504 | 1.139 | 1.246 | 1.371 | 1.229 |
| full_mainland | low_moderate_high | 15,691 | 3.187e-06 | 1.176 | 1.495 | 1.136 | 1.226 | 1.207 | 1.228 |
| narrow_mainland | moderate_high | 10,832 | 4.616e-06 | 1.152 | 1.328 | 1.102 | 1.215 | 1.239 | 1.204 |
| intermediate_mainland | moderate_high | 11,436 | 4.372e-06 | 1.154 | 1.376 | 1.094 | 1.172 | 1.209 | 1.162 |
| full_mainland | moderate_high | 12,378 | 4.039e-06 | 1.181 | 1.44 | 1.111 | 1.164 | 1.173 | 1.161 |

## How to read it

1. Compare the two curves within a cell first; that comparison holds everything but the statistic fixed.
2. Then read across a row (adding controls) and down a column (loosening the impact filter).
3. Never compare heights between cells: each has its own gene family and threshold.

## What this figure does *not* establish

- It is not replication; the cohorts are nested and share nearly all cases.
- It does not establish any association; the platform confound is in every curve.
- λ_GC over non-independent, discrete gene tests is a coarse instrument.

## Model

```
per cell: rvtest CMC burden and SKAT-O, logistic, sex + ancestry PCs, on the cell's own snpEff map; band = Beta(i, n − i + 1) pointwise 95 % interval
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
