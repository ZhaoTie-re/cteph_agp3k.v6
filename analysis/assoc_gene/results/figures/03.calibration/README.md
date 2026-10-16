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

λ_GC spans 1.16–1.24 over the 6 cells. Read across a row: the cohorts are nested and differ almost entirely in controls, so a trend as the sample set widens is a statement about the controls added, not a replication. Read down a column: the variant sets are nested, and the wider set's extra genes carry only LOW-impact variants. Within a cell the two curves share genes, samples and denominator, so their difference is the statistic. λ is reported, not targeted: the statistic is a count over a handful of minor alleles, and its median-based λ is inflated by the two-variant genes (table column lambda_gc_nvar2 against lambda_gc_nvar_ge3); the residual 1.1–1.2 over larger sets is consistent with the platform confound and cannot be separated from it here (METHODS §4.4).

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | 3 |
| variant sets | 2 |
| cells drawn | 6 |
| lambda_GC range | 1.16–1.24 |

## Full statistics

**Denominator and calibration per cell**

| cohort | stratum | n_genes_mapped | bonferroni_threshold | lambda_gc_cmc | lambda_gc_nvar2_cmc | lambda_gc_nvar_ge3_cmc | lambda_gc_skato | lambda_gc_nvar2_skato | lambda_gc_nvar_ge3_skato |
|---|---|---|---|---|---|---|---|---|---|
| narrow_mainland | low_moderate_high | 14,848 | 3.367e-06 | 1.176 | 1.353 | 1.15 | 1.232 | 1.276 | 1.22 |
| intermediate_mainland | low_moderate_high | 15,219 | 3.285e-06 | 1.189 | 1.468 | 1.151 | 1.243 | 1.333 | 1.229 |
| full_mainland | low_moderate_high | 15,829 | 3.159e-06 | 1.176 | 1.432 | 1.142 | 1.221 | 1.174 | 1.231 |
| narrow_mainland | moderate_high | 11,192 | 4.467e-06 | 1.165 | 1.339 | 1.119 | 1.22 | 1.232 | 1.215 |
| intermediate_mainland | moderate_high | 11,720 | 4.266e-06 | 1.163 | 1.384 | 1.103 | 1.184 | 1.215 | 1.169 |
| full_mainland | moderate_high | 12,554 | 3.983e-06 | 1.167 | 1.407 | 1.096 | 1.172 | 1.147 | 1.181 |

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
