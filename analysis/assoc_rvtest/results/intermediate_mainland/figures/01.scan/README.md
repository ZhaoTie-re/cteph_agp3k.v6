# Gene-based scans — intermediate_mainland

**Figures:** `01.scan/scan.<peak>.png` — 6 locus figure(s): 0 genome-wide, 0 suggestive.

One document covers the whole family: the panels, the reading and the limits are properties of the figure, identical for every locus. Only the numbers differ, and they are tabulated below, one row per locus.

## The question this figure answers

Is each scan calibrated, and which genes does it single out?

## Panels

**(a) Gene-based association across the genome**

One point per gene, at its midpoint, banded by chromosome. The two lines are the P-values that Benjamini-Hochberg implies for FDR 0.05 and 0.10 IN THIS SCAN. They are not shared constants: BH adapts to the number of genes tested, and that ranges from ~160 under HIGH to ~13,000 under LOW+MODERATE+HIGH, so the implied P differs by roughly an order of magnitude between strata. Each figure therefore prints its own.

**(b) Quantile-quantile plot**

The same gene P-values against a uniform null, with a Beta(i, n-i+1) 95% band and the gene-level genomic-control inflation factor. Read over GENES, not variants — it measures whether the gene-level null is calibrated.

## Interpretation

The scans are the shape of the data, not a list of results. What a reader should take from them is the calibration in (b) and the handful of genes that stand clear of the BH line in (a); the per-gene figures are where the carriers behind each of those become visible.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| loci | 6 |
| genome-wide loci | 0 |
| suggestive loci | 0 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-6 |

## Full statistics

**Per-locus values**

| peak | cohort | impact stratum | test method | genes tested | \lambda_{\mathrm{GC}} | smallest P | top gene | P at FDR 0.05 | P at FDR 0.10 | genes at FDR < 0.05 | genes at FDR < 0.10 | genes named |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| high.cmc | intermediate_mainland | high | cmc | 201 | 0.9793 | 0.004674 | CCDC178 | — | — | 0 | 0 | 0 |
| high.skato | intermediate_mainland | high | skato | 201 | 1.018 | 0.001859 | ABCC12 | — | — | 0 | 0 | 0 |
| low_moderate_high.cmc | intermediate_mainland | low_moderate_high | cmc | 13,326 | 1.132 | 2.759e-06 | SCAMP3 | 2.759e-06 | 2.825e-05 | 1 | 4 | 4 |
| low_moderate_high.skato | intermediate_mainland | low_moderate_high | skato | 13,326 | 1.236 | 3.483e-07 | PIK3C2G | 3.045e-05 | 1.117e-04 | 9 | 15 | 12 |
| moderate_high.cmc | intermediate_mainland | moderate_high | cmc | 8,962 | 1.088 | 7.197e-06 | STBD1 | — | 7.197e-06 | 0 | 1 | 1 |
| moderate_high.skato | intermediate_mainland | moderate_high | skato | 8,962 | 1.167 | 5.641e-07 | PIK3C2G | 2.791e-05 | 1.601e-04 | 7 | 15 | 12 |

## How to read it

1. Read lambda first. A scan whose gene-level lambda is far from 1 has an inflated null, and its tier counts are correspondingly optimistic.
2. Compare strata down the table below rather than across the figures: the same tier label means a different P in each stratum, and the table states both.
3. A gene named in red cleared FDR 0.05; grey cleared only 0.10.

## What this figure does *not* establish

- It cannot separate real polygenic signal from residual confounding. Under this design cases and controls differ in sequencing platform and depth, and the minAC inherited from tuning.rv removes the genome-wide burden artefact, not every per-gene one.
- A tier is not a finding. There is no independent cohort for this phenotype.

## Symbols

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
