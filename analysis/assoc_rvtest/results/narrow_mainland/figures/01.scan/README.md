# Gene-based scans — narrow_mainland

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
| cohort | narrow_mainland |
| loci | 6 |
| genome-wide loci | 0 |
| suggestive loci | 0 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-6 |

## Full statistics

**Per-locus values**

| peak | cohort | impact stratum | test method | genes tested | \lambda_{\mathrm{GC}} | smallest P | top gene | P at FDR 0.05 | P at FDR 0.10 | genes at FDR < 0.05 | genes at FDR < 0.10 | genes named |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| high.cmc | narrow_mainland | high | cmc | 153 | 0.9194 | 0.006728 | COQ8A | — | — | 0 | 0 | 0 |
| high.skato | narrow_mainland | high | skato | 153 | 1.054 | 0.002591 | COQ8A | — | — | 0 | 0 | 0 |
| low_moderate_high.cmc | narrow_mainland | low_moderate_high | cmc | 12,698 | 1.123 | 4.102e-06 | STBD1 | — | 4.102e-06 | 0 | 1 | 1 |
| low_moderate_high.skato | narrow_mainland | low_moderate_high | skato | 12,698 | 1.236 | 8.458e-07 | PIK3C2G | 3.773e-05 | 9.156e-05 | 11 | 14 | 12 |
| moderate_high.cmc | narrow_mainland | moderate_high | cmc | 8,292 | 1.105 | 1.316e-05 | STBD1 | — | — | 0 | 0 | 0 |
| moderate_high.skato | narrow_mainland | moderate_high | skato | 8,292 | 1.216 | 4.132e-07 | PIK3C2G | 2.132e-05 | 7.193e-05 | 8 | 9 | 9 |

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
