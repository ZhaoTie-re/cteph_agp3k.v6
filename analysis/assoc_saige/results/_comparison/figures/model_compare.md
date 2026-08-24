# Reported model against its calibration probe

**Figure file:** `model_compare.png`

## The question this figure answers

How much of the population structure does the GRM absorb on its own, how much do the principal components still remove on top of it, and does any signal depend on them?

## Panels

**(a) Calibration**

lambda_GC = median(chi2)/0.4549 for each model in each cohort, taken from that cohort's own scan QC table rather than recomputed. Right: the quantile-quantile curves behind those numbers, one panel per cohort with both models overlaid and every panel on one shared limit, so the curves are comparable between cohorts as well as between models. The curves are thinned for drawing only — the whole tail is kept and lambda is computed from every P-value.

**(b) Evidence, model against probe**

-log10 P under the reported model on x and the probe on y, one point per locus per cohort. The loci are EXACTLY those cohort_compare draws: diamonds are its panel (b) (genome-wide in at least one cohort), circles its panel (c) (suggestive in all three); colour is the cohort. Both figures read the same `panel` column of lead_crosscohort.tsv, so a locus can be followed between them. The red line is equality: a point on it is a locus the principal components do not touch. Dotted lines mark the genome-wide threshold on each axis, so whether a locus crosses under one model and not the other is readable directly.

**(c) Effect, model against probe**

The same leads by odds ratio, log axes, with dotted lines at OR = 1. Effects that sit on the identity line are unchanged by the covariates even where the evidence moves, which is the usual pattern when the difference is in the standard error rather than in the estimate.

## Interpretation

The two models fit the SAME samples and the SAME full genetic relationship matrix and differ only in whether the principal components are among the covariates. The pair is therefore a direct read-out of what the PCs contribute once relatedness and broad structure are already in the model. A probe whose lambda_GC is already near 1 says the GRM has done the work by itself; a probe that is inflated while the reported model is not says the PCs are removing structure the GRM cannot reach. In panel (b) a lead whose interval is stable across the pair is not an artefact of the covariate choice, while one that moves is a statement about the covariates rather than about the locus. The probe defines no peaks, receives no fine-mapping, conditional analysis or regional plot, and must not be reported as a second result.

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | narrow_mainland, intermediate_mainland, full_mainland |
| models | main, detect |
| peak model | main |
| lambda_GC full_mainland / detect | 1.017 |
| lambda_GC full_mainland / main | 1.021 |
| lambda_GC intermediate_mainland / detect | 1.024 |
| lambda_GC intermediate_mainland / main | 1.038 |
| lambda_GC narrow_mainland / detect | 1.04 |
| lambda_GC narrow_mainland / main | 1.039 |
| leads compared | 16 |
| rows in the forest | 48 |

## Full statistics

**Every peak lead under both models**

| cohort | panel | variant_id | chrom | pos | OR_main | L95_main | U95_main | P_main | OR_detect | L95_detect | U95_detect | P_detect |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| narrow_mainland | b | chr16:53882967:G:A | 16 | 53,882,967 | 1.671 | 1.372 | 2.036 | 1.358e-07 | 1.671 | 1.375 | 2.031 | 9.390e-08 |
| narrow_mainland | c | chr2:232438125:T:C | 2 | 232,438,125 | 0.504 | 0.3698 | 0.687 | 5.539e-06 | 0.5011 | 0.3684 | 0.6816 | 3.836e-06 |
| narrow_mainland | c | chr2:233068979:G:A | 2 | 233,068,979 | 4.04 | 2.113 | 7.723 | 5.195e-06 | 4.043 | 2.123 | 7.7 | 3.850e-06 |
| narrow_mainland | c | chr3:94020684:T:C | 3 | 94,020,684 | 5.248 | 2.587 | 10.65 | 5.773e-07 | 4.355 | 2.215 | 8.565 | 3.596e-06 |
| narrow_mainland | c | chr3:100783010:T:C | 3 | 100,783,010 | 1.906 | 1.445 | 2.515 | 1.970e-06 | 1.939 | 1.475 | 2.549 | 7.373e-07 |
| narrow_mainland | c | chr3:154069965:A:G | 3 | 154,069,965 | 3.762 | 2.076 | 6.814 | 2.820e-06 | 3.878 | 2.137 | 7.036 | 1.537e-06 |
| narrow_mainland | c | chr4:171998307:T:C | 4 | 171,998,307 | 4.386 | 2.276 | 8.453 | 1.734e-06 | 3.947 | 2.079 | 7.491 | 5.868e-06 |
| narrow_mainland | c | chr6:32632861:A:G | 6 | 32,632,861 | 0.6678 | 0.5603 | 0.7958 | 3.173e-06 | 0.6514 | 0.5482 | 0.774 | 5.339e-07 |
| narrow_mainland | c | chr8:132794261:A:T | 8 | 132,794,261 | 3.611 | 1.998 | 6.526 | 4.495e-06 | 3.167 | 1.784 | 5.624 | 2.434e-05 |
| narrow_mainland | c | chr9:104559024:C:T | 9 | 104,559,024 | 1.62 | 1.31 | 2.005 | 4.011e-06 | 1.627 | 1.318 | 2.008 | 2.650e-06 |
| narrow_mainland | c | chr10:1747167:CGTG:C | 10 | 1,747,167 | 0.3379 | 0.2042 | 0.559 | 5.033e-06 | 0.3511 | 0.2127 | 0.5797 | 9.433e-06 |
| narrow_mainland | c | chr14:88019710:G:A | 14 | 88,019,710 | 1.499 | 1.265 | 1.776 | 1.421e-06 | 1.506 | 1.275 | 1.779 | 6.853e-07 |
| narrow_mainland | c | chr17:13528059:G:A | 17 | 13,528,059 | 2.305 | 1.569 | 3.387 | 7.918e-06 | 2.535 | 1.717 | 3.745 | 7.560e-07 |
| narrow_mainland | c | chr17:76739260:T:C | 17 | 76,739,260 | 2.847 | 1.773 | 4.571 | 4.416e-06 | 3.1 | 1.921 | 5.005 | 7.439e-07 |
| narrow_mainland | c | chr19:10630195:G:C | 19 | 10,630,195 | 1.5 | 1.256 | 1.792 | 3.763e-06 | 1.49 | 1.251 | 1.776 | 3.940e-06 |
| narrow_mainland | c | chr20:60915853:C:T | 20 | 60,915,853 | 3.091 | 1.899 | 5.031 | 1.341e-06 | 3.05 | 1.877 | 4.957 | 1.466e-06 |
| intermediate_mainland | b | chr16:53882967:G:A | 16 | 53,882,967 | 1.679 | 1.386 | 2.033 | 3.836e-08 | 1.664 | 1.375 | 2.014 | 5.968e-08 |
| intermediate_mainland | c | chr2:232438125:T:C | 2 | 232,438,125 | 0.499 | 0.3683 | 0.676 | 2.386e-06 | 0.4978 | 0.3674 | 0.6745 | 2.181e-06 |
| intermediate_mainland | c | chr2:233068979:G:A | 2 | 233,068,979 | 4.102 | 2.147 | 7.834 | 3.133e-06 | 4.131 | 2.162 | 7.893 | 2.621e-06 |
| intermediate_mainland | c | chr3:94020684:T:C | 3 | 94,020,684 | 6.078 | 2.938 | 12.57 | 6.839e-08 | 4.818 | 2.413 | 9.619 | 9.732e-07 |
| intermediate_mainland | c | chr3:100783010:T:C | 3 | 100,783,010 | 1.857 | 1.421 | 2.426 | 2.095e-06 | 1.837 | 1.409 | 2.394 | 2.538e-06 |
| intermediate_mainland | c | chr3:154069965:A:G | 3 | 154,069,965 | 3.904 | 2.148 | 7.097 | 1.214e-06 | 3.853 | 2.125 | 6.988 | 1.442e-06 |
| intermediate_mainland | c | chr4:171998307:T:C | 4 | 171,998,307 | 4.484 | 2.369 | 8.488 | 4.801e-07 | 4.26 | 2.265 | 8.011 | 9.633e-07 |
| intermediate_mainland | c | chr6:32632861:A:G | 6 | 32,632,861 | 0.6638 | 0.561 | 0.7853 | 8.570e-07 | 0.6574 | 0.556 | 0.7773 | 4.397e-07 |
| intermediate_mainland | c | chr8:132794261:A:T | 8 | 132,794,261 | 3.78 | 2.082 | 6.862 | 1.915e-06 | 3.332 | 1.866 | 5.951 | 1.110e-05 |
| intermediate_mainland | c | chr9:104559024:C:T | 9 | 104,559,024 | 1.608 | 1.31 | 1.975 | 2.418e-06 | 1.613 | 1.314 | 1.981 | 2.114e-06 |
| intermediate_mainland | c | chr10:1747167:CGTG:C | 10 | 1,747,167 | 0.3522 | 0.2163 | 0.5737 | 5.274e-06 | 0.3561 | 0.218 | 0.5818 | 7.361e-06 |
| intermediate_mainland | c | chr14:88019710:G:A | 14 | 88,019,710 | 1.503 | 1.277 | 1.768 | 4.406e-07 | 1.504 | 1.279 | 1.769 | 3.832e-07 |
| intermediate_mainland | c | chr17:13528059:G:A | 17 | 13,528,059 | 2.548 | 1.735 | 3.74 | 4.175e-07 | 2.648 | 1.797 | 3.902 | 1.622e-07 |
| intermediate_mainland | c | chr17:76739260:T:C | 17 | 76,739,260 | 2.976 | 1.862 | 4.758 | 1.076e-06 | 3.008 | 1.88 | 4.813 | 8.443e-07 |
| intermediate_mainland | c | chr19:10630195:G:C | 19 | 10,630,195 | 1.492 | 1.257 | 1.771 | 2.312e-06 | 1.485 | 1.251 | 1.762 | 2.935e-06 |
| intermediate_mainland | c | chr20:60915853:C:T | 20 | 60,915,853 | 2.908 | 1.802 | 4.69 | 2.744e-06 | 2.811 | 1.748 | 4.52 | 5.009e-06 |
| full_mainland | b | chr16:53882967:G:A | 16 | 53,882,967 | 1.611 | 1.338 | 1.941 | 1.832e-07 | 1.576 | 1.311 | 1.895 | 5.141e-07 |
| full_mainland | c | chr2:232438125:T:C | 2 | 232,438,125 | 0.5113 | 0.3786 | 0.6906 | 3.914e-06 | 0.5228 | 0.3887 | 0.7031 | 5.877e-06 |
| full_mainland | c | chr2:233068979:G:A | 2 | 233,068,979 | 4.311 | 2.261 | 8.221 | 1.017e-06 | 4.69 | 2.432 | 9.044 | 2.910e-07 |
| full_mainland | c | chr3:94020684:T:C | 3 | 94,020,684 | 6.034 | 2.911 | 12.51 | 6.422e-08 | 5 | 2.487 | 10.05 | 5.097e-07 |
| full_mainland | c | chr3:100783010:T:C | 3 | 100,783,010 | 1.828 | 1.407 | 2.375 | 2.198e-06 | 1.802 | 1.392 | 2.333 | 2.738e-06 |
| full_mainland | c | chr3:154069965:A:G | 3 | 154,069,965 | 3.888 | 2.17 | 6.964 | 5.734e-07 | 3.949 | 2.204 | 7.077 | 4.135e-07 |
| full_mainland | c | chr4:171998307:T:C | 4 | 171,998,307 | 4.166 | 2.238 | 7.756 | 7.583e-07 | 4.155 | 2.236 | 7.723 | 6.814e-07 |
| full_mainland | c | chr6:32632861:A:G | 6 | 32,632,861 | 0.6682 | 0.5669 | 0.7877 | 7.319e-07 | 0.6665 | 0.5666 | 0.7839 | 4.520e-07 |
| … | 8 further rows in the TSV |  |  |  |  |  |  |  |  |  |  |  |

## How to read it

1. Read panel (a) as pairs, not as a series: the two markers above one cohort are the comparison, and the distance between them is what the PCs are worth there.
2. A probe already near 1 means the GRM absorbed the structure by itself.
3. In panels (b) and (c), read distance from the red identity line. A cloud lying on it means the covariates are doing nothing for these loci.
4. A point far below the line in (b) but on it in (c) is a locus whose ESTIMATE is unchanged and whose EVIDENCE weakened — a standard-error effect, not a signal that depended on the covariates.
5. lambda is reported and never applied, here as everywhere else in this component.

## What this figure does *not* establish

- It cannot say which model is correct. It reports the difference between them; the reported model is chosen on design grounds, not on which lambda is smaller.
- The probe is not an independent analysis. It shares every sample and the whole GRM with the reported model, so agreement between them carries no replication weight.
- Only the peak model defines loci, so panels (b) and (c) can only show leads that the reported model already found. A locus the probe alone would have flagged does not appear.
- Suggestive leads are shown to give the comparison a population to be read over. They remain a description of the scan, not findings.

## Symbols

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

## Model

```
lambda_GC = median(chi^2),/,chi^2_1,0.5,   chi^2_1,0.5 = 0.4549
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
