# HLA allele scan: fixed-effects against full-GRM model

**Figure file:** `model_compare.allele.png`

## The question this figure answers

Does putting a full genetic relationship matrix in the model change which HLA allele markers are associated, or by how much — and is either fit calibrated?

## Panels

**(a) Evidence, model against model**

-log10 P under fixed effects on x and full GRM on y, one point per allele marker, one sub-panel per cohort. All sub-panels share one limit, so a marker can be compared between cohorts as well as between models. The grey line is equality; the dashed red lines mark the significance threshold on each axis, so whether a marker crosses under one model and not the other is readable directly. Points clearing the threshold under either model are drawn in red at a larger size.

**(b) Effect concordance**

log odds ratio under the two models, pooled over cohorts and coloured by cohort, with 95% confidence whiskers drawn on the significant markers only — every marker would put a few hundred crossed intervals over a cloud whose shape is the whole message. Dotted lines mark log OR = 0 (no effect) on both axes. The annotated r is the Pearson correlation over every marker fitted under both models.

**(c) Calibration**

lambda_GC = median(chi2)/0.4549 for each cohort under each model, computed from the P-values of markers with a usable fit (ERRCODE = "."). The bars are drawn as deviations from 1 rather than from 0, because a bar from zero renders every plausible lambda at visually the same height. An HLA scan is a few hundred tests over one region under strong LD, so this lambda is a much noisier and more structure-dependent quantity than a genome-wide one; it is read as a coarse read-out of whether one model is systematically shifted relative to the other, not as an inflation estimate.

**(d) The hits under both models**

Every marker that clears the threshold under EITHER model, in every cohort, with its odds ratio and 95% CI under BOTH. One row group per marker; within a group one row per cohort, and within a row the two models are offset from each other, so the pair that must be compared is adjacent by construction. Colour is the cohort; the filled marker is fixed effects. The right-hand column carries both P-values in the order fixed effects / full GRM.

## Interpretation

The two models fit the SAME samples with the SAME phenotype and covariates and differ only in whether a full genetic relationship matrix is included as a random effect. The difference between them is therefore what relatedness and fine-scale structure were contributing to the fixed-effects result, and nothing else. The HLA region is where that matters most in this study: it is the most polymorphic region of the genome, its allele frequencies differ between Japanese subpopulations, and the cases and the controls of this study were sequenced on different platforms. A marker sitting on the identity line in panels (a) and (b) is one the GRM does not touch. A marker that loses evidence under the GRM was, in part, being carried by relatedness; panel (d) separates the two ways that can happen, because an odds ratio that holds while its interval widens is a standard-error effect, while an odds ratio that moves toward 1 is the estimate itself deflating. Neither model is declared correct here and lambda_GC is never applied as a correction — it is a calibration read-out, as everywhere else in this study.

## Values in this rendering

| quantity | value |
|---|---|
| marker class | allele |
| cohorts | narrow_mainland, intermediate_mainland, full_mainland |
| models | fixed (fixed effects), random (full GRM) |
| markers with a usable fit | 170 |
| fits excluded on ERRCODE | 0 |
| significance threshold | 2.941e-04 |
| threshold source | Bonferroni 0.05 / 170 |
| markers significant under either model | 2 |
| markers drawn in the forest | 2 |
| forest rows | 6 |
| pooled log-OR correlation | 0.9772 |
| lambda_GC narrow_mainland / fixed | 1.753 |
| lambda_GC narrow_mainland / random | 1.706 |
| lambda_GC intermediate_mainland / fixed | 1.901 |
| lambda_GC intermediate_mainland / random | 1.811 |
| lambda_GC full_mainland / fixed | 1.602 |
| lambda_GC full_mainland / random | 1.474 |
| markers significant, narrow_mainland / fixed | 2 |
| markers significant, narrow_mainland / random | 1 |
| markers significant, intermediate_mainland / fixed | 1 |
| markers significant, intermediate_mainland / random | 1 |
| markers significant, full_mainland / fixed | 1 |
| markers significant, full_mainland / random | 1 |

## Full statistics

**Every significant marker, in every cohort, under both models**

| marker | gene | position | cohort | model | OR | L95 | U95 | log_or | SE | P | A1_FREQ | OBS_CT | significant |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F |  | narrow_mainland | fixed | 2.906 | 1.663 | 5.078 | 1.067 | 0.2847 | 1.787e-04 | 0.006458 | 2,168 | 1 |
| HLA_F*01 | F |  | narrow_mainland | random | 2.886 | 1.585 | 5.253 | 1.06 | 0.3056 | 2.625e-04 | 0.006378 | 2,195 | 1 |
| HLA_F*01 | F |  | intermediate_mainland | fixed | 3.069 | 1.774 | 5.308 | 1.121 | 0.2796 | 6.066e-05 | 0.005645 | 2,480 | 1 |
| HLA_F*01 | F |  | intermediate_mainland | random | 3.043 | 1.678 | 5.52 | 1.113 | 0.3038 | 1.247e-04 | 0.005582 | 2,508 | 1 |
| HLA_F*01 | F |  | full_mainland | fixed | 3.002 | 1.782 | 5.055 | 1.099 | 0.266 | 3.583e-05 | 0.004884 | 3,071 | 1 |
| HLA_F*01 | F |  | full_mainland | random | 2.995 | 1.651 | 5.434 | 1.097 | 0.3039 | 1.536e-04 | 0.004837 | 3,101 | 1 |
| HLA_B*40:02 | B |  | narrow_mainland | fixed | 1.701 | 1.279 | 2.263 | 0.5313 | 0.1457 | 2.654e-04 | 0.05835 | 2,168 | 1 |
| HLA_B*40:02 | B |  | narrow_mainland | random | 1.721 | 1.255 | 2.359 | 0.5428 | 0.1609 | 3.705e-04 | 0.05809 | 2,195 | 0 |
| HLA_B*40:02 | B |  | intermediate_mainland | fixed | 1.55 | 1.187 | 2.024 | 0.4382 | 0.1363 | 0.0013 | 0.06149 | 2,480 | 0 |
| HLA_B*40:02 | B |  | intermediate_mainland | random | 1.557 | 1.153 | 2.103 | 0.443 | 0.1533 | 0.001928 | 0.0614 | 2,508 | 0 |
| HLA_B*40:02 | B |  | full_mainland | fixed | 1.446 | 1.117 | 1.871 | 0.3685 | 0.1315 | 0.005069 | 0.06529 | 3,071 | 0 |
| HLA_B*40:02 | B |  | full_mainland | random | 1.456 | 1.071 | 1.979 | 0.3755 | 0.1568 | 0.008301 | 0.06514 | 3,101 | 0 |

## How to read it

1. Read panel (a) as one cloud per cohort against the identity line, not as three separate scans: the limits are shared.
2. A red point off the identity line in (a) is a marker whose EVIDENCE depends on the model. Go to (d) and read its two intervals to see whether the estimate moved or only the standard error did.
3. In (c), read the pair above each cohort, not the series: the distance between the two bars is what the GRM is worth there.
4. In (d), read across a row for the model comparison and down a group for the cohort comparison. The cohorts are nested, so the three rows of a group share samples.
5. Nothing here ranks the two models. Which one is reported is a design decision, stated in METHODS, not the one with the smaller lambda.

## What this figure does *not* establish

- It is not replication. The two models fit the same samples, so agreement between them says the result does not depend on the random effect — not that it is real.
- lambda_GC over a few hundred markers in one strongly-LD region is a coarse statistic. It cannot separate confounding from a genuine concentration of signal, and in the HLA the latter is expected.
- Only the allele class is drawn. The other marker class is a separate analysis with its own multiple-testing burden and its own figure.
- A marker absent from one model's output does not appear in that model's panels. Rows where one model produced no fit are annotated rather than dropped in (d).
- The threshold is a Bonferroni over the markers of this class and treats them as independent tests, which HLA markers are not — they are in strong LD, so this is conservative in the number of tests and says nothing about the correlation between them.

## Symbols

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

- **ERRCODE** — plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, SEPARATION) did not fit cleanly and are excluded from lambda_GC and from the hit list rather than silently carried.

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

## Model

```
lambda_GC = median(chi^2),/,chi^2_1,0.5,   chi^2_1,0.5 = 0.4549
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
