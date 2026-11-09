# CTEPH with DVT vs CTEPH without DVT (cases only): fixed against random model, marker by marker

**Covers:** `models.png` in `03.models/`

## The question this figure answers

Does accounting for relatedness and the genome-wide polygenic background (random model) change any marker's result compared with the fixed model?

## Panels

**(models.png a–c) −log10 P, fixed (x) against random (y), one panel per cohort**

Each point is one marker that both models returned; green = alleles, blue = residues. Grey line: identity. Dotted lines: each class's LD-adjusted (GEC) threshold.

**(models.png d–f) β, fixed (x) against random (y)**

Fixed: log OR from plink2. Random: SAIGE's one-step score estimate Tstat / var.

## Interpretation

r(−log10 P) ranges 0.993–0.999 and r(β) 0.993–0.999. The models differ in samples (the fixed model drops relatives), in the treatment of missing calls (complete case against mean imputation) and in the estimator, so points scatter around the identity; what matters is whether any marker crosses a threshold in one model and not the other (model_concordance.tsv).

## Values in this rendering

| quantity | value |
|---|---|
| models: allele narrow: r(-log10 P) / r(beta) | 0.997 / 0.998 |
| models: allele intermediate: r(-log10 P) / r(beta) | 0.997 / 0.997 |
| models: allele full: r(-log10 P) / r(beta) | 0.993 / 0.993 |
| models: residue narrow: r(-log10 P) / r(beta) | 0.999 / 0.999 |
| models: residue intermediate: r(-log10 P) / r(beta) | 0.999 / 0.999 |
| models: residue full: r(-log10 P) / r(beta) | 0.999 / 0.996 |

## Full statistics

**Fixed against random, per cohort × class**

| class | cohort | markers_both | pearson_neglog10p | spearman_p | pearson_beta | median_ratio_neglog10p | gec_fixed | gec_random | gec_both |
|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | 132 | 0.9969 | 0.9982 | 0.9975 | 0.9848 | 0 | 0 | 0 |
| allele | intermediate | 133 | 0.997 | 0.998 | 0.9974 | 0.9885 | 0 | 0 | 0 |
| allele | full | 135 | 0.9932 | 0.9984 | 0.9925 | 0.9891 | 0 | 0 | 0 |
| residue | narrow | 663 | 0.9994 | 0.9989 | 0.9985 | 0.9812 | 0 | 0 | 0 |
| residue | intermediate | 672 | 0.9994 | 0.9992 | 0.9985 | 0.9817 | 0 | 0 | 0 |
| residue | full | 672 | 0.9991 | 0.9989 | 0.9963 | 0.9819 | 0 | 0 | 0 |

## How to read it

1. Points on the identity line: the two models agree.
2. Points between the dotted lines: called by one model only.

## What this figure does *not* establish

- Which model is right for a given marker: they answer the same question under different assumptions; a marker significant under only one is reported as such.
- The three cohorts share almost all their cases and this set has no controls, so the three cohorts are nearly the same analysis; agreement across them is expected, not a replication.
- Both groups are cases, so the case–control platform split does not enter. But the sequencing platform also differs between cases with and without DVT (docs/STUDY_NOTES.md), so a platform-dependent typing difference would still look like an association; platform is not a covariate here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
