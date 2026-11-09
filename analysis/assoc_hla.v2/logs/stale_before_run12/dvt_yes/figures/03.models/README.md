# CTEPH with DVT vs controls: fixed against random model, marker by marker

**Covers:** `models.png` in `03.models/`

## The question this figure answers

Does accounting for relatedness and the genome-wide polygenic background (random model) change any marker's result compared with the fixed model?

## Panels

**(models.png a–c) −log10 P, fixed (x) against random (y), one panel per cohort**

Each point is one marker that both models returned; green = alleles, blue = residues. Grey line: identity. Dotted lines: each class's LD-adjusted (GEC) threshold.

**(models.png d–f) β, fixed (x) against random (y)**

Fixed: log OR from plink2. Random: SAIGE's one-step score estimate Tstat / var.

## Interpretation

r(−log10 P) ranges 0.998–0.999 and r(β) 0.962–0.975. The models differ in samples (the fixed model drops relatives), in the treatment of missing calls (complete case against mean imputation) and in the estimator, so points scatter around the identity; what matters is whether any marker crosses a threshold in one model and not the other (model_concordance.tsv).

## Values in this rendering

| quantity | value |
|---|---|
| models: allele narrow: r(-log10 P) / r(beta) | 0.998 / 0.971 |
| models: allele intermediate: r(-log10 P) / r(beta) | 0.998 / 0.972 |
| models: allele full: r(-log10 P) / r(beta) | 0.998 / 0.962 |
| models: residue narrow: r(-log10 P) / r(beta) | 0.998 / 0.973 |
| models: residue intermediate: r(-log10 P) / r(beta) | 0.999 / 0.975 |
| models: residue full: r(-log10 P) / r(beta) | 0.999 / 0.967 |

## Full statistics

**Fixed against random, per cohort × class**

| class | cohort | markers_both | pearson_neglog10p | spearman_p | pearson_beta | median_ratio_neglog10p | gec_fixed | gec_random | gec_both |
|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | 134 | 0.9983 | 0.9957 | 0.9706 | 0.9908 | 0 | 0 | 0 |
| allele | intermediate | 135 | 0.9983 | 0.997 | 0.9717 | 0.9959 | 0 | 0 | 0 |
| allele | full | 137 | 0.9977 | 0.9953 | 0.9619 | 0.9926 | 0 | 0 | 0 |
| residue | narrow | 663 | 0.998 | 0.9979 | 0.973 | 0.9618 | 0 | 0 | 0 |
| residue | intermediate | 672 | 0.9985 | 0.9974 | 0.9747 | 0.9694 | 0 | 0 | 0 |
| residue | full | 672 | 0.9985 | 0.9977 | 0.9674 | 0.9606 | 0 | 0 | 0 |

## How to read it

1. Points on the identity line: the two models agree.
2. Points between the dotted lines: called by one model only.

## What this figure does *not* establish

- Which model is right for a given marker: they answer the same question under different assumptions; a marker significant under only one is reported as such.
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
