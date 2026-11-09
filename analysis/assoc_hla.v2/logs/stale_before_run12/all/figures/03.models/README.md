# CTEPH (all cases) vs controls: fixed against random model, marker by marker

**Covers:** `models.png` in `03.models/`

## The question this figure answers

Does accounting for relatedness and the genome-wide polygenic background (random model) change any marker's result compared with the fixed model?

## Panels

**(models.png a–c) −log10 P, fixed (x) against random (y), one panel per cohort**

Each point is one marker that both models returned; green = alleles, blue = residues. Grey line: identity. Dotted lines: each class's LD-adjusted (GEC) threshold.

**(models.png d–f) β, fixed (x) against random (y)**

Fixed: log OR from plink2. Random: SAIGE's one-step score estimate Tstat / var.

## Interpretation

r(−log10 P) ranges 0.995–0.999 and r(β) 0.969–0.995. The models differ in samples (the fixed model drops relatives), in the treatment of missing calls (complete case against mean imputation) and in the estimator, so points scatter around the identity; what matters is whether any marker crosses a threshold in one model and not the other (model_concordance.tsv).

## Values in this rendering

| quantity | value |
|---|---|
| models: allele narrow: r(-log10 P) / r(beta) | 0.997 / 0.976 |
| models: allele intermediate: r(-log10 P) / r(beta) | 0.996 / 0.972 |
| models: allele full: r(-log10 P) / r(beta) | 0.995 / 0.969 |
| models: residue narrow: r(-log10 P) / r(beta) | 0.999 / 0.995 |
| models: residue intermediate: r(-log10 P) / r(beta) | 0.999 / 0.995 |
| models: residue full: r(-log10 P) / r(beta) | 0.999 / 0.995 |

## Full statistics

**Fixed against random, per cohort × class**

| class | cohort | markers_both | pearson_neglog10p | spearman_p | pearson_beta | median_ratio_neglog10p | gec_fixed | gec_random | gec_both |
|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | 143 | 0.9972 | 0.9947 | 0.9759 | 0.9571 | 1 | 0 | 0 |
| allele | intermediate | 144 | 0.996 | 0.9957 | 0.972 | 0.9597 | 0 | 0 | 0 |
| allele | full | 142 | 0.9951 | 0.9939 | 0.9694 | 0.9377 | 0 | 0 | 0 |
| residue | narrow | 672 | 0.9989 | 0.9972 | 0.9951 | 0.9367 | 1 | 1 | 1 |
| residue | intermediate | 681 | 0.999 | 0.9974 | 0.9954 | 0.9422 | 1 | 1 | 1 |
| residue | full | 679 | 0.999 | 0.9975 | 0.9953 | 0.91 | 3 | 1 | 1 |

## How to read it

1. Points on the identity line: the two models agree.
2. Points between the dotted lines: called by one model only.

## What this figure does *not* establish

- Which model is right for a given marker: they answer the same question under different assumptions; a marker significant under only one is reported as such.
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
