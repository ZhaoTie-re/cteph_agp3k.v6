# Robust genes — tiers, effect sizes and carriers

**Figure file:** `robust_genes.png`

## The question this figure answers

Which genes survive a change of statistic and the addition of controls, how strong is their burden effect, and how many people carry it?

## Panels

**(a) Evidence in every cell, and the tier per variant set**

Rows are gene x cohort: every gene reaching Tier 1-3 is repeated once per cohort, so that (b) and (c) can give each cohort its own interval and its own denominator. Within each variant set (left LOW+MODERATE+HIGH, right MODERATE+HIGH) the two cells are the two statistics. Four states: a full red circle is above that cell's own Bonferroni threshold (and therefore also BH-called); a lighter red circle is BH-called only; a small blue dot is tested and called by neither; a grey cross is NOT tested -- the gene has fewer than MinNumVar variants in that cohort's map. The key shows only the states this run actually draws, so a mark a reader cannot find is never advertised; every cell here was tested, so no cross appears and the key has three states. The boxed number after each block is the tier that variant set awards; it is a property of the gene, not of a cohort, so it is written once per block: 1 = called in every cohort and by both statistics in at least one; 2 = called in every cohort; 3 = called in at least two cohorts; 0 = tested but called in fewer than two; -- = not tested in any cohort. The row band and the bar in the right gutter carry the gene's BEST tier over the variant sets, so a gene can sit in a Tier 1 band and carry a 0 in the set that did not call it. The tier channel is grey throughout: the red and blue belong to the cell states.

**(b) CMC burden odds ratio**

exp(b) with exp(b +/- 1.96 SE) from rvtest --burden cmcWald, log axis, dashed line at OR = 1. One cohort per row, so every row carries its OWN full 95 % interval; colour follows the cohort on the ordered ramp, darkest narrow to lightest full, which restates the row label rather than adding a variable. Two markers per row, one per variant set (filled MODERATE+HIGH, hollow LOW+MODERATE+HIGH). A cell whose Wald fit rvtest reset -- no minor allele in one group, which sends the interval to 1e-56 .. 1e66 -- has no usable estimate. It is marked n.e. in the margin, excluded from the axis limits, and reads "not estimable" in the table. n.e. is deliberately NOT the × of (a): × means the gene is not in that cohort's map and was never tested, whereas an n.e. cell WAS tested and keeps its state in (a). P2RY2 is the case here -- no control carries a MODERATE+HIGH allele in narrow or intermediate, which (c) shows as a carrier count of 0, yet the gene is BH-called in intermediate. Both facts are true, and the two marks now let the figure say both without contradicting itself.

**(c) Carrier frequency**

The fraction of controls (blue) and of cases (red) carrying at least one qualifying minor allele of the gene, for MODERATE+HIGH where the gene is mapped there and the wider set otherwise (hollow markers). Each row is scored against ITS OWN cohort's N, and the counts at the right are that cohort's. Carrier counts need not grow with the cohort: each cohort sets its own minAC, so a gene's mapped variant set is decided per cohort and a variant that passes in a smaller cohort can be dropped from a larger one. PES1 is the clearest case here -- 4 MODERATE+HIGH variants in narrow and intermediate, 3 in full -- and the n_var_map column of the tables below carries the number for every row. A carrier holding several sites counts once here and once per site in the cumulative MAC of the summary table.

## Interpretation

11 gene(s) reach a tier. Tier 1 is the strongest reading this design offers; Tier 2 drops the both-statistics condition; Tier 3 asks only for two cohorts. Two caveats outrank the tiers. The cohorts are nested and differ almost entirely in controls, so agreement across them shows the signal survives adding controls — a sample-selection sensitivity analysis, not replication. And the variant sets are nested, so a gene tiered in both is one finding tested twice on overlapping sets. What the figure genuinely separates is the statistic: a gene called by SKAT-O and not by CMC (JAK2-like rows, with an OR near 1) carries a signal in a subset of its variants or in opposite directions, which a collapse cannot see. Cases and controls were sequenced on different platforms, so no gene here is established as a CTEPH gene.

## Values in this rendering

| quantity | value |
|---|---|
| genes reaching a tier | 11 |
| genes drawn | 11 |
| Tier 1 genes | 2 |
| Tier 2 genes | 3 |
| Tier 3 genes | 6 |
| report cohort | full_mainland |
| N case / N control | 438 / 2,630 |
| cells not in a cohort's map | 0 |
| cells with no OR estimate | 2 of 66 |
| omitted for space | — |

## Full statistics

**Tiers (from tiers.tsv)**

| gene_symbol | stratum | tier | n_cohorts_called | cohorts_called | n_cohorts_both_methods | called_bonferroni_anywhere | cohorts_not_in_map | min_p_cmc | min_p_skato |
|---|---|---|---|---|---|---|---|---|---|
| STBD1 | moderate_high | 1 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 2 | 1 | — | 7.22605e-08 | 1.41108e-07 |
| OSBP | moderate_high | 1 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 1 | 1 | — | 3.54112e-07 | 3.21117e-07 |
| STBD1 | low_moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 4.10249e-06 | 8.88365e-07 |
| PIK3C2G | low_moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 0.00047545 | 3.48259e-07 |
| RPH3A | low_moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 0.000479202 | 5.44525e-07 |
| JAK2 | moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 0.17174 | 1.31438e-06 |
| PIK3C2G | moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 0.000203921 | 4.1323e-07 |
| RPH3A | moderate_high | 2 | 3 | narrow_mainland,intermediate_mainland,full_mainland | 0 | 1 | — | 0.000110518 | 3.908e-07 |
| SCAMP3 | low_moderate_high | 3 | 2 | intermediate_mainland,full_mainland | 1 | 1 | — | 2.75881e-06 | 4.15802e-06 |
| CD248 | low_moderate_high | 3 | 2 | narrow_mainland,intermediate_mainland | 0 | 1 | — | 0.0100787 | 2.02765e-06 |
| SCFD1 | low_moderate_high | 3 | 2 | narrow_mainland,full_mainland | 0 | 0 | — | 0.0381582 | 1.02818e-05 |
| PES1 | low_moderate_high | 3 | 2 | narrow_mainland,intermediate_mainland | 0 | 0 | — | 0.000147683 | 1.09234e-05 |
| NKD2 | moderate_high | 3 | 2 | narrow_mainland,intermediate_mainland | 0 | 0 | — | 5.91929e-05 | 6.18644e-06 |
| P2RY2 | moderate_high | 3 | 2 | intermediate_mainland,full_mainland | 1 | 0 | — | 1.30515e-05 | 1.60768e-05 |
| SCFD1 | moderate_high | 3 | 2 | narrow_mainland,full_mainland | 0 | 0 | — | 0.0183765 | 1.83712e-05 |
| PES1 | moderate_high | 3 | 2 | narrow_mainland,intermediate_mainland | 0 | 0 | — | 4.71219e-05 | 8.65008e-06 |

**Odds ratios drawn in (b)**

| gene_symbol | stratum | cohort | or | or_l95 | or_u95 | cmc_state | cmc_p | n_var_map | skato_p |
|---|---|---|---|---|---|---|---|---|---|
| STBD1 | low_moderate_high | narrow_mainland | 6.946 | 2.729 | 17.68 | not_significant | 4.10249e-06 | 5 | 6.72217e-06 |
| STBD1 | moderate_high | narrow_mainland | 9.978 | 2.931 | 33.96 | not_significant | 1.31634e-05 | 3 | 2.13224e-05 |
| STBD1 | low_moderate_high | intermediate_mainland | 4.767 | 2.178 | 10.44 | not_significant | 2.03095e-05 | 6 | 1.32989e-05 |
| STBD1 | moderate_high | intermediate_mainland | 7.349 | 2.705 | 19.96 | called_bh_only | 7.1974e-06 | 3 | 1.20885e-05 |
| STBD1 | low_moderate_high | full_mainland | 4.059 | 1.984 | 8.304 | not_significant | 3.66104e-05 | 8 | 8.88365e-07 |
| STBD1 | moderate_high | full_mainland | 9.601 | 3.565 | 25.86 | called_bonferroni | 7.22605e-08 | 3 | 1.41108e-07 |
| OSBP | low_moderate_high | narrow_mainland | 4.054 | 1.692 | 9.715 | not_significant | 0.000782652 | 6 | 0.000216417 |
| OSBP | moderate_high | narrow_mainland | 16.23 | 3.374 | 78.1 | not_significant | 4.70736e-06 | 2 | 5.14437e-06 |
| OSBP | low_moderate_high | intermediate_mainland | 4.245 | 1.887 | 9.548 | not_significant | 0.000163422 | 6 | 6.70093e-05 |
| OSBP | moderate_high | intermediate_mainland | 19.39 | 4.047 | 92.87 | called_bonferroni | 3.54112e-07 | 2 | 3.21117e-07 |
| OSBP | low_moderate_high | full_mainland | 2.767 | 1.349 | 5.672 | not_significant | 0.0039365 | 6 | 0.00323611 |
| OSBP | moderate_high | full_mainland | 7.105 | 2.53 | 19.95 | not_significant | 1.93294e-05 | 2 | 2.08318e-05 |
| PIK3C2G | low_moderate_high | narrow_mainland | 1.968 | 1.078 | 3.591 | not_significant | 0.0251171 | 15 | 8.45814e-07 |
| PIK3C2G | moderate_high | narrow_mainland | 2.648 | 1.361 | 5.153 | not_significant | 0.00305426 | 10 | 4.1323e-07 |
| PIK3C2G | low_moderate_high | intermediate_mainland | 2.377 | 1.37 | 4.124 | not_significant | 0.00155852 | 17 | 3.48259e-07 |
| PIK3C2G | moderate_high | intermediate_mainland | 2.674 | 1.437 | 4.979 | not_significant | 0.00130649 | 11 | 5.64051e-07 |
| PIK3C2G | low_moderate_high | full_mainland | 2.252 | 1.413 | 3.591 | not_significant | 0.00047545 | 21 | 5.7749e-06 |
| PIK3C2G | moderate_high | full_mainland | 2.566 | 1.535 | 4.289 | not_significant | 0.000203921 | 12 | 1.11519e-05 |
| RPH3A | low_moderate_high | narrow_mainland | 2.199 | 1.307 | 3.701 | not_significant | 0.00244435 | 11 | 1.80208e-05 |
| RPH3A | moderate_high | narrow_mainland | 2.615 | 1.475 | 4.637 | not_significant | 0.000693429 | 7 | 1.97181e-05 |
| RPH3A | low_moderate_high | intermediate_mainland | 2.314 | 1.405 | 3.809 | not_significant | 0.000723272 | 11 | 1.13239e-06 |
| RPH3A | moderate_high | intermediate_mainland | 2.826 | 1.634 | 4.888 | not_significant | 0.000113055 | 7 | 7.95035e-07 |
| RPH3A | low_moderate_high | full_mainland | 2.314 | 1.428 | 3.751 | not_significant | 0.000479202 | 11 | 5.44525e-07 |
| RPH3A | moderate_high | full_mainland | 2.715 | 1.606 | 4.589 | not_significant | 0.000110518 | 7 | 3.908e-07 |
| JAK2 | low_moderate_high | narrow_mainland | 1.175 | 0.685 | 2.017 | not_significant | 0.557152 | 9 | 0.00086482 |
| JAK2 | moderate_high | narrow_mainland | 1.707 | 0.7868 | 3.702 | not_significant | 0.17174 | 6 | 5.98467e-06 |
| JAK2 | low_moderate_high | intermediate_mainland | 1.045 | 0.6221 | 1.756 | not_significant | 0.867549 | 10 | 0.000543947 |
| JAK2 | moderate_high | intermediate_mainland | 1.425 | 0.6864 | 2.957 | not_significant | 0.339901 | 6 | 5.60992e-06 |
| JAK2 | low_moderate_high | full_mainland | 1.048 | 0.6315 | 1.739 | not_significant | 0.856374 | 11 | 0.000481945 |
| JAK2 | moderate_high | full_mainland | 1.564 | 0.7657 | 3.194 | not_significant | 0.21622 | 7 | 1.31438e-06 |
| SCAMP3 | low_moderate_high | narrow_mainland | 2.534 | 1.572 | 4.085 | not_significant | 8.66446e-05 | 8 | 0.000139187 |
| SCAMP3 | moderate_high | narrow_mainland | 2.099 | 1.094 | 4.027 | not_significant | 0.0230561 | 5 | 0.0536292 |
| SCAMP3 | low_moderate_high | intermediate_mainland | 2.9 | 1.826 | 4.606 | called_bonferroni | 2.75881e-06 | 8 | 4.15802e-06 |
| SCAMP3 | moderate_high | intermediate_mainland | 2.363 | 1.262 | 4.428 | not_significant | 0.0058082 | 5 | 0.0144526 |
| SCAMP3 | low_moderate_high | full_mainland | 2.522 | 1.633 | 3.896 | not_significant | 1.71334e-05 | 10 | 2.05815e-05 |
| SCAMP3 | moderate_high | full_mainland | 2.234 | 1.223 | 4.081 | not_significant | 0.00737174 | 5 | 0.018589 |
| CD248 | low_moderate_high | narrow_mainland | 2.512 | 1.152 | 5.477 | not_significant | 0.0172045 | 6 | 2.02765e-06 |
| CD248 | moderate_high | narrow_mainland | 0.6235 | 0.07425 | 5.235 | not_significant | 0.660605 | 2 | 0.771239 |
| CD248 | low_moderate_high | intermediate_mainland | 2.597 | 1.225 | 5.507 | not_significant | 0.0100787 | 6 | 4.22616e-06 |
| CD248 | moderate_high | intermediate_mainland | 0.7135 | 0.08599 | 5.92 | not_significant | 0.753406 | 2 | 0.848813 |
| CD248 | low_moderate_high | full_mainland | 2.151 | 1.057 | 4.379 | not_significant | 0.0307607 | 6 | 0.000300587 |
| CD248 | moderate_high | full_mainland | 0.5714 | 0.07202 | 4.534 | not_significant | 0.591708 | 2 | 0.719939 |
| SCFD1 | low_moderate_high | narrow_mainland | 2.273 | 0.9913 | 5.213 | not_significant | 0.0471149 | 4 | 1.94136e-05 |
| SCFD1 | moderate_high | narrow_mainland | 2.679 | 1.149 | 6.248 | not_significant | 0.0183765 | 3 | 1.83712e-05 |
| SCFD1 | low_moderate_high | intermediate_mainland | 1.869 | 0.8853 | 3.946 | not_significant | 0.0960194 | 5 | 0.000170785 |
| SCFD1 | moderate_high | intermediate_mainland | 2.164 | 1.012 | 4.627 | not_significant | 0.0417087 | 4 | 0.00016009 |
| SCFD1 | low_moderate_high | full_mainland | 1.947 | 1.026 | 3.696 | not_significant | 0.0381582 | 6 | 1.02818e-05 |
| SCFD1 | moderate_high | full_mainland | 1.991 | 0.9624 | 4.119 | not_significant | 0.0586716 | 4 | 2.2038e-05 |
| PES1 | low_moderate_high | narrow_mainland | 2.236 | 1.454 | 3.437 | not_significant | 0.000183236 | 8 | 1.09234e-05 |
| PES1 | moderate_high | narrow_mainland | 2.789 | 1.67 | 4.656 | not_significant | 5.01841e-05 | 4 | 8.65008e-06 |
| PES1 | low_moderate_high | intermediate_mainland | 2.192 | 1.449 | 3.317 | not_significant | 0.000147683 | 8 | 1.4951e-05 |
| PES1 | moderate_high | intermediate_mainland | 2.671 | 1.637 | 4.356 | not_significant | 4.71219e-05 | 4 | 1.07988e-05 |
| PES1 | low_moderate_high | full_mainland | 1.079 | 0.7076 | 1.645 | not_significant | 0.723959 | 9 | 0.167861 |
| PES1 | moderate_high | full_mainland | 1.189 | 0.3963 | 3.569 | not_significant | 0.756855 | 3 | 0.670125 |
| NKD2 | low_moderate_high | narrow_mainland | 5.007 | 2.121 | 11.82 | not_significant | 5.91929e-05 | 3 | 6.18644e-06 |
| NKD2 | moderate_high | narrow_mainland | 5.007 | 2.121 | 11.82 | not_significant | 5.91929e-05 | 3 | 6.18644e-06 |
| NKD2 | low_moderate_high | intermediate_mainland | 3.896 | 1.747 | 8.687 | not_significant | 0.000385884 | 3 | 2.79127e-05 |
| NKD2 | moderate_high | intermediate_mainland | 3.896 | 1.747 | 8.687 | not_significant | 0.000385884 | 3 | 2.79127e-05 |
| NKD2 | low_moderate_high | full_mainland | 3.426 | 1.62 | 7.243 | not_significant | 0.000657816 | 3 | 5.40276e-05 |
| NKD2 | moderate_high | full_mainland | 3.426 | 1.62 | 7.243 | not_significant | 0.000657816 | 3 | 5.40276e-05 |
| P2RY2 | low_moderate_high | narrow_mainland | 2.841 | 1.238 | 6.522 | not_significant | 0.0108683 | 3 | 0.00959046 |
| P2RY2 | moderate_high | narrow_mainland | not estimable | NA | NA | not_significant | 0.000102812 | 2 | 0.000140924 |
| P2RY2 | low_moderate_high | intermediate_mainland | 3.117 | 1.413 | 6.876 | not_significant | 0.00319724 | 3 | 0.00242745 |
| P2RY2 | moderate_high | intermediate_mainland | not estimable | NA | NA | called_bh_only | 1.30515e-05 | 2 | 1.60768e-05 |
| P2RY2 | low_moderate_high | full_mainland | 3.275 | 1.582 | 6.778 | not_significant | 0.000767725 | 5 | 0.000634595 |
| P2RY2 | moderate_high | full_mainland | 9.416 | 2.541 | 34.9 | not_significant | 6.49379e-05 | 4 | 2.08325e-05 |

**Carriers drawn in (c)**

| gene_symbol | stratum | cohort | carriers_case | n_case | pct_case | carriers_control | n_control | pct_control | n_var_map |
|---|---|---|---|---|---|---|---|---|---|
| STBD1 | moderate_high | narrow_mainland | 9 | 418 | 2.153 | 4 | 1,747 | 0.229 | 3 |
| STBD1 | moderate_high | intermediate_mainland | 10 | 428 | 2.336 | 7 | 2,049 | 0.342 | 3 |
| STBD1 | moderate_high | full_mainland | 10 | 438 | 2.283 | 7 | 2,630 | 0.266 | 3 |
| OSBP | moderate_high | narrow_mainland | 8 | 418 | 1.914 | 2 | 1,747 | 0.114 | 2 |
| OSBP | moderate_high | intermediate_mainland | 8 | 428 | 1.869 | 2 | 2,049 | 0.098 | 2 |
| OSBP | moderate_high | full_mainland | 8 | 438 | 1.826 | 8 | 2,630 | 0.304 | 2 |
| PIK3C2G | moderate_high | narrow_mainland | 15 | 418 | 3.589 | 26 | 1,747 | 1.488 | 10 |
| PIK3C2G | moderate_high | intermediate_mainland | 16 | 428 | 3.738 | 33 | 2,049 | 1.611 | 11 |
| PIK3C2G | moderate_high | full_mainland | 22 | 438 | 5.023 | 60 | 2,630 | 2.281 | 12 |
| RPH3A | moderate_high | narrow_mainland | 21 | 418 | 5.024 | 35 | 1,747 | 2.003 | 7 |
| RPH3A | moderate_high | intermediate_mainland | 22 | 428 | 5.14 | 39 | 2,049 | 1.903 | 7 |
| RPH3A | moderate_high | full_mainland | 22 | 438 | 5.023 | 51 | 2,630 | 1.939 | 7 |
| JAK2 | moderate_high | narrow_mainland | 10 | 418 | 2.392 | 23 | 1,747 | 1.317 | 6 |
| JAK2 | moderate_high | intermediate_mainland | 10 | 428 | 2.336 | 32 | 2,049 | 1.562 | 6 |
| JAK2 | moderate_high | full_mainland | 10 | 438 | 2.283 | 39 | 2,630 | 1.483 | 7 |
| SCAMP3 | moderate_high | narrow_mainland | 15 | 418 | 3.589 | 29 | 1,747 | 1.66 | 5 |
| SCAMP3 | moderate_high | intermediate_mainland | 16 | 428 | 3.738 | 30 | 2,049 | 1.464 | 5 |
| SCAMP3 | moderate_high | full_mainland | 16 | 438 | 3.653 | 38 | 2,630 | 1.445 | 5 |
| CD248 | moderate_high | narrow_mainland | 1 | 418 | 0.239 | 7 | 1,747 | 0.401 | 2 |
| CD248 | moderate_high | intermediate_mainland | 1 | 428 | 0.234 | 7 | 2,049 | 0.342 | 2 |
| CD248 | moderate_high | full_mainland | 1 | 438 | 0.228 | 10 | 2,630 | 0.38 | 2 |
| SCFD1 | moderate_high | narrow_mainland | 9 | 418 | 2.153 | 17 | 1,747 | 0.973 | 3 |
| SCFD1 | moderate_high | intermediate_mainland | 10 | 428 | 2.336 | 25 | 2,049 | 1.22 | 4 |
| SCFD1 | moderate_high | full_mainland | 10 | 438 | 2.283 | 36 | 2,630 | 1.369 | 4 |
| PES1 | moderate_high | narrow_mainland | 27 | 418 | 6.459 | 44 | 1,747 | 2.519 | 4 |
| PES1 | moderate_high | intermediate_mainland | 27 | 428 | 6.308 | 52 | 2,049 | 2.538 | 4 |
| PES1 | moderate_high | full_mainland | 4 | 438 | 0.913 | 19 | 2,630 | 0.722 | 3 |
| NKD2 | moderate_high | narrow_mainland | 11 | 418 | 2.632 | 12 | 1,747 | 0.687 | 3 |
| NKD2 | moderate_high | intermediate_mainland | 11 | 428 | 2.57 | 16 | 2,049 | 0.781 | 3 |
| NKD2 | moderate_high | full_mainland | 11 | 438 | 2.511 | 24 | 2,630 | 0.913 | 3 |
| P2RY2 | moderate_high | narrow_mainland | 5 | 418 | 1.196 | 0 | 1,747 | 0 | 2 |
| P2RY2 | moderate_high | intermediate_mainland | 5 | 428 | 1.168 | 0 | 2,049 | 0 | 2 |
| P2RY2 | moderate_high | full_mainland | 6 | 438 | 1.37 | 4 | 2,630 | 0.152 | 4 |

**Denominator per cohort × variant set**

| cohort | stratum | n_genes_mapped | bonferroni_threshold |
|---|---|---|---|
| intermediate_mainland | low_moderate_high | 15009 | 3.33133453261e-06 |
| intermediate_mainland | moderate_high | 11436 | 4.37215809724e-06 |
| full_mainland | low_moderate_high | 15691 | 3.18654005481e-06 |
| full_mainland | moderate_high | 12378 | 4.03942478591e-06 |
| narrow_mainland | low_moderate_high | 14557 | 3.43477364842e-06 |
| narrow_mainland | moderate_high | 10832 | 4.61595273264e-06 |

## How to read it

1. Read a row across: the tier numbers summarise (a); (b) and (c) are what the summary table prints for that gene, drawn.
2. No cell of (a) is a cross: every gene drawn here is in every cohort's map, which is expected -- a gene absent from a cohort's map cannot be called there and so rarely reaches a tier at all.
3. Compare the CMC and SKAT-O cells within one cohort block; that pair holds everything but the statistic fixed.
4. A wide interval in (b) with few carriers in (c) is a statement about a handful of people.

## What this figure does *not* establish

- Not replication: the cohorts are nested and share nearly all cases.
- Not a CTEPH gene: the platform confound is in every cell.
- The OR is the collapsed-carrier burden; it says nothing about any single variant.

## Model

```
called = BH q < alpha over the cell's mapped genes; tiers from robust_genes.tier_of; OR = exp(beta) from rvtest --burden cmcWald
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
