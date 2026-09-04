# Robust genes — full_mainland

3,068 samples (438 cases, 2,630 controls). 5 gene(s) at Tier 1 or Tier 2 in at least one stratum.

## LOW + MODERATE + HIGH

| Gene | Chr | Tier | Variants | Carriers case (%) | Carriers control (%) | MAC case | MAC control | MAF case | MAF control | SKAT-O rho | SKAT-O P | CMC OR (95% CI) | CMC P | Multi-site carriers case / control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| *STBD1* | 4 | 2 | 8 | 13 (2.97%) | 21 (0.80%) | 13 | 21 | 0.015 | 0.004 | 0.1 | 8.88E-07 | 4.06 (1.98-8.30) | 3.66E-05 | 0 / 0 |
| *OSBP* | 11 | 0 | 6 | 12 (2.74%) | 25 (0.95%) | 12 | 25 | 0.014 | 0.005 | 0.2 | 3.24E-03 | 2.77 (1.35-5.67) | 3.94E-03 | 0 / 0 |
| *JAK2* | 9 | 0 | 11 | 19 (4.34%) | 107 (4.07%) | 21 | 108 | 0.024 | 0.021 | 0 | 4.82E-04 | 1.05 (0.63-1.74) | 8.56E-01 | 2 / 1 |
| *PIK3C2G* | 12 | 2 | 21 | 26 (5.94%) | 80 (3.04%) | 26 | 82 | 0.030 | 0.016 | 0.1 | 5.77E-06 | 2.25 (1.41-3.59) | 4.75E-04 | 0 / 2 |
| *RPH3A* | 12 | 2 | 11 | 25 (5.71%) | 66 (2.51%) | 25 | 66 | 0.029 | 0.013 | 0 | 5.45E-07 | 2.31 (1.43-3.75) | 4.79E-04 | 0 / 0 |

## MODERATE + HIGH

| Gene | Chr | Tier | Variants | Carriers case (%) | Carriers control (%) | MAC case | MAC control | MAF case | MAF control | SKAT-O rho | SKAT-O P | CMC OR (95% CI) | CMC P | Multi-site carriers case / control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| *STBD1* | 4 | 1 | 3 | 10 (2.28%) | 7 (0.27%) | 10 | 7 | 0.011 | 0.001 | 0.3 | 1.41E-07 | 9.60 (3.56-25.86) | 7.23E-08 | 0 / 0 |
| *OSBP* | 11 | 1 | 2 | 8 (1.83%) | 8 (0.30%) | 8 | 8 | 0.009 | 0.002 | 1 | 2.08E-05 | 7.10 (2.53-19.95) | 1.93E-05 | 0 / 0 |
| *JAK2* | 9 | 2 | 7 | 10 (2.28%) | 39 (1.48%) | 11 | 39 | 0.013 | 0.007 | 0 | 1.31E-06 | 1.56 (0.77-3.19) | 2.16E-01 | 1 / 0 |
| *PIK3C2G* | 12 | 2 | 12 | 22 (5.02%) | 60 (2.28%) | 22 | 62 | 0.025 | 0.012 | 0.1 | 1.12E-05 | 2.57 (1.54-4.29) | 2.04E-04 | 0 / 2 |
| *RPH3A* | 12 | 2 | 7 | 22 (5.02%) | 51 (1.94%) | 22 | 51 | 0.025 | 0.010 | 0 | 3.91E-07 | 2.72 (1.61-4.59) | 1.11E-04 | 0 / 0 |

**Definitions.** Variants = variants of the stratum mapped to the gene in this cohort (>= MinNumVar). Carriers = samples with >= 1 minor allele at any of them (%: of the group). MAC = cumulative minor-allele count over the set. MAF = MAC / (2 x N of the group), a fixed denominator. SKAT-O rho = the mixing parameter rvtest selected; P = SKAT-O p. CMC OR (95% CI) = exp(beta) from rvtest --burden cmcWald with a Wald normal CI; CMC P = the CMC score-test p (the Wald p is not shown). Multi-site carriers = samples carrying a minor allele at >= 2 distinct sites of the gene (case / control); CMC counts each once, MAC counts every site. — = the gene is not in this cohort's map for that stratum, or the value is undefined.

**Tiers.** A gene is *called* in a cohort x stratum x method when Benjamini-Hochberg q < 0.05 over that cell's mapped genes. Tier 1 = called in every cohort and by both methods in at least one; Tier 2 = called in every cohort by at least one method; Tier 3 = called in at least two cohorts; 0 = tested in this stratum but called in fewer than two cohorts. The tier shown is the stratum's own; a gene is listed because it reaches Tier 1 or 2 in at least one stratum.

**Caveat.** The three cohorts are nested and differ almost entirely in controls; agreement across them shows the signal survives adding controls and is a sample-selection sensitivity analysis, not replication. Cases and controls were sequenced on different platforms, so no gene here is established as a CTEPH gene.
