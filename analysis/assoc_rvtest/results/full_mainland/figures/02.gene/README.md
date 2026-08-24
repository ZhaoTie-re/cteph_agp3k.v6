# Gene-based hits — full_mainland

**Figures:** `02.gene/<tier>/gene.<peak>.png` — 49 locus figure(s): 0 genome-wide, 31 suggestive.

One document covers the whole family: the panels, the reading and the limits are properties of the figure, identical for every locus. Only the numbers differ, and they are tabulated below, one row per locus.

## The question this figure answers

For each gene the burden test flagged, how many carriers actually produced the signal, and is the difference between cases and controls or between the cohort and the population?

## Panels

**(a) Case against control, per variant**

Alternate-allele frequency of every variant the gene-based test collapsed, in cases and in controls. A gene whose signal rests on ONE variant is a different claim from one where several variants agree, and the gene-level P cannot distinguish them.

**(b) Control against the population reference**

The same variants, control frequency against ToMMo 60KJPN. This is the artefact check that matters under this design: cases and controls were sequenced on different platforms at different depth, so a control frequency that departs from the population reference is a genotyping difference rather than a disease association. Agreement here is what makes the case excess in (a) worth reading.

**(c) The collapsed burden**

Cumulative minor-allele frequency over the gene in each group, labelled with the allele counts the test actually saw.

## Interpretation

A gene-based hit under this design is a candidate, not a confirmed association. tuning.rv fixed the minor-allele-count threshold at the point where the adjusted apparent group effect reaches a statistical null, which removes the genome-wide burden artefact; it does not make any individual gene real. What these figures add to the P-value is the number of people behind it and an external frequency check.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| loci | 49 |
| genome-wide loci | 0 |
| suggestive loci | 31 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-6 |

## Full statistics

**Per-locus values**

| peak | tier | gene | impact stratum | test method | variants in gene | cumulative MAC, case | cumulative MAC, control | cumulative MAF, case | cumulative MAF, control | P | FDR | burden OR | burden OR 95% CI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sug001_CD14 | suggestive | CD14 | low_moderate_high | cmc | 3 | 5 | 2 | 0.005708 | 3.800e-04 | — | — | 17.95 | [3.377, 95.39] |
| sug001_CEACAM7 | suggestive | CEACAM7 | low_moderate_high | skato | 5 | 5 | 8 | 0.005708 | 0.001521 | — | — | 3.951 | [1.255, 12.43] |
| sug001_NKD2 | suggestive | NKD2 | moderate_high | skato | 3 | 11 | 24 | 0.01256 | 0.004563 | — | — | 3.426 | [1.62, 7.243] |
| sug002_NFATC3 | suggestive | NFATC3 | moderate_high | skato | 3 | 5 | 2 | 0.005708 | 3.800e-04 | — | — | 13.72 | [2.624, 71.73] |
| sug002_RNF150 | suggestive | RNF150 | low_moderate_high | skato | 6 | 19 | 48 | 0.02169 | 0.009125 | — | — | 2.719 | [1.56, 4.74] |
| sug003_NKD2 | suggestive | NKD2 | low_moderate_high | skato | 3 | 11 | 24 | 0.01256 | 0.004563 | — | — | 3.426 | [1.62, 7.243] |
| sug003_THBS4 | suggestive | THBS4 | moderate_high | skato | 11 | 20 | 35 | 0.02283 | 0.006654 | — | — | 2.918 | [1.601, 5.318] |
| sug004_MMP20 | suggestive | MMP20 | moderate_high | skato | 5 | 4 | 10 | 0.004566 | 0.001901 | — | — | 2.298 | [0.7037, 7.507] |
| sug004_YTHDC2 | suggestive | YTHDC2 | low_moderate_high | skato | 11 | 16 | 35 | 0.01826 | 0.006654 | — | — | 2.753 | [1.484, 5.106] |
| sug005_ARHGEF10 | suggestive | ARHGEF10 | low_moderate_high | skato | 36 | 34 | 172 | 0.03881 | 0.0327 | — | — | 1.196 | [0.8063, 1.774] |
| sug005_NETO1 | suggestive | NETO1 | moderate_high | skato | 5 | 5 | 6 | 0.005708 | 0.001141 | — | — | 5.49 | [1.632, 18.47] |
| sug006_PNPT1 | suggestive | PNPT1 | moderate_high | skato | 6 | 12 | 25 | 0.0137 | 0.004753 | — | — | 2.836 | [1.365, 5.891] |
| sug006_TBC1D14 | suggestive | TBC1D14 | low_moderate_high | skato | 8 | 10 | 23 | 0.01142 | 0.004373 | — | — | 3.384 | [1.539, 7.44] |
| sug007_NFATC3 | suggestive | NFATC3 | low_moderate_high | skato | 5 | 7 | 6 | 0.007991 | 0.001141 | — | — | 7.05 | [2.306, 21.55] |
| sug007_TMEM131L | suggestive | TMEM131L | moderate_high | skato | 14 | 15 | 59 | 0.01712 | 0.01122 | — | — | 1.567 | [0.8704, 2.822] |
| sug008_NETO1 | suggestive | NETO1 | low_moderate_high | skato | 8 | 9 | 11 | 0.01027 | 0.002091 | — | — | 4.975 | [1.916, 12.92] |
| sug008_PELP1 | suggestive | PELP1 | moderate_high | skato | 5 | 23 | 55 | 0.02626 | 0.01046 | — | — | 2.448 | [1.471, 4.074] |
| sug009_COPS2 | suggestive | COPS2 | low_moderate_high | skato | 5 | 19 | 34 | 0.02169 | 0.006464 | — | — | 3.12 | [1.746, 5.574] |
| sug009_OR9A4 | suggestive | OR9A4 | moderate_high | skato | 3 | 5 | 2 | 0.005708 | 3.800e-04 | — | — | 13.69 | [2.605, 71.97] |
| sug010_DPP6 | suggestive | DPP6 | low_moderate_high | skato | 372 | 553 | 2,954 | 0.6313 | 0.5616 | — | — | 1.071 | [0.8698, 1.32] |
| sug010_GTPBP8 | suggestive | GTPBP8 | moderate_high | skato | 3 | 14 | 26 | 0.01598 | 0.004943 | — | — | 3.565 | [1.806, 7.035] |
| sug011_FRS3 | suggestive | FRS3 | moderate_high | skato | 7 | 14 | 30 | 0.01598 | 0.005703 | — | — | 2.874 | [1.484, 5.566] |
| sug011_OR9A4 | suggestive | OR9A4 | low_moderate_high | skato | 3 | 5 | 2 | 0.005708 | 3.800e-04 | — | — | 13.69 | [2.605, 71.97] |
| sug012_EFCAB6 | suggestive | EFCAB6 | moderate_high | skato | 17 | 49 | 154 | 0.05594 | 0.02928 | — | — | 2.005 | [1.411, 2.848] |
| sug012_KDM1A | suggestive | KDM1A | low_moderate_high | skato | 7 | 6 | 22 | 0.006849 | 0.004183 | — | — | 1.619 | [0.645, 4.065] |
| sug013_MCM10 | suggestive | MCM10 | moderate_high | skato | 11 | 19 | 45 | 0.02169 | 0.008555 | — | — | 2.768 | [1.583, 4.838] |
| sug013_TADA2A | suggestive | TADA2A | low_moderate_high | skato | 5 | 9 | 22 | 0.01027 | 0.004183 | — | — | 2.704 | [1.22, 5.992] |
| sug014_SMG1 | suggestive | SMG1 | moderate_high | skato | 11 | 25 | 102 | 0.02854 | 0.01939 | — | — | 1.385 | [0.8639, 2.222] |
| sug015_KCNAB2 | suggestive | KCNAB2 | moderate_high | skato | 8 | 22 | 61 | 0.02511 | 0.0116 | — | — | 2.348 | [1.41, 3.912] |
| sug016_DYSF | suggestive | DYSF | moderate_high | skato | 24 | 30 | 167 | 0.03425 | 0.03175 | — | — | 1.031 | [0.6769, 1.569] |
| sug017_FAM47E-STBD1 | suggestive | FAM47E-STBD1 | moderate_high | skato | 6 | 13 | 26 | 0.01484 | 0.004943 | — | — | 3.234 | [1.623, 6.445] |
| sig001_FAM217B | significant | FAM217B | moderate_high | skato | 6 | 9 | 13 | 0.01027 | 0.002471 | — | — | 4.217 | [1.754, 10.14] |
| sig001_RPH3A | significant | RPH3A | low_moderate_high | skato | 11 | 25 | 66 | 0.02854 | 0.01255 | — | — | 2.314 | [1.428, 3.751] |
| sig001_STBD1 | significant | STBD1 | moderate_high | cmc | 3 | 10 | 7 | 0.01142 | 0.001331 | — | — | 9.601 | [3.565, 25.86] |
| sig002_STBD1 | significant | STBD1 | low_moderate_high | skato | 8 | 13 | 21 | 0.01484 | 0.003992 | — | — | 4.059 | [1.984, 8.304] |
| sig002_STBD1 | significant | STBD1 | moderate_high | skato | 3 | 10 | 7 | 0.01142 | 0.001331 | — | — | 9.601 | [3.565, 25.86] |
| sig003_FAS | significant | FAS | low_moderate_high | skato | 3 | 4 | 2 | 0.004566 | 3.800e-04 | — | — | 15.38 | [2.74, 86.32] |
| sig003_RPH3A | significant | RPH3A | moderate_high | skato | 7 | 22 | 51 | 0.02511 | 0.009696 | — | — | 2.715 | [1.606, 4.589] |
| sig004_FAM217B | significant | FAM217B | low_moderate_high | skato | 9 | 10 | 25 | 0.01142 | 0.004753 | — | — | 2.495 | [1.169, 5.322] |
| sig004_JAK2 | significant | JAK2 | moderate_high | skato | 7 | 11 | 39 | 0.01256 | 0.007414 | — | — | 1.564 | [0.7656, 3.194] |
| sig005_PIK3C2G | significant | PIK3C2G | low_moderate_high | skato | 21 | 26 | 82 | 0.02968 | 0.01559 | — | — | 2.252 | [1.413, 3.591] |
| sig005_PIK3C2G | significant | PIK3C2G | moderate_high | skato | 12 | 22 | 62 | 0.02511 | 0.01179 | — | — | 2.566 | [1.535, 4.289] |
| sig006_CD14 | significant | CD14 | low_moderate_high | skato | 3 | 5 | 2 | 0.005708 | 3.800e-04 | — | — | 17.95 | [3.377, 95.39] |
| sig006_CEACAM7 | significant | CEACAM7 | moderate_high | skato | 4 | 5 | 5 | 0.005708 | 9.510e-04 | — | — | 6.179 | [1.723, 22.15] |
| sig007_P2RY2 | significant | P2RY2 | moderate_high | skato | 4 | 6 | 4 | 0.006849 | 7.600e-04 | — | — | 9.416 | [2.541, 34.9] |
| sig007_SCFD1 | significant | SCFD1 | low_moderate_high | skato | 6 | 14 | 45 | 0.01598 | 0.008555 | — | — | 1.947 | [1.026, 3.696] |
| sig008_SCAMP3 | significant | SCAMP3 | low_moderate_high | skato | 10 | 32 | 77 | 0.03653 | 0.01464 | — | — | 2.522 | [1.633, 3.896] |
| sig008_SCFD1 | significant | SCFD1 | moderate_high | skato | 4 | 10 | 36 | 0.01142 | 0.006844 | — | — | 1.991 | [0.9624, 4.119] |
| sig009_STARD10 | significant | STARD10 | moderate_high | skato | 3 | 5 | 3 | 0.005708 | 5.700e-04 | — | — | 11.81 | [2.688, 51.9] |

## How to read it

1. Read the allele counts in (c) first. A gene-level P of 1e-7 computed from a dozen carriers is a statement about a dozen people, and the odds ratio that comes with it has an interval spanning an order of magnitude.
2. Check (b) before believing (a). Under a fully confounded design the control side is where a technical artefact shows, and ToMMo is the only external reference available.
3. In (a), see whether the excess is spread across variants or carried by one. Burden tests weight every qualifying variant equally, so a single variant can carry a whole gene.

## What this figure does *not* establish

- It cannot establish that a variant is causal, or that the gene is. A burden test reports that qualifying variants are unequally distributed between groups.
- ToMMo agreement rules out one artefact class, not all of them. Residual platform chemistry is not fixable from within this data set.
- The odds ratio is fitted by logistic regression on the same covariates as the scan, but at these allele counts it is a description of the carriers observed, not a population estimate.

## Symbols

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **EAF** — frequency of the EFFECT allele (plink2 A1) among called samples of that group. Case and control EAF are the two numbers the odds ratio is computed from, so plotting them directly shows what drives an effect estimate and whether one group carries the whole difference.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
