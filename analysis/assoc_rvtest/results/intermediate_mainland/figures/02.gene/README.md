# Gene-based hits — intermediate_mainland

**Figures:** `02.gene/<tier>/gene.<peak>.png` — 35 locus figure(s): 0 genome-wide, 18 suggestive.

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
| cohort | intermediate_mainland |
| loci | 35 |
| genome-wide loci | 0 |
| suggestive loci | 18 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-6 |

## Full statistics

**Per-locus values**

| peak | tier | gene | impact stratum | test method | variants in gene | cumulative MAC, case | cumulative MAC, control | cumulative MAF, case | cumulative MAF, control | P | FDR | burden OR | burden OR 95% CI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sug001_CDYL2 | suggestive | CDYL2 | moderate_high | skato | 3 | 9 | 7 | 0.01051 | 0.001708 | — | — | 4.96 | [1.686, 14.59] |
| sug001_EFCAB6 | suggestive | EFCAB6 | low_moderate_high | cmc | 18 | 51 | 119 | 0.05958 | 0.02904 | — | — | 2.176 | [1.52, 3.115] |
| sug001_FAS | suggestive | FAS | low_moderate_high | skato | 3 | 4 | 2 | 0.004673 | 4.880e-04 | — | — | 12.5 | [2.23, 70.01] |
| sug001_STBD1 | suggestive | STBD1 | moderate_high | cmc | 3 | 10 | 7 | 0.01168 | 0.001708 | — | — | 7.349 | [2.705, 19.96] |
| sug002_EFCAB6 | suggestive | EFCAB6 | moderate_high | skato | 14 | 48 | 113 | 0.05607 | 0.02757 | — | — | 2.14 | [1.481, 3.093] |
| sug002_STBD1 | suggestive | STBD1 | low_moderate_high | cmc | 6 | 13 | 14 | 0.01519 | 0.003416 | — | — | 4.767 | [2.178, 10.44] |
| sug002_TH | suggestive | TH | low_moderate_high | skato | 8 | 18 | 38 | 0.02103 | 0.009273 | — | — | 2.398 | [1.341, 4.289] |
| sug003_FAM217B | suggestive | FAM217B | moderate_high | skato | 5 | 8 | 9 | 0.009346 | 0.002196 | — | — | 3.95 | [1.488, 10.49] |
| sug003_OSBP | suggestive | OSBP | low_moderate_high | skato | 6 | 12 | 13 | 0.01402 | 0.003172 | — | — | 4.245 | [1.887, 9.548] |
| sug003_TBC1D14 | suggestive | TBC1D14 | low_moderate_high | cmc | 4 | 9 | 8 | 0.01051 | 0.001952 | — | — | 6.659 | [2.448, 18.11] |
| sug004_KDM1A | suggestive | KDM1A | low_moderate_high | skato | 5 | 6 | 12 | 0.007009 | 0.002928 | — | — | 2.47 | [0.9098, 6.703] |
| sug004_PNPT1 | suggestive | PNPT1 | moderate_high | skato | 6 | 12 | 22 | 0.01402 | 0.005368 | — | — | 2.828 | [1.337, 5.982] |
| sug005_COPS2 | suggestive | COPS2 | low_moderate_high | skato | 4 | 19 | 28 | 0.0222 | 0.006833 | — | — | 3.265 | [1.787, 5.966] |
| sug005_TH | suggestive | TH | moderate_high | skato | 4 | 12 | 27 | 0.01402 | 0.006589 | — | — | 2.273 | [1.128, 4.583] |
| sug006_CEACAM7 | suggestive | CEACAM7 | moderate_high | skato | 3 | 5 | 3 | 0.005841 | 7.320e-04 | — | — | 8.459 | [1.965, 36.41] |
| sug006_YTHDC2 | suggestive | YTHDC2 | low_moderate_high | skato | 8 | 15 | 22 | 0.01752 | 0.005368 | — | — | 3.107 | [1.575, 6.128] |
| sug007_THBS4 | suggestive | THBS4 | moderate_high | skato | 10 | 20 | 28 | 0.02336 | 0.006833 | — | — | 3.025 | [1.618, 5.653] |
| sug008_SCFD1 | suggestive | SCFD1 | moderate_high | skato | 4 | 10 | 25 | 0.01168 | 0.006101 | — | — | 2.164 | [1.012, 4.627] |
| sig001_PIK3C2G | significant | PIK3C2G | low_moderate_high | skato | 17 | 20 | 46 | 0.02336 | 0.01123 | — | — | 2.377 | [1.37, 4.124] |
| sig001_PIK3C2G | significant | PIK3C2G | moderate_high | skato | 11 | 16 | 34 | 0.01869 | 0.008297 | — | — | 2.674 | [1.437, 4.979] |
| sig001_SCAMP3 | significant | SCAMP3 | low_moderate_high | cmc | 8 | 32 | 53 | 0.03738 | 0.01293 | — | — | 2.9 | [1.826, 4.606] |
| sig002_RPH3A | significant | RPH3A | low_moderate_high | skato | 11 | 25 | 53 | 0.02921 | 0.01293 | — | — | 2.314 | [1.405, 3.809] |
| sig002_RPH3A | significant | RPH3A | moderate_high | skato | 7 | 22 | 39 | 0.0257 | 0.009517 | — | — | 2.826 | [1.634, 4.888] |
| sig003_JAK2 | significant | JAK2 | moderate_high | skato | 6 | 11 | 32 | 0.01285 | 0.007809 | — | — | 1.425 | [0.6864, 2.957] |
| sig003_TBC1D14 | significant | TBC1D14 | low_moderate_high | skato | 4 | 9 | 8 | 0.01051 | 0.001952 | — | — | 6.659 | [2.448, 18.11] |
| sig004_PES1 | significant | PES1 | moderate_high | skato | 4 | 28 | 52 | 0.03271 | 0.01269 | — | — | 2.671 | [1.637, 4.356] |
| sig004_SCAMP3 | significant | SCAMP3 | low_moderate_high | skato | 8 | 32 | 53 | 0.03738 | 0.01293 | — | — | 2.9 | [1.826, 4.606] |
| sig005_CD248 | significant | CD248 | low_moderate_high | skato | 6 | 11 | 22 | 0.01285 | 0.005368 | — | — | 2.597 | [1.225, 5.507] |
| sig005_STBD1 | significant | STBD1 | moderate_high | skato | 3 | 10 | 7 | 0.01168 | 0.001708 | — | — | 7.349 | [2.705, 19.96] |
| sig006_OTOP3 | significant | OTOP3 | moderate_high | skato | 8 | 26 | 57 | 0.03037 | 0.01391 | — | — | 2.285 | [1.406, 3.712] |
| sig006_STBD1 | significant | STBD1 | low_moderate_high | skato | 6 | 13 | 14 | 0.01519 | 0.003416 | — | — | 4.767 | [2.178, 10.44] |
| sig007_NKD2 | significant | NKD2 | moderate_high | skato | 3 | 11 | 16 | 0.01285 | 0.003904 | — | — | 3.896 | [1.747, 8.687] |
| sig007_PES1 | significant | PES1 | low_moderate_high | skato | 9 | 37 | 88 | 0.04322 | 0.02147 | — | — | 2.143 | [1.419, 3.237] |
| sig008_NKD2 | significant | NKD2 | low_moderate_high | skato | 3 | 11 | 16 | 0.01285 | 0.003904 | — | — | 3.896 | [1.747, 8.687] |
| sig009_EFCAB6 | significant | EFCAB6 | low_moderate_high | skato | 18 | 51 | 119 | 0.05958 | 0.02904 | — | — | 2.176 | [1.52, 3.115] |

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
