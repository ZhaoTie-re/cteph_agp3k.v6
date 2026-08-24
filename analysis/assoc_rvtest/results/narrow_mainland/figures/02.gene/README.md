# Gene-based hits — narrow_mainland

**Figures:** `02.gene/<tier>/gene.<peak>.png` — 24 locus figure(s): 0 genome-wide, 5 suggestive.

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
| cohort | narrow_mainland |
| loci | 24 |
| genome-wide loci | 0 |
| suggestive loci | 5 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-6 |

## Full statistics

**Per-locus values**

| peak | tier | gene | impact stratum | test method | variants in gene | cumulative MAC, case | cumulative MAC, control | cumulative MAF, case | cumulative MAF, control | P | FDR | burden OR | burden OR 95% CI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sug001_CHRND | suggestive | CHRND | low_moderate_high | skato | 5 | 17 | 27 | 0.02033 | 0.007728 | — | — | 2.72 | [1.437, 5.149] |
| sug001_FRS3 | suggestive | FRS3 | moderate_high | skato | 7 | 14 | 19 | 0.01675 | 0.005438 | — | — | 3.188 | [1.554, 6.539] |
| sug001_STBD1 | suggestive | STBD1 | low_moderate_high | cmc | 5 | 12 | 8 | 0.01435 | 0.00229 | — | — | 6.946 | [2.729, 17.68] |
| sug002_TH | suggestive | TH | low_moderate_high | skato | 7 | 18 | 30 | 0.02153 | 0.008586 | — | — | 2.604 | [1.412, 4.802] |
| sug003_MYH13 | suggestive | MYH13 | low_moderate_high | skato | 26 | 37 | 112 | 0.04426 | 0.03206 | — | — | 1.397 | [0.9297, 2.101] |
| sig001_PIK3C2G | significant | PIK3C2G | low_moderate_high | skato | 15 | 17 | 39 | 0.02033 | 0.01116 | — | — | 1.968 | [1.078, 3.591] |
| sig001_PIK3C2G | significant | PIK3C2G | moderate_high | skato | 10 | 15 | 27 | 0.01794 | 0.007728 | — | — | 2.648 | [1.361, 5.153] |
| sig002_CD248 | significant | CD248 | low_moderate_high | skato | 6 | 11 | 20 | 0.01316 | 0.005724 | — | — | 2.512 | [1.152, 5.477] |
| sig002_JAK2 | significant | JAK2 | moderate_high | skato | 6 | 11 | 23 | 0.01316 | 0.006583 | — | — | 1.707 | [0.7868, 3.702] |
| sig003_NKD2 | significant | NKD2 | low_moderate_high | skato | 3 | 11 | 12 | 0.01316 | 0.003434 | — | — | 5.007 | [2.121, 11.82] |
| sig003_NKD2 | significant | NKD2 | moderate_high | skato | 3 | 11 | 12 | 0.01316 | 0.003434 | — | — | 5.007 | [2.121, 11.82] |
| sig004_PES1 | significant | PES1 | moderate_high | skato | 4 | 28 | 44 | 0.03349 | 0.01259 | — | — | 2.789 | [1.67, 4.656] |
| sig004_STBD1 | significant | STBD1 | low_moderate_high | skato | 5 | 12 | 8 | 0.01435 | 0.00229 | — | — | 6.946 | [2.729, 17.68] |
| sig005_SCFD1 | significant | SCFD1 | moderate_high | skato | 3 | 9 | 17 | 0.01077 | 0.004865 | — | — | 2.679 | [1.149, 6.249] |
| sig005_ZKSCAN7 | significant | ZKSCAN7 | low_moderate_high | skato | 5 | 28 | 52 | 0.03349 | 0.01488 | — | — | 2.377 | [1.448, 3.903] |
| sig006_PES1 | significant | PES1 | low_moderate_high | skato | 9 | 37 | 75 | 0.04426 | 0.02147 | — | — | 2.171 | [1.415, 3.329] |
| sig006_SMG1 | significant | SMG1 | moderate_high | skato | 8 | 25 | 74 | 0.0299 | 0.02118 | — | — | 1.52 | [0.9243, 2.5] |
| sig007_RPH3A | significant | RPH3A | low_moderate_high | skato | 11 | 24 | 47 | 0.02871 | 0.01345 | — | — | 2.199 | [1.307, 3.701] |
| sig007_RPH3A | significant | RPH3A | moderate_high | skato | 7 | 21 | 35 | 0.02512 | 0.01002 | — | — | 2.615 | [1.475, 4.637] |
| sig008_SCFD1 | significant | SCFD1 | low_moderate_high | skato | 4 | 9 | 19 | 0.01077 | 0.005438 | — | — | 2.273 | [0.9913, 5.213] |
| sig008_STBD1 | significant | STBD1 | moderate_high | skato | 3 | 9 | 4 | 0.01077 | 0.001145 | — | — | 9.978 | [2.931, 33.96] |
| sig009_MYOM2 | significant | MYOM2 | low_moderate_high | skato | 42 | 75 | 189 | 0.08971 | 0.05409 | — | — | 1.764 | [1.296, 2.402] |
| sig010_EFCAB6 | significant | EFCAB6 | low_moderate_high | skato | 15 | 44 | 84 | 0.05263 | 0.02404 | — | — | 2.369 | [1.585, 3.541] |
| sig011_SMG1 | significant | SMG1 | low_moderate_high | skato | 17 | 33 | 116 | 0.03947 | 0.0332 | — | — | 1.327 | [0.866, 2.035] |

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
