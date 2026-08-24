# Conditional analysis — narrow_mainland

**Figures:** `04.conditional/<tier>/conditional.<peak>.png` — 10 locus figure(s): 0 genome-wide, 10 suggestive.

One document covers the whole family: the panels, the reading and the limits are properties of the figure, identical for every locus. Only the numbers differ, and they are tabulated below, one row per locus.

## The question this figure answers

Does each peak carry one association signal, or more than one?

## Panels

**(a) The peak, by conditioning round**

Round 0 is the unconditioned fit restricted to the peak window. Each later round adds the previous round's top variant to the covariate set and re-fits every variant in the window. Same samples, same model and the same covariate set as the genome-wide scan, so the rounds are comparable to it and to each other.

**(b) The stepwise decision**

The statistic the procedure actually acted on: the smallest P remaining at each round, against the threshold of the peak's own tier. The loop stops at the first round in which nothing in the window clears it.

## Interpretation

Collapse of the whole locus after conditioning on the lead means every significant variant there was tagging one underlying signal — the expected outcome inside a single LD block. A residual peak would mark a second, independent signal. Conditioning is performed on the same samples, model and covariates as the genome-wide scan, so the rounds are comparable to it and to each other.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | narrow_mainland |
| loci | 10 |
| genome-wide loci | 0 |
| suggestive loci | 10 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-5 |

## Full statistics

**Per-locus values**

| peak | tier | gene | lead variant | rsID | decision threshold | rounds run | independent signals | max variants conditioned on |
|---|---|---|---|---|---|---|---|---|
| sg004_3_94020684 | suggestive | ARL13B | chr3:94020684:T:C | rs143483852 | 1e-5 | 2 | 1 | 1 |
| sg005_3_100783010 | suggestive | ABI3BP | chr3:100783010:T:C | rs7648469 | 1e-5 | 2 | 1 | 1 |
| sg007_3_154069965 | suggestive | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | 1e-5 | 2 | 1 | 1 |
| sg008_4_171998307 | suggestive | GALNTL6 | chr4:171998307:T:C | rs75214171 | 1e-5 | 2 | 1 | 1 |
| sg015_8_138575524 | suggestive | COL22A1 | chr8:138575524:T:C | rs62530012 | 1e-5 | 2 | 1 | 1 |
| sg022_14_88019710 | suggestive | LINC01146 | chr14:88019710:G:A | rs4899932 | 1e-5 | 2 | 1 | 1 |
| sg023_16_53882967 | suggestive | FTO | chr16:53882967:G:A | rs9934504 | 1e-5 | 2 | 1 | 1 |
| sg026_19_10631611 | suggestive | SLC44A2 | chr19:10631611:C:T | rs1560711 | 1e-5 | 2 | 1 | 1 |
| sg028_20_7553182 | suggestive | RN7SL547P | chr20:7553182:C:T | rs117662178 | 1e-5 | 2 | 1 | 1 |
| sg029_20_60915853 | suggestive | CDH4 | chr20:60915853:C:T | rs79313991 | 1e-5 | 2 | 1 | 1 |

## How to read it

1. If the whole peak collapses after conditioning on the lead, every significant variant there was tagging one underlying signal — the expected result inside a single LD block.
2. A residual peak that still clears the threshold is a second, independent signal, and the procedure will have added it to the conditioning set and continued.
3. Panel (b) is the audit trail: it shows the number the stopping rule compared, round by round.
4. Both tier lines are drawn; the one the peak was judged against — 5e-8 for a genome-wide peak, 1e-5 for a suggestive one — is drawn heavier and marked "decision".

## What this figure does *not* establish

- It cannot separate two causal variants in near-perfect LD; conditioning on one removes both. That limit is a property of the sample, not of the method.
- Absence of a secondary signal at this N is weak evidence — the conditional test faces the same detection floor as the primary scan, applied to a residual effect.
- A suggestive peak is judged against the suggestive threshold, so its signal count is not comparable to a genome-wide peak's without accounting for the different bar.

## Symbols

- **conditional analysis** — the locus re-fitted with the lead variant (then each further signal) added as a covariate. A signal that survives conditioning is independent of the lead; one that vanishes was tagging it.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^10gamma_k,PC_k,i      (g_i = genotype under the stated model; PCs = bbj_mainland)
```

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
