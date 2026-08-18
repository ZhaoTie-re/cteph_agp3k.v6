# Conditional analysis — full_mainland, peak gw003_2_11891017

**Figure file:** `conditional.gw003_2_11891017.png`

## The question this figure answers

Does this peak carry one association signal, or more than one?

## Panels

**(a) The peak, by conditioning round**

Round 0 is the unconditioned fit restricted to the peak window. Each later round adds the previous round's top variant to the covariate set and re-fits every variant in the window. Same samples, same model and the same covariate set as the genome-wide scan, so the rounds are comparable to it and to each other.

**(b) The stepwise decision**

The statistic the procedure actually acted on: the smallest P remaining at each round, against the same genome-wide threshold. The loop stops at the first round in which nothing in the window clears it.

## Interpretation

Collapse of the whole locus after conditioning on the lead means every significant variant there was tagging one underlying signal — the expected outcome inside a single LD block. A residual peak would mark a second, independent signal. Conditioning is performed on the same samples, model and covariates as the genome-wide scan, so the rounds are comparable to it and to each other.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| peak | gw003_2_11891017 |
| lead variant | chr2:11891017:G:A |
| rounds run | 2 |
| independent signals | 1 |
| max variants conditioned on | 1 |

## How to read it

1. If the whole peak collapses after conditioning on the lead, every significant variant there was tagging one underlying signal — the expected result inside a single LD block.
2. A residual peak that still clears the threshold is a second, independent signal, and the procedure will have added it to the conditioning set and continued.
3. Panel (b) is the audit trail: it shows the number the stopping rule compared, round by round.

## What this figure does *not* establish

- It cannot separate two causal variants in near-perfect LD; conditioning on one removes both. That limit is a property of the sample, not of the method.
- Absence of a secondary signal at this N is weak evidence — the conditional test faces the same detection floor as the primary scan, applied to a residual effect.

## Symbols

- **conditional analysis** — the locus re-fitted with the lead variant (then each further signal) added as a covariate. A signal that survives conditioning is independent of the lead; one that vanishes was tagging it.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^10gamma_k,PC_k,i      (g_i = genotype under the stated model; PCs = bbj_mainland)
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
