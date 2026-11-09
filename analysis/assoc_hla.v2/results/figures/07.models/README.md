# Supplementary: fixed against random model

**Supplementary figure.** **Covers:** `models.png` in `07.models/`

## The question this figure answers

Does accounting for relatedness and the polygenic background (random model) change any test compared with the fixed model?

## Panels

**(models.png a-c) −log10 P, fixed (x) against random (y)**

One panel per cohort; one point per test, all four case sets (colour). Grey: identity.

**(models.png d-f) β, fixed against random**

Fixed: log OR from plink2. Random: SAIGE's one-step estimate Tstat / var.

## Interpretation

The two models differ in samples (the fixed model drops relatives), missing calls (complete case against mean imputation) and estimator (Wald against score), so points scatter around the identity; r(−log10 P) is 0.993-0.999 per case set, cohort and class. Tests significant at the LD-adjusted threshold in one model only: narrow all cases B*40:02; narrow DVT no DQB1:-27:A; narrow DVT no DQB1:-5:S; narrow DVT no DQB1:-9:I; full all cases B:114:D; full all cases B:114:N. The random model is primary; calls made by the fixed model alone are not counted.

## Values in this rendering

| quantity | value |
|---|---|
| models: narrow all cases: r(-log10 P) | 0.9985 |
| models: narrow DVT yes: r(-log10 P) | 0.9979 |
| models: narrow DVT no: r(-log10 P) | 0.9988 |
| models: narrow DVT yes vs no, cases only: r(-log10 P) | 0.9991 |
| models: intermediate all cases: r(-log10 P) | 0.9985 |
| models: intermediate DVT yes: r(-log10 P) | 0.9983 |
| models: intermediate DVT no: r(-log10 P) | 0.999 |
| models: intermediate DVT yes vs no, cases only: r(-log10 P) | 0.9991 |
| models: full all cases: r(-log10 P) | 0.9983 |
| models: full DVT yes: r(-log10 P) | 0.9982 |
| models: full DVT no: r(-log10 P) | 0.9991 |
| models: full DVT yes vs no, cases only: r(-log10 P) | 0.9982 |
| models: tests significant at GEC in one model only (cohort, case set, marker) | narrow all cases B*40:02; narrow DVT no DQB1:-27:A; narrow DVT no DQB1:-5:S; narrow DVT no DQB1:-9:I; full all cases B:114:D; full all cases B:114:N |

## Full statistics

**Per case set, cohort and class**

| case set | cohort | class | tests | r(-log10 P) | r(beta) | GEC fixed / random / both |
|---|---|---|---|---|---|---|
| all cases | narrow | allele | 143 | 0.9972 | 0.9759 | 1 / 0 / 0 |
| all cases | narrow | residue | 494 | 0.9988 | 0.9949 | 1 / 1 / 1 |
| all cases | intermediate | allele | 144 | 0.996 | 0.972 | 0 / 0 / 0 |
| all cases | intermediate | residue | 503 | 0.9989 | 0.9953 | 1 / 1 / 1 |
| all cases | full | allele | 142 | 0.9951 | 0.9694 | 0 / 0 / 0 |
| all cases | full | residue | 512 | 0.999 | 0.9955 | 3 / 1 / 1 |
| DVT yes | narrow | allele | 134 | 0.9983 | 0.9706 | 0 / 0 / 0 |
| DVT yes | narrow | residue | 486 | 0.9978 | 0.9708 | 0 / 0 / 0 |
| DVT yes | intermediate | allele | 135 | 0.9983 | 0.9717 | 0 / 0 / 0 |
| DVT yes | intermediate | residue | 495 | 0.9984 | 0.973 | 0 / 0 / 0 |
| DVT yes | full | allele | 137 | 0.9977 | 0.9619 | 0 / 0 / 0 |
| DVT yes | full | residue | 505 | 0.9985 | 0.9681 | 0 / 0 / 0 |
| DVT no | narrow | allele | 141 | 0.9969 | 0.973 | 1 / 1 / 1 |
| DVT no | narrow | residue | 491 | 0.999 | 0.9832 | 30 / 33 / 30 |
| DVT no | intermediate | allele | 142 | 0.9968 | 0.976 | 0 / 0 / 0 |
| DVT no | intermediate | residue | 500 | 0.9992 | 0.985 | 0 / 0 / 0 |
| DVT no | full | allele | 140 | 0.9976 | 0.9756 | 0 / 0 / 0 |
| DVT no | full | residue | 510 | 0.9993 | 0.9876 | 0 / 0 / 0 |
| DVT yes vs no, cases only | narrow | allele | 132 | 0.9969 | 0.9975 | 0 / 0 / 0 |
| DVT yes vs no, cases only | narrow | residue | 451 | 0.9994 | 0.9986 | 0 / 0 / 0 |
| DVT yes vs no, cases only | intermediate | allele | 133 | 0.997 | 0.9974 | 0 / 0 / 0 |
| DVT yes vs no, cases only | intermediate | residue | 457 | 0.9994 | 0.9986 | 0 / 0 / 0 |
| DVT yes vs no, cases only | full | allele | 135 | 0.9932 | 0.9925 | 0 / 0 / 0 |
| DVT yes vs no, cases only | full | residue | 457 | 0.9992 | 0.9969 | 0 / 0 / 0 |

## How to read it

1. Points on the identity: the models agree.
2. The table and the interpretation list the calls that differ.

## What this figure does *not* establish

- Which model is right for a marker: they answer the same question under different assumptions.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
