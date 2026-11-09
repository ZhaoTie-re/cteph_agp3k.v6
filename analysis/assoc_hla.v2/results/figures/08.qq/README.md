# Supplementary: P values against the null

**Supplementary figure.** **Covers:** `qq.png` in `08.qq/`

## The question this figure answers

How far does each scan depart from the uniform null?

## Panels

**(qq.png a-h) observed against expected −log10 P**

Rows case sets, columns classes; one curve per cohort (random model, one point per test); grey: pointwise 95 % band of the uniform null, Beta(i, n − i + 1); corner: λGC per cohort.

## Interpretation

λGC spans 1.12-3.05 over the 24 scans (random model), 1.12-1.53 in the case-only test. All markers come from one region of chromosome 6 and are strongly correlated, so a few real or artefactual differences between the groups lift the whole curve; in the genome-wide fixed scan the MHC SNPs themselves have λGC 1.89 against 1.13 elsewhere (docs/STUDY_NOTES.md §8). λGC here describes the region, not the calibration of the test, and is never applied.

## Values in this rendering

| quantity | value |
|---|---|
| qq: all cases allele narrow: tests, lambda_GC | 143, 1.747 |
| qq: all cases allele intermediate: tests, lambda_GC | 144, 1.728 |
| qq: all cases allele full: tests, lambda_GC | 142, 1.374 |
| qq: all cases residue narrow: tests, lambda_GC | 494, 2.090 |
| qq: all cases residue intermediate: tests, lambda_GC | 503, 2.060 |
| qq: all cases residue full: tests, lambda_GC | 512, 2.199 |
| qq: DVT yes allele narrow: tests, lambda_GC | 143, 1.464 |
| qq: DVT yes allele intermediate: tests, lambda_GC | 144, 1.699 |
| qq: DVT yes allele full: tests, lambda_GC | 142, 1.414 |
| qq: DVT yes residue narrow: tests, lambda_GC | 492, 1.953 |
| qq: DVT yes residue intermediate: tests, lambda_GC | 501, 2.028 |
| qq: DVT yes residue full: tests, lambda_GC | 511, 1.681 |
| qq: DVT no allele narrow: tests, lambda_GC | 143, 1.662 |
| qq: DVT no allele intermediate: tests, lambda_GC | 144, 1.564 |
| qq: DVT no allele full: tests, lambda_GC | 142, 1.549 |
| qq: DVT no residue narrow: tests, lambda_GC | 491, 3.047 |
| qq: DVT no residue intermediate: tests, lambda_GC | 500, 2.632 |
| qq: DVT no residue full: tests, lambda_GC | 510, 2.239 |
| qq: DVT yes vs no, cases only allele narrow: tests, lambda_GC | 143, 1.142 |
| qq: DVT yes vs no, cases only allele intermediate: tests, lambda_GC | 144, 1.153 |
| qq: DVT yes vs no, cases only allele full: tests, lambda_GC | 142, 1.120 |
| qq: DVT yes vs no, cases only residue narrow: tests, lambda_GC | 457, 1.430 |
| qq: DVT yes vs no, cases only residue intermediate: tests, lambda_GC | 463, 1.530 |
| qq: DVT yes vs no, cases only residue full: tests, lambda_GC | 462, 1.438 |

## Full statistics

**λGC per scan**

| scan | tests, lambda_GC |
|---|---|
| all cases allele narrow: tests, lambda_GC | 143, 1.747 |
| all cases allele intermediate: tests, lambda_GC | 144, 1.728 |
| all cases allele full: tests, lambda_GC | 142, 1.374 |
| all cases residue narrow: tests, lambda_GC | 494, 2.090 |
| all cases residue intermediate: tests, lambda_GC | 503, 2.060 |
| all cases residue full: tests, lambda_GC | 512, 2.199 |
| DVT yes allele narrow: tests, lambda_GC | 143, 1.464 |
| DVT yes allele intermediate: tests, lambda_GC | 144, 1.699 |
| DVT yes allele full: tests, lambda_GC | 142, 1.414 |
| DVT yes residue narrow: tests, lambda_GC | 492, 1.953 |
| DVT yes residue intermediate: tests, lambda_GC | 501, 2.028 |
| DVT yes residue full: tests, lambda_GC | 511, 1.681 |
| DVT no allele narrow: tests, lambda_GC | 143, 1.662 |
| DVT no allele intermediate: tests, lambda_GC | 144, 1.564 |
| DVT no allele full: tests, lambda_GC | 142, 1.549 |
| DVT no residue narrow: tests, lambda_GC | 491, 3.047 |
| DVT no residue intermediate: tests, lambda_GC | 500, 2.632 |
| DVT no residue full: tests, lambda_GC | 510, 2.239 |
| DVT yes vs no, cases only allele narrow: tests, lambda_GC | 143, 1.142 |
| DVT yes vs no, cases only allele intermediate: tests, lambda_GC | 144, 1.153 |
| DVT yes vs no, cases only allele full: tests, lambda_GC | 142, 1.120 |
| DVT yes vs no, cases only residue narrow: tests, lambda_GC | 457, 1.430 |
| DVT yes vs no, cases only residue intermediate: tests, lambda_GC | 463, 1.530 |
| DVT yes vs no, cases only residue full: tests, lambda_GC | 462, 1.438 |

## How to read it

1. Compare the three curves in a panel: the cohorts should behave alike.
2. Do not read λGC as genome-wide inflation.

## What this figure does *not* establish

- Whether the departure is biology or artefact: see the platform limit below.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
