# Cases with and without DVT

**Main figure.** **Covers:** `dvt.png` in `05.dvt/`

## The question this figure answers

Does any HLA allele or residue associate differently with CTEPH in cases with a history of DVT than in cases without?

## Panels

**(dvt.png a-f) log OR with DVT against without DVT**

Rows alleles, residues; columns narrow, intermediate, full mainland. x: log OR of cases without DVT against all controls; y: cases with DVT against the same controls; random model, one point per test. Blue / red: significant at the LD-adjusted threshold in DVT no / DVT yes; ring: significant in the case-only test. Grey line: identity. Signal leads are labelled.

## Interpretation

Random model, one point per test. r(log OR) between the two sets is 0.13-0.24 over the six panels. Points far from the identity line differ between cases with and without DVT, but the two sets share their controls, so the scatter is not a test; the case-only test (rings) is. The DVT sets are small (140-159 cases), so a marker significant in one set and not the other is expected by chance for some markers.

## Values in this rendering

| quantity | value |
|---|---|
| dvt: allele narrow: tests, r(log OR) | 143, 0.132 |
| dvt: allele narrow: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 1 / 0 |
| dvt: allele intermediate: tests, r(log OR) | 144, 0.149 |
| dvt: allele intermediate: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 0 / 0 |
| dvt: allele full: tests, r(log OR) | 142, 0.239 |
| dvt: allele full: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 0 / 0 |
| dvt: residue narrow: tests, r(log OR) | 506, 0.220 |
| dvt: residue narrow: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 34 / 0 |
| dvt: residue intermediate: tests, r(log OR) | 511, 0.192 |
| dvt: residue intermediate: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 0 / 0 |
| dvt: residue full: tests, r(log OR) | 510, 0.170 |
| dvt: residue full: points drawn GEC-significant in DVT yes / DVT no / case-only | 0 / 0 / 0 |

## Full statistics

**Signal leads in every comparison (random model)**

| marker | cohort | OR all cases vs controls | P all cases vs controls | OR DVT yes vs controls | P DVT yes vs controls | OR DVT no vs controls | P DVT no vs controls | OR DVT yes vs DVT no | P DVT yes vs DVT no |
|---|---|---|---|---|---|---|---|---|---|
| DPB1*05:01 | narrow | 1.274 | 0.004197 | 1.277 | 0.05819 | 1.629 | 2.815e-04 | 0.773 | 0.142 |
| DPB1*05:01 | intermediate | 1.203 | 0.02306 | 1.162 | 0.2298 | 1.554 | 8.435e-04 | 0.748 | 0.09191 |
| DPB1*05:01 | full | 1.174 | 0.04341 | 1.115 | 0.3736 | 1.507 | 0.001415 | 0.75 | 0.08819 |
| DQB1 57 V | narrow | 0.621 | 4.570e-05 | 0.644 | 0.0112 | 0.521 | 3.630e-04 | 1.469 | 0.2174 |
| DQB1 57 V | intermediate | 0.652 | 1.847e-04 | 0.672 | 0.02043 | 0.549 | 0.001027 | 1.434 | 0.2382 |
| DQB1 57 V | full | 0.656 | 1.741e-04 | 0.688 | 0.02701 | 0.55 | 7.137e-04 | 1.441 | 0.2196 |
| DQB1 87 Y | narrow | 0.721 | 0.001152 | 0.808 | 0.1585 | 0.543 | 1.190e-04 | 1.909 | 0.0123 |
| DQB1 87 Y | intermediate | 0.737 | 0.001736 | 0.819 | 0.1756 | 0.569 | 3.290e-04 | 1.765 | 0.02277 |
| DQB1 87 Y | full | 0.742 | 0.001852 | 0.838 | 0.2241 | 0.577 | 3.149e-04 | 1.704 | 0.02864 |

## How to read it

1. Points on the identity line: the same effect with and without DVT.
2. A coloured point off the line: significant in one set with a different estimate in the other.
3. A ring: the case-only test supports a difference.

## What this figure does *not* establish

- DVT status is recorded only for cases; about 30 % of cases have unknown status and are in neither set.
- A difference read off the scatter alone: the two axes share controls; use the case-only test.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
