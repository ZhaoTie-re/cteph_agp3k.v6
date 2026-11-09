# Supplementary: the three thresholds

**Supplementary figure.** **Covers:** `thresholds.png` in `06.thresholds/`

## The question this figure answers

How strict is each threshold, and how close does each come to a family-wise error rate of 0.05?

## Panels

**(thresholds.png a) alleles: tests and effective numbers**

Per case set and cohort: the number of tests (conventional Bonferroni divides by it), the Me implied by the exact max-T threshold (0.05 / threshold), and GEC v0.2's Me (the LD-adjusted threshold divides by it).

**(thresholds.png b) residues: the same**

A residue mirror pair counts as one test.

**(thresholds.png c) family-wise error rate actually achieved**

For each of the 24 families, the probability under the max-T simulation (400,000 multivariate-normal draws on the covariate-adjusted genotype correlation) that any test passes the GEC (filled) or the conventional Bonferroni (open) threshold; dashed: the nominal 0.05.

## Interpretation

Over the 24 (case set, cohort, class) families the GEC threshold's family-wise error under max-T is 0.052-0.053 for alleles and 0.045-0.047 for residues; GEC v0.2 estimates Me within blocks of at most 100 markers and ignores LD between blocks. Conventional Bonferroni counts every test as independent: 0.041-0.043 for alleles and 0.019-0.021 for residues, whose mirror pairs and haplotypes make them strongly correlated.

## Values in this rendering

| quantity | value |
|---|---|
| thresholds: allele: tests / Me max-T / Me GEC (range over families) | 142–144 / 114.6–123.1 / 110.5–115.3 |
| thresholds: residue: tests / Me max-T / Me GEC (range over families) | 457–512 / 169.9–183.2 / 185.2–201.8 |
| thresholds: allele: family-wise error under max-T, GEC / Bonferroni | 0.052–0.053 / 0.041–0.043 |
| thresholds: residue: family-wise error under max-T, GEC / Bonferroni | 0.045–0.047 / 0.019–0.021 |

## Full statistics

**Every family**

| case_set | cohort | class | markers | mirror_pairs | tests | me_gec | me_maxt | p_gec | p_maxt | p_bonf | fwer_gec | fwer_bonf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DVT yes vs no, cases only | full | allele | 142 | 0 | 142 | 110.5 | 114.6 | 4.523e-04 | 4.362e-04 | 3.521e-04 | 0.052 | 0.041 |
| all cases | intermediate | allele | 144 | 0 | 144 | 115 | 120.6 | 4.349e-04 | 4.145e-04 | 3.472e-04 | 0.052 | 0.042 |
| DVT yes | narrow | residue | 672 | 180 | 492 | 200.2 | 180.4 | 2.497e-04 | 2.772e-04 | 1.016e-04 | 0.045 | 0.02 |
| DVT yes vs no, cases only | narrow | allele | 143 | 0 | 143 | 111.3 | 116 | 4.490e-04 | 4.309e-04 | 3.497e-04 | 0.052 | 0.041 |
| all cases | full | residue | 679 | 167 | 512 | 200.9 | 182.5 | 2.489e-04 | 2.740e-04 | 9.766e-05 | 0.046 | 0.02 |
| DVT yes | intermediate | residue | 681 | 180 | 501 | 200.5 | 182.2 | 2.494e-04 | 2.744e-04 | 9.980e-05 | 0.046 | 0.02 |
| DVT yes | intermediate | allele | 144 | 0 | 144 | 115.3 | 123.1 | 4.337e-04 | 4.063e-04 | 3.472e-04 | 0.053 | 0.043 |
| DVT yes vs no, cases only | intermediate | allele | 144 | 0 | 144 | 112.3 | 117.4 | 4.454e-04 | 4.260e-04 | 3.472e-04 | 0.052 | 0.041 |
| DVT no | narrow | residue | 672 | 181 | 491 | 200.1 | 180.5 | 2.499e-04 | 2.770e-04 | 1.018e-04 | 0.046 | 0.02 |
| DVT yes | narrow | allele | 143 | 0 | 143 | 114.3 | 120 | 4.374e-04 | 4.166e-04 | 3.497e-04 | 0.052 | 0.042 |
| all cases | intermediate | residue | 681 | 178 | 503 | 200.8 | 183.2 | 2.491e-04 | 2.730e-04 | 9.940e-05 | 0.046 | 0.02 |
| DVT yes | full | allele | 142 | 0 | 142 | 113.8 | 120.6 | 4.394e-04 | 4.146e-04 | 3.521e-04 | 0.053 | 0.043 |
| all cases | narrow | allele | 143 | 0 | 143 | 114 | 120.7 | 4.388e-04 | 4.144e-04 | 3.497e-04 | 0.053 | 0.042 |
| all cases | full | allele | 142 | 0 | 142 | 113.6 | 120.9 | 4.402e-04 | 4.134e-04 | 3.521e-04 | 0.053 | 0.043 |
| DVT no | intermediate | residue | 681 | 181 | 500 | 200.3 | 183.1 | 2.497e-04 | 2.731e-04 | 1.000e-04 | 0.046 | 0.02 |
| DVT no | full | residue | 679 | 169 | 510 | 200.4 | 180.7 | 2.495e-04 | 2.767e-04 | 9.804e-05 | 0.045 | 0.019 |
| DVT yes vs no, cases only | narrow | residue | 672 | 215 | 457 | 185.2 | 170.9 | 2.700e-04 | 2.925e-04 | 1.094e-04 | 0.047 | 0.021 |
| DVT no | narrow | allele | 143 | 0 | 143 | 114 | 119.7 | 4.385e-04 | 4.179e-04 | 3.497e-04 | 0.052 | 0.042 |
| DVT yes | full | residue | 679 | 168 | 511 | 200.6 | 181.3 | 2.493e-04 | 2.758e-04 | 9.785e-05 | 0.046 | 0.019 |
| DVT no | intermediate | allele | 144 | 0 | 144 | 115 | 121.1 | 4.346e-04 | 4.129e-04 | 3.472e-04 | 0.052 | 0.042 |
| DVT yes vs no, cases only | intermediate | residue | 681 | 218 | 463 | 186.6 | 169.9 | 2.680e-04 | 2.943e-04 | 1.080e-04 | 0.046 | 0.02 |
| DVT no | full | allele | 142 | 0 | 142 | 113.7 | 120.8 | 4.399e-04 | 4.139e-04 | 3.521e-04 | 0.053 | 0.043 |
| all cases | narrow | residue | 672 | 178 | 494 | 201.8 | 180.2 | 2.477e-04 | 2.775e-04 | 1.012e-04 | 0.045 | 0.02 |
| DVT yes vs no, cases only | full | residue | 679 | 217 | 462 | 186.5 | 172.3 | 2.681e-04 | 2.902e-04 | 1.082e-04 | 0.047 | 0.02 |

## How to read it

1. Compare the red dot (GEC) with the grey diamond (exact): the closer, the better calibrated.
2. Panel c: points on the dashed line hold the error rate exactly; below it, the threshold is stricter than needed.

## What this figure does *not* establish

- The max-T simulation assumes the test statistics are jointly normal with the genotype correlation; it is a calibration, not a permutation of phenotypes.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
