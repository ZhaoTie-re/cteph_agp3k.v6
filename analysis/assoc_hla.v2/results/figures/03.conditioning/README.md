# Conditional analyses: independent signals and their carriers

**Main figure.** **Covers:** `carriers.png`, `conditioning.png` in `03.conditioning/`

## The question this figure answers

Does one marker explain each set of significant markers, and which alleles carry each residue signal?

## Panels

**(conditioning.png) each walk before and after conditioning**

One row per (case set, class) walk that entered a signal, random model. Left: the marginal scan (step 0); right: the scan conditional on every entered lead (last step). Bars as in 01.landscape (largest to smallest P over the cohorts, coloured by n/3 at GEC). The triangle on the right marks the conditioned lead.

**(carriers.png) the alleles that carry each residue signal**

One panel per residue signal: every tested allele of the gene, OR per copy with 95 % CI in each cohort (random model, the signal's own comparison), carriers of the residue first. Carriers are read from hla.typing.v2: allele_match_map (the reference allele a 2-field allele was matched to) and residue_reference (that allele's residue at the position).

## Interpretation

Random model. Signals: DQB1 57 V (all cases, residue: GEC 3/3, max-T 3/3, Bonferroni 1/3); DQB1 87 Y (DVT no, residue: GEC 1/3, max-T 1/3, Bonferroni 0/3); DPB1*05:01 (DVT no, allele: GEC 1/3, max-T 1/3, Bonferroni 1/3). After conditioning on the lead(s), the number of tests still significant at the LD-adjusted threshold in some cohort is all cases residue: 0, DVT no allele: 0, DVT no residue: 0. A walk that enters one signal and leaves nothing explains all its significant markers with one lead. The carrier panels show which tested alleles carry each residue signal: a residue-level explanation predicts that the carriers share the residue's direction and the other alleles sit near OR 1.

## Values in this rendering

| quantity | value |
|---|---|
| carriers: DQB1 57 V (all cases): tested carrier alleles | DQB1*05:01, DQB1*06:04 |
| carriers: DQB1 57 V (all cases): other tested alleles | 10 |
| carriers: DQB1 87 Y (DVT no): tested carrier alleles | DQB1*05:01, DQB1*05:02, DQB1*05:03, DQB1*06:04 |
| carriers: DQB1 87 Y (DVT no): other tested alleles | 8 |
| conditioning: all cases residue step 0: drawn tests 3/3 / 2/3 / 1/3 at GEC | 1 / 0 / 0 |
| conditioning: all cases residue step 1: drawn tests 3/3 / 2/3 / 1/3 at GEC | 0 / 0 / 0 |
| conditioning: DVT no residue step 0: drawn tests 3/3 / 2/3 / 1/3 at GEC | 0 / 0 / 34 |
| conditioning: DVT no residue step 1: drawn tests 3/3 / 2/3 / 1/3 at GEC | 0 / 0 / 0 |
| conditioning: DVT no allele step 0: drawn tests 3/3 / 2/3 / 1/3 at GEC | 0 / 0 / 1 |
| conditioning: DVT no allele step 1: drawn tests 3/3 / 2/3 / 1/3 at GEC | 0 / 0 / 0 |

## Full statistics

**Every step of every walk (random model)**

| case_set | class | step | conditioned_on | lead | entered | reason | gec | n_sig_gec_any | p_narrow | p_intermediate | p_full |
|---|---|---|---|---|---|---|---|---|---|---|---|
| all cases | allele | 0 | - | A*31:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.003143 | 0.003684 | 0.002306 |
| all cases | residue | 0 | - | DQB1:57:V | 1 | entered | 3/3 | 1 | 4.570e-05 | 1.847e-04 | 1.741e-04 |
| all cases | residue | 1 | DQB1:57:V | DPB1:65:I | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.01445 | 0.01161 | 0.02384 |
| DVT yes | allele | 0 | - | C*03:04 | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.001569 | 0.003032 | 0.005926 |
| DVT yes | residue | 0 | - | B:114:D | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.004402 | 0.003551 | 0.001722 |
| DVT no | allele | 0 | - | DPB1*05:01 | 1 | entered | 1/3 | 1 | 2.815e-04 | 8.435e-04 | 0.001415 |
| DVT no | allele | 1 | DPB1*05:01 | DQB1*04:02 | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.001818 | 0.002593 | 0.001209 |
| DVT no | residue | 0 | - | DQB1:87:Y | 1 | entered | 1/3 | 51 | 1.190e-04 | 3.290e-04 | 3.149e-04 |
| DVT no | residue | 1 | DQB1:87:Y | A:99:Y | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.004708 | 0.006592 | 0.006294 |
| DVT yes vs no, cases only | allele | 0 | - | DQB1*05:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 0.008098 | 0.004794 | 0.007008 |
| DVT yes vs no, cases only | residue | 0 | - | DQB1:-27:A | 0 | no marker significant at GEC in any cohort | 0/3 | 0 | 3.592e-04 | 4.607e-04 | 0.001131 |

**Alleles of each residue signal's gene**

| signal | case set | allele | carries the residue | OR N / I / F | P N / I / F | control frequency (full) |
|---|---|---|---|---|---|---|
| DQB1 57 V | all cases | DQB1*03:01 | 0 | 0.97 / 0.92 / 0.93 | 8.2e-01 / 5.2e-01 / 5.5e-01 | 0.106 |
| DQB1 57 V | all cases | DQB1*03:02 | 0 | 1.34 / 1.32 / 1.29 | 3.6e-02 / 4.0e-02 / 5.6e-02 | 0.0879 |
| DQB1 57 V | all cases | DQB1*03:03 | 0 | 1.10 / 1.08 / 1.03 | 4.2e-01 / 4.9e-01 / 7.5e-01 | 0.1475 |
| DQB1 57 V | all cases | DQB1*04:01 | 0 | 1.08 / 1.09 / 1.12 | 5.4e-01 / 4.6e-01 / 3.3e-01 | 0.1225 |
| DQB1 57 V | all cases | DQB1*04:02 | 0 | 1.52 / 1.52 / 1.61 | 5.5e-02 / 4.5e-02 / 2.1e-02 | 0.0331 |
| DQB1 57 V | all cases | DQB1*05:01 | 1 | 0.61 / 0.65 / 0.66 | 3.7e-03 / 1.0e-02 / 9.6e-03 | 0.0695 |
| DQB1 57 V | all cases | DQB1*05:02 | 0 | 1.60 / 1.53 / 1.65 | 8.7e-02 / 1.0e-01 / 4.8e-02 | 0.0224 |
| DQB1 57 V | all cases | DQB1*05:03 | 0 | 0.86 / 0.81 / 0.78 | 4.8e-01 / 3.0e-01 / 2.1e-01 | 0.0413 |
| DQB1 57 V | all cases | DQB1*06:01 | 0 | 0.93 / 0.94 / 0.96 | 5.1e-01 / 5.5e-01 / 6.5e-01 | 0.1952 |
| DQB1 57 V | all cases | DQB1*06:02 | 0 | 1.12 / 1.16 / 1.14 | 4.8e-01 / 3.5e-01 / 4.0e-01 | 0.0658 |
| DQB1 57 V | all cases | DQB1*06:03 | 0 | 1.11 / 1.19 / - | 8.5e-01 / 7.3e-01 / - | — |
| DQB1 57 V | all cases | DQB1*06:04 | 1 | 0.74 / 0.76 / 0.77 | 6.6e-02 / 7.9e-02 / 8.6e-02 | 0.0737 |
| DQB1 87 Y | DVT no | DQB1*03:01 | 0 | 1.32 / 1.22 / 1.22 | 2.1e-01 / 3.5e-01 / 3.3e-01 | 0.106 |
| DQB1 87 Y | DVT no | DQB1*03:02 | 0 | 1.30 / 1.28 / 1.30 | 2.4e-01 / 2.6e-01 / 2.3e-01 | 0.0879 |
| DQB1 87 Y | DVT no | DQB1*03:03 | 0 | 1.23 / 1.23 / 1.17 | 2.5e-01 / 2.5e-01 / 3.7e-01 | 0.1475 |
| DQB1 87 Y | DVT no | DQB1*04:01 | 0 | 1.06 / 1.07 / 1.06 | 7.7e-01 / 7.3e-01 / 7.4e-01 | 0.1225 |
| DQB1 87 Y | DVT no | DQB1*04:02 | 0 | 2.81 / 2.64 / 2.83 | 2.8e-03 / 3.9e-03 / 1.7e-03 | 0.0331 |
| DQB1 87 Y | DVT no | DQB1*05:01 | 1 | 0.42 / 0.43 / 0.45 | 8.4e-04 / 1.4e-03 / 1.8e-03 | 0.0695 |
| DQB1 87 Y | DVT no | DQB1*05:02 | 1 | 1.23 / 1.18 / 1.35 | 6.5e-01 / 7.1e-01 / 4.7e-01 | 0.0224 |
| DQB1 87 Y | DVT no | DQB1*05:03 | 1 | 0.51 / 0.55 / 0.54 | 5.2e-02 / 6.9e-02 / 5.6e-02 | 0.0413 |
| DQB1 87 Y | DVT no | DQB1*06:01 | 0 | 0.84 / 0.83 / 0.86 | 3.1e-01 / 2.5e-01 / 3.7e-01 | 0.1952 |
| DQB1 87 Y | DVT no | DQB1*06:02 | 0 | 1.16 / 1.20 / 1.14 | 5.7e-01 / 4.9e-01 / 5.9e-01 | 0.0658 |
| DQB1 87 Y | DVT no | DQB1*06:03 | 0 | 1.62 / 1.94 / - | 5.3e-01 / 4.0e-01 / - | — |
| DQB1 87 Y | DVT no | DQB1*06:04 | 1 | 0.76 / 0.80 / 0.78 | 2.8e-01 / 3.8e-01 / 3.0e-01 | 0.0737 |

## How to read it

1. conditioning.png: compare left and right in each row; the coloured bars should turn grey.
2. carriers.png: do the carriers share the residue's direction, and do the non-carriers sit near 1?

## What this figure does *not* establish

- Statistical, not functional, separation: markers in near-perfect LD with the lead cannot be told apart from it.
- Untested alleles: only alleles that passed the typing QC (MAF ≥ 0.5 %) are shown.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
