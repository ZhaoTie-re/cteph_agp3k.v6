# Conditional analyses: independent HLA signals and the GWAS MHC signal

**Covers:** `ld_lead_hla.png`, `lead_given_hla.png`, `hla_given_lead.png`; tables in this directory (start with `signals_all.tsv` and `consistency.tsv`)

## The question this figure answers

Within each marker class, how many independent signals are there, and in how many cohorts does each hold; and is the genome-wide scan's MHC signal (chr6:32632861:A:G) the same signal as an HLA allele or amino-acid residue?

## Panels

**(ld_lead_hla.png) LD of the lead with every HLA marker**

r² between the lead's ALT-allele dosage and each allele / residue dosage, case set all, one row per cohort; markers in the gene order of the region figures.

**(lead_given_hla.png) the lead conditioned on HLA markers**

Rows case sets, columns cohorts. Each point is the lead's −log10 P after adding ONE HLA marker (the 10 in highest LD with it, the 10 with the smallest marginal P, and the independent-signal leads) to the model; x is that marker's |r| with the lead. The dashed line is the lead's marginal −log10 P in the same samples and model. A point far below the line is a marker that explains the lead.

**(hla_given_lead.png) HLA markers conditioned on the lead (random model)**

Each point is one HLA marker: marginal −log10 P (x) against −log10 P with the lead in the model (y). Points that drop below the identity line lose their signal to the lead; points that stay on it carry signal the lead does not.

## Interpretation

Random model, all cases: residue signals DQB1:57:V (3/3); allele signals - (consistency at the LD-adjusted threshold, conditional on the earlier signals). The GWAS lead's marginal P is 7.0e-06/1.8e-06/1.5e-06 (narrow / intermediate / full); given the first residue signal's lead it is 6.2e-03/9.8e-04/8.0e-04, and that lead given the GWAS lead is 4.8e-02/1.5e-01/1.4e-01. The single marker that explains most of the GWAS lead is DQB1:57:V/DQB1:57:V/DQB1:57:V; its best LD tag (DRB1:10:Q/DRB1:-25:K/DRB1:-16:A) leaves it at 1.5e-06/7.0e-06/3.4e-06. Signals are selected jointly over the three cohorts: at each step the candidate markers are those significant at the LD-adjusted (GEC) threshold in at least one cohort, and the lead is the candidate whose largest P over the cohorts is smallest. The cohorts share almost all cases, so 3/3 is a sample-selection sensitivity analysis, not replication.

## Full statistics

**Independent signals (all case sets, both models; random is primary)**

| case_set | model | class | signal | lead | gec | maxt | bonf | gec_marginal | p_narrow | p_intermediate | p_full | or_narrow | or_intermediate | or_full | n_members | carrier_alleles |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | random | residue | 1 | DQB1:57:V | 3/3 | 3/3 | 1/3 | 3/3 | 4.570e-05 | 1.847e-04 | 1.741e-04 | 0.6205 | 0.6521 | 0.6559 | 1 | DQB1*05:01, DQB1*06:04 |
| all | fixed | allele | 1 | B*40:02 | 1/3 | 1/3 | 0/3 | 1/3 | 3.743e-04 | 0.001778 | 0.007445 | 1.684 | 1.534 | 1.425 | 1 | — |
| all | fixed | residue | 1 | DQB1:57:V | 3/3 | 3/3 | 1/3 | 3/3 | 3.540e-05 | 1.482e-04 | 1.076e-04 | 0.591 | 0.6262 | 0.6273 | 3 | DQB1*05:01, DQB1*06:04 |
| dvt_no | random | allele | 1 | DPB1*05:01 | 1/3 | 1/3 | 1/3 | 1/3 | 2.815e-04 | 8.435e-04 | 0.001415 | 1.629 | 1.554 | 1.507 | 1 | — |
| dvt_no | random | residue | 1 | DQB1:87:Y | 1/3 | 1/3 | 0/3 | 1/3 | 1.190e-04 | 3.290e-04 | 3.149e-04 | 0.5427 | 0.5692 | 0.5774 | 51 | DQB1*05:01, DQB1*05:02, DQB1*05:03, DQB1*06:04 |
| dvt_no | fixed | allele | 1 | DPB1*05:01 | 1/3 | 1/3 | 1/3 | 1/3 | 1.908e-04 | 6.024e-04 | 0.001031 | 1.61 | 1.538 | 1.496 | 1 | — |
| dvt_no | fixed | residue | 1 | DQB1:87:Y | 1/3 | 1/3 | 0/3 | 1/3 | 1.267e-04 | 3.231e-04 | 3.104e-04 | 0.4768 | 0.5096 | 0.5184 | 46 | DQB1*05:01, DQB1*05:02, DQB1*05:03, DQB1*06:04 |

**Per case set × model (values narrow / intermediate / full)**

| case_set | model | allele_signals | residue_signals | gwas_marginal_p<1e-5_all | gwas_marginal_p | gwas_p_given_residue_signal1 | gwas_residue_signal1_p_given_lead | gwas_most_explaining | gwas_p_given_most_explaining | gwas_best_ld_tag | gwas_p_given_best_ld_tag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| all | random | - | DQB1:57:V (3/3) | 1 | 7.0e-06/1.8e-06/1.5e-06 | 6.2e-03/9.8e-04/8.0e-04 | 4.8e-02/1.5e-01/1.4e-01 | DQB1:57:V/DQB1:57:V/DQB1:57:V | 6.2e-03/9.8e-04/8.0e-04 | DRB1:10:Q/DRB1:-25:K/DRB1:-16:A | 1.5e-06/7.0e-06/3.4e-06 |
| all | fixed | B*40:02 (1/3) | DQB1:57:V (3/3) | 1 | 3.6e-06/8.1e-07/4.3e-07 | 7.6e-03/1.0e-03/6.0e-04 | 2.2e-02/9.6e-02/1.1e-01 | DQB1:57:V/DQA1:175:E/DQA1:175:E | 7.6e-03/1.1e-03/6.3e-04 | DRB1:10:Q/DRB1:-25:K/DRB1:-16:A | 8.9e-07/1.6e-06/1.3e-06 |
| dvt_yes | random | - | - | 0 | 1.4e-03/3.8e-04/3.5e-04 | -/-/- | -/-/- | DQA1:175:E/DRB1:-25:R/DQA1:175:E | 3.9e-02/1.7e-02/1.3e-02 | DQA1:175:E/DRB1:-25:R/DRB1:-25:R | 3.9e-02/1.7e-02/8.7e-03 |
| dvt_yes | fixed | - | - | 0 | 1.3e-03/3.3e-04/2.7e-04 | -/-/- | -/-/- | DQA1:175:E/DQA1:175:E/DQA1:175:E | 6.3e-02/2.2e-02/1.8e-02 | DQA1:175:E/DRB1:-25:R/DRB1:-25:R | 6.3e-02/4.0e-03/2.1e-03 |
| dvt_no | random | DPB1*05:01 (1/3) | DQB1:87:Y (1/3) | 0 | 3.9e-04/4.5e-04/5.1e-04 | 1.3e-01/9.5e-02/1.0e-01 | 3.3e-02/6.7e-02/5.9e-02 | DQB1:87:Y/DQB1:87:Y/DQB1:87:Y | 1.3e-01/9.5e-02/1.0e-01 | DQA1:175:E/DRB1:-25:K/DRB1:-25:R | 4.3e-02/4.4e-06/6.7e-06 |
| dvt_no | fixed | DPB1*05:01 (1/3) | DQB1:87:Y (1/3) | 0 | 3.5e-04/3.8e-04/4.6e-04 | 1.9e-01/1.3e-01/1.4e-01 | 1.0e-02/2.8e-02/2.9e-02 | DQB1:87:Y/DQB1:87:Y/DQB1:87:Y | 1.9e-01/1.3e-01/1.4e-01 | DQA1:175:E/DRB1:-25:K/DRB1:-25:R | 5.5e-02/1.2e-05/3.2e-05 |
| dvt_case_only | random | - | - | 0 | 5.5e-01/7.4e-01/8.6e-01 | -/-/- | -/-/- | DQB1:116:I/DQB1:116:V/DQB1:116:V | 4.7e-01/3.0e-01/3.0e-01 | DRB1:10:Q/DRB1:12:K/DRB1:149:H | 6.6e-02/6.5e-02/1.2e-01 |
| dvt_case_only | fixed | - | - | 0 | 5.1e-01/7.1e-01/8.2e-01 | -/-/- | -/-/- | DQB1:116:I/DQB1:116:I/DQB1:116:I | 4.5e-01/2.8e-01/2.8e-01 | DRB1:10:Q/DRB1:12:K/DRB1:149:H | 9.4e-02/9.3e-02/2.5e-01 |

**GWAS lead against the HLA markers, per case set × cohort × model**

| case_set | cohort | model | residue_signal1 | gwas_marginal_p | gwas_p_given_residue_signal1 | gwas_p_given_all_signals | gwas_most_explaining | gwas_p_given_most_explaining | gwas_best_ld_tag | gwas_best_ld_tag_r2 | gwas_p_given_best_ld_tag | gwas_residue_signal1_p_given_lead |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | narrow | random | DQB1:57:V | 7.027e-06 | 0.006237 | 0.006237 | DQB1:57:V | 0.006237 | DRB1:10:Q | 0.4691 | 1.490e-06 | 0.0479 |
| all | narrow | fixed | DQB1:57:V | 3.577e-06 | 0.007596 | 0.006635 | DQB1:57:V | 0.007596 | DRB1:10:Q | 0.4691 | 8.908e-07 | 0.02248 |
| all | intermediate | random | DQB1:57:V | 1.839e-06 | 9.773e-04 | 9.773e-04 | DQB1:57:V | 9.773e-04 | DRB1:-25:K | 0.4863 | 7.033e-06 | 0.1478 |
| all | intermediate | fixed | DQB1:57:V | 8.131e-07 | 0.001024 | 0.001232 | DQA1:175:E | 0.001104 | DRB1:-25:K | 0.4863 | 1.558e-06 | 0.09559 |
| all | full | random | DQB1:57:V | 1.455e-06 | 7.968e-04 | 7.968e-04 | DQB1:57:V | 7.968e-04 | DRB1:-16:A | 0.4934 | 3.404e-06 | 0.1442 |
| all | full | fixed | DQB1:57:V | 4.258e-07 | 5.979e-04 | 5.966e-04 | DQA1:175:E | 6.253e-04 | DRB1:-16:A | 0.4934 | 1.265e-06 | 0.1057 |
| dvt_yes | narrow | random |  | 0.001416 | — | — | DQA1:175:E | 0.0386 | DQA1:175:E | 0.4747 | 0.0386 | — |
| dvt_yes | narrow | fixed |  | 0.00129 | — | — | DQA1:175:E | 0.063 | DQA1:175:E | 0.4747 | 0.063 | — |
| dvt_yes | intermediate | random |  | 3.794e-04 | — | — | DRB1:-25:R | 0.01716 | DRB1:-25:R | 0.4811 | 0.01716 | — |
| dvt_yes | intermediate | fixed |  | 3.335e-04 | — | — | DQA1:175:E | 0.02206 | DRB1:-25:R | 0.4811 | 0.004036 | — |
| dvt_yes | full | random |  | 3.528e-04 | — | — | DQA1:175:E | 0.01294 | DRB1:-25:R | 0.4891 | 0.008658 | — |
| dvt_yes | full | fixed |  | 2.731e-04 | — | — | DQA1:175:E | 0.01784 | DRB1:-25:R | 0.4891 | 0.002058 | — |
| dvt_no | narrow | random | DQB1:87:Y | 3.946e-04 | 0.1297 | 0.1245 | DQB1:87:Y | 0.1297 | DQA1:175:E | 0.4684 | 0.04276 | 0.03298 |
| dvt_no | narrow | fixed | DQB1:87:Y | 3.465e-04 | 0.1879 | 0.1728 | DQB1:87:Y | 0.1879 | DQA1:175:E | 0.4684 | 0.05496 | 0.01015 |
| dvt_no | intermediate | random | DQB1:87:Y | 4.457e-04 | 0.0946 | 0.0876 | DQB1:87:Y | 0.0946 | DRB1:-25:K | 0.4858 | 4.417e-06 | 0.06672 |
| dvt_no | intermediate | fixed | DQB1:87:Y | 3.820e-04 | 0.1321 | 0.1158 | DQB1:87:Y | 0.1321 | DRB1:-25:K | 0.4858 | 1.205e-05 | 0.02796 |
| dvt_no | full | random | DQB1:87:Y | 5.096e-04 | 0.1029 | 0.09245 | DQB1:87:Y | 0.1029 | DRB1:-25:R | 0.4939 | 6.655e-06 | 0.05918 |
| dvt_no | full | fixed | DQB1:87:Y | 4.565e-04 | 0.136 | 0.1193 | DQB1:87:Y | 0.136 | DRB1:-25:R | 0.4939 | 3.236e-05 | 0.02889 |
| dvt_case_only | narrow | random |  | 0.5495 | — | — | DQB1:116:I | 0.4707 | DRB1:10:Q | 0.491 | 0.06609 | — |
| dvt_case_only | narrow | fixed |  | 0.512 | — | — | DQB1:116:I | 0.4453 | DRB1:10:Q | 0.491 | 0.09433 | — |
| dvt_case_only | intermediate | random |  | 0.7432 | — | — | DQB1:116:V | 0.3046 | DRB1:12:K | 0.4964 | 0.06531 | — |
| dvt_case_only | intermediate | fixed |  | 0.7065 | — | — | DQB1:116:I | 0.2822 | DRB1:12:K | 0.4964 | 0.0929 | — |
| dvt_case_only | full | random |  | 0.8582 | — | — | DQB1:116:V | 0.3 | DRB1:149:H | 0.4974 | 0.1233 | — |
| dvt_case_only | full | fixed |  | 0.8238 | — | — | DQB1:116:I | 0.2819 | DRB1:149:H | 0.4974 | 0.2476 | — |

**Every stepwise step**

| case_set | model | class | step | conditioned_on | lead | entered | reason | gec | p_narrow | p_intermediate | p_full | n_sig_gec_any |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | random | allele | 0 | - | A*31:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.003143 | 0.003684 | 0.002306 | 0 |
| all | random | residue | 0 | - | DQB1:57:V | 1 | entered | 3/3 | 4.570e-05 | 1.847e-04 | 1.741e-04 | 1 |
| all | random | residue | 1 | DQB1:57:V | DPB1:65:I | 0 | no marker significant at GEC in any cohort | 0/3 | 0.01445 | 0.01161 | 0.02384 | 0 |
| all | fixed | allele | 0 | - | B*40:02 | 1 | entered | 1/3 | 3.743e-04 | 0.001778 | 0.007445 | 1 |
| all | fixed | allele | 1 | B*40:02 | A*31:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.001067 | 0.001463 | 7.449e-04 | 0 |
| all | fixed | residue | 0 | - | DQB1:57:V | 1 | entered | 3/3 | 3.540e-05 | 1.482e-04 | 1.076e-04 | 3 |
| all | fixed | residue | 1 | DQB1:57:V | DPB1:65:I | 0 | no marker significant at GEC in any cohort | 0/3 | 0.01512 | 0.01144 | 0.02054 | 0 |
| dvt_yes | random | allele | 0 | - | C*03:04 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.001569 | 0.003032 | 0.005926 | 0 |
| dvt_yes | random | residue | 0 | - | B:114:D | 0 | no marker significant at GEC in any cohort | 0/3 | 0.004402 | 0.003551 | 0.001722 | 0 |
| dvt_yes | fixed | allele | 0 | - | C*03:04 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.001568 | 0.003107 | 0.004921 | 0 |
| dvt_yes | fixed | residue | 0 | - | B:114:D | 0 | no marker significant at GEC in any cohort | 0/3 | 0.003206 | 0.002638 | 0.001279 | 0 |
| dvt_no | random | allele | 0 | - | DPB1*05:01 | 1 | entered | 1/3 | 2.815e-04 | 8.435e-04 | 0.001415 | 1 |
| dvt_no | random | allele | 1 | DPB1*05:01 | DQB1*04:02 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.001818 | 0.002593 | 0.001209 | 0 |
| dvt_no | random | residue | 0 | - | DQB1:87:Y | 1 | entered | 1/3 | 1.190e-04 | 3.290e-04 | 3.149e-04 | 51 |
| dvt_no | random | residue | 1 | DQB1:87:Y | A:99:Y | 0 | no marker significant at GEC in any cohort | 0/3 | 0.004708 | 0.006592 | 0.006294 | 0 |
| dvt_no | fixed | allele | 0 | - | DPB1*05:01 | 1 | entered | 1/3 | 1.908e-04 | 6.024e-04 | 0.001031 | 1 |
| dvt_no | fixed | allele | 1 | DPB1*05:01 | DQB1*04:02 | 0 | no marker significant at GEC in any cohort | 0/3 | 9.897e-04 | 0.001311 | 7.483e-04 | 0 |
| dvt_no | fixed | residue | 0 | - | DQB1:87:Y | 1 | entered | 1/3 | 1.267e-04 | 3.231e-04 | 3.104e-04 | 46 |
| dvt_no | fixed | residue | 1 | DQB1:87:Y | A:99:Y | 0 | no marker significant at GEC in any cohort | 0/3 | 0.003358 | 0.004843 | 0.005224 | 0 |
| dvt_case_only | random | allele | 0 | - | DQB1*05:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.008098 | 0.004794 | 0.007008 | 0 |
| dvt_case_only | random | residue | 0 | - | DQB1:-27:A | 0 | no marker significant at GEC in any cohort | 0/3 | 3.592e-04 | 4.607e-04 | 0.001131 | 0 |
| dvt_case_only | fixed | allele | 0 | - | DQB1*05:01 | 0 | no marker significant at GEC in any cohort | 0/3 | 0.00994 | 0.006312 | 0.00845 | 0 |
| dvt_case_only | fixed | residue | 0 | - | DQB1:-27:A | 0 | no marker significant at GEC in any cohort | 0/3 | 4.452e-04 | 5.460e-04 | 0.001246 | 0 |

## How to read it

1. Start with the signals table: each row is one independent signal, its consistency (n/3) at the three thresholds, and the alleles that carry a residue signal.
2. lead_given_hla.png: which HLA marker pulls the GWAS lead down, in every cohort?
3. hla_given_lead.png: does the HLA signal survive the GWAS lead (points on the line)?

## What this figure does *not* establish

- Causality: conditioning separates statistically, not functionally; markers in near-perfect LD cannot be told apart (collinear rows have no conditional P), so a signal's lead stands for its members.
- Independence of the cohorts: they share almost all cases.
- Platform: all cases DNBSeq / NovaSeq, all controls HiSeq X; conditioning does not remove a platform effect shared by the lead and the HLA markers.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
