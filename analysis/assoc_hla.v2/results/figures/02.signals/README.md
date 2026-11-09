# The independent signals

**Main figure.** **Covers:** `signals.png` in `02.signals/`

## The question this figure answers

For each independent signal, how large is the effect in every case-set comparison and cohort, and at which thresholds is it significant?

## Panels

**(signals.png a) odds ratio in every comparison and cohort**

Rows grouped by signal (bold: the lead, and the case set it was selected in); within a group, one row per comparison: all cases, DVT yes and DVT no against all controls, and DVT yes against DVT no within cases. Three estimates per row, one per cohort (dark narrow → light full), OR per copy with 95 % CI, random model. Filled: significant at that comparison's LD-adjusted threshold.

**(signals.png b) significance at each threshold**

For the same rows: significance in each cohort (N, I, F) at the LD-adjusted (GEC), exact (max-T) and conventional Bonferroni thresholds of that comparison, cohort and class (filled = significant), and the number of cohorts with the same direction (n/3).

## Interpretation

3 independent signal(s) in the random model: DQB1 57 V (all cases, residue: GEC 3/3, max-T 3/3, Bonferroni 1/3); DQB1 87 Y (DVT no, residue: GEC 1/3, max-T 1/3, Bonferroni 0/3); DPB1*05:01 (DVT no, allele: GEC 1/3, max-T 1/3, Bonferroni 1/3). Each is shown in every comparison, so a signal selected in one case set can be read in the others. An effect of the same size and direction in DVT yes and DVT no, with a case-only OR near 1, is a shared effect that the smaller DVT sets lack the power to call; a case-only OR away from 1, significant, would be the evidence of a difference. The DVT-yes and DVT-no comparisons share their controls, so their side-by-side estimates are not a test.

## Values in this rendering

| quantity | value |
|---|---|
| signals: DQB1:57:V (all cases signal), all cases vs controls: OR N / I / F | 0.62 / 0.65 / 0.66 |
| signals: DQB1:57:V (all cases signal), DVT yes vs controls: OR N / I / F | 0.64 / 0.67 / 0.69 |
| signals: DQB1:57:V (all cases signal), DVT no vs controls: OR N / I / F | 0.52 / 0.55 / 0.55 |
| signals: DQB1:57:V (all cases signal), DVT yes vs DVT no: OR N / I / F | 1.47 / 1.43 / 1.44 |
| signals: DQB1:87:Y (DVT no signal), all cases vs controls: OR N / I / F | 0.72 / 0.74 / 0.74 |
| signals: DQB1:87:Y (DVT no signal), DVT yes vs controls: OR N / I / F | 0.81 / 0.82 / 0.84 |
| signals: DQB1:87:Y (DVT no signal), DVT no vs controls: OR N / I / F | 0.54 / 0.57 / 0.58 |
| signals: DQB1:87:Y (DVT no signal), DVT yes vs DVT no: OR N / I / F | 1.91 / 1.76 / 1.70 |
| signals: DPB1*05:01 (DVT no signal), all cases vs controls: OR N / I / F | 1.27 / 1.20 / 1.17 |
| signals: DPB1*05:01 (DVT no signal), DVT yes vs controls: OR N / I / F | 1.28 / 1.16 / 1.11 |
| signals: DPB1*05:01 (DVT no signal), DVT no vs controls: OR N / I / F | 1.63 / 1.55 / 1.51 |
| signals: DPB1*05:01 (DVT no signal), DVT yes vs DVT no: OR N / I / F | 0.77 / 0.75 / 0.75 |

## Full statistics

**Independent signals (random model)**

| case_set | class | signal | lead | gec | maxt | bonf | n_members | carrier_alleles |
|---|---|---|---|---|---|---|---|---|
| all cases | residue | 1 | DQB1:57:V | 3/3 | 3/3 | 1/3 | 1 | DQB1*05:01, DQB1*06:04 |
| DVT no | residue | 1 | DQB1:87:Y | 1/3 | 1/3 | 0/3 | 51 | DQB1*05:01, DQB1*05:02, DQB1*05:03, DQB1*06:04 |
| DVT no | allele | 1 | DPB1*05:01 | 1/3 | 1/3 | 1/3 | 1 | — |

**Every signal in every comparison**

| signal | selected in | comparison | OR N / I / F | P N / I / F | gec | maxt | bonf |
|---|---|---|---|---|---|---|---|
| DQB1 57 V | all cases | all cases vs controls | 0.62 / 0.65 / 0.66 | 4.6e-05 / 1.8e-04 / 1.7e-04 | 3/3 | 3/3 | 1/3 |
| DQB1 57 V | all cases | DVT yes vs controls | 0.64 / 0.67 / 0.69 | 1.1e-02 / 2.0e-02 / 2.7e-02 | 0/3 | 0/3 | 0/3 |
| DQB1 57 V | all cases | DVT no vs controls | 0.52 / 0.55 / 0.55 | 3.6e-04 / 1.0e-03 / 7.1e-04 | 0/3 | 0/3 | 0/3 |
| DQB1 57 V | all cases | DVT yes vs DVT no | 1.47 / 1.43 / 1.44 | 2.2e-01 / 2.4e-01 / 2.2e-01 | 0/3 | 0/3 | 0/3 |
| DQB1 87 Y | DVT no | all cases vs controls | 0.72 / 0.74 / 0.74 | 1.2e-03 / 1.7e-03 / 1.9e-03 | 0/3 | 0/3 | 0/3 |
| DQB1 87 Y | DVT no | DVT yes vs controls | 0.81 / 0.82 / 0.84 | 1.6e-01 / 1.8e-01 / 2.2e-01 | 0/3 | 0/3 | 0/3 |
| DQB1 87 Y | DVT no | DVT no vs controls | 0.54 / 0.57 / 0.58 | 1.2e-04 / 3.3e-04 / 3.1e-04 | 1/3 | 1/3 | 0/3 |
| DQB1 87 Y | DVT no | DVT yes vs DVT no | 1.91 / 1.76 / 1.70 | 1.2e-02 / 2.3e-02 / 2.9e-02 | 0/3 | 0/3 | 0/3 |
| DPB1*05:01 | DVT no | all cases vs controls | 1.27 / 1.20 / 1.17 | 4.2e-03 / 2.3e-02 / 4.3e-02 | 0/3 | 0/3 | 0/3 |
| DPB1*05:01 | DVT no | DVT yes vs controls | 1.28 / 1.16 / 1.11 | 5.8e-02 / 2.3e-01 / 3.7e-01 | 0/3 | 0/3 | 0/3 |
| DPB1*05:01 | DVT no | DVT no vs controls | 1.63 / 1.55 / 1.51 | 2.8e-04 / 8.4e-04 / 1.4e-03 | 1/3 | 1/3 | 1/3 |
| DPB1*05:01 | DVT no | DVT yes vs DVT no | 0.77 / 0.75 / 0.75 | 1.4e-01 / 9.2e-02 / 8.8e-02 | 0/3 | 0/3 | 0/3 |

## How to read it

1. Read the first row of each group: the comparison the signal was selected in.
2. Read down the group: is the OR the same in DVT yes and DVT no, and is the case-only OR near 1?
3. Read panel b across: does the call hold at max-T and at conventional Bonferroni?

## What this figure does *not* establish

- Fine-mapping: the lead stands for the markers it explains (03.conditioning, members).
- Random-model ORs are SAIGE one-step estimates; read them for direction and rough size.
- Replication: the three cohorts are nested and share almost all cases; n/3 is a sample-selection sensitivity analysis, not independent replication.
- Platform: all cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; a typing difference between platforms would look like an association, and no covariate separates the two. Within cases the platform mix also differs by DVT status.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
