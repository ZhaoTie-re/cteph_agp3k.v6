# 3 · Main analysis in the three cohorts

## Slide summary

**Take-home: the main analysis (random model: SAIGE, SEX + 10 PCs + full GRM) is well calibrated in all three nested cohorts (λGC 1.039 / 1.038 / 1.021). *FTO* (rs9934504): OR 1.67 / 1.68 / 1.61, P = 1.4 × 10⁻⁷ / 3.8 × 10⁻⁸ / 1.8 × 10⁻⁷; genome-wide in intermediate mainland. 15 further loci are suggestive in all three cohorts, each with the same direction of effect in every cohort.**

- Cohorts (cases / controls): Narrow mainland 419 / 1,776; Intermediate mainland 429 / 2,079; Full mainland 439 / 2,662.
- Loci, genome-wide + suggestive: narrow 0 + 29; intermediate 1 + 27; full 0 + 32.
- Suggestive in all three cohorts: *ARL13B*, *HS3ST3A1*, *LINC01146*, *GALNTL6*, *ARHGEF26-AS1*, *HLA-DQA1*, *INPP5D*, *MFSD11*, *CDH4*, *SLC44A2*, *PHF20L1*, *ABI3BP*, *ALPI*, *OR13C8*, *ADARB2*.
- The cohorts are nested and share every case: agreement between them shows robustness to the ancestry filter, not replication.

## Figures

- **`main_manhattan.png`** — one row per cohort: Manhattan and QQ of the main analysis on one genomic axis and one −log₁₀P limit (the layout of `cohort_manhattan.png`, without its caption).
- **`main_summary.png`** — (a) cases and controls per cohort with λGC and the number of genome-wide (GW) and suggestive (SG) loci; (b) OR (95% CI) of every locus genome-wide in any cohort or suggestive in all three, the cohorts side by side, with the smallest P in its cohort's colour.

## Detailed description

### Question
What does the main analysis show in the three nested cohorts?

### Design
SAIGE logistic mixed model, SEX + PC1–10 + full GRM, related samples retained; score test with the normal
approximation. Loci: upstream peak calls (variants at P < 10⁻⁵ merged within ±250 kb; genome-wide = lead
P < 5 × 10⁻⁸). The loci of panel b: genome-wide in at least one cohort, or suggestive in all three cohorts,
each at one representative variant shared by the three cohorts.

### How to read the figures
- **`main_manhattan.png`** — as `cohort_manhattan.png` of `analysis/assoc_saige`: rows are cohorts on one
  genomic axis and one −log₁₀P limit; red line 5 × 10⁻⁸, grey dashed 10⁻⁵; the dotted guide and the label
  mark the locus genome-wide in any cohort, the diamond the cohort that called it.
- **`main_summary.png`** — (a) composition and calibration per cohort; (b) one row per locus, the three cohorts
  side by side (diamond = genome-wide in that cohort), OR on a log scale with the line at 1; right, the
  smallest P of the three, in the colour of its cohort.

### Results
| Cohort | cases / controls | N_eff | variants | λGC | GW loci | SG loci | smallest P |
|---|---|---|---|---|---|---|---|
| Narrow mainland | 419 / 1,776 | 1,356 | 5,113,626 | 1.039 | 0 | 29 | 1.4 × 10⁻⁷ |
| Intermediate mainland | 429 / 2,079 | 1,422 | 5,113,982 | 1.038 | 1 | 27 | 3.8 × 10⁻⁸ |
| Full mainland | 439 / 2,662 | 1,507 | 5,115,116 | 1.021 | 0 | 32 | 6.4 × 10⁻⁸ |

| Group | Gene | Variant | rsID | OR (95% CI) narrow | OR (95% CI) intermediate | OR (95% CI) full | smallest P |
|---|---|---|---|---|---|---|---|
| GW in ≥1 | FTO | chr16:53882967:G:A | rs9934504 | 1.67 (1.37–2.04) | 1.68 (1.39–2.03) | 1.61 (1.34–1.94) | 3.8 × 10⁻⁸ (intermediate) |
| SG in all 3 | ARL13B | chr3:94020684:T:C | rs143483852 | 5.25 (2.59–10.65) | 6.08 (2.94–12.57) | 6.03 (2.91–12.51) | 6.4 × 10⁻⁸ (full) |
| SG in all 3 | HS3ST3A1 | chr17:13528059:G:A | rs34804183 | 2.31 (1.57–3.39) | 2.55 (1.74–3.74) | 2.53 (1.73–3.70) | 3.4 × 10⁻⁷ (full) |
| SG in all 3 | LINC01146 | chr14:88019710:G:A | rs4899932 | 1.50 (1.27–1.78) | 1.50 (1.28–1.77) | 1.43 (1.22–1.68) | 4.4 × 10⁻⁷ (intermediate) |
| SG in all 3 | GALNTL6 | chr4:171998307:T:C | rs75214171 | 4.39 (2.28–8.45) | 4.48 (2.37–8.49) | 4.17 (2.24–7.76) | 4.8 × 10⁻⁷ (intermediate) |
| SG in all 3 | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | 3.76 (2.08–6.81) | 3.90 (2.15–7.10) | 3.89 (2.17–6.96) | 5.7 × 10⁻⁷ (full) |
| SG in all 3 | HLA-DQA1 | chr6:32632861:A:G | rs9272130 | 0.67 (0.56–0.80) | 0.66 (0.56–0.79) | 0.67 (0.57–0.79) | 7.3 × 10⁻⁷ (full) |
| SG in all 3 | INPP5D | chr2:233068979:G:A | rs72982244 | 4.04 (2.11–7.72) | 4.10 (2.15–7.83) | 4.31 (2.26–8.22) | 1.0 × 10⁻⁶ (full) |
| SG in all 3 | MFSD11 | chr17:76739260:T:C | rs9897202 | 2.85 (1.77–4.57) | 2.98 (1.86–4.76) | 2.68 (1.70–4.23) | 1.1 × 10⁻⁶ (intermediate) |
| SG in all 3 | CDH4 | chr20:60915853:C:T | rs79313991 | 3.09 (1.90–5.03) | 2.91 (1.80–4.69) | 2.87 (1.78–4.62) | 1.3 × 10⁻⁶ (narrow) |
| SG in all 3 | SLC44A2 | chr19:10630195:G:C | rs4548995 | 1.50 (1.26–1.79) | 1.49 (1.26–1.77) | 1.49 (1.26–1.76) | 1.5 × 10⁻⁶ (full) |
| SG in all 3 | PHF20L1 | chr8:132794261:A:T | rs150804769 | 3.61 (2.00–6.53) | 3.78 (2.08–6.86) | 3.57 (2.00–6.39) | 1.9 × 10⁻⁶ (intermediate) |
| SG in all 3 | ABI3BP | chr3:100783010:T:C | rs7648469 | 1.91 (1.44–2.51) | 1.86 (1.42–2.43) | 1.83 (1.41–2.37) | 2.0 × 10⁻⁶ (narrow) |
| SG in all 3 | ALPI | chr2:232438125:T:C | rs790023 | 0.50 (0.37–0.69) | 0.50 (0.37–0.68) | 0.51 (0.38–0.69) | 2.4 × 10⁻⁶ (intermediate) |
| SG in all 3 | OR13C8 | chr9:104559024:C:T | rs7036847 | 1.62 (1.31–2.00) | 1.61 (1.31–1.97) | 1.54 (1.26–1.89) | 2.4 × 10⁻⁶ (intermediate) |
| SG in all 3 | ADARB2 | chr10:1747167:CGTG:C | — | 0.34 (0.20–0.56) | 0.35 (0.22–0.57) | 0.36 (0.22–0.59) | 5.0 × 10⁻⁶ (narrow) |

All 16 of 16 loci have the same direction of effect in the three cohorts.

### Interpretation
- Calibration is good in every cohort; the relaxed ancestry filter (narrow → full) adds controls without
  adding inflation.
- *FTO* is the one genome-wide locus; the suggestive loci shared by all three cohorts carry concordant
  effects of similar size.
- The cohorts are nested — narrow ⊂ intermediate ⊂ full — and share every case, so concordance is expected
  and is not independent replication.

### Methods
As `01.fixed_vs_random`. Loci and their representative variants: `analysis/assoc_saige`
`_comparison/tables/lead_crosscohort.tsv`; the per-variant details are in `integration/saige_signals`.
