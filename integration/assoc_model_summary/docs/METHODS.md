# assoc_model_summary — methods

## 1. Scans compared

| name here | upstream | model | covariates | samples |
|---|---|---|---|---|
| `fixed` | `assoc_plink2` `additive` | PLINK 2 logistic regression, Wald test, no Firth | SEX + PC1–10 | `fixed_model` genotypes: relatives (PI_HAT > 0.2) removed |
| `grm_pc` | `assoc_saige` `main` | SAIGE logistic mixed model, full GRM, LOCO, score test with the normal approximation (`p.value.NA`) | SEX + PC1–10 | `random_model` genotypes: relatives retained |
| `grm_only` | `assoc_saige` `detect` | as `grm_pc` | SEX | as `grm_pc` |

All three read the same variant set (same IDs, same order; MAF ≥ 0.01 in the fixed-model controls). The PCs are
the projections onto the BBJ mainland-Japanese PCA (`PC1_AVG`–`PC10_AVG`).

## 2. Scan reduction and calibration (`reduce_scan.py`)

- Usable variant: ERRCODE `.`, 0 < P ≤ 1, finite OR > 0, finite SE > 0. This is the rule of the upstream
  `call_peaks.py` / `plot_manhattan_qq.py`.
- λ<sub>GC</sub> = median(χ²₁ of P) / 0.4549364. It is computed overall and within study-wide MAF bins:
  < 0.01, 0.01–0.02, 0.02–0.05, 0.05–0.10, 0.10–0.20, 0.20–0.50.
  - The < 0.01 bin holds variants that pass the control-MAF ≥ 0.01 filter but are rarer overall. It is
    reported in the tables and not plotted.
- Counts of variants at P < 10⁻⁵ and < 5 × 10⁻⁸, against n × α.
- λ<sub>GC</sub>, the P < 10⁻⁵ / 5 × 10⁻⁸ counts and the number analysed are asserted equal to the upstream
  `scan_qc.tsv` (\|Δλ\| < 10⁻⁴).

## 3. Figures

- **Style.** The shared `plot_style.py` at its `paper` scale, one point larger: the figures are about 10 in
  wide, for a 16:9 slide. Palette, threshold lines, chromosome bands, marker sizes and the QQ band are those of
  `cohort_manhattan.png`. Panel letters only, no caption block.
- **Manhattan** (`plot_paired.py`).
  - One cumulative genomic axis and one −log₁₀P limit, max(8.4, 1.08 × top), shared by every panel.
  - Lines at 5 × 10⁻⁸ (red) and 10⁻⁵ (grey dashed).
  - Dotted guides at every locus genome-wide in any panel, named once above the top row.
  - Diamonds at a panel's own genome-wide leads, taken from the upstream peak calls.
  - Every variant at P ≤ 10⁻³ is drawn. Below that, one point per cell of a 1200 × 450 grid over the panel:
    about 2 px at 600 dpi against a marker about 15 px across. The image is unchanged, and the drawing is
    about 30 × faster.
- **QQ.**
  - The models are overlaid in their colours, each with its λ<sub>GC</sub>.
  - Pointwise 95% band: 2.5–97.5% quantiles of Beta(k, n − k + 1).
  - Drawn at the 20,000 smallest P plus a log-spaced sample of 30,000 ranks.
  - One limit for all rows.
- **Calibration** (`plot_calibration.py`):
  - **a**, λ<sub>GC</sub> per cohort;
  - **b**, λ<sub>GC</sub> by MAF bin from 0.01 (model = colour, cohort = marker);
  - **c**, variants at P < 10⁻⁵ against n × 10⁻⁵.
- **Main summary** (`plot_main_summary.py`).
  - **a**, cases / controls per cohort with λ<sub>GC</sub> and the genome-wide / suggestive locus counts.
  - **b**, OR (95% CI) of every locus genome-wide in any cohort (panel b of `lead_crosscohort.tsv`) or
    suggestive in all three (panel c), at its representative variant (`locus_rep == 1`). The three cohorts
    sit side by side in the slide colours; a diamond marks a genome-wide call in that cohort. The right
    column gives the smallest P in its cohort's colour.

## 4. Question-specific statistics

- **Fixed-model genome-wide leads** (`write_readme.py`): the random-model P is the smallest within ±250 kb of
  the lead. This makes the comparison per locus, because the two models may elect different lead variants of
  one locus.
- **PCs and case status** (`pc_assoc.py`):
  - logistic regression (`statsmodels` Logit) of case status on SEX vs SEX + PC1–10, on the mixed-model samples;
  - likelihood-ratio χ²(10) and its P, Wald P per PC, McFadden pseudo-R² gain.
- **Model-specific peaks** (`pc_effect_tables.py`):
  - a peak of one model with no peak of the other within 500 kb in the same cohort;
  - upstream peak calls (variants at P < 10⁻⁵ merged within ±250 kb).
  - The extended MHC is taken as chr6:28.4–33.5 Mb.
- **Effect stability**: Δ ln OR of the 16 cross-cohort leads × 3 cohorts, from `lead_bymodel.tsv`.

## 5. Text

`write_readme.py` writes every README from the run's tables and JSON; no number is typed in.
- Each README has a **Slide summary** (a take-home line and bullets) and a **Detailed description**
  (question, design, how to read, results, interpretation, limits, methods).
- `02.pc_effect` states facts and leaves its conclusion open on purpose.
- `results/README.md` repeats the three slide summaries.

## 6. Checks

- `reduce_scan.py` asserts agreement with the upstream `scan_qc.tsv` (§2).
- `verify.sh`:
  - syntax and site rules;
  - every script flag the `.nf` passes is accepted;
  - the upstream files are present;
  - task cardinality;
  - the 9 scans against the upstream `scan_qc_all.tsv`;
  - 5 PNGs and 4 READMEs, with slide summaries and no NaN; question 2 still marked pending.
