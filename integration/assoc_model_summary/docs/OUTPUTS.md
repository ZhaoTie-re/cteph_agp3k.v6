# assoc_model_summary — outputs

```
results/
  README.md                         index + the three slide summaries
  00.scans/                         <cohort>.<model>.stats.json, <cohort>.<model>.top.tsv   (model: fixed, grm_pc, grm_only)
  01.fixed_vs_random/               fixed_vs_random.png, calibration.png, README.md,
                                    tail_counts.tsv, lambda_by_maf.tsv, gw_leads.tsv
  02.pc_effect/                     grm_pc_vs_grm.png, README.md, pc_case_association.tsv, pc_effect.json,
                                    peak_changes.tsv, lead_or_stability.tsv, tail_counts.tsv, lambda_by_maf.tsv, gw_leads.tsv
  03.main_results/                  main_manhattan.png, main_summary.png, README.md, main_loci.tsv
  _run_info/                        run_manifest.json, trace.txt, report.html, timeline.html, dag.html
```

## Figures (≈ 10 in wide, 600 dpi, no caption block)

| file | panels |
|---|---|
| `fixed_vs_random.png` | rows a–c = cohorts; columns: fixed-model Manhattan, random-model Manhattan, QQ of both |
| `calibration.png` | a λ<sub>GC</sub> per cohort; b λ<sub>GC</sub> by MAF bin; c variants at P < 10⁻⁵ vs chance |
| `grm_pc_vs_grm.png` | as `fixed_vs_random.png`, random model without vs with PCs |
| `main_manhattan.png` | rows a–c = cohorts; Manhattan + QQ of the main analysis |
| `main_summary.png` | a cases / controls, λ<sub>GC</sub>, loci per cohort; b OR (95% CI) per locus × cohort |

## Tables

`00.scans/<cohort>.<model>.stats.json` — `n_all`, `n_excluded`, `n`, `lambda_gc`, `lambda_by_maf` (bin, lo, hi,
n, lambda_gc), `n_p_lt_suggestive`, `n_p_lt_genomewide`, their expectations, `min_p`, `n_case`, `n_ctrl`,
`n_eff`, `upstream_lambda_gc`.

`00.scans/<cohort>.<model>.top.tsv` — `ID CHROM POS A1_FREQ OR P` of every variant at P < 10⁻³.

`tail_counts.tsv` — per cohort × model:
- `n_case`, `n_ctrl`, `n_eff`, `n_variants`, `lambda_gc`;
- `n_p_lt_1e5`, `expected_p_lt_1e5`, `ratio_p_lt_1e5`, `n_p_lt_5e8`, `min_p`;
- `gw_loci`, `sg_loci` (upstream peaks) and `gw_leads` (labels).

`lambda_by_maf.tsv` — `cohort model maf_bin maf_lo maf_hi n lambda_gc`.

`gw_leads.tsv` — every genome-wide lead: `model cohort Gene rsID variant_id P OR`.

`pc_case_association.tsv` — per cohort:
- `n`, `n_case`, `n_ctrl`, `n_pcs`;
- `lr_chi2`, `lr_df`, `lr_p`;
- `pseudo_r2_sex`, `pseudo_r2_sex_pcs`;
- `n_pcs_p_lt_0_05`, `p_PC1` … `p_PC10`.

`peak_changes.tsv` — peaks of one model with no peak of the other within 500 kb:
- `cohort`, `only_in`, `only_in_label`;
- `gene`, `chrom`, `pos`, `variant_id`, `rsID`;
- `tier`, `n_sig_variants`, `P`, `OR`.

`lead_or_stability.tsv` — `lead_bymodel.tsv` of `assoc_saige` plus `dlnOR` (ln OR with PCs − without) and
`dnlp` (−log₁₀P with − without).

`pc_effect.json` — per cohort:
- λ<sub>GC</sub> and P < 10⁻⁵ counts with and without PCs;
- peaks per model and model-specific peaks;
- the PC / case-status test;
- OR shift summary.

`main_loci.tsv` — `group` (b = genome-wide in ≥ 1 cohort, c = suggestive in all three), `locus_id`, `gene`,
`variant_id`, `rsID`, `cohort`, `OR`, `L95`, `U95`, `P`, `called_peak`.
