# tuning.rv — Outputs & data dictionary

Single source of truth for the output tree and every column emitted. (Replaces
the former per-cohort `METRICS_DESCRIPTION.md` and per-minAC `methods.md`.)

---

## Output tree

```
results/
├── README.md                              # TOP-LEVEL index — where to start; what each dir holds
├── <cohort>/                              # one per --Cohorts entry
│   ├── README.md                          # run manifest (params, in→out, recommended minAC)
│   ├── figures/                           # the four publication PNGs (600 dpi), gathered in one place
│   ├── 00.precheck/
│   │   ├── metrics/                        # base metrics on RAW set (permanent storeDir cache)
│   │   │   ├── sample_metrics.txt.gz
│   │   │   ├── variant_metrics.txt.gz(.tbi)
│   │   │   └── calc_metrics.log
│   │   ├── 01_sample_qc/
│   │   │   ├── sample_outliers.remove.txt          # PLINK2 --remove (FID<TAB>IID)
│   │   │   ├── sample_qc_rate_residual.png         # rate-residual vs depth (group/depth/platform)
│   │   │   ├── rate_distribution.png               # rate before/after QC + by depth/platform/group
│   │   │   └── *.log
│   │   └── 02_variant_qc/
│   │       ├── variant_metrics.txt.gz(.tbi)        # RECOMPUTED on sample-cleaned set (symlink)
│   │       ├── variants_depthdiff.exclude.txt      # PLINK2 --exclude (VariantID)
│   │       ├── protect_variants.txt                # candidate variants never excluded
│   │       ├── depthdiff_variants_cdf.png
│   │       ├── depthdiff_stats.tsv
│   │       └── *.log
│   ├── 01.sample_filter/                   # sample-cleaned genotype (symlink) + filter_samples.log
│   ├── 02.callset_filter/                  # fully-filtered genotype (symlink) + filter_variants.log
│   ├── 03.qc_metrics/minac*/               # per-minAC metrics (symlink) + calc_metrics.log
│   ├── 04.qc_summary/minac*/               # per-minAC stats TSV + log
│   └── 05.qc_collect/
│       ├── all_qc_summary_stats.tsv        # merged per-minAC stats (+ Cohort column)
│       ├── minac_recommendation.tsv        # the recommended minAC + rationale
│       ├── minac_sweep.png                 # HEADLINE figure
│       └── cohort_summary.tsv              # one-row per-cohort summary
├── _comparison/                           # cross-cohort (cohort_comparison.png/.tsv, decision_axis.png/.tsv, all_cohorts_*.tsv)
└── _run_info/                             # Nextflow trace / report / timeline / dag
```

Bulky genotype and metric files are **symlinked** into `results/` (they live in
`work/`); text records and figures are **copied** (durable).

---

## Data dictionaries

### `sample_metrics.txt.gz` (one row per sample)
| Column | Meaning |
| --- | --- |
| `SampleID` | sample ID (`ID_JHRPv6`) |
| `Group` | `Case` (PH) / `Control` |
| `TargetDP` | target depth (15x / 30x) |
| `MeanDP` | observed mean depth `D_mean` |
| `SNumHomRef` / `SNumHet` / `SNumHomAlt` | genotype counts (callable = their sum) |
| `SMissCount` | missing genotypes |
| `SMinAC` | minor-allele burden = `SNumHet + 2·SNumHomAlt` |

Derived: `rate = SMinAC / (SNumHomRef+SNumHet+SNumHomAlt)`, reported ×1000.

### `variant_metrics.txt.gz` (one row per variant)
| Column | Meaning |
| --- | --- |
| `#CHROM`, `POS`, `VariantID`, `REF`, `ALT` | variant identity (`VariantID` = `chr<N>:<pos>:<ref>:<alt>`) |
| `RefAC` / `AltAC` / `MinAC` | allele counts (all samples) |
| `*_15X` / `*_30X` | the same, within the 15x / 30x depth strata |
| `VMISS_15X` / `VMISS_30X` | missing rate within each stratum (drives variant QC) |
| `VNumHomRef` / `VNumHet` / `VNumHomAlt` / `VMissCount` / `VMISS` | per-variant genotype tallies + overall missing rate |

### `all_qc_summary_stats.tsv` (one row per minAC)
| Column | Meaning |
| --- | --- |
| `Cohort`, `MinAC_Threshold` | keys |
| `Variant_Count`, `Sample_Count` | retained counts at this minAC |
| `Rate_Mean` | mean burden rate (×1000) |
| `Beta_group_noPC`, `_L`, `_U`, `P_group_noPC` | OLS `β_g` (decision), 95% CI, p — no-PC primary model (`rate ~ group + depth`) |
| `Beta_depth_noPC`, `_L`, `_U`, `P_depth_noPC` | OLS `β_d` (depth diagnostic) — no-PC model |
| `Beta_group_pc`, `_L`, `_U`, `P_group_pc` | OLS `β_g` — +ancestry-PC sensitivity model |
| `Beta_depth_pc`, `_L`, `_U`, `P_depth_pc` | OLS `β_d` — +PC model |
| `PC_source`, `N_PC` | ancestry PC label + count (e.g. `bbj_mainland`, 10) |
| `SE_type` | `classical` or `HC3` |
| `MWU_P_TargetDP`, `Wasserstein_Dist_TargetDP` | depth-distribution divergence across target depth |

### `minac_recommendation.tsv` (key/value)
| Key | Meaning |
| --- | --- |
| `recommended_minac` | the recommended threshold `κ` (V-trough) |
| `preferred_model` | `noPC` (primary) or `pc` |
| `calibrated` | `true` if a non-significant trough was found, else `false` (fell back to `argmin|β_g|`) |
| `rule` | human-readable selection rule |
| `alpha`, `min_floor`, `pc_source` | rule parameters + ancestry PC label |
| `trough_minac`, `nonsig_band` | `argmin|β_g|` and the list of non-significant minACs |
| `beta_group_at_rec`, `ci_low_at_rec`, `ci_high_at_rec`, `p_group_at_rec` | `β_g` at the recommendation |
| `variant_count_at_rec` | variants retained at the recommendation (power) |
| `sensitivity_model`, `sensitivity_rec`, `sensitivity_nonsig_band` | the +PC (ancestry) sensitivity model's recommendation |
| `models_agree_at_rec` | whether the +PC model is also non-significant at the recommended minAC |

### `depthdiff_stats.tsv` (key/value)
`knee`, `n_total`, `n_valid`, `n_removed`, `removed_pct`, `retained_pct`,
`n_protected_requested`, `n_protected_present`, `n_protected_rescued`.

### `cohort_summary.tsv` / `_comparison/cohort_comparison.tsv`
Per-cohort: samples removed / variants excluded / retained counts, depth-diff
`knee` + `retained_pct`, `recommended_minac` + `calibrated`, `β_group` / `p_group` /
variant count at the recommendation, and the +PC sensitivity (`pc_sensitivity_rec` /
`models_agree_at_rec`).

### `_comparison/decision_axis.tsv` (section / key / value)
The numbers behind `decision_axis.png`, so the model-choice argument in METHODS §5 is
auditable rather than asserted. Three sections:

| Section | Keys |
| --- | --- |
| `premise` | `platforms_total`, `platforms_shared_by_both_groups` (0 = zero overlap) |
| `beta_choice` | per cohort: `argmin_abs_beta_group`, `argmin_abs_beta_depth`, `beta_group_at_argmin_beta_depth`, `p_group_at_argmin_beta_depth`, `amplification_if_selected_on_beta_depth` |
| `model_choice` | per cohort: `beta_group_noPC_at_k2`, `beta_group_pc_at_k2`, `shift_at_k2`, `shift_over_abs_beta`, `sign_flipped_at_k2`, `r2_case_given_<N>_pcs` |

---

## Figures

Every figure is a **600-dpi PNG**, self-contained and overlap-free, styled after a
print journal. Bold panel letters (**a**, **b**, …) sit in the title slot at each
panel's top-left, so they never overlap the data, ticks, or axis labels. A thin
rule separates the panels from a **Nature/Science-style caption**: a **bold
declarative title** (the finding), then a bold-letter run-in describing each panel
(**a**, …; **b**, …), a plain sentence of objective interpretation (探讨), and a
bold **Data** (source, marker/CI definitions, statistics) and **Model** (estimator /
equation in professional notation) footer. Colour keys sit in an empty panel corner
or above the panels — never over data. They are copied into each cohort's `figures/`
folder.

| File | What it shows |
| --- | --- |
| `sample_qc_rate_residual.png` | depth-adjusted rate residual per sample; outliers ringed; coloured by group / depth / platform |
| `rate_distribution.png` | rate before vs after sample QC (skew/kurtosis) + rate by depth / platform / group |
| `depthdiff_variants_cdf.png` | signed & absolute ΔVMISS distributions + the Δ CDF with the Kneedle knee |
| `minac_sweep.png` | **headline** — 2×2: `β_g` full range, `β_g` trough zoom, `β_d` (depth diagnostic), power cost; recommended minAC marked |
| `cohort_comparison.png` | 2×3 — rows are the two models, columns are `β_g` / `β_d` / (power cost, decision summary). Panel **f** gives `β_g` ± 95 % CI at each cohort's own recommended minAC under both models |
| `_comparison/decision_axis.png` | why the decision is read from `β_g` in the no-PC model: platform × case/control overlap, the two coefficients' differing minima, and what PC adjustment does to the reading |
