# tuning.rv — run summary (intermediate_mainland)

Rare-variant depth/platform-confounding QC + minAC sweep. This run QUANTIFIES and
MITIGATES the confounding and recommends a minAC; it does not confirm associations
(the design is fully confounded: cases 30x, controls 15x, no platform overlap).
Methodology: `../../docs/METHODS.md` · Output dictionary: `../../docs/OUTPUTS.md`.

## Resolved parameters
| Parameter | Value |
| --- | --- |
| Cohort | `intermediate_mainland` |
| Rare call set | `/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/intermediate_mainland/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold` |
| Info file | `/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx` |
| Ancestry PC (+PC model) | `/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/intermediate_mainland/covariates/bbj_mainland_pc.sex.tsv` (bbj_mainland, 10 PCs) |
| minAC sweep | 0-20 |
| OLS SE | classical |
| Robust-Z threshold | 5.0 |
| Kneedle (S, wx, wy) | 1.0, 1.0, 2.0 |
| minAC recommendation | V-trough; primary model = nopc (alpha 0.05, floor 1) |
| Conda env | cteph_geno_pro |

## Run
| Field | Value |
| --- | --- |
| Run name | gloomy_stallman |
| Started | 2026-08-08T16:26:30.636607705+09:00 |
| Command | `nextflow run tuning.rv.nf -resume` |
| Revision | n/a |
| Nextflow | 26.04.4 |
| Timing | see `../_run_info/report.html` |

## Protected variants (never excluded by variant QC)
- `chr4:76306990:G:C`
- `chr4:76309487:C:G`
- `chr4:76309559:G:T`

## Output tree
- `figures/` - the four publication figures gathered in one place (sample QC, rate distribution, variant QC, minAC sweep)
- `00.precheck/metrics/` - base metrics (permanent storeDir cache) + `calc_metrics.log`
- `00.precheck/01_sample_qc/` - `sample_outliers.remove.txt`, `sample_qc_rate_residual.png`, `rate_distribution.png`, logs
- `00.precheck/02_variant_qc/` - recomputed `variant_metrics.txt.gz` (symlink), `variants_depthdiff.exclude.txt`, `depthdiff_variants_cdf.png`, `depthdiff_stats.tsv`, logs
- `01.sample_filter/` - sample-cleaned genotype (symlink) + `filter_samples.log`
- `02.callset_filter/` - fully-filtered genotype (symlink) + `filter_variants.log`
- `03.qc_metrics/minac*/` - per-minAC metrics (symlink) + `calc_metrics.log`
- `04.qc_summary/minac*/` - per-minAC stats TSV + log
- `05.qc_collect/` - `all_qc_summary_stats.tsv`, `minac_recommendation.tsv`, `minac_sweep.png`, `cohort_summary.tsv`
- `../_comparison/` - cross-cohort comparison (figure + table)
- `../_run_info/` - Nextflow trace / report / timeline / dag

## Result (in → out, and the recommended minAC)

| Metric | Value |
| --- | --- |
| Samples removed (rate-residual outliers) | 3 |
| Variants excluded (depth-diff) | 976407 |
| Samples retained | 2477 |
| Variants retained (minAC 0) | 16508040 |
| Depth-diff knee | 0.03 (retained 94.4156%) |
| Protected variants (never excluded) | 3 present, 1 rescued from exclusion |
| **Recommended minAC** (V-trough) | **2**  (calibrated: true) |
| Primary model / non-sig band | noPC / {2} |
| beta_group @ rec (noPC) | 0.00224664  (p = 0.751042) |
| Variants kept @ rec | 7620757 |
| +PC sensitivity | recommends minAC 2; agrees at rec: true |
