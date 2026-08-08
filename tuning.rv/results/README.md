# tuning.rv — results

Rare-variant depth/platform-confounding QC + a null-calibrated minAC recommendation.
This does **not** confirm associations (the design is fully confounded: cases 30x,
controls 15x, no platform overlap) — it quantifies and mitigates the confounding.

## Where to start
1. **Per cohort** → `<cohort>/README.md` (params, in→out counts, recommended minAC,
   β_g at the recommendation) and `<cohort>/figures/` (the four figures).
2. **Across cohorts** → `_comparison/cohort_comparison.png` + `.tsv`.
3. **Methods & data dictionary** → `../docs/METHODS.md`, `../docs/OUTPUTS.md`.

## Cohorts
`narrow_mainland` ⊂ `intermediate_mainland` ⊂ `full_mainland` — nested sample-selection variants.

## Per-cohort layout
- `figures/` — the four figures gathered in one place (sample QC, rate distribution,
  variant QC, minAC sweep). Each is also in its stage dir.
- `00.precheck/metrics/` — base metrics (permanent cache; the large `variant_metrics.txt.gz`).
- `00.precheck/01_sample_qc/`, `02_variant_qc/` — detection outputs.
- `01.sample_filter/`, `02.callset_filter/` — filtered genotypes (symlinks into work/).
- `03.qc_metrics/minac*/`, `04.qc_summary/minac*/` — per-minAC provenance.
- `05.qc_collect/` — merged stats, `minac_recommendation.tsv`, `minac_sweep.png`.

## Shared
- `_comparison/` — cross-cohort figure + table. · `_run_info/` — Nextflow trace / report / timeline / dag.

## Reading the sweep
**β_g** (Case−Control) is the DECISION axis; **β_d** (depth) is a DIAGNOSTIC only —
do NOT pick minAC to minimise it. Choose the minAC at the β_g non-significant trough;
confirm downstream with the genome-wide inflation factor λ (PC-adjusted rv test).
