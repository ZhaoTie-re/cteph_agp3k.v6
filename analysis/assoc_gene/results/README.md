# assoc_gene results

| layer | where |
|---|---|
| deliverable: robust genes | `_comparison/robust_genes/summary_table.md` (+ `summary_table.<cohort>.tsv`, `tiers.tsv`, `evidence.tsv`, `variant_detail.tsv`, `multi_carrier_detail.tsv`) |
| main figure | `figures/04.robust_genes/robust_genes.png` |
| per-gene figures (Tier 1 and 2) | `figures/05.genes/<GENE>.png`, `index.tsv` |
| other figure families | `figures/01.variant_sets/`, `figures/02.gene_scan/`, `figures/03.calibration/` — one `README.md` each |
| cross-cohort tables | `_comparison/tables/` |
| per-cohort stages | `<cohort>/00.prep … 05.signals` |
| run record | `_run_info/run_manifest.json`, `trace.txt`, `report.html`, `timeline.html`, `dag.html` |

Cohorts: narrow_mainland, intermediate_mainland, full_mainland. Variant sets: low_moderate_high, moderate_high. Statistics: cmc, skato. Report cohort: full_mainland. Column dictionaries: `docs/OUTPUTS.md`; the figure system: `docs/FIGURES.md`.
