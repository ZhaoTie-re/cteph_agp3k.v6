# assoc_hla.v2 results

| layer | where |
|---|---|
| per case set (`all`, `dvt_yes`, `dvt_no`, `dvt_case_only`) | `<case_set>/` — same layout in each (below) |
| all results, one schema | `<case_set>/_comparison/assoc_all.tsv` |
| markers significant in any cohort, with n/3 per threshold | `<case_set>/_comparison/hits.tsv` (all markers: `marker_summary.tsv`) |
| tests and thresholds (GEC, max-T, Bonferroni) | `<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv`, `tests.<class>.tsv`, `gec.<class>.log` |
| benchmark: exact max-T by parametric bootstrap, fixed-model samples | `<case_set>/<cohort>/01.markers/<class>/maxt_boot.<class>.tsv` (+ `.minp.tsv.gz` per replicate, `.observed.tsv`) |
| scan summaries | `<case_set>/_comparison/scan_qc_all.tsv`; fixed vs random: `model_concordance.tsv` |
| figures | `figures/` — main `01.landscape/`, `02.signals/`, `03.conditioning/`, `04.gwas_link/`, `05.dvt/`; supplementary `06.thresholds/`, `07.models/`, `08.qq/`; one `README.md` each |
| per-cohort stages | `<case_set>/<cohort>/00.prep`, `01.markers/<class>`, `02.assoc/<model>/<class>`, `03.tables`, `04.conditional` |
| DVT yes against DVT no | `_case_sets/` — `dvt_contrast.tsv`, `dvt_hits.tsv`, `scan_qc_case_sets.tsv`; the case-only test itself is `dvt_case_only/`; figures `figures/02.signals/`, `05.dvt/` |
| conditional analyses | `_conditional/` — `signals_all.tsv` (independent signals), `consistency.tsv`, `summary.tsv`, `steps_all.tsv`, `hla_given_lead_all.tsv`, `lead_given_hla_all.tsv`, `ld_lead_hla.tsv`, `gwas_leads.tsv`; figures `figures/03.conditioning/`, `04.gwas_link/`; per walk: `<case_set>/_conditional/<model>/`; per GWAS-lead run: `<case_set>/<cohort>/04.conditional/<model>/` (`hla_given_lead.tsv`, `lead_given_hla.tsv`; every HLA marker given each signal lead and whole position: `hla_given_signal.tsv`) |
| omnibus tests of residue positions (a main analysis; fixed model, Bonferroni) | `_omnibus/` — `omnibus_units.tsv` (each unit beside its best single markers), `omnibus_signals_all.tsv`, `omnibus_all.tsv`, `omnibus_steps_all.tsv`; per case set: `<case_set>/_omnibus/` (positions given each GWAS lead and residue signal lead: `omnibus_given.residue.tsv`) |
| reference tables | `_reference/hla_gene_coords.tsv`, `_reference/case_sets.tsv` |
| run record | `_run_info/run_manifest.json`, `trace.txt`, `report.html`, `timeline.html`, `dag.html` |

Case sets: all, dvt_yes, dvt_no, dvt_case_only (the first three against all controls; dvt_case_only = DVT yes vs DVT no within cases). Cohorts: narrow_mainland, intermediate_mainland, full_mainland. Classes: allele, residue. Models: random (SAIGE, primary), fixed (plink2). Thresholds: LD-adjusted GEC (primary), exact max-T (validation), conventional Bonferroni (sensitivity); consistency = significant in n/3 cohorts. Column dictionaries: `../docs/OUTPUTS.md`; methods: `../docs/METHODS.md`.
