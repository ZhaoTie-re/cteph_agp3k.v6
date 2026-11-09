# assoc_hla.v2 — outputs

Every path below is under `results/`. `<case_set>` is `all`, `dvt_yes`, `dvt_no` or `dvt_case_only`, and each has the same layout.

In `dvt_case_only` there are no controls. `PHENO1` is 1 for cases with DVT and 0 for cases without, so every `*_case` column below means **cases with DVT** and every `*_control` column means **cases without DVT**.

## Layout

| path | content |
|---|---|
| `_reference/hla_gene_coords.tsv` | GRCh38 span and chromosomal order of the 16 typed genes, from RefSeq |
| `_reference/case_sets.tsv` | `IID`, `dvt` (`yes` / `no` / `unknown`) and `dvt_raw`, from the v6 sample sheet |
| `<case_set>/<cohort>/00.prep/` | `pheno_covar.fixed.tsv`, `pheno_covar.random.tsv`, `fixed.keep`, `pheno_coding.json` (with `cases_dropped_by_dvt`); `null_model.json`, `null_model.samples.tsv`, `null_check.json`. `all` only: `validation.json`. DVT sets only: `null_model/` (`null.rda`, `null.varianceRatio.txt`, `null.fit.log` from `FIT_NULL`) |
| `<case_set>/<cohort>/01.markers/<class>/` | `<class>.bed/.bim/.fam` (the case set's samples; the set both engines read), `<class>.log`, `marker_info.tsv`; `effective_tests.<class>.tsv` and `tests.<class>.tsv` (the tests and three thresholds, below), `gec.<class>.log` and `gec.<class>.block.txt` (GEC v0.2's own output); `maxt_boot.<class>.tsv`, `.minp.tsv.gz` and `.observed.tsv` (the benchmark, below) |
| `<case_set>/<cohort>/02.assoc/fixed/<class>/` | plink2 output `fixed.<class>.PHENO1.glm.logistic` and `fixed.<class>.log` |
| `<case_set>/<cohort>/02.assoc/random/<class>/` | SAIGE output `random.<class>.saige.txt` and `random.<class>.log` |
| `<case_set>/<cohort>/03.tables/` | `assoc.<model>.<class>.tsv` and `scan_qc.<model>.<class>.tsv`, in the component's own schema |
| `<case_set>/<cohort>/04.conditional/pooled/` | `pooled.bed/.bim/.fam` (SENTINEL + alleles + residues + GWAS lead[s]), `pooled.log`, `pooled_markers.tsv` (marker, class, position, `file_order`, `r_<lead>`, `r_group0_<lead>`) |
| `<case_set>/<cohort>/04.conditional/<model>/` | `hla_given_lead.tsv`, `lead_given_hla.tsv`, `conditional.json` (below) |
| `<case_set>/_conditional/<model>/` | the independent-signal walk of each class, the three cohorts jointly: `steps.<class>.tsv`, `stepwise.<class>.tsv`, `signals.<class>.tsv`, `leads.<class>.txt` (below) |
| `<case_set>/_comparison/` | the case set's 12 scans stacked and compared (below) |
| `_case_sets/` | the case sets side by side (below) |
| `_conditional/` | every conditional analysis stacked and read across cohorts (below) |
| `figures/0N.<family>/` | the eight figure families, PNGs plus one `README.md` each (see `FIGURES.md`) |
| `_run_info/` | `run_manifest.json`, `trace.txt`, `report.html`, `timeline.html`, `dag.html` |

## `assoc.<model>.<class>.tsv` and `_comparison/assoc_all.tsv`

There is one row per tested marker and scan. Counts are those of the samples the model analysed: the `fixed` or the `random` sample set.

| column | meaning |
|---|---|
| `case_set`, `cohort`, `class`, `model` | scan identity; `class` is `allele` or `residue`; `model` is `fixed` or `random` |
| `marker` | matrix column name: `A*02:01` (allele) or `DQB1:57:V` (gene : IMGT position : residue) |
| `gene`, `position`, `residue` | parsed from the name; `position` and `residue` are blank for alleles |
| `chrom`, `bp`, `rank` | PLINK placement: chromosome 6, gene start + rank within gene. A label, not a physical position |
| `n_case`, `n_control` | samples in the model's set |
| `called_case`, `called_control` | samples with a determined call for this marker |
| `copies_case`, `copies_control` | copies carrying the marker among called samples |
| `freq_case`, `freq_control`, `freq` | carrier-copy frequency = copies / (2 × called) |
| `mac` | minor allele count over called copies (≥ 20 in every scan) |
| `missing_case`, `missing_control` | 1 − called / n (≤ 4.7 %) |
| `obs` | samples entering the test. Fixed: complete case for this marker. Random: all samples, missing calls mean-imputed |
| `or`, `l95`, `u95` | odds ratio per copy carrying the marker, and its 95 % CI |
| `beta`, `se` | fixed: log OR and its Wald SE. Random: SAIGE Tstat / var and 1 / √var |
| `p` | fixed: Wald P. Random: SAIGE `p.value.NA` (unadjusted score P). **Primary P** |
| `p_spa` | random only: SAIGE `p.value` (saddle-point); blank for fixed |
| `errcode` | plink2 `ERRCODE` (`.` = none); `NO_RESULT` if SAIGE returned no P; `MONOMORPHIC` if the marker has no minor allele in the analysed samples, so no engine tested it (case-only set only) |
| `test_id` | the test the marker belongs to: itself, or for the second residue of a mirror pair the first (blank if the marker is not a test) |
| `mirror_of` | residues only: the other residue of a mirror pair in this case set and cohort |
| `sig_gec`, `sig_maxt`, `sig_bonf` | `p` below the LD-adjusted (GEC, primary), exact (max-T) and conventional Bonferroni threshold of this (case set, cohort, class); `False` without a P |

## `scan_qc.<model>.<class>.tsv` and `_comparison/scan_qc_all.tsv`

There is one row per scan.

| column | meaning |
|---|---|
| `samples`, `cases`, `controls` | the model's sample set, from `pheno_coding.json` |
| `markers_tested` | markers in the matrix, equal to the typing export's tested count |
| `positions` | residues only: distinct (gene, IMGT position) |
| `with_p`, `errcodes` | markers with a P; error codes as `CODE:n` (`.` = none) |
| `lambda_gc` | median χ²(1) / 0.455 over the scan; reported, never applied |
| `tests`, `mirror_pairs` | tests (mirror pairs counted once) and mirror pairs |
| `me_gec` | effective number of independent tests from GEC v0.2 |
| `p_gec`, `p_maxt`, `p_bonf` | the three thresholds (from `effective_tests.<class>.tsv`) |
| `n_gec`, `n_maxt`, `n_bonf` | markers significant at each threshold |
| `tests_gec`, `tests_maxt`, `tests_bonf` | tests significant at each threshold (a mirror pair once) |
| `min_p`, `top_marker`, `top_or` | the scan's smallest P and its marker |
| `engine_cases`, `engine_controls` | what the engine reports it analysed. `COLLECT` stops the run if these differ from `cases` / `controls` |
| `engine_covariates` | plink2: covariates loaded, from its log. Random: the null model's covariate count |
| `impute_method` | random: `mean`, read from the SAIGE log |
| `variance_ratio` | random: as SAIGE reports it |
| `n_spa` | random: markers where SAIGE also computed the saddle-point P |

## `_comparison/model_concordance.tsv`

There is one row per cohort × class, comparing the fixed and random models.

| column | meaning |
|---|---|
| `markers_both` | markers with a P in both models |
| `pearson_neglog10p`, `spearman_p`, `pearson_beta` | agreement of −log10 P, of the P ranks and of β |
| `median_ratio_neglog10p` | median of −log10 P(random) / −log10 P(fixed) |
| `gec_fixed`, `gec_random`, `gec_both` | markers significant at the LD-adjusted (GEC) threshold per model, and in both |

## `<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv` and `tests.<class>.tsv`

`effective_tests` has one row: the tests and the three thresholds of the (case set, cohort, class), fixed from the genotypes before any P (`EFFECTIVE_TESTS`, METHODS §5).

| column | meaning |
|---|---|
| `alpha` | family-wise error rate of all three thresholds (0.05) |
| `markers`, `polymorphic`, `monomorphic` | markers in the PLINK set; those varying among the case set's random-model samples, and those not |
| `mirror_pairs`, `tests` | residue mirror pairs; tests = polymorphic markers with each mirror pair counted once |
| `gec_samples`, `gec_markers`, `gec_added` | GEC input: the case set's fixed-model (unrelated) samples; tests entered; tests monomorphic there, added to Me as independent |
| `me_gec_raw`, `me_gec`, `gec_blocks` | GEC v0.2's effective number; plus `gec_added`; its blocks |
| `p_gec` | **primary**: 0.05 / `me_gec` |
| `maxt_samples`, `maxt_draws`, `maxt_seed` | max-T simulation: random-model samples, draws, seed |
| `p_maxt`, `p_maxt_lo`, `p_maxt_hi`, `me_maxt` | exact family-wise threshold, its 95 % Monte Carlo interval, and 0.05 / `p_maxt` |
| `p_bonf` | conventional Bonferroni: 0.05 / `tests` |
| `fwer_gec`, `fwer_bonf` | family-wise error rate of `p_gec` and `p_bonf` under the max-T simulation |
| `gec_version`, `gec_sha256` | the GEC release and its jar's sha256 |

`tests` has one row per marker: `polymorphic` (0/1), `test_id`, `mirror_of`, `in_gec` (0/1: entered GEC; the second residue of a mirror pair does not).

## Omnibus tests of residue positions (a main analysis): `<case_set>/_omnibus/` and `_omnibus/`

Per case set (`OMNIBUS`, METHODS §8), residues only (no gene-level tests; `class` is always `residue`):

| file | content |
|---|---|
| `omnibus.residue.tsv` | marginal scan: every residue position (`unit` = `gene:position`) × cohort with `reference`, `categories` (tested), `pooled` (the categories in `rare`), `n`, `cases`, `controls`, `df`, `chisq`, `p`, `threshold` (Bonferroni), `sig_bonf`, `ors` (each category's OR against the reference), `errcode` (`.`, `SEPARATION`, `UNFINISHED`, `COLLINEAR`, `NOT_TESTED`) |
| `omnibus_stepwise.residue.tsv` | every position × cohort × step of the stepwise walk |
| `omnibus_steps.residue.tsv` | one row per step: the lead, whether it entered, its P per cohort, consistency |
| `omnibus_signals.residue.tsv` | one row per entered position: consistency at entry and marginally, P and df per cohort, members |
| `omnibus_given.residue.tsv` | every position × cohort × conditioning marker (`OMNIBUS_GIVEN`, METHODS §7.5): `conditioned_on` (a GWAS lead or a residue signal lead), `why` (`gwas_lead`, `signal_lead`), `unit`, `df`, `chisq`, `p` (given the marker), `p_marginal` (= `omnibus.residue.tsv`), `threshold` (Bonferroni, as the marginal scan), `sig_bonf`, `missing_imputed` (the marker's missing calls, mean-imputed), `errcode` |

Across case sets (`OMNIBUS_SUMMARY`), `_omnibus/`:

| file | content |
|---|---|
| `omnibus_all.tsv`, `omnibus_steps_all.tsv` | the per-case-set tables stacked |
| `omnibus_units.tsv` | one row per (case set, class, position): `df_`, `p_`, `sig_<cohort>`, `bonf` (n/3); `best_fixed`, `best_random` (the position's single residue whose largest P over the cohorts is smallest, per model) with their P per cohort; `random_signal_leads` (random-model single-marker signal leads at the position) |
| `omnibus_signals_all.tsv` | the signals stacked, with `random_signal_leads_in_unit` and `random_signal_leads_in_members` |

## `<case_set>/<cohort>/01.markers/<class>/maxt_boot.<class>.*` (benchmark)

The exact family-wise threshold by parametric bootstrap on the case set's fixed-model samples (`MAXT_BOOT`, METHODS §5.2). No call depends on it.

`maxt_boot.<class>.tsv`, one row:

| column | meaning |
|---|---|
| `samples`, `cases`, `controls`, `tests` | the fixed-model samples and the tests (one per mirror pair) |
| `draws`, `seed`, `null_model` | the replicates and the null logistic model they are drawn from |
| `cases_mean`, `cases_sd` | cases per replicate |
| `p_gec`, `p_bonf` | the scan's thresholds, from `effective_tests.<class>.tsv` |
| `p_maxt_<s>` [`_lo`, `_hi`], `me_maxt_<s>` | the exact threshold (0.05 quantile of the smallest P), 95 % interval from 1,000 resamples, and 0.05 / threshold; `s` = `wald` (plink2, the fixed model's test) or `score` (score test without SPA, blanks mean-imputed) |
| `fwer_gec_<s>`, `fwer_bonf_<s>` [`_lo`, `_hi`] | the share of replicates whose smallest P is below the threshold, Wilson 95 % interval |
| `replicates_with_na_<s>`, `tests_ok_min_<s>` | replicates in which a test returned no P (plink2: non-convergence or separation), and the fewest tests with a P |
| `p_maxt_mvn`, `fwer_gec_mvn`, `fwer_bonf_mvn` | the multivariate-normal max-T of `effective_tests.<class>.tsv`, for comparison |
| `observed_min_p_wald`, `observed_min_p_score` | the smallest P of the observed phenotype rerun |

`maxt_boot.<class>.minp.tsv.gz` holds one row per replicate (`cases`, `tests_ok_<s>`, `min_p_<s>`); `maxt_boot.<class>.observed.tsv` the observed phenotype's P per test (`p_wald`, equal to the fixed model's; `p_score`).

## `_comparison/marker_summary.tsv` and `hits.tsv`

There is one row per (class, marker), with the 6 cells (3 cohorts × 2 models) in the order narrow / intermediate / full mainland, fixed / random. Consistency is read in the random model, the primary model. `hits.tsv` keeps the rows significant under any threshold in at least one cohort (random model), or at the LD-adjusted threshold in at least one cohort in the fixed model.

| column | meaning |
|---|---|
| `p_<cohort>_<model>`, `or_<cohort>_<model>` | per cell; blank if the cohort did not test the marker |
| `mirror_of` | residues only: the other residue of a mirror pair in all three cohorts (one test with reciprocal ORs) |
| `cohorts_tested` | cohorts in which the marker was tested |
| `gec`, `maxt`, `bonf` | cross-cohort consistency at each threshold, random model: `n/3` = significant in n of the 3 cohorts with the same direction of effect (1/3 = "Single cohort") |
| `gec_n`, `maxt_n`, `bonf_n` | the same n as a number |
| `gec_cohorts`, `maxt_cohorts`, `bonf_cohorts` | the cohorts counted (`-` if none) |
| `gec_direction_conflict`, … | significant cohorts disagree in direction (the count is then the majority's) |
| `direction` | `risk` (OR > 1) or `protective` in the cohorts significant at GEC; blank if none |
| `gec_fixed`, `gec_n_fixed` | consistency at the LD-adjusted threshold in the fixed model, for model concordance |
| `min_p`, `max_p` | smallest and largest random-model P over the cohorts |
| `same_direction` | every random-model OR is on the same side of 1 |

## `_case_sets/`

| file | content |
|---|---|
| `scan_qc_case_sets.tsv` | the 48 scan summaries (4 case sets × 3 cohorts × 2 classes × 2 models), with the `scan_qc` columns above |
| `dvt_contrast.tsv` | one row per (class, cohort, model, marker), with columns `<x>_all`, `<x>_dvt_yes` and `<x>_dvt_no` for `x` in `n_case`, `freq_case`, `or`, `l95`, `u95`, `p`, `sig_gec`, `sig_maxt` and `sig_bonf`, and the same for `<x>_dvt_case_only`, which is the case-only test (`n_case` / `freq_case` there = DVT yes, plus `freq_control_dvt_case_only` = DVT no); also `gene`, `position`, `residue`, `freq_control` (the controls are shared), `called_dvt` (significant at any threshold in the DVT-yes or DVT-no set), `called_case_only` (significant at any threshold in the case-only test) and `min_p_dvt` (smallest P over the DVT-yes and DVT-no sets) |
| `dvt_hits.tsv` | every `dvt_contrast` row of a marker significant at any threshold in the DVT-yes set, the DVT-no set or the case-only test, in at least one cohort or model |

## Conditional analyses

### Independent signals: `<case_set>/_conditional/<model>/`

One walk per (case set, model, class), the three cohorts jointly (`STEPWISE`, METHODS §7.4). Columns `p_<cohort>` / `or_<cohort>` use the short cohort names (narrow, intermediate, full).

| file | one row per | columns |
|---|---|---|
| `steps.<class>.tsv` | step | `step`, `conditioned_on` (`;`-joined), `lead`, `entered`, `reason`, `lead_max_p` (the lead's largest P over the cohorts), `n_sig_gec_any` (markers significant at GEC in some cohort, not in the model), `n_sig_gec_<cohort>`, `p_<cohort>` / `or_<cohort>` of the lead, `gec` / `maxt` / `bonf` (the lead's consistency at that step) |
| `stepwise.<class>.tsv` | step × cohort × marker | `p`, `or`, `l95`, `u95`, `errcode` (`CONDITIONED` for markers in the model; `COLLINEAR` or a plink2 code when determined by them; `NO_RESULT` if the engine skipped a monomorphic marker), `sig_gec` |
| `signals.<class>.tsv` | entered signal | `signal`, `lead`, `step`, `conditioned_on`, `gec` / `maxt` / `bonf` (consistency at entry, conditional on the earlier signals), `gec_marginal`, `p_<cohort>` / `or_<cohort>` (at entry), `p_marginal_<cohort>` / `or_marginal_<cohort>`, `n_members` / `members` (markers significant before the lead entered and in no cohort after, the lead included), `carrier_alleles` (residue signals: tested alleles of the gene carrying the residue) |
| `leads.<class>.txt` | entered signal | the lead, one per line (read by `CONDITIONAL`) |

### GWAS lead against the HLA markers: `<case_set>/<cohort>/04.conditional/<model>/`

| file | one row per | columns |
|---|---|---|
| `hla_given_lead.tsv` | lead × HLA marker | `r` with the lead, `p_marginal`, `or_marginal`, `p_given_lead`, `or_given_lead`, `errcode` |
| `lead_given_hla.tsv` | lead × conditioning set | `conditioned_on`, `why` (`ld_tag`, `top_p`, `signal_lead`, `signal_leads_jointly`, or `position_jointly:` with its source, `residue_lead` and/or `omnibus_lead`), `position` and `position_reference` (whole-position sets only: the position, and the most frequent residue left out), `r`, `lead_p_marginal`, `lead_or_marginal`, `lead_p_given`, `lead_or_given`, `marker_p_marginal`, `marker_p_given_lead` |
| `conditional.json` | run | HLA markers, engine runs, the signal leads used, the GWAS lead(s) |
| `hla_given_signal.tsv` | conditioning set × marker (HLA markers and the GWAS lead) | `HLA_GIVEN_SIGNAL` (METHODS §7.5): `conditioned_on` (a signal lead, or a whole position's columns joined by `;`), `why` (`signal_lead` or `position_jointly:…`), `position`, `marker`, `class` (`allele`, `residue`, `gwas_lead`), `p_marginal`, `or_marginal`, `p_given`, `or_given`, `errcode` (`CONDITIONED` for the conditioning markers themselves); header only when the case set and model have no signal lead |

### Across runs: `_conditional/`

| file | content |
|---|---|
| `gwas_leads.tsv` | every MHC peak of the genome-wide scan: per-cohort P, OR, A1 and A1 frequency, `status` (`kept` / `excluded`) and `reason` |
| `signals_all.tsv`, `steps_all.tsv` | the per-walk tables stacked, each with `case_set`, `model` and `class` |
| `hla_given_lead_all.tsv`, `lead_given_hla_all.tsv` | the per-run tables stacked, each with `case_set`, `cohort` and `model` |
| `ld_lead_hla.tsv` | lead × cohort × HLA marker: `r`, `r2`, `r_controls` (case set `all`) |
| `summary.tsv` | one row per (case set, cohort, model): the first residue and allele signal leads; then, prefixed with the GWAS lead's ID, its marginal P; its P given each first signal lead and given all signal leads; its P given the whole position of the first residue signal lead (`position1`, `p_given_position1`) and of the first omnibus position lead (`omnibus_position1`, `p_given_omnibus_position1`); the HLA marker that explains it most and the P given it; its best LD tag, that tag's r² and the P given it; the HLA minimum P given the lead; and the first residue signal lead's P given the GWAS lead |
| `consistency.tsv` | one row per (case set, model): the signals of each class with their consistency, and the GWAS-lead quantities written narrow / intermediate / full |
