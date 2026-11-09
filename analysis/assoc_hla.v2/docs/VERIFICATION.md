# assoc_hla.v2 — verification

## Gates inside the run

A gate that fails stops the run, and its reason is written to disk.

| where | gate |
|---|---|
| workflow start | Every cohort name is valid. `ImputeMethod` is `mean`. The marker counts are read from the typing export's `summary.json` |
| `CASE_SETS` | Every `DVT_Presence` value is `Yes`, `No`, `Unknown` or blank; anything else stops the run. Sample IDs are unique |
| `VALIDATE_INPUTS` | 30 checks per cohort, written to `all/<cohort>/00.prep/validation.json` (listed below) |
| `CHECK_NULL` | Every case set's null model, reused (`all`) or fitted (`FIT_NULL`, DVT sets), is binary and LOCO, has exactly the case set's samples and the same label for every sample, and its covariates are SEX, PC1_AVG … PC10_AVG. Written to `null_check.json`; `ASSOC_RANDOM` runs only on a checked null |
| `BUILD_MARKERS` | Every case-set sample is a matrix row, kept in phenotype-file order. Every gene has coordinates. Positions are unique and increasing |
| `EFFECTIVE_TESTS` | The GEC jar's sha256 equals `params.GecSha256` (`tools/gec-0.2/SOURCE.md`). plink2's exported dosages are the `.bim` markers in order. GEC reads and uses every input test (none dropped by a filter). The covariate file covers exactly the PLINK set's samples |
| `MAXT_BOOT` | `fixed.keep` and the fixed-model phenotype file hold the same samples, with no missing phenotype or covariate and a 0/1 phenotype. The null logistic model converges. The tests equal `effective_tests.<class>.tsv`'s count, and plink2 exports exactly those tests and samples. plink2 writes a result for every phenotype column, and every replicate has a P for at least one test |
| `OMNIBUS` | `fixed.keep` and the fixed-model phenotype file hold the same samples, with no missing phenotype or covariate and a 0/1 phenotype. The exported markers equal the `.bim`. Every unit's samples are blank in all of its columns or in none, and no sample carries more than 2 listed copies |
| `EXTRACT_LEAD` | The lead VCF holds exactly the kept lead(s) |
| `STEPWISE`, `CONDITIONAL` (`engine.py`) | Every run analyses the case set's cases and controls, and the effect allele is the ALT on every plink2 row. An engine's first SAIGE run is unconditional, and every SAIGE run with `--condition` reproduces every marker's marginal statistics from it, which guards against the SAIGE first-marker defect the sentinel absorbs. `STEPWISE`: every per-cohort input belongs to its case set, cohort and class |
| `COLLECT` | The engine analysed the expected cases and controls; plink2 loaded 11 covariates (plink2's own log); SAIGE ran with `mean` imputation and LOCO on chromosome 6 (SAIGE's own log). Every marker with a minor allele in the analysed samples has exactly one result, and the effect allele is `P`. Only a monomorphic marker may lack a result, and it keeps its row as `MONOMORPHIC`. The thresholds are those of the same case set, cohort and class, and the tests table lists the same markers |

The 30 `VALIDATE_INPUTS` checks are:

- the md5 of each matrix equals `summary.json`;
- the columns equal the `tested` rows of `markers.tsv`;
- the rows equal the phenotype samples, in order;
- matrix values are in {blank, 0, 1, 2};
- the phenotype is 0/1 with no missing values, and its counts match `summary.json`;
- all 11 covariates are present and numeric for every sample;
- the fixed set is a subset of the cohort;
- the reused SAIGE null model is binary and LOCO, has the same samples and the same case/control label for every sample, and its covariates are SEX, PC1_AVG … PC10_AVG;
- every case has a row in the sample sheet, and no control carries a DVT status.

## `verify.sh` (read-only, after the run)

`bash verify.sh` returns exit 0 only if every check passes. It recomputes from the inputs rather than reading the pipeline's summaries back. §2–§9 and §12 run for every case set × cohort.

| § | check |
|---|---|
| 0 | No file of this component refers to a path in `analysis/assoc_hla/`. No `cpus` or `memory` directive is used, and every `--rsc` asks for exactly cores × 4,571 MB. The sha256 of `tools/gec-0.2/gec.jar` equals `params.GecSha256` and the one in `SOURCE.md` |
| 1 | `validation.json`: 30 checks per cohort, none failed |
| 12 | Case set recomputed from `case_sets.tsv` and the phenotype file. The random-model samples are all controls plus the set's cases, or, for `dvt_case_only`, the cases of known status with `PHENO1` recoded (1 = DVT yes). Both the samples and `PHENO1` are compared in phenotype-file order, and the counts equal `pheno_coding.json` |
| 2 | Each PLINK `.bed`, decoded directly, equals the typing matrix rows of the case set cell by cell, blanks included. The `.fam` equals those rows. The bim alleles are `P` (column 5) and `A` (column 6) |
| 3 | plink2 log: the analysed cases and controls equal `pheno_coding.json`; 11 covariates; `--1`; `no-firth`. SAIGE: `N_case` and `N_ctrl` on every row, `mean` imputation, LOCO on chromosome 6, and the variance ratio equals the null model's |
| 4 | Every tested marker has exactly one result row per scan, including `MONOMORPHIC` rows |
| 5 | `marker_info.tsv` (n, called and copies per group, both sample sets) recomputed from the matrices |
| 6 | Per case set × cohort × class: the tests recomputed from the matrix (polymorphic markers, mirror pairs as complements in every sample) equal `tests.<class>.tsv` and the `tests` count. GEC v0.2, re-run here on an input rebuilt from the PLINK set, the fixed-model samples and the recomputed tests, reproduces Me (1e-3), and its log shows v0.2. The thresholds follow: `p_bonf` = 0.05 / tests, `p_gec` = 0.05 / Me with 1 ≤ Me ≤ tests, `p_maxt` inside its Monte Carlo interval and ≥ `p_bonf`. The GEC threshold's family-wise error under max-T lies within 0.035–0.065. Per scan: every `sig_gec` / `sig_maxt` / `sig_bonf` flag, the per-threshold marker and test counts, and `test_id` recomputed |
| 7 | **Every** fixed-model marker with a P is refitted exactly (statsmodels `Logit`, Newton, tol 1e-10) on the same samples and covariates. An `INFO` line counts the markers plink2 misses by more than 1e-3 in log10 P (METHODS §3, numerical note). PASS requires that the ORs agree within 1e-3 and that no call at any of the three thresholds changes when the exact P replaces plink2's |
| 8 | SAIGE rows: BETA = Tstat / var, `p.value.NA` = 2Φ(−\|Tstat\| / √var), and in the table `or` = exp(BETA) and `p` = `p.value.NA` |
| 9 | `null_check.json` has no failure, and the null model's samples equal the case set's random-model samples |
| 10 | Per case set: `_comparison/assoc_all.tsv` equals the 12 per-cohort tables, and `scan_qc_all.tsv` has 12 scans. The n/3 consistency of every marker at every threshold is recomputed from the random-model rows (same direction, majority), and `hits.tsv` holds exactly the markers significant somewhere. `_case_sets/scan_qc_case_sets.tsv` has 48 distinct scans |
| 13 | Conditional analyses. SENTINEL is the first record of every pooled set, and the samples equal the case set. The pooled HLA columns equal the matrices cell by cell, and the lead equals the `random_model` genotypes cell by cell. For every (case set, model, class) walk: step 0 equals the class scans (random model: every P identical; fixed model: an `INFO` line counts plink2's file-dependent deviations), and the joint selection rule is recomputed step by step from `stepwise.<class>.tsv` (candidates significant at GEC in some cohort, the lead's largest P over the cohorts the smallest, entry, at most 5 signals, members). Every signal lead a cohort tested is among that cohort's lead \| HLA conditioning sets. For case set `all`, SAIGE with `best_guess` on the pooled set reproduces the genome-wide scan's P for each lead in every cohort. Every lead | whole-position set holds the position's residue columns but its most frequent (recomputed from the pooled genotypes), for the residue signal leads' and the omnibus leads' positions |
| 14 | Benchmark (`MAXT_BOOT`) per case set × cohort × class: 10,000 replicates on the fixed-model samples, with the scan's cases, controls, tests and thresholds. The observed phenotype's Wald rerun equals `ASSOC_FIXED`'s P on every test, and its score P is recomputed for 10 tests by weighted least squares. Every exact threshold and family-wise error rate is recomputed from the replicates |
| 15 | Omnibus tests of residue positions per case set: no gene-level (allele) omnibus output exists. For 5 positions per cohort the categories, reference and df are rebuilt from the matrix by the documented rule, and the likelihood-ratio test is refitted with statsmodels. The Bonferroni threshold (0.05 / units tested) and every call are recomputed. The stepwise rule (candidates, lead, members) is recomputed step by step. `_omnibus/omnibus_all.tsv` equals the per-case-set tables |
| 16 | HLA \| signal and positions \| one marker (`HLA_GIVEN_SIGNAL`, `OMNIBUS_GIVEN`). Per (case set, cohort, model): every marker's marginal P equals `CONDITIONAL`'s; the GWAS lead given each conditioning set equals `lead_given_hla.tsv`; the residues given a residue signal lead equal `STEPWISE`'s step 1 for that lead. Per case set: the positions' marginal P in `omnibus_given.residue.tsv` equal `OMNIBUS`'s, with one row per position for every cohort and conditioning marker. The scripts stop on their own if the recomputed marginal scan differs |
| 11 | `trace.txt` has one completed task per expected process instance. Expected counts: GENE_COORDS 1, CASE_SETS 1, NULL_INFO_ALL 3, VALIDATE_INPUTS 3, PREP 12, FIT_NULL 9, NULL_INFO_FIT 9, CHECK_NULL 12, BUILD_MARKERS 24, EFFECTIVE_TESTS 24, ASSOC_FIXED 24, ASSOC_RANDOM 24, COLLECT 48, COMPARE 4, FIGURES 8, CASE_SET_COMPARE 1, GWAS_LEADS 1, EXTRACT_LEAD 3, BUILD_POOLED 12, STEPWISE 16, CONDITIONAL 24, CONDITIONAL_COMPARE 1, WRITE_RUN_MANIFEST 1 |

Last run, 2026-10-26 (run 17, with HLA_GIVEN_SIGNAL and OMNIBUS_GIVEN): **814 passed, 0 failed**.
