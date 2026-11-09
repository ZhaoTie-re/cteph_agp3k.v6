# assoc_hla.v2 — HLA allele and amino-acid residue association, CTEPH

This component tests the HLA markers produced by `hla.typing.v2` for association with CTEPH:

- **case definitions:** four. Three are against all controls: all cases (`all`), cases with a history of DVT (`dvt_yes`) and cases without one (`dvt_no`). The fourth, `dvt_case_only`, is a case-only test of cases with DVT against cases without DVT;
- **marker classes:** 2-field alleles and amino-acid residues at IMGT positions, as 0/1/2 copies;
- **cohorts:** three nested cohorts (narrow, intermediate and full mainland);
- **models:** two, both adjusted for sex and 10 BBJ-mainland PCs;
- **levels of testing:** two main analyses. Single markers (each allele or residue) are tested in both models, with the random model primary. Whole residue positions are tested by omnibus tests in the fixed model; genes are not tested as a whole.

| model | engine | samples | missing calls | P |
|---|---|---|---|---|
| fixed | plink2 `--glm` logistic, `no-firth` | fixed_model set (relatives removed) within the case set | complete case | Wald |
| random | SAIGE step 2. `all` uses the genome-wide scan's **reused** null model; each DVT set uses its own null model, fitted here with the same settings and GRM markers (full GRM, LOCO) | all samples of the case set | mean-imputed | unadjusted score P (`p.value.NA`); the saddle-point P is kept |

**Significance.** Thresholds are fixed per (case set, cohort, class) from the genotypes, before any P, at a family-wise 0.05 over the class's tests (a two-residue mirror pair is one test):

| role | threshold |
|---|---|
| primary | **LD-adjusted (GEC)**: 0.05 / Me, Me from the original GEC v0.2 (archived in `tools/gec-0.2/`) |
| validation | **Exact (max-T)**: the family-wise threshold under the joint null of the covariate-adjusted genotype correlation |
| sensitivity | **Conventional Bonferroni**: 0.05 / tests |

Benjamini–Hochberg is not used, because LD makes its calls depend on how many correlated markers move together. **For single markers, the random model is primary.** A marker is reported as significant in n/3 cohorts, with the same direction of effect, per threshold, and alleles and residues are counted separately; 1/3 is "Single cohort". λGC is reported and never applied.

**Omnibus tests of residue positions (a main analysis).** Each residue position is also tested as a whole, by logistic regression on the fixed-model samples. The model is PHENO1 ~ 1 + SEX + PC1–10 + one copy column per residue (plus `unlisted` and a pooled `rare`), with the most frequent residue left out. It is compared with the model without the position by a likelihood-ratio test, df = residues − 1, at Bonferroni over the positions tested (0.05 / 290–296). Independent positions are selected stepwise, conditioning on all residues of each entered position. Gene-level omnibus tests over alleles are not run (METHODS §8). See [docs/METHODS.md](docs/METHODS.md) §8.

This component was designed and run independently of `analysis/assoc_hla`: none of that component's files, scripts, null models or results is read.

## What it establishes, and what it does not

- **What it gives.** For every tested marker in every case set, cohort and model: OR, 95 % CI and P, with counts in cases and controls. It also gives the fixed-against-random agreement, each marker's cross-cohort consistency at the three thresholds, the independent signals within each class, and the DVT-yes and DVT-no estimates side by side.
- **The DVT difference is tested in `dvt_case_only`, not by setting the two DVT sets side by side.** Those two sets share their controls. The case-only test avoids the case–control platform split, but the platform mix differs between cases with and without DVT (METHODS §2.1), and platform is not a covariate.
- **No proof of an HLA risk factor.** Cases and controls were sequenced on different platforms, and HLA typing can depend on read quality.
- **No replication.** The cohorts are nested and share almost all cases. Cross-cohort consistency is a sample-selection sensitivity analysis: it shows that a signal survives adding controls.
- **Conditional analyses separate signals statistically, not causally.** Independent signals are selected within each class, jointly over the three cohorts. The genome-wide scan's MHC lead is conditioned on HLA markers and the reverse (`results/_conditional/`).

See [docs/METHODS.md](docs/METHODS.md) §6–8.

## Run

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_hla.v2
source activate dsl2
nextflow lint assoc_hla.v2.nf                  # syntax check; never `run --help` / -preview
nextflow run assoc_hla.v2.nf -resume -ansi-log false > logs/run.log 2>&1
bash verify.sh                                 # read-only; exit 0 = all checks passed
```

The tests are small: SAIGE step 2 needs about 5 s per set. The one heavy step is `FIT_NULL`, nine full-GRM null fits (three DVT sets × three cohorts) of about 1 h each, run in parallel.

## Where things are

| | |
|---|---|
| results, one schema | `results/<case_set>/_comparison/assoc_all.tsv` |
| significant markers, n/3 per threshold | `results/<case_set>/_comparison/hits.tsv` |
| tests and thresholds | `results/<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv` |
| scan summaries | `results/<case_set>/_comparison/scan_qc_all.tsv`, `model_concordance.tsv` |
| figures | `results/figures/`: main `01.landscape` … `05.dvt`, supplementary `06.thresholds` … `08.qq`, each with a `README.md` ([docs/FIGURES.md](docs/FIGURES.md)) |
| omnibus tests of residue positions | `results/_omnibus/`: start with `omnibus_units.tsv` (each position beside its best single residues) and `omnibus_signals_all.tsv`; per case set `results/<case_set>/_omnibus/` |
| independent signals; GWAS MHC lead vs HLA | `results/_conditional/`: start with `signals_all.tsv` and `consistency.tsv`; figures `results/figures/03.conditioning`, `04.gwas_link` |
| every HLA marker given a signal lead or whole position; positions given the GWAS lead or a residue lead | `results/<case_set>/<cohort>/04.conditional/<model>/hla_given_signal.tsv`; `results/<case_set>/_omnibus/omnibus_given.residue.tsv` |
| DVT yes against DVT no | the test: `results/dvt_case_only/`; side by side: `results/_case_sets/`: `dvt_contrast.tsv`, `dvt_hits.tsv`, `scan_qc_case_sets.tsv`; figures `results/figures/02.signals`, `05.dvt` |
| measured findings | [docs/STUDY_NOTES.md](docs/STUDY_NOTES.md) |
| methods | [docs/METHODS.md](docs/METHODS.md) |
| column dictionary | [docs/OUTPUTS.md](docs/OUTPUTS.md) |
| figure system | [docs/FIGURES.md](docs/FIGURES.md) |
| checks | [docs/VERIFICATION.md](docs/VERIFICATION.md) |

## Files

| file | role |
|---|---|
| `assoc_hla.v2.nf`, `nextflow.config` | pipeline; SLURM resources via `--rsc` only |
| `tools/gec-0.2/` | GEC v0.2 (`gec.jar`, the original zip and manual), with `SOURCE.md`: where it came from, sha256, how it is run |
| `scripts/gene_coords.py` | gene spans from RefSeq hg38, used for order and placement only |
| `scripts/case_sets.py` | DVT status (yes / no / unknown) of every sample, from the v6 sample sheet |
| `scripts/null_model_info.R` | reads back a SAIGE null model, reused or fitted |
| `scripts/check_null.py` | proves each case set's null model matches its samples, labels and covariates |
| `scripts/validate_inputs.py` | 30 input checks per cohort, run before any test |
| `scripts/prep_pheno_cov.py` | phenotype and covariates for both models, from one join |
| `scripts/build_markers.py` | matrix → VCF hard calls → PLINK, plus per-marker counts |
| `scripts/effective_tests.py` | tests (mirror pairs once) and the three thresholds: GEC v0.2, max-T, Bonferroni |
| `scripts/collect.py` | engine output → one schema, significance at the three thresholds, scan summary, gates |
| `scripts/compare.py` | per case set: stacks the 12 scans, compares models, counts each marker's cohorts (n/3) per threshold |
| `scripts/case_set_compare.py` | the four case sets side by side (tables) |
| `scripts/gwas_leads.py` | the genome-wide random-model scan's MHC leads that pass the three-cohort rule |
| `scripts/build_pooled.py` | alleles + residues + GWAS lead(s) in one PLINK set (sentinel first), with LD to the lead |
| `scripts/engine.py` | one plink2 / SAIGE run on the pooled set, optionally conditioned, with its gates |
| `scripts/stepwise.py` | independent signals of one class, the three cohorts jointly |
| `scripts/conditional.py` | HLA given the GWAS lead and the lead given HLA markers and the signal leads |
| `scripts/conditional_compare.py` | conditional results across case sets, cohorts and models (tables) |
| `scripts/omnibus.py`, `scripts/omnibus_summary.py` | omnibus tests of residue positions, the three cohorts jointly, with stepwise selection; their summary beside the single-marker scans |
| `scripts/figures.py`, `scripts/catalogue.py`, `scripts/vocab.py` | the eight figure families, one README per family, shared vocabulary and encodings |
| `scripts/verify_checks.py`, `verify.sh` | post-run checks |
