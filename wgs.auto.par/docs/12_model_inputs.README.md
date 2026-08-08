# `results/12_model_inputs/` — association model inputs

Per **cohort mainland sample set**, this stage prepares the genotypes, covariates and
phenotype needed by the downstream (external) association analyses, plus a consolidated
per-set manifest. Sets are nested: **`narrow_mainland` (2,195) ⊂
`intermediate_mainland` (2,508) ⊂ `full_mainland` (3,101)**. The BBJ reference
(`reference_full_mainland`) is the PCA base, not a cohort set, and never appears here.

Produced by `select.auto.par.v6.nf` (processes `POPGMM_SUBSET_AND_PLOT_BBJ_PROJECTION`,
`POPGMM_PIHAT_INTERSECTION_PROJECTION`, `PREPARE_FIXED_MODEL_GENOTYPE`,
`PREPARE_RANDOM_MODEL_GENOTYPE`, `PREPARE_POPGMM_COV_PHENO_FILES`, `SUMMARIZE_MODEL_INPUTS`).

## Layout

```
12_model_inputs/<set>/
  genotype/
    base/          <pre>.popgmm.{bed,bim,fam}            all set samples, monomorphic-removed
    fixed_model/   <pre>.maf_ge_threshold.{bed,bim,fam}  COMMON (MAF>=thr)  -> single-variant tests
                   <pre>.maf_lt_threshold.{bed,bim,fam}  RARE   (MAF<thr)   -> burden/aggregate tests
                   <pre>.maf_ge_threshold.variants.txt, <pre>.maf_ref.afreq, <pre>.fixed_ready.*
    random_model/  <pre>.random_model.{bed,bim,fam}      COMMON variants, ALL samples -> GRM/random effect
  covariates/      geno_pc.sex.tsv          geno_pc.sex_age.tsv
                   bbj_pc.sex.tsv           bbj_pc.sex_age.tsv
                   bbj_mainland_pc.sex.tsv  bbj_mainland_pc.sex_age.tsv
  phenotype/       pheno.tsv
  pca/
    bbj/           full-BBJ projection sscore + variance_summary.png + pc_pairs.pdf
    insample/      relatedness-aware in-sample PCA sscore + figures + eigen/base artifacts
    bbj_mainland/  <set>.mainland.sscore + <set>.bbj_mainland.variance_summary.png + .pc_pairs.pdf
  logs/            *_prep.log.txt, cov_pheno.log.txt, pihat_*exclude.fid_iid, age_na.fid_iid,
                   sex_na.fid_iid, *.done.txt, subset.log.txt
  MANIFEST.tsv     one row per prepared artifact (category, file, n_samples, n_variants_or_pcs, note)
  prep_summary.log human-readable review block with consistency checks
```

## Model genotypes

- **Fixed model** (fixed-effects test on unrelated samples): remove PI_HAT-related samples →
  remove `mac=0` monomorphic variants (`plink2 --mac 1`) → split by MAF threshold (group
  `ctrl`, thr `0.01`, estimated in that group):
  - `maf_ge_threshold` = **common** (MAF≥thr) → single-variant tests
  - `maf_lt_threshold` = **rare** (MAF<thr) → burden/aggregate tests
- **Random model** (GRM / random effect): **no sample removal**; extract exactly the
  fixed-model **common** variant list from all set samples → `random_model.{bed,bim,fam}`.

## Covariates & phenotype

Three PC sources, each in two variants (`--n-pcs 20`, PC columns are `PC1_AVG..PC20_AVG`):

| label             | PC source                                                        |
|-------------------|-----------------------------------------------------------------|
| `geno_pc`         | relatedness-aware **in-sample** PCA of the set (project all back)|
| `bbj_pc`          | projection onto the **full-BBJ** PC space                        |
| `bbj_mainland_pc` | projection onto the **BBJ-mainland** PC space                    |

- `<label>.sex.tsv`      : `#FID IID SEX PC1_AVG..PC20_AVG`
- `<label>.sex_age.tsv`  : `#FID IID SEX AGE AGE_Z PC1_AVG..PC20_AVG` (age-NA samples dropped)
- `pheno.tsv`            : `#FID IID PHENO1` (PLINK 1=control / 2=case)

Covariate/phenotype files are the **full set superset** (all set samples with valid sex).
Downstream tools intersect by ID: the fixed model uses `genotype/fixed_model/*` (unrelated
samples), the random model uses all samples. `logs/sex_na.fid_iid` and `logs/age_na.fid_iid`
record samples dropped for missing sex/age. Sex map: male→1, female→2 (`Sex` F/M); age from
`Age_at_DNA_Collection`.

## Reproducibility & logs

Every prep step writes a structured `*_prep.log.txt` (counts, deltas, thresholds, tool
version). `MANIFEST.tsv` + `prep_summary.log` consolidate one set's artifacts with sanity
checks (e.g. `base − related_removed == fixed`, `random_variants == fixed_common`). The
in-sample and BBJ-mainland PCAs use `plink2 --pca approx --seed 42` for reproducibility.
