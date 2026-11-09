# assoc_hla.v2 — methods

HLA allele and amino-acid residue association with CTEPH. The markers come from `hla.typing.v2`.

The design has four axes:

- **Four case definitions.** Three compare cases against all of the cohort's controls: all cases, cases with a history of DVT, and cases without one. The fourth, a case-only test, compares cases with DVT against cases without DVT.
- **Three nested cohorts.**
- **Two marker classes:** alleles and amino-acid residues.
- **Two models:** a fixed-effect logistic regression and a random-effect logistic mixed model. Both adjust for sex and 10 BBJ-mainland principal components.

This component was designed and run without reading any file, script, null model or result of `analysis/assoc_hla`.

## 1. Markers

| | |
|---|---|
| source | `hla.typing.v2/results/06.assoc_inputs/<cohort>/{allele,residue}_dosage.tsv`, read unchanged |
| allele column | copies (0/1/2) of one 2-field allele in a sample's two copies of the gene |
| residue column | copies (0/1/2) carrying one amino acid at one IMGT position |
| blank cell | the copy was undetermined at typing; analysed as missing |
| tested per cohort | alleles 143 / 144 / 142, residues 672 / 681 / 679 (narrow / intermediate / full mainland) |

The marker set is the typing component's, after its QC:

- the marker varies in the cohort;
- DRB3/4/5 are left out;
- call rate is at least 95 % in cases and at least 95 % in controls;
- MAF is at least 0.5 %.

Nothing is added, dropped or re-filtered here. `VALIDATE_INPUTS` proves the copy read is the verified export. It checks the md5 against `summary.json`, the columns against the `tested` rows of `markers.tsv`, and the rows against the phenotype file in the same order.

**Coding.** `BUILD_MARKERS` writes each matrix as VCF hard calls:

| matrix cell | genotype |
|---|---|
| 0 | 0/0 |
| 1 | 0/1 |
| 2 | 1/1 |
| blank | ./. |

REF is `A`, meaning the copy does not carry the marker; ALT is `P`, meaning it carries it. `plink2 --vcf --double-id --make-bed` turns the VCF into one PLINK set that both engines read. Every effect is therefore per copy carrying the marker. `COLLECT` asserts this from each engine's output: the plink2 `A1` and the SAIGE `Allele2` are `P` on every row.

**Positions.** Markers sit on chromosome 6 at the gene start plus the marker's rank within its gene. Gene starts are GRCh38, taken from the union of the gene's chr6 transcripts in RefSeq `hg38_refGene.txt`. Within a gene, residues are ranked by IMGT position and then residue letter; alleles by their numeric fields. These positions are unique labels in chromosomal gene order, not the physical position of a residue. The only thing the random model's LOCO step needs from them is the chromosome.

**Mirror columns.** At a position with exactly two residues, the two residue columns are complements in every sample (blanks in the same samples). They return the same P and reciprocal ORs, and are counted as one test (§5). Residues on one haplotype are also strongly correlated with each other and with the alleles that carry them. The number of significant markers is therefore not the number of independent signals (§7).

## 2. Samples, phenotype, covariates

| case set `all` | narrow | intermediate | full mainland |
|---|---|---|---|
| random model: samples (cases / controls) | 2,195 (419 / 1,776) | 2,508 (429 / 2,079) | 3,101 (439 / 2,662) |
| fixed model: samples (cases / controls) | 2,168 (419 / 1,749) | 2,480 (429 / 2,051) | 3,071 (439 / 2,632) |

**Sources.**

- Phenotype: `wgs.auto.par/results/12_model_inputs/<cohort>/phenotype/pheno.tsv`. `PHENO1` is coded 0 = control, 1 = case.
- Covariates: `…/covariates/bbj_mainland_pc.sex.tsv`. The model uses `SEX` and `PC1_AVG`–`PC10_AVG`, the 10 BBJ-mainland projected PCs; this choice is documented with the genome-wide scans.

### 2.1 Case sets

| case set | cases | controls |
|---|---|---|
| `all` | every case | all |
| `dvt_yes` | cases whose `DVT_Presence` is `Yes` | all |
| `dvt_no` | cases whose `DVT_Presence` is `No` | all |
| `dvt_case_only` | cases whose `DVT_Presence` is `Yes`, coded 1 | none; the comparison group is cases whose `DVT_Presence` is `No`, coded 0 |

DVT status comes from `info/cteph_agp3k.v6.20260507.xlsx`. The `DVT_Presence` column is joined to the samples on `ID_JHRPv6` (`CASE_SETS`, written to `results/_reference/case_sets.tsv`). The column is recorded for cases only. `Unknown` and blank are both "unknown", and those cases enter `all` only; any other value stops the run. `VALIDATE_INPUTS` checks that every case has a row in the sheet and that no control carries a status.

| cases (random model) | narrow | intermediate | full mainland |
|---|---|---|---|
| all | 419 | 429 | 439 |
| DVT yes | 150 | 155 | 159 |
| DVT no | 140 | 141 | 145 |
| unknown (in `all` only) | 129 | 133 | 135 |

`dvt_case_only` holds 290 / 296 / 304 cases. It recodes `PHENO1` within the cases (1 = DVT yes, 0 = DVT no), and its tables keep the generic column names: `*_case` means DVT yes and `*_control` means DVT no. It is the direct test of whether a marker's frequency differs between the two kinds of case. The three cohorts share almost all their cases and this set has no controls, so the case-only test is nearly the same analysis three times.

**Platform differs with DVT status.** Among cases (full mainland), DNBSeq-T7 30× makes up 49 % of DVT-no cases (71 / 145) but 77 % of DVT-yes cases (123 / 159). DNBSeq-G400RS 30× is 31 / 145 against 8 / 159, and NovaSeq 30× is 29 / 145 against 16 / 159 (χ² = 30.8, 3 df, P = 9e-7; narrow mainland is the same). The case-only test therefore avoids the case–control platform split, but it is not free of platform effects. Platform is not a covariate in any model here.

The marker set is the cohort's, whatever the case set: the typing component's QC ran on the whole cohort. A marker can therefore be rarer in a DVT set's cases than its cohort MAF suggests.

### 2.2 Sample sets of the two models

Each model uses the same sample set as its genome-wide scan, restricted to the case set:

- **Fixed model:** the cohort's `genotype/fixed_model/*.maf_ge_threshold.fam`, which has relatives removed upstream.
- **Random model:** every sample in the phenotype file. The GRM in the null model accounts for relatedness.

A case–control set keeps all controls and only its own cases in both sample sets; `dvt_case_only` keeps the cases of known DVT status only. `PREP` writes both tables from one join, so the two models cannot disagree about labels or covariates. `BUILD_MARKERS` writes each PLINK set with exactly the case set's random-model samples.

## 3. Fixed model

```
plink2 --bfile <class> --keep fixed.keep
       --pheno pheno_covar.fixed.tsv --pheno-name PHENO1 --1
       --covar pheno_covar.fixed.tsv --covar-name SEX,PC1_AVG-PC10_AVG
       --glm omit-ref no-firth hide-covar cols=+a1countcc,+totallelecc --ci 0.95
```

- **Engine and settings:** logistic regression, with the same binary (`/home/b/b37974/plink2`, v2.0.0-a.5.17LM) and the same `no-firth` setting as the genome-wide fixed scan (`analysis/assoc_plink2`).
- **Missing calls:** a sample with a blank call is left out of that marker only (complete case). The number analysed is reported as `obs`.
- **Phenotype coding:** `--1` makes plink2 read 0/1 as control/case. Without it, 0 would be read as missing.
- **Gate:** `COLLECT` stops the run unless plink2's log reports the expected cases, controls and 11 covariates.
- **Outputs:** OR, Wald 95 % CI and Wald P from plink2. `ERRCODE` is kept as `errcode`, where `.` means no error.

**Numerical note: plink2 can stop its logistic fit short.** plink2's logistic fit stops short of the maximum-likelihood estimate for a few markers, and which markers depends on the other variants in the file. The same happens in v2.0.0-a5.17 and a6.20. The OR is barely affected, but the Wald SE, and so P, can be.

- **Example (case-only set):** DQA1:-7:M and its mirror DQA1:-7:V must have identical P. From the conditional analysis's pooled file plink2 gave P 0.0068 and 0.0036, with an SE of 0.433 instead of 0.466. From the residue file it gave 0.0068 for both, the exact value.
- **Extent:** every fixed-model test was refitted exactly (statsmodels, Newton, tol 1e-10). About 30 of 9,736 deviated by more than 1e-3 in log10 P, the largest 0.22 (DOA\*01:01, case-only set, P 0.19 against an exact 0.31).
- **Impact:** no call at any of the three thresholds (§5) changes when the exact P replaces plink2's.

plink2 is kept as the engine, the same one as the genome-wide fixed scan. `verify.sh` §7 refits every marker, reports the deviations, and fails if any call would change.

## 4. Random model

```
step2_SPAtests.R --bedFile/--bimFile/--famFile=<class> --AlleleOrder=alt-first --chrom=6
       --GMMATmodelFile=null.rda --varianceRatioFile=null.varianceRatio.txt
       --is_Firth_beta=FALSE --LOCO=TRUE --is_output_moreDetails=TRUE
       --minMAC=0.5 --maxMissing=0.15 --impute_method=mean
```

**Null model for `all`, reused and not refitted.** The model is `analysis/assoc_saige/results/<cohort>/01.assoc/main/null.rda`, the genome-wide SAIGE scan's own null model. It is a logistic mixed model with SEX + PC1–PC10, a full GRM from LD-pruned genome-wide markers, and LOCO. Variance ratios are 0.928 / 0.931 / 0.918.

**Null models for the DVT sets, fitted here.** A null model is specific to its phenotype and samples, so the DVT sets, including `dvt_case_only`, cannot reuse the genome-wide one. `FIT_NULL` fits one per (DVT set, cohort) with the genome-wide scan's settings and GRM markers:

- **Genotypes:** `12_model_inputs/<cohort>/genotype/random_model`, restricted to the scan's published `00.prep/grm.prune.in` (LD-pruned, high-LD regions excluded; about 378 k markers) and to the case set's samples.
- **SAIGE step 1 settings:** `--traitType=binary --isDiagofKinSetAsOne=True --numRandomMarkerforVarianceRatio=200 --useSparseGRMtoFitNULL=FALSE --LOCO=TRUE --isCovariateOffset=FALSE`, with covariates SEX, PC1_AVG … PC10_AVG.
- **Outputs:** `results/<dvt set>/<cohort>/00.prep/null_model/`.

For every null model, reused or fitted, `NULL_MODEL_INFO` reads back the model's sample IDs, case labels, covariates and LOCO flag. `CHECK_NULL` then requires the same things of each case set's null model before any random-model test runs. For the reused models, `VALIDATE_INPUTS` additionally requires, before anything else, that:

- the sample set equals the phenotype file;
- every case/control label agrees;
- the covariates are exactly SEX, PC1_AVG … PC10_AVG.

The variance ratio was estimated from common genome-wide markers. It is applied unchanged to the HLA markers, all of which have MAF ≥ 0.5 %.

**LOCO.** The chromosome-6-left-out predictor is used, so the HLA region's own effect is not absorbed into the random effect. The SAIGE log line `Leave chromosome 6 out will be applied` is asserted in `COLLECT`.

**Missing calls.** Missing calls are mean-imputed: `--impute_method=mean`, asserted from the log. SAIGE's default, best_guess, would set each undetermined copy to the most common genotype. For a rare marker that means "absent", which would turn every undetermined copy into a confident non-carrier.

**Primary P.** The primary P is `p.value.NA`, the unadjusted score test. This matches the genome-wide SAIGE scan. The saddle-point P (`p.value`) is kept as `p_spa`. All markers have MAF ≥ 0.5 % and MAC ≥ 20 in the cohort, where the normal approximation of the score statistic holds for `all`. The DVT-yes and DVT-no sets are more unbalanced, with about 140–160 cases against 1,776–2,662 controls. There, a marker that is rare among the set's cases is where the unadjusted P can be anti-conservative, which is exactly what SPA corrects. For any DVT-set hit, compare `p` with `p_spa`.

**Effect size.** `beta` is SAIGE's one-step estimate, Tstat / var. `se` is 1 / √var, the standard error that matches `p.value.NA`; `verify.sh` §8 checks this identity on every row. OR is exp(beta), and the 95 % CI is exp(beta ± 1.96·se). These are approximations to the mixed-model MLE; SAIGE's Firth-corrected estimate was not requested, as in the genome-wide scan. Read them for direction and rough size.

## 5. Significance

### 5.1 Tests

Thresholds are set per (case set, cohort, class) by `EFFECTIVE_TESTS`, from the genotypes alone, before any association P exists. Both models use the same thresholds.

- **Tests** are the class's markers that vary among the case set's random-model samples.
- **A mirror pair counts once.** Two residues at one IMGT position whose dosages are complements in every sample are one test reported twice (`test_id` names the pair's first residue in file order, `mirror_of` the other).
- **Alleles are counted as they are.** No two allele columns are complements.

| case set `all` | narrow / intermediate / full mainland |
|---|---|
| alleles: markers = tests | 143 / 144 / 142 |
| residues: markers | 672 / 681 / 679 |
| residues: mirror pairs | 178 / 178 / 167 |
| residues: tests | 494 / 503 / 512 |

In the DVT-yes and DVT-no sets the counts are within a few of these. `dvt_case_only`, with about 300 cases, has more complement pairs (215–218) and 457 / 463 / 462 residue tests. Every case set's counts are in `results/<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv`.

### 5.2 Three thresholds

All three control the family-wise error rate at 0.05 over the class's tests. They differ in how they treat the strong LD among HLA markers.

| role | name in the tables | rule |
|---|---|---|
| **primary** | LD-adjusted (GEC), `gec` | P < 0.05 / Me |
| validation | Exact (max-T), `maxt` | P < the exact family-wise threshold under the joint null |
| sensitivity | Conventional Bonferroni, `bonf` | P < 0.05 / tests |

**LD-adjusted threshold (GEC).** Me is the effective number of independent tests of Li et al. (Hum Genet 2012;131:747–756), computed with the original standalone GEC v0.2. This is the release used for the same purpose in HLA fine-mapping of leprosy (Dallmann-Sauer et al., PLoS Pathog 2020). The jar is archived in `tools/gec-0.2/` with its source and sha256, which `EFFECTIVE_TESTS` and `verify.sh` check. The GEC function of KGGSEE, the site's later release, is not used. The run:

```
java -Xmx4g -jar gec.jar --effect-number --plink-binary <tests> --genome --maf 0 --no-web --out gec
```

- **Input:** one marker per test, from the case set's fixed-model samples. These are unrelated, as GEC requires.
- **Alleles:** `P` is written as `T`, because GEC reads nucleotide alleles only.
- **`--maf 0`:** v0.2 drops markers with MAF < 0.01 by default, and every tested marker (MAF ≥ 0.5 %) must count.
- **Blocks:** v0.2 estimates Me within blocks of at most 100 consecutive markers and sums over the blocks. This is not settable, and LD between blocks is ignored, which makes Me conservative.
- **Monomorphic tests:** a test that is monomorphic among the unrelated samples cannot enter GEC. It is added to Me as one independent test (`gec_added`).

**Exact threshold (max-T).** This is the 0.05 quantile of the smallest P over the tests under the joint null (Westfall–Young max-T). It is computed as max |Z| from a multivariate normal (Conneely & Boehnke, AJHG 2007;81:1158–1168) whose correlation is that of the class's dosages:

- after regressing out sex and the 10 PCs;
- with blanks mean-imputed, as in the random model;
- over the case set's random-model samples.

The simulation uses 400,000 draws (seed 20261023) and gives a 95 % Monte Carlo interval (`p_maxt_lo`, `p_maxt_hi`). It also gives the family-wise error rate that the GEC and conventional-Bonferroni thresholds actually have (`fwer_gec`, `fwer_bonf`). The max-T threshold is reported to show how far the primary threshold is from exact. It is not used to make calls on its own.

**Benchmark by parametric bootstrap (`MAXT_BOOT`).** No call depends on it. The multivariate-normal max-T assumes normal test statistics, so it cannot see how the P values of rare markers behave when cases are few.

- **Null phenotypes:** the null logistic model (PHENO1 ~ 1 + SEX + PC1–10) is fitted on the case set's unrelated samples, and 10,000 null phenotypes are drawn from it (y* ∼ Bernoulli(μ̂), seed 20261023). The covariates keep their effect; the genotypes have none. Case labels are not permuted, because the PCs predict case status.
- **Two statistics on the same replicates:**
  - plink2's Wald test with the fixed model's settings (the observed phenotype is rerun too, and equals `ASSOC_FIXED`);
  - the score test without SPA, blanks mean-imputed, i.e. the random model's statistic without the GRM.
- **Outputs:** the exact threshold is the 0.05 quantile of the smallest P, and each threshold's family-wise error rate is the share of replicates whose smallest P is below it.

Family-wise error rate of the GEC threshold by case set (range over cohorts):

| case set (case fraction) | alleles: Wald / score / normal | residues: Wald / score / normal |
|---|---|---|
| all (17 %) | 0.045–0.058 / 0.065–0.075 / 0.053 | 0.041–0.047 / 0.049–0.055 / 0.046 |
| DVT yes, DVT no (6–7 %) | 0.070–0.080 / 0.115–0.126 / 0.053 | 0.049–0.052 / 0.073–0.077 / 0.046 |
| case-only (52 %) | 0.014–0.016 / 0.031–0.036 / 0.052 | 0.027–0.033 / 0.039–0.049 / 0.047 |

With few cases, the score statistic without SPA makes rare markers' small P values too frequent: the GEC threshold's real error rate then exceeds 0.05. This is the per-test calibration that SPA corrects, not a failure of the LD adjustment. The calls hold with the SPA P (`p_spa`):

- all cases: the three DQB1 57 V calls are unchanged;
- DVT no: 31 of the 34 test × cohort calls remain, DQB1 87 Y and DPB1\*05:01 among them.

Their leads are common markers, which SPA barely changes.

**Conventional Bonferroni** ignores LD, so it is conservative. It is reported as a sensitivity analysis.

| case set `all` | narrow | intermediate | full mainland |
|---|---|---|---|
| alleles: Me (GEC v0.2) | 114.0 | 115.0 | 113.6 |
| alleles: **LD-adjusted (GEC)** | **4.39e-04** | **4.35e-04** | **4.40e-04** |
| alleles: exact (max-T) [95 % MC interval] | 4.14e-04 [4.09e-04–4.19e-04] | 4.15e-04 [4.08e-04–4.18e-04] | 4.13e-04 [4.07e-04–4.19e-04] |
| alleles: conventional Bonferroni | 3.50e-04 | 3.47e-04 | 3.52e-04 |
| alleles: family-wise error of GEC / Bonferroni under max-T | 0.053 / 0.042 | 0.052 / 0.042 | 0.053 / 0.043 |
| residues: Me (GEC v0.2) | 201.8 | 200.8 | 200.9 |
| residues: **LD-adjusted (GEC)** | **2.48e-04** | **2.49e-04** | **2.49e-04** |
| residues: exact (max-T) [95 % MC interval] | 2.77e-04 [2.74e-04–2.81e-04] | 2.73e-04 [2.70e-04–2.77e-04] | 2.74e-04 [2.70e-04–2.77e-04] |
| residues: conventional Bonferroni | 1.01e-04 | 9.94e-05 | 9.77e-05 |
| residues: family-wise error of GEC / Bonferroni under max-T | 0.045 / 0.020 | 0.046 / 0.020 | 0.046 / 0.020 |

Across all 24 (case set, cohort, class) families, the GEC threshold's family-wise error rate under max-T is 0.045–0.053: slightly above 0.05 for alleles (0.052–0.053) and slightly below for residues (0.045–0.047), where v0.2's blocks of at most 100 markers ignore LD between blocks. Conventional Bonferroni gives 0.019–0.021 for residues and 0.041–0.043 for alleles.

**Why not Benjamini–Hochberg.** In this region the tests are strongly correlated: mirror pairs, residues on one haplotype, and alleles and the residues they carry. A BH call then depends on how many correlated markers move together, not on the strength of any one signal. In this data (random model, intermediate mainland, residue tests), a DQ haplotype in DVT-no cases reached q = 0.007 through 59 correlated tests, while DQB1 57 V in all cases, with a smaller P, had q = 0.09. BH is therefore not used.

**Not corrected across case sets, cohorts or models.** These test the same markers on nested or overlapping samples, not independent hypotheses. A marker with no P is significant under none of the thresholds.

**λGC** = median χ²(1) / 0.455 is reported per scan and never applied. All markers lie in one region and are strongly correlated, so a few true or artefactual case–control differences lift the whole distribution. λGC here describes the region, not the calibration of the test. The calibration reference for these models is the genome-wide λ of the same scans (`analysis/assoc_plink2`, `analysis/assoc_saige`).

### 5.3 Cross-cohort consistency

The random model is the primary model, as in the genome-wide scan; the fixed and random models barely differ here (§6). For each case set, class and threshold, a marker is reported by the number of the three nested cohorts in which its random-model P is below that cohort's threshold, with the same direction of effect:

| label | meaning |
|---|---|
| 3/3 | significant in all three cohorts |
| 2/3 | significant in two |
| 1/3 | significant in one: "Single cohort", as in the gene-based analysis |

**Direction.** If significant cohorts disagree in direction, the count is that of the majority direction, and `<threshold>_direction_conflict` is set. This did not occur in any case set.

**Reading.**

- Alleles and residues are counted separately.
- The primary statement uses the LD-adjusted threshold. The exact and conventional thresholds are reported beside it.
- The cohorts are nested and share almost all cases, so cross-cohort consistency is a sample-selection sensitivity analysis, not replication.
- The fixed model is reported for its effect sizes and for model concordance (`gec_fixed` in the marker tables), not for consistency.

## 6. Comparison across cohorts, models and case sets

`COMPARE` stacks the 12 scans of a case set.

**Model concordance** (`model_concordance.tsv`, per cohort × class) gives:

- the Pearson correlation of −log10 P and of β between the fixed and random models;
- Spearman's ρ of P;
- the overlap of the markers each model finds significant at the LD-adjusted threshold.

**Marker summary** (`marker_summary.tsv`) gives, for each marker:

- its P and OR in all 6 (cohort, model) cells;
- its cross-cohort consistency under each threshold: `gec`, `maxt` and `bonf` as n/3, `*_n` as a count, and `*_cohorts` naming the cohorts;
- the direction of effect, and whether every random-model OR is on the same side of 1 (`same_direction`);
- `mirror_of`, set when the marker and its partner are a mirror pair in all three cohorts.

`hits.tsv` holds the markers significant under any threshold in at least one cohort (random model), or at the LD-adjusted threshold in the fixed model.

**Reading agreement across models.** The two models differ in three ways: relatives (dropped by the fixed model), missing calls (complete case vs mean imputation) and the estimator (Wald vs score). Agreement between them tests the model.

**Case sets.** `CASE_SET_COMPARE` (`results/_case_sets/`) sets the four case sets side by side:

- `scan_qc_case_sets.tsv` holds the 48 scan summaries.
- `dvt_contrast.tsv` has one row per (class, cohort, model, marker), giving OR, P and significance at each threshold for all cases, DVT yes, DVT no and the case-only test.
- `dvt_hits.tsv` lists every marker significant at some threshold in a DVT set or in the case-only test.
- The figures are `results/figures/02.signals/` (each signal in every comparison) and `05.dvt/` (log OR with DVT against without DVT, every test).

The two DVT sets share all their controls, so their estimates are positively correlated even under the null, and the gap between them is not a test. The test of a difference is the case-only set `dvt_case_only` (§2.1), whose calls `dvt_contrast.tsv` carries beside the two sets and the figure rings.

## 7. Conditional analyses

### 7.1 The GWAS MHC signal

The genome-wide random-model scan is `analysis/assoc_saige`, model `main`: SAIGE with SEX + 10 BBJ-mainland PCs, P = `p.value.NA`.

**How the lead is chosen (`GWAS_LEADS`).** The rule is the project's three-cohort rule at the scan's suggestive line:

- every peak overlapping the extended MHC (chr6:28,510,120–33,480,577, GRCh38) is a candidate;
- its lead is kept if its P < 1e-5 in all three cohorts, read from each cohort's own `sumstats.tsv`.

`results/_conditional/gwas_leads.tsv` lists every candidate with its per-cohort P and the reason it was kept or excluded. In this run one lead is kept:

| | narrow | intermediate | full mainland |
|---|---|---|---|
| chr6:32632861:A:G, effect allele G | P 3.2e-6 | 8.6e-7 | 7.3e-7 |

The lead lies between HLA-DRB1 and HLA-DQA1.

One peak is excluded: chr6:29770977:A:G is suggestive in narrow only (P 5.2e-6), with 2.7e-5 and 2.1e-5 in the other two cohorts.

**Where the lead's genotypes come from (`EXTRACT_LEAD`).** They are taken from each cohort's `random_model` genotypes, the set the genome-wide scan tested, keeping its REF/ALT.

**Imputation, and why the lead's P here differs from the scan's.** The scan imputed missing calls with SAIGE's default (`best_guess`). Here every SAIGE run uses `mean`, as for the HLA markers (§4). About 4 % of the lead's calls are missing, so its marginal random-model P here is somewhat weaker than the scan's: 7.0e-6 against 3.2e-6 in narrow. `verify.sh` §13 shows the difference is the imputation and nothing else. With `best_guess`, SAIGE on this component's pooled set reproduces the scan's P for the lead exactly.

### 7.2 The pooled set

`BUILD_POOLED` writes one PLINK set per (case set, cohort) with the cohort's alleles and residues (coded as in §1), the lead(s), and a first record `SENTINEL`.

**Why the sentinel exists.** SAIGE 1.1.9 step 2 has a defect when run with `--condition`: the output row of the first marker in the file carries another marker's statistics. In testing, the first residue took the allele count and P of a DQB1 residue, while every other row was correct. The sentinel, a copy of the first HLA marker, absorbs this, and every reader discards its row.

**How the defect is guarded against.** `engine.py` also requires that every conditional SAIGE run reproduces every marker's marginal statistics (allele count, `p.value.NA`) from the unconditional run.

**Conditioning order.** SAIGE expects the conditioning markers in file order, and they are passed that way.

**LD.** The set also records the Pearson r between each HLA marker's dosage and the lead's ALT dosage. It is computed over all samples of the set and over the PHENO1 = 0 group.

### 7.3 Engines

Every conditional run uses the scan's own engine (`engine.py`):

- **Fixed model:** plink2 `--condition-list`, with the same samples and covariates as §3. The numerical note of §3 applies: in the pooled file, a few step-0 P values differ slightly from the class scans'.
- **Random model:** SAIGE step 2 `--condition` on the case set's checked null model. P = `p.value.NA_c`, SE = 1/√`var_c`.

A marker that the conditioning set (nearly) determines gets no conditional P. In plink2 this shows as an error code. In SAIGE it is marked `COLLINEAR`, applied when `var_c` < 1e-6 × the marginal `var`.

### 7.4 Independent signals (`STEPWISE`)

Signals are selected per (case set, model, class), jointly over the three cohorts. Alleles and residues are walked separately, like the consistency counts; the GWAS lead is never a candidate.

1. **Step 0** is the marginal scan in each cohort. Its P values equal the class scans' (`verify.sh` §13).
2. **Step k** tests every class marker in each cohort, conditional on the leads entered so far. A cohort is conditioned only on the entered leads it tested.
3. **Candidates:** the markers significant at the cohort's LD-adjusted (GEC) threshold in at least one cohort.
4. **Lead:** the candidate whose largest P over the three cohorts is smallest. A marker a cohort did not test counts as P = 1 there. This picks the marker most consistently associated, not the one with the single smallest P.
5. **Entry:** the lead enters as an independent signal. The walk stops when no marker is a candidate, after at most 5 signals.

For each signal, `signals_all.tsv` reports:

- **Consistency:** at entry (conditional on the earlier signals) and marginally, under each threshold.
- **Members:** the markers significant in some cohort before its lead entered and in none after, the lead included. These are the markers the signal explains.
- **Carrier alleles:** for a residue signal, the tested 2-field alleles of the gene that carry the residue. They are read from `hla.typing.v2/results/05.exports/allele_match_map.tsv` (the reference allele each 2-field allele was matched to) and `residue_reference/<gene>.tsv` (that allele's residue at the position).

Every step is in `steps_all.tsv`, and every marker at every step in `<case_set>/_conditional/<model>/stepwise.<class>.tsv`. `verify.sh` §13 recomputes the selection rule step by step.

### 7.5 The GWAS lead and the HLA markers (`CONDITIONAL`)

These runs are per (case set, cohort, model):

1. **HLA | lead.** Every HLA marker is retested with the GWAS lead in the model.
2. **Lead | HLA.** The lead is retested with one HLA marker at a time in the model. The markers used are:
   - the 10 in highest |r| with it;
   - the 10 with the smallest marginal P;
   - the signal leads of both classes (§7.4).

   A separate run puts all signal leads in the model together.
3. **Lead | whole position.** The lead is retested with every residue column of a position in the model, except the position's most frequent residue (a position's copies sum to 2, so the last column would be redundant). The positions are those of the residue signal leads (§7.4) and the omnibus position leads of the same case set (§8). The reference residue is chosen from the pooled set's samples, so both models condition on the same columns. This asks whether the position as a whole explains the GWAS lead, not one residue only. `OMNIBUS` therefore runs before `CONDITIONAL`.

4. **HLA | signal (`HLA_GIVEN_SIGNAL`).** Every HLA marker, and the GWAS lead, is retested with each signal lead (§7.4) in the model, and with each whole position of item 3. The stepwise selection (§7.4) steps only within a class that has a signal, so it never retests, for example, the alleles given the residue lead; this run does. Same engine, samples and pooled set as items 1–3.
5. **Positions | one marker (`OMNIBUS_GIVEN`).** Every residue position is retested with its omnibus test (§8) with one marker added as a covariate: each GWAS lead (its dosage from the pooled set) and each residue signal lead of either model (from the class's set). Missing calls of that marker are mean-imputed, so every position keeps its samples. The scan without the marker is recomputed with the same code and must equal `OMNIBUS`'s, so the marker is the only difference. The omnibus stepwise walk (§8) conditions positions on whole positions only; this run gives positions given the GWAS lead and given a single residue.

**Reading the results.** `results/_conditional/signals_all.tsv` lists the independent signals. `consistency.tsv` reads the signals and the GWAS-lead results across the three cohorts (narrow / intermediate / full), with the random model first as the primary model.

## 8. Omnibus tests of residue positions (`OMNIBUS`; a main analysis)

The single-marker scans test one residue or allele at a time. The omnibus tests ask the second main question: is an amino-acid position associated as a whole? Both are main analyses, each with its own model and threshold: single markers in the random model at the GEC threshold, positions in the fixed model at Bonferroni. Each residue signal is read at both levels.

**No gene-level omnibus tests.** All alleles of a gene are not tested jointly, for two reasons:

- It is not standard in HLA fine-mapping: Hirata et al. (Nat Genet 2019) test amino-acid positions only.
- It depends on the `unlisted` category, the copies of alleles outside the tested matrix. Its share is lower in cases than in controls in every HLA gene (STUDY_NOTES §10), as a platform difference in typing would produce.

Gene-level tests ran once (run 14) and were removed on 2026-10-25; their outputs are kept in `logs/stale_before_run15/`.

**Model.** For a position with residues k = 1 … K, of which `ref` is the most frequent:

```
logit Pr(y_i = 1) = x_iᵀ β + Σ_{k ≠ ref} g_ik γ_k
H0: γ_k = 0 for every k ≠ ref          LRT = 2 (ℓ̂_1 − ℓ̂_0) ~ χ²_df,   df = K − 1
```

- y_i: case (1) or control (0);
- x_i: intercept, SEX, PC1–10, and in the stepwise walk the residue columns of every entered position;
- g_ik: copies (0, 1, 2) of residue k carried by sample i. The reference column is left out because a sample's copies over all residues of a position sum to 2;
- γ_k: log OR per copy of residue k against the reference;
- ℓ̂_0, ℓ̂_1: maximized log-likelihoods without and with the position's columns, on the same samples;
- df: the rank the position's columns add, K − 1 unless a column is collinear with the model.

There is no random effect. As in `ASSOC_FIXED`, relatives are removed rather than modelled.

**Why the fixed model.** It is the standard for HLA omnibus tests (logistic regression with principal components and a likelihood-ratio test; Hirata et al., Nat Genet 2019; Sakaue et al., Nat Protoc 2023), and SAIGE has no multi-degree-of-freedom test. The samples are the case set's unrelated samples, as in `ASSOC_FIXED`, which differ from the random-model samples only by 27–30 related controls. The covariates are sex and 10 PCs.

**Categories.**

- A position's categories are its residues.
- Copies that carry a residue the matrix does not list (2 minus the listed copies) form one more category, `unlisted`: at most 0.6 % of copies.
- Within the cohort's samples called at the position, a category with copy frequency < 1 % joins the pool `rare`. A pool still below 0.5 % (the single-marker floor) joins the reference.
- The reference is the most frequent category, and df = categories − 1. A position with no category left besides the reference is not tested.
- A sample is blank in all of a position's columns or in none, so each position is tested on its called samples (complete case).

**Fitting.** Both models are fitted by maximum likelihood (Newton–Raphson with step halving) on the samples called at the position and at every conditioning position. Copies are additive (0/1/2), as in the single-marker scans. A fit that separates (|log OR| > 15) or does not converge keeps its P but is never significant.

**Threshold.** Bonferroni: 0.05 / the positions tested in the cohort's marginal scan (290–296). GEC applies to single-degree-of-freedom tests only, and no published LD-adjusted method covers multi-degree-of-freedom tests, so the conservative Bonferroni correction is used.

**Stepwise selection** follows the rule of `STEPWISE` (§7.4), with positions in place of markers:

- candidates are the positions significant in at least one cohort;
- the lead is the candidate whose largest P over the three cohorts is smallest;
- the lead's every residue enters the model as covariates;
- the walk stops when no position is a candidate, after at most 5 signals.

A signal's members are the positions significant before its lead entered and in no cohort after. Consistency is the number of cohorts in which a position is significant. An omnibus test has no single direction, so no direction is required.

**Comparison with the single-marker scans.** `OMNIBUS_SUMMARY` sets each position beside two residues: the best single residue it contains in the fixed model, which has the same model and samples so only the unit of testing differs; and the best single residue in the random model, the primary scan. It also records which position holds each random-model signal lead, or holds it among its members.

## 9. Known limits

- **Platform.** All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X (15×). Sequencing platform and depth can affect HLA typing. A marker whose typing differs by platform would look like an association here, and no covariate can separate the two. A platform-sensitivity check of the hits is planned separately; it is not part of this component.
- **Conditioning separates statistically, not functionally.** Markers in near-perfect LD cannot be told apart (collinear rows have no conditional P). A platform effect shared by the lead and the HLA markers is not removed by conditioning.
- **Omnibus tests use the fixed model.** The single-marker scans test each residue on its own, against all other residues at that position. The per-position omnibus tests (§8), the second main analysis, use the fixed model on the unrelated samples and a conservative Bonferroni threshold.
- **Platform within cases.** The platform mix differs between cases with and without DVT (§2.1). This affects the case-only test and the DVT-yes and DVT-no sets alike.
- **Monomorphic markers.** A marker with no minor allele among the analysed samples, possible only in `dvt_case_only`, keeps its row with `errcode` `MONOMORPHIC` and no P. It is not a test (§5.1) unless it varies among the random-model samples.
- **DVT sets are small.** They have about 140–160 cases each, so power is lower than for all cases, and some apparent differences between the sets are expected by chance. DVT status is unknown for about 30 % of cases.
- **Missing-call rules differ between the models.** At call rates of at least 95 % this difference is small; `model_concordance.tsv` quantifies it.
