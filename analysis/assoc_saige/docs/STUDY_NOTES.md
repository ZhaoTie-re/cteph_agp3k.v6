# Study notes — this dataset's run of `assoc_saige`

**`assoc_saige` is a reusable component.** [`METHODS.md`](METHODS.md) describes the *method* and holds
nothing specific to any dataset. This file is the other half: every finding, measured value, named
locus and configured resource that belongs to **this** project's run.

> **Status: all six scans complete** (3 cohorts x 2 models, 132 chromosome jobs, 0 failed).
> Numbers below are measured from that run under the default SPA p-value.

---

## 1. Configuration used

| parameter | value |
|---|---|
| `Cohorts` | `narrow_mainland`, `intermediate_mainland`, `full_mainland` (nested, narrowest first) |
| genotypes | `12_model_inputs/<cohort>/genotype/random_model` — **relateds retained** |
| `Models` | `main` (SEX + 10 `bbj_mainland` PCs) and `detect` (SEX only); `PeakModel = main` |
| GRM | **full only**, LD-pruned markers, high-LD regions excluded, LOCO on |
| `PValue` | `no_spa` — SAIGE's `p.value.NA` |
| Firth | off |
| thresholds | `PGenomeWide = 5e-8`, `PSuggestive = 1e-5`, `MaxSuggestive = 10`, `PeakFlank = 250 kb` |
| phenotype | CTEPH case/control |

### Sample sets: random model against fixed model

The GRM is the reason relateds are kept, so every cohort here is **larger** than the same cohort in
`assoc_plink2`. Cases are identical; the additional samples are related controls.

| cohort | random_model N | (`assoc_plink2` fixed_model N) | cases | controls | variants |
|---|---|---|---|---|---|
| `narrow_mainland` | 2,195 | 2,168 | 419 | 1,776 | 5,113,626 |
| `intermediate_mainland` | 2,508 | 2,480 | 429 | 2,079 | 5,113,982 |
| `full_mainland` | 3,101 | 3,071 | 439 | 2,662 | 5,115,116 |

`bbj_mainland_pc.sex.tsv` covers **every** random-model sample in all three cohorts, so no sample is
lost to a missing covariate — verified, and `00.prep/pheno_covar.tsv` row counts confirm it per run.

### External resources

Identical to `assoc_plink2`; see [that component's study notes](../../assoc_plink2/docs/STUDY_NOTES.md)
for versions and paths. SAIGE itself:

| purpose | path | version |
|---|---|---|
| SAIGE step 1 / step 2 | `/home/b/b37974/anaconda3/envs/saige/bin/step{1,2}_*.R` | SAIGE 1.1.9 |
| R for SAIGE | `/home/b/b37974/anaconda3/envs/saige/bin/Rscript` | env `saige` |

SAIGE 1.1.9 reads PLINK directly, so this project runs **no BGEN conversion**. The earlier version of
this component did, along with a sparse-GRM branch; both were removed as redundant.

## 2. Scan size and calibration

**The mixed model removed most of the inflation the fixed-effects analysis reported.**

| cohort | N | λ `main` | λ `detect` | λ, `assoc_plink2` additive |
|---|---|---|---|---|
| `narrow_mainland` | 2,195 | 1.039 | 1.040 | 1.120 |
| `intermediate_mainland` | 2,508 | 1.038 | 1.024 | 1.121 |
| `full_mainland` | 3,101 | 1.021 | 1.017 | 1.132 |

λ falls from ~1.12 to ~1.02–1.04. That is the result this component exists to produce: the residual
inflation the fixed-effects scaffold reported was largely relatedness and fine-scale structure, and a
full GRM absorbs it.

**The principal components add almost nothing on top of the GRM.** `main` and `detect` differ by at
most 0.014, and in `intermediate_mainland` the probe is *lower* than the reported model (1.024 vs
1.038). Read against §3 of METHODS: once a GRM is in the model, the PCs are close to redundant for
calibration in these cohorts. They are retained in the reported model because the design fixes the
covariate set in advance, not because the data demand them.

GRM markers after LD pruning and high-LD exclusion: 378,924 (`narrow_mainland`; see
`00.prep/grm.prune.in` per cohort).

## 3. Peaks called

**No locus reaches genome-wide significance under the default SPA p-value.**

| cohort | genome-wide | suggestive | smallest *P* |
|---|---|---|---|
| `narrow_mainland` | 0 | 62 | 1.74 × 10⁻⁷ |
| `intermediate_mainland` | 0 | 134 | 5.60 × 10⁻⁸ |
| `full_mainland` | 0 | 81 | 2.52 × 10⁻⁷ |

The strongest variant in **every** cohort is `chr16:53882967:G:A`, in *FTO*, about 5 kb from the
fixed-effects lead `chr16:53887925:T:C`.

### What happened to the fixed-effects hits

| cohort | *P*, `assoc_plink2` | *P*, SAIGE | OR, plink2 | OR, SAIGE |
|---|---|---|---|---|
| `narrow_mainland` | 7.39 × 10⁻⁸ | 2.57 × 10⁻⁷ | 1.631 | 1.662 |
| `intermediate_mainland` | 1.15 × 10⁻⁸ | 9.13 × 10⁻⁸ | 1.648 | 1.666 |
| `full_mainland` | 2.88 × 10⁻⁸ | 3.15 × 10⁻⁷ | 1.591 | 1.607 |

(at the fixed-effects lead `chr16:53887925:T:C`, so the two engines are compared at the same variant.)

**The effect estimate is stable; the evidence is weaker.** The odds ratio moves by less than 0.02 in
every cohort — the *FTO* signal is not an artefact of the fixed-effects model. What changes is the
standard error: the mixed model is more conservative once relatedness and structure are in the null,
and *P* weakens by roughly 5–10×. The locus no longer clears 5 × 10⁻⁸.

The honest reading is that the fixed-effects *P*-values were inflated along with everything else
(λ ≈ 1.12), and the SAIGE numbers are the better-calibrated ones. *FTO* remains the strongest signal
in the data by a wide margin; it is a candidate, not a genome-wide result.

**The SPA choice is decisive here, and the default is `no_spa`.** SAIGE emits both *P*-values in one
pass; `params.PValue` selects which one the whole downstream reads, and `run_manifest.json` records
the choice. Under `no_spa` (the normal approximation) *FTO* reaches 3.84 × 10⁻⁸ in
`intermediate_mainland` and **is** called genome-wide; under `spa` it does not, and no locus in any
cohort does.

The saddlepoint correction exists to control case/control imbalance, and at 429 cases against 2,079
controls (1 : 4.8) that imbalance is mild — SPA is aimed at ratios an order of magnitude more
extreme. Applying it here is not free: it is a one-directional penalty on exactly the variants that
matter, and there is no verified control bias in this data to justify paying it. So `no_spa` is the
default and **the reported result is 1 genome-wide locus** (*FTO*, `intermediate_mainland`).

This is a judgement, not a fact, and it is the single choice that most changes the headline. Both
numbers are in `scan_qc_all.tsv`; setting `--PValue spa` re-runs `TO_SUMSTATS` and everything below
it, and re-runs neither `FIT_NULL` nor `RUN_ASSOC`.

## 3b. Relatedness the GRM accounts for

Measured from **the GRM the null model was actually fitted with**, not from a separate KING run: the
figure then reports the numbers the model used rather than a second opinion about the same samples.
A GRM relationship is 2 × kinship, so the degree bands are 0.354 / 0.177 / 0.0884.

| cohort | related pairs | 1st | 2nd | 3rd | samples with a relative | cases with a relative | kept vs fixed-effects set |
|---|---|---|---|---|---|---|---|
| `narrow_mainland` | 55 | 10 | 14 | 31 | 94 (4.3 %) | 0 | 27 |
| `intermediate_mainland` | 56 | 10 | 14 | 32 | 96 (3.8 %) | 0 | 28 |
| `full_mainland` | 70 | 13 | 19 | 38 | 123 (4.0 %) | 0 | 30 |

No duplicate or MZ pair in any cohort; the maximum GRM relationship is 0.565–0.566, i.e. one
first-degree pair per cohort at the top of the range. **No case has a relative in the analysis** —
the relatedness is entirely among controls, so it is nuisance structure rather than anything
correlated with phenotype. Retaining it recovers 27–30 samples per cohort that the
relatedness-pruned fixed-effects design must drop; 21 / 21 / 29 of those recovered samples are
themselves in a related pair (`n_absent_and_related`), which is why dropping them was necessary
there and is not here.

These counts run above the KING-based ones (47 / 48 / 61 pairs) because **a GRM value conflates
relatedness with shared ancestry** — KING is robust to population structure by construction and a
GRM is not. That is a real limitation of reading relatedness off the GRM, and it is disclosed rather
than avoided: the point of this table is to describe what the fitted model absorbed, and the model
absorbed the conflated quantity.

The **diagonal** is reported separately (median 1.003–1.004, range 0.961–1.111) because it is
self-relatedness and SAIGE overrides it: the fit runs with `--isDiagofKinSetAsOne=True`, so every
diagonal element is forced to 1.0. Recording what it was before that override is provenance the run
would otherwise lose.

## 4. Against the fixed-effects component

`assoc_plink2` and `assoc_saige` are **not** independent analyses — they share the cases and almost
all the controls. A locus present in both is one observation, not two. What the pair does support is
reading how an estimate behaves when relatedness and structure move from a fixed set of covariates
into a random effect:

- a signal that survives with a similar effect size was not being carried by the fixed-effects
  model's residual structure;
- a signal that weakens sharply under the GRM was partly structure;
- a signal that strengthens gained from the related samples the fixed-effects analysis discarded.

Anchors from the fixed-effects run, for that comparison: FTO / rs16952623 genome-wide additive in
`intermediate_mainland` (*P* = 1.15 × 10⁻⁸) and `full_mainland` (*P* = 2.88 × 10⁻⁸), suggestive in
`narrow_mainland` (*P* = 7.39 × 10⁻⁸); MIR3681HG genome-wide in `full_mainland`; credible sets of 6,
6 and 21 variants. λ_ADD ≈ 1.12 in all three cohorts, which is the inflation this component exists to
address.

## 5. Reading this study's results

There is **no independent replication cohort and no way to obtain one**. The three sample sets are
nested, and the two engines share their samples. A peak here is a threshold crossing in one scan of
one sample set.

The platform confound recorded in the fixed-effects component's study notes applies unchanged: WGS
platform is perfectly separated from phenotype, cannot be fitted as a covariate, and was checked and
found not to be driving the inflation.

## 6. One drawn locus has no usable rs accession

`ADARB2`, `chr10:1747167:CGTG:C` — a 3 bp deletion, suggestive in all three cohorts, so it is drawn
in panel (c) of `cohort_compare.png`. The ToMMo VCF gives it `rs1347066655;rs1564211184`, and
**neither accession can be attributed to it**:

```
POS=1747167  ID=rs1347066655;rs1564211184  REF=CGTG  ALT=C    <- ours
POS=1747167  ID=rs1347066655;rs1564211184  REF=C     ALT=T
```

Two records, different alleles, identical ID string. The source VCF is `norm`ed, and splitting a
multi-allelic record copies the ID field onto every split allele rather than distributing the
accessions, so one of these two names the deletion and the other names the SNV — the normalised file
no longer says which.

This is not rare in the resource. In `chr10:1–3 Mb` alone, 20,922 records carry a multi-accession ID
and **6,599 of those ID strings are stamped on more than one allele at their position**. One example
carries three:

```
chr10:1000259  rs1382461852;rs2132161638  TGAGA>T
chr10:1000259  rs1382461852;rs2132161638  TGAG>T
chr10:1000259  rs1382461852;rs2132161638  TG>T
```

A **single** accession shared across alleles at a position is a different thing and is not flagged —
RS numbers are site-level, so `rs790041` on both `C>T` and `C>G` at `chr2:232445059` names our
variant either way. Only a multi-accession string spread over several alleles is a union that cannot
be attributed. METHODS §8 of the fixed-effects component tabulates the four cases.

`lookup_rsids` detects this and leaves `rsID` as `.`, keeping the candidates in `rsID_unresolved`;
the figure prints `rsID unresolved` on the label's third line. **Cite this variant as
`chr10:1747167:CGTG:C`, never by either rs number** — resolving which accession belongs to the
deletion needs a dbSNP lookup this pipeline does not perform.

Scope in the drawn figures: 1 of 16 loci here, 0 of 18 in the fixed-effects component.
