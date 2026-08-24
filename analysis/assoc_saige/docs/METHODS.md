# Methods — `assoc_saige`

What the component does and why, independent of any dataset. Everything measured on a particular
dataset — counts, λ values, named loci, configured resources — belongs in
[`STUDY_NOTES.md`](STUDY_NOTES.md), not here. Values written as `params.X` are configuration.

The peak rule, the annotation, the LD design, the fine-mapping and every figure are **identical to
[`assoc_plink2`](../../assoc_plink2/docs/METHODS.md)** and are implemented by the same code in
`_shared/scripts`. That is by construction, not coincidence: it is what makes a mixed-model scan and a
fixed-effects scan of the same cohort comparable row for row. This document states the method in full,
and marks each section that is shared.

---

## 1. What this component is

The **random-effect** analysis. `assoc_plink2` is a fixed-effects scaffold whose own calibration
read-out points here: residual λ > 1 after global PCs is fine-scale structure that a fixed set of
components cannot absorb, and a genetic relationship matrix is the standard way to absorb it.

Two things follow from using a GRM, and they are the whole difference from the fixed-effects
component:

1. **Relateds are retained.** The GRM accounts for relatedness, so removing related samples first
   would discard information for no reason. The random-model genotypes are therefore used, and the
   analysed sample sets are larger than the fixed-model ones.
2. **A calibration probe becomes possible.** With a GRM in the model one can ask what the principal
   components are still doing, by fitting the same samples and the same GRM without them (§4).

## 2. Association model

SAIGE fits a logistic mixed model,

```
logit Pr(case_i) = β₀ + β g_i + Σ_j γ_j X_{j,i} + b_i ,     b ~ N(0, τ K)
```

where `K` is the **full** genetic relationship matrix built from LD-pruned markers and `X` is the
covariate set of the model being fitted (§4). `g_i` is the additive allele count — this component
fits **additive only**; the dominant and recessive codings are the fixed-effects component's
business.

```
step1_fitNULLGLMM.R  --plinkFile=<pruned markers>  --traitType=binary
                     --covarColList=<the model's covariates>
                     --useSparseGRMtoFitNULL=FALSE  --LOCO=TRUE
                     --isDiagofKinSetAsOne=True
                     --numRandomMarkerforVarianceRatio=<params.VarianceRatioMarkers>

step2_SPAtests.R     --bedFile/--bimFile/--famFile  --AlleleOrder=ref-first  --chrom=<N>
                     --GMMATmodelFile  --varianceRatioFile
                     --is_Firth_beta=FALSE  --LOCO=TRUE  --is_output_moreDetails=TRUE
```

**Full GRM only.** A sparse GRM approximates `K` by discarding relationships below a cutoff. It is
the right trade at biobank scale, where a dense `K` is intractable; at a few thousand samples the
dense matrix is affordable and the approximation buys nothing but a threshold to defend.

**No Firth.** Same reasoning as the fixed-effects component: Firth penalisation changes the estimand,
so a penalised and an unpenalised effect are not on one scale and cannot share a table.

**LOCO in both steps.** Leave-one-chromosome-out excludes the tested chromosome from the GRM, so a
variant does not contribute to the random effect it is being tested against. Without it a real signal
is partly absorbed by the very matrix meant to absorb only structure. Step 2 requires `--chrom` when
LOCO is on, which the per-chromosome layout provides anyway.

**GRM markers** are LD-pruned (`params.Prune*`) with the known high-LD regions excluded
(`params.HighLdFile`). The GRM describes relatedness and broad structure, so correlated blocks and
inversion-prone regions would let a handful of loci dominate it.

## 2b. Reproducibility

**Nothing in this component samples without a seed**, so a re-run on the same inputs reproduces the
same fit:

| step | is it stochastic? |
|---|---|
| GRM marker selection (`plink2 --indep-pairwise`) | no — deterministic given the same genotypes |
| the GRM itself | no — built from **every** pruned marker; nothing is subsampled |
| variance-ratio markers | yes, but `step1_fitNULLGLMM.R` calls `set.seed()` before drawing them |
| step 2 score tests, SPA | no |

SAIGE exposes no option to override its internal seed, so `params.RandomSeed` **records** the value
rather than setting it; it is written to `run_manifest.json` for provenance. The marker list the GRM
was built from is published as `00.prep/grm.prune.in`, so the matrix is reconstructible exactly.

## 2c. The relatedness the GRM accounts for

The design claim of this component is that related samples are **kept** and the GRM accounts for
them. That claim is measured rather than asserted: `RELATEDNESS` computes the **variance-standardized
relationship matrix** (`plink2 --make-rel`) **on exactly the LD-pruned markers the GRM is built
from**. That is the same quantity SAIGE builds, so the read-out reports the values the model used
rather than a proxy for them.

Two parts are reported separately because they mean different things:

- **off-diagonal** — relatedness between two samples. Degree bands are the expected relationships:
  0.5 (1st), 0.25 (2nd), 0.125 (3rd), ≈ 1.0 (duplicate/MZ), i.e. twice the corresponding KING kinship
  boundaries because a GRM relationship is 2 × kinship. `params.MinRelationship` sets what counts as
  related.
- **diagonal** — a sample's relationship with itself, 1 + inbreeding. It is reported because SAIGE
  **overrides** it (`--isDiagofKinSetAsOne=True`), so this is the only record of what the data said.

**Limitation, disclosed rather than avoided:** a GRM value conflates relatedness with shared ancestry.
These cohorts are ancestry-filtered subsets of one another, so the off-diagonal bulk shifts with the
filter and the degree bands are indicative rather than diagnostic. A structure-insensitive estimator
such as KING-robust would classify degrees more reliably; what is gained in exchange is that this
reports the matrix the model actually fitted.

The read-out also states the cost the alternative design pays: how many analysed samples are absent
from `params.CompareGenotypeDir` (the relatedness-pruned set the fixed-effects component uses).
A duplicate or MZ pair appearing here is a sample-handling question, not a biological one.

## 3. Which p-value, and why the choice is free

SAIGE reports two p-values for every variant in a single pass:

| column | what it is |
|---|---|
| `p.value` | saddlepoint-approximated where `Is.SPA` is true, else the normal approximation |
| `p.value.NA` | the normal approximation, always |

Score tests on unbalanced case/control data have an inflated tail under the normal approximation, and
the saddlepoint approximation corrects it. **The default here is nevertheless `no_spa`**
(`params.PValue = 'no_spa'`, i.e. `p.value.NA`): the correction earns its cost when the imbalance is
severe, and this design's is not judged to be. The alternative remains one parameter away.

**The choice is consequential, not cosmetic.** It decides whether the strongest locus clears
5 × 10⁻⁸, so it is recorded in `run_manifest.json` and stated in the study notes rather than left to
be inferred. Anyone reading a result from this component should know which p-value produced it.

Because both are produced by the same run, **the choice is made downstream** — `TO_SUMSTATS` selects
the column when it writes the canonical schema (§5). Switching `params.PValue` re-runs the adapter and
the figures and never re-runs SAIGE; this is verified, not assumed. Both columns are always kept in `saige_extra.tsv`, so the choice
is reversible at no cost. Which one a run used is recorded in `_run_info/run_manifest.json`
(`p_value_source`, `p_value_column`); it is deliberately not printed on the figures.

## 4. Two covariate models, one of which is a result

`params.ModelCovars` defines the models; `params.PeakModel` names the one that is reported.

| model | covariates | what it is for |
|---|---|---|
| `main` | `SEX` + `params.NPcs` PCs | the reported analysis; defines every peak |
| `detect` | `SEX` only | a calibration probe |

The probe fits the **same samples** and the **same GRM** with the principal components removed. The
difference between the pair is therefore a direct read-out of what the PCs contribute *on top of* a
GRM — not a second analysis, and not a model-selection exercise:

- if the probe's λ is already near 1, the GRM has absorbed the structure by itself;
- if the probe is inflated while the reported model is not, the PCs are removing structure the GRM
  cannot reach;
- a lead whose effect is stable across the pair is not an artefact of the covariate choice.

**The probe defines no peaks and nothing fans out from it** — no fine-mapping, no conditional
analysis, no regional plot — exactly the treatment the fixed-effects component gives its non-primary
genetic models. Reporting it as a second result would be reporting the same data twice.

`params.PeakModel` is the only model that defines the fan-out, so the whole downstream is unchanged
from the fixed-effects component, which already implements "peaks for every model, follow-up from one".

## 5. The canonical summary-statistic schema

**The canonical schema is plink2's `--glm` schema**, and `scripts/saige_to_sumstats.py` adapts SAIGE
to it. Everything downstream then reads one format regardless of engine, and the two components'
results are directly comparable.

| canonical | from SAIGE |
|---|---|
| `#CHROM` `POS` | `CHR` `POS` |
| `ID` | **reconstructed** `chr:pos:REF:ALT` from `CHR`/`POS`/`Allele1`/`Allele2` |
| `REF` `ALT` `A1` | `Allele1` `Allele2` `Allele2` — SAIGE's `BETA` is per `Allele2` |
| `A1_FREQ` | `AF_Allele2` |
| `OBS_CT` | `(N_case + N_ctrl) × (1 − MissingRate)`, rounded |
| `OR` `LOG(OR)_SE` | `exp(BETA)` `SE` |
| `L95` `U95` | `exp(BETA ∓ 1.96 · SE)` |
| `Z_STAT` | `BETA / SE` |
| `P` | `p.value` or `p.value.NA` (§3) |
| `ERRCODE` | `.`, or `SAIGE_NO_FIT` where the row carries no usable estimate |

Two details are load-bearing:

- **`OBS_CT` is a SAMPLE count.** plink2 uses that name for a non-missing *sample* count under
  `--glm` and for an *allele* count under `--freq`. The downstream reads the `--glm` convention.
- **The variant ID is reconstructed, not taken from `MarkerID`.** It is the join key for every
  downstream lookup — rsID, the snpEff index, the LD panels, the credible sets — and what SAIGE puts
  in `MarkerID` depends on how the genotypes were exported. The adapter rebuilds it from the alleles
  and cross-checks the result against the `.bim`, reporting the match rate, so a convention mismatch
  is loud rather than silent.

Everything SAIGE knows that the canonical schema has no slot for — `Is.SPA`, `p.value.NA`, `AF_case`,
`AF_ctrl`, `N_case`, `N_ctrl`, the per-group genotype counts — is written beside it in
`saige_extra.tsv` and never dropped.

## 6. Which rows count as a result  *(shared)*

`usable_mask()` in `_shared/scripts/call_peaks.py`, unchanged: `ERRCODE == '.'`, `OR` and `SE` finite
and positive, `0 < P ≤ 1`. The adapter sets `ERRCODE` so that anything failing the numeric guard is
labelled rather than silently dropped.

## 7. Peaks — two tiers, one merge  *(shared)*

Called **once** at `params.PSuggestive` and then tiered on `params.PGenomeWide`, because the
thresholds are nested. Merging is by **distance** (`params.PeakFlank`), not LD clumping, so the peak
definition needs no reference panel and cannot shift with the choice of LD sample. Lead = smallest
*P* in the peak under the peak model.

**Both tiers receive the full per-peak follow-up** — LD, fine-mapping and conditional analysis — with
outputs separated by tier (`03.peaks/<tier>/`, `figures/0{2,3,4}.*/<tier>/`) so the two are never
mixed. The conditional analysis is judged against the threshold that **defined** the peak: a
suggestive peak tested against `PGenomeWide` would fail at round 0 every time and report zero signals
everywhere, which is meaningless rather than empty.

What separates the tiers is now interpretation, not machinery. The suggestive tier still describes
the shape of the scan rather than listing findings; the follow-up exists so that shape can be
inspected, not so that every suggestive locus can be reported as a result.

## 8. Calibration  *(shared rule, extra read-out here)*

λ_GC is reported and **never applied**, as everywhere else in this project. What is new here is that
the probe (§4) turns λ from a single number into a comparison: the pair (reported, probe) says how the
structure is being absorbed, not merely how much of it is left.

A mixed model does not guarantee λ = 1. Residual inflation after a full GRM is either polygenicity or
structure the GRM cannot express, and this component does not attempt to distinguish them — that
needs LDSC, which is unusable at a few thousand effective samples.

## 9. Lead annotation, gene models, per-peak follow-up  *(shared)*

Identical to the fixed-effects component and implemented by the same code: the lead table with `OR`,
`L95`, `U95` as three numeric columns; per-group genotype counts, EAF, missingness and HWE from two
plink2 calls; rsID from `params.rsid_vcf`, attributed only when the accession is resolved to the
variant's own allele and otherwise held in `rsID_unresolved` (see the fixed-effects component's
METHODS §8 — one drawn locus here is affected, recorded in [STUDY_NOTES.md](STUDY_NOTES.md) §6);
the gene rule in `_shared/scripts/gene_utils.R`; SuSiE
with an in-sample LD matrix; LD from three sources on one binned scale, where a reporting floor is a
**bound**, not missing data.

**Conditional analysis is the one piece that differs.** The loop is the same — round 0 unconditioned,
each later round adds the previous round's top variant, stop when nothing in the window clears the
threshold — but it runs SAIGE step 2 with `--condition` and **re-uses the fitted null model and its
variance ratio**. The GRM is never re-fitted for a conditioning round: it would cost hours per round,
and it would change the null the rounds are compared against, making them incomparable. The variant
fed to `--condition` is the `MarkerID` SAIGE itself printed, so the identifier convention is
self-consistent by construction.

## 10. Multiplicity, cross-cohort reporting, figures  *(shared)*

`params.PGenomeWide` applied identically to every cohort, with no further adjustment across nested
sample sets. Every peak lead reported in every cohort, with `called_peak` stating explicitly what each
cohort made of it rather than leaving a blank. **Locus identity** — proximity grouping split only where
a cohort's conditional analysis reported two independent signals, and one representative variant per
locus — is documented in [the fixed-effects component's METHODS, §12a](../../assoc_plink2/docs/METHODS.md#12a-locus-identity--what-counts-as-the-same-signal); the code is shared, and here it merges 52 lead
variants into 47 loci. Figures follow
[`_shared/docs/FIGURES.md`](../../_shared/docs/FIGURES.md); every PNG has a companion `.md`.

One figure is specific to this component: **`model_compare.png`**, the reported model against its
probe — λ and the QQ curves in one row, and in the other the OR under both models for **exactly the
loci `cohort_compare` draws**: panels (b) and (c) of that figure, read from the shared `panel` column
so the two figures can never drift apart. Here that is 16 loci — 1 genome-wide somewhere, 15
suggestive in all three cohorts — not the genome-wide leads alone.

## 11. What is hash-bearing

`FIT_NULL` and `RUN_ASSOC` dominate the cost. Their task hashes are computed from their script blocks
and their input files, so their scripts interpolate **only** what they genuinely need: the covariate
list, the staged input names, the SAIGE entry points, and the thread count. Everything presentational
— labels, thresholds, the p-value choice, figure parameters — is consumed downstream of
`TO_SUMSTATS`, which is why none of it can trigger a re-scan.

**No comment may be added inside a `script:` block.** Nextflow hashes that text, and `//` is not a
shell comment — it becomes part of the command and both changes the hash and breaks the script.
