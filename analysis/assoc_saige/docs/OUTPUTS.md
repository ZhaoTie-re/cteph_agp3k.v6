# Outputs — `assoc_saige`

Every path below is under `results/`. Column definitions that are shared with the fixed-effects
component are the same columns, produced by the same code in `_shared/scripts`.

```
<cohort>/
  00.prep/
    pheno_covar.tsv         phenotype + every covariate, on the genotype samples
    grm.prune.in            LD-pruned marker list the GRM is built from
    relatedness_summary.tsv one row: pair counts by degree, samples with a relative
    relatedness_pairs.tsv   every pair at kinship >= params.MinKinship
    relatedness_hist.tsv    the GRM off-diagonal distribution, binned
    relatedness_diag.tsv    the GRM diagonal, binned (SAIGE overrides it with 1.0)
  01.assoc/<model>/
    null.rda                SAIGE step-1 null model (the fitted GRM)
    null.varianceRatio.txt  variance ratio step 2 needs
    null.fit.log            step-1 log, kept: it carries the fitted variance components
    saige.assoc.tsv         SAIGE step-2 output, all chromosomes, RAW
    sumstats.tsv            the canonical schema — what everything downstream reads
    saige_extra.tsv         every SAIGE column the canonical schema has no slot for
  02.scan/
    scan_qc.tsv             per model: N, lambda_GC, counts, min P
  03.peaks/
    peaks.tsv               peak-model peaks, two tiers
    model_peaks.tsv         EVERY model's peaks, two tiers
    lead_variants.tsv       one lead per peak — the fan-out key
    lead_annotation.tsv     annotated peak-model leads
    model_peaks_annotation.tsv   annotated peaks of every model
    <tier>/<peak_id>/       per-peak follow-up, peak model only
                            <tier> = genome_wide | suggestive
  figures/
    01.scan/   scan.<model>.png + .md
    02.regional/<tier>/ 03.finemap/<tier>/ 04.conditional/<tier>/
_comparison/
  tables/    scan_qc_all.tsv · peaks_all.tsv · lead_annotation_all.tsv
             model_peaks_all.tsv · cs_variants_all.tsv
             lead_crosscohort.tsv   every lead, in every cohort
             lead_bymodel.tsv       every genome-wide lead, under every model
             relatedness_all.tsv    the per-cohort relatedness summary, stacked
  figures/   cohort_compare.png · cohort_manhattan.png · model_compare.png
             relatedness.png        (+ .md each)
_run_info/   trace · report · timeline · dag · run_manifest.json
```

**Every PNG has a companion `.md`** of the same name, written by `_shared/scripts/figure_doc.py`.

## `01.assoc/<model>/sumstats.tsv` — the canonical schema

The one file the rest of the pipeline reads. Column names and meanings are **plink2's `--glm`
schema**, so a scan from either engine is interchangeable here.

| column | meaning |
|---|---|
| `#CHROM` `POS` | position on the analysis build |
| `ID` | `chr:pos:REF:ALT`, **reconstructed** from the alleles — the join key for every lookup |
| `REF` `ALT` | reference and alternate allele |
| `A1` | the **effect** allele; SAIGE's `BETA` is per `Allele2`, so `A1 = ALT` |
| `A1_FREQ` | frequency of `A1` among analysed samples |
| `TEST` | `ADD` — this component fits additive only |
| `OBS_CT` | **non-missing SAMPLE count**, not an allele count (plink2 uses the name both ways) |
| `OR` `LOG(OR)_SE` | `exp(BETA)` and `SE` |
| `L95` `U95` | 95 % CI, as **two numeric columns**; never one packed string |
| `Z_STAT` | `BETA / SE` |
| `P` | `p.value` or `p.value.NA`, per `params.PValue` — see METHODS §3 |
| `ERRCODE` | `.`, or `SAIGE_NO_FIT` where the row carries no usable estimate |

## `01.assoc/<model>/saige_extra.tsv` — everything else SAIGE reported

Keyed by the same `ID`, so it joins to `sumstats.tsv` directly. Nothing SAIGE produced is discarded.

| column | why it is kept |
|---|---|
| `MarkerID_saige` | SAIGE's own name for the variant — what `--condition` must be given |
| `Is.SPA` | whether the saddlepoint correction was applied to this variant |
| `p.value.NA` | the other p-value, so `params.PValue` is reversible without re-running SAIGE |
| `AF_case` `AF_ctrl` | per-group frequencies straight from SAIGE |
| `N_case` `N_ctrl` | the counts `OBS_CT` is derived from |
| `N_case_hom` `N_case_het` `N_ctrl_hom` `N_ctrl_het` | per-group genotype counts — a free cross-check on the ones the annotator computes independently |
| `Tstat` `var` | the score statistic and its variance |
| `MissingRate` `AC_Allele2` | per-variant missingness and allele count |

## `02.scan/scan_qc.tsv`

One row per model: sample sizes, `lambda_gc`, the counts behind it, and the smallest *P*. λ is
**reported, never applied**.

## `_comparison/tables/lead_bymodel.tsv`

One row per (**locus representative** × **cohort** × **model**) — the table behind
`model_compare.png`. The loci are panels (b) and (c) of `cohort_compare`, selected on
`lead_crosscohort.tsv`'s `panel` column, so both figures cover the same set.

| column group | columns |
|---|---|
| identity | `label` (gene or variant), `variant_id`, `chrom`, `pos` |
| context | `cohort`, `model` |
| effect | `OR`, `L95`, `U95`, `P` |

The leads are those the **peak model** called; the probe contributes its estimate at the same
variants, never its own peak list.

## `00.prep/relatedness_*.tsv`

`relatedness_summary.tsv` is one row per cohort: sample and pair counts, the kinship median and
maximum, pairs in each KING degree class, how many samples carry at least one relative (and how many
of those are cases), and — when `params.CompareGenotypeDir` is set — how many analysed samples are
**absent** from the relatedness-pruned set, which is the cost a fixed-effects design pays.

`relatedness_pairs.tsv` is one row per related pair (`id1`, `id2`, `kinship`, `degree`).
`relatedness_hist.tsv` bins the whole pairwise distribution so the figure never carries the matrix.

Computed on the LD-pruned GRM markers, so these describe the matrix the null model fitted.

## `_run_info/run_manifest.json`

Every parameter the run resolved, including the two that cannot be read off the figures:

- `p_value_source` / `p_value_column` — which SAIGE p-value the whole downstream used;
- `model_covariates` — the exact covariate list fitted for each model.

Also records that `sparse_grm` is `false` and `loco` is `true`, so a later reader does not have to
infer either from the code.
