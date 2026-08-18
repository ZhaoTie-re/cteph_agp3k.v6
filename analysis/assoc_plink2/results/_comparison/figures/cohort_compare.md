# Additive scan across three nested sample sets

**Figure file:** `cohort_compare.png`

## The question this figure answers

How do calibration and effect size move as the ancestry filter is relaxed from narrow to full?

## Panels

**(a) Calibration against effective size**

lambda_GC (y) against N_eff (x), one marker per scan: colour = cohort, shape = model. A grey line joins each model's three cohorts in nesting order. This replaces a dual-axis stem plot, where lambda sat on the left axis and N_eff on the right and the reader had to align two series by eye; plotting one against the other makes the trade-off a direction on the page, and leaves room for all nine scans rather than the additive three.

**(b) Effect concordance across all three cohorts**

A forest plot over (genome-wide lead variant) x (cohort). Every variant that reached genome-wide significance in *any* cohort is shown in *all three*, because a nested design always has an estimate in the larger sets. Three states are distinguished, not two: filled diamond = genome-wide in that cohort, filled circle = suggestive there, open circle = not significant there. A row with no estimate at all is annotated "not in this cohort's call set". A blank would have been indistinguishable from a missing estimate, which is why the earliest version of this panel was misleading.

## Interpretation

The cohorts are nested — narrow within intermediate within full — so they share most of their samples and are *not* replication of one another; agreement between them is expected and carries little independent information, while a signal present only in the narrowest set is a candidate for an ancestry-driven artefact. Panel (b) is therefore a description of how one estimate behaves as the sample grows, not a replication test: an interval that narrows while the point estimate holds is a variant gaining sample size, and an estimate that moves toward 1 as the filter relaxes is what a structure-driven signal does. Panel (a) shows why narrow is the most exposed set: relaxing the filter from narrow to full adds 20 cases but 883 controls, so the narrow cohort keeps 95% of the cases against only 67% of the controls. No cohort is designated correct here; a fixed-effects scan cannot settle that, and the GRM-based random-effect follow-up is what will. Dominant and recessive calibration is not in this figure — it is in scan_qc_all.tsv and on each scan's own figure.

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | narrow, intermediate, full |
| lambda_GC additive, narrow | 1.12 |
| lambda_GC additive, intermediate | 1.121 |
| lambda_GC additive, full | 1.132 |
| N_eff, narrow | 1,352 |
| N_eff, intermediate | 1,419 |
| N_eff, full | 1,505 |
| N cases / controls, narrow | 419 / 1,749 |
| N cases / controls, intermediate | 429 / 2,051 |
| N cases / controls, full | 439 / 2,632 |
| genome-wide additive peaks, narrow | 0 |
| genome-wide additive peaks, intermediate | 1 |
| genome-wide additive peaks, full | 2 |
| suggestive additive peaks, narrow | 40 |
| suggestive additive peaks, intermediate | 37 |
| suggestive additive peaks, full | 51 |
| distinct genome-wide lead variants | 2 |
| forest rows | 6 |
| lead variants reported in every cohort | 84 |
| rows in lead_crosscohort.tsv | 252 |

## Full statistics

**Calibration and size, all nine scans**

| cohort | model | n_case | n_ctrl | n_eff | n_analysed | lambda_gc | n_genomewide | n_suggestive |
|---|---|---|---|---|---|---|---|---|
| narrow_mainland | additive | 419 | 1,749 | 1,352 | 5,113,526 | 1.12 | 0 | 130 |
| narrow_mainland | dominant | 419 | 1,749 | 1,352 | 5,092,286 | 1.059 | 0 | 131 |
| narrow_mainland | recessive | 419 | 1,749 | 1,352 | 4,675,710 | 0.6699 | 0 | 162 |
| intermediate_mainland | additive | 429 | 2,051 | 1,419 | 5,113,902 | 1.121 | 4 | 216 |
| intermediate_mainland | dominant | 429 | 2,051 | 1,419 | 5,095,406 | 1.06 | 1 | 191 |
| intermediate_mainland | recessive | 429 | 2,051 | 1,419 | 4,725,886 | 0.665 | 1 | 198 |
| full_mainland | additive | 439 | 2,632 | 1,505 | 5,115,079 | 1.132 | 3 | 204 |
| full_mainland | dominant | 439 | 2,632 | 1,505 | 5,100,619 | 1.065 | 0 | 217 |
| full_mainland | recessive | 439 | 2,632 | 1,505 | 4,802,496 | 0.6543 | 7 | 169 |

**Genome-wide leads — complete statistics in all three cohorts**

| variant_id | rsID | Gene | cohort | called_peak | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Case_Missing_Rate | Case_HWE_P | Control_Genotype_Distribution | Control_EAF | Control_Missing_Rate | Control_HWE_P | A1_FREQ | OBS_CT | N_case | N_ctrl |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2:11891017:G:A | rs66898022 | MIR3681HG | narrow_mainland | not_a_peak | A | G | 2.691 | 1.669 | 4.338 | 4.846e-05 | 378/32/0 | 0.03902 | 0.02148 | 1 | 1656/48/0 | 0.01408 | 0.02573 | 1 | 0.01892 | 2,114 | 419 | 1,749 |
| chr2:11891017:G:A | rs66898022 | MIR3681HG | intermediate_mainland | suggestive | A | G | 2.893 | 1.833 | 4.564 | 4.993e-06 | 387/33/0 | 0.03929 | 0.02098 | 1 | 1945/55/0 | 0.01375 | 0.02487 | 1 | 0.01818 | 2,420 | 429 | 2,051 |
| chr2:11891017:G:A | rs66898022 | MIR3681HG | full_mainland | genome_wide | A | G | 3.31 | 2.157 | 5.08 | 4.342e-08 | 396/33/1 | 0.0407 | 0.0205 | 0.5137 | 2501/67/0 | 0.01305 | 0.02432 | 1 | 0.01701 | 2,998 | 439 | 2,632 |
| chr16:53887925:T:C | rs16952623 | FTO | narrow_mainland | suggestive | C | T | 1.631 | 1.365 | 1.95 | 7.388e-08 | 207/160/39 | 0.2931 | 0.03103 | 0.3381 | 1067/549/72 | 0.2053 | 0.03488 | 0.8815 | 0.2223 | 2,094 | 419 | 1,749 |
| chr16:53887925:T:C | rs16952623 | FTO | intermediate_mainland | genome_wide | C | T | 1.648 | 1.388 | 1.957 | 1.154e-08 | 209/166/40 | 0.2964 | 0.03263 | 0.4106 | 1245/656/84 | 0.2076 | 0.03218 | 0.8915 | 0.2229 | 2,400 | 429 | 2,051 |
| chr16:53887925:T:C | rs16952623 | FTO | full_mainland | genome_wide | C | T | 1.591 | 1.35 | 1.874 | 2.877e-08 | 215/169/41 | 0.2953 | 0.03189 | 0.3526 | 1582/850/121 | 0.2139 | 0.03002 | 0.6376 | 0.2255 | 2,978 | 439 | 2,632 |

## How to read it

1. Do not read agreement between cohorts as replication. They share samples by construction; narrow is a subset of intermediate, which is a subset of full.
2. In (a), read the vertical spread first: the three models separate far more than the three cohorts do, so model choice dominates calibration here.
3. In (b), read down each variant's three rows. An interval that narrows while the point estimate holds is a variant gaining sample size. An estimate that drifts toward 1 as the filter is relaxed is what a structure-driven signal looks like.
4. A peak that appears only in narrow deserves suspicion — the narrow filter retains 95% of cases but only 67% of controls, so it is the set most exposed to ancestry confounding.
5. For every lead of both tiers in every cohort, not just the genome-wide ones, read `_comparison/tables/lead_crosscohort.tsv`.

## What this figure does *not* establish

- It cannot adjudicate between the cohorts. Choosing one requires a model that absorbs fine-scale structure, which this fixed-effects scan is not.
- It is not a replication analysis and no cohort here is independent of another. The component has no external cohort available.
- Panel (b) covers the additive model only. Dominant and recessive peaks are in each cohort's `03.peaks/model_peaks_annotation.tsv` and on their own scan figures.
- lambda_GC in (a) cannot separate confounding from polygenicity; see METHODS §7.

## Symbols

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

- **N_eff** — 4/(1/N_case+1/N_ctrl), the balanced-design equivalent sample size. Reported as a cohort descriptor — it is the size any later meta-analysis or mixed-model run quotes — and it collapses toward 4x the smaller arm as the design becomes unbalanced.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **called_peak** — whether the variant was a peak *in that cohort* (genome_wide / suggestive) or not (not_a_peak). An estimate is reported either way — a cohort that called no peak still has an odds ratio there, and a blank cell would be indistinguishable from a missing one.

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

## Model

```
lambda_GC = median(chi^2),/,chi^2_1,0.5,   chi^2_1,0.5 = 0.4549
N_eff = 4,/,(1/N_case + 1/N_ctrl)
```

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)
