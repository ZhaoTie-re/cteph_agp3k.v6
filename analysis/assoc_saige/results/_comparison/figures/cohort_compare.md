# Additive scan across three nested sample sets

**Figure file:** `cohort_compare.png`

## The question this figure answers

How do calibration and effect size move as the ancestry filter is relaxed from narrow to full?

## Panels

**(a) Calibration against effective size**

lambda_GC (y) against N_eff (x), one marker per scan: colour = cohort, shape = model. A grey line joins each model's three cohorts in nesting order. This replaces a dual-axis stem plot, where lambda sat on the left axis and N_eff on the right and the reader had to align two series by eye; plotting one against the other makes the trade-off a direction on the page, and leaves room for all nine scans rather than the one model's scans alone.

**(b) Effect concordance across all three cohorts**

A forest plot over (genome-wide lead variant) x (cohort). Every variant that reached genome-wide significance in *any* cohort is shown in *all three*, because a nested design always has an estimate in the larger sets. Three states are distinguished, not two: filled diamond = genome-wide in that cohort, filled circle = suggestive there, open circle = not significant there. A row with no estimate at all is annotated "not in this cohort's call set". A blank would have been indistinguishable from a missing estimate, which is why the earliest version of this panel was misleading.

**(c) Leads that hold across every cohort**

The same forest, for leads called SUGGESTIVE in all three cohorts and genome-wide in none. Membership of (b) and (c) is assigned once, in cross_cohort.py, and written to the `panel` column of lead_crosscohort.tsv, so the two panels cannot share a variant and model_compare can draw exactly the same loci. There is no top-N cap: the definition itself bounds the panel, because the variant has to clear the suggestive threshold in every cohort to appear. Each group is labelled with its gene symbol, the representative variant as CHROM:POS:REF:ALT, and that variant's rsID — the three rows are one variant measured three times, and the label says which one. What this panel does NOT show is replication. The cohorts are nested — narrow within intermediate within full — and share every case, so a lead present in all three has survived the ancestry filter, not been confirmed by independent data. An earlier version of this figure carried a fourth panel counting cohorts per lead; it is gone because this panel answers the same question directly and its axis label invited exactly the reading the design forbids.

## Interpretation

The cohorts are nested — narrow within intermediate within full — so they share most of their samples and are *not* replication of one another; agreement between them is expected and carries little independent information, while a signal present only in the narrowest set is a candidate for an ancestry-driven artefact. Panel (b) is therefore a description of how one estimate behaves as the sample grows, not a replication test: an interval that narrows while the point estimate holds is a variant gaining sample size, and an estimate that moves toward 1 as the filter relaxes is what a structure-driven signal does. Panel (a) shows why narrow is the most exposed set: relaxing the filter from narrow to full adds 20 cases but 883 controls, so the narrow cohort keeps 95% of the cases against only 67% of the controls. No cohort is designated correct here; a fixed-effects scan cannot settle that, and the GRM-based random-effect follow-up is what will. Dominant and recessive calibration is not in this figure — it is in scan_qc_all.tsv and on each scan's own figure.

## Values in this rendering

| quantity | value |
|---|---|
| cohorts | narrow, intermediate, full |
| lambda_GC main, narrow | 1.039 |
| lambda_GC main, intermediate | 1.038 |
| lambda_GC main, full | 1.021 |
| N_eff, narrow | 1,356 |
| N_eff, intermediate | 1,422 |
| N_eff, full | 1,507 |
| N cases / controls, narrow | 419 / 1,776 |
| N cases / controls, intermediate | 429 / 2,079 |
| N cases / controls, full | 439 / 2,662 |
| genome-wide main peaks, narrow | 0 |
| genome-wide main peaks, intermediate | 1 |
| genome-wide main peaks, full | 0 |
| suggestive main peaks, narrow | 29 |
| suggestive main peaks, intermediate | 27 |
| suggestive main peaks, full | 32 |
| distinct genome-wide lead variants (b) | 1 |
| leads suggestive in every cohort (c) | 15 |
| forest rows in (b) | 3 |
| lead variants reported in every cohort | 47 |
| rows in lead_crosscohort.tsv | 141 |

## Full statistics

**Calibration and size, all nine scans**

| cohort | model | n_case | n_ctrl | n_eff | n_analysed | lambda_gc | n_genomewide | n_suggestive |
|---|---|---|---|---|---|---|---|---|
| narrow_mainland | main | 419 | 1,776 | 1,356 | 5,113,626 | 1.039 | 0 | 87 |
| narrow_mainland | detect | 419 | 1,776 | 1,356 | 5,113,626 | 1.04 | 0 | 180 |
| intermediate_mainland | main | 429 | 2,079 | 1,422 | 5,113,982 | 1.038 | 2 | 160 |
| intermediate_mainland | detect | 429 | 2,079 | 1,422 | 5,113,982 | 1.024 | 0 | 168 |
| full_mainland | main | 439 | 2,662 | 1,507 | 5,115,116 | 1.021 | 0 | 111 |
| full_mainland | detect | 439 | 2,662 | 1,507 | 5,115,116 | 1.017 | 0 | 121 |

**Genome-wide leads — complete statistics in all three cohorts**

| variant_id | rsID | Gene | cohort | called_peak | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Case_Missing_Rate | Case_HWE_P | Control_Genotype_Distribution | Control_EAF | Control_Missing_Rate | Control_HWE_P | A1_FREQ | OBS_CT | N_case | N_ctrl |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr16:53882967:G:A | rs9934504 | FTO | narrow_mainland | suggestive | A | G | 1.671 | 1.372 | 2.036 | 1.358e-07 | 214/157/40 | 0.2883 | 0.01909 | 0.1854 | 1082/548/74 | 0.2042 | 0.04054 | 0.655 | 0.2125 | 2,115 | 419 | 1,776 |
| chr16:53882967:G:A | rs9934504 | FTO | intermediate_mainland | genome_wide | A | G | 1.679 | 1.386 | 2.033 | 3.836e-08 | 217/163/41 | 0.291 | 0.01865 | 0.2368 | 1259/651/87 | 0.2066 | 0.03944 | 0.7851 | 0.2133 | 2,418 | 429 | 2,079 |
| chr16:53882967:G:A | rs9934504 | FTO | full_mainland | suggestive | A | G | 1.611 | 1.338 | 1.941 | 1.832e-07 | 223/166/42 | 0.29 | 0.01822 | 0.1977 | 1595/845/124 | 0.2131 | 0.03681 | 0.3772 | 0.2165 | 2,995 | 439 | 2,662 |

## How to read it

1. Do not read agreement between cohorts as replication. They share samples by construction; narrow is a subset of intermediate, which is a subset of full.
2. In (a), read the vertical spread first: the three models separate far more than the three cohorts do, so model choice dominates calibration here.
3. In (b), read down each variant's three rows. An interval that narrows while the point estimate holds is a variant gaining sample size. An estimate that drifts toward 1 as the filter is relaxed is what a structure-driven signal looks like.
4. A peak that appears only in narrow deserves suspicion — the narrow filter retains 95% of cases but only 67% of controls, so it is the set most exposed to ancestry confounding.
5. For every lead of both tiers in every cohort, not just the genome-wide ones, read `_comparison/tables/lead_crosscohort.tsv`.

## What this figure does *not* establish

- It cannot adjudicate between the cohorts. Choosing one requires a model that absorbs fine-scale structure, which this fixed-effects scan is not.
- It is not a replication analysis and no cohort here is independent of another. The component has no external cohort available.
- Panel (b) covers the main model only. Any other model's peaks are in each cohort's `03.peaks/model_peaks_annotation.tsv` and on their own scan figures.
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
