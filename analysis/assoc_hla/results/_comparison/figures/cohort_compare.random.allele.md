# HLA allele scan across 3 nested sample sets

**Figure file:** `cohort_compare.random.allele.png`

## The question this figure answers

How does the evidence for each HLA allele marker move as the ancestry filter is relaxed from the narrowest sample set to the widest?

## Panels

**(a) Evidence on one marker axis**

-log10 P for every allele marker with a usable fit, one series per cohort on ONE shared marker axis. The axis is ordered by gene, then by the coordinate build_hla_markers.py mints within the gene; alternating bands separate the genes and the tick labels name them. HLA markers have no meaningful base-pair spacing, so this is an ordered categorical axis and never a genomic distance. The solid line is the significance threshold and the dashed line nominal P = 0.05; both are drawn on every -log10 P axis in this study. Diamonds mark markers clearing the threshold in that cohort — shape carries the call, colour carries the cohort, so the two encodings never compete. At most 8 markers are named above the data.

**(b) Effect between cohort pairs**

log odds ratio in one cohort against another, one sub-panel per pair, over the markers fitted in both. The grey line is equality; dotted lines mark no effect. Points significant in either member of the pair are red. The annotated r is the Pearson correlation over every shared marker. These correlations are NOT replication statistics: the smaller cohort is a subset of the larger, so a high r is what the sample overlap alone would produce.

**(c) The hits in every cohort**

A forest over (marker) x (cohort). Every marker significant in ANY cohort is drawn in ALL of them, because a nested design always has an estimate in the larger sets. Three states are distinguished, not two: filled diamond = significant in that cohort, open circle = fitted but not significant, and a row with no fit at all is annotated in words. A blank would be indistinguishable from a missing estimate. The group label names the marker and its gene; the row ticks name the cohort, so the marker id is written once rather than three times.

**(d) Shared and cohort-specific**

Per cohort, the number of markers clearing the threshold, split into those also significant somewhere else (solid) and those significant in this cohort and no other (hatched), with the cohort-specific markers named to the right. Because the cohorts are nested, "shared" here means the same samples were counted twice — it is a bookkeeping statement about the three analyses, not evidence.

## Interpretation

The cohorts are NESTED: each is a subset of the next, and they share every case and most controls. Nothing in this figure is replication, and the three-row groups of panel (c) are the shape a reader is most likely to misread as such. What the figure does report is the behaviour of one estimate as the sample grows: an interval that narrows while the point estimate holds is a marker gaining sample size, and an estimate that drifts toward OR = 1 as the filter is relaxed is what a structure-driven signal does. The HLA region makes this sharper than elsewhere in the study, because class I and II allele frequencies differ between Japanese subpopulations, so the ancestry filter that defines the three sets acts directly on the markers being tested. A marker significant only in the narrowest set is therefore the case that most needs explaining, and panel (d) is what makes those countable instead of something the reader has to find by eye.

## Values in this rendering

| quantity | value |
|---|---|
| marker class | allele |
| model | random |
| cohorts (nesting order) | narrow_mainland, intermediate_mainland, full_mainland |
| markers with a usable fit | 170 |
| fits excluded on ERRCODE | 0 |
| significance threshold | 2.941e-04 |
| threshold source | Bonferroni 0.05 / 170 |
| markers significant in at least one cohort | 1 |
| markers significant in every cohort | 1 |
| markers drawn in the forest | 1 |
| significant, narrow_mainland | 1 |
| significant, intermediate_mainland | 1 |
| significant, full_mainland | 1 |
| significant only in narrow_mainland | 0 |
| significant only in intermediate_mainland | 0 |
| significant only in full_mainland | 0 |
| log-OR correlation narrow_mainland vs intermediate_mainland | 0.9802 |
| log-OR correlation narrow_mainland vs full_mainland | 0.9422 |
| log-OR correlation intermediate_mainland vs full_mainland | 0.9667 |

## Full statistics

**Every significant marker, in every cohort**

| marker | gene | position | cohort | model | OR | L95 | U95 | log_or | SE | P | A1_FREQ | OBS_CT | significant | n_cohorts_significant | n_cohorts_fitted | shared_by_all_fitted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HLA_F*01 | F |  | narrow_mainland | random | 2.886 | 1.585 | 5.253 | 1.06 | 0.3056 | 2.625e-04 | 0.006378 | 2,195 | 1 | 3 | 3 | 1 |
| HLA_F*01 | F |  | intermediate_mainland | random | 3.043 | 1.678 | 5.52 | 1.113 | 0.3038 | 1.247e-04 | 0.005582 | 2,508 | 1 | 3 | 3 | 1 |
| HLA_F*01 | F |  | full_mainland | random | 2.995 | 1.651 | 5.434 | 1.097 | 0.3039 | 1.536e-04 | 0.004837 | 3,101 | 1 | 3 | 3 | 1 |

## How to read it

1. Read panel (a) down, not across: a marker is at the same x in all three series, so the vertical spread at one x is what the ancestry filter did to that marker.
2. Do not read agreement between cohorts as replication anywhere in this figure. They are nested and share every case.
3. In (c), read down a group. An interval that narrows while the point estimate holds is a marker gaining sample size; an estimate drifting toward 1 as the filter opens is what a structure-driven signal does.
4. In (d), a hatched segment on the narrowest cohort is the finding that most needs explaining, because that set is the one most exposed to ancestry confounding.
5. For every marker in every cohort, not only the significant ones, read the TSV beside this figure.

## What this figure does *not* establish

- It is not a replication analysis and cannot be one: the cohorts are nested and no independent HLA-typed cohort is available to this study.
- It cannot adjudicate between the cohorts. Which sample set is reported is a design decision, and this figure is one of its inputs, not its verdict.
- It covers the allele class under the random model only. The other marker class and the other model are separate analyses with their own figures and their own multiple-testing burdens.
- The threshold treats the markers of this class as independent tests. HLA markers are in strong LD, so it is conservative in the number of tests and says nothing about the correlation between them.
- Panel (a) names at most the strongest few markers. An unnamed point is not a point the figure is claiming does not exist.

## Symbols

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **ERRCODE** — plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, SEPARATION) did not fit cleanly and are excluded from lambda_GC and from the hit list rather than silently carried.

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^Kgamma_k,PC_k,i      (g_i = genotype under the stated model)
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
