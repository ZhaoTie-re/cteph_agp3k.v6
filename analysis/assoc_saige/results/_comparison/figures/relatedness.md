# Relatedness among the analysed samples

**Figure file:** `relatedness.png`

## The question this figure answers

What relatedness does the GRM encode among the samples this component analyses, and what does retaining it buy over a design that prunes it away?

## Panels

**(a) GRM off-diagonal, whole distribution**

The variance-standardized relationship between every pair of analysed samples, binned, on a log count axis so the tail is visible against a bulk that is orders of magnitude larger. This is the matrix the null model fitted, computed on the same LD-pruned markers. Dotted vertical lines are the expected values per degree — 0.0884, 0.177, 0.354, 0.708 — which are twice the corresponding KING kinship thresholds, because a GRM relationship is 2 x kinship. Negative values are expected: the estimator is not bounded below at zero and returns small negative values for pairs less related than the sample average.

**(b) GRM diagonal**

Each sample's relationship with itself, which is 1 + inbreeding. It is reported because SAIGE replaces it with exactly 1.0 when fitting the null (--isDiagofKinSetAsOne=True), so this is the only record of what the data actually said. A diagonal far from 1 indicates inbreeding or a genotyping problem.

**(c) Related pairs by degree**

Counts of pairs in each degree class, per cohort. A pair is assigned the strongest class whose threshold it clears.

**(d) What retaining relatedness is worth**

Filled bars: samples with at least one relative in the analysis. Open bars: samples present in this component's sample set but absent from the fixed-effects one, which prunes relateds. The open bar is the direct cost the fixed-effects design pays and this one does not.

## Interpretation

The mixed model exists so that related samples can be KEPT: the genetic relationship matrix accounts for the correlation between them instead of the design discarding one of every related pair. This figure measures that. Panel (a) shows the whole pairwise distribution, whose bulk sits near zero — the overwhelming majority of pairs are unrelated, and the informative content is the tail beyond the third-degree boundary. Panel (b) resolves that tail into degree classes. Panel (c) is the design claim as a number: how many samples carry a relative, and how many are present here but absent from the relatedness-pruned set the fixed-effects component analyses. Relatedness is accounted for, not removed, so these samples contribute to every estimate. Note that kinship is computed on the LD-pruned GRM markers, which is the matrix the null model actually fitted, not a separately ascertained marker set.

## Values in this rendering

| quantity | value |
|---|---|
| narrow_mainland: n_samples | 2,195 |
| narrow_mainland: n_related_pairs | 55 |
| narrow_mainland: n_samples_with_relative | 94 |
| narrow_mainland: n_cases_with_relative | 0 |
| narrow_mainland: n_absent_from_fixed_model | 27 |
| intermediate_mainland: n_samples | 2,508 |
| intermediate_mainland: n_related_pairs | 56 |
| intermediate_mainland: n_samples_with_relative | 96 |
| intermediate_mainland: n_cases_with_relative | 0 |
| intermediate_mainland: n_absent_from_fixed_model | 28 |
| full_mainland: n_samples | 3,101 |
| full_mainland: n_related_pairs | 70 |
| full_mainland: n_samples_with_relative | 123 |
| full_mainland: n_cases_with_relative | 0 |
| full_mainland: n_absent_from_fixed_model | 30 |

## Full statistics

**Every related pair**

| cohort | id1 | id2 | grm | degree |
|---|---|---|---|---|
| narrow_mainland | NAG16954 | NAG10123 | 0.5648 | first |
| narrow_mainland | NAG15145 | NAG11677 | 0.5339 | first |
| narrow_mainland | NAG19002 | NAG16040 | 0.5289 | first |
| narrow_mainland | AYUM0119 | AYUM0101 | 0.5258 | first |
| narrow_mainland | HVUM0478 | AYUM0643 | 0.5187 | first |
| narrow_mainland | NAG19095 | NAG11456 | 0.5006 | first |
| narrow_mainland | AYUM0047 | AYUM0008 | 0.4981 | first |
| narrow_mainland | NAG18985 | NAG16557 | 0.4955 | first |
| narrow_mainland | AYUM0655 | AYUM0285 | 0.4858 | first |
| narrow_mainland | AYUM0626 | AYUM0580 | 0.4477 | first |
| narrow_mainland | NAG13819 | NAG12423 | 0.2985 | second |
| narrow_mainland | NAG17626 | NAG12647 | 0.2853 | second |
| narrow_mainland | NAG16040 | NAG12313 | 0.2835 | second |
| narrow_mainland | NAG16581 | NAG15233 | 0.2756 | second |
| narrow_mainland | NAG16239 | NAG13402 | 0.2745 | second |
| narrow_mainland | NAG17496 | NAG12259 | 0.2727 | second |
| narrow_mainland | NAG19002 | NAG12313 | 0.2658 | second |
| narrow_mainland | NAG13381 | NAG10775 | 0.2549 | second |
| narrow_mainland | AYUM0119 | AYUM0118 | 0.2549 | second |
| narrow_mainland | NAG15050 | NAG13426 | 0.2462 | second |
| narrow_mainland | NAG19159 | NAG16629 | 0.2381 | second |
| narrow_mainland | NAG17523 | NAG15068 | 0.23 | second |
| narrow_mainland | HVUM0617 | HVUM0390 | 0.2275 | second |
| narrow_mainland | AYUM0118 | AYUM0101 | 0.2188 | second |
| narrow_mainland | NAG14708 | NAG11177 | 0.1757 | third |
| narrow_mainland | NAG12161 | NAG10536 | 0.167 | third |
| narrow_mainland | NAG17733 | NAG12161 | 0.1587 | third |
| narrow_mainland | NAG17815 | NAG14246 | 0.1523 | third |
| narrow_mainland | NAG15059 | NAG12765 | 0.1505 | third |
| narrow_mainland | NAG12313 | NAG10418 | 0.1491 | third |
| narrow_mainland | NAG11848 | NAG10323 | 0.1468 | third |
| narrow_mainland | NAG17733 | NAG10536 | 0.1434 | third |
| narrow_mainland | NAG17466 | NAG16231 | 0.139 | third |
| narrow_mainland | NAG19141 | NAG14686 | 0.1372 | third |
| narrow_mainland | NAG17298 | NAG16260 | 0.1349 | third |
| narrow_mainland | NAG16120 | NAG13997 | 0.1347 | third |
| narrow_mainland | NAG12496 | NAG12313 | 0.1321 | third |
| narrow_mainland | NAG15737 | NAG15589 | 0.1319 | third |
| narrow_mainland | NAG16260 | NAG13364 | 0.1234 | third |
| narrow_mainland | NAG17298 | NAG13364 | 0.1224 | third |
| … | 141 further rows in the TSV |  |  |  |

## How to read it

1. Read panel (a) on the log axis: the bulk near zero is every unrelated pair, and the question is only how far the right tail extends.
2. Check panel (b) against 1.0. A systematic shift there is a property of the marker set or the samples, not of any one pair, and it is invisible in the fitted model because SAIGE overwrites the diagonal.
3. A pair beyond the 1st-degree line is a parent-offspring or full-sibling pair. A duplicate or MZ twin would sit near 1.0 and would be worth investigating as a sample-handling question, not a biological one.
4. In panel (d), compare the two bars within a cohort, not across cohorts.
5. If cases carry relatives, the relatedness is not independent of phenotype and the GRM is doing more than nuisance correction — check n_cases_with_relative in the table.

## What this figure does *not* establish

- A GRM value conflates relatedness with SHARED ANCESTRY. These cohorts are ancestry-filtered subsets of one another, so the off-diagonal bulk shifts with the filter and the degree bands are indicative rather than diagnostic. A kinship estimator such as KING-robust is designed to be insensitive to structure and would classify degrees more reliably; what this figure gains in exchange is that it reports the values the model actually used.
- It says nothing about whether the GRM has absorbed the relatedness adequately. That question is answered by the calibration read-out, not here.
- The diagonal shown is the computed one. The fitted model does not use it.
- Degree classes are threshold calls on a noisy estimator. A pair near a boundary may belong to either class, and the counts should be read as approximate.

## Symbols

- **N_eff** — 4/(1/N_case+1/N_ctrl), the balanced-design equivalent sample size. Reported as a cohort descriptor — it is the size any later meta-analysis or mixed-model run quotes — and it collapses toward 4x the smaller arm as the design becomes unbalanced.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
