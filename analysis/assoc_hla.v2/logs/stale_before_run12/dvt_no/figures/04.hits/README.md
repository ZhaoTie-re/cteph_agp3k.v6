# CTEPH without DVT vs controls: markers significant in at least one cohort

**Covers:** `hits.allele.png`, `hits.residue.png` in `04.hits/`

## The question this figure answers

For each marker significant at any of the three thresholds in at least one cohort, how large is the effect, and in how many of the three cohorts does it hold?

## Panels

**(hits.allele.png) one row per test, six estimates per row**

Odds ratio per copy carrying the marker, with its 95 % CI, in each cohort (colour, dark narrow → light full) under each model (circle fixed, square random, the primary model). Filled = significant at the LD-adjusted (GEC) threshold in that cell. Rows in gene order, then position. A two-residue position whose two columns are complements in every cohort (a mirror pair, `mirror_of` in the tables) is one test: it is drawn once, as the residue with OR > 1 in the full-mainland random cell. At most 24 tests are drawn — those significant at GEC in most cohorts, then smallest P; the caption gives the totals.

**(hits.residue.png) one row per test, six estimates per row**

Odds ratio per copy carrying the marker, with its 95 % CI, in each cohort (colour, dark narrow → light full) under each model (circle fixed, square random, the primary model). Filled = significant at the LD-adjusted (GEC) threshold in that cell. Rows in gene order, then position. A two-residue position whose two columns are complements in every cohort (a mirror pair, `mirror_of` in the tables) is one test: it is drawn once, as the residue with OR > 1 in the full-mainland random cell. At most 24 tests are drawn — those significant at GEC in most cohorts, then smallest P; the caption gives the totals.

## Interpretation

alleles: 1 markers significant in at least one cohort under some threshold, 0 at the LD-adjusted threshold in 3/3 cohorts; amino-acid residues: 51 markers significant in at least one cohort under some threshold, 0 at the LD-adjusted threshold in 3/3 cohorts. Consistency is read in the random model with the same direction of effect in every significant cohort; the cohorts are nested, so 3/3 is a sample-selection sensitivity analysis, not replication.

## Values in this rendering

| quantity | value |
|---|---|
| hits.allele: markers significant in >= 1 cohort | 1 |
| hits.allele: mirror pairs merged | 0 |
| hits.allele: distinct tests | 1 |
| hits.allele: markers drawn | 1 |
| hits.allele: LD-adjusted (GEC), random model: markers 3/3 | 0 |
| hits.allele: LD-adjusted (GEC), random model: markers 2/3 | 0 |
| hits.allele: LD-adjusted (GEC), random model: markers 1/3 | 1 |
| hits.allele: Exact (max-T), random model: markers 3/3 | 0 |
| hits.allele: Exact (max-T), random model: markers 2/3 | 0 |
| hits.allele: Exact (max-T), random model: markers 1/3 | 1 |
| hits.allele: Conventional Bonferroni, random model: markers 3/3 | 0 |
| hits.allele: Conventional Bonferroni, random model: markers 2/3 | 0 |
| hits.allele: Conventional Bonferroni, random model: markers 1/3 | 1 |
| hits.residue: markers significant in >= 1 cohort | 51 |
| hits.residue: mirror pairs merged | 17 |
| hits.residue: distinct tests | 34 |
| hits.residue: markers drawn | 24 |
| hits.residue: LD-adjusted (GEC), random model: markers 3/3 | 0 |
| hits.residue: LD-adjusted (GEC), random model: markers 2/3 | 0 |
| hits.residue: LD-adjusted (GEC), random model: markers 1/3 | 51 |
| hits.residue: Exact (max-T), random model: markers 3/3 | 0 |
| hits.residue: Exact (max-T), random model: markers 2/3 | 0 |
| hits.residue: Exact (max-T), random model: markers 1/3 | 51 |
| hits.residue: Conventional Bonferroni, random model: markers 3/3 | 0 |
| hits.residue: Conventional Bonferroni, random model: markers 2/3 | 0 |
| hits.residue: Conventional Bonferroni, random model: markers 1/3 | 0 |

## Full statistics

**Every marker significant in at least one cohort (also _comparison/hits.tsv)**

| class | marker | mirror_of | GEC (random) | max-T (random) | Bonferroni (random) | GEC (fixed) | direction | min_p | max_p | or_narrow_fixed | or_narrow_random | or_intermediate_fixed | or_intermediate_random | or_full_fixed | or_full_random |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | DPB1*05:01 | — | 1/3 | 1/3 | 1/3 | 1/3 | risk | 2.815e-04 | 0.001415 | 1.61 | 1.629 | 1.538 | 1.554 | 1.496 | 1.507 |
| residue | DQB1:87:Y | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.190e-04 | 3.290e-04 | 0.4768 | 0.5427 | 0.5096 | 0.5692 | 0.5184 | 0.5774 |
| residue | DQA1:52:R | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.519e-04 | 4.514e-04 | 1.679 | 1.658 | 1.617 | 1.602 | 1.585 | 1.567 |
| residue | DQA1:11:C | DQA1:11:Y | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:11:Y | DQA1:11:C | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:18:F | DQA1:18:S | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:18:S | DQA1:18:F | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:45:A | DQA1:45:V | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:45:V | DQA1:45:A | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:47:R | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:48:L | DQA1:48:W | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:48:W | DQA1:48:L | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:50:E | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:52:S | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:53:K | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:55:G | DQA1:55:R | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:55:R | DQA1:55:G | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:61:F | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:64:R | DQA1:64:T | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:64:T | DQA1:64:R | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:66:I | DQA1:66:M | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:66:M | DQA1:66:I | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:69:A | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:76:M | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:80:S | DQA1:80:Y | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.253e-04 | 1.678 | 1.654 | 1.619 | 1.6 | 1.591 | 1.57 |
| residue | DQA1:80:Y | DQA1:80:S | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.253e-04 | 0.5959 | 0.6046 | 0.6178 | 0.6251 | 0.6285 | 0.637 |
| residue | DQA1:218:Q | DQA1:218:R | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.565e-04 | 4.168e-04 | 0.5959 | 0.6046 | 0.6174 | 0.6245 | 0.6281 | 0.6365 |
| residue | DQA1:218:R | DQA1:218:Q | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.565e-04 | 4.168e-04 | 1.678 | 1.654 | 1.62 | 1.601 | 1.592 | 1.571 |
| residue | DQA1:61:G | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.631e-04 | 4.360e-04 | 0.5965 | 0.6053 | 0.6183 | 0.6256 | 0.6289 | 0.6374 |
| residue | DQB1:53:L | DQB1:53:Q | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.754e-04 | 4.760e-04 | 1.671 | 1.645 | 1.613 | 1.593 | 1.583 | 1.562 |
| residue | DQB1:53:Q | DQB1:53:L | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.754e-04 | 4.760e-04 | 0.5986 | 0.6081 | 0.6201 | 0.6277 | 0.6317 | 0.6404 |
| residue | DQB1:85:L | DQB1:85:V | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.797e-04 | 4.807e-04 | 1.67 | 1.644 | 1.612 | 1.593 | 1.583 | 1.561 |
| residue | DQB1:85:V | DQB1:85:L | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.797e-04 | 4.807e-04 | 0.5989 | 0.6083 | 0.6203 | 0.6279 | 0.6318 | 0.6405 |
| residue | DQB1:86:E | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.797e-04 | 4.807e-04 | 1.67 | 1.644 | 1.612 | 1.593 | 1.583 | 1.561 |
| residue | DQB1:87:L | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.797e-04 | 4.807e-04 | 1.67 | 1.644 | 1.612 | 1.593 | 1.583 | 1.561 |
| residue | DQB1:89:G | DQB1:89:T | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.797e-04 | 4.807e-04 | 0.5989 | 0.6083 | 0.6203 | 0.6279 | 0.6318 | 0.6405 |
| residue | DQB1:89:T | DQB1:89:G | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.797e-04 | 4.807e-04 | 1.67 | 1.644 | 1.612 | 1.593 | 1.583 | 1.561 |
| residue | DQB1:90:I | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.797e-04 | 4.876e-04 | 0.5989 | 0.6083 | 0.6207 | 0.6283 | 0.632 | 0.6407 |
| residue | DQB1:90:T | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.797e-04 | 4.807e-04 | 1.67 | 1.644 | 1.612 | 1.593 | 1.583 | 1.561 |
| residue | DQB1:84:E | DQB1:84:Q | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.860e-04 | 4.894e-04 | 0.5994 | 0.6089 | 0.6207 | 0.6283 | 0.6321 | 0.6408 |
| residue | DQB1:84:Q | DQB1:84:E | 1/3 | 1/3 | 0/3 | 1/3 | risk | 1.860e-04 | 4.894e-04 | 1.668 | 1.642 | 1.611 | 1.592 | 1.582 | 1.56 |
| residue | DQA1:175:Q | — | 1/3 | 1/3 | 0/3 | 1/3 | protective | 1.915e-04 | 5.034e-04 | 0.6011 | 0.6094 | 0.6246 | 0.6314 | 0.633 | 0.6412 |
| residue | DQB1:220:H | DQB1:220:R | 1/3 | 1/3 | 0/3 | 1/3 | risk | 2.151e-04 | 5.433e-04 | 1.656 | 1.639 | 1.601 | 1.59 | 1.58 | 1.559 |
| residue | DQB1:220:R | DQB1:220:H | 1/3 | 1/3 | 0/3 | 1/3 | protective | 2.151e-04 | 5.433e-04 | 0.6037 | 0.6101 | 0.6247 | 0.6289 | 0.6328 | 0.6413 |
| residue | DQB1:221:H | DQB1:221:Q | 1/3 | 1/3 | 0/3 | 1/3 | risk | 2.151e-04 | 5.433e-04 | 1.656 | 1.639 | 1.601 | 1.59 | 1.58 | 1.559 |
| residue | DQB1:221:Q | DQB1:221:H | 1/3 | 1/3 | 0/3 | 1/3 | protective | 2.151e-04 | 5.433e-04 | 0.6037 | 0.6101 | 0.6247 | 0.6289 | 0.6328 | 0.6413 |
| residue | DQB1:125:A | — | 1/3 | 1/3 | 0/3 | 1/3 | risk | 2.174e-04 | 5.728e-04 | 1.658 | 1.639 | 1.601 | 1.587 | 1.58 | 1.556 |
| residue | DQB1:-27:A | DQB1:-27:S | 1/3 | 1/3 | 0/3 | 0/3 | risk | 2.308e-04 | 6.192e-04 | 2.601 | 2.048 | 2.445 | 1.966 | 2.293 | 1.885 |
| residue | DQB1:-27:S | DQB1:-27:A | 1/3 | 1/3 | 0/3 | 0/3 | protective | 2.308e-04 | 6.192e-04 | 0.3845 | 0.4883 | 0.4091 | 0.5085 | 0.436 | 0.5304 |
| residue | DQB1:-9:I | DQB1:-9:M | 1/3 | 1/3 | 0/3 | 0/3 | protective | 2.308e-04 | 6.192e-04 | 0.3845 | 0.4883 | 0.4091 | 0.5085 | 0.436 | 0.5304 |
| residue | DQB1:-9:M | DQB1:-9:I | 1/3 | 1/3 | 0/3 | 0/3 | risk | 2.308e-04 | 6.192e-04 | 2.601 | 2.048 | 2.445 | 1.966 | 2.293 | 1.885 |
| residue | DQB1:-5:S | — | 1/3 | 1/3 | 0/3 | 0/3 | protective | 2.308e-04 | 6.192e-04 | 0.3845 | 0.4883 | 0.4091 | 0.5085 | 0.436 | 0.5304 |

## How to read it

1. Read the consistency columns first: GEC (primary), then max-T and Bonferroni.
2. In the figure, check that the six estimates lie on the same side of 1.
3. In the table, a row with `mirror_of` set and its partner are one test (OR and 1/OR, same P).

## What this figure does *not* establish

- Independence of the markers drawn: many are on one haplotype; see the conditional analyses (results/_conditional/signals_all.tsv).
- Exact effect sizes in the random model: they are SAIGE's one-step score estimates (Tstat / var), an approximation (SAIGE's Firth-corrected estimate was not requested, as in the genome-wide scan); read them for direction and rough size.
- The three cohorts are nested and share almost all cases, so agreement across cohorts is a check of the control selection, not a replication.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; an association that reflects a typing difference between platforms would look the same as a true one here.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
