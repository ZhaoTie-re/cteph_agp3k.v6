# CTEPH (all cases) vs controls: markers significant in at least one cohort

**Covers:** `hits.allele.png`, `hits.residue.png` in `04.hits/`

## The question this figure answers

For each marker significant at any of the three thresholds in at least one cohort, how large is the effect, and in how many of the three cohorts does it hold?

## Panels

**(hits.allele.png) one row per test, six estimates per row**

Odds ratio per copy carrying the marker, with its 95 % CI, in each cohort (colour, dark narrow → light full) under each model (circle fixed, square random, the primary model). Filled = significant at the LD-adjusted (GEC) threshold in that cell. Rows in gene order, then position. A two-residue position whose two columns are complements in every cohort (a mirror pair, `mirror_of` in the tables) is one test: it is drawn once, as the residue with OR > 1 in the full-mainland random cell. At most 24 tests are drawn — those significant at GEC in most cohorts, then smallest P; the caption gives the totals.

**(hits.residue.png) one row per test, six estimates per row**

Odds ratio per copy carrying the marker, with its 95 % CI, in each cohort (colour, dark narrow → light full) under each model (circle fixed, square random, the primary model). Filled = significant at the LD-adjusted (GEC) threshold in that cell. Rows in gene order, then position. A two-residue position whose two columns are complements in every cohort (a mirror pair, `mirror_of` in the tables) is one test: it is drawn once, as the residue with OR > 1 in the full-mainland random cell. At most 24 tests are drawn — those significant at GEC in most cohorts, then smallest P; the caption gives the totals.

## Interpretation

alleles: 1 markers significant in at least one cohort under some threshold, 0 at the LD-adjusted threshold in 3/3 cohorts; amino-acid residues: 3 markers significant in at least one cohort under some threshold, 1 at the LD-adjusted threshold in 3/3 cohorts. Consistency is read in the random model with the same direction of effect in every significant cohort; the cohorts are nested, so 3/3 is a sample-selection sensitivity analysis, not replication.

## Values in this rendering

| quantity | value |
|---|---|
| hits.allele: markers significant in >= 1 cohort | 1 |
| hits.allele: mirror pairs merged | 0 |
| hits.allele: distinct tests | 1 |
| hits.allele: markers drawn | 1 |
| hits.allele: LD-adjusted (GEC), random model: markers 3/3 | 0 |
| hits.allele: LD-adjusted (GEC), random model: markers 2/3 | 0 |
| hits.allele: LD-adjusted (GEC), random model: markers 1/3 | 0 |
| hits.allele: Exact (max-T), random model: markers 3/3 | 0 |
| hits.allele: Exact (max-T), random model: markers 2/3 | 0 |
| hits.allele: Exact (max-T), random model: markers 1/3 | 0 |
| hits.allele: Conventional Bonferroni, random model: markers 3/3 | 0 |
| hits.allele: Conventional Bonferroni, random model: markers 2/3 | 0 |
| hits.allele: Conventional Bonferroni, random model: markers 1/3 | 0 |
| hits.residue: markers significant in >= 1 cohort | 3 |
| hits.residue: mirror pairs merged | 0 |
| hits.residue: distinct tests | 3 |
| hits.residue: markers drawn | 3 |
| hits.residue: LD-adjusted (GEC), random model: markers 3/3 | 1 |
| hits.residue: LD-adjusted (GEC), random model: markers 2/3 | 0 |
| hits.residue: LD-adjusted (GEC), random model: markers 1/3 | 0 |
| hits.residue: Exact (max-T), random model: markers 3/3 | 1 |
| hits.residue: Exact (max-T), random model: markers 2/3 | 0 |
| hits.residue: Exact (max-T), random model: markers 1/3 | 0 |
| hits.residue: Conventional Bonferroni, random model: markers 3/3 | 0 |
| hits.residue: Conventional Bonferroni, random model: markers 2/3 | 0 |
| hits.residue: Conventional Bonferroni, random model: markers 1/3 | 1 |

## Full statistics

**Every marker significant in at least one cohort (also _comparison/hits.tsv)**

| class | marker | mirror_of | GEC (random) | max-T (random) | Bonferroni (random) | GEC (fixed) | direction | min_p | max_p | or_narrow_fixed | or_narrow_random | or_intermediate_fixed | or_intermediate_random | or_full_fixed | or_full_random |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | B*40:02 | — | 0/3 | 0/3 | 0/3 | 1/3 | — | 5.088e-04 | 0.01133 | 1.684 | 1.786 | 1.534 | 1.599 | 1.425 | 1.469 |
| residue | DQB1:57:V | — | 3/3 | 3/3 | 1/3 | 3/3 | protective | 4.570e-05 | 1.847e-04 | 0.591 | 0.6205 | 0.6262 | 0.6521 | 0.6273 | 0.6559 |
| residue | B:114:N | — | 0/3 | 0/3 | 0/3 | 1/3 | — | 3.607e-04 | 0.001069 | 1.346 | 1.329 | 1.325 | 1.309 | 1.352 | 1.335 |
| residue | B:114:D | — | 0/3 | 0/3 | 0/3 | 1/3 | — | 3.819e-04 | 0.001051 | 0.7399 | 0.7495 | 0.7539 | 0.7635 | 0.7403 | 0.7498 |

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
