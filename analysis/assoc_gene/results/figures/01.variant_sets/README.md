# Variant sets — what the map hands to the tests

**Figures:** `01.variant_sets/` — 3 figure(s): `variant_sets.full_mainland`, `variant_sets.intermediate_mainland`, `variant_sets.narrow_mainland`. One document covers the family; the numbers of each member are one row of the table below.

## The question this figure answers

Of the variants that pass the cohort's minor-allele-count floor, which enter a gene set, under which consequence, and how large are the sets each stratum produces?

## Panels

**(a) Call set by snpEff impact**

Every minAC-passing variant of the cohort, once, at its most severe annotation. MODIFIER (intergenic, deep intronic) is the vast majority and never enters a gene set; the three classes that can are HIGH, MODERATE and LOW.

**(b) Consequences by impact class**

The most frequent snpEff consequence terms within HIGH, MODERATE and LOW, each term counted once per variant at its most severe annotation; the remaining terms of each class are folded into an "Other" row whose count of distinct terms is given. Compound terms are shown by their first component with the number of additional components in brackets; the full strings are in annotation_summary.tsv.

**(c) Gene sets and variants per stratum**

For each impact stratum, the number of gene sets the map assigns (any size) beside the number that are tested (at least MinNumVar mapped variants), and the number of distinct tested variants. The tested count is the Bonferroni denominator of that cohort x stratum.

**(d) Set-size distribution**

The fraction of tested sets with at least x variants, per stratum, on a log axis. Most tested sets are small; the dashed guide is MinNumVar, below which a set is not tested and not counted.

## Interpretation

The figure shows the experiment that was designed, before any association was run: which variants qualify, how they distribute over consequence classes, and how many genes each stratum can test. Because MinNumVar is applied before the denominator is fixed, a set absent here is absent from the family the thresholds correct over. The three cohorts differ almost entirely in controls, so their maps differ mainly by which rare variants clear the minor-allele-count floor.

## Values in this rendering

| quantity | value |
|---|---|
| members | 3 |
| alpha | 0.05 |
| MinNumVar | 2 |

## Full statistics

**Per-cohort values**

| figure | cohort | variants in the call set | chromosomes | annotated | unannotated | MODIFIER-only | HIGH | MODERATE | LOW | MODIFIER | eligible (HIGH+MODERATE+LOW) | eligible percent | distinct consequence terms | multi-gene variants in the map | MinNumVar | legend rows in (a) | LOW+MODERATE+HIGH: sets produced | LOW+MODERATE+HIGH: sets tested | LOW+MODERATE+HIGH: variants tested | LOW+MODERATE+HIGH: Bonferroni | LOW+MODERATE+HIGH: median variants per set | MODERATE+HIGH: sets produced | MODERATE+HIGH: sets tested | MODERATE+HIGH: variants tested | MODERATE+HIGH: Bonferroni | MODERATE+HIGH: median variants per set |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| variant_sets.full_mainland | full_mainland | 8,832,209 | 22 | 8,832,209 | 0 | 8,657,032 | 4,283 | 74,808 | 99,916 | 8,653,202 | 179,007 | 2.027 | 48 | 2,067 | 2 | 6 | 19,398 | 15,691 | 171,641 | 3.18654005481e-06 | 6 | 15,815 | 12,378 | 75,787 | 4.03942478591e-06 | 4 |
| variant_sets.intermediate_mainland | intermediate_mainland | 7,620,757 | 22 | 7,620,757 | 0 | 7,470,269 | 3,604 | 63,894 | 86,335 | 7,466,924 | 153,833 | 2.019 | 48 | 1,762 | 2 | 6 | 18,857 | 15,009 | 146,820 | 3.33133453261e-06 | 6 | 15,207 | 11,436 | 63,867 | 4.37215809724e-06 | 4 |
| variant_sets.narrow_mainland | narrow_mainland | 6,857,863 | 22 | 6,857,863 | 0 | 6,722,024 | 3,260 | 57,512 | 78,080 | 6,719,011 | 138,852 | 2.025 | 48 | 1,612 | 2 | 6 | 18,468 | 14,557 | 132,096 | 3.43477364842e-06 | 5 | 14,801 | 10,832 | 56,931 | 4.61595273264e-06 | 4 |

## How to read it

1. Read (c) first: the tested-set count per stratum is the denominator every p-value in that cohort x stratum is corrected over.
2. Then (a) and (b): which consequences actually make up the sets — a stratum dominated by one term is a different experiment from a mixed one.
3. (d) says how small the sets are; a gene with two variants is tested, but its statistic is a count over a handful of alleles (METHODS §4.4).

## What this figure does *not* establish

- It shows annotation, not association: nothing here says which gene carries signal.
- Counts are per cohort and depend on that cohort's minor-allele-count floor; they are not comparable across cohorts as if the callset were fixed.

## Model

```
snpEff canonical map, per-gene most-severe consequence; sets need >= MinNumVar variants
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
