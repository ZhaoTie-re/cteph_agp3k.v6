# Fine-mapping — full_mainland

**Figures:** `03.finemap/<tier>/finemap.<peak>.png` — 12 locus figure(s): 2 genome-wide, 10 suggestive.

One document covers the whole family: the panels, the reading and the limits are properties of the figure, identical for every locus. Only the numbers differ, and they are tabulated below, one row per locus.

## The question this figure answers

Does the posterior concentrate on a few variants, or does LD spread it out?

## Panels

**(a) Association**

The peak window, coloured by in-sample r2 with the lead (purple diamond). Shown so the PIP panel can be read against the evidence that produced it. Both tier lines are drawn.

**(b) Posterior inclusion probability**

Per-variant PIP from susie_rss on the summary statistics and the in-sample signed-r LD matrix. Rings mark credible-set membership; the bracket spans a set's physical extent. A credible set is the smallest group of variants carrying 95% of the posterior mass for one signal.

**(c) Resolution**

Cumulative posterior mass against the PIP-ranked members of each set. A curve reaching 0.95 in a few steps means the signal is resolved to those variants; a slow curve means LD has spread the mass and the set cannot be narrowed at this sample size.

## Interpretation

Fine-mapping is conditional on the LD matrix being the one the statistics came from, which is why in-sample LD is used here. At this effective sample size the posterior is driven by a small number of strongly associated variants, so a wide credible set should be read as insufficient resolution rather than as evidence against a single causal variant.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| loci | 12 |
| genome-wide loci | 2 |
| suggestive loci | 10 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-5 |

## Full statistics

**Per-locus values**

| peak | tier | gene | lead variant | rsID | variants fine-mapped | GWAS n used | L (max signals) | coverage | purity filter min |r| | converged | credible sets | largest-PIP set size | top PIP | susieR version |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gw003_2_11891017 | genome_wide | MIR3681HG | chr2:11891017:G:A | rs66898022 | 1,581 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 21 | 0.394 | 0.14.2 |
| gw043_16_53887925 | genome_wide | FTO | chr16:53887925:T:C | rs16952623 | 1,125 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 6 | 0.4251 | 0.14.2 |
| sg006_2_232445059 | suggestive | ALPI | chr2:232445059:C:T | rs790041 | 982 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 8 | 0.3456 | 0.14.2 |
| sg007_2_233068979 | suggestive | INPP5D | chr2:233068979:G:A | rs72982244 | 1,284 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 17 | 0.1759 | 0.14.2 |
| sg009_3_94020684 | suggestive | ARL13B | chr3:94020684:T:C | rs143483852 | 400 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.9994 | 0.14.2 |
| sg010_3_100797923 | suggestive | ABI3BP | chr3:100797923:T:C | rs9881972 | 986 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 8 | 0.2722 | 0.14.2 |
| sg012_3_154069965 | suggestive | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | 886 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.9823 | 0.14.2 |
| sg014_4_171998307 | suggestive | GALNTL6 | chr4:171998307:T:C | rs75214171 | 890 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 2 | 0.621 | 0.14.2 |
| sg018_6_32632724 | suggestive | HLA-DQA1 | chr6:32632724:A:G | rs9272120 | 4,166 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 30 | 0.07871 | 0.14.2 |
| sg031_10_93984133 | suggestive | SLC35G1 | chr10:93984133:C:G | rs142698390 | 818 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.9589 | 0.14.2 |
| sg045_17_13528059 | suggestive | HS3ST3A1 | chr17:13528059:G:A | rs34804183 | 1,597 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 9 | 0.4094 | 0.14.2 |
| sg052_20_60915853 | suggestive | CDH4 | chr20:60915853:C:T | rs79313991 | 1,421 | 3,071 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.9722 | 0.14.2 |

## How to read it

1. A credible set of one or two variants is a resolved signal. A set of twenty means the data cannot distinguish among them, not that twenty variants are causal.
2. Check that the set contains the lead. If it does not, the lead is a tag and the posterior prefers a neighbour.
3. PIP is conditional on the LD matrix being the one the statistics came from, which is why in-sample LD is used here rather than a reference panel.
4. Read `credible sets` = 0 in the table as "SuSiE found no set meeting the purity and coverage requirements", not as "no signal".

## What this figure does *not* establish

- SuSiE assumes the causal variant is present in the data. A causal variant not genotyped or filtered out cannot appear, and its posterior mass will be distributed over its tags.
- It cannot rank the biological plausibility of set members — only their statistical compatibility with the observed association pattern.

## Symbols

- **PIP** — SuSiE posterior inclusion probability — the probability a variant is causal given the locus summary statistics and an in-sample LD matrix. Credible sets are the smallest variant groups covering 95% posterior mass.

- **r^2** — linkage disequilibrium with the lead variant, binned on one scale for all sources. A co-occurrence panel typically publishes only pairs with r^2>=0.2, so a variant present in that panel but absent from the query is bounded below 0.2 and sits in the lowest bin; only variants absent from the panel entirely are unknown (grey).

## Model

```
susie_rss(hatbeta, se(hatbeta), R, n),   R = in-sample LD (signed r), L = 10, credible-set coverage 0.95
```

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
