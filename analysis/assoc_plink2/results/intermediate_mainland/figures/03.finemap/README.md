# Fine-mapping — intermediate_mainland

**Figures:** `03.finemap/<tier>/finemap.<peak>.png` — 11 locus figure(s): 1 genome-wide, 10 suggestive.

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
| cohort | intermediate_mainland |
| loci | 11 |
| genome-wide loci | 1 |
| suggestive loci | 10 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-5 |

## Full statistics

**Per-locus values**

| peak | tier | gene | lead variant | rsID | variants fine-mapped | GWAS n used | L (max signals) | coverage | purity filter min |r| | converged | credible sets | largest-PIP set size | top PIP | susieR version |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gw031_16_53887925 | genome_wide | FTO | chr16:53887925:T:C | rs16952623 | 1,270 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 6 | 0.377 | 0.14.2 |
| sg010_3_94020684 | suggestive | ARL13B | chr3:94020684:T:C | rs143483852 | 402 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.999 | 0.14.2 |
| sg011_3_100797923 | suggestive | ABI3BP | chr3:100797923:T:C | rs9881972 | 989 | 2,480 | 10 | 0.95 | 0.5 | 1 | 0 | — | — | 0.14.2 |
| sg012_3_154069965 | suggestive | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | 879 | 2,480 | 10 | 0.95 | 0.5 | 1 | 2 | 2 | 0.9435 | 0.14.2 |
| sg014_4_171998307 | suggestive | GALNTL6 | chr4:171998307:T:C | rs75214171 | 894 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 2 | 0.6864 | 0.14.2 |
| sg015_5_66174098 | suggestive | SREK1 | chr5:66174098:A:G | rs4700094 | 843 | 2,480 | 10 | 0.95 | 0.5 | 1 | 0 | — | — | 0.14.2 |
| sg017_6_32632724 | suggestive | HLA-DQA1 | chr6:32632724:A:G | rs9272120 | 4,161 | 2,480 | 10 | 0.95 | 0.5 | 1 | 0 | — | — | 0.14.2 |
| sg028_14_87924989 | suggestive | GALC | chr14:87924989:G:A | rs12881373 | 1,642 | 2,480 | 10 | 0.95 | 0.5 | 1 | 0 | — | — | 0.14.2 |
| sg033_17_13528059 | suggestive | HS3ST3A1 | chr17:13528059:G:A | rs34804183 | 1,541 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 9 | 0.438 | 0.14.2 |
| sg034_17_76739260 | suggestive | MFSD11 | chr17:76739260:T:C | rs9897202 | 1,194 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 24 | 0.1091 | 0.14.2 |
| sg036_20_60915853 | suggestive | CDH4 | chr20:60915853:C:T | rs79313991 | 1,411 | 2,480 | 10 | 0.95 | 0.5 | 1 | 1 | 1 | 0.9593 | 0.14.2 |

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
