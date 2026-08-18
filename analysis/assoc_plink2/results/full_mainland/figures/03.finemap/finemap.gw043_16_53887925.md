# Fine-mapping — full_mainland, peak gw043_16_53887925

**Figure file:** `finemap.gw043_16_53887925.png`

## The question this figure answers

Does the posterior concentrate on a few variants, or does LD spread it out?

## Panels

**(a) Association**

The peak window, coloured by in-sample r2 with the lead (purple diamond). Shown so the PIP panel can be read against the evidence that produced it.

**(b) Posterior inclusion probability**

Per-variant PIP from susie_rss on the summary statistics and the in-sample signed-r LD matrix, L = 10. Rings mark credible-set membership; the bracket spans a set's physical extent. A credible set is the smallest group of variants carrying 95% of the posterior mass for one signal.

**(c) Resolution**

Cumulative posterior mass against the PIP-ranked members of each set. A curve reaching 0.95 in a few steps means the signal is resolved to those variants; a slow curve means LD has spread the mass and the set cannot be narrowed at this sample size.

## Interpretation

Fine-mapping is conditional on the LD matrix being the one the statistics came from, which is why in-sample LD is used here. At this effective sample size the posterior is driven by a small number of strongly associated variants, so a wide credible set should be read as insufficient resolution rather than as evidence against a single causal variant.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| peak | gw043_16_53887925 |
| lead variant | chr16:53887925:T:C |
| variants fine-mapped | 1,125 |
| GWAS n used | 3,071 |
| L (max signals) | 10 |
| coverage | 0.95 |
| purity filter min |r| | 0.5 |
| converged | 1 |
| credible sets | 1 |
| largest-PIP set size | 6 |
| top PIP | 0.4251 |
| susieR version | 0.14.2 |

## How to read it

1. A credible set of one or two variants is a resolved signal. A set of twenty means the data cannot distinguish among them, not that twenty variants are causal.
2. Check that the set contains the lead. If it does not, the lead is a tag and the posterior prefers a neighbour.
3. PIP is conditional on the LD matrix being the one the statistics came from, which is why in-sample LD is used here rather than a reference panel.

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

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
