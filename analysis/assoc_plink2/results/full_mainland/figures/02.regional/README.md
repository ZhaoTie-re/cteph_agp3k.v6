# Regional association — full_mainland

**Figures:** `02.regional/<tier>/regional.<peak>.png` — 12 locus figure(s): 2 genome-wide, 10 suggestive.

One document covers the whole family: the panels, the reading and the limits are properties of the figure, identical for every locus. Only the numbers differ, and they are tabulated below, one row per locus.

## The question this figure answers

Is the LD structure that produces each peak a property of the East Asian population, or an artefact of this particular sample?

## Panels

**(a) In-sample cohort LD**

The association statistics of the peak window, each variant coloured by r2 with the lead computed on the *same* samples the statistics came from. This is the LD matrix the SuSiE fine-mapping uses, so this panel shows exactly what fine-mapping saw.

**(b) ToMMo 54KJPN**

Identical statistics, recoloured by r2 from the ToMMo 54KJPN co-occurrence tables. That resource publishes only pairs with r2 >= 0.2, so a variant present in the panel but returning no pair with the lead is *bounded* below 0.2 — it belongs in the lowest colour bin, not in grey. Grey is reserved for variants absent from the panel entirely.

**(c) 1000 Genomes EAS (n = 504)**

Identical statistics again, coloured by r2 in 1000 Genomes EAS (n = 504) — an out-of-sample population reference with no relationship to this cohort.

**(d) Gene models**

Filled boxes are exons, the connecting line spans the introns and chevrons give the transcribed strand. One representative transcript per gene (the longest protein-coding one), because drawing every transcript of a multi-transcript gene would need many rows to say the same thing. Only genes with an official symbol appear; clone-accession models name a sequencing clone rather than a gene and are suppressed.

## Interpretation

Agreement between the three panels means the LD structure is a property of the population rather than of this sample; divergence localised to the in-sample panel would indicate that the credible set is being shaped by sampling noise at this N. Each panel prints its own coverage in the table below, so a sparse source cannot be mistaken for a low-LD region.

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

| peak | tier | gene | lead variant | rsID | window | variants in window | genes drawn | cohort: r2 measured | cohort: bounded r2 < 0.2 | cohort: not in panel | tommo: r2 measured | tommo: bounded r2 < 0.2 | tommo: not in panel | 1000g_eas: r2 measured | 1000g_eas: bounded r2 < 0.2 | 1000g_eas: not in panel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gw003_2_11891017 | genome_wide | MIR3681HG | chr2:11891017:G:A | rs66898022 | chr2:11,624,908-12,305,379 | 1,581 | 8 | 1,581 | 0 | 0 | 55 | 1,526 | 0 | 1,547 | 0 | 34 |
| gw043_16_53887925 | genome_wide | FTO | chr16:53887925:T:C | rs16952623 | chr16:53,633,096-54,140,373 | 1,125 | 2 | 1,125 | 0 | 0 | 59 | 1,066 | 0 | 1,109 | 0 | 16 |
| sg006_2_232445059 | suggestive | ALPI | chr2:232445059:C:T | rs790041 | chr2:232,156,434-232,695,009 | 982 | 13 | 982 | 0 | 0 | 17 | 965 | 0 | 966 | 0 | 16 |
| sg007_2_233068979 | suggestive | INPP5D | chr2:233068979:G:A | rs72982244 | chr2:232,816,555-233,318,834 | 1,284 | 11 | 1,284 | 0 | 0 | 32 | 1,252 | 0 | 1,272 | 0 | 12 |
| sg009_3_94020684 | suggestive | ARL13B | chr3:94020684:T:C | rs143483852 | chr3:93,773,800-94,270,655 | 400 | 8 | 400 | 0 | 0 | 3 | 397 | 0 | 390 | 0 | 10 |
| sg010_3_100797923 | suggestive | ABI3BP | chr3:100797923:T:C | rs9881972 | chr3:100,534,960-101,080,800 | 986 | 6 | 986 | 0 | 0 | 41 | 945 | 0 | 971 | 0 | 15 |
| sg012_3_154069965 | suggestive | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | chr3:153,820,389-154,318,917 | 886 | 4 | 886 | 0 | 0 | 5 | 881 | 0 | 870 | 0 | 16 |
| sg014_4_171998307 | suggestive | GALNTL6 | chr4:171998307:T:C | rs75214171 | chr4:171,749,341-172,250,104 | 890 | 1 | 890 | 0 | 0 | 2 | 888 | 0 | 877 | 0 | 13 |
| sg018_6_32632724 | suggestive | HLA-DQA1 | chr6:32632724:A:G | rs9272120 | chr6:32,378,114-32,883,298 | 4,166 | 20 | 4,166 | 0 | 0 | 329 | 3,837 | 0 | 4,098 | 0 | 68 |
| sg031_10_93984133 | suggestive | SLC35G1 | chr10:93984133:C:G | rs142698390 | chr10:93,734,414-94,231,138 | 818 | 5 | 818 | 0 | 0 | 6 | 812 | 0 | 809 | 0 | 9 |
| sg045_17_13528059 | suggestive | HS3ST3A1 | chr17:13528059:G:A | rs34804183 | chr17:13,252,557-13,795,921 | 1,597 | 3 | 1,597 | 0 | 0 | 27 | 1,570 | 0 | 1,580 | 0 | 17 |
| sg052_20_60915853 | suggestive | CDH4 | chr20:60915853:C:T | rs79313991 | chr20:60,667,637-61,164,570 | 1,421 | 0 | 1,421 | 0 | 0 | 1 | 1,420 | 0 | 1,396 | 0 | 25 |

## How to read it

1. Compare the three colour patterns. If they agree, the LD block is a population property and the credible set can be trusted to reflect real correlation structure.
2. If the in-sample panel alone shows tight LD, the correlation is being driven by sampling noise at this N and the credible set is correspondingly fragile.
3. Read the lead against the gene track: whether it sits in an exon, an intron or between genes constrains which mechanisms are plausible.
4. Both tier lines are drawn in every panel, so a suggestive locus can be read against the genome-wide line it did not reach.

## What this figure does *not* establish

- It does not identify a causal variant. LD colour is correlation with the lead, not evidence of function — that is what the fine-mapping figures address.
- The two external panels differ from each other and from the study samples in ancestry breadth and in genomic coverage. Disagreement between (b) and (c) may reflect that difference rather than an error in either.

## Symbols

- **r^2** — linkage disequilibrium with the lead variant, binned on one scale for all sources. A co-occurrence panel typically publishes only pairs with r^2>=0.2, so a variant present in that panel but absent from the query is bounded below 0.2 and sits in the lowest bin; only variants absent from the panel entirely are unknown (grey).

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

---

Methods and rationale: [`METHODS.md`](../../../../docs/METHODS.md)
