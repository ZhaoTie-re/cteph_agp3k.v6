# Regional association — narrow_mainland

**Figures:** `02.regional/<tier>/regional.<peak>.png` — 10 locus figure(s): 0 genome-wide, 10 suggestive.

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
| cohort | narrow_mainland |
| loci | 10 |
| genome-wide loci | 0 |
| suggestive loci | 10 |
| genome-wide threshold | 5e-8 |
| suggestive threshold | 1e-5 |

## Full statistics

**Per-locus values**

| peak | tier | gene | lead variant | rsID | window | variants in window | genes drawn | cohort: r2 measured | cohort: bounded r2 < 0.2 | cohort: not in panel | tommo: r2 measured | tommo: bounded r2 < 0.2 | tommo: not in panel | 1000g_eas: r2 measured | 1000g_eas: bounded r2 < 0.2 | 1000g_eas: not in panel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sg004_3_94020684 | suggestive | ARL13B | chr3:94020684:T:C | rs143483852 | chr3:93,773,800-94,270,655 | 396 | 8 | 396 | 0 | 0 | 3 | 393 | 0 | 388 | 0 | 8 |
| sg005_3_100783010 | suggestive | ABI3BP | chr3:100783010:T:C | rs7648469 | chr3:100,534,960-101,080,800 | 989 | 6 | 989 | 0 | 0 | 41 | 948 | 0 | 974 | 0 | 15 |
| sg007_3_154069965 | suggestive | ARHGEF26-AS1 | chr3:154069965:A:G | rs1021580972 | chr3:153,820,389-154,318,917 | 878 | 4 | 878 | 0 | 0 | 5 | 873 | 0 | 864 | 0 | 14 |
| sg008_4_171998307 | suggestive | GALNTL6 | chr4:171998307:T:C | rs75214171 | chr4:171,749,341-172,250,104 | 891 | 1 | 891 | 0 | 0 | 2 | 889 | 0 | 875 | 0 | 16 |
| sg015_8_138575524 | suggestive | COL22A1 | chr8:138575524:T:C | rs62530012 | chr8:138,246,096-138,824,057 | 1,311 | 2 | 1,311 | 0 | 0 | 46 | 1,265 | 0 | 1,293 | 0 | 18 |
| sg022_14_88019710 | suggestive | LINC01146 | chr14:88019710:G:A | rs4899932 | chr14:87,675,490-88,273,183 | 1,621 | 5 | 1,621 | 0 | 0 | 315 | 1,306 | 0 | 1,603 | 0 | 18 |
| sg023_16_53882967 | suggestive | FTO | chr16:53882967:G:A | rs9934504 | chr16:53,633,096-54,140,373 | 1,117 | 2 | 1,117 | 0 | 0 | 59 | 1,058 | 0 | 1,102 | 0 | 15 |
| sg026_19_10631611 | suggestive | SLC44A2 | chr19:10631611:C:T | rs1560711 | chr19:10,378,495-10,881,151 | 536 | 23 | 536 | 0 | 0 | 94 | 442 | 0 | 527 | 0 | 9 |
| sg028_20_7553182 | suggestive | RN7SL547P | chr20:7553182:C:T | rs117662178 | chr20:7,305,040-7,802,912 | 1,123 | 2 | 1,123 | 0 | 0 | 2 | 1,121 | 0 | 1,112 | 0 | 11 |
| sg029_20_60915853 | suggestive | CDH4 | chr20:60915853:C:T | rs79313991 | chr20:60,667,637-61,164,570 | 1,404 | 0 | 1,404 | 0 | 0 | 1 | 1,403 | 0 | 1,385 | 0 | 19 |

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
