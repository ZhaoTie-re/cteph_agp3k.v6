# Regional association — intermediate_mainland, peak gw031_16_53887925

**Figure file:** `regional.gw031_16_53887925.png`

## The question this figure answers

Is the LD structure that produces this peak a property of the East Asian population, or an artefact of this particular sample?

## Panels

**(a) In-sample cohort LD**

The association statistics of the peak window, each variant coloured by r2 with the lead computed on the *same* samples the statistics came from. This is the LD matrix the SuSiE fine-mapping uses, so this panel shows exactly what fine-mapping saw.

**(b) ToMMo 54KJPN**

Identical statistics, recoloured by r2 from the ToMMo 54KJPN co-occurrence tables. That resource publishes only pairs with r2 >= 0.2, so a variant present in the panel but returning no pair with the lead is *bounded* below 0.2 — it belongs in the lowest colour bin, not in grey. Grey is reserved for variants absent from the panel entirely. Treating the bound as missing data is what made this source look uninformative in the earlier implementation.

**(c) 1000 Genomes EAS (n = 504)**

Identical statistics again, coloured by r2 in 1000 Genomes EAS (n = 504) — an out-of-sample population reference with no relationship to this cohort.

**(d) Gene models**

Ensembl 86 / GRCh38. Filled boxes are exons, the connecting line spans the introns and chevrons give the transcribed strand. One representative transcript per gene (the longest protein-coding one), because drawing every transcript of a multi-transcript gene would need 22 rows to say the same thing. Only genes with an official symbol appear; clone-accession models such as AC007347.1 or RP11-357N13.3 name a sequencing clone rather than a gene and are suppressed.

## Interpretation

Agreement between the three panels means the LD structure is a property of the population rather than of this sample; divergence localised to the in-sample panel would indicate that the credible set is being shaped by sampling noise at this N. Each panel prints its own coverage, so a sparse source cannot be mistaken for a low-LD region.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| peak | gw031_16_53887925 |
| lead variant | chr16:53887925:T:C |
| window | chr16:53,633,096-54,197,338 |
| variants in window | 1,270 |
| genes drawn | 2 |
| cohort: r2 measured | 1,270 |
| cohort: bounded r2 < 0.2 | 0 |
| cohort: not in panel | 0 |
| tommo: r2 measured | 59 |
| tommo: bounded r2 < 0.2 | 1,211 |
| tommo: not in panel | 0 |
| 1000g_eas: r2 measured | 1,251 |
| 1000g_eas: bounded r2 < 0.2 | 0 |
| 1000g_eas: not in panel | 19 |

## How to read it

1. Compare the three colour patterns. If they agree, the LD block is a population property and the credible set can be trusted to reflect real correlation structure.
2. If the in-sample panel alone shows tight LD, the correlation is being driven by sampling noise at this N and the credible set is correspondingly fragile.
3. Read the lead against the gene track: whether it sits in an exon, an intron or between genes constrains which mechanisms are plausible.

## What this figure does *not* establish

- It does not identify a causal variant. LD colour is correlation with the lead, not evidence of function — that is what the fine-mapping figure addresses.
- The two external panels differ from each other and from the study samples in ancestry breadth and in genomic coverage. Disagreement between (b) and (c) may reflect that difference rather than an error in either.

## Symbols

- **r^2** — linkage disequilibrium with the lead variant, binned on one scale for all sources. A co-occurrence panel typically publishes only pairs with r^2>=0.2, so a variant present in that panel but absent from the query is bounded below 0.2 and sits in the lowest bin; only variants absent from the panel entirely are unknown (grey).

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
