# Per-gene detail — Tier 1 and Tier 2 genes

**Figures:** `<GENE>.png` — 5 gene figure(s), one per Tier 1 or Tier 2 gene (STBD1, OSBP, PIK3C2G, RPH3A, JAK2). One document covers the family; each gene is one row of the tables below.

## The question this figure answers

Which variants make up each robust gene's set, where are its carriers, and does the burden effect hold in every cell of the design?

## Panels

**(a) Variants and carrier frequency, one strip per cohort**

Every variant of the widest set the gene is mapped in, in full (the report cohort, which decides WHICH set the figure is about), at its genomic position in Mb. The three strips are the three cohorts, drawn on ONE shared y scale so the same variant is the same height in each. The y axis is the percentage of that cohort's own group carrying at least one minor allele of the variant -- cases up, controls down. It is a percentage and not a count because the groups are of very different size: on a count axis a variant carried by 0.9 % of 2,630 controls out-stems one carried by 2.1 % of 438 cases, and the panel reads as control-dominated when the case excess is the finding. Marker shape is the snpEff impact class (diamond HIGH, circle MODERATE, square LOW). A variant absent from a narrower cohort's map draws nothing in that strip, which is not the same as a drawn 0 %: each cohort sets its own minAC, so its map is its own, and reading down the strips shows exactly which variants the larger sample adds. The x axis is shared with (b), so a marker sits exactly above the exon it falls in.

**(b) Gene models**

The exon/intron structure of the tested gene and of any gene overlapping its span, from Ensembl 86 (GRCh38) — one representative transcript per gene, the one with the greatest total exonic length. Filled boxes are exons, the connecting line spans the introns and the chevrons give the transcribed strand. Genes that merely lie in the window without overlapping the tested gene are not drawn: the figure is about one gene, and its neighbours would cost rows without adding to it. Where Ensembl 86 has no model under the tested symbol the panel says so rather than drawing a neighbour in its place.

**(c) CMC burden odds ratio in every cell, with both P-values**

exp(β) and exp(β ± 1.96 SE) from rvtest --burden cmcWald per variant set and cohort, log axis, dashed line at OR = 1. Marker colour is the CMC state of that cell (full red Bonferroni + BH, light red BH only, blue tested and not called); a grey cross is a cell in which the gene is not in the map; "OR not estimable" is a Wald fit rvtest reset. To the right, the CMC score P and the SKAT-O P of the same cell, EACH COLOURED BY ITS OWN STATISTIC'S STATE IN THAT CELL: dark red = above that cell's Bonferroni threshold (and therefore BH-called), light red = BH-called only, blue = tested and called by neither, — = the gene is not in that cohort's map. So a row can carry a red SKAT-O P beside a blue CMC P: the two statistics disagree, and the marker shows what the burden estimate was.

## Interpretation

Read (a) for how the signal is distributed: one tall stem is a single-variant finding wearing a gene-set P; many short stems is a burden. Read (b) for whether the effect is stable as controls are added (down a block) and as the set widens (between blocks). The cohorts are nested, so stability across them is not replication. Where a P is coloured but the OR sits near 1, the SKAT-O kernel is seeing structure the collapse cannot (a subset of variants, or opposite directions). Where samples carry two or more sites of the gene the caption says so; CMC counts such a sample once while the cumulative MAC counts every site.

## Values in this rendering

| quantity | value |
|---|---|
| genes drawn | 5 |
| tiers drawn | 1, 2 |
| report cohort | full_mainland |
| N case / N control | 438 / 2,630 |

## Full statistics

**Per-gene values**

| gene_symbol | set_name | chrom | pos_min | pos_max | tiers | stratum_drawn | n_variants_drawn | n_variants_per_cohort | n_genes_in_model_panel | n_multi_site_carriers_case | n_multi_site_carriers_control | png |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| STBD1 | STBD1 | 4 | 76,306,990 | 76,309,559 | low_moderate_high:2;moderate_high:1 | low_moderate_high | 8 | narrow:5;intermediate:6;full:8 | 2 | 0 | 0 | STBD1.png |
| OSBP | OSBP | 11 | 59,576,952 | 59,594,182 | moderate_high:1 | moderate_high | 2 | narrow:2;intermediate:2;full:2 | 2 | 0 | 0 | OSBP.png |
| PIK3C2G | PIK3C2G | 12 | 18,282,395 | 18,647,868 | low_moderate_high:2;moderate_high:2 | low_moderate_high | 21 | narrow:15;intermediate:17;full:21 | 2 | 0 | 2 | PIK3C2G.png |
| RPH3A | RPH3A | 12 | 112,866,794 | 112,896,647 | low_moderate_high:2;moderate_high:2 | low_moderate_high | 11 | narrow:11;intermediate:11;full:11 | 1 | 0 | 0 | RPH3A.png |
| JAK2 | JAK2 | 9 | 5,044,432 | 5,080,294 | moderate_high:2 | moderate_high | 7 | narrow:6;intermediate:6;full:7 | 1 | 1 | 0 | JAK2.png |

**Cells drawn in (b) (from evidence.tsv)**

| gene_symbol | stratum | cohort | cmc_state | skato_state | or | or_l95 | or_u95 | cmc_p | skato_p |
|---|---|---|---|---|---|---|---|---|---|
| STBD1 | low_moderate_high | narrow_mainland | not_significant | called_bh_only | 6.946 | 2.729 | 17.68 | 4.102e-06 | 6.722e-06 |
| STBD1 | low_moderate_high | intermediate_mainland | not_significant | called_bh_only | 4.767 | 2.178 | 10.44 | 2.031e-05 | 1.330e-05 |
| STBD1 | low_moderate_high | full_mainland | not_significant | called_bonferroni | 4.059 | 1.984 | 8.304 | 3.661e-05 | 8.884e-07 |
| STBD1 | moderate_high | narrow_mainland | not_significant | called_bh_only | 9.978 | 2.931 | 33.96 | 1.316e-05 | 2.132e-05 |
| STBD1 | moderate_high | intermediate_mainland | called_bh_only | called_bh_only | 7.349 | 2.705 | 19.96 | 7.197e-06 | 1.209e-05 |
| STBD1 | moderate_high | full_mainland | called_bonferroni | called_bonferroni | 9.601 | 3.565 | 25.86 | 7.226e-08 | 1.411e-07 |
| OSBP | low_moderate_high | narrow_mainland | not_significant | not_significant | 4.054 | 1.692 | 9.715 | 7.827e-04 | 2.164e-04 |
| OSBP | low_moderate_high | intermediate_mainland | not_significant | not_significant | 4.245 | 1.887 | 9.548 | 1.634e-04 | 6.701e-05 |
| OSBP | low_moderate_high | full_mainland | not_significant | not_significant | 2.767 | 1.349 | 5.672 | 0.003936 | 0.003236 |
| OSBP | moderate_high | narrow_mainland | not_significant | called_bh_only | 16.23 | 3.374 | 78.1 | 4.707e-06 | 5.144e-06 |
| OSBP | moderate_high | intermediate_mainland | called_bonferroni | called_bonferroni | 19.39 | 4.047 | 92.87 | 3.541e-07 | 3.211e-07 |
| OSBP | moderate_high | full_mainland | not_significant | called_bh_only | 7.105 | 2.53 | 19.95 | 1.933e-05 | 2.083e-05 |
| PIK3C2G | low_moderate_high | narrow_mainland | not_significant | called_bonferroni | 1.968 | 1.078 | 3.591 | 0.02512 | 8.458e-07 |
| PIK3C2G | low_moderate_high | intermediate_mainland | not_significant | called_bonferroni | 2.377 | 1.37 | 4.124 | 0.001559 | 3.483e-07 |
| PIK3C2G | low_moderate_high | full_mainland | not_significant | called_bh_only | 2.252 | 1.413 | 3.591 | 4.754e-04 | 5.775e-06 |
| PIK3C2G | moderate_high | narrow_mainland | not_significant | called_bonferroni | 2.648 | 1.361 | 5.153 | 0.003054 | 4.132e-07 |
| PIK3C2G | moderate_high | intermediate_mainland | not_significant | called_bonferroni | 2.674 | 1.437 | 4.979 | 0.001306 | 5.641e-07 |
| PIK3C2G | moderate_high | full_mainland | not_significant | called_bh_only | 2.566 | 1.535 | 4.289 | 2.039e-04 | 1.115e-05 |
| RPH3A | low_moderate_high | narrow_mainland | not_significant | called_bh_only | 2.199 | 1.307 | 3.701 | 0.002444 | 1.802e-05 |
| RPH3A | low_moderate_high | intermediate_mainland | not_significant | called_bonferroni | 2.314 | 1.405 | 3.809 | 7.233e-04 | 1.132e-06 |
| RPH3A | low_moderate_high | full_mainland | not_significant | called_bonferroni | 2.314 | 1.428 | 3.751 | 4.792e-04 | 5.445e-07 |
| RPH3A | moderate_high | narrow_mainland | not_significant | called_bh_only | 2.615 | 1.475 | 4.637 | 6.934e-04 | 1.972e-05 |
| RPH3A | moderate_high | intermediate_mainland | not_significant | called_bonferroni | 2.826 | 1.634 | 4.888 | 1.131e-04 | 7.950e-07 |
| RPH3A | moderate_high | full_mainland | not_significant | called_bonferroni | 2.715 | 1.606 | 4.589 | 1.105e-04 | 3.908e-07 |
| JAK2 | low_moderate_high | narrow_mainland | not_significant | not_significant | 1.175 | 0.685 | 2.017 | 0.5572 | 8.648e-04 |
| JAK2 | low_moderate_high | intermediate_mainland | not_significant | not_significant | 1.045 | 0.6221 | 1.756 | 0.8675 | 5.439e-04 |
| JAK2 | low_moderate_high | full_mainland | not_significant | not_significant | 1.048 | 0.6315 | 1.739 | 0.8564 | 4.819e-04 |
| JAK2 | moderate_high | narrow_mainland | not_significant | called_bh_only | 1.707 | 0.7868 | 3.702 | 0.1717 | 5.985e-06 |
| JAK2 | moderate_high | intermediate_mainland | not_significant | called_bh_only | 1.425 | 0.6864 | 2.957 | 0.3399 | 5.610e-06 |
| JAK2 | moderate_high | full_mainland | not_significant | called_bonferroni | 1.564 | 0.7657 | 3.194 | 0.2162 | 1.314e-06 |

**Variants drawn in (a) (from variant_detail.tsv)**

| variant_id | impact | effect | n_carrier_case | n_carrier_control | mac_case | mac_control | n_hom_case | n_hom_control | gene_symbol |
|---|---|---|---|---|---|---|---|---|---|
| chr4:76306990:G:C | HIGH | splice_donor_variant&intron_variant | 4 | 3 | 4 | 3 | 0 | 0 | STBD1 |
| chr4:76307143:C:T | LOW | sequence_feature | 0 | 2 | 0 | 2 | 0 | 0 | STBD1 |
| chr4:76307188:G:C | LOW | sequence_feature | 0 | 5 | 0 | 5 | 0 | 0 | STBD1 |
| chr4:76307249:G:T | LOW | sequence_feature | 0 | 2 | 0 | 2 | 0 | 0 | STBD1 |
| chr4:76307817:C:A | LOW | sequence_feature | 1 | 3 | 1 | 3 | 0 | 0 | STBD1 |
| chr4:76308731:T:C | LOW | sequence_feature | 2 | 2 | 2 | 2 | 0 | 0 | STBD1 |
| chr4:76309487:C:G | MODERATE | missense_variant | 2 | 3 | 2 | 3 | 0 | 0 | STBD1 |
| chr4:76309559:G:T | MODERATE | missense_variant | 4 | 1 | 4 | 1 | 0 | 0 | STBD1 |
| chr11:59576952:G:A | MODERATE | missense_variant | 5 | 5 | 5 | 5 | 0 | 0 | OSBP |
| chr11:59594182:C:T | MODERATE | missense_variant | 3 | 3 | 3 | 3 | 0 | 0 | OSBP |
| chr12:18282395:C:T | MODERATE | missense_variant | 0 | 4 | 0 | 4 | 0 | 0 | PIK3C2G |
| chr12:18282639:C:A | MODERATE | missense_variant | 9 | 10 | 9 | 10 | 0 | 0 | PIK3C2G |
| chr12:18286861:C:G | LOW | synonymous_variant | 2 | 0 | 2 | 0 | 0 | 0 | PIK3C2G |
| chr12:18290996:T:C | LOW | synonymous_variant | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18314059:C:T | HIGH | stop_gained | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18362773:G:A | LOW | synonymous_variant | 0 | 3 | 0 | 3 | 0 | 0 | PIK3C2G |
| chr12:18362887:G:A | HIGH | splice_donor_variant&intron_variant | 1 | 10 | 1 | 10 | 0 | 0 | PIK3C2G |
| chr12:18391228:T:C | MODERATE | missense_variant | 6 | 15 | 6 | 15 | 0 | 0 | PIK3C2G |
| chr12:18391243:C:T | MODERATE | missense_variant | 1 | 1 | 1 | 1 | 0 | 0 | PIK3C2G |
| chr12:18399674:G:T | MODERATE | missense_variant | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18424028:T:C | LOW | synonymous_variant | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18496162:A:G | LOW | splice_region_variant&intron_variant | 2 | 2 | 2 | 2 | 0 | 0 | PIK3C2G |
| chr12:18503281:G:T | MODERATE | missense_variant&splice_region_variant | 0 | 4 | 0 | 4 | 0 | 0 | PIK3C2G |
| chr12:18538166:C:T | MODERATE | missense_variant | 1 | 1 | 1 | 1 | 0 | 0 | PIK3C2G |
| chr12:18546418:A:C | LOW | synonymous_variant | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18566953:C:A | MODERATE | missense_variant | 1 | 4 | 1 | 4 | 0 | 0 | PIK3C2G |
| chr12:18567027:G:A | LOW | synonymous_variant | 0 | 2 | 0 | 2 | 0 | 0 | PIK3C2G |
| chr12:18609627:A:C | MODERATE | missense_variant&splice_region_variant | 2 | 2 | 2 | 2 | 0 | 0 | PIK3C2G |
| chr12:18640421:C:T | LOW | splice_region_variant&intron_variant | 0 | 4 | 0 | 4 | 0 | 0 | PIK3C2G |
| chr12:18640498:C:T | MODERATE | missense_variant | 1 | 7 | 1 | 7 | 0 | 0 | PIK3C2G |
| chr12:18647868:G:A | LOW | splice_region_variant&intron_variant | 0 | 3 | 0 | 3 | 0 | 0 | PIK3C2G |
| chr12:112866794:G:A | MODERATE | missense_variant | 3 | 9 | 3 | 9 | 0 | 0 | RPH3A |
| chr12:112869778:G:A | LOW | synonymous_variant | 0 | 5 | 0 | 5 | 0 | 0 | RPH3A |
| chr12:112875092:C:T | MODERATE | missense_variant | 1 | 3 | 1 | 3 | 0 | 0 | RPH3A |
| chr12:112875132:C:T | MODERATE | missense_variant | 0 | 10 | 0 | 10 | 0 | 0 | RPH3A |
| chr12:112875720:C:T | MODERATE | missense_variant | 0 | 2 | 0 | 2 | 0 | 0 | RPH3A |
| chr12:112876865:A:G | LOW | splice_region_variant&synonymous_variant | 1 | 5 | 1 | 5 | 0 | 0 | RPH3A |
| chr12:112883362:C:T | MODERATE | missense_variant | 1 | 3 | 1 | 3 | 0 | 0 | RPH3A |
| chr12:112887868:A:G | MODERATE | missense_variant | 1 | 1 | 1 | 1 | 0 | 0 | RPH3A |
| chr12:112890858:C:T | MODERATE | missense_variant | 16 | 23 | 16 | 23 | 0 | 0 | RPH3A |
| chr12:112894582:C:T | LOW | synonymous_variant | 0 | 2 | 0 | 2 | 0 | 0 | RPH3A |
| chr12:112896647:C:T | LOW | splice_region_variant&intron_variant | 2 | 3 | 2 | 3 | 0 | 0 | RPH3A |
| chr9:5044432:G:A | MODERATE | missense_variant | 1 | 23 | 1 | 23 | 0 | 0 | JAK2 |
| chr9:5065000:G:A | MODERATE | missense_variant | 0 | 3 | 0 | 3 | 0 | 0 | JAK2 |
| chr9:5065040:C:T | MODERATE | missense_variant&splice_region_variant | 0 | 3 | 0 | 3 | 0 | 0 | JAK2 |
| chr9:5073743:C:T | HIGH | structural_interaction_variant | 1 | 2 | 1 | 2 | 0 | 0 | JAK2 |
| chr9:5073770:G:T | MODERATE | missense_variant | 9 | 0 | 9 | 0 | 0 | 0 | JAK2 |
| chr9:5080268:T:C | MODERATE | missense_variant | 0 | 5 | 0 | 5 | 0 | 0 | JAK2 |
| chr9:5080294:G:T | HIGH | structural_interaction_variant | 0 | 3 | 0 | 3 | 0 | 0 | JAK2 |

## How to read it

1. (a) first: is this one variant or many?
2. (b) against (a): which exon does each stem fall in, and is the signal in coding sequence at all?
3. (c) down a block: does adding controls move the OR?
4. (c) across blocks: does widening the set dilute it?
5. The P columns say which statistic called the cell, in that statistic's own state colour; the marker says what the burden estimate is.

## What this figure does *not* establish

- Not replication: the cohorts are nested.
- Not a CTEPH gene: the platform confound is in every cell.
- The OR is the collapsed-carrier burden; it says nothing about any single variant.

## Model

```
rvtest --burden cmc,cmcWald and --kernel skato, logistic, sex + ancestry PCs; gene models from Ensembl 86 (GRCh38) via the shared region_tracks.R
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
