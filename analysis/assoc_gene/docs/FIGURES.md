# Figures — `assoc_gene`

The figure system of this component: five families, read in the order the analysis is read, every
family documented once. The shared `plot_style` / `figure_doc` library is used as a toolkit (colour
layers, the caption block, the collision helpers); the conventions below are this component's own.

## Where they are

```
results/figures/
├── 01.variant_sets/   variant_sets.<cohort>.png ×3 · tables/<cohort>/   README.md
├── 02.gene_scan/      gene_scan.<cohort>.<stratum>.png ×6       README.md
├── 03.calibration/    calibration.png                            README.md
├── 04.robust_genes/   robust_genes.png                           README.md      ← main figure
└── 05.genes/          <GENE>.png (Tier 1 and 2) · index.tsv      README.md
```

**One `README.md` per family.** A family whose members differ only in their numbers (01, 02, 05)
gets one document with one table row per member; a single-figure family's README is its document.
No other Markdown sits beside a PNG. Members of 01 and 02 drop a `<name>.stats.json` that the
`CATALOGUE` step tabulates; 05 writes its README and `index.tsv` itself.

Suggested manuscript roles: **04** and **05** main figures; **01–03** supplementary. The directory
numbers are the reading order, not manuscript numbers.

## The plot block is a fixed physical size

Every figure is authored at journal double-column width, 7.20 in, with the `paper` type scale
(8 pt labels, 7 pt ticks, panel letters 9 pt, nothing below 6.5 pt). `caption_block(plot_h=…)`
fixes the plot block; the canvas grows downward for the three-line caption (finding · one line per
panel · `Data:` line). Crop `plot_h` inches off the top and you have the publication figure.

| family | layout | plot block (in) |
|---|---|---|
| `variant_sets` | 3 × 2: (a) donut, (c) gene sets per stratum, (d) set-size survival down the left; (b) consequences by impact spanning the right column | 7.20 × 6.6 |
| `gene_scan` | (a) CMC and (b) SKAT-O Manhattans full width on one shared genomic axis and one y-limit; (c) square QQ with the null band; (d) square CMC-against-SKAT-O | 7.20 × **computed**: both Manhattan data boxes are fixed at 1.15 in and the canvas grows by the measured gene-label bands (7.4 in with short bands, 8.2 in with the longest) |
| `calibration` | 2 variant-set rows × 3 cohort columns of square QQs | 7.20 × 5.0 |
| `robust_genes` | three row-aligned panels on **one row per (gene × cohort)**: (a) evidence matrix + tier columns, (b) OR forest, (c) carrier dumbbells. Every row carries its own 95 % CI and its own denominator | 7.20 × **computed**: 0.21 in per row × 3 cohorts per gene + header, legend and axis slots (11 genes = 33 rows ≈ 8.3 in). Capped at `MAX_PLOT_H`; the overflow is reported in the caption |
| `genes/<GENE>` | (a) **three strips, one per cohort**, of per-variant carrier **percentage** on one shared y scale (case up, control down) and (b) Ensembl 86 gene models, all four sharing one genomic axis, (c) a table-style OR forest with two P-value columns; a header row names the gene and its tiers inside the crop | 7.20 × 7.0 + 0.24 header |

## Vocabulary

Defined once in `scripts/vocab.py` and imported by every figure:

| thing | how it is written |
|---|---|
| statistics | `CMC burden` (`SERIES[0]`, blue), `SKAT-O` (`SERIES[1]`, red) |
| variant sets | `LOW+MODERATE+HIGH`, `MODERATE+HIGH` — never abbreviated; where both appear in one panel the narrower set is the filled marker and the wider the hollow one |
| cohorts | `narrow`, `intermediate`, `full` — the run's tags with their shared suffix stripped (`plot_style.shorten`), never typed by hand. The cohort is always given **position**: a row in `robust_genes`, a strip in `genes/<GENE>` (a). Where hue is also free it is reinforced by **`plot_style.COHORT_RAMP`**, narrowest darkest; it is never carried by marker size, and never by a channel that already means something else |
| cell states | `Bonferroni + BH` (● `ACCENT`), `BH only` (● `ACCENT_LT`), `tested, not called` (· `DATA`), `not in map` (× `REFERENCE`). **× means `not in map` in every panel of a figure**, and the key lists only the states that figure actually draws — `not_in_map` is usually absent from `robust_genes`, because a gene missing from a cohort's map cannot be called there and so rarely reaches a tier; the README carries the count either way. A cell that WAS tested but has no effect estimate — one group holds no minor allele, so rvtest's Wald fit is unbounded — is `n.e.`, never ×: it keeps its state in the evidence matrix, and only the effect-size panel has nothing to draw |
| groups | cases `ACCENT` (red), controls `DATA` (blue), wherever the two are contrasted |
| tier | the **robustness tier** 1 / 2 / 3 from `robust_genes.tier_of` (0 = tested, called in fewer than two cohorts), one grey channel throughout: the boxed number, the bar in the matrix's right gutter and the row band all deepen with the tier (`TIER_FILL`, `TIER_BAND`). `plot_style.TIER_STYLE` and `tier_lines` are the variant components' genome-wide / suggestive tiers and are never used here |
| impact | snpEff's own colours, in the variant-set figure only |
| λ | `λ_GC` is the only calibration number drawn; its split by set size (`n_var = 2` against `≥ 3`) is in `scan_qc.tsv` and the README tables |
| gene names on a scan | the label set **is** the rejection set: every gene with `significant_bh == 1` and no other, capped by `MaxLabelGenes`; a panel that called nothing names nothing and says so where the BH line would be. Always written vertically (`rotation=90`), so every panel of the family is labelled the same way |
| P-value text | coloured by that statistic's own state in the cell, the same four colours as the cell states; stated in every caption that uses it |
| gene models | Ensembl 86 (GRCh38), one representative transcript per gene, drawn by the shared `plot_regional.draw_gene_track` from tables `region_tracks.R` writes |

## What each figure answers

| family | the one question |
|---|---|
| `variant_sets` | of the minAC-passing call set, what enters a gene set, under which consequence, and how many sets each stratum can test — the **denominator**, before any test |
| `gene_scan` | where does each scan show association, are both statistics calibrated, and do they agree about the genes they single out |
| `calibration` | is the null calibrated in every cohort × variant-set cell, and does that change as controls are added or the filter loosens |
| `robust_genes` | which genes survive a change of statistic and the addition of controls, how large is their burden effect, and how many people carry it |
| `genes/<GENE>` | which variants make up a robust gene's set, where are its carriers, and does the effect hold in every cell |

## Collisions are prevented by measurement

Every script follows one order: draw → `despine` → limits → `caption_block` → then only the helpers
that measure the rendered figure (`gene_labels`, `equalise_row_heights`, `legend_above` /
`place_legend`, `value_labels`, `spread_labels`, `panel_tag(pad=…)`, `thin_tick_labels` last) →
`savefig` without a tight bbox. Legends are anchored to an axes and placed where they cannot cover
data (above the panel, or in a corner that is empty by construction — the lower right of a QQ, the
upper right of a survival curve). Gene labels live in a strip above the Manhattan, never inside it;
the right margin of the scan figure is widened so the last label cannot clip. The main figure's
height is computed from its row count rather than crushing rows into a fixed canvas. In the
variant-set figure the impact table's width is a hard constraint measured before the ring is sized,
the ring drops below the table when the panel is too narrow to hold both side by side, and the
column gap is re-measured after the layout is final so the consequence panel's long tick labels
cannot reach the panel beside it, and two compound consequences that would shorten to the same
label are spelled out to their first differing component. In the per-gene figures the lollipop
and the gene model panel are one nested block whose internal gap holds only a panel letter,
so the model panel is not stranded in the white space the forest below it needs.

## Deliberate omissions

* No `−log10 P`-against-set-size panel: Spearman ρ is a number, tabulated in the gene-scan README.
* No threshold line on a QQ: a threshold is a decision about genes and belongs on the scan.
* No λ split in figures (see Vocabulary).
* No per-gene sidecars, no per-gene figures for Tier 3: the evidence matrix shows every tiered gene;
  the detail figures are for the genes the summary table reports.
* No "cells calling / cells tested" bars and no denominator bar panel: both are columns of
  `tiers.tsv` and `denominators.tsv`, tabulated in the robust-genes README.
* No stratum-membership fill encoding in the lollipop: membership is a function of impact by
  construction, and the impact marker already says it.
* No neighbouring genes in the gene model panel: only the tested gene and genes overlapping it. The
  figure is about one gene, and neighbours cost rows without adding to it.
* No gene names in the scan figure's (d): they are in (a) and (b) against the genome and the
  thresholds; (d) marks the called genes in red and counts them.
