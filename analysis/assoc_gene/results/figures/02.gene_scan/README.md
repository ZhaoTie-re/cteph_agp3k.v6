# Gene-based scans — CMC burden and SKAT-O, every cohort x stratum

**Figures:** `02.gene_scan/` — 6 figure(s): `gene_scan.full_mainland.low_moderate_high`, `gene_scan.full_mainland.moderate_high`, `gene_scan.intermediate_mainland.low_moderate_high`, `gene_scan.intermediate_mainland.moderate_high`, `gene_scan.narrow_mainland.low_moderate_high`, `gene_scan.narrow_mainland.moderate_high`. One document covers the family; the numbers of each member are one row of the table below.

## The question this figure answers

Where does each scan show association, are the two statistics calibrated, and do they agree about the genes they single out?

## Panels

**(a) CMC burden across the genome**

One point per gene at the midpoint of its mapped variant span, chromosomes banded, y = -log10 of the CMC score-test P (logistic model, sex + ancestry PCs, minor alleles collapsed to carrier / non-carrier). The solid line is the Bonferroni threshold alpha / n_genes_mapped fixed at the map barrier; the dashed line is the Benjamini-Hochberg cut of THIS statistic — the largest P it rejected — over the same gene family; where BH rejected nothing it has no cut, and the panel says so instead of leaving the reader to notice a missing line. Genes called by both rules are drawn in the full accent, genes called by BH only in the lighter one, and EVERY called gene -- and only a called gene -- is named above the panel in the same two tones, written vertically so that every panel of the family is labelled the same way. A scan that called nothing names nothing.

**(b) SKAT-O across the genome**

The same construction over the SKAT-O P-values (same model, same sets, same samples, mean-imputed dosages for the kernel). The panel shares (a)'s x mapping, y-limit and Bonferroni line, so a taller point is a stronger point; the BH cut is this statistic's own. Both Manhattan data boxes are the same fixed height in every member of the family, and the y limit is the largest -log10 P this variant set reached in ANY cohort (scan_scale.tsv), so a point may be compared across the three cohorts of one set; the threshold lines stay each cohort's own.

**(c) Calibration**

Observed against expected -log10 P for both statistics over the genes each returned, with the pointwise Beta(i, n-i+1) 95 % band of the null. The legend carries the genomic-control inflation factor lambda_GC of each statistic, read from scan_qc.tsv. The split of lambda by set size (n_var = 2 against >= 3) is in the calibration table below, not in the figure: for two-variant genes the normal approximation of a discrete statistic is anti-conservative near P = 0.5, which inflates a median-based lambda without touching the tail where calls are made (METHODS §4.4).

**(d) CMC against SKAT-O, gene by gene**

Genes returned by both statistics, matched on set name. A gene is accented when BH calls it under either statistic; the solid lines are the Bonferroni threshold on both axes, the dashed lines each statistic's BH cut. A gene far above the diagonal carries a signal the kernel sees and the collapse does not (a subset of variants, or opposite directions); a gene below it is one where every minor allele points the same way. No gene is named here -- the names are in (a) and (b), against the genome and the thresholds; this panel marks the called genes in red and counts them.

## Interpretation

Each cohort x stratum is one scan judged by one threshold, alpha / n_genes_mapped, fixed before any P existed, and by BH q < alpha over the same family. A gene called by one statistic and not the other is therefore a statement about the statistics, not about two thresholds. Calibration bounds everything: a curve leaving the band early is producing small P-values everywhere. None of this makes a gene a CTEPH gene — cases and controls were sequenced on different platforms with no overlap, so every P carries that confound and a called gene is a candidate. Which candidates survive the change of statistic and the addition of controls is the question of the robust-genes figure.

## Values in this rendering

| quantity | value |
|---|---|
| members (cohort x stratum) | 6 |
| alpha | 0.05 |
| MinNumVar | 2 |

## Full statistics

**Per-scan values**

| figure | cohort | stratum | genes mapped | MinNumVar | Bonferroni threshold | BH cut, CMC | BH cut, SKAT-O | called (BH), CMC | of which Bonferroni, CMC | called (BH), SKAT-O | of which Bonferroni, SKAT-O | called by both | strongest, CMC | strongest, SKAT-O | lambda_GC, CMC | lambda_GC n_var=2, CMC | lambda_GC n_var>=3, CMC | lambda_GC, SKAT-O | lambda_GC n_var=2, SKAT-O | lambda_GC n_var>=3, SKAT-O | Spearman rho set size vs signal, CMC | Spearman rho set size vs signal, SKAT-O | genes named in (a) | genes named in (b) | label rule | name cap reached | genes called by either statistic (d) | gene tests with no P |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gene_scan.full_mainland.low_moderate_high | full_mainland | low_moderate_high | 15,691 | 2 | 3.187e-06 | — | 3.165e-05 | 0 | 0 | 10 | 3 | 0 | CD14 (5e-6) | RPH3A (5e-7) | 1.176 | 1.495 | 1.136 | 1.226 | 1.207 | 1.228 | -0.01375 | 0.02977 |  | SCAMP3, STBD1, CD14, RP11-350J20.5, FAS, PIK3C2G, RPH3A, SCFD1, GABRA5, FAM217B | every called gene, nothing else | 0 | 10 | 0 |
| gene_scan.full_mainland.moderate_high | full_mainland | moderate_high | 12,378 | 2 | 4.039e-06 | 7.226e-08 | 2.433e-05 | 1 | 1 | 11 | 4 | 1 | STBD1 (7e-8) | FAM217B (1e-7) | 1.181 | 1.44 | 1.111 | 1.164 | 1.173 | 1.161 | -0.02136 | 0.01048 | STBD1 | STBD1, JAK2, OSBP, STARD10, P2RY2, PIK3C2G, RPH3A, SCFD1, CA12, CEACAM7, FAM217B | every called gene, nothing else | 0 | 11 | 0 |
| gene_scan.intermediate_mainland.low_moderate_high | intermediate_mainland | low_moderate_high | 15,009 | 2 | 3.331e-06 | 2.759e-06 | 1.495e-05 | 1 | 1 | 7 | 3 | 1 | SCAMP3 (3e-6) | PIK3C2G (3e-7) | 1.185 | 1.504 | 1.139 | 1.246 | 1.371 | 1.229 | -0.01591 | 0.0211 | SCAMP3 | SCAMP3, TBC1D14, STBD1, CD248, PIK3C2G, RPH3A, PES1 | every called gene, nothing else | 0 | 7 | 0 |
| gene_scan.intermediate_mainland.moderate_high | intermediate_mainland | moderate_high | 11,436 | 2 | 4.372e-06 | 1.305e-05 | 2.791e-05 | 3 | 1 | 9 | 3 | 3 | OSBP (4e-7) | OSBP (3e-7) | 1.154 | 1.376 | 1.094 | 1.172 | 1.209 | 1.162 | -0.01246 | 0.00978 | STBD1, OSBP, P2RY2 | STBD1, NKD2, JAK2, OSBP, P2RY2, PIK3C2G, RPH3A, OTOP3, PES1 | every called gene, nothing else | 0 | 9 | 0 |
| gene_scan.narrow_mainland.low_moderate_high | narrow_mainland | low_moderate_high | 14,557 | 2 | 3.435e-06 | — | 3.773e-05 | 0 | 0 | 11 | 2 | 0 | STBD1 (4e-6) | PIK3C2G (8e-7) | 1.164 | 1.367 | 1.126 | 1.237 | 1.275 | 1.227 | -0.003058 | 0.03046 |  | ZKSCAN7, STBD1, NKD2, MYOM2, CD248, PIK3C2G, RPH3A, SCFD1, SMG1, PES1, EFCAB6 | every called gene, nothing else | 0 | 11 | 0 |
| gene_scan.narrow_mainland.moderate_high | narrow_mainland | moderate_high | 10,832 | 2 | 4.616e-06 | — | 2.132e-05 | 0 | 0 | 9 | 1 | 0 | OSBP (5e-6) | PIK3C2G (4e-7) | 1.152 | 1.328 | 1.102 | 1.215 | 1.239 | 1.204 | -0.004302 | 0.01498 |  | STBD1, NKD2, JAK2, OSBP, PIK3C2G, RPH3A, SCFD1, SMG1, PES1 | every called gene, nothing else | 0 | 9 | 0 |

## How to read it

1. Read (a) and (b) against the same solid line, then against each other.
2. Go to (c) before believing either: lambda_GC is reported, not targeted, and the table gives its set-size split.
3. Use (d) to see which genes the two statistics disagree about; those are the genes whose variant-level structure matters.
4. A gene-set P says nothing about which variant carries it; the per-gene figures do.

## What this figure does *not* establish

- It does not establish that any gene is associated with CTEPH (platform confound).
- Point HEIGHTS are comparable across the three cohorts of one variant set, whose y axis is shared (scan_scale.tsv), and across (a) and (b) within a figure. They are NOT comparable between variant sets, and no DECISION is comparable anywhere: each scan has its own denominator, so the threshold lines sit at different heights in different scans.
- The diagonal in (d) is not a null expectation; two statistics on the same data are correlated by construction.

## Model

```
rvtest --burden cmc,cmcWald and --kernel skato, logistic, sex + ancestry PCs, --impute mean; threshold = alpha / n_genes_mapped per cohort x stratum; called = BH q < alpha over the same family
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
