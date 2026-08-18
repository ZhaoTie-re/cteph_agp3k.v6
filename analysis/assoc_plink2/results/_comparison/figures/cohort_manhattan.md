# Additive scan across three nested sample sets — narrow, intermediate, full

**Figure file:** `cohort_manhattan.png`

## The question this figure answers

Which associations are shared across the three nested sample sets, and which appear only as the ancestry filter is relaxed?

## Panels

**(a) narrow_mainland — Manhattan and QQ**

LEFT — every analysed variant of the narrow_mainland additive scan at its genomic position, y = -log10 P, on the shared axis; drawn without thinning and rasterised. The box is 4.08 x 1.77 in, the same box scan.additive.png gives its Manhattan, so a peak has the same shape on both. RIGHT — the same P-values ranked against a uniform null, with the pointwise 95% concentration band from the order-statistic Beta(i, n-i+1) distribution and this cohort's lambda_GC. All three QQs carry one limit, for the same reason the Manhattans do. The header counts and lambda_GC come from this cohort's own 02.scan/scan_qc.tsv rather than being recomputed, so they agree with scan.additive.png and cohort_compare.png exactly. Chromosome numbers are on the bottom row, which this row shares.

**(b) intermediate_mainland — Manhattan and QQ**

LEFT — every analysed variant of the intermediate_mainland additive scan at its genomic position, y = -log10 P, on the shared axis; drawn without thinning and rasterised. The box is 4.08 x 1.77 in, the same box scan.additive.png gives its Manhattan, so a peak has the same shape on both. RIGHT — the same P-values ranked against a uniform null, with the pointwise 95% concentration band from the order-statistic Beta(i, n-i+1) distribution and this cohort's lambda_GC. All three QQs carry one limit, for the same reason the Manhattans do. The header counts and lambda_GC come from this cohort's own 02.scan/scan_qc.tsv rather than being recomputed, so they agree with scan.additive.png and cohort_compare.png exactly. Chromosome numbers are on the bottom row, which this row shares.

**(c) full_mainland — Manhattan and QQ**

LEFT — every analysed variant of the full_mainland additive scan at its genomic position, y = -log10 P, on the shared axis; drawn without thinning and rasterised. The box is 4.08 x 1.77 in, the same box scan.additive.png gives its Manhattan, so a peak has the same shape on both. RIGHT — the same P-values ranked against a uniform null, with the pointwise 95% concentration band from the order-statistic Beta(i, n-i+1) distribution and this cohort's lambda_GC. All three QQs carry one limit, for the same reason the Manhattans do. The header counts and lambda_GC come from this cohort's own 02.scan/scan_qc.tsv rather than being recomputed, so they agree with scan.additive.png and cohort_compare.png exactly. Chromosome numbers are on this row, shared by all three.

## Interpretation

The three sample sets are NESTED — narrow within intermediate within full — so they share most of their samples and agreement between rows is expected rather than informative. This figure is not replication and must not be read as such. What it does show is where a locus sits in every set at once: a dotted guide marks each locus that is genome-wide in any cohort, and a diamond marks it only in the cohorts that actually called it, so a column that rises across the rows is a signal gaining sample size rather than a new finding. All three rows are drawn on one genomic axis, one lambda-independent y-limit and one threshold pair, and on equal data heights, so peak heights are directly comparable between rows.

## Values in this rendering

| quantity | value |
|---|---|
| model | additive |
| cohorts | narrow_mainland, intermediate_mainland, full_mainland |
| narrow_mainland: cases | 419 |
| intermediate_mainland: cases | 429 |
| full_mainland: cases | 439 |
| narrow_mainland: controls | 1,749 |
| intermediate_mainland: controls | 2,051 |
| full_mainland: controls | 2,632 |
| narrow_mainland: N_eff | 1,352 |
| intermediate_mainland: N_eff | 1,419 |
| full_mainland: N_eff | 1,505 |
| narrow_mainland: lambda_GC | 1.12 |
| intermediate_mainland: lambda_GC | 1.121 |
| full_mainland: lambda_GC | 1.132 |
| narrow_mainland: variants analysed | 5,113,526 |
| intermediate_mainland: variants analysed | 5,113,902 |
| full_mainland: variants analysed | 5,115,079 |
| narrow_mainland: excluded on ERRCODE / degenerate fit | 100 |
| intermediate_mainland: excluded on ERRCODE / degenerate fit | 80 |
| full_mainland: excluded on ERRCODE / degenerate fit | 37 |
| narrow_mainland: genome-wide peaks | 0 |
| intermediate_mainland: genome-wide peaks | 1 |
| full_mainland: genome-wide peaks | 2 |
| shared -log10 P ceiling | 8.573 |
| loci named (union, genome-wide in any cohort) | MIR3681HG, FTO |
| annotation style | auto |
| repel force | 0.03 |

## Full statistics

**Loci genome-wide in at least one cohort**

| label | chrom | pos | P | smallest_P_in |
|---|---|---|---|---|
| MIR3681HG | 2 | 11,891,017 | 4.342e-08 | full_mainland |
| FTO | 16 | 53,887,925 | 1.154e-08 | intermediate_mainland |

## How to read it

1. Read DOWN a guide, not across a row. The guide is the same genomic position in all three panels, so it shows what one locus does as the ancestry filter is relaxed.
2. A diamond means that cohort called the locus genome-wide. A guide with no diamond means the locus is there but below the threshold in that cohort — MIR3681HG in narrow_mainland is the case here, at P = 4.85e-05.
3. Compare peak heights between rows directly: the y-limit and the axes height are shared, so a taller column is a smaller P and not a rescaled panel.
4. Read the header numbers together with the row, and the QQ beside it. lambda_GC spans 1.120-1.132 across the 3 sets; if that span is narrow and the QQ curves have the same shape, the calibration does not change as the filter is relaxed, so a column that grows down the figure is growing with N and not with inflation.

## What this figure does *not* establish

- The three sample sets are nested, so this is not replication and no independent cohort is available to this study. A locus present in all three rows has been seen once.
- Only the additive model is drawn. Peaks called under the other genetic models are on their own scan figures and are absent here by construction.
- Suggestive peaks are not marked. The dashed line shows the tier; naming them across 3 rows would put 128 labels on the figure.
- Variants plink2 could not fit cleanly are absent from every row; their counts and error codes are in each cohort's 02.scan/scan_qc.tsv.

## Symbols

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **N_eff** — 4/(1/N_case+1/N_ctrl), the balanced-design equivalent sample size. Reported as a cohort descriptor — it is the size any later meta-analysis or mixed-model run quotes — and it collapses toward 4x the smaller arm as the design becomes unbalanced.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^10gamma_k,PC_k,i      (g_i = genotype under the stated model; PCs = bbj_mainland)
lambda_GC = median(chi^2),/,chi^2_1,0.5,   chi^2_1,0.5 = 0.4549
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
