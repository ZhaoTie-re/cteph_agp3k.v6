# Genome-wide scan — intermediate_mainland, main model

**Figure file:** `scan.main.png`

## The question this figure answers

Where does this scan show association, what are its peaks, and is it calibrated well enough for those signals to be believed?

## Panels

**(a) Manhattan plot**

Every analysed variant at its genomic position, y = -log10 P. All points are drawn — no thinning — and the scatter is rasterised so the file stays small while axes and text remain vector. Thinning would distort the QQ plot in (c), which is read from the same P-values. Solid red line is genome-wide significance (5e-8); dashed grey is the suggestive threshold (1e-5). The x tick labels are on (b), which shares this axis.

**(b) Peaks of this scan**

One marker per peak, at its lead variant, on the same x-axis as (a). Peaks are formed by merging variants that pass P < 1e-5 and lie within 250 kb; the lead is the smallest-P variant in the merged window. Distance merging is used rather than LD clumping so the peak definition needs no reference panel and cannot shift with the choice of LD sample. Labels: every genome-wide lead, plus the N smallest-P suggestive leads (N is --label-suggestive, default 10).

**(c) Quantile-quantile plot**

The same P-values ranked against their expectation under a uniform null. The grey band is the pointwise 95% concentration from the order-statistic Beta(i, n-i+1) distribution, so departures are judged against a stated reference rather than by eye. lambda_GC = median(chi2)/0.4549 is printed in the panel.

**(d) Effect against frequency**

MAF of each lead on x; |effect| as an odds ratio, exp|log OR|, on a log y axis. Taking the absolute value puts protective and risk alleles on one scale. The upward slope toward low MAF is a property of threshold selection, not of biology.

## Interpretation

0 of 5,113,982 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the peak model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on its peaks alone — so every other model is reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| model | main |
| variants in file | 5,113,982 |
| excluded on ERRCODE / degenerate fit | 0 |
| variants analysed | 5,113,982 |
| lambda_GC | 1.038 |
| variants P < 5e-8 | 2 |
| genome-wide peaks | 1 |
| suggestive peaks | 27 |
| suggestive peaks labelled | 10 |
| genome-wide loci | FTO |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 3.836e-08 |

## Full statistics

**Genome-wide peaks of this scan**

| peak_id | rsID | Gene | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Control_Genotype_Distribution | Control_EAF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mai024_16_53882967 | rs9934504 | FTO | A | G | 1.679 | 1.386 | 2.033 | 3.836e-08 | 217/163/41 | 0.291 | 1259/651/87 | 0.2066 |

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| mai008_3_94020684 | rs143483852 | ARL13B | 6.078 | 2.938 | 12.57 | 6.839e-08 | 0.01495 |
| mai025_17_13528059 | rs34804183 | HS3ST3A1 | 2.548 | 1.735 | 3.74 | 4.175e-07 | 0.04526 |
| mai023_14_88019710 | rs4899932 | LINC01146 | 1.503 | 1.277 | 1.768 | 4.406e-07 | 0.3985 |
| mai011_4_171998307 | rs75214171 | GALNTL6 | 4.484 | 2.369 | 8.488 | 4.801e-07 | 0.01715 |
| mai012_5_66174098 | rs4700094 | SREK1 | 0.6248 | 0.5159 | 0.7568 | 5.964e-07 | 0.2239 |
| mai013_6_32632861 | rs9272130 | HLA-DQA1 | 0.6638 | 0.561 | 0.7853 | 8.570e-07 | 0.373 |
| mai026_17_76739260 | rs9897202 | MFSD11 | 2.976 | 1.861 | 4.758 | 1.076e-06 | 0.0309 |
| mai010_3_154069965 | rs1021580972 | ARHGEF26-AS1 | 3.904 | 2.148 | 7.097 | 1.214e-06 | 0.01854 |
| mai015_8_132794261 | rs150804769 | PHF20L1 | 3.78 | 2.082 | 6.862 | 1.915e-06 | 0.02113 |
| mai027_19_10631611 | rs1560711 | SLC44A2 | 1.494 | 1.259 | 1.773 | 2.029e-06 | 0.3449 |
| mai009_3_100783010 | rs7648469 | ABI3BP | 1.857 | 1.421 | 2.426 | 2.095e-06 | 0.09769 |
| mai005_2_232438125 | rs790023 | ALPI | 0.499 | 0.3683 | 0.6761 | 2.386e-06 | 0.07516 |
| mai017_9_104559024 | rs7036847 | OR13C8 | 1.608 | 1.31 | 1.975 | 2.418e-06 | 0.188 |
| mai028_20_60915853 | rs79313991 | CDH4 | 2.908 | 1.802 | 4.69 | 2.744e-06 | 0.02911 |
| mai016_9_71722660 | rs201776247 | TMEM2 | 2.225 | 1.56 | 3.175 | 2.859e-06 | 0.05403 |

## How to read it

1. Check lambda_GC in (b) first. Near 1 means the bulk of the distribution is calibrated; a scan far from 1 in either direction is discussed in METHODS §7.
2. In (b), a curve that lifts off the band only at the extreme right is what a true signal looks like. A curve that leaves the band along its whole length indicates residual structure rather than a handful of associations.
3. Read (c) against (a): a peak marker sits directly below the Manhattan column it was called from, so the two can be compared without re-reading positions.
4. Read the peaks as belonging to *this* model only. A locus shown on the recessive figure is a recessive result and is not carried into fine-mapping or conditional analysis, both of which run on the peak model alone.
5. For a deflated scan, treat the peak count as an upper bound: the variants plink2 could fit are a biased subset of the call set.

## What this figure does *not* establish

- It cannot separate confounding from polygenicity. That needs LDSC, which is not usable at this effective sample size, or a random-effect model — see METHODS §7.
- Variants plink2 could not fit cleanly are absent from every panel, so the figure says nothing about them; their counts and error codes are in `02.scan/scan_qc.tsv`.
- A peak here is a threshold crossing in one scan of one nested sample set. It is not a replicated finding, and no independent cohort is available to this study.
- Panel (d) compares effects across frequencies but against no detection reference; this component makes no power claim.

## Symbols

- **genetic model** — the encoding of the genotype in the GLM: ADD counts alt alleles (0/1/2), DOM contrasts carriers against non-carriers, REC contrasts alt-homozygotes against the rest. Three separate genome-wide scans, not one joint test; ADD is the primary.

- **lambda_GC** — genomic-control inflation factor = median chi^2 / 0.4549. lambda>1 is inflation, lambda<1 deflation; neither is corrected for here — it is reported as a calibration read-out, and its interpretation is in METHODS §7.

- **genome-wide significance** — P<5x10^-8, applied identically to every cohort and model. The three cohorts are nested and the three models correlated, so these are not independent tests and no further multiplicity adjustment is made — stated, not silently assumed.

- **OR** — odds ratio per copy of the A1 (effect) allele from logistic regression, with its 95% CI. Compared on the log scale, so a protective and a risk allele of equal strength are equally far from 1.

- **ERRCODE** — plink2's per-variant fit diagnostic. Variants with a non-'.' code (e.g. VIF_INFINITE, SEPARATION) did not fit cleanly and are excluded from lambda_GC and from the hit list rather than silently carried.

## Model

```
logit,Pr(case_i) = beta_0 + beta,g_i + gamma_sex,SEX_i + sum_k=1^10gamma_k,PC_k,i      (g_i = genotype under the stated model; PCs = bbj_mainland)
lambda_GC = median(chi^2),/,chi^2_1,0.5,   chi^2_1,0.5 = 0.4549
```

---

Methods and rationale: [`METHODS.md`](../../../docs/METHODS.md)
