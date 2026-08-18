# Genome-wide scan — intermediate_mainland, dominant (DOM) model

**Figure file:** `scan.dominant.png`

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

18,576 of 5,113,982 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. lambda_GC = 1.060 is above 1 across the whole distribution, not just at the tail, which is residual structure rather than a handful of true signals; no genomic-control correction is applied to any reported statistic. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the additive model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on additive peaks alone — so dominant and recessive are reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| model | dominant |
| variants in file | 5,113,982 |
| excluded on ERRCODE / degenerate fit | 18,576 |
| variants analysed | 5,095,406 |
| lambda_GC | 1.06 |
| variants P < 5e-8 | 1 |
| genome-wide peaks | 1 |
| suggestive peaks | 46 |
| suggestive peaks labelled | 10 |
| genome-wide loci | SREK1 |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 4.579e-08 |

## Full statistics

**Genome-wide peaks of this scan**

| peak_id | rsID | Gene | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Control_Genotype_Distribution | Control_EAF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dom020_5_66174098 | rs4700094 | SREK1 | G | A | 0.5214 | 0.4129 | 0.6586 | 4.579e-08 | 297/104/15 | 0.1611 | 1122/736/119 | 0.2463 |

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| dom011_3_94020684 | rs143483852 | ARL13B | 3.964 | 2.39 | 6.575 | 9.555e-08 | 0.0155 |
| dom039_17_13528059 | rs34804183 | HS3ST3A1 | 2.314 | 1.689 | 3.17 | 1.770e-07 | 0.04642 |
| dom040_17_76739260 | rs9897202 | MFSD11 | 2.574 | 1.791 | 3.701 | 3.284e-07 | 0.03134 |
| dom017_4_171998307 | rs75214171 | GALNTL6 | 3.305 | 2.084 | 5.242 | 3.752e-07 | 0.01766 |
| dom035_14_88030553 | rs10147393 | LINC01146 | 0.5747 | 0.4637 | 0.7122 | 4.174e-07 | 0.3604 |
| dom037_16_53887925 | rs16952623 | FTO | 1.729 | 1.391 | 2.149 | 7.805e-07 | 0.2229 |
| dom045_20_60915853 | rs79313991 | CDH4 | 2.556 | 1.755 | 3.723 | 9.889e-07 | 0.02922 |
| dom026_9_71722660 | rs201776247 | TMEM2 | 2.109 | 1.562 | 2.848 | 1.092e-06 | 0.05583 |
| dom012_3_100797923 | rs9881972 | ABI3BP | 1.84 | 1.436 | 2.356 | 1.378e-06 | 0.1022 |
| dom014_3_154069965 | rs1021580972 | ARHGEF26-AS1 | 3.019 | 1.927 | 4.73 | 1.409e-06 | 0.01904 |
| dom004_2_76843478 | rs62173094 | LRRTM4 | 1.718 | 1.376 | 2.144 | 1.703e-06 | 0.1685 |
| dom043_18_61226084 | rs1845750 | CDH20 | 0.1906 | 0.09645 | 0.3765 | 1.829e-06 | 0.1164 |
| dom027_10_24375907 | rs146165269 | KIAA1217 | 2.681 | 1.78 | 4.037 | 2.362e-06 | 0.02527 |
| dom007_2_233068979 | rs72982244 | INPP5D | 3.118 | 1.942 | 5.005 | 2.491e-06 | 0.01687 |
| dom025_8_132794261 | rs150804769 | PHF20L1 | 2.78 | 1.808 | 4.276 | 3.233e-06 | 0.02231 |

## How to read it

1. Check lambda_GC in (b) first. Near 1 means the bulk of the distribution is calibrated; the additive scans run ~1.12 and the recessive ~0.65, which is discussed in METHODS §7.
2. In (b), a curve that lifts off the band only at the extreme right is what a true signal looks like. A curve that leaves the band along its whole length indicates residual structure rather than a handful of associations.
3. Read (c) against (a): a peak marker sits directly below the Manhattan column it was called from, so the two can be compared without re-reading positions.
4. Read the peaks as belonging to *this* model only. A locus shown on the recessive figure is a recessive result and is not carried into fine-mapping or conditional analysis, both of which run on the additive peaks alone.
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
