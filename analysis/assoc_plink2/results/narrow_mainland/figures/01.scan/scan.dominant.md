# Genome-wide scan — narrow_mainland, dominant (DOM) model

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

21,340 of 5,113,626 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. lambda_GC = 1.059 is above 1 across the whole distribution, not just at the tail, which is residual structure rather than a handful of true signals; no genomic-control correction is applied to any reported statistic. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the additive model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on additive peaks alone — so dominant and recessive are reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | narrow_mainland |
| model | dominant |
| variants in file | 5,113,626 |
| excluded on ERRCODE / degenerate fit | 21,340 |
| variants analysed | 5,092,286 |
| lambda_GC | 1.059 |
| variants P < 5e-8 | 0 |
| genome-wide peaks | 0 |
| suggestive peaks | 39 |
| suggestive peaks labelled | 10 |
| genome-wide loci | none |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 3.591e-07 |

## Full statistics

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| dom012_5_66174098 | rs4700094 | SREK1 | 0.5358 | 0.4213 | 0.6813 | 3.591e-07 | 0.227 |
| dom038_20_60915853 | rs79313991 | CDH4 | 2.776 | 1.871 | 4.118 | 3.905e-07 | 0.02923 |
| dom035_18_61226084 | rs1845750 | CDH20 | 0.145 | 0.06821 | 0.3083 | 5.219e-07 | 0.1152 |
| dom005_3_94020684 | rs143483852 | ARL13B | 3.739 | 2.22 | 6.299 | 7.184e-07 | 0.0168 |
| dom029_14_88030553 | rs10147393 | LINC01146 | 0.5773 | 0.4624 | 0.7207 | 1.212e-06 | 0.3539 |
| dom033_17_76739260 | rs9897202 | MFSD11 | 2.561 | 1.75 | 3.748 | 1.280e-06 | 0.03163 |
| dom010_5_10794370 | rs55744610 | DAP | 0.5761 | 0.4606 | 0.7205 | 1.351e-06 | 0.3643 |
| dom009_4_171998307 | rs75214171 | GALNTL6 | 3.324 | 2.036 | 5.426 | 1.563e-06 | 0.01786 |
| dom019_8_138573836 | rs56779193 | COL22A1 | 3.099 | 1.95 | 4.925 | 1.705e-06 | 0.02083 |
| dom037_20_7553182 | rs117662178 | RN7SL547P | 3.551 | 2.11 | 5.973 | 1.807e-06 | 0.01563 |
| dom021_9_71722660 | rs201776247 | TMEM2 | 2.132 | 1.561 | 2.911 | 1.931e-06 | 0.05664 |
| dom006_3_100783010 | rs7648469 | ABI3BP | 1.87 | 1.443 | 2.423 | 2.169e-06 | 0.09986 |
| dom028_14_34300615 | rs2383582 | EGLN3 | 0.3353 | 0.2132 | 0.5272 | 2.222e-06 | 0.2076 |
| dom017_8_105989633 | rs117799920 | ZFPM2-AS1 | 2.482 | 1.703 | 3.619 | 2.292e-06 | 0.0342 |
| dom008_3_154069965 | rs1021580972 | ARHGEF26-AS1 | 3.06 | 1.915 | 4.888 | 2.883e-06 | 0.01991 |

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
