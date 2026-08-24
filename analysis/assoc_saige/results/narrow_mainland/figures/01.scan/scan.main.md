# Genome-wide scan — narrow_mainland, main model

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

0 of 5,113,626 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the peak model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on its peaks alone — so every other model is reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | narrow_mainland |
| model | main |
| variants in file | 5,113,626 |
| excluded on ERRCODE / degenerate fit | 0 |
| variants analysed | 5,113,626 |
| lambda_GC | 1.039 |
| variants P < 5e-8 | 0 |
| genome-wide peaks | 0 |
| suggestive peaks | 29 |
| suggestive peaks labelled | 10 |
| genome-wide loci | none |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 1.358e-07 |

## Full statistics

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| mai023_16_53882967 | rs9934504 | FTO | 1.671 | 1.372 | 2.036 | 1.358e-07 | 0.2125 |
| mai004_3_94020684 | rs143483852 | ARL13B | 5.248 | 2.587 | 10.65 | 5.773e-07 | 0.01617 |
| mai029_20_60915853 | rs79313991 | CDH4 | 3.091 | 1.899 | 5.031 | 1.341e-06 | 0.02916 |
| mai022_14_88019710 | rs4899932 | LINC01146 | 1.499 | 1.265 | 1.776 | 1.421e-06 | 0.3936 |
| mai008_4_171998307 | rs75214171 | GALNTL6 | 4.386 | 2.276 | 8.453 | 1.734e-06 | 0.01731 |
| mai005_3_100783010 | rs7648469 | ABI3BP | 1.906 | 1.445 | 2.515 | 1.970e-06 | 0.0959 |
| mai015_8_138575524 | rs62530012 | COL22A1 | 3.965 | 2.145 | 7.329 | 2.143e-06 | 0.02005 |
| mai026_19_10631611 | rs1560711 | SLC44A2 | 1.514 | 1.268 | 1.808 | 2.161e-06 | 0.3447 |
| mai028_20_7553182 | rs117662178 | RN7SL547P | 4.434 | 2.268 | 8.665 | 2.611e-06 | 0.01481 |
| mai007_3_154069965 | rs1021580972 | ARHGEF26-AS1 | 3.762 | 2.076 | 6.814 | 2.820e-06 | 0.01936 |
| mai011_6_32632861 | rs9272130 | HLA-DQA1 | 0.6678 | 0.5603 | 0.7958 | 3.173e-06 | 0.3706 |
| mai017_9_71722660 | rs201776247 | TMEM2 | 2.254 | 1.569 | 3.238 | 3.447e-06 | 0.05467 |
| mai018_9_104559024 | rs7036847 | OR13C8 | 1.62 | 1.31 | 2.005 | 4.011e-06 | 0.1893 |
| mai021_11_119708698 | rs117569646 | NECTIN1 | 3.588 | 2.005 | 6.42 | 4.378e-06 | 0.02027 |
| mai025_17_76739260 | rs9897202 | MFSD11 | 2.847 | 1.773 | 4.571 | 4.416e-06 | 0.03121 |

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
