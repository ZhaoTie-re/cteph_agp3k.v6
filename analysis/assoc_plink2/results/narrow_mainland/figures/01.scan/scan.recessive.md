# Genome-wide scan — narrow_mainland, recessive (REC) model

**Figure file:** `scan.recessive.png`

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

437,916 of 5,113,626 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. This scan is deflated (lambda_GC = 0.670): the recessive coding leaves most variants with too few alt-homozygotes to fit, plink2 drops them on VIF_INFINITE, and what remains is tested against a null that is too narrow. Its peaks are not comparable with the additive ones. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the additive model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on additive peaks alone — so dominant and recessive are reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | narrow_mainland |
| model | recessive |
| variants in file | 5,113,626 |
| excluded on ERRCODE / degenerate fit | 437,916 |
| variants analysed | 4,675,710 |
| lambda_GC | 0.6699 |
| variants P < 5e-8 | 0 |
| genome-wide peaks | 0 |
| suggestive peaks | 31 |
| suggestive peaks labelled | 10 |
| genome-wide loci | none |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 1.716e-07 |

## Full statistics

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| rec023_14_87977321 | rs12587119 | GALC | 1.814 | 1.451 | 2.267 | 1.716e-07 | 0.3849 |
| rec031_18_61231020 | rs17067980 | CDH20 | 7.167 | 3.381 | 15.19 | 2.770e-07 | 0.1108 |
| rec005_2_144967963 | rs11674165 | TEX41 | 3.001 | 1.949 | 4.622 | 6.073e-07 | 0.201 |
| rec007_3_101923016 | rs1707610 | Y_RNA | 0.573 | 0.4575 | 0.7177 | 1.248e-06 | 0.2565 |
| rec019_12_129771981 | rs117390528 | TMEM132D | 6.064 | 2.9 | 12.68 | 1.668e-06 | 0.1144 |
| rec028_17_76746460 | rs237051 | MFSD11 | 0.3963 | 0.2714 | 0.5789 | 1.682e-06 | 0.03187 |
| rec025_15_45252090 | rs2413775 | SLC28A2 | 4.111 | 2.282 | 7.407 | 2.515e-06 | 0.1568 |
| rec002_1_61169799 | rs17121764 | NFIA | 4.367 | 2.361 | 8.078 | 2.631e-06 | 0.1449 |
| rec014_8_41305281 | rs79465702 | SFRP1 | 3.053 | 1.915 | 4.865 | 2.689e-06 | 0.2126 |
| rec013_6_165529675 | rs576853 | PDE10A | 1.719 | 1.37 | 2.157 | 2.926e-06 | 0.2771 |
| rec015_9_87981804 | rs584280 | CDK20 | 2.05 | 1.517 | 2.77 | 3.021e-06 | 0.3588 |
| rec020_13_100936730 | rs1570898 | NALCN-AS1 | 1.804 | 1.407 | 2.312 | 3.260e-06 | 0.4754 |
| rec021_14_34275011 | rs1594548173 | EGLN3 | 3.114 | 1.928 | 5.032 | 3.464e-06 | 0.1904 |
| rec022_14_81136689 | rs12891336 | TSHR | 3.735 | 2.136 | 6.531 | 3.805e-06 | 0.1544 |
| rec003_1_161421885 | rs61298075 | RNU6-481P | 5.27 | 2.585 | 10.74 | 4.801e-06 | 0.1215 |

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
