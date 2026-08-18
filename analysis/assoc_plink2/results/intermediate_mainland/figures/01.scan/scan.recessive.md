# Genome-wide scan — intermediate_mainland, recessive (REC) model

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

388,096 of 5,113,982 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. This scan is deflated (lambda_GC = 0.665): the recessive coding leaves most variants with too few alt-homozygotes to fit, plink2 drops them on VIF_INFINITE, and what remains is tested against a null that is too narrow. Its peaks are not comparable with the additive ones. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the additive model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on additive peaks alone — so dominant and recessive are reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | intermediate_mainland |
| model | recessive |
| variants in file | 5,113,982 |
| excluded on ERRCODE / degenerate fit | 388,096 |
| variants analysed | 4,725,886 |
| lambda_GC | 0.665 |
| variants P < 5e-8 | 1 |
| genome-wide peaks | 1 |
| suggestive peaks | 33 |
| suggestive peaks labelled | 10 |
| genome-wide loci | GALC |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 4.108e-08 |

## Full statistics

**Genome-wide peaks of this scan**

| peak_id | rsID | Gene | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Control_Genotype_Distribution | Control_EAF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rec029_14_87977321 | rs12587119 | GALC | G | A | 1.827 | 1.473 | 2.266 | 4.108e-08 | 54/156/212 | 0.6872 | 338/943/706 | 0.5926 |

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| rec008_2_144967963 | rs11674165 | TEX41 | 2.932 | 1.956 | 4.394 | 1.873e-07 | 0.1984 |
| rec004_1_161421885 | rs61298075 | RNU6-481P | 6.254 | 3.121 | 12.53 | 2.356e-07 | 0.1166 |
| rec033_17_76746460 | rs237051 | MFSD11 | 0.3937 | 0.2743 | 0.5653 | 4.371e-07 | 0.03149 |
| rec017_6_74501830 | rs6904268 | COL12A1 | 2.597 | 1.78 | 3.788 | 7.319e-07 | 0.2394 |
| rec014_4_172000448 | rs10520213 | GALNTL6 | 0.3183 | 0.2016 | 0.5027 | 9.119e-07 | 0.01823 |
| rec034_18_61231020 | rs17067980 | CDH20 | 5.463 | 2.77 | 10.77 | 9.578e-07 | 0.1115 |
| rec018_6_165529675 | rs576853 | PDE10A | 1.722 | 1.383 | 2.144 | 1.178e-06 | 0.2794 |
| rec020_9_87981804 | rs584280 | CDK20 | 2.023 | 1.519 | 2.693 | 1.408e-06 | 0.3604 |
| rec007_2_127146116 | rs2118507 | CYP27C1 | 0.5894 | 0.475 | 0.7313 | 1.570e-06 | 0.2412 |
| rec021_9_133268885 | rs597974 | ABO | 0.5227 | 0.4006 | 0.682 | 1.766e-06 | 0.4672 |
| rec011_3_101923016 | rs1707610 | Y_RNA | 0.5924 | 0.4769 | 0.736 | 2.247e-06 | 0.2602 |
| rec009_2_232445059 | rs790041 | ALPI | 0.5239 | 0.4004 | 0.6855 | 2.451e-06 | 0.0791 |
| rec028_14_81136689 | rs12891336 | TSHR | 3.534 | 2.084 | 5.994 | 2.799e-06 | 0.1541 |
| rec027_14_34275011 | rs1594548173 | EGLN3 | 2.984 | 1.879 | 4.739 | 3.602e-06 | 0.1874 |
| rec023_12_129771981 | rs117390528 | TMEM132D | 5.116 | 2.564 | 10.21 | 3.625e-06 | 0.1134 |

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
