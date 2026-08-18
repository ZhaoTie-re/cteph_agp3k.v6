# Genome-wide scan — full_mainland, recessive (REC) model

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

312,620 of 5,115,116 variants carried a non-null ERRCODE or a degenerate fit and are excluded from every panel and from lambda_GC. This scan is deflated (lambda_GC = 0.654): the recessive coding leaves most variants with too few alt-homozygotes to fit, plink2 drops them on VIF_INFINITE, and what remains is tested against a null that is too narrow. Its peaks are not comparable with the additive ones. The suggestive tier describes the shape of this scan and is not a list of findings; it receives no per-peak follow-up. Only the additive model defines the downstream fan-out — fine-mapping, conditional analysis and regional plots run on additive peaks alone — so dominant and recessive are reported here in full and go no further.

## Values in this rendering

| quantity | value |
|---|---|
| cohort | full_mainland |
| model | recessive |
| variants in file | 5,115,116 |
| excluded on ERRCODE / degenerate fit | 312,620 |
| variants analysed | 4,802,496 |
| lambda_GC | 0.6543 |
| variants P < 5e-8 | 7 |
| genome-wide peaks | 1 |
| suggestive peaks | 28 |
| suggestive peaks labelled | 10 |
| genome-wide loci | TEX41 |
| annotation style | auto |
| repel force | 0.03 |
| smallest P | 2.015e-08 |

## Full statistics

**Genome-wide peaks of this scan**

| peak_id | rsID | Gene | EA | OA | OR | L95 | U95 | P | Case_Genotype_Distribution | Case_EAF | Control_Genotype_Distribution | Control_EAF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rec005_2_144936354 | rs11682379 | TEX41 | T | C | 3.468 | 2.246 | 5.355 | 2.015e-08 | 268/124/35 | 0.2272 | 1795/719/68 | 0.1656 |

**Suggestive peaks (smallest P first)**

| peak_id | rsID | Gene | OR | L95 | U95 | P | MAF |
|---|---|---|---|---|---|---|---|
| rec014_6_74404472 | rs9352110 | CD109 | 2.085 | 1.575 | 2.759 | 2.802e-07 | 0.337 |
| rec011_4_172000448 | rs10520213 | GALNTL6 | 0.3219 | 0.2085 | 0.4969 | 3.089e-07 | 0.01773 |
| rec028_18_61231020 | rs17067980 | CDH20 | 5.012 | 2.699 | 9.309 | 3.343e-07 | 0.1121 |
| rec024_14_87977321 | rs12587119 | GALC | 1.701 | 1.381 | 2.094 | 5.776e-07 | 0.3849 |
| rec002_1_161421885 | rs61298075 | RNU6-481P | 4.747 | 2.56 | 8.805 | 7.729e-07 | 0.115 |
| rec023_14_81130209 | rs11848926 | TSHR | 4.321 | 2.406 | 7.762 | 9.692e-07 | 0.1314 |
| rec027_17_76738241 | rs237054 | MFSD11 | 0.4202 | 0.2962 | 0.5961 | 1.176e-06 | 0.03157 |
| rec006_2_232445059 | rs790041 | ALPI | 0.5237 | 0.4031 | 0.6804 | 1.271e-06 | 0.07801 |
| rec022_13_49895953 | rs9535399 | RNY4P30 | 3.962 | 2.266 | 6.927 | 1.374e-06 | 0.1388 |
| rec012_4_177869169 | rs1037489 | LINC01098 | 0.5882 | 0.4726 | 0.732 | 1.983e-06 | 0.1424 |
| rec019_11_72256820 | rs7942001 | PHOX2A | 2.932 | 1.867 | 4.604 | 3.014e-06 | 0.1744 |
| rec009_3_101923016 | rs1707610 | Y_RNA | 0.6061 | 0.4912 | 0.7479 | 3.033e-06 | 0.2586 |
| rec010_3_125122636 | rs12492332 | SLC12A8 | 1.789 | 1.401 | 2.284 | 3.111e-06 | 0.4286 |
| rec017_9_87981804 | rs584280 | CDK20 | 1.919 | 1.457 | 2.526 | 3.472e-06 | 0.3608 |
| rec018_9_133271745 | rs554833 | ABO | 0.5421 | 0.4184 | 0.7023 | 3.571e-06 | 0.4673 |

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
