# 2 · Random model with and without PCs

## Slide summary

**Conclusion: pending — to be agreed. The facts:**

- With the full GRM in the model, adding 10 PCs leaves λGC nearly unchanged: 1.040 / 1.024 / 1.017 without → 1.039 / 1.038 / 1.021 with PCs (narrow / intermediate / full mainland).
- The PCs differ between cases and controls in every cohort: likelihood-ratio test of PC1–10 on case status, P = 5.4 × 10⁻¹⁰ / 0.0081 / 6.0 × 10⁻⁴.
- Variants at P < 10⁻⁵ (≈ 51 expected): 180 / 168 / 121 without PCs, 87 / 160 / 111 with PCs.
- Peaks found only without PCs: 8 / 2 / 3, among them in the extended MHC (chr6 28.4–33.5 Mb): *XXbac-BPG308K3.5* (narrow mainland, 1 variant), *DDR1* (narrow mainland, 27 variants), *HCG17* (intermediate mainland, 1 variant); only with PCs: 6 / 3 / 6.
- Effect sizes of the 48 cross-cohort lead × cohort pairs barely move: median |Δ ln OR| = 0.014 (max 0.23).

## Figures

- **`grm_pc_vs_grm.png`** — one row per cohort: Manhattan of the random model without PCs (left) and with PCs (middle, the main analysis) on one axis and scale, and the QQ plot of both (right).

## Detailed description

### Question
Within the random (mixed) model with a full GRM, what do the 10 ancestry PCs add as fixed covariates?

### Design
Both scans are SAIGE logistic mixed models on the same samples (relatives retained), the same full GRM and the
same variants; the only difference is the covariates: SEX + full GRM, against SEX + PC1–10 + full GRM (the main
analysis). P values: score test, normal approximation.

### How to read `grm_pc_vs_grm.png`
As `01.fixed_vs_random`: rows are cohorts, both Manhattans share one axis and one −log₁₀P limit, red line
5 × 10⁻⁸, grey dashed 10⁻⁵, diamonds = a panel's own genome-wide leads, QQ with the 95% band.

### Results
Calibration and peaks (peaks: upstream calls, ±250 kb; "only" = no peak of the other model within 500 kb):

| Cohort | λGC without PCs | λGC with PCs | P<10⁻⁵ without | P<10⁻⁵ with | peaks without / with | only without | only with |
|---|---|---|---|---|---|---|---|
| Narrow mainland | 1.040 | 1.039 | 180 | 87 | 31 / 29 | 8 | 6 |
| Intermediate mainland | 1.024 | 1.038 | 168 | 160 | 27 / 28 | 2 | 3 |
| Full mainland | 1.017 | 1.021 | 121 | 111 | 29 / 32 | 3 | 6 |

Do the PCs differ between cases and controls? Logistic regression of case status on SEX vs SEX + PC1–10, samples
of the mixed-model scans:

| Cohort | cases / controls | LR χ²(10) | P | PCs with P < 0.05 | pseudo-R² gain |
|---|---|---|---|---|---|
| Narrow mainland | 419/1,776 | 64.4 | 5.4 × 10⁻¹⁰ | PC1, PC2, PC6, PC10 | 0.0301 |
| Intermediate mainland | 429/2,079 | 23.8 | 0.0081 | PC2, PC6, PC10 | 0.0104 |
| Full mainland | 439/2,662 | 31.0 | 6.0 × 10⁻⁴ | PC3, PC4, PC10 | 0.0122 |

Peaks called under only one of the two models:

| Cohort | only in | Gene | Position | rsID | variants | P |
|---|---|---|---|---|---|---|
| Narrow mainland | with PCs | XPR1 | chr1:180,866,004 | rs77630739 | 6 | 4.8 × 10⁻⁶ |
| Narrow mainland | with PCs | MAGI2 | chr7:78,457,650 | rs73150307 | 1 | 8.1 × 10⁻⁶ |
| Narrow mainland | with PCs | PHF20L1 | chr8:132,794,261 | rs150804769 | 1 | 4.5 × 10⁻⁶ |
| Narrow mainland | with PCs | SPINK4 | chr9:33,217,214 | rs830406 | 1 | 9.3 × 10⁻⁶ |
| Narrow mainland | with PCs | INSC | chr11:15,225,340 | rs1714337 | 1 | 8.1 × 10⁻⁶ |
| Narrow mainland | with PCs | NECTIN1 | chr11:119,708,698 | rs117569646 | 1 | 4.4 × 10⁻⁶ |
| Narrow mainland | without PCs | HDAC4 | chr2:239,515,147 | rs7602325 | 3 | 5.7 × 10⁻⁶ |
| Narrow mainland | without PCs | GAP43 | chr3:115,430,183 | rs9859487 | 3 | 5.1 × 10⁻⁶ |
| Narrow mainland | without PCs | XXbac-BPG308K3.5 | chr6:28,788,591 | rs117473648 | 1 | 9.2 × 10⁻⁶ |
| Narrow mainland | without PCs | DDR1 | chr6:30,886,019 | rs3131034 | 27 | 3.6 × 10⁻⁶ |
| Narrow mainland | without PCs | MIR4289 | chr9:88,825,028 | rs143267648 | 1 | 6.4 × 10⁻⁶ |
| Narrow mainland | without PCs | NAALAD2 | chr11:90,195,615 | rs11603498 | 1 | 9.4 × 10⁻⁶ |
| Narrow mainland | without PCs | LINC00598 | chr13:40,072,903 | rs60080071 | 1 | 9.8 × 10⁻⁶ |
| Narrow mainland | without PCs | PSTPIP2 | chr18:46,008,989 | rs79902549 | 1 | 7.5 × 10⁻⁶ |
| Intermediate mainland | with PCs | COL12A1 | chr6:74,508,980 | rs9360810 | 2 | 8.8 × 10⁻⁶ |
| Intermediate mainland | with PCs | PHF20L1 | chr8:132,794,261 | rs150804769 | 1 | 1.9 × 10⁻⁶ |
| Intermediate mainland | with PCs | LINC01309 | chr13:103,722,812 | rs9284206 | 3 | 6.1 × 10⁻⁶ |
| Intermediate mainland | without PCs | HCG17 | chr6:30,228,997 | rs117586454 | 1 | 8.1 × 10⁻⁶ |
| Intermediate mainland | without PCs | RN7SL547P | chr20:7,553,182 | rs117662178 | 1 | 7.9 × 10⁻⁶ |
| Full mainland | with PCs | MTPN | chr7:135,996,912 | rs143188925 | 1 | 9.8 × 10⁻⁶ |
| Full mainland | with PCs | PHF20L1 | chr8:132,794,261 | rs150804769 | 1 | 2.4 × 10⁻⁶ |
| Full mainland | with PCs | COL22A1 | chr8:138,575,524 | rs62530012 | 1 | 5.6 × 10⁻⁶ |
| Full mainland | with PCs | OR13C8 | chr9:104,559,024 | rs7036847 | 1 | 9.7 × 10⁻⁶ |
| Full mainland | with PCs | NECTIN1 | chr11:119,708,698 | rs117569646 | 1 | 9.8 × 10⁻⁶ |
| Full mainland | with PCs | MCTP2 | chr15:94,387,061 | rs41444346 | 1 | 9.6 × 10⁻⁶ |
| Full mainland | without PCs | PAX7 | chr1:18,618,717 | rs4027744 | 1 | 6.2 × 10⁻⁶ |
| Full mainland | without PCs | SREK1 | chr5:66,174,098 | rs4700094 | 1 | 7.6 × 10⁻⁶ |
| Full mainland | without PCs | BUB3 | chr10:123,414,776 | rs545536036 | 2 | 3.9 × 10⁻⁶ |

Lead variants of the 16 cross-cohort loci (genome-wide in any cohort or suggestive in all three) under both
models: median |Δ ln OR| = 0.014, largest 0.23; per variant in `lead_or_stability.tsv`.

λGC by minor-allele frequency: `lambda_by_maf.tsv`.

### Open question (to be agreed)
Whether, and on which of the facts above, the main analysis should include the PCs on top of the GRM. This
section is left without a conclusion on purpose.

### Methods
As `01.fixed_vs_random` (λGC, QQ band, drawing). PC test: `statsmodels` Logit, likelihood-ratio
χ²(10) of SEX + PC1–10 against SEX alone, per-PC Wald P and McFadden pseudo-R². Upstream peak calls of the
`main` (with PCs) and `detect` (without PCs) scans of `analysis/assoc_saige`.
