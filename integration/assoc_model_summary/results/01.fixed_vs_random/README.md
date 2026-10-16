# 1 · Fixed vs random model

## Slide summary

**Take-home: a fixed-effect logistic model remains inflated even after adjusting for 10 PCs and removing close relatives; the logistic mixed model with a full GRM is well calibrated. The main analysis therefore uses the random (mixed) model.**

- Fixed model (PLINK2; SEX + 10 PCs; relatives removed): λGC = 1.120 / 1.121 / 1.132 (narrow / intermediate / full mainland).
- Random model (SAIGE; SEX + 10 PCs + full GRM; relatives retained): λGC = 1.039 / 1.038 / 1.021.
- The excess is genome-wide: in every allele-frequency bin (MAF 0.01–0.50) the fixed model has λGC 1.03–1.15, the random model 0.95–1.06.
- Variants at P < 10⁻⁵ (≈ 51 expected by chance): fixed 130 / 216 / 204; random 87 / 160 / 111.
- Genome-wide loci: fixed 0 / 1 / 2, random 0 / 1 / 0; genome-wide only under the fixed model: *MIR3681HG* in full mainland (P 4.3 × 10⁻⁸ → 1.0 × 10⁻⁷); *FTO* in full mainland (P 2.9 × 10⁻⁸ → 1.8 × 10⁻⁷).

## Figures

- **`fixed_vs_random.png`** — one row per cohort: Manhattan of the fixed model (left) and of the random model (middle) on one genomic axis and one −log₁₀P scale, and the QQ plot of both (right, λGC in the model colour).
- **`calibration.png`** — (a) λGC per cohort; (b) λGC by minor-allele frequency; (c) variants at P < 10⁻⁵ against the number expected by chance.

## Detailed description

### Question
Does a fixed-effect logistic regression control population structure and relatedness in the three nested
cohorts, or is a mixed model with a genetic relationship matrix (GRM) needed?

### Design
|  | Fixed model | Random model |
|---|---|---|
| Software | PLINK 2 logistic regression (`--glm`, no Firth) | SAIGE logistic mixed model |
| Covariates | SEX + PC1–10 (fixed effects) | SEX + PC1–10 (fixed) + polygenic random effect, full GRM |
| Samples | relatives removed (PI_HAT > 0.2, vertex cover) | relatives retained; the GRM models them |
| Cases / controls | 419/1,749 · 429/2,051 · 439/2,632 | 419/1,776 · 429/2,079 · 439/2,662 |
| Variants | 5,113,526 (narrow) | 5,113,626 (narrow) |
| P value | Wald test | score test, normal approximation |

The two scans share the variant set (the same IDs in the same order) and the covariates (SEX + 10 PCs). They
differ in the polygenic random effect and in the 27 / 28 / 30 related controls that the fixed model removes
and the random model keeps (narrow / intermediate / full mainland).

### How to read the figures
- **`fixed_vs_random.png`** — rows a–c are the cohorts. Every Manhattan shares one genomic axis and one
  −log₁₀P limit, so heights are comparable between models and cohorts. Red line: 5 × 10⁻⁸; grey dashed line:
  10⁻⁵. Dotted guides mark loci genome-wide in any panel (named above the top row); a diamond marks a
  panel's own genome-wide lead. In the QQ panel, both models are drawn against the uniform null, with the
  pointwise 95% band (grey); a curve that leaves the diagonal from the bottom, not only in the tail, is
  inflation rather than signal.
- **`calibration.png`** — (a) λGC per cohort (red line = 1). (b) λGC within minor-allele
  frequency bins: inflation that persists in every bin is genome-wide; a few true loci cannot raise the
  median of millions of variants. (c) Variants at P < 10⁻⁵ against the number expected by chance
  (dashed line, n × 10⁻⁵).

### Results
| Cohort | λGC fixed | λGC random | P<10⁻⁵ fixed | P<10⁻⁵ random | expected | GW loci fixed | GW loci random | SG loci fixed / random |
|---|---|---|---|---|---|---|---|---|
| Narrow mainland | 1.120 | 1.039 | 130 | 87 | 51 | 0 | 0 | 40 / 29 |
| Intermediate mainland | 1.121 | 1.038 | 216 | 160 | 51 | 1 (FTO) | 1 (FTO) | 37 / 27 |
| Full mainland | 1.132 | 1.021 | 204 | 111 | 51 | 2 (MIR3681HG, FTO) | 0 | 51 / 32 |

λGC by minor-allele frequency:

| MAF | variants (random, full mainland) | fixed narrow | fixed intermediate | fixed full | random narrow | random intermediate | random full |
|---|---|---|---|---|---|---|---|
| <0.01 | 42,690 | 2.758 | 2.844 | 2.805 | 2.007 | 1.619 | 1.458 |
| 0.01–0.02 | 629,103 | 1.035 | 1.033 | 1.049 | 0.958 | 0.964 | 0.954 |
| 0.02–0.05 | 802,916 | 1.103 | 1.103 | 1.120 | 1.016 | 1.010 | 1.003 |
| 0.05–0.10 | 715,906 | 1.097 | 1.094 | 1.105 | 1.012 | 1.011 | 0.991 |
| 0.10–0.20 | 967,707 | 1.131 | 1.139 | 1.146 | 1.053 | 1.059 | 1.038 |
| 0.20–0.50 | 1,956,794 | 1.139 | 1.139 | 1.152 | 1.063 | 1.062 | 1.046 |

Genome-wide leads of the fixed model, with the smallest random-model P within ±250 kb:

| Cohort | Gene | Variant | rsID | P fixed (lead) | P random (locus) |
|---|---|---|---|---|---|
| Intermediate mainland | FTO | chr16:53887925:T:C | rs16952623 | 1.2 × 10⁻⁸ | 3.8 × 10⁻⁸ |
| Full mainland | MIR3681HG | chr2:11891017:G:A | rs66898022 | 4.3 × 10⁻⁸ | 1.0 × 10⁻⁷ |
| Full mainland | FTO | chr16:53887925:T:C | rs16952623 | 2.9 × 10⁻⁸ | 1.8 × 10⁻⁷ |

### Interpretation
- Removing close relatives and adjusting for 10 PCs does not calibrate the fixed model: λGC
  1.12–1.13 in all three cohorts, raised in every frequency bin, and 130–216 variants
  at P < 10⁻⁵ where about 51 are expected. A genome-wide excess of this kind is the signature of
  structure that the fixed covariates do not capture: distant relatedness below the PI_HAT cut-off and
  fine-scale ancestry beyond 10 PCs.
- The random model, on the same variants with the same covariates, absorbs that genome-wide sharing through
  the GRM: λGC 1.021–1.039, and 0.95–1.06 in every frequency bin from MAF 0.01. It also keeps the
  related samples instead of discarding them.
- Loci: the inflated fixed model calls more suggestive signals (40 / 37 / 51
  against 29 / 27 / 32) and genome-wide calls that the random model
  does not reproduce (*MIR3681HG* in full mainland (P 4.3 × 10⁻⁸ → 1.0 × 10⁻⁷); *FTO* in full mainland (P 2.9 × 10⁻⁸ → 1.8 × 10⁻⁷)). *FTO* stays the strongest locus under both models.
- Hence the main analysis is the random (mixed) model.

### Limits
- The sample sets differ by the removed relatives (27 / 28 / 30 controls). Those samples are absent from the
  fixed model, so they cannot be the cause of its inflation.
- λGC mixes polygenic signal and confounding and does not separate them. The two models are compared on
  the same phenotype and variants, so polygenic signal cannot explain the difference between them.
- Sequencing platform is confounded with case status (cases and controls were sequenced on different
  platforms). Part of a genome-wide excess can be technical; the mixed model absorbs genome-wide components
  of either origin, but single signals still need the locus-level checks.
- Variants with overall MAF < 0.01 (41,060 / 47,560 / 42,690 in the random-model scans, kept
  because the upstream MAF filter is applied to the controls) are inflated under both models (λGC
  2.76 / 2.84 / 2.80 fixed, 2.01 / 1.62 / 1.46 random). They are less than 1% of
  the variants, are not plotted in calibration panel b, and are listed for follow-up.
- P values: PLINK 2 Wald test and SAIGE score test with the normal approximation (no saddle-point
  correction); both are asymptotic.

### Methods
- λGC = median(χ²₁ of P) / 0.4549 over every usable variant: ERRCODE '.', 0 < P ≤ 1, finite OR
  and SE. λGC and the counts are recomputed here and equal the upstream `scan_qc.tsv` (asserted).
- QQ band: pointwise 2.5–97.5% quantiles of Beta(k, n − k + 1) for the k-th smallest of n uniform P. The
  QQ draws every one of the 20,000 smallest P and a log-spaced sample of the rest.
- Manhattan: every variant at P ≤ 10⁻³ is drawn; below that, overlapping points (closer than ~2 px at
  600 dpi) are drawn once, which leaves the image unchanged.
- Sources: `analysis/assoc_plink2` (additive scans) and `analysis/assoc_saige` (`main` scans), with their
  peak calls (±250 kb).
