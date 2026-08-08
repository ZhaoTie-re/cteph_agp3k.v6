# tuning.rv — Methods

Single source of truth for the methodology. The pipeline, the per-cohort run
manifests, and the figures all refer here; nothing is duplicated elsewhere.

---

## The whole pipeline on one screen

Three thresholds are chosen, and two modelling choices decide how the evidence for
them is read. Numbers are `narrow_mainland`; the other cohorts are in
`_comparison/cohort_comparison.tsv`.

| | Decides | Rule | Threshold from | Cost | Evidence |
| --- | --- | --- | --- | --- | --- |
| **§3 Sample QC** | which samples to drop | \|robust-`Z`\| > 5 on the Theil–Sen depth residual | **fixed parameter** `RobustZThreshold` | 3 samples | `sample_qc_rate_residual.png`, `rate_distribution.png` |
| **§4 Variant QC** | which variants to drop | `Δ` > Kneedle knee of `F(Δ)` | **data-driven** (`S`, `wx`, `wy`) | 1,352,762 variants (8.4 %) | `depthdiff_variants_cdf.png` |
| **§6 minAC** | the threshold `κ` | smallest `κ ≥ floor` with `p ≥ α` **and** `0 ∈ 95 % CI` | **data-driven**, falls back to `argmin \|β_group\|` | 53 % of variants | `minac_sweep.png` |
| **§2 Metric** | *how burden is measured* | rate `r = S_minAC / N_callable`, not the raw count | design choice | — | `rate_distribution.png` |
| **§5 Model** | *which statistic the decision is read from* | `β_group`, from the **no-PC** model | design choice, evidenced | — | `decision_axis.png` |

The last two rows remove nothing, but they fix how every reading above is interpreted —
see §5, which is where the two most common misreadings of this pipeline are addressed
(selecting on `β_depth`, and adjusting the tuning model on ancestry).

**What this pipeline does not decide:** whether any association is real. See §8.

---

## Symbols & notation

The same symbols are used in the code, the figures, and this document. Every figure
now prints the definitions of the symbols it plots in its own caption, drawn from
`SYMBOL_DEFS` in `scripts/plot_style.py` — that dict and this table are the two faces
of one definition and must be kept in step.

| Symbol | Meaning **in this analysis** |
| --- | --- |
| `r_i` | per-sample minor-allele burden **rate** = (n_het + 2·n_hom-alt) / N_call, ×10³ — the response variable of every model here |
| `N_call` | non-missing genotypes for a sample (hom-ref + het + hom-alt) |
| `Case_i` | case indicator ∈ {0,1} (1 = Case/PH) |
| `D_i`, `z(D_i)` | sample mean depth; `z(D)` = standardized (z-scored) mean depth |
| **`β_group`** | OLS coefficient on `Case` = the **apparent** Case−Control difference in `r`. Cases are 30× and controls 15× with no platform overlap, so at low `κ` this measures the depth/platform **artifact**, not phenotype. The **decision** axis: `κ` is taken at its non-significant trough |
| **`β_depth`** | OLS coefficient on `z(D)`, per +1 SD depth = residual depth sensitivity of `r` given group (and PCs). A **diagnostic** only — smallest exactly where `β_group` is largest, so it must never drive `κ` |
| `β_k`, `PC_k` | ancestry-PC coefficients / bbj_mainland principal components (the +PC **sensitivity** model; no-PC is primary) |
| `ε_i` | OLS residual; also the Theil–Sen depth residual `r_i − r̂(D_i)` in sample QC |
| robust-`Z_i` | 0.6745(ε_i − median ε)/MAD(ε) — median/MAD sample-outlier score, so extreme samples cannot inflate their own threshold |
| `κ` | minor-allele-count (minAC) threshold being swept — the quantity being tuned |
| `Δ_v` | per-variant depth-stratified missingness gap = \|VMISS_30× − VMISS_15×\| — how differentially a variant drops out with depth |
| knee | Kneedle elbow of the CDF `F(Δ)`; variants above it are excluded (protected candidates exempt) |
| `VMISS_{d,v}` | missing rate of variant *v* within depth stratum *d* (15× or 30×) |
| `α`, CI | a coefficient is **non-significant** only when *p* ≥ `α` **and** 0 ∈ [L, U] of its 95% CI — the same rule that selects the recommended `κ` |
| `λ` | genome-wide inflation factor from the real rv test — the final calibration check |

---

## 1. The problem: a fully confounded design

In CTEPH/AGP3K v6 the sequencing design is **completely confounded** with
phenotype:

| Group | Target depth | Platform(s) |
| --- | --- | --- |
| Cases (PH) | 30x | DNBSeq-T7, NovaSeq, DNBSeq-G400RS |
| Controls (AGP3K) | 15x | HiSeqX |

There is **zero platform overlap** between cases and controls. Any raw
rare-variant burden difference between cases and controls is therefore
inseparable from a depth/platform difference. No statistical adjustment can undo
this: depth and platform are *collinear* with phenotype.

Consequently `tuning.rv` does **not** try to confirm associations. Its job is to:

1. **quantify** the depth/platform inflation of the rare-variant burden test,
2. **mitigate** it layer by layer (sample QC → variant QC → minAC threshold →
   ancestry adjustment),
3. **recommend** a minAC threshold that brings the *adjusted* apparent phenotype
   effect to a statistical null while preserving as much power as possible,
4. run this across three nested cohorts as a **sample-selection sensitivity**
   analysis, and
5. hand off the residual (platform-chemistry) axis to downstream **depth-matched
   down-sampling** — the internal replacement for external replication in a rare
   disease — plus read-level checks and a genome-wide λ scan.

**Guiding principle — one layer, one confounder:**

| Confounder | Handled by | Stage |
| --- | --- | --- |
| variant-level depth-differential missingness | depth-diff VMISS knee | variant QC |
| sample-level burden anomaly | rate-residual robust-Z | sample QC |
| singleton / low-count calling artifact | minAC threshold | sweep |
| ancestry | bbj_mainland PCs — *measured* in the +PC sensitivity model, deliberately **not** adjusted on in the primary one (§5) | sweep regression |
| **platform chemistry (residual)** | **not fixable internally** → down-sampling + stated limit | downstream |

---

## 2. The metric: burden **rate**, not raw count

For each sample define the minor-allele burden **rate**

```
r = S_minAC / N_callable            (reported ×1000, i.e. per 1,000 callable sites)
```

- `S_minAC` = het + 2·hom-alt minor-allele count over the retained variants,
- `N_callable` = het + hom-ref + hom-alt = **non-missing** genotypes.

Using the rate (rather than the raw `S_minAC`) **decouples burden from
differential missingness**: 15x controls miss more genotypes, which deflates the
raw count in a depth-dependent way. Dividing by the callable denominator removes
that first-order artifact so the residual signal reflects genuine burden, not
coverage.

---

## 3. Sample QC (FIRST) — depth-adjusted rate-residual outliers

Sample QC runs **before** variant QC, because variant metrics are computed
*across samples* — anomalous samples must be removed before per-variant missing
rates are trusted.

1. Robustly regress rate on depth: Theil–Sen fit `r ~ D_mean` (median-based;
   insensitive to outliers).
2. Residual `ε = r − (slope·D_mean + intercept)`.
3. Robust-Z of the residual: `Z = 0.6745·(ε − median(ε)) / MAD(ε)`.
4. Flag `|Z| > Z_thr` (default `RobustZThreshold = 5`) → written to a PLINK2
   `--remove` list.

The **Sample-QC figure** (`sample_qc_rate_residual.png`) shows the residual vs
depth coloured three ways (group / target-depth / platform); a residual that
tracks platform/depth rather than phenotype is the confounding, made visible.
The **rate-distribution figure** (`rate_distribution.png`) shows the rate before
vs after removal (skewness / excess-kurtosis quantify the tail that QC trims) and
the same rate split by depth / platform / group — panels that should mirror each
other under this design.

Variant metrics (`VMISS_15X`, `VMISS_30X`, …) are then **recomputed on the
sample-cleaned set** (process `CALC_VARIANT_METRICS`).

---

## 4. Variant QC — depth-differential missingness (Kneedle knee)

For every variant compute the depth-stratified missing rates within the 15x and
30x strata and their absolute gap

```
Δ = | VMISS_30X − VMISS_15X |
```

> **Two different depth variables — do not conflate them.** These strata are built from
> the **nominal `TargetDP`** in the info file (`run_calc_metrics.sh`, step 4), i.e. the
> 15x/30x *label* a sample was sequenced under. The OLS in §5 instead uses the
> **observed `MeanDP`**. Both are correct for their purpose — the strata define who is
> compared with whom, the regression conditions on what depth was actually achieved —
> but they are not the same variable and the achieved depths overlap heavily between
> the two labels.

Build the CDF of Δ and locate its **Kneedle knee** (concave/increasing;
sensitivity `S`, `weight_x`, `weight_y` tunable; knee rounded to 2 decimals).
Variants with `Δ > knee` have depth-dependent missingness and are written to a
PLINK2 `--exclude` list (`variants_depthdiff.exclude.txt`). Figure:
`depthdiff_variants_cdf.png`.

**Protected variants.** VariantIDs listed in `params.ProtectVariants`
(candidate signals; bim format `chr<N>:<pos>:<ref>:<alt>`) are **never** placed
on the exclude list — they are dropped from it even if `Δ > knee`, so they always
survive into the filtered call set. The count rescued this way is logged and
recorded in `depthdiff_stats.tsv` (`n_protected_rescued`).

The sample-`--remove` and variant-`--exclude` lists are applied in two PLINK2
steps: `FILTER_SAMPLES` → `FILTER_VARIANTS` → fully-filtered call set.

---

## 5. The minAC sweep — quantifying and calibrating the residual effect

On the fully-filtered set, for `minAC = κ = 0 … MaxAC`, recompute metrics and fit
**OLS on the rate** per threshold, in **two declared models** (one fit loop;
classical SE by default; `--robust` → HC3):

```
r_i  =  β0  +  β_group·Case_i  +  β_depth·z(D_i)                          (noPC — primary, depth-adjusted)
r_i  =  β0  +  β_group·Case_i  +  β_depth·z(D_i)  +  Σ_k β_k·PC_{k,i}      (pc   — +ancestry sensitivity)
```

### 5.1 Two coefficients, two roles — and why `β_group` is the decision axis

**`β_group` is a DETECTOR, not an effect estimate.** This one idea decides both of the
modelling questions in this section. The design gives cases and controls zero platform
overlap, so the artifact *is* the case/control contrast; `β_group` is the instrument
that reports how much of it survives at each κ. A detector is chosen for **sensitivity
to its target**, not for the unbiasedness of an estimate — which is why the minAC is
read from `β_group`, and read from the model that leaves it intact (§5.4).
- **`β_depth` (depth) is a DIAGNOSTIC only.** It reports the residual depth–burden
  coupling (this replaced the earlier Spearman ρ, and comes from the *same* OLS).
  **Do NOT choose minAC to minimise `β_depth`:** the two coefficients simply do not
  bottom out in the same place, and `β_depth`'s minimum sits well inside the range
  where the apparent case/control effect is already large. Measured
  (`_comparison/decision_axis.tsv`, panel **b** of `decision_axis.png`):

  | cohort | `argmin \|β_group\|` | `argmin \|β_depth\|` | `β_group` there | p | vs. its own minimum |
  | --- | --- | --- | --- | --- | --- |
  | narrow_mainland | κ = 2 | κ = 6 | 0.1262 | 4.5e−24 | **17×** larger |
  | intermediate_mainland | κ = 2 | κ = 8 | 0.1678 | 3.2e−37 | **75×** larger |
  | full_mainland | κ = 3 | κ = 12 | 0.1807 | 6.1e−38 | **107×** larger |

  So selecting on `β_d` would maximise, not minimise, the test bias.

### 5.2 Why depth is a covariate

The two models are:

| Model | Formula | Role |
| --- | --- | --- |
| `noPC` | `rate ~ group + z(D)` | **primary** — depth-adjusted apparent effect (drives the minAC choice) |
| `pc` | `rate ~ group + z(D) + PC` | **sensitivity** — does ancestry change the picture? |

Depth is kept in the primary model for three reasons: (i) it is a *precision
covariate* — observed depth predicts the burden rate (`β_depth` is strongly
significant), so conditioning on it shrinks the residual variance and tightens
`β_group`'s CI; (ii) it is only **weakly correlated with group** (corr ≈ 0.17–0.27,
R² ≈ 0.03–0.07 — cases are 93% 30x but with a 7% 15x overlap, and observed depths
overlap heavily), so including it does **not** steal group's variance (VIF ≈
1.03–1.08) — the textbook "include it" case; (iii) it is the *correct specification* —
omitting a real rate-predictor that is (weakly) group-correlated would bias `β_group`.

**Honest limit — the depth term is not "the depth control."** Observed `z(D)` is only a
weak proxy for the true confounder (platform / target-depth chemistry), which is
**collinear with group and therefore unadjustable**. So it is a *minor
precision/diagnostic* covariate. The real depth/platform confounding is handled by the
metric (rate), the depth-diff variant QC, the minAC threshold, and downstream λ /
down-sampling — not by this coefficient.

### 5.3 Why the PC source is bbj_mainland, and geno-PC is excluded

**geno-PC is deliberately NOT used** for adjustment: in this design the genotype PCs
capture platform, so adjusting on them would regress out the very effect under study.
Ancestry is measured with **bbj_mainland** reference PCs (default 10) — axes defined by
an external reference, so they cannot be driven by this cohort's own batch structure.
geno-PC may be added only as an optional negative control.

### 5.4 Why the primary model is no-PC

Given §5.1 — `β_group` is a detector — a covariate that carries case/control contrast
will absorb part of what the detector is meant to see, by an amount that cannot be
recovered. The `bbj_mainland` PCs do carry some: they span
**R²(case/control) = 0.009–0.028** depending on the cohort. The consequence is measured
in §6 and in `decision_axis.png` panel **c**: at the decision point the PC-adjusted
reading moves by **1.9–2.6× the size of the reading itself**, in cohort-dependent
directions. That is why the tuning model omits them, and why the downstream association
test — which needs an estimator rather than a detector — keeps them.

### 5.5 Empirical finding — `β_group` is V-shaped in minAC, not monotone

In all three cohorts it is significantly inflated at minAC 0–1, dips to a
**non-significant trough at minAC = 2** (the decoupling point), then **climbs steeply
again** for minAC ≥ 4 (β up to ~0.5–0.75, p → 1e−100 or smaller). Raising minAC past
the trough does *not* help — the higher-count variants are exactly where 30x-vs-15x
calling differs most, so they re-introduce (worse) platform inflation. The metric
choice does not change this: raw `S_minAC` gives the same V with the same trough.

> **Poisson was dropped.** `S_minAC` is a large near-normal sum (~10⁴); a Poisson
> GLM was massively over-dispersed and produced meaningless p-values (~1e-300).
> OLS on the rate is the correct, well-calibrated model here.

### 5.6 Figures

`minac_sweep.png` — 2×2 vs minAC: (a) `β_g` full range (the V), (b) `β_g` zoomed on the
trough (the decision region), (c) `β_d` — the residual depth effect (diagnostic),
(d) variants retained (power cost). no-PC and +PC series carry 95% CI ribbons with
significant/non-significant markers; the recommended minAC is marked on every panel.

`_comparison/decision_axis.png` — why the decision is read from `β_group` in the no-PC
model: (a) platform × case/control (zero overlap), (b) the two coefficients' minima,
(c) what PC adjustment does to the reading. Numbers in `decision_axis.tsv`.

---

## 6. minAC recommendation (V-trough, power-preserving)

Because `β_group` is V-shaped, the recommendation is the **trough**, not a
"raise-until-flat" point. `select_minac.py` recommends the **smallest** minAC
(≥ `MinACFloor`, default 1) at which `β_group` is statistically indistinguishable
from zero:

- `p_group ≥ alpha` (default 0.05, non-significant), **and**
- the 95% CI for `β_group` contains 0.

This is the **left edge of the non-significant band** — it decouples the burden
from depth/platform while keeping the most variants (maximal power). If nothing is
non-significant, it falls back to `argmin |β_group|` (the least-confounded
threshold) and flags `calibrated = false`. (The earlier "stable for N following
thresholds" rule was removed: it is wrong for a V, where the band is narrow and
the right arm is significant by design.) Output: `minac_recommendation.tsv`.

### Primary model = no-PC; +PC is a sensitivity check

**`β_group` here is a detector, not an effect estimate.** It is the instrument that
reports how much technical artifact survives at each κ, and a detector is chosen for
**sensitivity to its target**, not for the unbiasedness of an estimate. That settles
the model choice:

1. The design gives cases and controls **zero platform overlap**, so the artifact *is*
   the case/control contrast — anything that differs between the groups is entangled
   with the technical axis.
2. The ancestry PCs carry part of that contrast: the 10 `bbj_mainland` PCs span
   **R²(case/control) = 0.009–0.028** depending on the cohort.
3. So conditioning on them perturbs the detector, at the exact κ where the decision is
   made:

   | cohort | `β_group` no-PC | `β_group` +PC | shift | vs. \|β\| |
   | --- | --- | --- | --- | --- |
   | narrow_mainland | −0.00759 | **+0.00669** | +0.01428 | 1.9× — sign flip |
   | intermediate_mainland | +0.00225 | **−0.00364** | −0.00589 | 2.6× — sign flip |
   | full_mainland | −0.00824 | −0.02531 | −0.01706 | 2.1× |

4. The shift is **larger than the quantity being measured**, and its **direction is
   cohort-dependent** (+0.014 in one, −0.017 in another) — so it is not a bias that
   could be corrected for, and the adjusted β cannot locate the trough.

`select_minac.py --prefer nopc` (the pipeline default) therefore reads the trough from
the no-PC model, and additionally reports the +PC recommendation and whether the two
**agree at the recommended minAC** (`models_agree_at_rec`). Empirically no-PC gives
**minAC = 2 in all three cohorts**; +PC agrees in narrow/intermediate and recommends
κ = 1 in `full_mainland` — that divergence is this same perturbation surfacing in the
recommendation. Evidence: `_comparison/decision_axis.png` (panel **c**) and
`decision_axis.tsv`.

**Why the downstream `rvtest` still includes PCs.** Because it is doing a different
job. Tuning needs a **detector** (maximum sensitivity to the artifact); the association
test needs an **estimator** (ancestry adjusted). Different purpose, different model —
there is no inconsistency between the two.

The genome-wide inflation factor **λ** from the real burden scan (assoc_rvtest
layer, PC-adjusted) is the gold-standard follow-up that confirms the chosen minAC.

> **Note on the candidate (protected) variants.** They are rare (MinAC 3–7) and
> sit just above the minAC = 2 trough, so they survive the recommended threshold
> (MinAC ≥ 2) — but they live in the low-count regime the threshold is designed to
> treat with caution. They must be validated by the targeted route (depth-matched
> down-sampling + read-level review), not read off the genome-wide calibrated scan.

---

## 7. Cross-cohort comparison

The three nested cohorts (`narrow_mainland` ⊂ `intermediate_mainland` ⊂
`full_mainland`) differ only in **sample selection**. Overlaying their sweeps
(`cohort_comparison.png`) tests whether the `β_group` trajectory and the
recommended minAC are **robust to sample selection**. This does **not** break the
platform confounding — cases are 30x and controls 15x in every cohort alike.

---

## 8. Scope and honest limits

`tuning.rv` delivers: a quantified confounding profile, layered mitigation, a
recommended minAC, qualified candidate variants, and explicit limits — **not**
confirmation of association. The residual **platform-chemistry** axis is not
fixable from these data; it is handed to depth-matched **down-sampling**
(internal replication surrogate, depth axis only), read-level inspection, and the
genome-wide λ scan.
