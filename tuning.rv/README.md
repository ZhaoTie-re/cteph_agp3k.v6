# tuning.rv

Rare-variant depth/platform-confounding QC and a null-calibrated minor-allele-count
(minAC) recommendation for CTEPH/AGP3K v6, run across three nested cohorts with a
cross-cohort comparison.

## The problem, in one line

The design is **fully confounded**: cases (PH) are 30x (DNBSeq/NovaSeq), controls
(AGP3K) are 15x (HiSeqX), with **zero platform overlap**. A raw rare-variant burden
difference cannot be separated from a depth/platform difference. So this pipeline
does **not** confirm associations — it **quantifies** the confounding, **mitigates**
it one layer at a time, and **recommends** a minAC that brings the *adjusted*
apparent phenotype effect to a statistical null while preserving power.

→ Full methodology: **[docs/METHODS.md](docs/METHODS.md)** ·
Output tree & data dictionary: **[docs/OUTPUTS.md](docs/OUTPUTS.md)**

## What gets decided, and on what evidence

Three thresholds, plus two choices that fix how the evidence for them is read. Costs
are `narrow_mainland`. Full argument in [docs/METHODS.md](docs/METHODS.md).

| | Decides | Rule | Threshold from | Cost | Evidence |
| --- | --- | --- | --- | --- | --- |
| Sample QC | which samples to drop | \|robust-`Z`\| > 5 on the depth residual | **fixed parameter** | 3 samples | `sample_qc_rate_residual.png` |
| Variant QC | which variants to drop | `Δ` > Kneedle knee | **data-driven** | 1,352,762 (8.4 %) | `depthdiff_variants_cdf.png` |
| minAC | the threshold `κ` | smallest `κ` with `p ≥ α` **and** `0 ∈ CI` | **data-driven** + fallback | 53 % of variants | `minac_sweep.png` |
| Metric | *how* burden is measured | rate, not raw count | design | — | `rate_distribution.png` |
| Model | *which statistic* the decision reads | `β_group`, from the **no-PC** model | design, evidenced | — | `decision_axis.png` |

Note the asymmetry: the sample-QC cutoff is a fixed parameter, the other two are
data-driven. **Not decided here:** whether any association is real.

## Method at a glance (one layer, one confounder)

```
CALC_METRICS_BASE (raw, minAC 0; permanent storeDir cache)
  └─ DETECT_OUTLIER_SAMPLES   rate = S_minAC/N_callable; Theil–Sen r~D_mean; robust-Z   → remove
       ├─ PLOT_RATE_DISTRIBUTION   rate before/after QC + by depth/platform/group
       └─ FILTER_SAMPLES  (--remove)   → sample-cleaned set
            └─ CALC_VARIANT_METRICS   recompute VMISS_15x/30x on the cleaned set
                 └─ DETECT_DEPTHDIFF_VARIANTS   Kneedle knee of |ΔVMISS| CDF   → exclude
                      └─ FILTER_VARIANTS  (--exclude)   → fully-filtered set
                           └─ CALC_METRICS (minAC 0..Max)
                                └─ QC_SUMMARY   OLS  r ~ group + Z(D_mean) [+bbj_mainland PC]
                                     └─ MERGE → SELECT_MINAC (+ PLOT_SWEEP) → WRITE_RUN_MANIFEST
...then COMPARE_COHORTS across cohorts.
```

- **Metric = rate** (`S_minAC / N_callable`, ×1000) — removes depth-driven
  differential missingness. **Model = OLS** (Poisson dropped — overdispersed).
- **Sample QC before variant QC** (variant metrics recomputed on the cleaned set).
- **Ancestry** via **bbj_mainland** PCs (default 10); **geno-PC is NOT used** (it
  captures platform here).
- **Recommended minAC** = the **V-trough**: `β_group` is V-shaped in minAC
  (inflated low, non-significant trough at minAC≈2, climbs again); we pick the
  smallest minAC where `β_group` is non-significant (0 in CI). **Primary model =
  no-PC** (`rate ~ group + depth`; isolates the technical artifact); **+PC is a
  sensitivity check** — PCs are for the downstream rvtest, not for choosing minAC.
  Depth is a *precision/diagnostic* covariate (weakly group-correlated), not the
  depth "control"; see [docs/METHODS.md](docs/METHODS.md) §5.
- **Residual platform chemistry is not fixable internally** → downstream
  depth-matched down-sampling + read-level checks + genome-wide λ.

## Run it

```bash
source activate dsl2          # provides nextflow 26.x
cd tuning.rv
nextflow run tuning.rv.nf                       # all three cohorts + comparison
nextflow run tuning.rv.nf --PrecheckOnly true   # stop after detection (fast tuning loop)
nextflow run tuning.rv.nf --Cohorts narrow_mainland,full_mainland
nextflow run tuning.rv.nf -preview              # build the DAG without executing
```

Analysis scripts run under the `cteph_geno_pro` conda env (activated inside each
process). Base metrics are cached in a permanent `storeDir`
(`00.precheck/metrics/`) — recomputed only if the raw call set changes, surviving
`work/` wipes with no `-resume` needed.

### Key parameters (`params.*`, override on the CLI)

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Cohorts` | all three | cohorts to tune (comma list on CLI) |
| `MaxAC` | 20 | top of the minAC sweep |
| `RobustZThreshold` | 5 | sample rate-residual outlier cutoff |
| `KneeS`, `KneeWeightX`, `KneeWeightY` | 1, 1, 2 | Kneedle knee for variant QC |
| `PcTemplate`, `PcLabel`, `NPcs` | bbj_mainland, 10 | ancestry PCs for the +PC model |
| `Robust` | false | OLS SE: classical (default) or HC3 (`--Robust true`) |
| `Alpha`, `MinACFloor` | 0.05, 1 | V-trough non-significance level + floor |
| `SelPrefer` | `nopc` | primary model for the recommendation (`nopc`\|`pc`; pc = sensitivity) |
| `ProtectVariants` | 3 STBD1 IDs | VariantIDs that variant QC must **never** exclude |
| `PrecheckOnly` | false | stop after sample+variant detection |

**Protected variants.** `params.ProtectVariants` lists candidate signals
(`chr<N>:<pos>:<ref>:<alt>`) that are dropped from the depth-diff `--exclude` list
even if `|ΔVMISS| > knee`, so they always survive into the filtered call set.
Currently the three STBD1-region CTEPH candidates. Add new candidates here.

## Reading the figures

Every figure is a **600-dpi PNG**, self-contained and overlap-free, styled after a
print journal. Bold panel letters (**a**, **b**, …) sit at each panel's top-left (in
the title slot, never over data); colour keys sit in an empty panel corner or above
the panels. A thin rule separates the panels from a **Nature/Science-style caption**
— a **bold declarative finding**, then bold-letter run-in panel descriptions (**a**,
…; **b**, …), a sentence of interpretation, and a bold **Data** / **Model** footer
(professional notation). All four are gathered per cohort under `results/<cohort>/figures/`.

- **`sample_qc_rate_residual.png`** — each point is a sample; y = burden-rate
  residual after removing the depth trend. Rings = removed outliers. If the
  residual tracks platform/depth rather than phenotype, that's the confounding.
- **`rate_distribution.png`** — (a) QC trims the heavy tail (skew/kurtosis fall);
  (b–d) the same rate by depth / platform / group — under this design they mirror
  each other, so the group contrast **is** the depth/platform contrast.
- **`depthdiff_variants_cdf.png`** — the Δ = |VMISS₃₀ₓ − VMISS₁₅ₓ| CDF and its
  Kneedle knee; variants above the knee are excluded.
- **`minac_sweep.png`** (headline, 2×2) — (a) `β_g` full range (V-shaped),
  (b) `β_g` trough zoom (the decision region), (c) `β_d` = residual depth effect
  (**diagnostic only** — do *not* pick minAC to minimise it), (d) power cost.
  `β_g` is the decision axis; the recommended minAC (its non-significant trough)
  is marked on every panel. no-PC = primary, +PC = ancestry sensitivity.
- **`cohort_comparison.png`** (2×3) — rows are the two models, columns are `β_g` /
  `β_d` / (power cost, decision summary). Agreement across cohorts means the minAC
  choice is robust to sample selection (not to platform). Panel **f** is the payoff:
  `β_g` ± 95 % CI at each cohort's own recommended minAC, under both models.
- **`_comparison/decision_axis.png`** — *why the decision is read from `β_g` in the
  no-PC model*: (a) platform × case/control, zero overlap; (b) `|β_g|` and `|β_d|`
  bottom out at different minAC, and `β_d`'s minimum sits where the apparent effect is
  17–107× larger; (c) PC adjustment moves the reading by 1.9–2.6× its own size, in
  cohort-dependent directions. Numbers in `decision_axis.tsv`.

## Scope & limits

Quantify + mitigate + recommend minAC + qualified candidates + explicit limits —
**not** confirmation. The platform-chemistry residual is handed to depth-matched
down-sampling (the internal replication surrogate in a rare disease), read-level
inspection, and a genome-wide λ scan.

## Layout

```
tuning.rv/
├── README.md                 # this file
├── docs/METHODS.md           # methodology (single source of truth)
├── docs/OUTPUTS.md           # output tree + data dictionary
├── tuning.rv.nf              # the DSL2 pipeline
├── nextflow.config           # resources + run reports
├── scripts/                  # analysis + figure scripts (plot_style.py = shared aesthetic)
└── results/                  # per-cohort outputs, _comparison/, _run_info/
```
