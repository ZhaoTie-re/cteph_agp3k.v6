# assoc_gene — gene-based rare-variant association, CTEPH

Three nested cohorts × two impact strata × two statistics (CMC burden, SKAT-O) under one fixed-effect
logistic model (sex + 10 ancestry PCs), with the multiple-testing denominator fixed before any test, a
robustness tier per gene across cohorts and statistics, a publication summary table and a main figure
for the robust genes.

## What it establishes — and what it does not

* Per (cohort, stratum, method): which genes are **called** (BH q < 0.05 over the mapped genes) and which
  also clear Bonferroni (α / `n_genes_mapped`).
* Per (gene, stratum): a **tier** — Tier 1 called in every cohort and by both statistics somewhere,
  Tier 2 called in every cohort, Tier 3 called in ≥ 2 cohorts.
* For Tier 1/2 genes: variants, carriers, cumulative MAC and MAF per group, SKAT-O ρ and p, CMC OR
  (95 % CI) and p, and the number of samples carrying more than one site of the gene.

It does **not** establish a CTEPH gene. Cases and controls were sequenced on different platforms with
zero overlap; that confound is in every p-value and no covariate removes it. And the cohorts are nested
(they differ almost entirely in controls: 1,747 → 2,049 → 2,630, against 418 → 428 → 438 cases), so
agreement across them shows a signal survives adding controls — a sample-selection sensitivity
analysis, **not replication**. See [docs/METHODS.md](docs/METHODS.md) §1, §3.

## Run

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_gene
./verify.sh                       # pre-flight; must be all green
source activate dsl2
nextflow run assoc_gene.nf        # 397 tasks; SLURM queue gr10478b
./verify.sh                       # post-run sections now fire
cat results/_comparison/robust_genes/summary_table.md
```

`nextflow lint assoc_gene.nf` is the syntax check (`nextflow run --help` launches the pipeline).
Never follow `-preview` with a bare `-resume`.

## Deliverables

| what | where |
|---|---|
| summary table (Tier 1 + 2 genes, per stratum block) | `results/_comparison/robust_genes/summary_table.md` (+ `.tsv` per cohort) |
| main figure: tiers, evidence, OR, carriers | `results/figures/04.robust_genes/robust_genes.png` |
| per-gene detail figures | `results/figures/05.genes/<GENE>.png` (+ `README.md`, `index.tsv`) |
| evidence and tiers | `results/_comparison/robust_genes/{evidence,tiers}.tsv` |
| every other figure family | `results/figures/01.variant_sets`, `02.gene_scan`, `03.calibration` |

`results/README.md`, written by the run, indexes the tree.

## Inputs (contract with `tuning.rv`)

| input | path |
|---|---|
| genotypes | `tuning.rv/results/<cohort>/02.callset_filter/filtered.{bed,bim,fam}` |
| minAC | `tuning.rv/results/<cohort>/05.qc_collect/minac_recommendation.tsv` |
| phenotype, covariates | `model_inputs/<cohort>/{phenotype/pheno.tsv, covariates/bbj_mainland_pc.sex.tsv}` |
| annotation | JHRP snpEff index (`params.SnpEffIndex`) |

Variant IDs must be `chr<N>:<POS>:<REF>:<ALT>`; the pipeline asserts it.

## Key parameters (`assoc_gene.nf`)

| param | default | effect |
|---|---|---|
| `Cohorts` | narrow, intermediate, full `_mainland` | tier rule requires a call in every listed cohort |
| `Strata` | `low_moderate_high`, `moderate_high` | snpEff impact classes per stratum; widest first |
| `MinNumVar` | 2 | genes with fewer mapped variants are dropped **before** the denominator |
| `Alpha` | 0.05 | BH q and Bonferroni α |
| `NPcs` / `CovarName` | 10 / `sex,pc1_avg…pc10_avg` | the model's covariates |
| `ImputeMethod` | `mean` | rvtest `--impute` (mean / hwe / drop only) |
| `ReportCohort` | `full_mainland` | whose numbers fill `summary_table.md` and panels (b)(c) of the main figure |
| `MaxLabelGenes` | 24 | cap on the names in one Manhattan. Only a called gene is ever named; a scan that calls nothing names nothing |
| `MaxRobustGenesPlotted` | 40 | rows in `robust_genes.png` (the figure grows with the count) |
| `StrictIdentity` | true | `SCAN` aborts on a set whose rvtest variant count is unexplained |

Resources are set only through `clusterOptions --rsc` in `nextflow.config`; scripts read
`task.ext.threads`.

## Pipeline

```
PREP_PHENO_COV ─┐
SPLIT_BIM ──────┼─ MAP_CHUNK ×22 ─ MERGE_MAP  [B1: denominator fixed] ─ PLOT_ANNOTATION ─┐
MAKE_MINAC_KEEP ┘        │                                                                CATALOGUE
                 EMIT_TEST_INPUTS  (set files + per-gene counts from ONE dict; bed decode cross-checked)
                         │
                 RVTEST ×22 ×2 methods  (5 post-flight gates on rvtest's own log)
                         │
                 MERGE_RVTEST [B2] ─ SCAN [B3: map identity + BH] ─ SIGNIFICANCE
                                            └ SCAN_SCALE ─ PLOT_SCAN ─ CATALOGUE
                 COLLECT_ALL [B4] ─ ROBUSTNESS ─ PLOT_ROBUST_GENES
                                  │             └ GENE_MODELS ─ PLOT_GENE_DETAIL
                                  └ PLOT_GRID
```

Four barriers: **B1** the denominator is decided from the map alone; **B2** the 22 parts of one
(cohort, stratum, method) become one table; **B3** the tested sets are proven to be the mapped sets and
both decision rules are applied over the fixed family; **B4** every cell exists before anything is
compared across cells. See [docs/METHODS.md](docs/METHODS.md).

## Outputs

`results/<cohort>/{00.prep, 01.map, 02.test_inputs, 03.assoc, 04.scan, 05.signals}`,
`results/_comparison/{tables, robust_genes}` and `results/figures/{01.variant_sets, 02.gene_scan,
03.calibration, 04.robust_genes, 05.genes}`. Column dictionaries in [docs/OUTPUTS.md](docs/OUTPUTS.md);
the figure system in [docs/FIGURES.md](docs/FIGURES.md).

## Verification

`./verify.sh` — 21 sections, explained in [docs/VERIFICATION.md](docs/VERIFICATION.md): syntax and site
rules, flag/schema consistency between the `.nf` and the scripts, a live rvtest smoke test that pins the
output names and the cmcWald block height, a vocabulary check, a fixture for the tier rule, and after a
run the DAG cardinality, the span-overlap predictor, the count invariants, the effect-size arithmetic,
the summary-table formats, evidence completeness, the variant-ID invariant, the one-README-per-family
rule and the figure geometry.

## Documents

* [docs/METHODS.md](docs/METHODS.md) — samples, map, model, statistics, significance, counts, tiers.
* [docs/OUTPUTS.md](docs/OUTPUTS.md) — every file and column.
* [docs/FIGURES.md](docs/FIGURES.md) — the figure system: families, sizes, vocabulary, omissions.
* [docs/STUDY_NOTES.md](docs/STUDY_NOTES.md) — the measured facts the design rests on, and limitations.
* [docs/VERIFICATION.md](docs/VERIFICATION.md) — what each check proves.

`scripts/` belongs to this component (`vocab.py` holds its figure vocabulary); `../_shared/scripts/`
(`plot_style`, `figure_doc`) is a shared toolkit and is never edited here.
