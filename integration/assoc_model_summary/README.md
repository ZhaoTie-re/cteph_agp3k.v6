# `assoc_model_summary` — the GWAS model choice and the main results, as slide figures and text

An integrative component (`cteph_agp3k.v6/integration/`). It reads two finished association components and
answers three questions with figures and a Markdown text per question (a brief **slide summary** and a
**detailed description**):

| # | Question | Compared scans |
|---|---|---|
| 1 | Is a fixed-effect model enough, or is a mixed (random) model needed? | `assoc_plink2` additive (SEX + 10 PCs, relatives removed) vs `assoc_saige` `main` (SEX + 10 PCs + full GRM, relatives retained) |
| 2 | What do the PCs add on top of the full GRM? *(facts only; the conclusion is to be agreed)* | `assoc_saige` `detect` (SEX + full GRM) vs `main` |
| 3 | What does the main analysis show in the three nested cohorts? | `assoc_saige` `main` |

No association is re-run. The scans are reduced (positions, −log₁₀P, MAF), their λ<sub>GC</sub> and tail counts
are recomputed and asserted equal to the upstream `scan_qc.tsv`, and everything is re-plotted in the style of
`analysis/assoc_saige/results/_comparison/figures/cohort_manhattan.png` (one genomic axis, one −log₁₀P limit,
the same thresholds and QQ band), sized for a 16:9 slide and without a caption block.

## Run

```bash
source activate dsl2
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/assoc_model_summary
./verify.sh                                  # pre-flight
nextflow run assoc_model_summary.nf          # 9 + 4 SLURM jobs and 2 local tasks, ~10 min
./verify.sh                                  # post-run
```

Resume with `nextflow run assoc_model_summary.nf -resume <session>`. Editing a figure or text script re-runs
only the figure tasks; the scan reductions stay cached. `nextflow lint assoc_model_summary.nf nextflow.config`
is the syntax check.

## Deliverables

| path | content |
|---|---|
| `results/README.md` | index of the three questions with the three slide summaries, ready to copy |
| `results/01.fixed_vs_random/` | `fixed_vs_random.png`, `calibration.png`, `README.md`, tables |
| `results/02.pc_effect/` | `grm_pc_vs_grm.png`, `README.md` (facts, conclusion pending), evidence tables |
| `results/03.main_results/` | `main_manhattan.png`, `main_summary.png`, `README.md`, `main_loci.tsv` |
| `results/00.scans/` | per cohort × model: calibration statistics (`stats.json`) and the variants at P < 10⁻³ |

Every file and column: [docs/OUTPUTS.md](docs/OUTPUTS.md).

## Inputs

| input | source |
|---|---|
| fixed-model scans | `analysis/assoc_plink2/results/<cohort>/01.assoc/additive/*.glm.logistic` (symlinks into its `work/`: keep that `work/`) |
| mixed-model scans | `analysis/assoc_saige/results/<cohort>/01.assoc/{main,detect}/sumstats.tsv` |
| scan QC, peaks | each component's `<cohort>/02.scan/scan_qc.tsv` and `<cohort>/03.peaks/model_peaks_annotation.tsv` |
| case status, PCs | `analysis/assoc_saige/results/<cohort>/00.prep/pheno_covar.tsv` |
| cross-cohort loci | `analysis/assoc_saige/results/_comparison/tables/{lead_crosscohort,lead_bymodel}.tsv` |
| style | `analysis/_shared/scripts/plot_style.py` |

The upstream guard reads both `run_manifest.json` files and stops before submitting anything unless:
- the components and models are the expected ones;
- the covariates are SEX + 10 PCs (with and without PCs on the SAIGE side);
- the cohorts are present upstream.

## Key parameters

| parameter | default | meaning |
|---|---|---|
| `saige_dir`, `plink2_dir` | the two `analysis/*/results` | upstream runs |
| `Cohorts` | `narrow_mainland,intermediate_mainland,full_mainland` | upstream order |
| `CohortLabels` | `Narrow mainland,Intermediate mainland,Full mainland` | display names |
| `CohortColors` | `0571B0,008837,B2182B` | the project's slide colours (main summary) |
| `FixedModel`, `GrmPcModel`, `GrmOnlyModel` | `additive`, `main`, `detect` | upstream model names |
| `NPcs` | `10` | PCs in the compared models and in the PC / case-status test |
| `python` | `cteph_geno_pro` interpreter | absolute path: `source activate` cannot be relied on here |

## Pipeline

```
upstream guard (both manifests)
REDUCE_SCAN × 9 (cohort × fixed / grm_pc / grm_only) ─┬─► FIG_FIXED_RANDOM ─┐
PC_ASSOC ─────────────────────────────────────────────┼─► FIG_PC_EFFECT ────┼─► INDEX ─► WRITE_RUN_MANIFEST
                                                      └─► FIG_MAIN ─────────┘
```

## Documents

[docs/METHODS.md](docs/METHODS.md) — statistics, figure rules and text generation ·
[docs/OUTPUTS.md](docs/OUTPUTS.md) — files and columns ·
[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md) — the numbers of the current run.
