# assoc_hla

HLA **allele** and **amino-acid residue** association over three nested cohorts, under both a
fixed-effects model and a full-GRM mixed model.

→ Method and rationale: **[docs/METHODS.md](docs/METHODS.md)** ·
Output tree and data dictionary: **[docs/OUTPUTS.md](docs/OUTPUTS.md)** ·
This study's numbers: **[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md)** ·
Shared visual grammar: **[../_shared/docs/FIGURES.md](../_shared/docs/FIGURES.md)**

## What this component does *not* establish

**Read [`../../hla.typing/docs/OPEN_QUESTIONS.md`](../../hla.typing/docs/OPEN_QUESTIONS.md)
before using the output.** Two things there change how a result should be read:

- **§3** — 4.6 % of typed chromosomes carry an allele that a 61,424-person Japanese reference
  panel has never observed, and it is **differential**: controls 4.74 % against cases 3.57 %,
  a gap that survives matching on measured depth (Fisher *P* = 9.6 × 10⁻¹⁰). It tracks
  platform-and-cohort, which is perfectly confounded with phenotype. A common allele depleted
  more in controls reads as enrichment in cases — **a false risk association from the typing
  alone**. Every allele tested here carries `in_reference`, and any hit with `in_reference = 0`
  is flagged on the figure and in the tables.
- **§1** — DRB3, DRB4 and DRB5 encode hemizygotes as homozygotes, so their dosages are wrong.
  They are **excluded from the primary analysis** here for that and two further reasons.

Neither is a reason to distrust the *allele calls*. Both are reasons to read an association at
these loci as provisional until the typing question of `OPEN_QUESTIONS` §2 is settled.

## The design

| | |
|---|---|
| cohorts | `narrow_mainland` 2,195 · `intermediate_mainland` 2,508 · `full_mainland` 3,101 — **100 % HLA-typed**, so nothing is lost to the join |
| cases / controls | 419/1,776 · 429/2,079 · 439/2,662 |
| models | `fixed` — plink2 `--glm`, relatives removed · `random` — SAIGE full GRM, relatives kept |
| covariates | `SEX` + 10 `bbj_mainland` PCs, both models |
| marker classes | `allele` (2-field, e.g. `HLA_A*24:02`) · `residue` (e.g. `AA_DQB1_57_D`) |
| tested | MAC ≥ 20, call rate ≥ 0.95 — about 150 alleles and 970 residues per cohort |
| significance | Li & Ji effective tests (primary), Bonferroni, and 5 × 10⁻⁸, all reported |

## The four layers

1. **Single allele** — logistic on 0/1/2 dosage.
2. **Single residue** — the same model, residue dosage.
3. **Amino-acid omnibus** — at each position, an *m*−1 df likelihood-ratio test of all its
   residues jointly. This is the layer that makes the analysis fine mapping rather than a
   marker scan.
4. **Forward stepwise conditioning** — how many *independent* signals there are.

## The steps

| # | process | writes |
|---|---|---|
| 1 | `GENE_COORDS` | `00.prep/hla_gene_coords.tsv` |
| 2 | `BUILD_MARKERS` | `<cohort>/00.prep/marker_map.tsv` + the chr6 bfile |
| 3 | `MARKER_QC` | `<cohort>/00.prep/marker_qc.tsv` |
| 4 | `PREP_PHENO_COV` | `<cohort>/00.prep/pheno_covar.tsv` |
| 5 | `PRUNE_MARKERS` | `<cohort>/00.prep/grm.prune.in` |
| 6 | `FIT_NULL` | `<cohort>/01.assoc/random/null.rda` |
| 7 | `ASSOC_FIXED` | (work dir) |
| 8 | `ASSOC_RANDOM` | (work dir) |
| 9 | `TO_SUMSTATS` | `<cohort>/01.assoc/<model>/<class>/sumstats.tsv` |
| 10 | `EXPORT_TESTED` | (work dir) |
| 11 | `OMNIBUS` | `<cohort>/02.omnibus/omnibus.tsv` |
| 12 | `SIGNIFICANCE` | `<cohort>/03.signals/significance.<model>.<class>.tsv` |
| 12b | `OMNIBUS_SIGNIFICANCE` | `<cohort>/03.signals/significance.omnibus.tsv` |
| 13 | `CONDITIONAL` | `<cohort>/03.signals/cond.*` |
| 14–16 | `PLOT_HLA_SCAN` · `PLOT_OMNIBUS` · `PLOT_FREQ_QC` | `<cohort>/figures/` |
| 17–18 | `COMPARE_MODELS` · `COMPARE_COHORTS` | `_comparison/` |
| 19 | `WRITE_RUN_MANIFEST` | `_run_info/run_manifest.json` |

## Running it

```bash
source activate dsl2
cd analysis/assoc_hla
nextflow run assoc_hla.nf
```

Wall clock is dominated by the three `FIT_NULL` fits (~1 h each, run in parallel); everything
else is minutes, because this scans about 1,100 markers rather than 5.1 M.

Useful overrides, none of which re-fit the null:

```bash
--MinMac 0 --MinMaf 0.01     # the MAF >= 1 % sensitivity arm
--DropGenes ''               # include DRB3/4/5 (see docs/METHODS.md before you do)
--PValue no_spa              # SAIGE's normal-approximation p-value instead of SPA
-resume                      # after an interruption; work/ must still exist
```

## Checks worth running after a run

```bash
# nothing failed
awk -F'\t' 'NR>1 && $4!="COMPLETED" && $4!="CACHED"' results/_run_info/trace.txt

# the effect allele is P everywhere — an inversion would flip every odds ratio
awk -F'\t' 'NR>1 && $6!="P"' results/*/01.assoc/*/*/sumstats.tsv | head

# the upstream components were not touched
md5sum -c /path/to/baseline.md5
```
