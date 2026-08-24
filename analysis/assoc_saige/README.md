# `assoc_saige` — mixed-model association scan with two-tier peak calling

The random-effect counterpart of [`assoc_plink2`](../assoc_plink2). Same design, same peak rule, same
tables and the same figures — a **full genetic relationship matrix** replaces the fixed set of
principal components as the way population structure and relatedness are absorbed.

Two consequences follow, and they are the whole difference:

- **The random-model genotypes are used**, not the fixed-model ones. Relateds are *retained*, because
  the GRM is what accounts for them. The sample sets are therefore larger than the fixed-model
  analysis of the same cohorts.
- **Two covariate models are fitted per cohort**, and only one of them is a result:

| model | covariates | GRM | role |
|---|---|---|---|
| `main` | `SEX` + `params.NPcs` PCs | full | **the reported analysis**; defines every peak |
| `detect` | `SEX` only | full | a **calibration probe** — how much the GRM absorbs on its own |

`detect` exists to answer one question: with a GRM already in the model, how much are the principal
components still doing? It is reported — its own scan figure, its own rows in `model_peaks.tsv`, and
the model-comparison figure — and **nothing fans out from it**, the same treatment `assoc_plink2`
gives its non-primary genetic models.

Peaks come in two tiers, treated differently on purpose:

| tier | threshold | gets |
|---|---|---|
| `genome_wide` | `params.PGenomeWide` | annotation · cross-cohort table · conditional analysis · SuSiE · three-source regional plot |
| `suggestive` | `params.PSuggestive` | annotation · cross-cohort table · scan figure; the `params.MaxSuggestive` **strongest** also get the full per-peak follow-up |

Everything from the summary statistics onward is the **shared** implementation in
[`../_shared/scripts`](../_shared/scripts), because it is about association *results*, not about which
engine produced them.

---

## Run

```bash
source activate dsl2
nextflow run assoc_saige.nf -resume
```

Everything environment- or project-specific lives in one **`SITE CONFIGURATION`** block at the top of
`assoc_saige.nf`. Override on the command line (`--key value`) or with `-params-file cfg.yaml`.

> **Cost note.** `FIT_NULL` and `RUN_ASSOC` are by far the most expensive processes — a full GRM is
> hours per cohort per model. Their task hashes are computed from their script blocks, so **any** edit
> there re-runs the scan, *including adding a comment*: Nextflow hashes that text and `//` is not a
> shell comment. Everything presentational is consumed downstream of `TO_SUMSTATS` and is free to
> change. After a resume, confirm the log lists no `FIT_NULL` or `RUN_ASSOC` process.

### Reproducibility

Nothing samples without a seed: the GRM uses **every** pruned marker, `plink2 --indep-pairwise` is
deterministic, and SAIGE's step 1 seeds itself before drawing variance-ratio markers.
`params.RandomSeed` records that value in `run_manifest.json`; SAIGE offers no way to override it.
`00.prep/grm.prune.in` publishes the exact marker list, so the GRM is reconstructible.

### The SPA switch is free

SAIGE emits **both** p-values in one pass — `p.value` (saddlepoint-corrected where `Is.SPA`) and
`p.value.NA` (normal approximation). `params.PValue` selects between them in `TO_SUMSTATS`, downstream
of every expensive stage. **The default is `no_spa`**; the saddlepoint correction earns its cost when
case/control imbalance is severe, and this design's is not judged to be.

```bash
nextflow run assoc_saige.nf -resume --PValue spa         # re-runs the adapter and the figures only
```

The choice decides whether the strongest locus clears 5 × 10⁻⁸, so treat it as part of the result.

Which one a run used is recorded in `results/_run_info/run_manifest.json` (`p_value_source`,
`p_value_column`). Both are always kept in `saige_extra.tsv`, so the choice is reversible without
re-running anything expensive.

## Inputs

Read from `params.model_inputs/<cohort>/`:

| item | path |
|---|---|
| genotype | `genotype/random_model/*.{bed,bim,fam}` — **relateds retained** |
| phenotype | `phenotype/pheno.tsv` (`#FID IID PHENO1`, 1 = control, 2 = case) |
| covariates | `covariates/${params.CovarFile}` (`SEX PC1_AVG … PCn_AVG`) |

`PREP_PHENO_COV` joins the two and restricts to the genotype samples, so the row count of
`00.prep/pheno_covar.tsv` is the N the null model actually fits — a silent sample loss shows up there
rather than as a puzzling N in a SAIGE log hours later.

## Design

| parameter | default | where |
|---|---|---|
| sample sets | ordered, narrowest first where nested | `params.Cohorts` |
| covariate models | `main`, `detect` | `params.Models`, `params.ModelCovars` |
| peak-calling model | `main` **only** | `params.PeakModel` |
| p-value | **non-SPA** (`p.value.NA`) | `params.PValue` |
| Firth | off (`is_Firth_beta=FALSE`) | fixed |
| GRM | full only; no sparse GRM | fixed |
| LOCO | on, in both SAIGE steps | fixed |
| thresholds | 5 × 10⁻⁸ / 1 × 10⁻⁵ | `params.PGenomeWide`, `params.PSuggestive` |
| suggestive follow-up | the 10 strongest | `params.MaxSuggestive` |
| peak merging | distance, ±250 kb, no LD clumping | `params.PeakFlank` |
| GRM markers | LD-pruned, high-LD regions excluded | `params.Prune*`, `params.HighLdFile` |
| fine-mapping | `susie_rss`, L = 10, coverage 0.95, in-sample LD | `params.Susie*` |

## Pipeline

```
PREP_PHENO_COV    per cohort              one table for SAIGE step 1
PRUNE_MARKERS     per cohort              LD-pruned GRM markers, high-LD excluded
RELATEDNESS       per cohort              KING kinship on those same markers
   └─► COMPARE_RELATEDNESS                what the GRM is accounting for
SPLIT_BY_CHROM    per cohort x chr        step-2 input, shared by BOTH models
   │
FIT_NULL          per cohort x model      SAIGE step 1, full GRM          <- expensive
RUN_ASSOC         per cohort x model x chr  SAIGE step 2                  <- expensive
MERGE_ASSOC       per cohort x model      concatenate, keep SAIGE raw
TO_SUMSTATS       per cohort x model      canonical schema + the P choice <- cheap
   │
   └─► SCAN_PEAKS → ANNOTATE_LEADS → PLOT_SCAN
          ├─► CROSS_COHORT → COMPARE_COHORTS · COMPARE_MANHATTAN
          ├─► COMPARE_MODELS                        the reported model vs its probe
          └─► per GENOME-WIDE peak of the peak model:
                 LD_SOURCES → SUSIE → CS_VARIANTS → COLLECT_CS
                 CONDITIONAL   SAIGE step 2 --condition, re-using the null model
                 PLOT_REGIONAL · PLOT_FINEMAP · PLOT_CONDITIONAL
```

SAIGE 1.1.9 reads PLINK directly (`--bedFile/--bimFile/--famFile`), so there is **no BGEN conversion
stage**; `SPLIT_BY_CHROM` produces exactly what step 2 needs and both models of a cohort share it.

Conditional analysis **never re-fits the null model** — conditioning changes which variants are
covariates in step 2 only. Re-fitting would cost hours per round and would change the null the rounds
are compared against.

## Layout

[docs/OUTPUTS.md](docs/OUTPUTS.md) has the full tree and every column definition.
[../_shared/docs/FIGURES.md](../_shared/docs/FIGURES.md) covers the visual grammar, shared with
`assoc_plink2`; **every PNG has a companion `.md`** with the numbers behind that rendering.

```
../_shared/scripts/        peak calling, annotation, LD, fine-mapping, every figure
scripts/                   engine-specific to this component
  saige_to_sumstats.py     SAIGE -> the canonical schema; applies params.PValue
  saige_conditional.py     stepwise conditioning via SAIGE step 2 --condition
  merge_pheno_cov.py       phenotype + covariates -> one SAIGE step-1 table
  plot_model_compare.py    the reported model against its calibration probe
  relatedness.py           KING kinship on the GRM markers, by degree
  plot_relatedness.py      what retaining related samples is worth
```

## Verification

```bash
# every task completed, and the scans were NOT re-run
awk -F'\t' 'NR>1{n++; if($4!="COMPLETED"&&$4!="CACHED")b++; if($3~/^(FIT_NULL|RUN_ASSOC)/)r[$4]++} \
  END{print n" tasks, "b+0" failed"; for(k in r) print "  "k": "r[k]}' results/_run_info/trace.txt

# the adapter's join key agrees with the genotypes (it prints the match rate)
grep -h 'id-match' results/*/01.assoc/*/../../../.command.out 2>/dev/null

# N per cohort is the RANDOM-model count, not the fixed-model one
column -t -s$'\t' results/_comparison/tables/scan_qc_all.tsv | cut -c1-120

# the probe drives nothing
find results -path '*03.peaks*' -name '*.cs.tsv' | wc -l   # peak model only

# every standalone figure has its sidecar; every fan-out family has its catalogue
find results -name '*.png' | while read p; do
  d=$(dirname "$p")
  case "$d" in */genome_wide|*/suggestive) [ -f "$d/../README.md" ] || echo "NO CATALOGUE $p" ;;
                *) [ -f "${p%.png}.md" ] || echo "MISSING $p" ;; esac
done
```

A project running this component should record its own sanity anchors — the loci, credible-set sizes
and λ values it expects to be unchanged by a refactor — in [docs/STUDY_NOTES.md](docs/STUDY_NOTES.md).
