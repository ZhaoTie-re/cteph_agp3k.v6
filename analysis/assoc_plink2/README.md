# `assoc_plink2` — association scan with two-tier peak calling

Genome-wide logistic association across **N sample sets × M genetic models**, with one model primary.
The primary model (`params.PeakModel`, additive by default) is the only one that defines peaks; the
others are reported as scans — Manhattan + QQ, calibration, and their own genome-wide loci with full
annotation — but nothing fans out from them.

The component is **dataset-independent**: cohorts, covariates, thresholds, reference resources and
every display label are parameters. [`docs/METHODS.md`](docs/METHODS.md) states the method;
[`docs/STUDY_NOTES.md`](docs/STUDY_NOTES.md) records what one project's run configured and found.

Peaks come in two tiers, treated differently on purpose:

| tier | threshold | gets |
|---|---|---|
| `genome_wide` | `params.PGenomeWide` | annotation · cross-cohort table · conditional analysis · SuSiE · three-source regional plot |
| `suggestive` | `params.PSuggestive` | annotation · cross-cohort table · one landscape figure — **no per-peak follow-up** |

The suggestive tier describes the shape of the scan; it is not a list of findings. Over *M* analysed
variants a threshold α yields ≈ *M*α crossings by chance alone, and an inflated scan yields a multiple
of that.

This component is a **fixed-effects scaffold**: it establishes the scale of the signal and produces
candidate peaks. The random-effect (GRM) follow-up its calibration calls for is SAIGE / REGENIE and is
out of scope here.

---

## Run

```bash
source activate dsl2
nextflow run assoc_plink2.nf -resume
```

Everything environment- or project-specific lives in one **`SITE CONFIGURATION`** block at the top of
`assoc_plink2.nf`. Override on the command line (`--key value`) or with `-params-file cfg.yaml`;
nothing outside that block needs to change to run the component elsewhere.

> **Cost note.** `RUN_ASSOC` is by far the most expensive process. Its task hash depends on the
> interpolated script text and the input file contents, so editing that process block, `cohortInput()`,
> or the **values** of the params they resolve (`plink2`, `conda_env`, `FirthMode`, `out_dir`,
> `covarNameArg`, `covarTag`) re-runs every scan. Everything else is cheap to change. `work/` must not
> be hand-deleted before a resume; `results/` may be deleted freely, as `publishDir` re-publishes on
> resume including for cached tasks.
>
> After a resume, confirm the log lists no `RUN_ASSOC` process and the scheduler shows no
> `nf-RUN_ASSOC` job. If it does, stop and diagnose before letting the scans re-run.

## Inputs

Read from `params.model_inputs/<cohort>/`:

| item | path |
|---|---|
| genotype | `genotype/fixed_model/*.maf_ge_threshold.{bed,bim,fam}` |
| phenotype | `phenotype/pheno.tsv` (`#FID IID PHENO1`, 1 = control, 2 = case) |
| covariates | `covariates/${params.CovarFile}` (`SEX PC1_AVG … PCn_AVG`) |

Phenotype and covariate files may be a full-cohort superset; plink2 intersects by ID, so the analysis N
is the genotype's.

`params.Cohorts` is an **ordered** list. Where the sets are nested, list them narrowest first — figures
use this order top to bottom, and the cross-cohort logic assumes an estimate exists in every set.

## Design

| parameter | default | where |
|---|---|---|
| peak-calling model | additive **only** | `params.PeakModel` |
| models scanned | additive, dominant, recessive | `params.Models` |
| covariates | `SEX` + `PC1_AVG–PCn_AVG` | `params.NPcs`, `params.PcLabel`, `params.CovarLabel` |
| Firth | `no-firth` | `params.FirthMode` |
| thresholds | 5 × 10⁻⁸ / 1 × 10⁻⁵ | `params.PGenomeWide`, `params.PSuggestive` |
| peak merging | distance, ±250 kb, no LD clumping | `params.PeakFlank` |
| fine-mapping | `susie_rss`, L = 10, coverage 0.95, in-sample LD | `params.Susie*` |
| LD panel names | shown on the regional figure | `params.LdPanelLabels` |

`params.CovarLabel` is the human-readable covariate set printed on every figure and sidecar. It
defaults to `SEX + <NPcs> <PcLabel> PCs`; set it explicitly when the covariates are not simply SEX+PCs,
so no figure can describe covariates the run did not fit.

## Pipeline

```
RUN_ASSOC ×(cohorts × models)          the expensive step
   │
   └─► SCAN_PEAKS ×cohorts   QC + N for every model; primary-model peaks (2 tiers);
          │                  peaks of EVERY model, 2 tiers (model_peaks.tsv)
          │
          ├─► ANNOTATE_LEADS       lead_annotation.tsv         (primary peak leads)
          │      │                 model_peaks_annotation.tsv  (every model's peaks)
          │      └─► PLOT_SCAN     one figure per scan: Manhattan + ITS OWN peaks
          │                        on a shared axis + QQ + effect/MAF
          │
          ├─► CROSS_COHORT ×1      every peak lead, in EVERY cohort (barrier)
          │
          └─► per GENOME-WIDE peak
                 ├─ LD_SOURCES   in-sample / population panel / reference panel,
                 │               + exon-level genes + recombination
                 ├─ SUSIE ─► CS_VARIANTS ─► COLLECT_CS ×1
                 ├─ CONDITIONAL
                 └─ PLOT_REGIONAL · PLOT_FINEMAP · PLOT_CONDITIONAL

          COMPARE_COHORTS ×1     calibration + composition + cross-cohort forest
          COMPARE_MANHATTAN ×1   every cohort's scan stacked on one genomic axis
```

The genome-wide fan-out is driven by `lead_variants.tsv`, filtered on `tier`. A cohort with no
genome-wide peak simply produces no follow-up tasks — no process needs to know about tiers.

`CROSS_COHORT` and `COMPARE_MANHATTAN` are **barriers**, deliberately: their inputs are unions or
shared axes over all cohorts, so none can be processed until every cohort has finished.

**There is no power / minimum-detectable-effect analysis.** Removed by decision — see
[docs/METHODS.md](docs/METHODS.md) §13.

## Layout

[docs/OUTPUTS.md](docs/OUTPUTS.md) has the full tree and every column definition.
[docs/FIGURES.md](docs/FIGURES.md) covers the shared visual grammar; **every PNG has a companion
`.md`** with the numbers behind that rendering and how to read it.

```
scripts/
  plot_style.py          shared style, LD colour scale, symbol glossary, layout helpers
  figure_doc.py          writes the sidecar .md beside every figure
  call_peaks.py          scan QC + peak calling, 2 tiers, every model
  annotate_leads.py      the lead-variant table + every model's peaks
  gene_utils.R           ONE definition of "an informative gene" + representative transcript
  gene_annotate.R        nearest/overlapping gene per lead
  region_tracks.R        exon structure + recombination for one window
  ld_sources.py          three LD sources + the signed r matrix SuSiE uses
  susie_finemap.R        SuSiE on summary statistics + in-sample LD
  conditional_stepwise.py
  variant_annot.py       ONE implementation of rsID / gene / snpEff / genotype lookups
  cs_variants.py         per-variant credible-set table
  cross_cohort.py        every peak lead, in every cohort
  plot_manhattan_qq.py   the per-scan figure (Manhattan + peaks + QQ + effect/MAF)
  plot_regional.py · plot_finemap.py · plot_conditional.py
  plot_cohort_compare.py · plot_cohort_manhattan.py
```

## Environments

| use | param |
|---|---|
| Nextflow | (activated before launching) |
| python + PLINK | `params.conda_env` |
| R (susieR, EnsDb, rtracklayer) | `params.conda_env_r`, invoked by **absolute path** |

`params.rscript` pins the R interpreter absolutely on purpose: `source activate <r env>` on top of an
already-active env leaves the first env's `Rscript` ahead on `PATH`, which runs the R scripts against a
library that has none of the required packages.

## Threads

Thread counts are set **explicitly**, and both the tool's thread argument and the scheduler request
come from one table in `nextflow.config`:

```groovy
params.threadsRunAssoc = 8
withName: 'RUN_ASSOC' {
    ext.threads    = params.threadsRunAssoc
    clusterOptions = "--rsc p=1:t=${params.threadsRunAssoc}:c=…:m=…"
}
```

and every tool reads `${task.ext.threads ?: 1}`. Two traps this avoids:

- **Do not use Nextflow's `cpus` directive** if the site's `sbatch` rejects `--cpus-per-task`
  (`error: forbidden option, cpus-per-task`). The resource request form in `clusterOptions` is
  site-specific; the one above is an example, not a portable default.
- **Do not fall back on `task.cpus`.** With `cpus` unset it reports Nextflow's default of **1**, not
  `null`, so `${task.cpus ?: 8}` silently runs single-threaded inside an 8-core allocation — a large
  slowdown with correct output and no warning.

Verify with `nextflow config`; `-preview` does not evaluate resource blocks and will not catch an
error there. `-preview` also opens a new, empty cache session, so a later bare `-resume` re-runs
everything.

## Verification

```bash
# every task completed, and the scans were NOT re-run
awk -F'\t' 'NR>1{n++; if($4!="COMPLETED"&&$4!="CACHED")b++; if($3~/^RUN_ASSOC/)r[$4]++} \
  END{print n" tasks, "b+0" failed"; for(k in r) print "  RUN_ASSOC "k": "r[k]}' results/_run_info/trace.txt

find results -name '*.glm.logistic' | wc -l                        # cohorts × models
grep -h "covariates loaded" results/*/01.assoc/*/*.log | sort -u   # one line: 1 + NPcs

# scan QC, including the degenerate-fit guard
column -t -s$'\t' results/_comparison/tables/scan_qc_all.tsv | cut -c1-140

# tiers
awk -F'\t' 'NR>1{n[$3]++} END{for(k in n) print k, n[k]}' results/_comparison/tables/peaks_all.tsv

# the annotation table
head -1 results/_comparison/tables/lead_annotation_all.tsv | tr '\t' '\n' | nl

# every figure has its sidecar, and none is stale
find results -name '*.png' | while read p; do [ -f "${p%.png}.md" ] || echo "MISSING $p"; done
```

A project running this component should record its own sanity anchors — the loci, credible-set sizes
and λ values it expects to be unchanged by a refactor — in `docs/STUDY_NOTES.md`.
