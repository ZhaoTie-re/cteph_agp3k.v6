# `saige_signals` — one academic table of the SAIGE signals and their credible sets

An integrative component (`cteph_agp3k.v6/integration/`): it reads a finished
[`analysis/assoc_saige`](../../analysis/assoc_saige) run and builds one table that a manuscript or a
project meeting can use as is — the genome-wide and cross-cohort suggestive loci, their SuSiE
credible-set variants, and for every variant × cohort the association statistics, the credible-set
status **with the reason when the variant is not in a set**, genotype counts, missingness and HWE, plus
snpEff / Ensembl annotation, ToMMo / 1000 Genomes / gnomAD effect-allele frequencies, and what the
Open Targets Platform holds on each variant (ClinVar, GWAS and molQTL credible sets — each molQTL classed cis
or trans by the variant's distance to the target gene's TSS — and enhancer-to-gene predictions; summarised in
the table, every record on a second sheet).

No association scan or fine-mapping is re-run. SuSiE is *refitted* on the stored inputs only to recover
what the upstream files do not record (why a cohort has no credible set); the refit must reproduce the
stored PIPs and sets exactly or the run stops.

## What it establishes — and what it does not

- It establishes which variants the three nested cohorts point to, how consistently, and what is known
  about them, from one source of numbers that is asserted equal to the upstream tables.
- It does **not** replicate anything: the cohorts are nested and share every case.
- It does **not** choose the P-value source. SPA vs normal approximation changes the peaks and credible
  sets themselves, so `--PValue` must equal the upstream run's `p_value_source`; otherwise the run
  aborts before submitting anything and prints the command that produces matching upstream results.

## Run

```bash
source activate dsl2
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/saige_signals
./verify.sh                            # pre-flight (lint, site rules, flags, upstream)
nextflow run saige_signals.nf          # ~60 short SLURM jobs + 1 local task, ~10 min
./verify.sh                            # post-run: DAG, checks, anchors, workbook
```

`OPENTARGETS` is the one local task: it queries the Open Targets GraphQL API (~170 queries, ~2 min)
from the host that runs `nextflow`, because compute nodes have no internet access. `-resume` reuses its
stored responses; `--OtQueryTag <date>` re-queries (e.g. after a new Open Targets release), and
`--OpenTargets false` leaves the Open Targets columns and sheet out.

Resume a later run with `nextflow run saige_signals.nf -resume <session>` (the session id from
`nextflow log`; never resume bare after a `-preview`). Only presentation changed? Editing
`scripts/write_table.py` re-runs only `WRITE_TABLE`. `nextflow lint saige_signals.nf nextflow.config`
is the syntax check.

**SPA instead of the normal approximation.** First produce SPA upstream results into their own folder,
then point this component at them:

```bash
cd ../../analysis/assoc_saige
nextflow run assoc_saige.nf -resume <assoc_saige session> --PValue spa \
    --out_dir $PWD/results.spa --run_info_dir $PWD/results.spa/_run_info
cd ../../integration/saige_signals
nextflow run saige_signals.nf --PValue spa --saige_dir ../../analysis/assoc_saige/results.spa \
    --out_dir $PWD/results.spa
```

## Deliverables

| file | content |
|---|---|
| `results/saige_signals.xlsx` | sheet `Common-variant signals`: the table, booktabs style, three header levels, lettered footnotes as short lines below; sheet `Common-variant Open Targets`: the same variants in the same order, each with every Open Targets record grouped (collapsible) below it |
| `results/saige_signals.long.tsv` | every number behind it, one row per variant × cohort, unformatted |
| `results/saige_signals.susie_sets.tsv` | one row per SuSiE fit × effect: prior variance, the effect's 95 % set, size, coverage, purity, kept or dropped |
| `results/checks.json` | counts and cross-checks of the run (selection, fits, CS status, reference coverage) |

Intermediate tables are under `results/00.selection` … `06.opentargets`; run records under
`results/_run_info`. Every column is defined in [docs/OUTPUTS.md](docs/OUTPUTS.md).

## Inputs

| input | source |
|---|---|
| loci, leads, called peaks | upstream `_comparison/tables/lead_crosscohort.tsv` (panels b, c) |
| credible-set members | upstream `_comparison/tables/cs_variants_all.tsv` |
| per-cohort P / OR / CI / PIP | upstream `<cohort>/03.peaks/<tier>/<peak>/<peak>.{pip,sumstat}.tsv`; genome-wide `<cohort>/01.assoc/<peak model>/sumstats.tsv` as fallback and cross-check |
| SuSiE refit inputs | the same peak folders: `.sumstat.tsv`, `.ld_matrix.tsv/.vars`, `.pip.tsv`, `.cs.tsv`, `.susie.json` |
| genotypes | `<model_inputs>/<cohort>/genotype/random_model/*.{bed,bim,fam}` — the set the scan ran on (path from the upstream manifest) |
| annotation | JHRP snpEff v6 index (as `analysis/assoc_gene`); Ensembl 86 via `_shared/scripts/gene_annotate.R` |
| frequencies | ToMMo 61KJPN (2025-06-16); 1000 Genomes 30x GRCh38 (NYGC) genotypes, 504 EAS / 104 JPT; gnomAD v4.1 genomes (ANNOVAR) |
| Open Targets | Platform GraphQL API (`ot_api`); release recorded in `06.opentargets/opentargets_meta.json` and footnote i |

## Key parameters

| parameter | default | meaning |
|---|---|---|
| `saige_dir` | `analysis/assoc_saige/results` | the upstream run |
| `PValue` | `no_spa` | must equal the upstream `p_value_source` |
| `Cohorts` | `narrow_mainland,intermediate_mainland,full_mainland` | upstream order |
| `CohortLabels` / `CohortCodes` | `Narrow mainland,Intermediate mainland,Full mainland` / `N,I,F` | display |
| `CohortColors` / `CohortTints` | `0571B0,008837,B2182B` / `EAF3F9,E9F5EE,FBEDEE` | the colours of the ancestry-refinement figure |
| `AfCohort` | `full_mainland` | cohort whose CTEPH / AGP3K EAF sit beside the references |
| `TableLabel` | `''` | optional label (e.g. `Table 2`) prefixed to the title and footer |
| `AfBars` | `true` | in-cell bars over the frequency columns |
| `OpenTargets` | `true` | query Open Targets; `false` leaves its columns and sheet out |
| `OtFocusTraits` | `thromb\|embol\|venous\|pulmonary hypertension\|…` | GWAS traits (regex) counted as thrombosis-related (not ranked) |
| `OtTopN` | `3` | records named per summary cell, by posterior probability |
| `OtQueryTag` | `''` | any new value re-queries the API on `-resume` |
| `python` / `rscript` | `cteph_geno_pro` / `r_work` interpreters | absolute paths: `source activate` cannot be relied on to put an env first on PATH here |
| `out_dir` | `results` | also moves `_run_info` |

## Pipeline

```
guard: upstream manifest, P source, cohorts          (Groovy, before any submission)
SELECT_SIGNALS ─┬─► ASSOC_STATS ──────────────────────────────┐
                ├─► SUSIE_REFIT × locus·cohort ─► CS_STATUS ───┤
                ├─► GENOTYPE_STATS × cohort ───────────────────┤
                ├─► ANNOTATE_FUNCTION ─────────────────────────┤
                ├─► OPENTARGETS (local) ───────────────────────┤  (records ─► WRITE_TABLE)
                └─► REF_TOMMO · REF_1KG · REF_GNOMAD ──────────┴─► ASSEMBLE ─► WRITE_TABLE ─► WRITE_RUN_MANIFEST
```

## Documents

[docs/METHODS.md](docs/METHODS.md) — selection, credible-set status, references, Open Targets, workbook, checks ·
[docs/OUTPUTS.md](docs/OUTPUTS.md) — files and columns ·
[docs/STUDY_NOTES.md](docs/STUDY_NOTES.md) — the numbers of the current run (sanity anchors) ·
`docs/anchors.no_spa.tsv` — the frozen reference values `verify.sh` compares against.
