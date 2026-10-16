# `gene_signals` — one workbook of the gene-based rare-variant signals

An integrative component (`cteph_agp3k.v6/integration/`). It reads two inputs:

- a finished [`analysis/assoc_gene`](../../analysis/assoc_gene) run;
- the [`tuning.rv`](../../tuning.rv) call sets and QC that the run tested.

From them it builds one workbook. The workbook follows the common-variant table
([`saige_signals`](../saige_signals)) wherever the two analyses share a question. Rare-variant content is
added only where the gene-based analysis needs it.

Each sheet comes after the sheets it refers to: input QC → what the tests were given → the genome-wide scans →
the genes → their variants → which variants carry the signal → external evidence on them → the individuals.

| sheet | one row per | content |
|---|---|---|
| Rare-variant QC | QC step | The rare-variant QC table: input → sample QC → variant QC → allele count, and the basis of the minor allele count threshold κ, recomputed from the tuning.rv outputs |
| Variant sets | step × cohort | Each step on one row, its variants and the genes they fall in side by side: analysed variants by highest snpEff impact; mapped to genes or not (and why); per variant set, qualifying → not tested / tested, tested genes with 2 and ≥ 3 variants |
| Consequence terms | impact × term | Every snpEff consequence term within each impact, with the impact's total: analysed variants, and variants of the tested genes per set, per cohort |
| Gene-based scans | cohort × variant set | Primary analysis (genes with ≥ 2 variants) and sensitivity analysis (genes with ≥ 3 variants) side by side: genes tested, Bonferroni threshold, and for SKAT-O and CMC the smallest P, significant genes (Bonferroni, BH), the largest BH-significant P and λ<sub>GC</sub>, as observed |
| Gene-based signals | gene × variant set × cohort | Per gene, the two variant sets stacked, each with its three cohorts: variants, carriers n (%), cumulative MAC / MAF, multi-site carriers, SKAT-O ρ / P / q, CMC OR (95% CI) / P / q, consistency tier; and the sensitivity analysis (≥ 3 variants): SKAT-O q, CMC q, tier. BH calls shaded |
| Gene-based variants | gene × variant | The columns of the common-variant table: rsID, Chr, Position, EA/OA, Ensembl 86 gene, snpEff consequence / impact / transcript / HGVS / other consequences, Open Targets (ClinVar, GWAS credible sets, L2G, molQTL credible sets, enhancer-to-gene), EAF in this study / ToMMo / 1KG / gnomAD. Per cohort: genotypes, EAF, missing, HWE P, and whether the cohort's gene test used the variant (and why not) |
| Leave-one-variant-out | gene × variant | Every gene in both variant sets side by side, the tiered sets marked: −log10 P of SKAT-O and CMC per cohort with all variants (grey gene row) and with each variant left out (rows below, in the order of Gene-based variants), shaded red for a fall and blue for a rise; carriers per cohort. 3 cohorts × 2 statistics × 2 sets. Descriptive: no multiple-testing correction |
| Gene-based Open Targets | variant → record | Every Open Targets record behind the summary columns, as the common-variant table's second sheet |
| Multi-site carriers | gene × sample | The samples behind the multi-site counts, by sample ID: the variants they carry (by impact, het / hom), the sets and cohorts counting them, project, platform, target and observed depth, sex, age at DNA collection, and for cases age at diagnosis, DVT and max mPAP |

The calls of the primary analysis use upstream's gene family, the genes with ≥ 2 qualifying variants. A
sensitivity analysis judges the same P values over the genes with ≥ 3 variants (`--SensMinNumVar 3`).

The only tests run here are the leave-one-variant-out re-tests. They use the upstream inputs and options, and
each gene's test on all its variants is asserted equal to the published one. Every count is recomputed from the
genotypes and maps the tests used, and asserted equal to the upstream tables before it is written. The workbook states what was observed and adds no
interpretation. No P value is genomic-control corrected.

## Run

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/gene_signals
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH ./verify.sh     # pre-flight
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH nextflow run gene_signals.nf
PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH ./verify.sh     # post-run: DAG, checks, anchors, workbook
```

- **Jobs.** 23 short SLURM jobs (6 of them the leave-one-variant-out re-tests, cohort × variant set) and 2
  local tasks (Open Targets, the run manifest).
- **Open Targets.** `OPENTARGETS` is the local task: it queries the Open Targets GraphQL API from the host
  that runs `nextflow`, because compute nodes have no internet access. `-resume` reuses the stored
  responses. `--OtQueryTag <date>` re-queries; `--OpenTargets false` leaves the Open Targets columns and
  sheet out.
- **Resume.** Use `nextflow run gene_signals.nf -resume <session>`, with the session id from
  `nextflow log`. Never resume bare after a `-preview`.
- **Syntax check.** `nextflow lint gene_signals.nf nextflow.config`. Do not use
  `nextflow run … --help`: it launches the pipeline.

## Files

```
gene_signals.nf        the pipeline (upstream guard, 14 processes)
nextflow.config        SLURM resources (--rsc only), trace / report
verify.sh              pre-flight and post-run checks
scripts/               one script per step (docs/METHODS.md)
docs/METHODS.md        what each step computes and asserts
docs/OUTPUTS.md        every output file and column
docs/STUDY_NOTES.md    decisions behind the design, and what the build showed
```
