# OUTPUTS — assoc_rvtest

The tree `results/` holds after a full run, and what every column means. Rationale is in
[METHODS.md](METHODS.md); this study's numbers are in [STUDY_NOTES.md](STUDY_NOTES.md).

```
results/
  <cohort>/
    00.callset/          bin<BIN>.split.log        plink2 log of the minAC subset
    01.prepare/          rvtest_prepare.log
    02.annotate/         snpEff logs + per-bin impact/effect stats
    03.info_filter/      *.summary.json            variants kept per stratum
    04.assoc/<stratum>/<method>/
        merged.assoc     the 12 bins concatenated — ONE table per scan
        gene_scan.tsv    every tested gene, with its tier
        gene_hits.tsv    only the genes reaching a tier (drives the fan-out)
        scan_qc.tsv      genes tested, gene-level lambda, tier counts, BH-implied P
    05.genes/<tier>/<gene_id>/
        <GENE>.<stratum>.<method>.summary.txt      the per-gene deep dive
        <GENE>.<stratum>.<method>.case.sample_details.tsv
        <GENE>.<stratum>.<method>.control.sample_details.tsv
    figures/
      01.scan/           scan.<stratum>.<method>.png  + .md sidecar
      02.gene/
        README.md        ONE catalogue for the whole family (see FIGURES.md)
        <tier>/          gene.<gene_id>.png
  _comparison/
    tables/  gene_crosscohort.tsv  scan_qc_all.tsv  gene_hits_all.tsv
    figures/ cohort_compare_genes.png + .md
  _run_info/ run_manifest.json  trace.txt  report.html  timeline.html  dag.html
```

`<tier>` is `significant` | `suggestive`. `<gene_id>` is `sig###_<GENE>` or `sug###_<GENE>`;
**the prefix is load-bearing** — `tierOf()` in `assoc_rvtest.nf` and `tier_of()` in
`_shared/scripts/figure_catalogue.py` both derive the tier from it.

## `gene_scan.tsv` / `gene_hits.tsv`

| column | meaning |
| --- | --- |
| `cohort`, `stratum`, `method` | which scan this row belongs to |
| `gene_id` | `sig###_<GENE>` / `sug###_<GENE>`; empty for a gene reaching no tier |
| `tier` | `significant` (FDR < 0.05) · `suggestive` (FDR < 0.10) · `not_a_hit` |
| `Gene` | gene symbol, as it appears in the refFlat gene model |
| `chrom`, `start`, `end` | position; parsed from `RANGE` when a method omits the columns |
| `NumVar` | qualifying variants the test collapsed |
| `NumPolyVar` | of those, the ones actually polymorphic in this cohort |
| `N_INFORMATIVE` | samples contributing to the test |
| `Pvalue` | the gene-level P, as rvtest reported it |
| `FDR_BH` | Benjamini-Hochberg over **this scan's** genes, recomputed here |

`FDR_BH` is recomputed rather than taken from `rvtest_post_process.py`, so the tier decision
has one owner; a disagreement above 1e-6 is reported as a warning.

## `scan_qc.tsv`

| column | meaning |
| --- | --- |
| `n_genes_read` / `n_genes_tested` | before and after dropping genes with no usable P |
| `n_genes_dropped_no_p` | tests that did not converge — excluded from the FDR denominator |
| `lambda_gc` | genomic-control inflation over the **gene-level** P |
| `min_p` | smallest gene P in the scan |
| `n_significant`, `n_suggestive` | tier counts |
| `p_at_significant`, `p_at_suggestive` | the raw P at which BH crossed each tier — what the scan figure draws |

## `<GENE>.<stratum>.<method>.summary.txt`

A flat `Key: Value` head followed by a `=== Variant Details ===` TSV block. The head carries
the gene-level statistics from both methods, the burden effect size (β, SE, OR, 95 % CI)
fitted by logistic regression on the pipeline's covariates, and the cumulative allele counts
and frequencies per group. The block carries one row per variant: genotype counts, case and
control frequency, ToMMo 60KJPN frequency, missingness and the snpEff impact/effect.

**Read the cumulative allele counts first.** A gene-level P computed from a dozen carriers is
a statement about a dozen people, and the odds ratio beside it has an interval spanning an
order of magnitude.

## `gene_crosscohort.tsv`

One row per (stratum, method, gene, cohort) for every gene reaching a tier in **at least one**
cohort.

| column | meaning |
| --- | --- |
| `best_tier` | the strongest tier this gene reached in any cohort |
| `called` | what THIS cohort made of it: `significant` · `suggestive` · `not_a_hit` · `not_tested` |

`not_a_hit` and `not_tested` are kept apart on purpose: the first means tested and null, the
second means the gene never entered that cohort's scan because too few qualifying variants
survived its minAC. A blank cell is never inferred.

**The cohorts are nested and share their cases**, so agreement across them is robustness to
the ancestry filter, not replication.

## `run_manifest.json`

Records the cohorts, the resolved minAC **and whether it was inherited or overridden**, the
chromosome bins, the strata and methods run, both FDR thresholds, `MinNumVar`, and the exact
phenotype and covariate column names passed to rvtest.

## Figure documents

A standalone figure (`01.scan/*.png`, `_comparison/figures/*.png`) has a companion `.md`.
The per-gene family does **not**: its prose is identical for every gene, so it is written
once as `figures/02.gene/README.md` with the genes as rows of one table. See
[../../_shared/docs/FIGURES.md](../../_shared/docs/FIGURES.md).
