# Outputs — `assoc_gene`

All paths are under `params.out_dir` (`results/`). `NA` means *not available from the source*;
`—` in the summary tables means *not defined for this cell* (gene not in the map, or a fit rvtest
reset). Three layers: the per-cohort analyst layer, the cross-cohort tables and the deliverable under
`_comparison/`, and every figure under `figures/`. `results/README.md` (written by the run) indexes
them.

```
results/
├── README.md
├── <cohort>/
│   ├── 00.prep/            pheno_rvt.tsv covar_rvt.tsv pheno_coding.json case.keep control.keep
│   ├── 01.map/             map.<s>.tsv.gz gene_index.<s>.tsv denominator.<s>.tsv
│   │                       map_collisions.<s>.tsv variants.<s>.txt anno_counts.tsv map_qc.tsv
│   ├── 02.test_inputs/<s>/ rvtest.set.chr*.txt set_stats.tsv variant_stats.tsv multi_carriers.tsv
│   │                       test_inputs_qc.tsv
│   ├── 03.assoc/<s>/<m>/   merged.tsv chr*.rvtest.log
│   ├── 04.scan/<s>/        gene_scan.tsv scan_qc.tsv map_identity.tsv
│   └── 05.signals/         significance.<m>.<s>.tsv
├── _comparison/
│   ├── tables/             *_all.tsv denominators.tsv cohort_samples.tsv gene_index_all.tsv annotation_all.tsv
│   └── robust_genes/       evidence.tsv tiers.tsv summary_table.<cohort>.tsv summary_table.md
│                           variant_detail.tsv multi_carrier_detail.tsv gene_windows.tsv
│                           cohort_samples.tsv
├── figures/
│   ├── 01.variant_sets/    variant_sets.<cohort>.png ×3 · tables/<cohort>/ · README.md
│   ├── 02.gene_scan/       gene_scan.<cohort>.<s>.png ×6 · README.md
│   ├── 03.calibration/     calibration.png · README.md
│   ├── 04.robust_genes/    robust_genes.png · README.md
│   └── 05.genes/           <GENE>.png · index.tsv · README.md
└── _run_info/              run_manifest.json trace.txt report.html timeline.html dag.html
```

`<s>` ∈ {`low_moderate_high`, `moderate_high`}; `<m>` ∈ {`cmc`, `skato`}.

---

## Per cohort

### `00.prep/`

| file | content |
|---|---|
| `pheno_rvt.tsv` | `fid iid fatid matid sex pheno1`, phenotype **1/2** (1 = control, 2 = case), `.fam` order |
| `covar_rvt.tsv` | `fid iid` + every covariate column of the model-inputs file, lower-cased headers |
| `pheno_coding.json` | `cohort n_samples n_case n_control input_coding rvtest_coding …` — the contract the RVTEST task checks rvtest's log against |
| `case.keep`, `control.keep` | `FID IID`, no header; the groups every count is taken over |

### `01.map/` (barrier B1)

| file | content |
|---|---|
| `map.<s>.tsv.gz` | the filtered map: `variant_id chrom pos ref alt impact effect biotype gene_symbol gene_id n_genes set_name`; one row per (variant, gene); only sets with ≥ `MinNumVar` variants |
| `gene_index.<s>.tsv` | `set_name gene_symbol chrom n_var_map pos_min pos_max span_bp n_high n_moderate n_low` — one row per tested set; **its row count is the denominator** |
| `denominator.<s>.tsv` | one row: `cohort stratum min_num_var n_genes_mapped n_genes_dropped_lt_minvar n_variants n_variants_multigene n_sets_chrom_split n_positions_colliding alpha bonferroni_threshold` |
| `map_collisions.<s>.tsv` | `chrom pos n_set_names set_names variant_ids` — positions shared by records of more than one set |
| `variants.<s>.txt` | the variant IDs handed to `plink2 --extract` |
| `anno_counts.tsv`, `map_qc.tsv` | annotation composition (feeds `figures/01.variant_sets`) and per-chromosome map QC |

### `02.test_inputs/<s>/`

| file | content |
|---|---|
| `rvtest.set.chr<N>.txt` | `<set_name>\t<chrom>:<pos>-<pos>,…` — one 1-bp range per mapped variant, contig string taken from the VCF |
| `set_stats.tsv` | per set: `cohort stratum set_name n_var mac mac_case mac_control n_carrier_case n_carrier_control n_case n_control n_multi_case n_multi_control max_var_per_sample_case max_var_per_sample_control n_flipped min_obs_ct` |
| `variant_stats.tsv` | per variant: `cohort stratum set_name variant_id chrom pos ref alt impact effect minor_allele flipped ac_case obs_case ac_control obs_control mac_case mac_control n_het_case n_hom_case n_carrier_case n_missing_case n_het_control n_hom_control n_carrier_control n_missing_control` |
| `multi_carriers.tsv` | one row per (set, sample) with a minor allele at ≥ 2 distinct sites: `cohort stratum set_name fid iid group n_variants_carried n_hom variant_ids genotypes` (genotypes 1 = het, 2 = hom-minor, `;`-separated in `variant_ids` order). **Contains sample IDs**; analyst file, not for the manuscript |
| `test_inputs_qc.tsv` | per chromosome: `cohort stratum chrom n_sets n_variants n_setfile_ranges vcf_chrom_string bim_chrom_string n_carrier_checks_passed all_assertions_passed` |

Invariants asserted in the task: `set(map) == set(bim)`; bim A1 == the ID's ALT; bed-decoded ALT
and observed-allele counts equal the `.acount` per variant and group; `carriers ≤ MAC ≤ 2 · n_var ·
carriers`; `n_multi ≤ carriers ≤ N`; `mac == mac_case + mac_control`.

### `03.assoc/<s>/<m>/` (barrier B2)

`merged.tsv` — the 22 rvtest parts on the `gene_scan` schema (27 columns, below) with the last four
columns `NA`; `chr<N>.rvtest.log` — rvtest's own log per chromosome, kept because the post-flight gates
read it.

### `04.scan/<s>/` (barrier B3)

**`gene_scan.tsv`** — one row per gene × method, genomic order:

| column | meaning |
|---|---|
| `cohort stratum method set_name gene_symbol chrom pos_min pos_max` | keys and position (span of the mapped variants) |
| `n_var_map` | variants in the map |
| `n_var_engine` | rvtest's `NumVar` (may exceed `n_var_map` by the span-overlap excess) |
| `n_informative` | rvtest's `N_INFORMATIVE` (asserted == N in the task) |
| `pvalue` | CMC **score** p, or SKAT-O p, verbatim |
| `rho` | SKAT-O ρ (NA on cmc rows) |
| `beta se or or_l95 or_u95` | cmcWald burden coefficient, SE, exp(β), exp(β ∓ 1.96 SE) (NA on skato rows, NA when rvtest reset the fit) |
| `mac mac_case mac_control n_carrier_case n_carrier_control` | from `set_stats.tsv`, identical on both method rows |
| `bonferroni_threshold` | copied from `denominator.<s>.tsv` |
| `significant_bonferroni` | 1 / 0, NA when there is no p |
| `fdr_bh` | BH q over `n_genes_mapped` (missing genes at p = 1) |
| `significant_bh` | 1 / 0 = **called**, NA when there is no p |

**`scan_qc.tsv`** — one row per method: `cohort stratum method n_genes_mapped n_genes_returned
n_genes_missing lambda_gc lambda_gc_nvar2 lambda_gc_nvar_ge3 min_p n_significant_bonferroni
n_significant_bh bh_threshold bonferroni_threshold min_num_var`. `bh_threshold` is the largest p BH
rejected (NA when none). The figures draw `lambda_gc` only; the two split columns are tabulated in
the figure READMEs.

**`map_identity.tsv`** — every set in the union of the gene index and rvtest's output: `cohort
stratum chrom set_name n_var_map n_var_engine span_overlap_excess ok reason`, `reason ∈ {ok,
explained_by_span_overlap, explained_by_collision, missing, unexplained}`. `SCAN` aborts on `missing` /
`unexplained` unless `--no-strict`.

### `05.signals/`

`significance.<m>.<s>.tsv` — three rows (`bonferroni`, `bh`, `nominal`): `threshold value basis
primary comparable_across_methods cohort stratum method n_genes_mapped n_genes_returned
n_genes_missing min_num_var lambda_gc min_p n_significant alpha`. `primary` = a hit can come from this
row; `comparable_across_methods` = the value is the same number for cmc and skato. `nominal` is
reference only and never a decision rule.

---

## Across cohorts — `_comparison/`

### `tables/`

Concatenations of the per-cohort files with their own keys (`gene_scan_all.tsv`, `scan_qc_all.tsv`,
`significance_all.tsv`, `map_identity_all.tsv`, `denominators.tsv`, `set_stats_all.tsv`,
`variant_stats_all.tsv`, `multi_carriers_all.tsv`), plus:

| file | content |
|---|---|
| `cohort_samples.tsv` | `cohort n_samples n_case n_control` from the three `pheno_coding.json` |
| `gene_index_all.tsv` | `cohort stratum` + the gene-index columns — **the designed experiment**, from which `not_in_map` is decided |
| `annotation_all.tsv` | `cohort` + `anno_counts.tsv` |
| `scan_scale.tsv` | `stratum y_data_max source_cohort source_method n_cohorts n_methods n_tests` — the −log10 P maximum each variant set reached in ANY cohort, so its three scan figures share a y axis. The threshold lines stay each cohort's own |

### `robust_genes/` — the deliverable

**`evidence.tsv`** — every (stratum, gene in the union of the three maps) × cohort × method:
`gene_symbol set_name chrom cohort stratum method state pvalue nlp fdr_bh n_var_map mac_case
mac_control n_carrier_case n_carrier_control rho or or_l95 or_u95 lambda_gc`;
`state ∈ {called_bonferroni, called_bh_only, not_significant, not_in_map}`. A mapped gene rvtest did not
return is `not_significant` (it entered the BH family at p = 1).

**`tiers.tsv`** — one row per (gene, stratum) with tier ≥ 1: `gene_symbol set_name chrom stratum tier
n_cohorts_called cohorts_called n_cohorts_called_cmc n_cohorts_called_skato n_cohorts_both_methods
called_bonferroni_anywhere n_cohorts_mapped cohorts_not_in_map min_p_cmc min_p_skato best_cohort
n_multi_site_carriers_case_max n_multi_site_carriers_control_max note`. `note` carries
`cmc_or_na:<cohorts>` and/or `multi_site_carriers`.

**`summary_table.<cohort>.tsv`** (one per cohort) and **`summary_table.md`** (`params.ReportCohort`) —
Tier 1 and 2 genes, one row per gene, one block per stratum. Block columns, prefixed `<stratum>.`:

| column | format | source |
|---|---|---|
| `tier` | `1` / `2` / `3` / `0` / `—` | the stratum's own tier; `0` = tested but called in < 2 cohorts; `—` = not in this cohort's map |
| `n_variants` | integer | `gene_index` `n_var_map` |
| `carriers_case`, `carriers_control` | `10 (2.28%)` = carriers (% of the group) | `set_stats` |
| `mac_case`, `mac_control` | integer | `set_stats` |
| `maf_case`, `maf_control` | `0.011` = MAC / (2 N<sub>group</sub>) | `set_stats`, `cohort_samples` |
| `skato_rho` | `0.3` | `gene_scan` skato row |
| `skato_p` | `1.41E-07` | `gene_scan` skato row |
| `cmc_or_95ci` | `9.60 (3.56-25.86)` | `gene_scan` cmc row (cmcWald) |
| `cmc_p` | `7.23E-08` | `gene_scan` cmc row (CMC score p) |
| `multi_site_carriers` | `4 / 22` = case / control samples with ≥ 2 sites | `set_stats` |

`—` fills a block when the gene is not in that cohort's map for the stratum. The `.md` carries the
definitions, the tier legend and the nested-cohort caveat under the tables.

**`gene_windows.tsv`** — `gene_symbol set_name chrom start end`, one row per Tier 1/2 gene: the
window whose Ensembl 86 gene models the per-gene figure draws (the variant span padded by
max(10 %, 20 kb)).

**`variant_detail.tsv`** — `variant_stats_all.tsv` rows for every tiered (gene, stratum), all cohorts.
**`multi_carrier_detail.tsv`** — `multi_carriers_all.tsv` rows for every tiered (gene, stratum); sample
IDs inside.

---

## `figures/` — see `docs/FIGURES.md`

Five families in reading order, each with one `README.md` (the eight-section document: question,
panels, interpretation, values, full statistics, how to read, what it does not establish, model).

| family | files | panels |
|---|---|---|
| `01.variant_sets/` | `variant_sets.<cohort>.png` ×3, `tables/<cohort>/{annotation_summary,stratum_summary}.tsv`, `README.md` (one row per cohort) | (a) call set by impact, (b) consequences by impact, (c) gene sets per stratum, (d) set-size distribution |
| `02.gene_scan/` | `gene_scan.<cohort>.<s>.png` ×6, `README.md` (one row per scan: thresholds, calls, λ overall / n_var = 2 / ≥ 3, Spearman ρ, named genes) | (a) CMC Manhattan, (b) SKAT-O Manhattan, (c) QQ with the null 95 % band and λ<sub>GC</sub>, (d) CMC against SKAT-O per gene |
| `03.calibration/` | `calibration.png`, `README.md` | 2 variant sets × 3 cohorts QQ grid, λ<sub>GC</sub> per statistic per cell |
| `04.robust_genes/` | `robust_genes.png`, `README.md` | (a) evidence matrix with the tier per variant set, (b) CMC OR (95 % CI) in the report cohort, (c) carrier frequency control → case with counts |
| `05.genes/` | `<GENE>.png` per Tier 1/2 gene, `index.tsv` (`gene_symbol set_name chrom pos_min pos_max tiers stratum_drawn n_variants_drawn n_variants_per_cohort n_genes_in_model_panel n_multi_site_carriers_case n_multi_site_carriers_control png`), `README.md` (one row per gene + the cells and variants drawn) | (a) three strips, one per cohort, of per-variant carrier percentage against that cohort's own N on one shared scale, (b) Ensembl 86 gene models on the same axis, (c) table-style OR forest with CMC P and SKAT-O P columns coloured by state. `n_variants_per_cohort` says how many of the drawn variants each cohort's own map holds — a narrower cohort can lack some, and that variant is then simply absent from its series |

The `*.stats.json` beside the members of 01 and 02 are build artefacts of the catalogue step.

### `_run_info/`

`run_manifest.json` records the cohorts, minAC, strata, model, methods, decision rules and tier rule of
the run; `trace.txt` is what `verify.sh` §9 reads.
