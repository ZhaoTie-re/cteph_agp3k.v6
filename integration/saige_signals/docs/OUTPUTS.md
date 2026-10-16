# saige_signals — outputs

```
results/
  saige_signals.xlsx               sheets "Common-variant signals" (the table) and "Common-variant Open Targets" (every record)
  saige_signals.long.tsv           one row per variant x cohort
  saige_signals.susie_sets.tsv     one row per SuSiE fit x effect
  checks.json                      counts and cross-checks
  00.selection/   loci.tsv variants.tsv sel.tsv cs_members.tsv pairs.tsv selection.json
  01.association/ assoc.tsv
  02.finemap/     <cohort>/<cohort>.<locus>.{effects.tsv,variants.tsv,fit.json}  assoc_status.tsv  cs_status.json
  03.genotypes/   geno.<cohort>.tsv
  04.annotation/  function.tsv
  05.frequencies/ tommo.tsv kg.tsv gnomad.tsv
  06.opentargets/ opentargets.tsv opentargets_detail.tsv opentargets_meta.json raw/<variant>.json
  _run_info/      run_manifest.json trace.txt report.html timeline.html dag.html
```

## The workbook

Sheet `Common-variant signals`: rows 1–2 title and descriptor, rows 3–5 the header (block spanners · sub-spanners ·
leaf headers), one row per variant, then lettered footnotes (a–i) and abbreviations below, one short line
per row in column A. Loci are separated by a thin
rule; the locus name and class appear on its first row; credible-set variants are an outline group under
their lead. Columns:

| block | columns |
|---|---|
| Variant | Locus · Class (GW / SG) · Role (Lead / Alt. lead / CS variant) · rsID · Chr · Position (GRCh38) · EA/OA · Gene (Ensembl 86) · CS in (cohort letters in cohort colours) |
| Functional annotation | Consequence · Impact · Transcript (snpEff gene · transcript) · HGVS (`c.… (p.…)`) · Other consequences |
| Open Targets | ClinVar (n · significance: condition) · GWAS credible sets (n (n thrombosis-related) · trait (PP)) · L2G (gene (median score)) · molQTL credible sets (n (k trans) · cis- / trans-type gene (tissue / cell type, PP); METHODS §6) · Enhancer-to-gene (n genes · gene score (best-scoring tissue / cell type, n tissues / cell types)); the number of records, then up to `OtTopN` (3) of them by PP (E2G by score, ClinVar alphabetical), full names, one unwrapped line per cell (select a cell to read text beyond the column); "—" = no record, across all five = not in the Open Targets variant index |
| Effect-allele frequency | This study: CTEPH · AGP3K (`AfCohort`); Japanese: ToMMo · 1KG JPT; East Asian: 1KG EAS · gnomAD EAS; European: gnomAD NFE |
| per cohort | Peak (GW / SG where the variant is the cohort's peak lead; – otherwise; blank on CS rows) · OR (95% CI) · P · PIP · CS (`CS1 (33)`, `Not in CS1 (PIP < 0.09489)`, `No CS: no signal`, `No CS: low LD (251 variants; min |r| 0.00)`; METHODS §3) · Genotypes CTEPH / AGP3K (OA/OA / OA/EA / EA/EA / missing) · EAF · Missing · HWE P |

Sheet `Common-variant Open Targets` (with `OpenTargets = true`): rows 1–2 title and descriptor, row 3 the header (note
letters as superscripts), then one row per variant of sheet `Common-variant signals`, same variants in
the same order (grey band; Role, Variant, rsID, then the number of records per source, molQTL with the
number in trans, "Not in the Open Targets variant index" or "no record"), each followed by its records,
under a header row per locus (locus,
number of variants by role, number in the Open Targets variant index). Outline, summary rows above: locus
header rows level 0, every variant row level 1 (Lead and Alt. lead included, so each variant's records fold
on their own), record rows level 2, collapsed when the file opens (buttons 1 loci, 2 variants, 3 records).
Records under a variant: ClinVar (alphabetical), GWAS and molQTL (by PP), E2G (by score); the first ones
are those named in `Common-variant signals`; thrombosis-related GWAS traits in bold. Autofilter on the header; record
rows repeat the variant ID (grey) for filtered views. Columns: Locus · Role · Variant · rsID · Source · Type (QTL type) · Trait,
condition or target gene · cis / trans (to target TSS) (molQTL only: `cis (15 kb)`, `trans (gene on chr2)`,
`neither (…)`, `unknown`; METHODS §6) · QTL tissue / cell type (of the molQTL study) · Study · Publication (first
author, year) · Study N · Cohorts (as Open Targets lists them; grey italic, the study's project for molQTL studies
listing none; — none) · Population (% of sample) (share of the study sample by gnomAD reference population, e.g.
`NFE 96%, EAS 2%, AFR 2%`; — not reported) · P · Beta · P / Beta from · Set lead variant · PP (the variant's posterior
probability in the set) · Fine-mapping (set) · L2G gene (set) · L2G score (set) · E2G score · E2G best tissue / cell type · E2G tissues / cell types (n, predicting
the link) · E2G interval · Clinical significance · Review status (ClinVar). Below the table, one blank row, then
lettered notes a–f and the abbreviations as short lines in column A, wrapped to the frozen columns A:D (as in
`Common-variant signals`): a records and their order, b study / cohorts / population, c P and Beta, d cis / trans,
e PP and the set-level columns, f E2G.

## `saige_signals.long.tsv`

| column | meaning |
|---|---|
| `order`, `locus_id`, `label`, `panel`, `tier` | locus order, upstream locus id, gene label, cross-cohort panel (b/c), class (GW/SG) |
| `role`, `role_rank` | Lead / Alt lead / CS member |
| `variant_id`, `chrom`, `pos`, `EA`, `OA` | chr:pos:REF:ALT; EA = ALT, OA = REF |
| `rsID`, `rsID_unresolved`, `rsID_source` | accession attributed to the allele, or the unattributable accession string; source: `ToMMo`, `dbSNP allele` (dbSNP lists it for this allele) or `dbSNP exclusion` (the only candidate dbSNP does not list for another allele) |
| `Gene_EnsDb`, `Gene_Biotype`, `Gene_Distance_bp` | Ensembl 86 gene, biotype, signed distance |
| `snpEff_*` | most severe annotation (Effect, Impact, Gene, FeatureID, Biotype, HGVS_c, HGVS_p), `LoF_NMD`, number of annotations, every effect term |
| `max_pip`, `cs_cohorts` | largest PIP over the cohorts' credible sets; cohort letters of the sets containing the variant |
| `cohort`, `peak_id`, `peak_tier`, `called_peak` | the cohort, its peak at the locus ('' = none), its tier, the variant's call there (leads only) |
| `fine_mapped`, `fm_reason` | False with `no_peak` / `outside_window` / `not_tested` when there is no PIP |
| `cohort_has_cs`, `n_cs`, `in_cs`, `cs_index`, `cs_size`, `cs_cut`, `cs_cuts_all` | credible-set membership; `cs_cut` = smallest PIP in the variant's set (or in the sets when outside) |
| `pip`, `pip_fitted` | PIP in that cohort; False when every PIP of the window is 0 |
| `cs_status`, `cs_reason`, `cs_reason_short` | status code (METHODS §3) and plain-text reason (lead rows / member rows) |
| `n_effects`, `n_sets_pre`, `n_sets_pure` | effects with V > 1e-9; 95 % sets before / after the purity filter |
| `dropped_set_size`, `dropped_set_min_r` | `no_cs_impure` fits: size and purity (min \|r\|) of the strongest effect's dropped 95% set (largest V; the only effect of a single-effect fit) |
| `alpha_pure`, `alpha_cut_pure` | the variant's α in, and the α cut of, the nearest surviving set |
| `P`, `OR`, `L95`, `U95`, `Beta`, `SE`, `A1_FREQ`, `OBS_CT` | SAIGE statistics (P as chosen upstream) |
| `n_case`, `n_ctrl`, `case_*`, `ctrl_*` | group sizes; genotype counts, `*_geno` string, EAF, missing rate, HWE P |
| `af_cteph`, `af_agp3k` | case / control EAF in `AfCohort` |
| `tommo_*`, `tommo60_af` | ToMMo 61KJPN AF / AN / AC / FILTER / ID / match mode; 60KJPN cross-check |
| `kg_*` | 1000 Genomes EAS and JPT AF / AN, FILTER, match mode, LD-panel cross-check |
| `gnomad_af`, `gnomad_eas_af`, `gnomad_nfe_af` | gnomAD v4.1 genomes |
| `external_af_note` | why an external frequency is missing or flagged |
| `ot_in_index`, `ot_rsids` | variant in the Open Targets variant index (empty when not queried); its rsIDs there |
| `ot_clinvar_n`, `ot_clinvar` | ClinVar records; "significance: condition" summary |
| `ot_gwas_n`, `ot_gwas_focus_n`, `ot_gwas_top`, `ot_l2g` | GWAS credible sets containing the variant, of which thrombosis-related; top traits (PP); most frequent top L2G gene (median score) |
| `ot_qtl_n`, `ot_qtl_trans_n`, `ot_qtl_top` | molQTL credible sets; of which trans; top records (`cis-` / `trans-` prefix) |
| `ot_e2g_n`, `ot_e2g_top` | enhancer-to-gene target genes; top records |

## `saige_signals.susie_sets.tsv`

One row per effect of every refit: `cohort`, `locus_id` (the peak id), `lead_id`, `n_gwas`, `n_variants`,
`L`, `niter`, `converged`, `n_effects`, `n_sets_pre`, `n_sets_pure`, `n_credible_sets_stored`,
`max_abs_dpip`, `effect`, `V`, `active`, `dedup_of`, `set_size`, `set_coverage`, `min_abs_corr`,
`mean_abs_corr`, `median_abs_corr`, `purity_pass`, `stored`, `contains_lead`, `top_variant`, `top_alpha`,
`n_table_variants`, `table_variants`.

## Intermediate tables

`pairs.tsv` (locus × cohort → peak, variant list) drives the SuSiE refits; `assoc.tsv` → `assoc_status.tsv`
add the credible-set status; `geno.<cohort>.tsv`, `function.tsv`, `tommo.tsv`, `kg.tsv`, `gnomad.tsv` are
per-source tables joined by `ASSEMBLE`. `selection.json` and `cs_status.json` are merged into `checks.json`.
`opentargets.tsv` (one row per variant) is joined by `ASSEMBLE`; `opentargets_detail.tsv` (one row per
record, the columns of sheet `Common-variant Open Targets` (credible-set records: `p_value` / `beta` as shown,
`p_source` / `beta_source` = variant | lead | '', `lead_variant`, `lead_rsid`, `lead_is_variant`, and the raw
`variant_p`, `variant_beta`, `lead_p`, `lead_beta`, the study's `cohorts` and `ld_population` (empty when Open
Targets lists none); molQTL records: `target_chrom`, `target_tss`,
`qtl_distance` (|variant − TSS|, empty on another chromosome), `qtl_class` (cis / trans / neither / unknown),
`ot_is_trans` (Open Targets' set-level flag)) plus `project`, `condition`, `in_95cs`, `distance_to_tss` (E2G),
`n_predictions`, `datasource`) and `opentargets_meta.json` (API / data release, query settings, counts)
go to `WRITE_TABLE`; `raw/` keeps the API responses.
