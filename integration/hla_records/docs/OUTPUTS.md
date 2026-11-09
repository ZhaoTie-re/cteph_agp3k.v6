# hla_records — outputs

```
results/
  hla_records.xlsx            the workbook (22 sheets)
  workbook.json               its layout: per sheet, the blocks with header row, first body row, columns, row keys
  workbook_cells.json.gz      what every body cell shows (the audit compares the file with it)
  audit.json                  the audit of the written workbook, with the md5 of the file it audited
  00.upstream/
    typing_verify.log         the typing component's verify.sh, run read only
    upstream.json             the upstream runs and their verification; every file read (sheets, path, date, md5)
  01.typing/                  the tables of sheets T1–T8 (one TSV per block) and typing_facts.json
  02.association/             the tables of sheets A1–A12 and assoc_facts.json
  03.decks/
    decks.json                every slide of both decks: number, key, section, title
    deck_check.json           agreement with the decks' data, per sheet, every mismatch listed
  04.checks/
    recheck.json              the independent recomputation, per check, every mismatch listed
  _run_info/                  run_manifest.json (written after the audit passed), trace, report, timeline, DAG
```

TSVs: tab-separated, a header row, `''` for missing; P-like columns as `%.6e`, other floats at ten significant
digits. Cohorts are `narrow_mainland`, `intermediate_mainland`, `full_mainland`; comparisons `all`,
`dvt_yes`, `dvt_no`, `dvt_case_only`; markers as upstream (`A*02:01`, `DQB1:57:V`, positions `DQB1:57`).

## The workbook

| sheet | built from | one row per |
|---|---|---|
| Contents | `03.decks/decks.json`, `assets/slide_map.tsv` | sheet; slide |
| T1 Typing run | `t1_run.tsv`, `t1_platforms.tsv` | setting; group × platform |
| T2 DRB3, DRB4, DRB5 | `t2_drb.tsv`, `t2_two_allele.tsv` | gene; expected copies |
| T3 Completeness | `t3_genes.tsv` | gene |
| T4 ToMMo comparison | `t4_harmonization.tsv`, `t4_panels.tsv`, `t4_outlier.tsv`, `t4_homozygosity.tsv`, `t4_confidence_vs_panel.tsv` | step; gene × panel; item; gene; confidence |
| T5 Marker QC | `t5_record.tsv`, `t5_genes.tsv`, `t5_floors.tsv`, `t5_completeness.tsv` | step; gene; cohort × class; gene |
| T6 P groups | `t6_split_counts.tsv`, `t6_split_groups.tsv`, `t6_dqa1.tsv` | cohort; P group |
| T7 Alleles (detail) | `t7_alleles.tsv` | allele column |
| T8 Residues (detail) | `t8_residues.tsv` | residue column |
| A1 Design | `a1_samples.tsv`, `a1_tests_by_gene.tsv`, `a1_methods.tsv` | comparison × cohort; gene; setting |
| A2 Thresholds | `a2_scans.tsv`, `a2_positions.tsv`, `a2_gec_check.tsv` | scan; class |
| A3 Q1 Association | `a3_scans.tsv` (all cases), `a3_top.tsv` | scan; marker |
| A4 Q2 Conditioning | `a4_steps.tsv` (all cases), `a4_after.tsv`, `assoc_facts.json` (B 114) | step; scan; cohort |
| A5 Q3 GWAS lead | `a5_leads.tsv`, `a5_link.tsv`, `a5_given_lead.tsv`, `a5_lead_given.tsv` | SNP; cohort; scan; conditioning set |
| A6 Inflation | `a6_lambda.tsv` | scan × conditioning |
| A7 Summary | `a7_v57.tsv`, `a7_genotypes.tsv`, `a7_residues.tsv`, `assoc_facts.json` | cohort; cohort × group; residue or allele |
| A8 DVT | `a8_strata.tsv`, `a3_scans.tsv` and `a4_steps.tsv` (DVT comparisons), `a8_nodvt_positions.tsv`, `a8_lrt_vs_score.tsv`, `a8_nodvt_dqb1.tsv`, `a8_57v_given_87y.tsv`, `a8_given57v.tsv`, `a8_spa_lost.tsv` | comparison × cohort; marker |
| A9 Reported alleles | `a9_published.tsv`, `a9_ours.tsv`, `a9_carriers.tsv`, `assoc_facts.json` | study × allele; comparison × cohort × allele |
| A10 Appendix | `a10_mirror.tsv`, `a10_concordance.tsv`, `a10_one_model_only.tsv`, `a10_platform.tsv`, `a10_platform_chi2.tsv`, `a10_pos57_pgroups.tsv` | scan; platform; P group |
| A11 Single markers (detail) | `a11_single_markers.tsv` | comparison × cohort × allele or residue column |
| A12 Positions (detail) | `a12_positions.tsv` | comparison × cohort × position |
| Sources and checks | `00.upstream/upstream.json`, `03.decks/deck_check.json`, `04.checks/recheck.json` | run; file; check |

Cell conventions: `—` missing, `–` not applicable, each with its reason on the sheet; shaded cells pass the
sheet's threshold (or, on T1, fall below the call-quality marks); P as the slides print it, numeric on the
detail sheets.

## Selected table columns

### `a11_single_markers.tsv` (comparison × cohort × column)

`comparison`, `cohort`, `class` (allele / residue), `marker`, `gene`, `position`, `residue`, `test_id` (the
column the test is reported under: a mirror pair shares one), `mirror_of`; `n_case`, `n_control`,
`called_*`, `copies_*`, `freq_*`, `mac`, `missing_*` (random-model samples); `r_or`, `r_l95`, `r_u95`,
`r_p`, `r_p_spa`, `r_errcode` (random model); `f_or`, `f_l95`, `f_u95`, `f_p`, `f_errcode`, `f_n_case`,
`f_n_control` (fixed model, unrelated samples); `gec`, `bonf` (the scan's thresholds); `r_sig_gec`,
`r_sig_bonf`, `f_sig_gec`; `n3_gec` (cohorts significant at GEC in this comparison).

### `a12_positions.tsv` (comparison × cohort × position)

`unit` (`DQB1:57`), `gene`, `position`, `reference` (most frequent residue), `categories` (the other residue
groups, `rare` = residues below 1% pooled), `n`, `cases`, `controls`, `df`, `chisq`, `p`, `threshold`
(Bonferroni over the positions tested), `sig_bonf`, `ors` (OR of each group against the reference),
`errcode` (`.`, `SEPARATION`, `NOT_TESTED`), `n3`.

### `a2_scans.tsv` (comparison × cohort × class)

`markers`, `polymorphic`, `mirror_pairs`, `tests`, `me_gec`, `p_gec`, `p_bonf`, `gec_blocks`, `me_share`,
`bonf_stricter`; the exact max-T benchmark `boot_*`: `draws`, `samples`, `cases`, `controls`, and for
`score` and `wald` the threshold `p_maxt_<s>` with its 95% interval (`_lo`, `_hi`), `me_maxt_<s>`, the FWER
of GEC and of Bonferroni with Wilson intervals, `replicates_with_na_<s>`, `beyond_*` (draws whose smallest P
passed).

### `t7_alleles.tsv` / `t8_residues.tsv` (column)

Per cohort with prefix `n_`, `i_`, `f_`: `status` (Tested, Below MAF 0.5%, Call rate below 95%, Gene not
tested, Not varying), `cr_cases`, `cr_controls`, `copies_cases`, `copies_controls`, `maf`. T7 adds
`p_group`, `first_release`, `match`, `matched_to`, `tommo_freq`, `kg_freq` (frequency of the P group in the
panel; missing when the P group is in neither our controls nor the panel) and `tommo_compared`,
`kg_compared`.

### `deck_check.json`, `recheck.json`, `audit.json`

`all_passed`; `checks`: per check `id`, `group`, `sheet`, `what`, `route`, `compared`, `mismatches` (audit:
`failures`), `passed`, `detail` (up to 50 mismatches: item, workbook value, deck or recomputed value).
`deck_check.json` adds `deck_files` (path and md5 of the deck data it read); `audit.json` adds `xlsx_md5`.
