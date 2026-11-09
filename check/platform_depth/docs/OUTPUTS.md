# Outputs — `check/platform_depth/results`

Final tables and figures are copies (not links into `work/`). `<label>` is the down-sampling run
(`pilot`, then `downsampling`); `<arm>` is `baseline` or `r01` … `rNN`; `<cohort>` is `narrow_mainland`,
`intermediate_mainland` or `full_mainland`. Which process writes which folder, and which gate or evidence-matrix
column each feeds: [PIPELINE.md](PIPELINE.md). Table conventions: [below](#table-conventions).

| Folder | File | Contents |
|---|---|---|
| `00.samples/` | `sample_table.tsv` | one row per genome of the three cohorts: status, `platform` (WGS_Platform verbatim), `platform_slug`, target and sheet depth, `deep` (down-sampled), resolved CRAM / index and their source, batch, surviving original GVCF folder; `in_*`, `fixed_*`, `rv_*` cohort and model sets |
| | `counts.tsv`, `sample_summary.json` | cohort × status × platform counts; deep-case and GVCF summary |
| `01.depth/` | `depth_per_sample.tsv` | effective depth (probe windows, duplicates excluded and kept), duplicate share, depth over the calling intervals, read length, the read length the sheet's `Observed_Depth` implies |
| | `platform_depth.tsv` | per cohort × platform: medians and IQR of the above (the regenerated read-length audit) |
| | `fractions.tsv`, `depth_summary.json` | the 89 deep cases' kept fraction f = D* / D_i; the target D* and its source |
| | `local_depth.tsv` | per calling interval: median depth by platform, the deep cases after scaling, ratio to the controls |
| | `depth_long.tsv.gz` | every genome × region measurement |
| `02.targets/` | `signal_variants.tsv` | the 168 GWAS signal variants and the genes' tested variants per cohort × variant set |
| | `hc.intervals`, `hc_intervals.tsv`, `extract.bed` | calling intervals (5-kb grid of the original intervals) and the read-extraction regions |
| | `fixed_sites.tsv`, `fixed_positions.intervals` | every variant re-genotyped (kind; membership per cohort fileset) |
| | `null_snps.tsv`, `targets_qc.json` | the null SNPs and indels by MAF bin (platform analysis); counts |
| `03.calls/<label>/` | `arm_plan.tsv`, `orig_gvcfs.tsv`, `deep_cases.tsv` | who runs in which arm, with fraction and seed; the original GVCFs joint-genotyped beside the arms |
| | `arms.qc.vcf.gz`, `per_chrom/` | every arm's genotypes at the fixed sites after the cohort's genotype QC (`FT` = the failed rule) |
| `04.genotypes/` | `original/<cohort>/gwas.*`, `gene.*` | the published genotypes at the fixed sites, rows copied byte for byte |
| | `<label>/<arm>/<cohort>/swap.*.json` | what each swap changed (untouched genomes identical) |
| `05.hla/` | `original/…`, `<label>/<arm>/06.assoc_inputs/`, `02.evidence/call_evidence.tsv` | hla.typing.v2's step-06 outputs and call table of each side run |
| | `<label>/<arm>/matrices/<cohort>/{allele,residue}_dosage.tsv` | the arm's matrices on the ORIGINAL tested columns |
| `06.retest/` | `original/<cohort>/`, `<label>/<arm>/<cohort>/` | `gwas.{best_guess,mean}.tsv` (β, SE = 1/√var, z, P no-SPA and SPA, allele frequencies), `gene.<set>.tsv` (P, OR, carriers, Bonferroni, BH), `assoc.random.<class>.tsv` (assoc_hla.v2 schema), `_omnibus/omnibus.residue.tsv` |
| `07.downsampling/<label>/` | `concordance/concordance.{summary,samples,replicates}.tsv`, `.long.tsv.gz` | genotype concordance by comparison, arm, sample, platform, variant class and MAF class |
| | `qc_reasons/qc_reasons.{summary,samples}.tsv` | share of genotypes each QC rule no-calls: arms, fidelity controls, the cohort's pre-QC genotypes by group (the summary counts the full-mainland members; per cohort: `summary/ds.qc_reasons.tsv`) |
| | `hla/hla.{calls_long,summary,read_support,dosage}.tsv` | HLA calls and tested dosages, arm vs archive and baseline; read support vs controls |
| | `summary/ds.*` | every signal across arms (GWAS leads and credible sets, genes, HLA markers, positions, DQB1:57:V); per cohort, the exposure of each signal in the deep cases (`ds.exposure.tsv`) and the share of genotypes each QC rule no-calls by group (`ds.qc_reasons.tsv`) |
| `08.subsets/<cohort>/` | `subsets/{pooled|only_<slug>|without_<slug>}.{random,fixed}.tsv`, `.keep`, `.rv_include.txt` | the platform subsets (all controls in each) |
| | `subsets/subset_counts.tsv`, `null/<subset>/null.{varianceRatio.txt,fit.log}` | counts per subset; the subset null fits (platform and control-source subsets) |
| | `control_sources/subsets/{ctrlonly,ctrlwithout}_<source>.*`, `subset_counts.tsv` | the control-source subsets: every case with one source's controls, or with all sources but one |
| `09.platform_assoc/` | `<cohort>/<subset>.saige.tsv`, `gene.<subset>.<set>.tsv`, `carriers_firth.tsv`; `_omnibus/omnibus.<subset>.tsv` | per-subset results |
| `10.platform/<cohort>/` | `platform.strata.tsv` | per marker × subset: β, SE, OR (95 % CI), P, counts |
| | `platform.heterogeneity.tsv` | per marker: Q, df, P (χ² and empirical), I² (95 % CI), FE and RE estimates, pooled, power per stratum |
| | `platform.lopo.tsv` | per marker × left-out platform: β, SE, P, z observed and expected, β_k − β_(−k) and its P |
| | `platform.case_only.tsv` | per marker: multinomial LRT within cases; genes: carriers by platform and permutation P |
| | `platform.rho.tsv`, `platform.null_calibration.tsv`, `platform.null_Q.tsv` | shared-control correlation (analytic vs null SNPs); per-subset λ and |z| > 3 on the null SNPs; null Q distribution |
| `11.array/` | `array.evidence.tsv` | per cohort × signal × marker × evidence: array marker, Minimac R², proxy r², genotype concordance with WGS in cases and controls (`conc_cases`, `conc_controls`; among non-reference genotypes `nonref_conc_*`), OR on WGS and on the array in the same people (per the lead's ALT allele, `or_allele`; DQB1:57:V per V), z on WGS and on the array (`z_wgs`, `z_array`), Δz, case frequency and r² in cases (HIBAG), verdict; gene drivers (full mainland only): WGS case carriers, confirmed on the imputed arrays, expected recovery and count (a diagnostic, not read), mean dosage, cases imputed as carriers without a WGS call, controls (panel members) (METHODS §7) |
| | `array.signals.tsv`, `array.summary.json`, `array.gate.json` | per cohort × signal: strongest evidence, verdict (consistent on arrays, contradicted, not checkable; genes: carriers counted, full mainland), case carriers recovered; per cohort the genome-wide numbers (people, concordance, λ, z correlation) and the verdict counts, the imputation, and the gene carrier totals; gate G10 |
| | `cohorts/<cohort>/{04.compare,05.rare_direct,06.imputed}/` | the per-cohort re-run of check/array_crosscheck in its own table layout (`array_cohorts.py`; counts and estimates only) |
| | `cohorts/<cohort>/06.imputed/concordance.tsv` | per imputed lead and driver × cases / controls: people called on both, best-guess genotype concordance with WGS, and among people non-reference on either side (not part of gate G10) |
| | `cohorts/<cohort>/07.hibag/summary.json` | DQB1:57:V imputed with HIBAG's published Asian model (`array_hibag.py`): DQB1 allele concordance, 57 V agreement, r² and frequencies in cases and controls, the typed-minus-imputed residual, and the association with the typed and the HIBAG dosage; the evidence row of DQB1:57:V reads it (the ridge prediction in 04.compare serves gate G10 only) |
| `12.controls/` | `controls.{counts,strata,heterogeneity,loo,null_lambda,genes,positions,signals}.tsv` | control-source sensitivity of every signal (METHODS §8): per source, across sources (Q, I²), leaving each out, calibration, genes, positions, the headline summary |
| `13.followup/` | `followup.gwas_af.tsv` | per cohort × GWAS lead × group (case platform × registration study, control source, ToMMo 61KJPN, 1KG JPT): ALT frequency (METHODS §9) |
| | `followup.gene_carriers.tsv`, `followup.gene_filters.tsv` | per cohort × tested variant and gene × group: carriers, possible carriers set missing, reference frequencies and the carriers ToMMo predicts; the re-applied genotype filters against the published genotypes |
| | `followup.hla_v57.tsv` | per cohort, every down-sampled case whose DQB1:57:V dosage differs between arms: dosage, DQB1 call, reads, HLA-HD confidence per arm |
| | `followup.lead_preqc.tsv`, `followup.lead_window.tsv` | per cohort × GWAS lead × group: pre-QC and published carriers, genotypes set missing (with ALT evidence); every pre-QC record within 50 bp of a lead with its carriers per group |
| `14.evidence/` | `evidence.long.tsv`, `evidence.loo.tsv`, `evidence.loo_calibration.tsv` | per signal × check × cohort: metric, value, reading (the gene sets' array check: one full-mainland row for the three cohorts); every leave-one-out z observed, expected and pooled, its conditional SD, P_diff and chance of a flag; per check × cohort the flags observed against those expected by chance (METHODS §10) |
| | `evidence.signals.tsv`, `evidence.rules.tsv` | one row per signal: discovery, the reading of every check over the three cohorts (one column per check of the registry, `scripts/checks.py`), the flags; the rule of every check |
| `figures/0N.<family>/` | `*.png`, `*.stats.json`, `README.md` | the figures, one README per family (`01.depth`, `02.downsampling`, `03.platform`, `04.array`, `05.controls`, `06.evidence`) |
| `_gates/` | `G1–G5*.json`, `<label>/G6,G8,G9.json`, `gates_G1-G4.tsv` | every gate and its numbers (G10: `11.array/array.gate.json`) |
| `_pilot/` | `R_decision.json` | the replicate count and the projection it rests on |
| `_run_info/<run_tag>/` | `trace.txt`, `report.html`, `timeline.html` | one record per run (never overwritten) |

## Table conventions

What holds in every report table (`07.downsampling/<label>/summary`, `10.platform`, `11.array`, `12.controls`,
`13.followup`, `14.evidence`; `verify.sh` §15 checks the first three points):

- Tab-separated, one header line of unique names, no index column; a missing value is an empty cell (never `NaN`,
  `NA` or `inf` text); numbers at full precision (rounding only in `metric` texts and figures).
- `cohort` takes the three cohort names; a table without it is single-cohort by its folder (`10.platform/<cohort>/`)
  or describes genomes, not a cohort (`07.downsampling/<label>/{concordance,hla}`: the deep cases belong to all
  three cohorts).
- Platforms by their `WGS_Platform` label in values (`DNBSeq-T7 30x`, … note `DNBseq-G400RS 15x`); slugs only in
  names: subsets `pooled`, `only_<slug>`, `without_<slug>` (`t7`, `novaseq`, `g400rs30`, `g400rs15`) and
  `ctrlonly_<source>`, `ctrlwithout_<source>` (`nagahama`, `bbj`, `htlv1`, `acc`, `lungca`); `left_out` holds the
  slug.
- Arms `original`, `baseline`, `r01` … ; wide tables append the arm (`p_r03`, `beta_original`); spread over the
  replicates `<quantity>_rep_{median,min,max}`; counts over replicates `reps_*`.
- Effects on the log-odds scale `beta` with `se` (SAIGE: 1/√var, METHODS §1), `or` with `or_l95` / `or_u95`; derived
  tests `p_<test>` (`p_Q`, `p_Q_empirical`, `p_diff`, `p_lrt`, `p_perm`).
- Variants `CHROM:POS:REF:ALT`; HLA markers as assoc_hla.v2 names them (`DQB1*05:01`, `DQB1:57:V`), positions
  `DQB1:57`; genes in their tiered sets `<gene> (<set>)` in `14.evidence`.
- Groups of `13.followup`: `case | <platform>` (GWAS leads also `case | <platform> | new study` / `old study`),
  `control | <source>`, `case | all`, `control | all`, `reference | ToMMo 61KJPN`, `reference | 1KG JPT`.
- Readings of `14.evidence`: the registry's vocabulary (`scripts/checks.py` `READINGS`; docs/PIPELINE.md).

Known exceptions (renaming them would change published columns; listed, not applied):
[REFACTOR_PROPOSAL.md](REFACTOR_PROPOSAL.md) §B.
