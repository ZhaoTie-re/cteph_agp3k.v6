# Pipeline — `check/platform_depth`

Written by `python3 scripts/checks.py markdown > docs/PIPELINE.md` from the registry in `scripts/checks.py` (do not edit by hand; `verify.sh` fails when it is out of date). Methods: [METHODS.md](METHODS.md); every file: [OUTPUTS.md](OUTPUTS.md). Paths are under `results/`; `<arm>` is `original` or `<label>/<arm>`.

## Stages

Run with `bash run.sh <stage>` (preflight, next run tag, detached launch, cache-miss guard; README).

| Stage | Arguments | Run tag | What |
|---|---|---|---|
| `inputs` | `--run inputs` | `inputs.N` | sample table, targets, depth and fractions (part of every run; cached after the first) |
| `gates` | `--run gates --run_label downsampling --Replicates 3` | `gates.N` | the ORIGINAL arm re-tested; gates G1-G4 |
| `pilot` | `--run gates,pilot,platform --run_label pilot --Replicates 3` | `pilot.N` | the pilot arms (costs for R) with the platform stage |
| `downsampling` | `--run gates,downsampling,platform --run_label downsampling --Replicates 3` | `full.N` | every down-sampling arm (R = 3; user, 2026-11-03) with the platform stage |
| `platform` | `--run gates,platform --run_label downsampling --Replicates 3` | `full.N` | the platform and control-source subsets |
| `report` | `--run report --run_label downsampling` | `report.N` | summaries, evidence matrix and figures from the published results (never cached but ARRAY_COHORTS) |
| `verify` | `bash verify.sh downsampling` | | read-only checks of every published result |

## Processes: inputs

| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |
|---|---|---|---|---|
| `SAMPLE_TABLE` | `sample_table.py` | info sheet; cohort pheno.tsv, fixed-model .fam; tuning.rv call set; CRAM paths | `00.samples/` | — |
| `BUILD_TARGETS` | `build_targets.py` | integration/saige_signals variants; assoc_gene tiers; original calling intervals; cohort .bim | `02.targets/` | — |
| `MAKE_CHUNKS` | (inline) | SAMPLE_TABLE | — | — |
| `DEPTH` | `depth_probe.py` | CRAMs; BUILD_TARGETS intervals | — | — |
| `DEPTH_SUMMARY` | `depth_summary.py` | DEPTH; SAMPLE_TABLE | `01.depth/` | — |

## Processes: gates

| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |
|---|---|---|---|---|
| `SUBSET_ORIGINAL` | `bedio.py` | random-model fileset; tuning.rv call set; BUILD_TARGETS fixed sites | `04.genotypes/original/<cohort>/` | — |
| `RETEST_GWAS` | `bedio.py`, `saige_collect.py` | an arm's GWAS fileset; assoc_saige null models | `06.retest/<arm>/<cohort>/gwas.*.tsv` | depth |
| `GATE_G1` | `saige_collect.py` | RETEST_GWAS (original, best_guess); assoc_saige saige.assoc.tsv | `_gates/` | G1 |
| `RETEST_GENE` | `gene_retest.py`, `bedio.py` | an arm's gene fileset; assoc_gene sets, pheno, scan | `06.retest/<arm>/<cohort>/gene.<set>.tsv`, `_gates/` | G2, depth |
| `HLA_CHAIN` | `hla_v2_wrapper.py` | hla.typing.v2 00.inputs, or an arm's re-typing (HLA_INPUTS) | `05.hla/<arm>/` | — |
| `GATE_G3` | `hla_v2_wrapper.py` | HLA_CHAIN (original); hla.typing.v2 06.assoc_inputs md5 | `_gates/` | G3 |
| `RETEST_HLA` | (assoc_hla.v2 scripts) | an arm's matrices; assoc_hla.v2 00.prep, 01.markers | `06.retest/<arm>/<cohort>/assoc.random.<class>.tsv` | depth |
| `RETEST_OMNIBUS` | (assoc_hla.v2 omnibus.py) | RETEST_HLA residues of the three cohorts | `06.retest/<arm>/_omnibus/` | depth |
| `GATE_G4` | `hla_gate.py` | RETEST_HLA, RETEST_OMNIBUS (original); assoc_hla.v2 published results | `_gates/` | G4 |
| `GATE_SUMMARY` | (inline) | the G1-G4 records | `_gates/gates_G1-G4.tsv` | G1-G4 |

## Processes: downsampling

The `pilot` and `downsampling` stages run these with their own `--run_label`; the gates-stage re-tests (RETEST_*, HLA_CHAIN) run on every arm.

| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |
|---|---|---|---|---|
| `T7_FRACTIONS` | `t7_fractions.py` | DEPTH_SUMMARY fractions and per-sample depth | `01.depth/fractions.all.tsv` | — |
| `ARM_PLAN` | `arm_plan.py` | T7_FRACTIONS (every down-sampled case); SAMPLE_TABLE | `03.calls/<label>/` | — |
| `EXTRACT_READS` | (samtools) | CRAMs; BUILD_TARGETS extraction regions | — | — |
| `SUBSAMPLE` | (samtools) | EXTRACT_READS; the arm's fraction and seed | — | G8 |
| `CALL_GVCF` | (GATK HaplotypeCaller) | SUBSAMPLE; calling intervals | — | — |
| `JOINT` | (GATK GenomicsDBImport, `GenotypeGVCFs)` | every arm's GVCF; the original GVCFs | — | — |
| `GENOTYPE_QC` | `genotype_qc.sh`, `expand_fixed.py` | JOINT; fixed sites | `03.calls/<label>/per_chrom/` | — |
| `CONCAT_QC` | (bcftools) | GENOTYPE_QC | `03.calls/<label>/arms.qc.vcf.gz` | — |
| `SWAP` | `bedio.py` | CONCAT_QC; SUBSET_ORIGINAL; ARM_PLAN deep cases | `04.genotypes/<label>/<arm>/<cohort>/swap.*.json` | G7 |
| `HLA_SUBSAMPLE` | (samtools) | hla.typing 01.reads; the arm's fraction and seed | — | — |
| `HLAHD` | (HLA-HD 1.7.1) | HLA_SUBSAMPLE | — | — |
| `HLA_INPUTS` | `hla_v2_wrapper.py` | HLAHD, HLA_SUBSAMPLE of one arm | — | — |
| `HLA_MATRICES` | `hla_matrices.py` | HLA_CHAIN exports of an arm | `05.hla/<label>/<arm>/matrices/` | — |
| `CONCORDANCE` | `concordance.py` | CONCAT_QC; SUBSET_ORIGINAL (full mainland); ARM_PLAN | `07.downsampling/<label>/concordance/` | G6 |
| `QC_REASONS` | `concordance.py` | CONCAT_QC; the cohort's pre-QC VCFs; ARM_PLAN | `07.downsampling/<label>/qc_reasons/` | — |
| `HLA_CONCORDANCE` | `hla_concordance.py` | HLA_CHAIN evidence; HLA_MATRICES; the archive | `07.downsampling/<label>/hla/` | G9 |
| `GATES_DS` | `gates_ds.py` | CONCORDANCE; SUBSAMPLE depth; HLA_CONCORDANCE | `_gates/<label>/` | G6, G8, G9 |

## Processes: platform

| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |
|---|---|---|---|---|
| `MAKE_SUBSETS` | `make_subsets.py` | SAMPLE_TABLE; cohort pheno, covariates; assoc_hla.v2 prep | `08.subsets/<cohort>/subsets/` | — |
| `MAKE_CTRL_SUBSETS` | `control_subsets.py` | MAKE_SUBSETS (pooled); info sheet Project | `08.subsets/<cohort>/control_sources/` | — |
| `FIT_NULL` | (SAIGE step 1) | a subset; random-model genotypes; assoc_saige GRM markers | `08.subsets/<cohort>/null/<subset>/` | G5 |
| `PLATFORM_GENO` | `bedio.py` | random-model fileset; fixed sites; null SNPs | — | — |
| `STEP2_SUBSET` | `bedio.py`, `saige_collect.py` | FIT_NULL; PLATFORM_GENO; RETEST_HLA (original) | `09.platform_assoc/<cohort>/<subset>.saige.tsv` | G5 |
| `PLATFORM_GENE` | `gene_retest.py`, `bedio.py` | subset include lists; SUBSET_ORIGINAL gene fileset | `09.platform_assoc/<cohort>/gene.<subset>.<set>.tsv` | lopo, sources_loo (genes) |
| `CARRIER_FIRTH` | `gene_retest.py`, `carrier_firth.py` | SUBSET_ORIGINAL gene fileset; MAKE_SUBSETS | `09.platform_assoc/<cohort>/carriers_firth.tsv` | — |
| `PLATFORM_OMNIBUS` | (assoc_hla.v2 omnibus.py) | RETEST_HLA (original) residues; a subset | `09.platform_assoc/_omnibus/` | lopo, sources_loo (position) |
| `PLATFORM_STATS` | `platform_stats.py` | STEP2_SUBSET (platform subsets); PLATFORM_GENO; CARRIER_FIRTH; null SNPs | `10.platform/<cohort>/` | platform_Q, platform_RE, lopo |
| `GATE_G5` | (inline) | FIT_NULL, STEP2_SUBSET (pooled); RETEST_GWAS (original) | `_gates/` | G5 |

## Processes: report

`cache false`: every report run re-runs them from the published results; editing their scripts never touches the cache.

| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |
|---|---|---|---|---|
| `DS_SUMMARY` | `ds_summary.py` | 06.retest; 04.genotypes/original/<cohort>; 03.calls deep cases; qc_reasons | `07.downsampling/<label>/summary/` | discovery, depth |
| `PLATFORM_REPORT` | `platform_report.py` | 10.platform; 08.subsets; 09.platform_assoc | `10.platform/_summary/` | — |
| `ARRAY_COHORTS` | `array_cohorts.py`, `array_hibag.py`, `array_hibag.R` | check/array_crosscheck inputs and results; each cohort's fixed-model members, phenotype, covariates; the QC'd and imputed arrays; the HIBAG model | `11.array/cohorts/` | G10 |
| `ARRAY_EVIDENCE` | `array_evidence.py` | 11.array/cohorts (ARRAY_COHORTS); integration/gene_signals LOVO (the gene drivers) | `11.array/array.*` | array |
| `CONTROL_REPORT` | `control_report.py` | 08.subsets; 09.platform_assoc | `12.controls/controls.*` | sources_Q, sources_loo |
| `HLA_REF_SHIFT` | `hla_reference_shift.py` | hla.typing.v2 05.exports; assoc_hla.v2 conditional scan; 06.retest/original | `12.controls/hla_reference_shift.*` | — |
| `FOLLOWUP` | `signal_followup.py` | 04.genotypes/original/<cohort>; the cohort pre-QC VCFs; 05.hla; integration frequencies | `13.followup/` | — |
| `EVIDENCE` | `evidence_matrix.py` | the tables above; 10.platform; 09.platform_assoc | `14.evidence/` | every check |
| `FIGURES` | `figures.py`, `catalogue.py` | every report table; 01.depth; 07.downsampling; 10.platform | `figures/` | — |

`choose_R.py` is run by hand after the pilot (`_pilot/R_decision.json`).

## Checks (evidence-matrix columns)

Readings: **consistent** the check agrees with the signal (arrays: same direction and |z WGS - z array| <= 3); **flag** the check disagrees; flags, never removes; **not testable** no left-out group has |z expected| >= 2.8; **not applicable** the check does not apply (Q and random effects at MAF < 5 %; per-stratum estimates of genes); **reported** a number without a rule (the residue position; the case carriers of a gene's drivers confirmed on the imputed arrays: one cell for the three cohorts, counted in full mainland, which contains the other two); **not checkable** arrays: not typed, not imputed with R2 >= 0.8, no LD proxy with r2 >= 0.8 (genes: no driver in the imputed arrays). Every check is read in each cohort and flags a signal when its rule holds in any one of the 3 (was two of three before 2026-11-04). Thresholds: leave-one-out |z expected| ≥ 2.8, |z observed| < 1.96; depth Δ|z| > 1; heterogeneity 0.05 / signals at MAF ≥ 0.05.

### `discovery` — Discovery

- **Artifact targeted:** none: the main analysis' own reading, shown beside the checks
- **Samples:** as published, three cohorts
- **Signals:** every signal
- **Readout:** P of the ORIGINAL arm (= the published result, gates G1-G4; `*_original` columns) against its threshold
- **Rule:** the main analysis: P against its pre-specified threshold (GWAS 5e-8, suggestive 1e-5; genes Bonferroni; HLA residues GEC; positions Bonferroni); a reading, never a flag
- **Tables:** `07.downsampling/<label>/summary/ds.*.tsv`
- **Figure:** —
- **METHODS:** §1

### `depth` — Read down-sampling

- **Artifact targeted:** read depth: every case platform but DNBseq-G400RS 15x against controls sequenced at 15x
- **Samples:** the 408 down-sampled cases at the controls' median depth (89 NovaSeq and DNBSeq-G400RS 30x cases each to it; 319 DNBSeq-T7 cases by one platform fraction), R = 3 replicates; three cohorts
- **Signals:** every signal
- **Readout:** replicate-median |z| against the original |z|; replicates significant
- **Rule:** read down-sampling of the down-sampled cases: flag if the replicate-median |z| is more than 1 below the original |z| (z from the two-sided P)
- **Tables:** `07.downsampling/<label>/summary/ds.{gwas,hla_dqb1_57v,positions,genes}.tsv`
- **Figure:** `02.downsampling/{gwas_leads,hla_dqb1_57v}`
- **METHODS:** §3, §5

### `platform_Q` — Platform heterogeneity

- **Artifact targeted:** a call artifact of one case platform
- **Samples:** each case platform's cases against all controls (shared); three cohorts
- **Signals:** GWAS leads, DQB1:57:V; genes not applicable
- **Readout:** Q across the four case platforms, its P against null SNPs of the same MAF bin, I²
- **Rule:** heterogeneity across the four case platforms (shared controls): flag if the empirical P of Q < 0.05 / number of SAIGE signals; applied at MAF >= 5% only
- **Tables:** `10.platform/<cohort>/platform.heterogeneity.tsv`
- **Figure:** `03.platform/forest.<cohort>`
- **METHODS:** §6a

### `platform_RE` — Random effects across platforms

- **Artifact targeted:** an effect carried by some case platforms only
- **Samples:** as platform_Q
- **Signals:** GWAS leads, DQB1:57:V; genes not applicable
- **Readout:** random-effects OR (REML, Hartung-Knapp) against the pooled OR
- **Rule:** random-effects estimate (REML, Hartung-Knapp) on the other side of the null from the pooled one; applied at MAF >= 5% only (as Q: SAIGE's one-step stratum beta is unreliable in small strata)
- **Tables:** `10.platform/<cohort>/platform.heterogeneity.tsv`
- **Figure:** `03.platform/forest.<cohort>`
- **METHODS:** §6b

### `lopo` — Leave one platform out

- **Artifact targeted:** a signal carried by one case platform
- **Samples:** all cases but one platform's, all controls; three cohorts
- **Signals:** GWAS leads, DQB1:57:V, genes (CMC burden β); the position DQB1 57 reported
- **Readout:** z observed without each platform against z expected = pooled β / SE without it
- **Rule:** leave one case platform out: flag if a left-out platform has |z expected| >= 2.8 (z expected = pooled beta / SE without it) and |z observed| < 1.96 or the opposite sign; genes: the CMC burden beta
- **Tables:** `10.platform/<cohort>/platform.lopo.tsv`; `09.platform_assoc/<cohort>/gene.without_*.tsv`; `14.evidence/evidence.loo.tsv`
- **Figure:** `06.evidence/z_observed_expected` (a-c)
- **METHODS:** §6c, §10

### `array` — SNP arrays

- **Artifact targeted:** a WGS calling error (a second technology in the same people)
- **Samples:** the unrelated genomes with arrays, DNBSeq-T7 cases and HiSeq X controls: narrow 1,398, intermediate 1,460, full 1,683; three cohorts
- **Signals:** GWAS leads, DQB1:57:V; genes described (one cell)
- **Readout:** one marker per signal: the lead typed, else imputed (R2 >= 0.8; DQB1:57:V with HIBAG), else an LD proxy (r2 >= 0.8); its association re-run in the same people against WGS; check/array_crosscheck re-run per cohort (gate G10). Genes: the WGS case carriers of each driver confirmed on the imputed arrays (dosage >= 0.5), counted once in full mainland, which contains the other two cohorts; a description, no rule
- **Rule:** SNP arrays (DNBSeq-T7 cases, HiSeq X controls): flag if the opposite direction or |z WGS - z array| > 3, in any cohort; the array P is not required; genes reported, never flagged
- **Tables:** `11.array/array.{signals,evidence}.tsv`; `11.array/cohorts/<cohort>/`
- **Figure:** `04.array/array_evidence`
- **METHODS:** §7

### `sources_Q` — Control-source heterogeneity

- **Artifact targeted:** a control-source batch or population difference
- **Samples:** every case against each control source's controls (5 sources); three cohorts
- **Signals:** GWAS leads, DQB1:57:V; genes not applicable
- **Readout:** Q across the five sources (shared cases), P, I²
- **Rule:** heterogeneity across the five control sources (shared cases): flag if the analytic (chi-square) P of Q < 0.05 / number of SAIGE signals; applied at MAF >= 5% only
- **Tables:** `12.controls/controls.heterogeneity.tsv`
- **Figure:** `05.controls/sources.<cohort>`
- **METHODS:** §8

### `sources_loo` — Leave one control source out

- **Artifact targeted:** a signal carried by one control source
- **Samples:** every case against all sources but one; three cohorts
- **Signals:** GWAS leads, DQB1:57:V, genes (CMC burden β); the position DQB1 57 reported
- **Readout:** z observed without each source against z expected
- **Rule:** leave one control source out: as lopo
- **Tables:** `12.controls/controls.{loo,genes,positions}.tsv`; `14.evidence/evidence.loo.tsv`
- **Figure:** `06.evidence/z_observed_expected` (d-f)
- **METHODS:** §8, §10

## Gates

| Gate | What | Rule | Process | Record |
|---|---|---|---|---|
| G1 | SAIGE step 2 on the ORIGINAL genotypes reproduces assoc_saige's saige.assoc.tsv | every row present; relative difference <= 1e-6 in every compared column | `GATE_G1` | `_gates/G1.<cohort>.json` |
| G2 | rvtest on the ORIGINAL genotypes reproduces the published gene P, OR and carriers | no mismatching row | `RETEST_GENE` | `_gates/G2.<cohort>.<set>.json` |
| G3 | the hla.typing.v2 wrapper reproduces the six 06.assoc_inputs matrices | md5 identical | `GATE_G3` | `_gates/G3.json` |
| G4 | the HLA re-test reproduces assoc_all.tsv and the residue-position omnibus | the same rows; relative difference <= 1e-6 | `GATE_G4` | `_gates/G4.json` |
| G5 | the pooled null re-fit reproduces the published null model | variance ratio and median z within 1 % | `GATE_G5` | `_gates/G5.<cohort>.json` |
| G6 | re-processing fidelity: baseline, re-genotyped original GVCFs and original agree | concordance >= 99.5 %, call rate within 0.5 percentage points | `GATES_DS` | `_gates/<label>/G6.json` |
| G7 | the swap leaves every genome but the deep cases' byte-identical | bedio.py swap aborts otherwise | `SWAP` | `04.genotypes/<label>/<arm>/<cohort>/swap.{gwas,gene}.json` |
| G8 | achieved depth: depth after / (f x depth before) over the calling intervals | median 0.95-1.05, every genome 0.85-1.15 | `GATES_DS` | `_gates/<label>/G8.json` |
| G9 | HLA re-typing at f = 1 reproduces the archived typing | every call identical | `GATES_DS` | `_gates/<label>/G9.json` |
| G10 | the array comparison re-run in full mainland reproduces check/array_crosscheck's tables | same rows and text; numbers within 1e-6 relative | `ARRAY_COHORTS` | `11.array/array.gate.json` |
