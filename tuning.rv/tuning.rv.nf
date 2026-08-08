nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════════════════════
// tuning.rv — rare-variant depth/platform-confounding QC + minor-allele-count
// sweep, across nested cohorts, with a cross-cohort comparison.
//
// The CTEPH/AGP3K v6 design is FULLY confounded (cases 30x DNBSeq/NovaSeq,
// controls 15x HiSeqX; zero platform overlap). This pipeline does not confirm
// associations — it QUANTIFIES and MITIGATES the confounding one layer at a
// time and recommends a null-calibrated, power-preserving minAC threshold.
//
// One `nextflow run` tunes ALL cohorts (default narrow/intermediate/full_mainland).
// Per cohort — SAMPLE QC BEFORE VARIANT QC (variant metrics are recomputed on
// the sample-cleaned set):
//
//   CALC_METRICS_BASE (minAC 0, RAW rare set; permanent storeDir cache)
//        └─ DETECT_OUTLIER_SAMPLES  (rate ~ MeanDP residual robust-Z) -> remove
//              ├─ PLOT_RATE_DISTRIBUTION (before/after rate distribution)
//              └─ FILTER_SAMPLES (--remove) -> sample-cleaned set
//                   └─ CALC_VARIANT_METRICS (VMISS_15x/30x on cleaned set)
//                        └─ DETECT_DEPTHDIFF_VARIANTS (Kneedle |ΔVMISS| CDF) -> exclude
//                             └─ FILTER_VARIANTS (--exclude) -> fully-filtered set
//                                  └─ CALC_METRICS (minAC 0..Max)
//                                       └─ QC_SUMMARY (rate + OLS, no-PC & +bbj_mainland-PC)
//                                            └─ MERGE_QC_SUMMARIES
//                                                 ├─ SELECT_MINAC (recommend)
//                                                 └─ PLOT_SWEEP (headline; marks rec. minAC)
//                                                      └─ WRITE_RUN_MANIFEST
//
//   ...then across cohorts:  COMPARE_COHORTS (overlay sweeps + summary table).
//
// Outputs live under results/<cohort>/…, the cross-cohort comparison under
// results/_comparison/, and Nextflow run reports under results/_run_info/.
// Full methodology: docs/METHODS.md ; output dictionary: docs/OUTPUTS.md.
// ═══════════════════════════════════════════════════════════════════════════

// ── Cohort selection ────────────────────────────────────────────────────────
// Default: tune all three. Override a subset with e.g. --Cohorts narrow_mainland,full_mainland
params.Cohorts      = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.ValidCohorts = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']

// ── Inputs / paths ───────────────────────────────────────────────────────────
params.GTPath     = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results'
// Rare (MAF < threshold) fixed-model call-set prefix, per cohort. @@COHORT@@ is
// substituted per cohort (.replace) — DRY template so the long path lives once.
params.BedTemplate = "${params.GTPath}/12_model_inputs/@@COHORT@@/genotype/fixed_model/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold"
params.InfoPath    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx'
// Authoritative sequencing-platform table (for the by-platform QC panels).
params.PlatformFile  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tmp/cram.v6/cram.v6.summary.csv'
params.PlatformIdCol = 'ID_JHRPv6'
params.PlatformCol   = 'WGS_Platform'
params.ResultsDir  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/results'
params.ScriptDir   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/scripts'

// Ancestry PC covariates for the +PC sweep model (bbj_mainland; per cohort).
// geno-PC is deliberately NOT used — it captures platform in this design.
params.PcTemplate = "${params.GTPath}/12_model_inputs/@@COHORT@@/covariates/bbj_mainland_pc.sex.tsv"
params.PcLabel    = 'bbj_mainland'
params.NPcs       = 10

// Conda env activated in every process; plink2 / tabix / bgzip / python come
// from PATH + this env (see the prelude at the top of each process script).
params.conda_env_activate = 'cteph_geno_pro'

// ── Info-file column mapping ──────────────────────────────────────────────────
params.IDCol     = 'ID_JHRPv6'
params.GroupCol  = 'Outcome'
params.CaseValue = 'PH'
params.TdpCol    = 'Target_Depth'
params.MdpCol    = 'Observed_Depth'

// ── Detection / model thresholds (frequently tuned) ──────────────────────────
params.RobustZThreshold = 5.0   // |robust z| cutoff for residual sample outliers
params.KneeS            = 1.0   // Kneedle sensitivity S
params.KneeWeightX      = 1.0   // Kneedle weight_x
params.KneeWeightY      = 2.0   // Kneedle weight_y
params.Robust           = false // OLS SE: classical (default) or HC3 robust (--robust)
params.Alpha            = 0.05  // non-significance level for the minAC V-trough

// ── Protected variants (candidate signals) — NEVER excluded by variant QC ─────
// These VariantIDs (bim/metrics format chr<N>:<pos>:<ref>:<alt>) are dropped from
// the depth-diff --exclude list even if their |ΔVMISS| exceeds the knee, so they
// always survive into the filtered call set. STBD1-region CTEPH candidates.
params.ProtectVariants = ['chr4:76306990:G:C', 'chr4:76309487:C:G', 'chr4:76309559:G:T']
params.MinACFloor       = 1     // smallest minAC eligible to be recommended
// PRIMARY model for the minAC recommendation. beta_group is V-SHAPED in minAC;
// the recommendation is the trough = smallest minAC where beta_group is non-sig
// (0 in 95% CI). Here beta_group is a DETECTOR for the technical artifact, not an
// effect estimate: with zero platform overlap the artifact IS the case/control
// contrast, and the ancestry PCs span 0.9-2.8% of that contrast, so adjusting on
// them moves the reading by 1.9-2.6x its own size and in cohort-dependent
// directions. 'nopc' is therefore primary; +PC is reported as a SENSITIVITY check
// (see docs/METHODS.md and results/_comparison/decision_axis.png). Downstream
// rvtest still uses PCs — it needs an estimator, not a detector.
params.SelPrefer        = 'nopc'  // 'nopc' (default) | 'pc'
// decision_axis.png illustrates the model choice with ONE cohort — the argument is
// about which statistic to read, not about cross-cohort spread, and overlaying all
// three forced cohort colours into panels that already encode the model. The
// per-cohort numbers are still written to decision_axis.tsv.
params.DecisionAxisCohort = 'intermediate_mainland'

// ── Fast threshold-tuning loop ────────────────────────────────────────────────
// Base metrics (CALC_METRICS_BASE) are threshold-INDEPENDENT and kept in a
// permanent storeDir (00.precheck/metrics) — computed once, never recomputed on
// re-runs (survives `work/` wipes, no -resume needed). --PrecheckOnly true stops
// after DETECTION (sample QC + variant QC), skipping the sweep + comparison.
params.PrecheckOnly = false

// ── Sweep range / resources ───────────────────────────────────────────────────
params.MaxAC       = 20
params.ThreadsCalc = 4

// Shared script prelude (PATH + conda) — a param, not a bare top-level var, so
// the strict DSL parser accepts it.
params.PRELUDE = """
    set -euo pipefail
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}
""".stripIndent()

// ═══════════════════════════════════════════════════════════════════════════
// CALC_METRICS_BASE — sample/variant metrics on the RAW rare set (minAC 0).
// Permanent storeDir cache (threshold-independent). Feeds sample detection.
// ═══════════════════════════════════════════════════════════════════════════
process CALC_METRICS_BASE {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    storeDir { "${params.ResultsDir}/${cohort}/00.precheck/metrics" }

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path("sample_metrics.txt.gz"), emit: sample_metrics
    tuple val(cohort), path("variant_metrics.txt.gz"), path("variant_metrics.txt.gz.tbi"), emit: variant_metrics
    tuple val(cohort), path("calc_metrics.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    bash ${params.ScriptDir}/run_calc_metrics.sh \\
        --bed-prefix "${bed.baseName}" \\
        --info-file "${params.InfoPath}" \\
        --out-sample "sample_metrics.txt" \\
        --out-variant "variant_metrics.txt" \\
        --script-dir "${params.ScriptDir}" \\
        --plink2 "\$(command -v plink2)" \\
        --tabix "\$(command -v tabix)" \\
        --id-col "${params.IDCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --tdp-col "${params.TdpCol}" \\
        --mdp-col "${params.MdpCol}" \\
        --min-ac 0 \\
        --threads ${params.ThreadsCalc}
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// DETECT_OUTLIER_SAMPLES — depth-confounded sample outliers (rate residual Z).
// ═══════════════════════════════════════════════════════════════════════════
process DETECT_OUTLIER_SAMPLES {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/01_sample_qc" }, mode: 'copy'
    publishDir { "${params.ResultsDir}/${cohort}/figures" }, mode: 'copy', pattern: "*.png"

    input:
    tuple val(cohort), path(sample_metrics)

    output:
    tuple val(cohort), path("sample_outliers.remove.txt"), emit: remove
    tuple val(cohort), path("sample_qc_rate_residual.png"), emit: pdf
    tuple val(cohort), path("detect_outlier_samples.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "DETECT_OUTLIER_SAMPLES[${cohort}]" "detect_outlier_samples.log"
    rv_kv "cohort" "${cohort}"
    rv_kv "z_threshold" "${params.RobustZThreshold}"

    python3 ${params.ScriptDir}/detect_outlier_samples.py \\
        --sample-metrics ${sample_metrics} \\
        --z-threshold ${params.RobustZThreshold} \\
        --platform-file "${params.PlatformFile}" \\
        --platform-id-col "${params.PlatformIdCol}" \\
        --platform-col "${params.PlatformCol}" \\
        --out-remove sample_outliers.remove.txt \\
        --out-png sample_qc_rate_residual.png 2>&1 | tee -a detect_outlier_samples.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// PLOT_RATE_DISTRIBUTION — burden-rate distribution before/after sample QC and
// its structure by depth / platform / group (the confounding, visualised).
// ═══════════════════════════════════════════════════════════════════════════
process PLOT_RATE_DISTRIBUTION {
    executor 'slurm'
    queue 'gr10478b'
    time '20m'
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/01_sample_qc" }, mode: 'copy'
    publishDir { "${params.ResultsDir}/${cohort}/figures" }, mode: 'copy', pattern: "*.png"

    input:
    tuple val(cohort), path(sample_metrics), path(sample_remove)

    output:
    tuple val(cohort), path("rate_distribution.png"), emit: pdf
    tuple val(cohort), path("plot_rate_distribution.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "PLOT_RATE_DISTRIBUTION[${cohort}]" "plot_rate_distribution.log"
    rv_kv "cohort" "${cohort}"

    python3 ${params.ScriptDir}/plot_rate_distribution.py \\
        --sample-metrics ${sample_metrics} \\
        --removed-file ${sample_remove} \\
        --platform-file "${params.PlatformFile}" \\
        --platform-id-col "${params.PlatformIdCol}" \\
        --platform-col "${params.PlatformCol}" \\
        --out-png rate_distribution.png 2>&1 | tee -a plot_rate_distribution.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// FILTER_SAMPLES — apply the sample outlier list to the RAW rare set (--remove).
// ═══════════════════════════════════════════════════════════════════════════
process FILTER_SAMPLES {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir { "${params.ResultsDir}/${cohort}/01.sample_filter" }, mode: 'symlink', pattern: "sample_cleaned.*"
    publishDir { "${params.ResultsDir}/${cohort}/01.sample_filter" }, mode: 'copy',    pattern: "filter_samples.log"

    input:
    tuple val(cohort), path(sample_remove)

    output:
    tuple val(cohort), path("sample_cleaned.bed"), path("sample_cleaned.bim"), path("sample_cleaned.fam"), emit: cleaned
    tuple val(cohort), path("filter_samples.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "FILTER_SAMPLES[${cohort}]" "filter_samples.log"
    BEDPREFIX="${params.BedTemplate.replace('@@COHORT@@', cohort)}"
    rv_kv "Input prefix" "\${BEDPREFIX}"
    rv_kv "Sample list" "${sample_remove}"
    rv_kv "Samples to remove" "\$(wc -l < "${sample_remove}")"

    plink2 --bfile "\${BEDPREFIX}" \\
        --remove "${sample_remove}" \\
        --make-bed --out sample_cleaned --threads ${params.ThreadsCalc} \\
        >> filter_samples.log 2>&1

    rv_kv "Samples retained" "\$(wc -l < sample_cleaned.fam)"
    rv_kv "Variants (unchanged)" "\$(wc -l < sample_cleaned.bim)"
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// CALC_VARIANT_METRICS — recompute depth-stratified variant metrics (VMISS_15X
// / VMISS_30X) on the SAMPLE-CLEANED set (sample QC before variant QC).
// ═══════════════════════════════════════════════════════════════════════════
process CALC_VARIANT_METRICS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'symlink', pattern: "variant_metrics.txt.gz*"
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'copy',    pattern: "calc_variant_metrics.log"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path("variant_metrics.txt.gz"), path("variant_metrics.txt.gz.tbi"), emit: variant_metrics
    tuple val(cohort), path("calc_variant_metrics.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    bash ${params.ScriptDir}/run_calc_metrics.sh \\
        --bed-prefix "${bed.baseName}" \\
        --info-file "${params.InfoPath}" \\
        --out-sample "sample_metrics.txt" \\
        --out-variant "variant_metrics.txt" \\
        --script-dir "${params.ScriptDir}" \\
        --plink2 "\$(command -v plink2)" \\
        --tabix "\$(command -v tabix)" \\
        --id-col "${params.IDCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --tdp-col "${params.TdpCol}" \\
        --mdp-col "${params.MdpCol}" \\
        --min-ac 0 \\
        --threads ${params.ThreadsCalc} 2>&1 | tee calc_variant_metrics.log
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// DETECT_DEPTHDIFF_VARIANTS — depth-differential variants (Kneedle |ΔVMISS|),
// computed on the sample-cleaned variant metrics.
// ═══════════════════════════════════════════════════════════════════════════
process DETECT_DEPTHDIFF_VARIANTS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'copy', pattern: "*.txt"
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'copy', pattern: "*.png"
    publishDir { "${params.ResultsDir}/${cohort}/figures" }, mode: 'copy', pattern: "*.png"
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'copy', pattern: "*.tsv"
    publishDir { "${params.ResultsDir}/${cohort}/00.precheck/02_variant_qc" }, mode: 'copy', pattern: "*.log"

    input:
    tuple val(cohort), path(variant_metrics), path(variant_metrics_tbi)

    output:
    tuple val(cohort), path("variants_depthdiff.exclude.txt"), emit: exclude
    tuple val(cohort), path("protect_variants.txt"), emit: protect
    tuple val(cohort), path("depthdiff_variants_cdf.png"), emit: pdf
    tuple val(cohort), path("depthdiff_stats.tsv"), emit: stats
    tuple val(cohort), path("detect_depthdiff_variants.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "DETECT_DEPTHDIFF_VARIANTS[${cohort}]" "detect_depthdiff_variants.log"
    rv_kv "cohort" "${cohort}"
    rv_kv "kneedle (S,wx,wy)" "${params.KneeS}, ${params.KneeWeightX}, ${params.KneeWeightY}"

    # Protected variants (candidate signals) — never excluded, even if |ΔVMISS| > knee.
    printf '%s\\n' ${params.ProtectVariants.collect { "'${it}'" }.join(' ')} > protect_variants.txt
    rv_kv "protected variants" "\$(grep -cve '^[[:space:]]*\$' protect_variants.txt)"

    python3 ${params.ScriptDir}/detect_depthdiff_variants.py \\
        --variant-metrics ${variant_metrics} \\
        --knee-s ${params.KneeS} \\
        --knee-weight-x ${params.KneeWeightX} \\
        --knee-weight-y ${params.KneeWeightY} \\
        --protect-file protect_variants.txt \\
        --out-exclude variants_depthdiff.exclude.txt \\
        --out-png depthdiff_variants_cdf.png \\
        --out-stats depthdiff_stats.tsv 2>&1 | tee -a detect_depthdiff_variants.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// FILTER_VARIANTS — apply the depth-diff exclude list to the sample-cleaned set.
// ═══════════════════════════════════════════════════════════════════════════
process FILTER_VARIANTS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir { "${params.ResultsDir}/${cohort}/02.callset_filter" }, mode: 'symlink', pattern: "filtered.*"
    publishDir { "${params.ResultsDir}/${cohort}/02.callset_filter" }, mode: 'copy',    pattern: "filter_variants.log"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), path(variant_exclude)

    output:
    tuple val(cohort), path("filtered.bed"), path("filtered.bim"), path("filtered.fam"), emit: filtered
    tuple val(cohort), path("filter_variants.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "FILTER_VARIANTS[${cohort}]" "filter_variants.log"
    rv_kv "Input prefix" "${bed.baseName}"
    rv_kv "Variant list" "${variant_exclude}"
    rv_kv "Variants to exclude" "\$(wc -l < "${variant_exclude}")"

    plink2 --bfile "${bed.baseName}" \\
        --exclude "${variant_exclude}" \\
        --make-bed --out filtered --threads ${params.ThreadsCalc} \\
        >> filter_variants.log 2>&1

    rv_kv "Samples retained" "\$(wc -l < filtered.fam)"
    rv_kv "Variants retained" "\$(wc -l < filtered.bim)"
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// CALC_METRICS — per-minAC sample/variant metrics on the FULLY-FILTERED set.
// ═══════════════════════════════════════════════════════════════════════════
process CALC_METRICS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag { "${cohort}:minac${min_ac}" }
    publishDir { "${params.ResultsDir}/${cohort}/03.qc_metrics/minac${min_ac}" }, mode: 'symlink', pattern: "*_metrics.txt.gz*"
    publishDir { "${params.ResultsDir}/${cohort}/03.qc_metrics/minac${min_ac}" }, mode: 'copy',    pattern: "calc_metrics.log"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), val(min_ac)

    output:
    tuple val(cohort), val(min_ac), path("sample_metrics.txt.gz"), path("variant_metrics.txt.gz"), path("variant_metrics.txt.gz.tbi"), emit: metrics
    tuple val(cohort), path("calc_metrics.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    bash ${params.ScriptDir}/run_calc_metrics.sh \\
        --bed-prefix "${bed.baseName}" \\
        --info-file "${params.InfoPath}" \\
        --out-sample "sample_metrics.txt" \\
        --out-variant "variant_metrics.txt" \\
        --script-dir "${params.ScriptDir}" \\
        --plink2 "\$(command -v plink2)" \\
        --tabix "\$(command -v tabix)" \\
        --id-col "${params.IDCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --tdp-col "${params.TdpCol}" \\
        --mdp-col "${params.MdpCol}" \\
        --min-ac ${min_ac} \\
        --threads ${params.ThreadsCalc}
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// QC_SUMMARY — per-minAC depth-confounding stats (rate + OLS, no-PC & +PC).
// ═══════════════════════════════════════════════════════════════════════════
process QC_SUMMARY {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag { "${cohort}:minac${min_ac}" }
    publishDir { "${params.ResultsDir}/${cohort}/04.qc_summary/minac${min_ac}" }, mode: 'copy'

    input:
    tuple val(cohort), val(min_ac), path(sample_metrics), path(variant_metrics), path(variant_metrics_tbi), path(pc_file)

    output:
    tuple val(cohort), val(min_ac), path("qc_summary_stats.minac${min_ac}.tsv"), emit: summary_stats
    tuple val(cohort), path("qc_summary.minac${min_ac}.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "QC_SUMMARY[${cohort}:minac${min_ac}]" "qc_summary.minac${min_ac}.log"
    rv_kv "cohort" "${cohort}"
    rv_kv "MinAC" "${min_ac}"
    rv_kv "PC source" "${params.PcLabel}"
    rv_kv "SE type" "${params.Robust ? 'HC3' : 'classical'}"

    python3 ${params.ScriptDir}/qc_summary.py \\
        --sample-metrics ${sample_metrics} \\
        --variant-metrics ${variant_metrics} \\
        --min-ac ${min_ac} \\
        --pc-file ${pc_file} \\
        --pc-label ${params.PcLabel} \\
        --n-pcs ${params.NPcs} \\
        ${params.Robust ? '--robust' : ''} \\
        --out-tsv qc_summary_stats.minac${min_ac}.tsv 2>&1 | tee -a qc_summary.minac${min_ac}.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// MERGE_QC_SUMMARIES — concatenate per-minAC rows; tag with a Cohort column.
// ═══════════════════════════════════════════════════════════════════════════
process MERGE_QC_SUMMARIES {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir { "${params.ResultsDir}/${cohort}/05.qc_collect" }, mode: 'copy'

    input:
    tuple val(cohort), path(tsvs)

    output:
    tuple val(cohort), path("all_qc_summary_stats.tsv"), emit: merged_stats
    tuple val(cohort), path("merge_qc_summaries.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "MERGE_QC_SUMMARIES[${cohort}]" "merge_qc_summaries.log"
    rv_kv "input_tables" "\$(echo ${tsvs} | wc -w)"

    awk 'NR==1 || FNR > 1' ${tsvs} > temp_merged.tsv
    python3 - <<'PY' 2>&1 | tee -a merge_qc_summaries.log
import pandas as pd

df = pd.read_csv('temp_merged.tsv', sep='\\t')
df = df.sort_values('MinAC_Threshold')
df.insert(0, 'Cohort', '${cohort}')
df.to_csv('all_qc_summary_stats.tsv', sep='\\t', index=False)
print(f'merged rows: {len(df)}  (minAC {int(df.MinAC_Threshold.min())}..{int(df.MinAC_Threshold.max())})')
PY
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// SELECT_MINAC — recommend the smallest null-calibrated, power-preserving minAC.
// ═══════════════════════════════════════════════════════════════════════════
process SELECT_MINAC {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir { "${params.ResultsDir}/${cohort}/05.qc_collect" }, mode: 'copy'

    input:
    tuple val(cohort), path(merged_stats)

    output:
    tuple val(cohort), path("minac_recommendation.tsv"), emit: rec
    tuple val(cohort), path("select_minac.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "SELECT_MINAC[${cohort}]" "select_minac.log"

    python3 ${params.ScriptDir}/select_minac.py \\
        --qc-stats ${merged_stats} \\
        --alpha ${params.Alpha} \\
        --min-floor ${params.MinACFloor} \\
        --prefer ${params.SelPrefer} \\
        --out-tsv minac_recommendation.tsv 2>&1 | tee -a select_minac.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// PLOT_SWEEP — headline minAC-sweep figure (marks the recommended minAC).
// ═══════════════════════════════════════════════════════════════════════════
process PLOT_SWEEP {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir { "${params.ResultsDir}/${cohort}/05.qc_collect" }, mode: 'copy'
    publishDir { "${params.ResultsDir}/${cohort}/figures" }, mode: 'copy', pattern: "*.png"

    input:
    tuple val(cohort), path(qc_stats), path(rec_tsv)

    output:
    tuple val(cohort), path("minac_sweep.png"), emit: pdf
    tuple val(cohort), path("plot_sweep.log"), emit: log

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "PLOT_SWEEP[${cohort}]" "plot_sweep.log"
    REC=\$(awk -F'\\t' '\$1=="recommended_minac"{print \$2}' ${rec_tsv})
    rv_kv "recommended_minac" "\${REC}"

    python3 ${params.ScriptDir}/plot_sweep.py \\
        --qc-stats ${qc_stats} \\
        --recommended-minac "\${REC}" \\
        --alpha ${params.Alpha} \\
        --out-png minac_sweep.png 2>&1 | tee -a plot_sweep.log
    rv_done
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// WRITE_RUN_MANIFEST — per-cohort durable README + machine-readable summary.
// ═══════════════════════════════════════════════════════════════════════════
process WRITE_RUN_MANIFEST {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir { "${params.ResultsDir}/${cohort}" }, mode: 'copy', pattern: 'README.md'
    publishDir { "${params.ResultsDir}/${cohort}/05.qc_collect" }, mode: 'copy', pattern: 'cohort_summary.tsv'

    input:
    tuple val(cohort), path(samp_log), path(var_log), path(depthdiff_stats), path(rec_tsv), path(sweep_pdf)
    val run_info

    output:
    tuple val(cohort), path("README.md"), emit: readme
    tuple val(cohort), path("cohort_summary.tsv"), emit: summary

    script:
    """
    ${params.PRELUDE}

    # in→out counts (filter logs) + knee (depthdiff sidecar) + recommendation (select_minac)
    SM_RM=\$(grep -m1 'Samples to remove'   ${samp_log} | sed 's/.*: *//')
    SR=\$(grep -m1 'Samples retained'       ${var_log}  | sed 's/.*: *//')
    VE=\$(grep -m1 'Variants to exclude'    ${var_log}  | sed 's/.*: *//')
    VR=\$(grep -m1 'Variants retained'      ${var_log}  | sed 's/.*: *//')
    KNEE=\$(awk -F'\\t' '\$1=="knee"{print \$2}'          ${depthdiff_stats})
    RETP=\$(awk -F'\\t' '\$1=="retained_pct"{print \$2}'  ${depthdiff_stats})
    PROT_PRES=\$(awk -F'\\t' '\$1=="n_protected_present"{print \$2}' ${depthdiff_stats})
    PROT_RESC=\$(awk -F'\\t' '\$1=="n_protected_rescued"{print \$2}' ${depthdiff_stats})
    REC=\$(awk -F'\\t'  '\$1=="recommended_minac"{print \$2}'  ${rec_tsv})
    CAL=\$(awk -F'\\t'  '\$1=="calibrated"{print \$2}'         ${rec_tsv})
    PMODEL=\$(awk -F'\\t' '\$1=="preferred_model"{print \$2}'  ${rec_tsv})
    BAND=\$(awk -F'\\t'  '\$1=="nonsig_band"{print \$2}'       ${rec_tsv})
    BREC=\$(awk -F'\\t' '\$1=="beta_group_at_rec"{print \$2}'  ${rec_tsv})
    PREC=\$(awk -F'\\t' '\$1=="p_group_at_rec"{print \$2}'     ${rec_tsv})
    VREC=\$(awk -F'\\t' '\$1=="variant_count_at_rec"{print \$2}' ${rec_tsv})
    AGREE=\$(awk -F'\\t' '\$1=="models_agree_at_rec"{print \$2}' ${rec_tsv})
    SENSREC=\$(awk -F'\\t' '\$1=="sensitivity_rec"{print \$2}'  ${rec_tsv})

    # machine-readable per-cohort summary (feeds the cross-cohort comparison)
    printf 'cohort\\tsamples_removed\\tvariants_excluded\\tsamples_retained\\tvariants_retained\\tknee\\tretained_pct\\trecommended_minac\\tpreferred_model\\tcalibrated\\tnonsig_band\\tbeta_group_at_rec\\tp_group_at_rec\\tvariant_count_at_rec\\tpc_sensitivity_rec\\tmodels_agree_at_rec\\n' > cohort_summary.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${cohort}" "\$SM_RM" "\$VE" "\$SR" "\$VR" "\$KNEE" "\$RETP" "\$REC" "\$PMODEL" "\$CAL" "\$BAND" "\$BREC" "\$PREC" "\$VREC" "\$SENSREC" "\$AGREE" >> cohort_summary.tsv

    # README header (Nextflow fills cohort/params/run_info; bash appends results)
    cat > README.md <<'RVEOF'
# tuning.rv — run summary (${cohort})

Rare-variant depth/platform-confounding QC + minAC sweep. This run QUANTIFIES and
MITIGATES the confounding and recommends a minAC; it does not confirm associations
(the design is fully confounded: cases 30x, controls 15x, no platform overlap).
Methodology: `../../docs/METHODS.md` · Output dictionary: `../../docs/OUTPUTS.md`.

## Resolved parameters
| Parameter | Value |
| --- | --- |
| Cohort | `${cohort}` |
| Rare call set | `${params.BedTemplate.replace('@@COHORT@@', cohort)}` |
| Info file | `${params.InfoPath}` |
| Ancestry PC (+PC model) | `${params.PcTemplate.replace('@@COHORT@@', cohort)}` (${params.PcLabel}, ${params.NPcs} PCs) |
| minAC sweep | 0-${params.MaxAC} |
| OLS SE | ${params.Robust ? 'HC3 (robust)' : 'classical'} |
| Robust-Z threshold | ${params.RobustZThreshold} |
| Kneedle (S, wx, wy) | ${params.KneeS}, ${params.KneeWeightX}, ${params.KneeWeightY} |
| minAC recommendation | V-trough; primary model = ${params.SelPrefer} (alpha ${params.Alpha}, floor ${params.MinACFloor}) |
| Conda env | ${params.conda_env_activate} |

## Run
| Field | Value |
| --- | --- |
${run_info}
| Timing | see `../_run_info/report.html` |

## Protected variants (never excluded by variant QC)
${params.ProtectVariants.collect { "- `${it}`" }.join('\n')}

## Output tree
- `figures/` - the four publication figures gathered in one place (sample QC, rate distribution, variant QC, minAC sweep)
- `00.precheck/metrics/` - base metrics (permanent storeDir cache) + `calc_metrics.log`
- `00.precheck/01_sample_qc/` - `sample_outliers.remove.txt`, `sample_qc_rate_residual.png`, `rate_distribution.png`, logs
- `00.precheck/02_variant_qc/` - recomputed `variant_metrics.txt.gz` (symlink), `variants_depthdiff.exclude.txt`, `depthdiff_variants_cdf.png`, `depthdiff_stats.tsv`, logs
- `01.sample_filter/` - sample-cleaned genotype (symlink) + `filter_samples.log`
- `02.callset_filter/` - fully-filtered genotype (symlink) + `filter_variants.log`
- `03.qc_metrics/minac*/` - per-minAC metrics (symlink) + `calc_metrics.log`
- `04.qc_summary/minac*/` - per-minAC stats TSV + log
- `05.qc_collect/` - `all_qc_summary_stats.tsv`, `minac_recommendation.tsv`, `minac_sweep.png`, `cohort_summary.tsv`
- `../_comparison/` - cross-cohort comparison (figure + table)
- `../_run_info/` - Nextflow trace / report / timeline / dag
RVEOF

    {
      echo ""
      echo "## Result (in → out, and the recommended minAC)"
      echo ""
      echo "| Metric | Value |"
      echo "| --- | --- |"
      printf '| Samples removed (rate-residual outliers) | %s |\\n' "\$SM_RM"
      printf '| Variants excluded (depth-diff) | %s |\\n'           "\$VE"
      printf '| Samples retained | %s |\\n'                          "\$SR"
      printf '| Variants retained (minAC 0) | %s |\\n'               "\$VR"
      printf '| Depth-diff knee | %s (retained %s%%) |\\n'           "\$KNEE" "\$RETP"
      printf '| Protected variants (never excluded) | %s present, %s rescued from exclusion |\\n' "\$PROT_PRES" "\$PROT_RESC"
      printf '| **Recommended minAC** (V-trough) | **%s**  (calibrated: %s) |\\n' "\$REC" "\$CAL"
      printf '| Primary model / non-sig band | %s / {%s} |\\n'       "\$PMODEL" "\$BAND"
      printf '| beta_group @ rec (%s) | %s  (p = %s) |\\n'          "\$PMODEL" "\$BREC" "\$PREC"
      printf '| Variants kept @ rec | %s |\\n'                       "\$VREC"
      printf '| +PC sensitivity | recommends minAC %s; agrees at rec: %s |\\n' "\$SENSREC" "\$AGREE"
    } >> README.md

    echo "[WRITE_RUN_MANIFEST] wrote README.md + cohort_summary.tsv for ${cohort} (rec minAC=\$REC)"
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// COMPARE_COHORTS — overlay the cohorts' sweeps + a comparison table.
// ═══════════════════════════════════════════════════════════════════════════
process COMPARE_COHORTS {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir "${params.ResultsDir}/_comparison", mode: 'copy', pattern: "{cohort_comparison.*,decision_axis.*,all_cohorts_*.tsv,compare_cohorts.log}"
    publishDir "${params.ResultsDir}", mode: 'copy', pattern: "results_README.md", saveAs: { 'README.md' }

    input:
    path stats_combined
    path summ_combined

    output:
    path "cohort_comparison.png", emit: pdf
    path "cohort_comparison.tsv", emit: table
    path "decision_axis.png", emit: axis_fig
    path "decision_axis.tsv", emit: axis_table
    path "all_cohorts_qc_stats.tsv", emit: combined_stats
    path "all_cohorts_summary.tsv", emit: combined_summary
    path "compare_cohorts.log", emit: log
    path "results_README.md", emit: index

    script:
    """
    ${params.PRELUDE}
    source ${params.ScriptDir}/rv_log.sh
    rv_log_init "COMPARE_COHORTS" "compare_cohorts.log"

    cp ${stats_combined} all_cohorts_qc_stats.tsv
    cp ${summ_combined}  all_cohorts_summary.tsv

    # Top-level results index (published to results/README.md).
    cat > results_README.md <<'IDXEOF'
# tuning.rv — results

Rare-variant depth/platform-confounding QC + a null-calibrated minAC recommendation.
This does **not** confirm associations (the design is fully confounded: cases 30x,
controls 15x, no platform overlap) — it quantifies and mitigates the confounding.

## Where to start
1. **Per cohort** → `<cohort>/README.md` (params, in→out counts, recommended minAC,
   β_g at the recommendation) and `<cohort>/figures/` (the four figures).
2. **Across cohorts** → `_comparison/cohort_comparison.png` + `.tsv`.
3. **Methods & data dictionary** → `../docs/METHODS.md`, `../docs/OUTPUTS.md`.

## Cohorts
`narrow_mainland` ⊂ `intermediate_mainland` ⊂ `full_mainland` — nested sample-selection variants.

## Per-cohort layout
- `figures/` — the four figures gathered in one place (sample QC, rate distribution,
  variant QC, minAC sweep). Each is also in its stage dir.
- `00.precheck/metrics/` — base metrics (permanent cache; the large `variant_metrics.txt.gz`).
- `00.precheck/01_sample_qc/`, `02_variant_qc/` — detection outputs.
- `01.sample_filter/`, `02.callset_filter/` — filtered genotypes (symlinks into work/).
- `03.qc_metrics/minac*/`, `04.qc_summary/minac*/` — per-minAC provenance.
- `05.qc_collect/` — merged stats, `minac_recommendation.tsv`, `minac_sweep.png`.

## Shared
- `_comparison/` — cross-cohort figure + table. · `_run_info/` — Nextflow trace / report / timeline / dag.

## Reading the sweep
**β_g** (Case−Control) is the DECISION axis; **β_d** (depth) is a DIAGNOSTIC only —
do NOT pick minAC to minimise it. Choose the minAC at the β_g non-significant trough;
confirm downstream with the genome-wide inflation factor λ (PC-adjusted rv test).
IDXEOF

    python3 ${params.ScriptDir}/compare_cohorts.py \\
        --stats-combined all_cohorts_qc_stats.tsv \\
        --summary-combined all_cohorts_summary.tsv \\
        --alpha ${params.Alpha} \\
        --out-png cohort_comparison.png \\
        --out-tsv cohort_comparison.tsv 2>&1 | tee -a compare_cohorts.log

    # Why the minAC decision is read from beta_group in the no-PC model. Run-level, so
    # it lives here rather than in a process of its own; needs no input the workflow
    # does not already resolve.
    python3 ${params.ScriptDir}/plot_decision_axis.py \\
        --stats-combined all_cohorts_qc_stats.tsv \\
        --platform-file "${params.PlatformFile}" \\
        --platform-id-col "${params.PlatformIdCol}" \\
        --platform-col "${params.PlatformCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --pc-template "${params.PcTemplate}" \\
        --n-pcs ${params.NPcs} \\
        --alpha ${params.Alpha} \\
        --reference-cohort "${params.DecisionAxisCohort}" \\
        --out-png decision_axis.png \\
        --out-tsv decision_axis.tsv 2>&1 | tee -a compare_cohorts.log
    rv_done
    """
}

workflow {
    // Normalize --Cohorts (list default, or CLI comma-string) + validate.
    def cohortList = (params.Cohorts instanceof List) ? params.Cohorts
                     : params.Cohorts.toString().split(',').collect { it.trim() }.findAll { it }
    cohortList.each { c ->
        if (!params.ValidCohorts.contains(c)) {
            error "Invalid cohort '${c}'. Must be one of: ${params.ValidCohorts.join(', ')}"
        }
    }
    log.info "[tuning.rv] cohorts: ${cohortList.join(', ')}  | minAC 0..${params.MaxAC}  | PrecheckOnly=${params.PrecheckOnly}"

    // Per-cohort raw rare call set (checked to exist up-front).
    cohorts_ch = channel.fromList(cohortList).map { c ->
        def bp = params.BedTemplate.replace('@@COHORT@@', c)
        tuple(c,
              file("${bp}.bed", checkIfExists: true),
              file("${bp}.bim", checkIfExists: true),
              file("${bp}.fam", checkIfExists: true))
    }

    // Per-cohort ancestry PC file (bbj_mainland) for the +PC sweep model.
    pc_ch = channel.fromList(cohortList).map { c ->
        tuple(c, file(params.PcTemplate.replace('@@COHORT@@', c), checkIfExists: true))
    }

    // ── Pre-check: base metrics -> SAMPLE QC (first) ─────────────────────────────
    base = CALC_METRICS_BASE(cohorts_ch)
    outl = DETECT_OUTLIER_SAMPLES(base.sample_metrics)
    PLOT_RATE_DISTRIBUTION(base.sample_metrics.join(outl.remove))

    // ── Sample filter -> recompute variant metrics -> VARIANT QC ─────────────────
    samp_filt = FILTER_SAMPLES(outl.remove)
    vmet      = CALC_VARIANT_METRICS(samp_filt.cleaned)
    depth     = DETECT_DEPTHDIFF_VARIANTS(vmet.variant_metrics)

    if (params.PrecheckOnly) {
        log.info "[tuning.rv] PrecheckOnly=true -> stopped after detection (sample + variant QC). See <cohort>/00.precheck/."
    } else {

    // ── Variant filter -> fully-filtered set ─────────────────────────────────────
    var_in = samp_filt.cleaned.join(depth.exclude)
    filt   = FILTER_VARIANTS(var_in)

    // ── minAC sweep on the filtered set (cohort × minAC) ─────────────────────────
    sweep_in = filt.filtered.combine(channel.from(0..params.MaxAC))
    calc     = CALC_METRICS(sweep_in)

    // Attach each cohort's PC file to every per-minAC metrics tuple, then summarise.
    qc_in = calc.metrics.combine(pc_ch, by: 0)
    qc    = QC_SUMMARY(qc_in)

    // Group per-minAC stats back per cohort for the merge.
    per_cohort_tsvs = qc.summary_stats.map { c, _minac, tsv -> tuple(c, tsv) }.groupTuple()
    merged = MERGE_QC_SUMMARIES(per_cohort_tsvs)

    // Recommend a minAC + draw the headline sweep (marks the recommendation).
    rec   = SELECT_MINAC(merged.merged_stats)
    sweep = PLOT_SWEEP(merged.merged_stats.join(rec.rec))

    // Run-level provenance (built once, broadcast to every cohort's manifest).
    run_info = "| Run name | ${workflow.runName} |\n| Started | ${workflow.start} |\n| Command | `${workflow.commandLine}` |\n| Revision | ${workflow.revision ?: 'n/a'} |\n| Nextflow | ${workflow.nextflow.version} |"

    // Per-cohort manifest: join the terminal artifacts by cohort key.
    man_in   = samp_filt.log.join(filt.log).join(depth.stats).join(rec.rec).join(sweep.pdf)
    manifest = WRITE_RUN_MANIFEST(man_in, run_info)

    // ── Cross-cohort comparison ─────────────────────────────────────────────────
    // Distinct staged names (≠ the copied outputs) so COMPARE_COHORTS' `cp`
    // doesn't hit a same-file collision.
    combined_stats = merged.merged_stats.map { _c, f -> f }
        .collectFile(name: 'combined_qc_stats.tsv', keepHeader: true, skip: 1)
    combined_summ  = manifest.summary.map { _c, f -> f }
        .collectFile(name: 'combined_summary.tsv', keepHeader: true, skip: 1)
    COMPARE_COHORTS(combined_stats, combined_summ)
    }
}
