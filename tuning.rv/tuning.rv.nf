nextflow.enable.dsl = 2

params.GTPath      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results'
params.InfoPath    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx'
params.OutDir      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/results'
params.ScriptDir   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/scripts'
params.Plink2      = '/home/b/b37974/plink2_alpha6/plink2'
params.Tabix       = '/home/b/b37974/htslib-1.9/tabix'
params.IDCol       = 'ID_JHRPv6'
params.GroupCol    = 'Outcome'
params.CaseValue   = 'PH'
params.TdpCol      = 'Target_Depth'
params.MdpCol      = 'Observed_Depth'
params.SampleRm    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/pre_check/results/01_outliers_z3/outliers_sminac_robustz.remove.txt'
params.ThreadsCalc = 4

// ─────────────────────────────────────────────────────────────────────────
// Remove outlier samples from the genotype file BEFORE metric calculation.
// Only executed when params.SampleRm resolves to an existing file.
// ─────────────────────────────────────────────────────────────────────────
process REMOVE_SAMPLES {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutDir}/00.sample_rm", mode: 'symlink'

    output:
    tuple path("filtered.bed"), path("filtered.bim"), path("filtered.fam"), emit: filtered
    path "remove_samples.log", emit: log

    script:
    def orig_prefix = "${params.GTPath}/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold"
    """
    set -euo pipefail

    echo "[REMOVE_SAMPLES] Removing samples listed in: ${params.SampleRm}" | tee remove_samples.log
    echo "[REMOVE_SAMPLES] Input prefix: ${orig_prefix}"               >> remove_samples.log

    n_rm=\$(wc -l < "${params.SampleRm}")
    echo "[REMOVE_SAMPLES] Samples to remove: \${n_rm}"                >> remove_samples.log

    ${params.Plink2} \\
        --bfile "${orig_prefix}" \\
        --remove "${params.SampleRm}" \\
        --make-bed \\
        --out filtered \\
        --threads ${params.ThreadsCalc} \\
        >> remove_samples.log 2>&1

    n_after=\$(wc -l < filtered.fam)
    echo "[REMOVE_SAMPLES] Samples retained: \${n_after}"             >> remove_samples.log
    """
}

process CALC_METRICS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "minac:${min_ac}"
    publishDir "${params.OutDir}/00.qc_metrics/minac${min_ac}", mode: 'symlink'

    input:
    val min_ac
    tuple path(bed), path(bim), path(fam)   // staged genotype files (original or sample-removed)

    output:
    tuple val(min_ac), path("sample_metrics.txt.gz"), path("variant_metrics.txt.gz"), path("variant_metrics.txt.gz.tbi"), emit: metrics
    path "calc_metrics.log", emit: log

    script:
    """
    set -euo pipefail

    bash ${params.ScriptDir}/run_calc_metrics.sh \\
        --bed-prefix "${bed.baseName}" \\
        --info-file "${params.InfoPath}" \\
        --out-sample "sample_metrics.txt" \\
        --out-variant "variant_metrics.txt" \\
        --script-dir "${params.ScriptDir}" \\
        --plink2 "${params.Plink2}" \\
        --tabix "${params.Tabix}" \\
        --id-col "${params.IDCol}" \\
        --group-col "${params.GroupCol}" \\
        --case-value "${params.CaseValue}" \\
        --tdp-col "${params.TdpCol}" \\
        --mdp-col "${params.MdpCol}" \\
        --min-ac ${min_ac} \\
        --threads ${params.ThreadsCalc}
    """
}

process QC_SUMMARY {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "minac:${min_ac}"
    publishDir "${params.OutDir}/01.qc_summary/minac${min_ac}", mode: 'symlink'

    input:
    tuple val(min_ac), path(sample_metrics), path(variant_metrics), path(variant_metrics_tbi)

    output:
    tuple val(min_ac), path("qc_summary_stats.minac${min_ac}.tsv"), emit: summary_stats
    path "qc_summary_plots.pdf", emit: summary_plot
    path "*.md", emit: summary_docs

    script:
    """
    set -euo pipefail

    python3 ${params.ScriptDir}/qc_summary.py \
        --sample-metrics ${sample_metrics} \
        --variant-metrics ${variant_metrics} \
        --min-ac ${min_ac} \
        --out-tsv qc_summary_stats.minac${min_ac}.tsv \
        --out-md qc_summary_methods.md

    python3 ${params.ScriptDir}/plot_qc.py \
        --sample-metrics ${sample_metrics} \
        --variant-metrics ${variant_metrics} \
        --qc-stats qc_summary_stats.minac${min_ac}.tsv \
        --min-ac ${min_ac} \
        --out-pdf qc_summary_plots.pdf \
        --hist-stat density
    """
}

process MERGE_ALL_SUMMARIES {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir "${params.OutDir}/02.qc_collect", mode: 'symlink'

    input:
    path tsvs

    output:
    path "all_qc_summary_stats.tsv", emit: merged_stats

    script:
    """
    set -euo pipefail

    awk 'NR==1 || FNR > 1' ${tsvs} > temp_merged.tsv
    python3 - <<'PY'
import pandas as pd

df = pd.read_csv('temp_merged.tsv', sep='\t')
df = df.sort_values('MinAC_Threshold')
df.to_csv('all_qc_summary_stats.tsv', sep='\t', index=False)
PY
    """
}

process PLOT_SUMMARY_TREND {
    executor 'slurm'
    queue 'gr10478b'
    time '10m'
    publishDir "${params.OutDir}/02.qc_collect", mode: 'symlink'

    input:
    path qc_stats

    output:
    path "qc_trend_plots.pdf", emit: trend_plot

    script:
    """
    set -euo pipefail

    N_SAMPLES=\$(wc -l < ${params.GTPath}/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold.fam)
    python3 ${params.ScriptDir}/plot_summary_trend.py \
        --qc-stats ${qc_stats} \
        --out-pdf qc_trend_plots.pdf \
        --sample-n \${N_SAMPLES}
    """
}

workflow {
    def orig_prefix = "${params.GTPath}/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold"
    def rm_file = file(params.SampleRm, checkIfExists: false)

    minac_ch = channel.from(0..20)

    // ── Conditional sample removal ────────────────────────────────────────────
    if (rm_file.exists()) {
        log.info "[workflow] SampleRm found: ${params.SampleRm} — running REMOVE_SAMPLES"
        bed_files_ch = REMOVE_SAMPLES().filtered
    } else {
        log.info "[workflow] SampleRm not found — using original genotype files"
        bed_files_ch = channel.value([
            file("${orig_prefix}.bed"),
            file("${orig_prefix}.bim"),
            file("${orig_prefix}.fam")
        ])
    }
    // ─────────────────────────────────────────────────────────────────────────

    calc_out = CALC_METRICS(minac_ch, bed_files_ch)
    qc_out = QC_SUMMARY(calc_out.metrics)

    summary_tsvs_ch = qc_out.summary_stats
        .map { _min_ac, tsv -> tsv }
        .collect()

    merged_out = MERGE_ALL_SUMMARIES(summary_tsvs_ch)
    PLOT_SUMMARY_TREND(merged_out.merged_stats)
}


