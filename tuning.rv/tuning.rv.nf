nextflow.enable.dsl = 2

// ── Inputs / paths ─────────────────────────────────────────────────────────
params.CohortName  = 'refined_core' // refined_core, expanded_core, or full_mainland
params.GTPath      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results'
// Single source of truth for the genotype prefix (reused by every process).
params.BedPrefix   = "${params.GTPath}/14_fixed_model_prep/${params.CohortName}/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold"
params.InfoPath    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/info/cteph_agp3k.v6.20260507.xlsx'
params.OutDir      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/results'
params.ScriptDir   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/scripts'

// Conda env activated in every process; plink2 / tabix / bgzip / python come
// from PATH + this env (see the prelude at the top of each process script).
params.conda_env_activate = 'cteph_geno_pro'

// ── Info-file column mapping ────────────────────────────────────────────────
params.IDCol       = 'ID_JHRPv6'
params.GroupCol    = 'Outcome'
params.CaseValue   = 'PH'
params.TdpCol      = 'Target_Depth'
params.MdpCol      = 'Observed_Depth'

// ── QC removal lists (both OPTIONAL; each applied only if the file exists) ───
//   SampleRm  : FID<TAB>IID per line   -> plink2 --remove   (pre_check step 1)
//   VariantRm : VariantID  per line    -> plink2 --exclude  (pre_check step 3)
params.SampleRm    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/pre_check/results/01_resid_outliers_z5/outliers_sminac_resid_robustz.remove.txt'
params.VariantRm   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/pre_check/results/02_vmiss_depthdiff_knee/variants_vmiss_depthdiff_gt0.0330.exclude.txt'

params.ThreadsCalc = 4

// ─────────────────────────────────────────────────────────────────────────
// FILTER_CALLSET — filter the genotype call set BEFORE metric calculation.
//   - removes outlier samples   (--remove,  if params.SampleRm exists)
//   - excludes flagged variants (--exclude, if params.VariantRm exists)
// At least one list must exist; the workflow only invokes this when so.
// ─────────────────────────────────────────────────────────────────────────
process FILTER_CALLSET {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutDir}/00.callset_filter", mode: 'symlink'

    output:
    tuple path("filtered.bed"), path("filtered.bim"), path("filtered.fam"), emit: filtered
    path "filter_callset.log", emit: log

    script:
    def has_sample   = params.SampleRm  && file(params.SampleRm).exists()
    def has_variant  = params.VariantRm && file(params.VariantRm).exists()
    def remove_flag  = has_sample  ? "--remove \"${params.SampleRm}\""   : ""
    def exclude_flag = has_variant ? "--exclude \"${params.VariantRm}\"" : ""
    """
    set -euo pipefail
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}

    {
      echo "[FILTER_CALLSET] Input prefix : ${params.BedPrefix}"
      echo "[FILTER_CALLSET] Sample list  : ${has_sample  ? params.SampleRm  : '(none)'}"
      echo "[FILTER_CALLSET] Variant list : ${has_variant ? params.VariantRm : '(none)'}"
    } | tee filter_callset.log

    if [[ -n "${remove_flag}" ]]; then
        echo "[FILTER_CALLSET] Samples to remove  : \$(wc -l < "${params.SampleRm}")"  >> filter_callset.log
    fi
    if [[ -n "${exclude_flag}" ]]; then
        echo "[FILTER_CALLSET] Variants to exclude: \$(wc -l < "${params.VariantRm}")" >> filter_callset.log
    fi

    plink2 --bfile "${params.BedPrefix}" ${remove_flag} ${exclude_flag} \\
        --make-bed --out filtered --threads ${params.ThreadsCalc} \\
        >> filter_callset.log 2>&1

    echo "[FILTER_CALLSET] Samples retained  : \$(wc -l < filtered.fam)" >> filter_callset.log
    echo "[FILTER_CALLSET] Variants retained : \$(wc -l < filtered.bim)" >> filter_callset.log
    """
}

// ─────────────────────────────────────────────────────────────────────────
// CALC_METRICS — sample/variant metrics via run_calc_metrics.sh.
// The script also emits depth-stratified missingness columns
// (VMISS, VMISS_15X, VMISS_30X) in variant_metrics; the CLI is unchanged.
// ─────────────────────────────────────────────────────────────────────────
process CALC_METRICS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "minac:${min_ac}"
    publishDir "${params.OutDir}/00.qc_metrics/minac${min_ac}", mode: 'symlink'

    input:
    val min_ac
    tuple path(bed), path(bim), path(fam)   // staged genotype files (original or filtered)

    output:
    tuple val(min_ac), path("sample_metrics.txt.gz"), path("variant_metrics.txt.gz"), path("variant_metrics.txt.gz.tbi"), emit: metrics
    path "calc_metrics.log", emit: log

    script:
    """
    set -euo pipefail
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}

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
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}

    python3 ${params.ScriptDir}/qc_summary.py \\
        --sample-metrics ${sample_metrics} \\
        --variant-metrics ${variant_metrics} \\
        --min-ac ${min_ac} \\
        --out-tsv qc_summary_stats.minac${min_ac}.tsv \\
        --out-md qc_summary_methods.md

    python3 ${params.ScriptDir}/plot_qc.py \\
        --sample-metrics ${sample_metrics} \\
        --variant-metrics ${variant_metrics} \\
        --qc-stats qc_summary_stats.minac${min_ac}.tsv \\
        --min-ac ${min_ac} \\
        --out-pdf qc_summary_plots.pdf \\
        --hist-stat density
    """
}

process MERGE_QC_SUMMARIES {
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
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}

    awk 'NR==1 || FNR > 1' ${tsvs} > temp_merged.tsv
    python3 - <<'PY'
import pandas as pd

df = pd.read_csv('temp_merged.tsv', sep='\\t')
df = df.sort_values('MinAC_Threshold')
df.to_csv('all_qc_summary_stats.tsv', sep='\\t', index=False)
PY
    """
}

process PLOT_QC_TREND {
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
    export PATH=/home/b/b37974/:/home/b/b37974/htslib-1.9/:\$PATH
    source activate ${params.conda_env_activate}

    N_SAMPLES=\$(wc -l < "${params.BedPrefix}.fam")
    python3 ${params.ScriptDir}/plot_summary_trend.py \\
        --qc-stats ${qc_stats} \\
        --out-pdf qc_trend_plots.pdf \\
        --sample-n \${N_SAMPLES}
    """
}

workflow {
    minac_ch = channel.from(0..20)

    // Resolve removal lists (both optional; applied only if the file exists).
    def sampleRmFile  = params.SampleRm  ? file(params.SampleRm,  checkIfExists: false) : null
    def variantRmFile = params.VariantRm ? file(params.VariantRm, checkIfExists: false) : null
    def hasSampleRm   = sampleRmFile  != null && sampleRmFile.exists()
    def hasVariantRm  = variantRmFile != null && variantRmFile.exists()

    // ── Conditional genotype filtering (samples and/or variants) ─────────────
    if (hasSampleRm || hasVariantRm) {
        log.info "[workflow] genotype filtering ON  (samples=${hasSampleRm}, variants=${hasVariantRm}) — running FILTER_CALLSET"
        genotypes_ch = FILTER_CALLSET().filtered
    } else {
        log.info "[workflow] no removal lists found — using original genotype files"
        genotypes_ch = channel.value([
            file("${params.BedPrefix}.bed"),
            file("${params.BedPrefix}.bim"),
            file("${params.BedPrefix}.fam")
        ])
    }
    // ─────────────────────────────────────────────────────────────────────────

    calc_out = CALC_METRICS(minac_ch, genotypes_ch)
    qc_out = QC_SUMMARY(calc_out.metrics)

    summary_tsvs_ch = qc_out.summary_stats
        .map { _min_ac, tsv -> tsv }
        .collect()

    merged_out = MERGE_QC_SUMMARIES(summary_tsvs_ch)
    PLOT_QC_TREND(merged_out.merged_stats)
}
