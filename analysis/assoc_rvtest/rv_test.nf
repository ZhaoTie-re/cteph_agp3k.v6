#!/usr/bin/env nextflow
/*
=========================================================================================
    rvtests-rvat  ·  Gene-based rare-variant association pipeline
=========================================================================================
    A reusable Nextflow (DSL2) workflow that runs gene-based rare-variant association
    testing with RVTESTS (SKAT-O and CMC burden) on a PLINK genotype set, sweeping over
    variant-impact strata × test methods. The workflow logic is project-agnostic; only
    the inputs in the PARAMETERS block are project-specific.

    Default configuration targets : CTEPH AGP3K v6 (WGS, hg38).

    Pipeline DAG
    ------------
        PLINK{bed,bim,fam}
              │
        [0] FILTER_GENOTYPE (optional)  ──(QC'd PLINK)──┐
              │ (skipped if --filterGenotype false)     │
              ▼                                          ▼
        PLINK + pheno + covar + refFlat ──────▶ [1] RVTEST_PREPARE
              │
        [1] ──(vcf)──▶ [2] SNPEFF_ANNOTATE ──(stats)──▶ [3] PLOT_SNPEFF_STATS
              │                │
              │           (annotated vcf)
              │                ▼
              │          [4] INFO_FILTER ──(impact-filtered vcf)──┐
              │                                                   ▼
              └──(pheno/covar + refFlat)──────────────▶ [5] RVTEST_RUN
                                                                  │
                                                               (assoc)
                                                                  ▼
                                                      [6] RVTEST_POST_PROCESS
                                                                  │
                                                             (fdr assoc)
                                                                  ▼
                                                      [7] RVTEST_VISUALIZATION

    Stage summary
    -------------
      0 FILTER_GENOTYPE       (optional) plink2 QC of the genotype matrix: remove outlier
                              samples/variants (Logic 1) + MAC filter (Logic 2)
      1 RVTEST_PREPARE        PLINK (bed/bim/fam) → bgzipped VCF + pheno/covar + refFlat
      2 SNPEFF_ANNOTATE       Functional annotation of the VCF with snpEff
      3 PLOT_SNPEFF_STATS     QC plot of the impact / effect distribution
      4 INFO_FILTER           Subset variants by predicted IMPACT stratum
      5 RVTEST_RUN            Gene-based association tests (impact strata × methods)
      6 RVTEST_POST_PROCESS   Min-NumVar filtering + Benjamini-Hochberg FDR
      7 RVTEST_VISUALIZATION  Manhattan / QQ plots per result

    Inputs   (PARAMETERS block — override on the CLI or with -params-file)
      inputPlinkPrefix  PLINK prefix, expects <prefix>.{bed,bim,fam}
      phenoFile         Phenotype table  (rvtest TSV)
      covarFile         Covariate table  (rvtest TSV)
      refFlatPath       refFlat gene model (.txt.gz)
      snpeffDir         snpEff index directory
      filterGenotype    Toggle the STEP 0 genotype QC filter (default true)

    Outputs  (published under params.outDir)
      00.filter_genotype · 01.rvtest_prepare · 02.snpeff_annotate · 03.info_filter
      04.rvtest_run · 05.post_process · 06.visualization

    Reuse on a new project
    ----------------------
      Change ONLY the "Required inputs" parameters (edit below, or pass on the CLI):
        nextflow run rv_test.nf \
            --inputPlinkPrefix <prefix> --phenoFile <tsv> --covarFile <tsv> \
            --refFlatPath <refFlat> --snpeffDir <index> --outDir <dir>
      Nothing in the process / workflow logic is hard-coded to a particular dataset.

    Requirements
    ------------
      Nextflow >=22.10 · SLURM scheduler · conda env (params.condaEnv) ·
      rvtest binary (on PATH) · helper scripts (params.scriptDir)

    Configuration split
    -------------------
      rv_test.nf        pipeline parameters + process/workflow logic (this file)
      nextflow.config   SLURM `--rsc` resource requests only

    Usage
    -----
      nextflow run rv_test.nf -resume
      nextflow run rv_test.nf --help
=========================================================================================
*/

nextflow.enable.dsl = 2

// =========================================================================================
//  PARAMETERS
// -----------------------------------------------------------------------------------------
//  Edit the values here, or override any of them on the command line, e.g.
//      nextflow run rv_test.nf --performNorm true --rvtestThreads 8
//  To port the pipeline to another project, only the "Required inputs" group needs editing.
// =========================================================================================

// ---- Project identity (label only; used for tags and the run banner) --------------------
params.project           = 'cteph_agp3k_v6'

// ----------------------------------------------------------------------------------------
//  Required inputs — the only project-specific values. Point these at your dataset.
// ----------------------------------------------------------------------------------------
// PLINK genotype set — the pipeline expects <inputPlinkPrefix>.{bed,bim,fam}
params.inputPlinkPrefix  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/14_fixed_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.maf_lt_threshold'

// Phenotype / covariate tables (rvtest format, tab-separated)
params.phenoFile         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/16_cov_pheno_prep/refined_core/popgmm_subset_on_bbj_pcs.pheno.tsv'
params.covarFile         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/16_cov_pheno_prep/refined_core/popgmm_subset_on_bbj_pcs.cov.sex.tsv'

// Reference resources
params.refFlatPath       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz'
params.snpeffDir         = '/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/output/snpEff.v6.index'

// ---- Output ----------------------------------------------------------------------------
params.outDir            = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_rvtest/results'
// Naming stem for all result files; defaults to the PLINK file basename.
params.outPrefix         = params.inputPlinkPrefix.tokenize('/').last()

// ---- Tools / environment ---------------------------------------------------------------
params.scriptDir         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_rvtest/scripts'
params.condaEnv          = 'cteph_geno_pro'

// ---- Analysis options ------------------------------------------------------------------
params.performNorm       = false        // bcftools normalisation during PREPARE
params.allowRefAltSwap   = false        // accept REF/ALT-swapped sites during annotation
params.num_var_threshold = 3            // min variants per gene to retain post-test
params.phenoName         = 'pheno1'                                       // pheno column
params.covarName         = 'sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg' // covar columns — MUST match the covarFile header (v6 PCs are named pc<N>_avg)

// ---- Genotype-matrix QC filter (STEP 0; disable with --filterGenotype false) -----------
//   Logic 1: drop outlier samples (--remove) and flagged variants (--exclude); each is
//            applied only if its file exists. Logic 2: drop variants with MAC < threshold.
//   Order is Logic 1 then Logic 2 (plink2 computes MAC after the sample/variant removal).
//   If neither Logic-1 file exists, only Logic 2 runs; if filterGenotype=false, both skip.
params.filterGenotype     = true
params.sampleRemoveFile   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/pre_check/results/01_resid_outliers_z5/outliers_sminac_resid_robustz.remove.txt'
params.variantExcludeFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tuning.rv/pre_check/results/02_vmiss_depthdiff_knee/variants_vmiss_depthdiff_gt0.0330.exclude.txt'
params.macThreshold       = 2

// ---- Per-process parallelism knobs -----------------------------------------------------
params.filterThreads     = 8
params.prepareThreads    = 32
params.snpeffWorkers     = 16
params.snpeffThreads     = 16
params.infoFilterThreads = 16
params.rvtestThreads     = 4

// ---- Test matrix: variant-impact strata × gene-based tests -----------------------------
//   Add/remove rows to extend the sweep; every (impact × method) pair is run independently.
params.impactFilters = [
    [ tag: 'low_moderate_high', values: 'LOW MODERATE HIGH' ],
    [ tag: 'moderate_high',     values: 'MODERATE HIGH'     ],
    [ tag: 'high',              values: 'HIGH'              ],
]
params.testMethods = [
    [ tag: 'skato',   opt: '--kernel skato'   ],   // SKAT-O (kernel)
    [ tag: 'cmc',     opt: '--burden cmc'     ],   // CMC burden
    [ tag: 'zeggini', opt: '--burden zeggini' ],   // Zeggini burden
]

// ---- Misc ------------------------------------------------------------------------------
params.help     = false

// =========================================================================================
//  HELPER FUNCTIONS
// =========================================================================================

/*
 * Print CLI help and the resolved default of every user-facing parameter.
 */
def helpMessage() {
    log.info """
    =====================================================================================
     rvtests-rvat  ·  gene-based rare-variant association pipeline
    =====================================================================================

    Usage:
      nextflow run rv_test.nf [-resume] [options]

    Key options:
      --outDir            Results directory                     [${params.outDir}]
      --inputPlinkPrefix  PLINK .bed/.bim/.fam prefix           [${params.inputPlinkPrefix}]
      --phenoFile         Phenotype TSV                         [${params.phenoFile}]
      --covarFile         Covariate TSV                         [${params.covarFile}]
      --refFlatPath       refFlat gene model (.txt.gz)          [${params.refFlatPath}]
      --snpeffDir         snpEff index directory                [${params.snpeffDir}]
      --performNorm       Run bcftools normalisation            [${params.performNorm}]
      --allowRefAltSwap   Allow REF/ALT swap during annotation  [${params.allowRefAltSwap}]
      --num_var_threshold Minimum NumVar per gene to keep       [${params.num_var_threshold}]
      --phenoName         Phenotype column name                 [${params.phenoName}]
      --covarName         Covariate column list                 [${params.covarName}]
      --help              Print this help and exit
    =====================================================================================
    """.stripIndent()
}

/*
 * Echo the resolved run configuration to the log at start-up.
 */
def startupSummary() {
    log.info """
=========================================================================================
 rvtests-rvat  |  ${params.project}
=========================================================================================
 PLINK prefix    : ${params.inputPlinkPrefix}
 Phenotype       : ${params.phenoFile}
 Covariates      : ${params.covarFile}
 refFlat         : ${params.refFlatPath}
 snpEff index    : ${params.snpeffDir}
 Output dir      : ${params.outDir}
 Filter genotype : ${params.filterGenotype}${params.filterGenotype ? " (remove samples/variants if present, then MAC >= ${params.macThreshold})" : ""}
 Normalise VCF   : ${params.performNorm}
 Impact strata   : ${params.impactFilters*.tag.join(', ')}
 Test methods    : ${params.testMethods*.tag.join(', ')}
=========================================================================================
""".stripIndent()
}

// =========================================================================================
//  PROCESSES
// -----------------------------------------------------------------------------------------
//  Each process declares its own executor/queue/wall-time; only the SLURM `--rsc` memory
//  and CPU requests are externalised to nextflow.config (matched by process name).
// =========================================================================================

/* ----------------------------------------------------------------------------------------
 * STEP 0 · FILTER_GENOTYPE   (optional; runs only when params.filterGenotype is true)
 *   QC-filter the PLINK genotype matrix with plink2:
 *     Logic 1 - remove outlier samples (--remove) and flagged variants (--exclude),
 *               each applied only if the corresponding file exists.
 *     Logic 2 - drop variants with MAC < params.macThreshold (--mac), computed on the
 *               sample set that remains after Logic 1 (so the order is Logic 1 -> Logic 2).
 *     in  : PLINK {bed,bim,fam} · sample-remove list · variant-exclude list
 *     out : filtered PLINK {bed,bim,fam}.gtqc
 *     via : plink2
 * ------------------------------------------------------------------------------------- */
process FILTER_GENOTYPE {
    executor 'slurm'
    queue    'gr10478b'
    time     '12h'
    tag      "${params.project}"

    publishDir "${params.outDir}/00.filter_genotype", mode: 'symlink'

    input:
    tuple path(bed), path(bim), path(fam)
    val sample_remove_file
    val variant_exclude_file

    output:
    tuple path("*.gtqc.bed"), path("*.gtqc.bim"), path("*.gtqc.fam"), emit: plink
    path("*.gtqc.log")

    script:
    def in_prefix  = bed.baseName
    def out_prefix = "${in_prefix}.gtqc"
    """
    export PATH=/home/b/b37974/plink2_alpha6:\$PATH

    # ---- Logic 1: outlier samples + flagged variants (each only if its file exists) ----
    REMOVE_ARG=""
    EXCLUDE_ARG=""
    if [ -s "${sample_remove_file}" ]; then
        REMOVE_ARG="--remove ${sample_remove_file}"
        echo "[FILTER_GENOTYPE] Logic1: removing samples from ${sample_remove_file}"
    else
        echo "[FILTER_GENOTYPE] Logic1: sample-remove file absent -> skipping --remove"
    fi
    if [ -s "${variant_exclude_file}" ]; then
        EXCLUDE_ARG="--exclude ${variant_exclude_file}"
        echo "[FILTER_GENOTYPE] Logic1: excluding variants from ${variant_exclude_file}"
    else
        echo "[FILTER_GENOTYPE] Logic1: variant-exclude file absent -> skipping --exclude"
    fi

    # ---- Logic 2: MAC filter (plink2 computes MAC after the Logic-1 sample removal) ----
    echo "[FILTER_GENOTYPE] Logic2: --mac ${params.macThreshold}"
    plink2 \\
        --bfile ${in_prefix} \\
        \$REMOVE_ARG \$EXCLUDE_ARG \\
        --mac ${params.macThreshold} \\
        --make-bed \\
        --threads ${params.filterThreads} \\
        --out ${out_prefix}
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 1 · RVTEST_PREPARE
 *   Convert the PLINK genotype set to a bgzip-compressed VCF, assemble the rvtest-format
 *   phenotype / covariate tables, and emit a chr-prefix-normalised refFlat gene model.
 *     in  : PLINK {bed,bim,fam} · pheno TSV · covar TSV · refFlat.txt.gz
 *     out : vcf (+tbi) · pheno/covar CSVs · refFlat (no-chr)
 *     via : scripts/rvtest_prepare_main.py
 * ------------------------------------------------------------------------------------- */
process RVTEST_PREPARE {
    executor 'slurm'
    queue    'gr10478b'
    time     '36h'
    tag      "${params.project}"

    publishDir "${params.outDir}/01.rvtest_prepare", mode: 'symlink'

    input:
    tuple path(bed), path(bim), path(fam)
    path pheno
    path covar
    path ref_flat

    output:
    path('*.log')
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"),                emit: vcf
    tuple path("*.pheno_rvt.tsv"), path("*.covar_rvt.tsv"),      emit: pheno_cov
    path("refFlat.hg38.nochr.txt.gz"),                           emit: refflat

    script:
    def bed_prefix = bed.baseName
    def norm_flag  = params.performNorm ? '--norm' : ''
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/rvtest_prepare_main.py \\
        --bed-prefix ${bed_prefix} \\
        --pheno-path ${pheno} \\
        --covar-path ${covar} \\
        --refflat-path ${ref_flat} \\
        --threads ${params.prepareThreads} \\
        --verbose \\
        --log-file rvtest_prepare.log \\
        ${norm_flag}
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 2 · SNPEFF_ANNOTATE
 *   Functionally annotate the prepared VCF with snpEff (adds the IMPACT INFO field used
 *   downstream for stratification).
 *     in  : vcf (+tbi) · snpEff index dir
 *     out : annotated vcf (+tbi) · per-impact/effect stats TSV
 *     via : scripts/snpeff_anno_main.py
 * ------------------------------------------------------------------------------------- */
process SNPEFF_ANNOTATE {
    executor 'slurm'
    queue    'gr10478b'
    time     '48h'
    tag      "${params.project}"

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    val snpeff_dir
    tuple path(vcf), path(vcf_tbi)

    output:
    path('*.log')
    path('*.tsv'),                                                emit: stats
    tuple path("*.snpeff.vcf.gz"), path("*.snpeff.vcf.gz.tbi"),   emit: annotated_vcf

    script:
    def swap_flag = params.allowRefAltSwap ? '--allow-ref-alt-swap' : ''
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/snpeff_anno_main.py \\
        --vcf-path ${vcf} \\
        --snpeff-dir ${snpeff_dir} \\
        --parallel \\
        --max-workers ${params.snpeffWorkers} \\
        --threads ${params.snpeffThreads} \\
        --keep-cache \\
        ${swap_flag}
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 3 · PLOT_SNPEFF_STATS
 *   QC visualisation of the snpEff impact / effect distribution.
 *     in  : stats TSV (from STEP 2)
 *     out : distribution plot (PDF)
 *     via : scripts/plot_snpeff_stats.py
 * ------------------------------------------------------------------------------------- */
process PLOT_SNPEFF_STATS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${params.project}"

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    path stats_file

    output:
    path("*.pdf")

    script:
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/plot_snpeff_stats.py \\
        --input ${stats_file} \\
        --output snpeff_impact_effect_dist
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 4 · INFO_FILTER
 *   Subset the annotated VCF to one predicted-IMPACT stratum (one task per filter row).
 *     in  : (filter_tag, filter_values) · annotated vcf (+tbi)
 *     out : (filter_tag, impact-filtered vcf (+tbi)) · filter stats JSON
 *     via : scripts/info_filter_main.py
 * ------------------------------------------------------------------------------------- */
process INFO_FILTER {
    executor 'slurm'
    queue    'gr10478b'
    time     '48h'
    tag      "${filter_tag}"

    publishDir "${params.outDir}/03.info_filter", mode: 'symlink'

    input:
    tuple val(filter_tag), val(filter_values), path(vcf), path(vcf_tbi)

    output:
    path('*.json')
    tuple val(filter_tag),
          path("${params.outPrefix}.${filter_tag}.vcf.gz"),
          path("${params.outPrefix}.${filter_tag}.vcf.gz.tbi"), emit: filtered_vcf

    script:
    def out_prefix = "${params.outPrefix}.${filter_tag}"
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/info_filter_main.py \\
        --input ${vcf} \\
        --info-key impact \\
        --values ${filter_values} \\
        --out-prefix ${out_prefix} \\
        --threads ${params.infoFilterThreads} \\
        --check-chr-prefix \\
        --no-keep-chr-prefix \\
        --verbose
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 5 · RVTEST_RUN
 *   Gene-based rare-variant association test. Fans out to one task per
 *   (impact stratum × test method) combination.
 *     in  : (method_tag, method_opt, filter_tag, vcf(+tbi), pheno, covar, refFlat)
 *     out : (filter_tag, method_tag, .assoc)
 *     via : rvtest (resolved from PATH)
 * ------------------------------------------------------------------------------------- */
process RVTEST_RUN {
    executor 'slurm'
    queue    'gr10478b'
    time     '72h'
    tag      "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/04.rvtest_run/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(method_tag), val(method_opt), val(filter_tag),
          path(vcf), path(vcf_tbi), path(pheno), path(covar), path(refflat)

    output:
    path("*.log")
    tuple val(filter_tag), val(method_tag), path("*.assoc"), emit: results

    script:
    def out_prefix = "${params.outPrefix}.${filter_tag}.${method_tag}"
    """
    export PATH=/home/b/b37974/rvtests/executable:\$PATH
    rvtest \\
        --inVcf ${vcf} \\
        --pheno ${pheno} --pheno-name ${params.phenoName} \\
        --covar ${covar} --covar-name ${params.covarName} \\
        --geneFile ${refflat} \\
        --out ${out_prefix} \\
        --noweb \\
        --numThread ${params.rvtestThreads} \\
        ${method_opt}
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 6 · RVTEST_POST_PROCESS
 *   Filter genes by minimum NumVar and apply Benjamini-Hochberg FDR correction.
 *     in  : (filter_tag, method_tag, .assoc)
 *     out : (filter_tag, method_tag, .fdr.assoc)
 *     via : scripts/rvtest_post_process.py
 * ------------------------------------------------------------------------------------- */
process RVTEST_POST_PROCESS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/05.post_process/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), path(assoc_file)

    output:
    tuple val(filter_tag), val(method_tag), path("*.fdr.assoc"), emit: final_assoc

    script:
    def out_file = "${assoc_file.baseName}.filtered.fdr.assoc"
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/rvtest_post_process.py \\
        --input ${assoc_file} \\
        --output ${out_file} \\
        --num-var-threshold ${params.num_var_threshold}
    """
}

/* ----------------------------------------------------------------------------------------
 * STEP 7 · RVTEST_VISUALIZATION
 *   Manhattan + QQ plots for each post-processed association result.
 *     in  : (filter_tag, method_tag, .fdr.assoc)
 *     out : (filter_tag, method_tag, Manhattan/QQ PNG + PDF)
 *     via : scripts/plot_manhattan_qq.py
 * ------------------------------------------------------------------------------------- */
process RVTEST_VISUALIZATION {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/06.visualization/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), path(assoc_file)

    output:
    tuple val(filter_tag), val(method_tag), path("*.png"), path("*.pdf")

    script:
    def out_prefix = "${assoc_file.baseName}"
    def title      = "RVTest: ${filter_tag} - ${method_tag}"
    """
    source activate ${params.condaEnv}
    python ${params.scriptDir}/plot_manhattan_qq.py \\
        --input ${assoc_file} \\
        --output-prefix ${out_prefix} \\
        --title "${title}"
    """
}

// =========================================================================================
//  WORKFLOW
// =========================================================================================
workflow {

    // --help short-circuits the whole pipeline.
    if (params.help) {
        helpMessage()
    }
    else {
        startupSummary()

        // -- Inputs -------------------------------------------------------------------
        //    `checkIfExists` fails the run immediately (with the offending path) if any
        //    required input is missing — no half-finished submissions on a bad config.
        ch_plink = channel.of(
            tuple(
                file("${params.inputPlinkPrefix}.bed", checkIfExists: true),
                file("${params.inputPlinkPrefix}.bim", checkIfExists: true),
                file("${params.inputPlinkPrefix}.fam", checkIfExists: true)
            )
        )
        pheno_file   = file(params.phenoFile,   checkIfExists: true)
        covar_file   = file(params.covarFile,   checkIfExists: true)
        refflat_file = file(params.refFlatPath, checkIfExists: true)
        ch_snpeff    = channel.value(params.snpeffDir)

        // -- Test matrix → channels ---------------------------------------------------
        ch_impact_filters = channel.fromList( params.impactFilters.collect { f -> [ f.tag, f.values ] } )
        ch_test_methods   = channel.fromList( params.testMethods.collect   { m -> [ m.tag, m.opt    ] } )

        // -- [0] Optional genotype-matrix QC filter -----------------------------------
        if (params.filterGenotype) {
            FILTER_GENOTYPE(ch_plink,
                            channel.value(params.sampleRemoveFile),
                            channel.value(params.variantExcludeFile))
            ch_genotypes = FILTER_GENOTYPE.out.plink
        }
        else {
            ch_genotypes = ch_plink
        }

        // -- [1] Prepare genotypes, phenotypes, covariates, gene model ----------------
        RVTEST_PREPARE(ch_genotypes, pheno_file, covar_file, refflat_file)

        // -- [2] Annotate, then [3] QC-plot the annotation ----------------------------
        SNPEFF_ANNOTATE(ch_snpeff, RVTEST_PREPARE.out.vcf)
        PLOT_SNPEFF_STATS(SNPEFF_ANNOTATE.out.stats)

        // -- [4] One impact-filtered VCF per stratum ----------------------------------
        ch_info_filter_in = ch_impact_filters.combine(SNPEFF_ANNOTATE.out.annotated_vcf)
        INFO_FILTER(ch_info_filter_in)

        // -- [5] Cross (methods) × (impact-filtered VCFs), join pheno/covar + refFlat --
        ch_rvtest_in = ch_test_methods
            .combine(INFO_FILTER.out.filtered_vcf)
            .combine(RVTEST_PREPARE.out.pheno_cov)
            .combine(RVTEST_PREPARE.out.refflat)
        RVTEST_RUN(ch_rvtest_in)

        // -- [6] FDR post-processing, then [7] visualisation --------------------------
        RVTEST_POST_PROCESS(RVTEST_RUN.out.results)
        RVTEST_VISUALIZATION(RVTEST_POST_PROCESS.out.final_assoc)
    }

    // -- Completion / failure handlers (must live inside the entry workflow) ----------
    workflow.onComplete {
        log.info """
=========================================================================================
 Pipeline ${workflow.success ? 'COMPLETED SUCCESSFULLY' : 'FAILED'}
-----------------------------------------------------------------------------------------
 Duration  : ${workflow.duration}
 Exit code : ${workflow.exitStatus}
 Results   : ${params.outDir}
=========================================================================================
""".stripIndent()
    }

    workflow.onError {
        log.error "Pipeline execution stopped: ${workflow.errorMessage}"
    }
}
