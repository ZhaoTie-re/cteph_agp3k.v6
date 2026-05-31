nextflow.enable.dsl=2

// =========================================================================================
// PIPELINE PARAMETERS
// =========================================================================================

// Paths
params.gtPath         = '/LARGE0/gr10478/b37974/igg4/raw/format_matching/results'
params.snpeffDir      = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index'
params.scriptDir      = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/scripts'
params.outDir         = '/LARGE0/gr10478/b37974/igg4/analysis/rv_test/results'
params.refFlatPath    = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz'

// Project Settings
params.project          = 'igg4'
params.genotype_prefix  = "${params.project}_genebase_test_genotype.mac_filtered"

// Input Files (Hardcoded based on previous script logic, can be parameterized)
params.inputPlinkPrefix = "${params.gtPath}/03b.mac_filtered_genotypes/${params.genotype_prefix}"

// Options
params.performNorm       = false
params.allowRefAltSwap   = false
params.num_var_threshold = 3
params.phenoName         = 'pheno1'
params.covarName         = 'sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10'

// =========================================================================================
// PROCESS DEFINITIONS
// =========================================================================================

process RVTEST_PREPARE {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    
    publishDir "${params.outDir}/01.rvtest_prepare", mode: 'symlink'

    input:
    tuple path(bed), path(bim), path(fam)
    path gt_path
    path ref_flat

    output:
    path('*.log')
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple path("*.pheno_df.csv"), path("*.cov_df.no_age.csv"), emit: pheno_cov
    path("refFlat.hg38.nochr.txt.gz"), emit: refflat

    script:
    def bed_prefix = bed.baseName
    def pheno = "${gt_path}/05.cov_pheno/${params.project}.pheno_df.csv"
    def covar = "${gt_path}/05.cov_pheno/${params.project}.cov_df.no_age.csv"
    def norm_flag = params.performNorm ? "--norm" : ""
    
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_prepare_main.py \
        --bed-prefix ${bed_prefix} \
        --pheno-path ${pheno} \
        --covar-path ${covar} \
        --refflat-path ${ref_flat} \
        --threads 32 \
        --verbose \
        --log-file rvtest_prepare.log \
        ${norm_flag}
    """
}

process SNPEFF_ANNOTATE {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    path snpeff_dir
    tuple path(vcf), path(vcf_tbi)

    output:
    path('*.log')
    path('*.tsv'), emit: stats
    tuple path("*.snpeff.vcf.gz"), path("*.snpeff.vcf.gz.tbi"), emit: annotated_vcf

    script:
    def swap_flag = params.allowRefAltSwap ? "--allow-ref-alt-swap" : ""
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/snpeff_anno_main.py \
        --vcf-path ${vcf} \
        --snpeff-dir ${snpeff_dir} \
        --parallel \
        --max-workers 16 \
        --threads 16 \
        --keep-cache \
        ${swap_flag}
    """
}

process PLOT_SNPEFF_STATS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    path stats_file

    output:
    path("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/plot_snpeff_stats.py \
        --input ${stats_file} \
        --output snpeff_impact_effect_dist
    """
}

process INFO_FILTER {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    tag "${filter_tag}"

    publishDir "${params.outDir}/03.info_filter", mode: 'symlink'

    input:
    tuple val(filter_tag), val(filter_values), path(vcf), path(vcf_tbi)

    output:
    path('*.json')
    tuple val(filter_tag), path("${params.genotype_prefix}.${filter_tag}.vcf.gz"), path("${params.genotype_prefix}.${filter_tag}.vcf.gz.tbi"), emit: filtered_vcf

    script:
    out_prefix="${params.genotype_prefix}.${filter_tag}"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/info_filter_main.py \
        --input ${vcf} \
        --info-key impact \
        --values ${filter_values} \
        --out-prefix ${out_prefix} \
        --threads 16 \
        --check-chr-prefix \
        --no-keep-chr-prefix \
        --verbose
    """
}

process RVTEST_RUN {
    executor 'slurm'
    queue 'gr10478b'
    time '72h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/04.rvtest_run/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(method_tag), val(method_opt), val(filter_tag), path(vcf), path(vcf_tbi), path(pheno), path(covar), path(refflat)

    output:
    path("*.log")
    tuple val(filter_tag), val(method_tag), path("*.assoc"), emit: results

    script:
    def out_prefix = "${params.genotype_prefix}.${filter_tag}.${method_tag}"
    """
    export PATH=/home/b/b37974/rvtests/executable/:\$PATH
    rvtest \
        --inVcf ${vcf} \
        --pheno ${pheno} --pheno-name ${params.phenoName} \
        --covar ${covar} --covar-name ${params.covarName} \
        --geneFile ${refflat} \
        --out ${out_prefix} \
        --noweb \
        --numThread 4 \
        ${method_opt}
    """
}

process RVTEST_POST_PROCESS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/05.post_process/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), path(assoc_file)

    output:
    tuple val(filter_tag), val(method_tag), path("*.fdr.assoc"), emit: final_assoc

    script:
    def out_file = "${assoc_file.baseName}.filtered.fdr.assoc"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_post_process.py \
        --input ${assoc_file} \
        --output ${out_file} \
        --num-var-threshold ${params.num_var_threshold}
    """
}

process RVTEST_VISUALIZATION {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/06.visualization/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), path(assoc_file)

    output:
    tuple val(filter_tag), val(method_tag), path("*.png"), path("*.pdf")

    script:
    def out_prefix = "${assoc_file.baseName}"
    def title = "RVTest: ${filter_tag} - ${method_tag}"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/plot_manhattan_qq.py \
        --input ${assoc_file} \
        --output-prefix ${out_prefix} \
        --title "${title}"
    """
}

// =========================================================================================
// WORKFLOW
// =========================================================================================

workflow {
    // 1. Inputs
    // PLINK Input Channels
    input_bed = file("${params.inputPlinkPrefix}.bed")
    input_bim = file("${params.inputPlinkPrefix}.bim")
    input_fam = file("${params.inputPlinkPrefix}.fam")
    plink_input_ch = channel.of([input_bed, input_bim, input_fam])

    // Reference files
    ref_flat_ch = file(params.refFlatPath)
    gt_path_ch  = channel.value(params.gtPath)
    snpeff_dir_ch = channel.value(params.snpeffDir)

    // Define Combinations
    impact_filters = channel.fromList([
        ["impact_modifier_low_moderate_high", "MODIFIER LOW MODERATE HIGH"],
        ["impact_low_moderate_high", "LOW MODERATE HIGH"],
        ["impact_moderate_high", "MODERATE HIGH"],
        ["impact_high", "HIGH"]
    ])

    test_methods = channel.fromList([
        ["skato", "--kernel skato"],
        ["burden", "--burden cmc"]
    ])

    // 2. Prepare VCF and Covariates
    RVTEST_PREPARE(plink_input_ch, gt_path_ch, ref_flat_ch)
    
    // 3. Annotate VCF
    SNPEFF_ANNOTATE(snpeff_dir_ch, RVTEST_PREPARE.out.vcf)
    
    // 4. Plot SnpEff Stats
    PLOT_SNPEFF_STATS(SNPEFF_ANNOTATE.out.stats)

    // 5. Filter VCF by Info (Impact)
    // Combine filters with annotated VCF
    info_filter_input = impact_filters.combine(SNPEFF_ANNOTATE.out.annotated_vcf)
    INFO_FILTER(info_filter_input)

    // 6. Run RVTest
    // Combine methods with filtered VCFs and other necessary files
    // Logic: Cross product of (Methods) x (Filtered VCFs) -> Combined with (Pheno/Cov) and (RefFlat)
    
    rvtest_run_input = test_methods
        .combine(INFO_FILTER.out.filtered_vcf)
        .combine(RVTEST_PREPARE.out.pheno_cov)
        .combine(RVTEST_PREPARE.out.refflat)

    RVTEST_RUN(rvtest_run_input)

    // 7. Post Process Results
    RVTEST_POST_PROCESS(RVTEST_RUN.out.results)

    // 8. Visualization
    RVTEST_VISUALIZATION(RVTEST_POST_PROCESS.out.final_assoc)
}