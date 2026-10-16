nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// gene_signals — one workbook of the gene-based rare-variant association signals
// across the three nested cohorts, with the variants behind them annotated as in
// the common-variant table (integration/saige_signals).
//
// Integrative: it reads a finished analysis/assoc_gene run (and the tuning.rv
// call sets that run tested) and does NOT re-run any test. What it adds:
//   · the genes worth a table: every gene with a cross-cohort consistency tier,
//     plus every gene below the Bonferroni threshold in any cell;
//   · per gene x variant set x cohort: the burden the test saw (carriers,
//     cumulative MAC / MAF, multi-site carriers) beside SKAT-O and CMC results;
//   · per variant: the annotation of the common-variant table (rsID, Ensembl 86
//     gene, snpEff, Open Targets, ToMMo / 1000 Genomes / gnomAD frequencies) and,
//     per cohort, genotypes, EAF, missingness, HWE P, and whether the cohort's
//     test used the variant — with the reason when it did not;
//   · per cohort: the analysed rare variants, their impact and consequence
//     terms, and the way to the variants of the tested genes;
//   · the rare-variant QC table, and the scan summaries (λGC as observed).
// Every count is asserted equal to the upstream tables before it is written.
//
//   upstream manifest guard (before anything is submitted)
//   RV_QC_STATS × cohort ─┬► VARIANT_SETS ─────────────────────────┐
//   MAP_RECOUNT × cohort ─┘                                         │
//   SELECT_GENES ─┬─► GENOTYPE_STATS × cohort ─────────────────────┤
//                 ├─► ANNOTATE_FUNCTION ───────────────────────────┤
//                 ├─► REF_TOMMO · REF_1KG · REF_ANNOVAR ───────────┤
//                 └─► OPENTARGETS (local task) ────────────────────┴─► ASSEMBLE ─► WRITE_TABLE
//                                                                                └─► WRITE_RUN_MANIFEST
//
// Full methodology: docs/METHODS.md. Output dictionary: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// =============================================================================
// SITE CONFIGURATION — everything below is environment- or project-specific.
// Override on the command line (`--key value`) or with `-params-file cfg.yaml`.
// out_dir and run_info_dir live in nextflow.config so the run records follow
// --out_dir.
// =============================================================================
// ---- Paths
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/integration/gene_signals"
params.script_dir        = "${params.base_dir}/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.gene_dir          = "${params.project_dir}/analysis/assoc_gene/results"         // a finished assoc_gene run
params.gene_script_dir   = "${params.project_dir}/analysis/assoc_gene/scripts"         // its build_gene_map.py
params.tuning_dir        = "${params.project_dir}/tuning.rv/results"                   // the call sets it tested
params.model_inputs      = "${params.project_dir}/wgs.auto.par/results/12_model_inputs" // rare / common split, PCs
params.sample_info       = "${params.project_dir}/tmp/cram.v6/cram.v6.summary.csv"     // ID_JHRPv6, WGS_Platform, Project
params.v6_info           = "${params.project_dir}/info/cteph_agp3k.v6.20260507.xlsx"   // sex, case ages, DVT, mPAP
params.agp3k_agesex      = '/LARGE1/gr10478/pub/JHRPv1/sample/AGP3K_agesex20211202rel.tsv' // control ages (release)
params.covariate_file    = 'bbj_mainland_pc.sex.tsv'   // the sex + PC table of the association model (per cohort)

// ---- Table
params.Cohorts      = 'narrow_mainland,intermediate_mainland,full_mainland'      // upstream order, nested
params.CohortCodes  = 'N,I,F'
params.CohortLabels = 'Narrow mainland,Intermediate mainland,Full mainland'
params.CohortColors = '0571B0,008837,B2182B'
params.CohortTints  = 'EAF3F9,E9F5EE,FBEDEE'
params.Strata       = 'moderate_high,low_moderate_high'   // the upstream variant sets (nested)
// their snpEff impacts, as assoc_gene's params.Strata defines them
params.StrataImpacts = 'moderate_high=HIGH,MODERATE low_moderate_high=HIGH,MODERATE,LOW'
params.MaxTier      = 3            // genes with a consistency tier <= MaxTier in either set
params.IncludeBonferroni = true    // + genes below the Bonferroni threshold in any cell
params.SensMinNumVar = 3           // sensitivity analysis: the multiple-testing family of the genes with >= this
                                   // many variants (the primary family is upstream's MinNumVar)
params.RareMaf      = 0.01         // the control-MAF cut of the rare / common split upstream
params.AfCohort     = 'full_mainland'   // the cohort whose CTEPH / AGP3K EAF sit beside the references
params.AfBars       = true
params.TableLabel   = ''           // e.g. 'Table 3': prefixed to every sheet title and footer

// ---- Tools (absolute: an activated env is not first on PATH at this site)
params.python    = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'  // pandas, statsmodels, xlsxwriter, openpyxl
params.rscript   = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'         // EnsDb.Hsapiens.v86
params.plink2    = '/home/b/b37974/plink2'
params.plink2_test = '/home/b/b37974/plink2_alpha6/plink2'   // the build assoc_gene exported its test VCF with (LOVO)
params.rvtest_bin  = '/home/b/b37974/rvtests/executable'     // assoc_gene's rvtest (LOVO)
params.tabix     = '/home/b/b37974/htslib-1.9/tabix'
params.bcftools  = '/home/b/b37974/bcftools/bcftools'
params.perl      = '/usr/bin/perl'
params.fasta     = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

// ---- Annotation resources (those of the common-variant table)
params.snpeff_index_dir     = null    // null = the index the upstream run mapped genes with
params.tommo_vcf            = '/LARGE0/gr10478/pub/ToMMo/20250512/src/tommo-61kjpn-20250616-GRCh38-snvindel-gf-autosome-chr@@CHROM@@.vcf.gz'
params.tommo_crosscheck_vcf = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz'
params.kg_gt_vcf            = '/LARGE0/gr10478/pub/KGP/genotype/vcf_format/30x-grch38/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr@@CHROM@@.recalibrated_variants.vcf.gz'
params.kg_sites_vcf         = '/LARGE0/gr10478/pub/KGP/genotype/vcf_format/30x-grch38/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr@@CHROM@@.recalibrated_variants.annotated.vcf.gz'
params.kg_eas_bfile         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all'
params.kg_jpt_list          = '/LARGE0/gr10478/pub/KGP/sample/G1K_samples_v5.txt.JPT'
params.kg_n_eas             = 504
params.kg_n_jpt             = 104
params.annovar_dir          = '/LARGE0/gr10478/bin/annovar'
params.annovar_gnomad       = 'gnomad41_genome'
params.annovar_dbsnp        = 'avsnp151'

// ---- Open Targets Platform (queried from the launch host: compute nodes have no internet)
params.OpenTargets   = true
params.ot_api        = 'https://api.platform.opentargets.org/api/v4/graphql'
params.OtFocusTraits = 'thromb|embol|venous|pulmonary hypertension|pulmonary arter|coagul|fibrin|haemosta|hemosta'
params.OtTopN        = 3              // records named per summary cell, by posterior probability
params.OtQueryTag    = ''             // any new value (e.g. a date) re-queries the API on -resume

// ---- Labels printed in the workbook (must describe the resources above)
params.TommoLabel  = '61KJPN (2025-06-16)'
params.KgLabel     = '30x GRCh38'
params.GnomadLabel = 'v4.1 genomes'

// ---- Helpers
// Scripts are staged as path inputs, never referenced as bare param paths, so
// that a changed script invalidates its task on -resume.
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

def shared_file(String name) {
    return file("${params.shared_script_dir}/${name}", checkIfExists: true)
}

def asList(v) {
    if (v instanceof List) {
        return v.collect { x -> x.toString() }
    }
    return v.toString().split(',').collect { x -> x.trim() }.findAll { x -> x }
}

def toJson(x) {
    return groovy.json.JsonOutput.toJson(x)
}

// The one file of a cohort's fixed-model genotype directory that matches a suffix.
def modelFile(String cohort, String suffix) {
    def fs = files("${params.model_inputs}/${cohort}/genotype/fixed_model/*${suffix}")
    if (fs.size() != 1) {
        error "gene_signals | expected one *${suffix} in ${params.model_inputs}/${cohort}/genotype/fixed_model, found ${fs.size()}"
    }
    return fs[0]
}

// Everything GENOTYPE_STATS needs for one cohort.
def cohortInput(String cohort, String code, int minac) {
    def T = "${params.tuning_dir}/${cohort}"
    return tuple(cohort, code, minac,
                 file("${T}/02.callset_filter/filtered.bed", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.bim", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.fam", checkIfExists: true),
                 file("${T}/00.precheck/02_variant_qc/variants_depthdiff.exclude.txt", checkIfExists: true),
                 modelFile(cohort, '.maf_ge_threshold.variants.txt'), modelFile(cohort, '.maf_ref.afreq'))
}

// Everything MAP_RECOUNT needs for one cohort: its call set and the published gene map.
def recountInput(String cohort, int minac) {
    def T = "${params.tuning_dir}/${cohort}"
    return tuple(cohort, minac,
                 file("${T}/02.callset_filter/filtered.bed", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.bim", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.fam", checkIfExists: true),
                 file("${params.gene_dir}/${cohort}/01.map", checkIfExists: true))
}

// Everything LOVO needs for one cohort x variant set: the call set and assoc_gene's test inputs.
def lovoInput(String cohort, String stratum) {
    def T = "${params.tuning_dir}/${cohort}"
    def U = "${params.gene_dir}/${cohort}"
    return tuple(cohort, stratum,
                 file("${T}/02.callset_filter/filtered.bed", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.bim", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.fam", checkIfExists: true),
                 file("${U}/01.map/variants.${stratum}.txt", checkIfExists: true),
                 file("${U}/01.map/map.${stratum}.tsv.gz", checkIfExists: true),
                 file("${U}/00.prep/pheno_rvt.tsv", checkIfExists: true),
                 file("${U}/00.prep/covar_rvt.tsv", checkIfExists: true),
                 file("${U}/00.prep/pheno_coding.json", checkIfExists: true))
}

// Everything RV_QC_STATS needs for one cohort (the tuning.rv files of its QC).
def qcInput(String cohort, int minac) {
    def T = "${params.tuning_dir}/${cohort}"
    def ks = ([1, minac] as Set).sort()
    return tuple(cohort, ks,
                 file("${T}/00.precheck/metrics/sample_metrics.txt.gz", checkIfExists: true),
                 file("${T}/00.precheck/01_sample_qc/detect_outlier_samples.log", checkIfExists: true),
                 file("${T}/00.precheck/01_sample_qc/sample_outliers.remove.txt", checkIfExists: true),
                 file("${T}/00.precheck/02_variant_qc/variant_metrics.txt.gz", checkIfExists: true),
                 file("${T}/00.precheck/02_variant_qc/depthdiff_stats.tsv", checkIfExists: true),
                 file("${T}/02.callset_filter/filtered.fam", checkIfExists: true),
                 file("${T}/05.qc_collect/minac_recommendation.tsv", checkIfExists: true),
                 ks.collect { k -> file("${T}/03.qc_metrics/minac${k}/sample_metrics.txt.gz", checkIfExists: true) },
                 file("${params.model_inputs}/${cohort}/covariates/${params.covariate_file}", checkIfExists: true))
}


// ═════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ═════════════════════════════════════════════════════════════════════════════

/* SELECT_GENES — the genes of the table and the variants each cohort tested them on. */
process SELECT_GENES {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag 'selection'
    publishDir "${params.out_dir}/00.selection", mode: 'copy'

    input:
    path tiers
    path gene_scan
    path gene_index
    path variant_stats
    val alpha
    val min_num_var
    path script
    path lib

    output:
    path 'genes.tsv', emit: genes
    path 'variants.tsv', emit: variants
    path 'selection.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tiers ${tiers} --gene-scan ${gene_scan} --gene-index ${gene_index} \\
        --variant-stats ${variant_stats} --cohorts ${asList(params.Cohorts).join(',')} \\
        --cohort-codes ${params.CohortCodes} --strata ${params.Strata} --max-tier ${params.MaxTier} \\
        --include-bonferroni ${params.IncludeBonferroni} --alpha ${alpha} --min-num-var ${min_num_var} \\
        --sens-min-num-var ${params.SensMinNumVar} \\
        --out-genes genes.tsv --out-variants variants.tsv --out-json selection.json
    """
}

/* GENOTYPE_STATS — per cohort: genotypes, EAF, missingness and HWE P on the call set
 * the tests used, the reason a variant was not used, and the carriers of every gene
 * x set (asserted equal to the upstream counts). geno / samples files carry sample
 * IDs (analyst files). */
process GENOTYPE_STATS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}"
    publishDir "${params.out_dir}/01.genotypes", mode: 'copy'

    input:
    tuple val(cohort), val(code), val(minac), path(bed), path(bim), path(fam), path(exclude), path(common),
          path(afreq)
    path variants
    path variant_stats
    path set_stats
    path sample_info
    path script
    path lib
    path annot_lib

    output:
    path "geno.${cohort}.tsv", emit: geno
    path "carriers.${cohort}.tsv", emit: carriers
    path "samples.${cohort}.tsv", emit: samples
    path "multisite.${cohort}.tsv", emit: multisite

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --cohort ${cohort} --cohort-code ${code} \\
        --strata ${params.Strata} --bfile ${bed.baseName} --exclude-list ${exclude} --common-list ${common} \\
        --control-afreq ${afreq} --rare-maf ${params.RareMaf} --minac ${minac} --sample-info ${sample_info} \\
        --variant-stats ${variant_stats} --set-stats ${set_stats} --plink2 ${params.plink2} \\
        --annot-lib ${annot_lib} --threads ${task.ext.threads ?: 1} \\
        --out-geno geno.${cohort}.tsv --out-carriers carriers.${cohort}.tsv --out-samples samples.${cohort}.tsv \\
        --out-multisite multisite.${cohort}.tsv
    """
}

/* RV_QC_STATS — per cohort: the numbers of the rare-variant QC table, recomputed from
 * the tuning.rv outputs and checked against them. */
process RV_QC_STATS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}"
    publishDir "${params.out_dir}/05.rv_qc", mode: 'copy'

    input:
    tuple val(cohort), val(ks), path(sample_metrics), path(outlier_log), path(outliers), path(variant_metrics),
          path(depthdiff), path(fam), path(recommendation), path(minac_metrics, stageAs: 'minac?/*'), path(covariates)
    path sweep
    val n_pcs
    path script
    path lib

    output:
    path "rv_qc.${cohort}.json", emit: json

    script:
    def mm = [ks, minac_metrics instanceof List ? minac_metrics : [minac_metrics]].transpose()
        .collect { k, f -> "${k}=${f}" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohort ${cohort} --sample-metrics ${sample_metrics} --outlier-log ${outlier_log} \\
        --outliers ${outliers} --variant-metrics ${variant_metrics} --depthdiff-stats ${depthdiff} --fam ${fam} \\
        --recommendation ${recommendation} --sweep ${sweep} --minac-metrics ${mm} --covariates ${covariates} \\
        --n-pcs ${n_pcs} --out rv_qc.${cohort}.json
    """
}

/* MAP_RECOUNT — per cohort: the gene map rebuilt before the >= MinNumVar filter with the
 * upstream build_gene_map.py, asserted equal to the published map; gives the qualifying
 * variants and the genes with 1, 2, >= 3 variants per variant set. */
process MAP_RECOUNT {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}"
    publishDir "${params.out_dir}/05.rv_qc", mode: 'copy'

    input:
    tuple val(cohort), val(minac), path(bed), path(bim), path(fam), path(map_dir, stageAs: 'published_map')
    val snpeff_index
    val min_num_var
    path build_script
    path script
    path lib

    output:
    path "map_recount.${cohort}.json", emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohort ${cohort} --bfile ${bed.baseName} --minac ${minac} \\
        --snpeff-index ${snpeff_index} --build-script ${build_script} --strata ${params.StrataImpacts} \\
        --min-num-var ${min_num_var} --map-dir ${map_dir} --plink2 ${params.plink2} \\
        --threads ${task.ext.threads ?: 1} --out map_recount.${cohort}.json
    """
}

/* VARIANT_SETS — per cohort: analysed variants by impact and consequence term, and the
 * way to the variants of the tested genes, from the upstream gene map. */
process VARIANT_SETS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'variant-sets'
    publishDir "${params.out_dir}/05.rv_qc", mode: 'copy'

    input:
    val map_args
    path map_dirs, stageAs: 'map?'
    path rv_qc
    path recount
    val min_num_var
    path script
    path lib

    output:
    path 'variant_sets.terms.tsv', emit: terms
    path 'variant_sets.flow.tsv', emit: flow
    path 'variant_sets.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohorts ${asList(params.Cohorts).join(',')} --strata ${params.Strata} \\
        ${map_args} --rv-qc ${rv_qc} --recount ${recount} --min-num-var ${min_num_var} \\
        --out-terms variant_sets.terms.tsv --out-flow variant_sets.flow.tsv --out-json variant_sets.json
    """
}

/* ANNOTATE_FUNCTION — snpEff for the gene each variant was tested in, and the Ensembl
 * 86 gene (gene_annotate.R sources gene_utils.R from its own directory: both staged). */
process ANNOTATE_FUNCTION {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'snpeff-ensembl'
    publishDir "${params.out_dir}/02.annotation", mode: 'copy'

    input:
    path variants
    val snpeff_index
    path script
    path lib
    path annot_lib
    path gene_script
    path gene_utils

    output:
    path 'function.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --tabix ${params.tabix} --snpeff-index ${snpeff_index} \\
        --rscript ${params.rscript} --gene-script ${gene_script} --annot-lib ${annot_lib} --out function.tsv
    """
}

/* REF_ANNOVAR — gnomAD v4.1 genomes and dbSNP 151 in one ANNOVAR run. */
process REF_ANNOVAR {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'annovar'
    publishDir "${params.out_dir}/03.frequencies", mode: 'copy'

    input:
    path variants
    path script
    path lib

    output:
    path 'annovar.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --annovar-dir ${params.annovar_dir} --perl ${params.perl} \\
        --gnomad ${params.annovar_gnomad} --dbsnp ${params.annovar_dbsnp} --out annovar.tsv
    """
}

/* REF_TOMMO — ToMMo 61KJPN ALT frequency and record ID. */
process REF_TOMMO {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'tommo'
    publishDir "${params.out_dir}/03.frequencies", mode: 'copy'

    input:
    path variants
    path script
    path lib

    output:
    path 'tommo.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} tommo --variants ${variants} --vcf-template '${params.tommo_vcf}' \\
        --crosscheck-vcf '${params.tommo_crosscheck_vcf}' --bcftools ${params.bcftools} --tabix ${params.tabix} \\
        --fasta ${params.fasta} --out tommo.tsv
    """
}

/* REF_1KG — 1000 Genomes 30x EAS and JPT ALT frequency from genotypes. */
process REF_1KG {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag '1kg'
    publishDir "${params.out_dir}/03.frequencies", mode: 'copy'

    input:
    path variants
    path script
    path lib

    output:
    path 'kg.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} kg --variants ${variants} --gt-template '${params.kg_gt_vcf}' \\
        --sites-template '${params.kg_sites_vcf}' --eas-bfile ${params.kg_eas_bfile} \\
        --jpt-list ${params.kg_jpt_list} --expect-eas ${params.kg_n_eas} --expect-jpt ${params.kg_n_jpt} \\
        --plink2 ${params.plink2} --threads ${task.ext.threads ?: 1} --plink-memory ${task.ext.plink_mb ?: 6000} \\
        --bcftools ${params.bcftools} --tabix ${params.tabix} --fasta ${params.fasta} --out kg.tsv
    """
}

/* OPENTARGETS — Open Targets Platform records of every variant (as the common-variant
 * table). A local task: compute nodes have no internet access. The query tag is part
 * of the task hash, so -resume reuses the stored responses until it changes. */
process OPENTARGETS {
    executor 'local'
    tag 'opentargets'
    publishDir "${params.out_dir}/04.opentargets", mode: 'copy'

    input:
    path variants
    val query_tag
    path script
    path lib

    output:
    path 'opentargets.tsv', emit: tsv
    path 'opentargets_detail.tsv', emit: detail
    path 'opentargets_meta.json', emit: meta
    path 'raw', emit: raw

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # query tag: ${query_tag}
    ${params.python} ${script} --variants ${variants} --api '${params.ot_api}' --enabled ${params.OpenTargets} \\
        --focus-traits '${params.OtFocusTraits}' --top-n ${params.OtTopN} \\
        --out opentargets.tsv --out-detail opentargets_detail.tsv --out-meta opentargets_meta.json
    """
}

/* LOVO — per cohort x variant set: the selected genes re-tested with each variant left out in turn
 * (descriptive). The inputs and options are assoc_gene's (EMIT_TEST_INPUTS + RVTEST): the call set
 * exported to the test VCF with the same plink2 build, the published phenotype / covariate files,
 * --impute mean, CMC (+ cmcWald for beta / se) and SKAT-O. The same gates as upstream, and the
 * baseline (no variant left out) is asserted equal to the published result. */
process LOVO {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${stratum}"
    publishDir "${params.out_dir}/06.lovo", mode: 'copy', pattern: '{lovo.*.*.tsv,lovo.*.log}'

    input:
    tuple val(cohort), val(stratum), path(bed), path(bim), path(fam), path(variants), path(gmap), path(pheno),
          path(covar), path(coding)
    path genes
    path gene_scan
    val covar_names
    path script
    path lib

    output:
    path "lovo.${cohort}.${stratum}.tsv", emit: tsv
    path "lovo.${cohort}.${stratum}.*.log"

    script:
    def nCovar = covar_names.toString().split(',').size()
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=${params.rvtest_bin}:\$PATH

    # the test VCF, as EMIT_TEST_INPUTS writes it (every sample; only the set's variants)
    ${params.plink2_test} --bfile ${bed.baseName} --extract ${variants} --make-bed --out tested \\
        --threads ${task.ext.threads ?: 1} --memory ${task.ext.plink_mb}
    ${params.plink2_test} --bfile tested --export vcf id-paste=iid bgz --out tested \\
        --threads ${task.ext.threads ?: 1} --memory ${task.ext.plink_mb}
    ${params.tabix} -f -p vcf tested.vcf.gz
    test \$(wc -l < tested.fam) -eq \$(wc -l < ${fam})

    ${params.python} ${script} sets --cohort ${cohort} --stratum ${stratum} --genes ${genes} --map ${gmap} \\
        --variants ${variants} --bim tested.bim --vcf tested.vcf.gz --out-set lovo.set.txt --out-key lovo.key.tsv

    # pre-flight: every named covariate exists (rvtest would exit 0 with an uncovaried model)
    head -1 ${covar} | tr '\\t ' '\\n' > .covar_header
    for c in \$(echo '${covar_names}' | tr ',' ' '); do
        grep -qx "\$c" .covar_header || { echo "ERROR: covariate '\$c' is absent from ${covar}" >&2; exit 1; }
    done
    pheno_name=\$(${params.python} -c "import json;print(json.load(open('${coding}'))['pheno_col_rvtest'])")
    exp_case=\$(${params.python} -c "import json;print(json.load(open('${coding}'))['n_case'])")
    exp_ctrl=\$(${params.python} -c "import json;print(json.load(open('${coding}'))['n_control'])")

    for m in cmc skato; do
        if [ "\$m" = cmc ]; then opt='--burden cmc,cmcWald'; else opt='--kernel skato'; fi
        log=lovo.${cohort}.${stratum}.\$m.log
        rvtest --inVcf tested.vcf.gz --pheno ${pheno} --pheno-name "\$pheno_name" \\
            --covar ${covar} --covar-name ${covar_names} --setFile lovo.set.txt --out lovo.\$m \\
            --impute mean --noweb --numThread ${task.ext.threads ?: 1} \$opt > "\$log" 2>&1
        # post-flight, as upstream: the case / control counts, the imputation, the SKAT-O weights
        line=\$(grep -a -m1 -E 'Loaded [0-9]+ cases' "\$log" || true)
        test -n "\$line" || { echo "ERROR: rvtest reported no case/control count" >&2; exit 1; }
        got_case=\$(echo "\$line" | sed -E 's/.*Loaded ([0-9]+) cases.*/\\1/')
        got_ctrl=\$(echo "\$line" | sed -E 's/.*, ([0-9]+) controls.*/\\1/')
        got_miss=\$(echo "\$line" | sed -E 's/.*and ([0-9]+) missing.*/\\1/')
        test "\$got_case" -eq "\$exp_case" && test "\$got_ctrl" -eq "\$exp_ctrl" && test "\$got_miss" -eq 0 || {
            echo "ERROR: rvtest loaded \$line; expected \$exp_case cases, \$exp_ctrl controls" >&2; exit 1; }
        if grep -aq 'Impute missing genotype to mean (by default)' "\$log"; then
            echo "ERROR: --impute never reached rvtest's parser" >&2; exit 1
        fi
        grep -aq 'Impute missing genotype to mean' "\$log" || { echo "ERROR: no imputation reported" >&2; exit 1; }
        if grep -aq 'Uninitialized consolidation methods' "\$log"; then echo "ERROR: --impute not recognised" >&2; exit 1; fi
        if [ "\$m" = skato ]; then
            grep -aq 'beta1 = 1.00, beta2 = 25.00' "\$log" || { echo "ERROR: SKAT-O weights are not Beta(1,25)" >&2; exit 1; }
        fi
    done

    ${params.python} ${script} collect --cohort ${cohort} --stratum ${stratum} --key lovo.key.tsv \\
        --cmc lovo.cmc.CMC.assoc --wald lovo.cmc.CMCWald.assoc --skato lovo.skato.SkatO.assoc \\
        --gene-scan ${gene_scan} --coding ${coding} --n-covar ${nCovar} --out lovo.${cohort}.${stratum}.tsv
    """
}

/* ASSEMBLE — the long tables and every cross-check against the upstream tables. */
process ASSEMBLE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'assemble'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path genes
    path variants
    path sel_json
    path geno
    path carriers
    path samples
    path multisite
    path func
    path tommo
    path kg
    path annovar
    path opentargets
    path tables, stageAs: 'up/*'
    path tiers
    path rv_qc
    path vs_json
    path sample_info
    path v6_info
    path agp3k_agesex
    path upstream, stageAs: 'upstream_manifest.json'
    path lovo
    path script
    path lib

    output:
    path 'gene_signals.genes.long.tsv', emit: genes_long
    path 'gene_signals.variants.long.tsv', emit: variants_long
    path 'gene_signals.carriers.tsv', emit: carriers
    path 'gene_signals.scans.tsv', emit: scans
    path 'gene_signals.multisite.tsv', emit: multisite
    path 'gene_signals.lovo.tsv', emit: lovo
    path 'checks.json', emit: checks

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --genes ${genes} --variants ${variants} --selection-json ${sel_json} \\
        --geno ${geno} --carriers ${carriers} --samples ${samples} --multisite ${multisite} \\
        --multi-upstream up/multi_carriers_all.tsv --sample-info ${sample_info} --v6-info ${v6_info} \\
        --agp3k-agesex ${agp3k_agesex} --function ${func} --tommo ${tommo} \\
        --kg ${kg} --annovar ${annovar} --opentargets ${opentargets} \\
        --gene-scan up/gene_scan_all.tsv --gene-index up/gene_index_all.tsv --scan-qc up/scan_qc_all.tsv \\
        --denominators up/denominators.tsv --cohort-samples up/cohort_samples.tsv --tiers ${tiers} \\
        --rv-qc ${rv_qc} --variant-sets ${vs_json} --upstream-manifest ${upstream} \\
        --cohorts ${asList(params.Cohorts).join(',')} --cohort-codes ${params.CohortCodes} --strata ${params.Strata} \\
        --af-cohort ${params.AfCohort} --sens-min-num-var ${params.SensMinNumVar} --lovo ${lovo} \\
        --out-genes gene_signals.genes.long.tsv --out-variants gene_signals.variants.long.tsv \\
        --out-carriers gene_signals.carriers.tsv --out-scans gene_signals.scans.tsv \\
        --out-multisite gene_signals.multisite.tsv --out-lovo gene_signals.lovo.tsv --out-checks checks.json
    """
}

/* WRITE_TABLE — the workbook (presentation only, re-opened and checked). */
process WRITE_TABLE {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag 'table'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path genes_long
    path variants_long
    path scans
    path genes
    path checks
    path ot_detail
    path ot_meta
    path rv_qc
    path terms
    path flow
    path multisite
    path samples
    path lovo
    path upstream, stageAs: 'upstream_manifest.json'
    path script
    path lib

    output:
    path 'gene_signals.xlsx', emit: xlsx
    path 'workbook.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --genes-long ${genes_long} --variants-long ${variants_long} --scans ${scans} \\
        --genes ${genes} --checks ${checks} --upstream-manifest ${upstream} \\
        --opentargets-detail ${ot_detail} --opentargets-meta ${ot_meta} --rv-qc ${rv_qc} \\
        --variant-terms ${terms} --variant-flow ${flow} --multisite ${multisite} --samples ${samples} --lovo ${lovo} \\
        --cohorts ${asList(params.Cohorts).join(',')} --cohort-codes ${params.CohortCodes} \\
        --cohort-labels '${params.CohortLabels}' --cohort-colors ${params.CohortColors} \\
        --cohort-tints ${params.CohortTints} --strata ${params.Strata} --af-cohort ${params.AfCohort} \\
        --af-bars ${params.AfBars} --rare-maf ${params.RareMaf} --sens-min-num-var ${params.SensMinNumVar} \\
        --table-label '${params.TableLabel}' \\
        --tommo-label '${params.TommoLabel}' --kg-label '${params.KgLabel}' \\
        --kg-n-eas ${params.kg_n_eas} --kg-n-jpt ${params.kg_n_jpt} --gnomad-label '${params.GnomadLabel}' \\
        --out gene_signals.xlsx --out-stats workbook.json
    """
}

/* WRITE_RUN_MANIFEST — provenance, written once the workbook exists. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.run_info_dir}", mode: 'copy'

    input:
    val _done
    val manifest

    output:
    path 'run_manifest.json'

    script:
    """
    cat > run_manifest.json <<'JSON'
${manifest}
JSON
    """
}


// ═════════════════════════════════════════════════════════════════════════════
// WORKFLOW
// ═════════════════════════════════════════════════════════════════════════════
workflow {
    // ── upstream guard: nothing is submitted unless it passes ─────────────────
    def cohortList = asList(params.Cohorts)
    def strataList = asList(params.Strata)
    def upFile = file("${params.gene_dir}/_run_info/run_manifest.json")
    if (!upFile.exists()) {
        error "gene_signals | no upstream manifest at ${upFile}\n  --gene_dir must be the results directory of an assoc_gene run."
    }
    def up = new groovy.json.JsonSlurper().parseText(upFile.text)
    def need = ['component', 'cohorts', 'minac', 'strata', 'min_num_var', 'covariates', 'methods', 'alpha',
                'snpeff_index']
    def lacking = need.findAll { k -> !up.containsKey(k) }
    if (lacking) {
        error "gene_signals | ${upFile} lacks: ${lacking.join(', ')}"
    }
    if (up.component != 'assoc_gene') {
        error "gene_signals | ${upFile} was written by '${up.component}', not assoc_gene"
    }
    def notUp = cohortList.findAll { c -> !(c in up.cohorts) }
    if (notUp) {
        error "gene_signals | cohort(s) not in the upstream run: ${notUp.join(', ')} (upstream: ${up.cohorts.join(', ')})"
    }
    if (cohortList != up.cohorts.findAll { c -> c in cohortList }) {
        error "gene_signals | --Cohorts must follow the upstream order ${up.cohorts}"
    }
    if (strataList.toSet() != up.strata.toSet()) {
        error "gene_signals | --Strata ${strataList} must be the upstream variant sets ${up.strata}"
    }
    if (up.methods.toSet() != ['cmc', 'skato'].toSet()) {
        error "gene_signals | the workbook is laid out for the statistics cmc and skato; upstream ran ${up.methods}"
    }
    ['CohortCodes', 'CohortLabels', 'CohortColors', 'CohortTints'].each { k ->
        if (asList(params[k]).size() != cohortList.size()) {
            error "gene_signals | --${k} needs ${cohortList.size()} comma-separated entries"
        }
    }
    if (!(params.AfCohort in cohortList)) {
        error "gene_signals | --AfCohort ${params.AfCohort} is not among ${cohortList}"
    }
    def C = "${params.gene_dir}/_comparison"
    def snpeffIndex = (params.snpeff_index_dir ?: up.snpeff_index).toString()
    if (!file(snpeffIndex).isDirectory()) {
        error "gene_signals | snpEff index ${snpeffIndex} is missing"
    }
    log.info "gene_signals | upstream ${params.gene_dir}"
    log.info "gene_signals | cohorts ${cohortList.join(', ')} | minAC ${cohortList.collect { c -> up.minac[c] }.join(', ')} | sets ${strataList.join(', ')}"
    log.info "gene_signals | genes: tier <= ${params.MaxTier}${params.IncludeBonferroni ? ' + Bonferroni in any cell' : ''} | Open Targets ${params.OpenTargets} | out ${params.out_dir}"

    // ── STAGE 1 — selection; the rare-variant QC and the variant sets ──────────
    def lib = script_file('gs_common.py')
    def annot = shared_file('variant_annot.py')
    def byName = { a, b -> a.name <=> b.name }
    def tiers = file("${C}/robust_genes/tiers.tsv", checkIfExists: true)
    def tableNames = ['gene_scan_all', 'gene_index_all', 'variant_stats_all', 'set_stats_all', 'scan_qc_all',
                      'denominators', 'cohort_samples', 'multi_carriers_all']
    def T = tableNames.collectEntries { n -> [n, file("${C}/tables/${n}.tsv", checkIfExists: true)] }

    if ((params.SensMinNumVar as int) <= (up.min_num_var as int)) {
        error "gene_signals | --SensMinNumVar ${params.SensMinNumVar} must exceed the upstream min_num_var ${up.min_num_var}"
    }
    SELECT_GENES(tiers, T.gene_scan_all, T.gene_index_all, T.variant_stats_all, up.alpha, up.min_num_var as int,
                 script_file('select_genes.py'), lib)
    def ch_variants = SELECT_GENES.out.variants      // a value channel: broadcast to every consumer

    def nPcs = up.covariates.toString().split(',').count { x -> x.trim().toLowerCase().startsWith('pc') }
    RV_QC_STATS(channel.fromList(cohortList).map { c -> qcInput(c, up.minac[c] as int) },
                file("${params.tuning_dir}/_comparison/all_cohorts_qc_stats.tsv", checkIfExists: true), nPcs,
                script_file('rv_qc_stats.py'), lib)
    def ch_rvqc = RV_QC_STATS.out.json.collect(sort: byName)
    // the gene map directory of every cohort, staged as map<k> (k = 1.. in cohort order), so the cohort
    // travels in the argument; every file the step reads must be there
    def mapNames = ['anno_counts.tsv', 'map_qc.tsv'] + strataList.collectMany { s ->
        ["map.${s}.tsv.gz", "gene_index.${s}.tsv", "denominator.${s}.tsv"] }
    def mapDirs = cohortList.collect { c ->
        mapNames.each { n -> file("${params.gene_dir}/${c}/01.map/${n}", checkIfExists: true) }
        file("${params.gene_dir}/${c}/01.map", checkIfExists: true) }
    def mapArgs = '--map-dir ' + cohortList.withIndex().collect { c, i -> "${c}=map${i + 1}" }.join(' ')
    MAP_RECOUNT(channel.fromList(cohortList).map { c -> recountInput(c, up.minac[c] as int) }, snpeffIndex,
                up.min_num_var as int, file("${params.gene_script_dir}/build_gene_map.py", checkIfExists: true),
                script_file('map_recount.py'), lib)
    VARIANT_SETS(mapArgs, mapDirs, ch_rvqc, MAP_RECOUNT.out.json.collect(sort: byName), up.min_num_var as int,
                 script_file('variant_sets.py'), lib)

    // ── STAGE 2 — genotypes and carriers ──────────────────────────────────────
    def codes = [cohortList, asList(params.CohortCodes)].transpose().collectEntries { c, k -> [c, k] }
    GENOTYPE_STATS(
        channel.fromList(cohortList).map { c -> cohortInput(c, codes[c], up.minac[c] as int) },
        ch_variants, T.variant_stats_all, T.set_stats_all, file(params.sample_info, checkIfExists: true),
        script_file('genotype_stats.py'), lib, annot)
    def ch_geno = GENOTYPE_STATS.out.geno.collect(sort: byName)
    def ch_samples = GENOTYPE_STATS.out.samples.collect(sort: byName)

    // ── STAGE 3 — annotation and reference frequencies ──────────────────────
    ANNOTATE_FUNCTION(ch_variants, snpeffIndex, script_file('annotate_function.py'), lib, annot,
                      shared_file('gene_annotate.R'), shared_file('gene_utils.R'))
    REF_ANNOVAR(ch_variants, script_file('annovar_annotations.py'), lib)
    REF_TOMMO(ch_variants, script_file('ref_frequencies.py'), lib)
    REF_1KG(ch_variants, script_file('ref_frequencies.py'), lib)
    OPENTARGETS(ch_variants, params.OtQueryTag.toString(), script_file('opentargets.py'), lib)

    // ── STAGE 3b — leave-one-variant-out re-tests of the selected genes ───────
    LOVO(channel.fromList(cohortList.collectMany { c -> strataList.collect { s -> lovoInput(c, s) } }),
         SELECT_GENES.out.genes, T.gene_scan_all, up.covariates.toString(), script_file('lovo.py'), lib)

    // ── STAGE 4 — the long tables and the workbook ─────────────────────────
    ASSEMBLE(SELECT_GENES.out.genes, ch_variants, SELECT_GENES.out.json, ch_geno,
             GENOTYPE_STATS.out.carriers.collect(sort: byName), ch_samples,
             GENOTYPE_STATS.out.multisite.collect(sort: byName),
             ANNOTATE_FUNCTION.out.tsv, REF_TOMMO.out.tsv, REF_1KG.out.tsv, REF_ANNOVAR.out.tsv, OPENTARGETS.out.tsv,
             tableNames.collect { n -> T[n] }, tiers, ch_rvqc, VARIANT_SETS.out.json,
             file(params.sample_info, checkIfExists: true), file(params.v6_info, checkIfExists: true),
             file(params.agp3k_agesex, checkIfExists: true), upFile, LOVO.out.tsv.collect(sort: byName),
             script_file('assemble_table.py'), lib)
    WRITE_TABLE(ASSEMBLE.out.genes_long, ASSEMBLE.out.variants_long, ASSEMBLE.out.scans, SELECT_GENES.out.genes,
                ASSEMBLE.out.checks, OPENTARGETS.out.detail, OPENTARGETS.out.meta, ch_rvqc, VARIANT_SETS.out.terms,
                VARIANT_SETS.out.flow, ASSEMBLE.out.multisite, ch_samples, ASSEMBLE.out.lovo, upFile,
                script_file('write_table.py'), lib)

    def manifest = """{
  "component": "gene_signals",
  "session": ${toJson(workflow.sessionId.toString())},
  "upstream_results": ${toJson(params.gene_dir.toString())},
  "upstream_minac": ${toJson(cohortList.collectEntries { c -> [c, up.minac[c]] })},
  "upstream_covariates": ${toJson(up.covariates)},
  "upstream_alpha": ${up.alpha},
  "upstream_min_num_var": ${up.min_num_var},
  "tuning_dir": ${toJson(params.tuning_dir.toString())},
  "model_inputs": ${toJson(params.model_inputs.toString())},
  "covariate_file": ${toJson(params.covariate_file)},
  "sample_info": ${toJson(params.sample_info.toString())},
  "v6_info": ${toJson(params.v6_info.toString())},
  "agp3k_agesex": ${toJson(params.agp3k_agesex.toString())},
  "cohorts": ${toJson(cohortList)},
  "cohort_codes": ${toJson(asList(params.CohortCodes))},
  "cohort_labels": ${toJson(asList(params.CohortLabels))},
  "cohort_colors": ${toJson(asList(params.CohortColors))},
  "cohort_tints": ${toJson(asList(params.CohortTints))},
  "strata": ${toJson(strataList)},
  "max_tier": ${params.MaxTier},
  "include_bonferroni": ${params.IncludeBonferroni},
  "sensitivity_min_num_var": ${params.SensMinNumVar},
  "rare_maf": ${params.RareMaf},
  "af_cohort": ${toJson(params.AfCohort)},
  "snpeff_index_dir": ${toJson(snpeffIndex)},
  "tommo_vcf": ${toJson(params.tommo_vcf)},
  "tommo_crosscheck_vcf": ${toJson(params.tommo_crosscheck_vcf)},
  "kg_gt_vcf": ${toJson(params.kg_gt_vcf)},
  "kg_sites_vcf": ${toJson(params.kg_sites_vcf)},
  "kg_eas_bfile": ${toJson(params.kg_eas_bfile)},
  "kg_jpt_list": ${toJson(params.kg_jpt_list)},
  "kg_n_eas": ${params.kg_n_eas},
  "kg_n_jpt": ${params.kg_n_jpt},
  "annovar_dir": ${toJson(params.annovar_dir)},
  "annovar_protocols": ${toJson([params.annovar_gnomad, params.annovar_dbsnp])},
  "strata_impacts": ${toJson(params.StrataImpacts)},
  "gene_map_script": ${toJson("${params.gene_script_dir}/build_gene_map.py".toString())},
  "opentargets": ${params.OpenTargets},
  "ot_api": ${toJson(params.ot_api)},
  "ot_focus_traits": ${toJson(params.OtFocusTraits)},
  "ot_top_n": ${params.OtTopN},
  "ot_query_tag": ${toJson(params.OtQueryTag.toString())},
  "fasta": ${toJson(params.fasta)},
  "python": ${toJson(params.python)},
  "rscript": ${toJson(params.rscript)},
  "plink2": ${toJson(params.plink2)},
  "plink2_test_vcf": ${toJson(params.plink2_test)},
  "rvtest": ${toJson(params.rvtest_bin)},
  "tabix": ${toJson(params.tabix)},
  "bcftools": ${toJson(params.bcftools)},
  "table_label": ${toJson(params.TableLabel)}
}"""
    WRITE_RUN_MANIFEST(WRITE_TABLE.out.xlsx.map { _f -> 'done' }, channel.value(manifest))
}
