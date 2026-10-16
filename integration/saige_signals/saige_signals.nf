nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// saige_signals — one academic table of the SAIGE association signals across
// the three nested cohorts, with the SuSiE credible sets they resolve to.
//
// Integrative: it reads a finished analysis/assoc_saige run and does NOT
// re-run any association scan or fine-mapping. What it adds:
//   · per-cohort statistics for every lead and credible-set variant, in every
//     cohort (not only where the variant was called);
//   · the credible-set status of every variant x cohort, WITH the reason when
//     the variant is not in a set — re-derived by refitting SuSiE on the stored
//     inputs (the upstream files record only whether a set survived);
//   · genotype counts / EAF / missingness / HWE per cases and controls;
//   · snpEff and Ensembl 86 annotation; ToMMo, 1000 Genomes and gnomAD
//     effect-allele frequencies;
//   · Open Targets Platform records per variant: ClinVar, GWAS and molQTL
//     credible sets, enhancer-to-gene predictions (queried from the launch host).
//
// The P-value source (SAIGE SPA or normal approximation) changes the peaks,
// the leads and the credible sets themselves, so it is not chosen here: the
// run aborts unless --PValue equals the upstream run's p_value_source, and
// prints how to produce matching upstream results.
//
//   upstream manifest guard (before anything is submitted)
//   SELECT_SIGNALS ─┬─► ASSOC_STATS ───────────────────────────────┐
//                   ├─► SUSIE_REFIT × locus·cohort ─► CS_STATUS ───┤
//                   ├─► GENOTYPE_STATS × cohort ───────────────────┤
//                   ├─► ANNOTATE_FUNCTION ─────────────────────────┤
//                   ├─► OPENTARGETS (local task) ──────────────────┤ (records ─► WRITE_TABLE)
//                   └─► REF_TOMMO · REF_1KG · REF_GNOMAD ──────────┴─► ASSEMBLE ─► WRITE_TABLE
//                                                                              └─► WRITE_RUN_MANIFEST
//
// Full methodology: docs/METHODS.md. Output dictionary: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// =============================================================================
// SITE CONFIGURATION — everything below is environment- or project-specific.
// Override on the command line (`--key value`) or with `-params-file cfg.yaml`;
// nothing outside this block needs to change to run the component elsewhere.
// out_dir and run_info_dir live in nextflow.config so the run records follow
// --out_dir.
// =============================================================================
// ---- Paths
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/integration/saige_signals"
params.script_dir        = "${params.base_dir}/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.saige_dir         = "${params.project_dir}/analysis/assoc_saige/results"   // a finished assoc_saige run

// ---- Design
params.PValue       = 'no_spa'       // 'no_spa' | 'spa' — must equal the upstream run's p_value_source
params.Cohorts      = 'narrow_mainland,intermediate_mainland,full_mainland'      // upstream order, nested
params.CohortCodes  = 'N,I,F'
params.CohortLabels = 'Narrow mainland,Intermediate mainland,Full mainland'   // as named in the project
// the cohort colours of the ancestry-refinement figure (hex without #)
params.CohortColors = '0571B0,008837,B2182B'
params.CohortTints  = 'EAF3F9,E9F5EE,FBEDEE'
params.AfCohort     = 'full_mainland'  // the cohort whose CTEPH/AGP3K EAF sit beside the references
params.GenotypeSet  = 'random_model'   // the genotype set the SAIGE scan ran on (relateds retained)
params.TableLabel   = ''                // e.g. 'Table 2': prefixed to the title and footer; '' = none
params.AfBars       = true
params.SusieTol     = 1e-3             // max |PIP refit - PIP stored| accepted
params.MinTommoFrac = 0.9

// ---- Tools
// Absolute interpreters, deliberately: `source activate` does not put an env
// first on PATH when another env's bin is exported ahead of it by the login
// profile, so a bare python3 / Rscript can silently be another env's.
params.python    = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'  // pandas, xlsxwriter, openpyxl
params.rscript   = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'         // susieR 0.14.2, EnsDb.Hsapiens.v86
params.plink2    = '/home/b/b37974/plink2'
params.tabix     = '/home/b/b37974/htslib-1.9/tabix'
params.bcftools  = '/home/b/b37974/bcftools/bcftools'
params.perl      = '/usr/bin/perl'
params.fasta     = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

// ---- Annotation and reference frequencies (@@CHROM@@ = 1..22)
params.snpeff_index_dir     = null    // null = the index the upstream run annotated with
params.tommo_vcf            = '/LARGE0/gr10478/pub/ToMMo/20250512/src/tommo-61kjpn-20250616-GRCh38-snvindel-gf-autosome-chr@@CHROM@@.vcf.gz'
params.tommo_crosscheck_vcf = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz'
params.kg_gt_vcf            = '/LARGE0/gr10478/pub/KGP/genotype/vcf_format/30x-grch38/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr@@CHROM@@.recalibrated_variants.vcf.gz'
params.kg_sites_vcf         = '/LARGE0/gr10478/pub/KGP/genotype/vcf_format/30x-grch38/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr@@CHROM@@.recalibrated_variants.annotated.vcf.gz'
params.kg_eas_bfile         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all'
params.kg_jpt_list          = '/LARGE0/gr10478/pub/KGP/sample/G1K_samples_v5.txt.JPT'
params.kg_n_eas             = 504
params.kg_n_jpt             = 104
params.annovar_dir          = '/LARGE0/gr10478/bin/annovar'
// allele-specific dbSNP (ANNOVAR avsnp, with its .idx): attributes the rsID where the ToMMo ID field lists
// several accessions at a multi-allelic position; '' skips the step
params.dbsnp_table          = "${params.annovar_dir}/humandb/hg38_avsnp151.txt"
// source descriptions printed in the table footnotes
params.TommoLabel  = '61KJPN (2025-06-16)'
params.KgLabel     = '30x GRCh38'
params.GnomadLabel = 'v4.1 genomes'

// ---- Open Targets Platform (GraphQL; queried by a local task on the launch host)
params.OpenTargets   = true           // false: no query; the Open Targets columns and sheet are left out
params.ot_api        = 'https://api.platform.opentargets.org/api/v4/graphql'
// GWAS traits matching this (case-insensitive) regex are counted as thrombosis-related
// (and set in bold on the Open Targets sheet); they are not ranked first
params.OtFocusTraits = 'thromb|embol|venous|pulmonary hypertension|pulmonary arter|coagul|fibrin|haemosta|hemosta'
params.OtTopN        = 3              // records named per summary cell, by posterior probability
params.OtQueryTag    = ''             // any new value (e.g. a date) re-queries the API on -resume

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

// A staged multi-file input is a single Path when it holds one file.
def asListOf(x) {
    return (x instanceof List) ? x : [x]
}

def toJson(x) {
    return groovy.json.JsonOutput.toJson(x)
}

// The one PLINK file set of a cohort's genotype directory.
def genoInput(String cohort, String modelInputs) {
    def beds = files("${modelInputs}/${cohort}/genotype/${params.GenotypeSet}/*.bed")
    if (beds.size() != 1) {
        error "saige_signals | expected one .bed in ${modelInputs}/${cohort}/genotype/${params.GenotypeSet}, found ${beds.size()}"
    }
    def pfx = beds[0].toString()[0..-5]
    return tuple(cohort, file("${pfx}.bed", checkIfExists: true), file("${pfx}.bim", checkIfExists: true),
                 file("${pfx}.fam", checkIfExists: true))
}

// The newest successful session of the upstream component (for the mismatch message only).
def lastOkSession(launch) {
    def h = file("${launch}/.nextflow/history")
    if (!h.exists()) {
        return null
    }
    def ok = h.readLines().collect { ln -> ln.split('\t') }.findAll { f -> f.size() > 5 && f[3] == 'OK' }
    return ok ? ok[-1][5] : null
}


// ═════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ═════════════════════════════════════════════════════════════════════════════

/* SELECT_SIGNALS — the loci of the cross-cohort panels, their variants, and the
 * peak each cohort called at each locus. */
process SELECT_SIGNALS {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag 'selection'
    publishDir "${params.out_dir}/00.selection", mode: 'copy'

    input:
    path cc
    path csall
    tuple val(cohorts), path(peaks, stageAs: 'peaks_*.tsv')
    path script
    path lib

    output:
    path 'loci.tsv', emit: loci
    path 'variants.tsv', emit: variants
    path 'sel.tsv', emit: sel
    path 'cs_members.tsv', emit: cs
    path 'pairs.tsv', emit: pairs
    path 'selection.json', emit: json

    script:
    def pk = [cohorts, asListOf(peaks)].transpose().collect { c, f -> "--peaks ${c}=${f}" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --lead-crosscohort ${cc} --cs-variants ${csall} ${pk} \\
        --cohorts ${cohorts.join(',')} --cohort-codes ${params.CohortCodes} \\
        --out-loci loci.tsv --out-variants variants.tsv --out-sel sel.tsv --out-cs cs_members.tsv \\
        --out-pairs pairs.tsv --out-json selection.json
    """
}

/* ASSOC_STATS — P / OR / CI / PIP per variant x cohort from the peak files, with
 * the genome-wide sumstats as fallback and cross-check. Peak ids repeat across
 * cohorts, so peak files are staged by index with a parallel key list. */
process ASSOC_STATS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'association'
    publishDir "${params.out_dir}/01.association", mode: 'copy'

    input:
    tuple val(keys), path(pips, stageAs: 'pip_*.tsv'), path(sss, stageAs: 'ss_*.tsv')
    tuple val(gcohorts), path(gwas, stageAs: 'gwas_*.tsv')
    path variants
    path pairs
    path cs
    path script
    path lib

    output:
    path 'assoc.tsv', emit: assoc

    script:
    def pk = [keys, asListOf(pips), asListOf(sss)].transpose().collect { k, p, s -> "--peak ${k} ${p} ${s}" }.join(' ')
    def gw = [gcohorts, asListOf(gwas)].transpose().collect { c, f -> "--gwas ${c}=${f}" }.join(' ')
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --pairs ${pairs} --cs-members ${cs} \\
        ${pk} ${gw} --cohorts ${gcohorts.join(',')} --out assoc.tsv
    """
}

/* SUSIE_REFIT — one per locus x cohort with a peak: refit SuSiE on the stored
 * inputs with the stored arguments, assert the stored PIPs and sets are
 * reproduced, and list every 95 % set before the purity filter. */
process SUSIE_REFIT {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}:${locus_id}"
    publishDir { "${params.out_dir}/02.finemap/${cohort}" }, mode: 'copy'

    input:
    tuple val(locus_id), val(cohort), val(peak_id), val(variants),
          path(sumstat), path(ld_matrix), path(ld_vars), path(pip), path(cs), path(susie_json)
    path script

    output:
    tuple val(locus_id), val(cohort), path("${cohort}.${locus_id}.effects.tsv"),
          path("${cohort}.${locus_id}.variants.tsv"), path("${cohort}.${locus_id}.fit.json"), emit: fit

    script:
    """
    set -euo pipefail
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    printf '%s\\n' ${variants.split(',').join(' ')} > variants.txt
    ${params.rscript} ${script} --sumstat ${sumstat} --ld-matrix ${ld_matrix} --ld-vars ${ld_vars} \\
        --pip ${pip} --cs ${cs} --json ${susie_json} --cohort ${cohort} --locus-id ${peak_id} \\
        --variants variants.txt --tol ${params.SusieTol} \\
        --out-effects ${cohort}.${locus_id}.effects.tsv --out-variants ${cohort}.${locus_id}.variants.tsv \\
        --out-fit ${cohort}.${locus_id}.fit.json
    """
}

/* CS_STATUS — every variant x cohort: in a set, below its cut, or why there is none. */
process CS_STATUS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'cs-status'
    publishDir "${params.out_dir}/02.finemap", mode: 'copy', pattern: '{assoc_status.tsv,cs_status.json}'
    publishDir "${params.out_dir}", mode: 'copy', pattern: 'saige_signals.susie_sets.tsv'

    input:
    path effects
    path fvars
    path fjson
    path assoc
    path pairs
    path cs
    path upstream, stageAs: 'upstream_manifest.json'
    path script
    path lib

    output:
    path 'assoc_status.tsv', emit: status
    path 'saige_signals.susie_sets.tsv', emit: sets
    path 'cs_status.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --effects ${effects} --fit-variants ${fvars} --fit-json ${fjson} \\
        --assoc ${assoc} --pairs ${pairs} --cs-members ${cs} --upstream-manifest ${upstream} \\
        --tol ${params.SusieTol} --out-status assoc_status.tsv --out-sets saige_signals.susie_sets.tsv \\
        --out-json cs_status.json
    """
}

/* GENOTYPE_STATS — per cohort: genotype counts, EAF, missingness, HWE for cases and controls. */
process GENOTYPE_STATS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}"
    publishDir "${params.out_dir}/03.genotypes", mode: 'copy'

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)
    path variants
    path script
    path lib
    path annot_lib

    output:
    tuple val(cohort), path("geno.${cohort}.tsv"), emit: geno

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --cohort ${cohort} --bfile ${bed.baseName} \\
        --plink2 ${params.plink2} --threads ${task.ext.threads ?: 1} --annot-lib ${annot_lib} \\
        --out geno.${cohort}.tsv
    """
}

/* ANNOTATE_FUNCTION — snpEff (JHRP v6 index) and Ensembl 86 gene. gene_annotate.R
 * sources gene_utils.R from its own directory, so both are staged side by side. */
process ANNOTATE_FUNCTION {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'annotation'
    publishDir "${params.out_dir}/04.annotation", mode: 'copy'

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

/* REF_TOMMO — ToMMo 61KJPN effect-allele frequency. */
process REF_TOMMO {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'tommo'
    publishDir "${params.out_dir}/05.frequencies", mode: 'copy'

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

/* REF_1KG — 1000 Genomes 30x EAS and JPT effect-allele frequency from genotypes. */
process REF_1KG {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag '1kg'
    publishDir "${params.out_dir}/05.frequencies", mode: 'copy'

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

/* REF_GNOMAD — gnomAD v4.1 genomes EAS and NFE through ANNOVAR (run in the task directory). */
process REF_GNOMAD {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'gnomad'
    publishDir "${params.out_dir}/05.frequencies", mode: 'copy'

    input:
    path variants
    path script
    path lib

    output:
    path 'gnomad.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} gnomad --variants ${variants} --annovar-dir ${params.annovar_dir} --perl ${params.perl} \\
        --out gnomad.tsv
    """
}

/* OPENTARGETS — Open Targets Platform records of every variant. A local task:
 * compute nodes have no internet access. The query tag is part of the task hash,
 * so -resume reuses the stored responses until it changes. */
process OPENTARGETS {
    executor 'local'
    tag 'opentargets'
    publishDir "${params.out_dir}/06.opentargets", mode: 'copy'

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

/* ASSEMBLE — the long table (one row per variant x cohort) and every cross-check
 * against the upstream tables. */
process ASSEMBLE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'assemble'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path variants
    path loci
    path sel
    path cs
    path sel_json
    path status
    path cs_json
    path geno
    path func
    path tommo
    path kg
    path gnomad
    path opentargets
    path upstream, stageAs: 'upstream_manifest.json'
    path script
    path lib

    output:
    path 'saige_signals.long.tsv', emit: long
    path 'checks.json', emit: checks

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --variants ${variants} --loci ${loci} --sel ${sel} --cs-members ${cs} \\
        --selection-json ${sel_json} --status ${status} --cs-json ${cs_json} --geno ${geno} \\
        --function ${func} --tommo ${tommo} --kg ${kg} --gnomad ${gnomad} --opentargets ${opentargets} \\
        --dbsnp '${params.dbsnp_table}' \\
        --cohorts ${asList(params.Cohorts).join(',')} --af-cohort ${params.AfCohort} \\
        --upstream-manifest ${upstream} --p-value ${params.PValue} --min-tommo-frac ${params.MinTommoFrac} \\
        --out-long saige_signals.long.tsv --out-checks checks.json
    """
}

/* WRITE_TABLE — the workbook (presentation only): sheet Signals, and sheet
 * Open Targets with every record behind the summary columns. */
process WRITE_TABLE {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag 'table'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path long_table
    path loci
    path checks
    path upstream, stageAs: 'upstream_manifest.json'
    path ot_detail
    path ot_meta
    path script
    path lib

    output:
    path 'saige_signals.xlsx', emit: xlsx

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --long ${long_table} --loci ${loci} --checks ${checks} --upstream-manifest ${upstream} \\
        --cohorts ${asList(params.Cohorts).join(',')} --cohort-codes ${params.CohortCodes} \\
        --cohort-labels '${params.CohortLabels}' --cohort-colors ${params.CohortColors} \\
        --cohort-tints ${params.CohortTints} --af-cohort ${params.AfCohort} \\
        --table-label '${params.TableLabel}' --af-bars ${params.AfBars} \\
        --tommo-label '${params.TommoLabel}' --kg-label '${params.KgLabel}' \\
        --kg-n-eas ${params.kg_n_eas} --kg-n-jpt ${params.kg_n_jpt} --gnomad-label '${params.GnomadLabel}' \\
        --opentargets-detail ${ot_detail} --opentargets-meta ${ot_meta} --out saige_signals.xlsx
    """
}

/* WRITE_RUN_MANIFEST — provenance, written once the table exists. */
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
    def upFile = file("${params.saige_dir}/_run_info/run_manifest.json")
    if (!upFile.exists()) {
        error "saige_signals | no upstream manifest at ${upFile}\n  --saige_dir must be the results directory of an assoc_saige run."
    }
    def up = new groovy.json.JsonSlurper().parseText(upFile.text)
    def need = ['component', 'cohorts', 'peak_model', 'p_value_source', 'p_value_column', 'p_genomewide',
                'p_suggestive', 'peak_flank_bp', 'n_pcs', 'loco', 'susie', 'model_inputs', 'snpeff_index_dir']
    def lacking = need.findAll { k -> !up.containsKey(k) }
    if (lacking) {
        error "saige_signals | ${upFile} lacks: ${lacking.join(', ')}"
    }
    if (up.component != 'assoc_saige') {
        error "saige_signals | ${upFile} was written by '${up.component}', not assoc_saige"
    }
    if (!(params.PValue in ['spa', 'no_spa'])) {
        error "saige_signals | --PValue must be 'spa' or 'no_spa', got '${params.PValue}'"
    }
    if (up.p_value_source != params.PValue) {
        def launch = file(params.saige_dir).parent
        def sess = lastOkSession(launch) ?: "<session id: run `nextflow log` in ${launch}>"
        def target = "${launch}/results.${params.PValue}"
        error([
            "saige_signals | P-value source mismatch: --PValue ${params.PValue}, but the upstream run used ${up.p_value_source} (${up.p_value_column}).",
            "  upstream manifest: ${upFile}",
            "  The P-value source changes the peaks, leads and credible sets themselves, so the table can only be built",
            "  from assoc_saige results made with the same --PValue. Nothing has been submitted.",
            "  Either re-run with --PValue ${up.p_value_source}, or first produce matching upstream results",
            "  (null fits and genome scans come from the cache; peaks, conditional analysis, LD, SuSiE and figures re-run):",
            "    cd ${launch}",
            "    nextflow run assoc_saige.nf -resume ${sess} --PValue ${params.PValue} \\",
            "        --out_dir ${target} --run_info_dir ${target}/_run_info",
            "  then, from ${projectDir}:",
            "    nextflow run saige_signals.nf --PValue ${params.PValue} --saige_dir ${target} \\",
            "        --out_dir ${projectDir}/results.${params.PValue}",
        ].join('\n'))
    }
    def notUp = cohortList.findAll { c -> !(c in up.cohorts) }
    if (notUp) {
        error "saige_signals | cohort(s) not in the upstream run: ${notUp.join(', ')} (upstream: ${up.cohorts.join(', ')})"
    }
    if (cohortList != up.cohorts.findAll { c -> c in cohortList }) {
        error "saige_signals | --Cohorts must follow the upstream order ${up.cohorts}"
    }
    ['CohortCodes', 'CohortLabels', 'CohortColors', 'CohortTints'].each { k ->
        if (asList(params[k]).size() != cohortList.size()) {
            error "saige_signals | --${k} needs ${cohortList.size()} comma-separated entries"
        }
    }
    if (!(params.AfCohort in cohortList)) {
        error "saige_signals | --AfCohort ${params.AfCohort} is not among ${cohortList}"
    }
    if (!file("${up.model_inputs}").isDirectory()) {
        error "saige_signals | upstream model_inputs ${up.model_inputs} is missing"
    }
    def snpeffIndex = (params.snpeff_index_dir ?: up.snpeff_index_dir).toString()
    log.info "saige_signals | upstream ${params.saige_dir}"
    log.info "saige_signals | P ${up.p_value_source} (${up.p_value_column}) | peak model ${up.peak_model} | cohorts ${cohortList.join(', ')}"
    log.info "saige_signals | SuSiE L=${up.susie.L} coverage=${up.susie.coverage} min|r|=${up.susie.min_abs_corr} | flank ${up.peak_flank_bp} bp"
    log.info "saige_signals | out ${params.out_dir}"

    // ── STAGE 1 — selection ─────────────────────────────────────────────────────
    def T = "${params.saige_dir}/_comparison/tables"
    def lib = script_file('sig_common.py')
    def annot = shared_file('variant_annot.py')
    def byName = { a, b -> a.name <=> b.name }

    SELECT_SIGNALS(
        file("${T}/lead_crosscohort.tsv", checkIfExists: true),
        file("${T}/cs_variants_all.tsv", checkIfExists: true),
        channel.value(tuple(cohortList, cohortList.collect { c ->
            file("${params.saige_dir}/${c}/03.peaks/peaks.tsv", checkIfExists: true) })),
        script_file('select_signals.py'), lib)
    def ch_variants = SELECT_SIGNALS.out.variants      // a value channel: broadcast to every consumer

    // pairs.tsv: locus_id label cohort peak_id('' = no peak) peak_tier variants
    def ch_pairs = SELECT_SIGNALS.out.pairs
        .splitCsv(header: true, sep: '\t')
        .filter { r -> r.peak_id }
        .map { r ->
            def p = "${params.saige_dir}/${r.cohort}/03.peaks/${r.peak_tier}/${r.peak_id}/${r.peak_id}"
            tuple(r.locus_id, r.cohort, r.peak_id, r.variants,
                  file("${p}.sumstat.tsv", checkIfExists: true), file("${p}.ld_matrix.tsv", checkIfExists: true),
                  file("${p}.ld_matrix.vars", checkIfExists: true), file("${p}.pip.tsv", checkIfExists: true),
                  file("${p}.cs.tsv", checkIfExists: true), file("${p}.susie.json", checkIfExists: true))
        }

    // ── STAGE 2 — per-cohort statistics and credible-set status ───────────────
    ASSOC_STATS(
        ch_pairs.map { _l, c, p, _v, ss, _ld, _lv, pip, _cs, _js -> tuple("${c}/${p}".toString(), pip, ss) }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows -> tuple(rows.collect { x -> x[0] }, rows.collect { x -> x[1] }, rows.collect { x -> x[2] }) },
        channel.value(tuple(cohortList, cohortList.collect { c ->
            file("${params.saige_dir}/${c}/01.assoc/${up.peak_model}/sumstats.tsv", checkIfExists: true) })),
        ch_variants, SELECT_SIGNALS.out.pairs, SELECT_SIGNALS.out.cs,
        script_file('assoc_stats.py'), lib)

    SUSIE_REFIT(ch_pairs, script_file('susie_cs_diagnostic.R'))

    CS_STATUS(
        SUSIE_REFIT.out.fit.map { _l, _c, e, _v, _j -> e }.collect(sort: byName),
        SUSIE_REFIT.out.fit.map { _l, _c, _e, v, _j -> v }.collect(sort: byName),
        SUSIE_REFIT.out.fit.map { _l, _c, _e, _v, j -> j }.collect(sort: byName),
        ASSOC_STATS.out.assoc, SELECT_SIGNALS.out.pairs, SELECT_SIGNALS.out.cs, upFile,
        script_file('cs_status.py'), lib)

    GENOTYPE_STATS(
        channel.fromList(cohortList).map { c -> genoInput(c, up.model_inputs.toString()) },
        ch_variants, script_file('genotype_stats.py'), lib, annot)

    // ── STAGE 3 — annotation and reference frequencies ──────────────────────
    ANNOTATE_FUNCTION(ch_variants, snpeffIndex, script_file('annotate_function.py'), lib, annot,
                      shared_file('gene_annotate.R'), shared_file('gene_utils.R'))
    REF_TOMMO(ch_variants, script_file('ref_frequencies.py'), lib)
    REF_1KG(ch_variants, script_file('ref_frequencies.py'), lib)
    REF_GNOMAD(ch_variants, script_file('ref_frequencies.py'), lib)
    OPENTARGETS(ch_variants, params.OtQueryTag.toString(), script_file('opentargets.py'), lib)

    // ── STAGE 4 — the long table and the workbook ──────────────────────────
    ASSEMBLE(ch_variants, SELECT_SIGNALS.out.loci, SELECT_SIGNALS.out.sel, SELECT_SIGNALS.out.cs,
             SELECT_SIGNALS.out.json, CS_STATUS.out.status, CS_STATUS.out.json,
             GENOTYPE_STATS.out.geno.map { _c, f -> f }.collect(sort: byName),
             ANNOTATE_FUNCTION.out.tsv, REF_TOMMO.out.tsv, REF_1KG.out.tsv, REF_GNOMAD.out.tsv,
             OPENTARGETS.out.tsv, upFile, script_file('assemble_table.py'), lib)
    WRITE_TABLE(ASSEMBLE.out.long, SELECT_SIGNALS.out.loci, ASSEMBLE.out.checks, upFile,
                OPENTARGETS.out.detail, OPENTARGETS.out.meta, script_file('write_table.py'), lib)

    def manifest = """{
  "component": "saige_signals",
  "session": ${toJson(workflow.sessionId.toString())},
  "upstream_results": ${toJson(params.saige_dir.toString())},
  "upstream_p_value_source": ${toJson(up.p_value_source)},
  "upstream_p_value_column": ${toJson(up.p_value_column)},
  "upstream_peak_model": ${toJson(up.peak_model)},
  "p_value_source": ${toJson(params.PValue)},
  "cohorts": ${toJson(cohortList)},
  "cohort_codes": ${toJson(asList(params.CohortCodes))},
  "cohort_labels": ${toJson(asList(params.CohortLabels))},
  "cohort_colors": ${toJson(asList(params.CohortColors))},
  "cohort_tints": ${toJson(asList(params.CohortTints))},
  "af_cohort": ${toJson(params.AfCohort)},
  "genotype_set": ${toJson(params.GenotypeSet)},
  "model_inputs": ${toJson(up.model_inputs.toString())},
  "susie_refit_tolerance": ${params.SusieTol},
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
  "dbsnp_table": ${toJson(params.dbsnp_table.toString())},
  "fasta": ${toJson(params.fasta)},
  "python": ${toJson(params.python)},
  "rscript": ${toJson(params.rscript)},
  "plink2": ${toJson(params.plink2)},
  "tabix": ${toJson(params.tabix)},
  "bcftools": ${toJson(params.bcftools)},
  "table_label": ${toJson(params.TableLabel)},
  "af_bars": ${params.AfBars},
  "opentargets": ${params.OpenTargets},
  "ot_api": ${toJson(params.ot_api)},
  "ot_focus_traits": ${toJson(params.OtFocusTraits)},
  "ot_top_n": ${params.OtTopN},
  "ot_query_tag": ${toJson(params.OtQueryTag.toString())},
  "ot_meta": ${toJson("${params.out_dir}/06.opentargets/opentargets_meta.json".toString())}
}"""
    WRITE_RUN_MANIFEST(WRITE_TABLE.out.xlsx.map { _f -> 'done' }, channel.value(manifest))
}
