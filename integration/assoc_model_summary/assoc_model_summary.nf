nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// assoc_model_summary — figures and text (for slides) on the choice of the GWAS
// model and on the main results, from two finished analysis components:
//   analysis/assoc_plink2  fixed-effect logistic regression (relatives removed)
//   analysis/assoc_saige   logistic mixed model, full GRM (relatives retained),
//                          with PCs (main) and without PCs (detect)
// Three questions, one result folder each:
//   01.fixed_vs_random  is a fixed-effect model enough?  (fixed vs GRM + PCs)
//   02.pc_effect        what do the PCs add on top of the GRM?  (facts only)
//   03.main_results     the main analysis in the three nested cohorts
// No association is re-run: the scans are read, reduced and re-plotted.
//
//   upstream guard (both manifests, before anything is submitted)
//   REDUCE_SCAN × cohort·model ─┬─► FIG_FIXED_RANDOM ─┐
//   PC_ASSOC ───────────────────┼─► FIG_PC_EFFECT ────┼─► INDEX ─► WRITE_RUN_MANIFEST
//                               └─► FIG_MAIN ─────────┘
//
// Methods: docs/METHODS.md. Files and columns: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// =============================================================================
// SITE CONFIGURATION — override with `--key value` or `-params-file cfg.yaml`.
// out_dir and run_info_dir live in nextflow.config so the run records follow
// --out_dir.
// =============================================================================
// ---- Paths
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/integration/assoc_model_summary"
params.script_dir        = "${params.base_dir}/scripts"
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.saige_dir         = "${params.project_dir}/analysis/assoc_saige/results"     // a finished assoc_saige run
params.plink2_dir        = "${params.project_dir}/analysis/assoc_plink2/results"    // a finished assoc_plink2 run

// ---- Design
params.Cohorts      = 'narrow_mainland,intermediate_mainland,full_mainland'      // upstream order, nested
params.CohortLabels = 'Narrow mainland,Intermediate mainland,Full mainland'
params.CohortColors = '0571B0,008837,B2182B'      // the project's slide colours (hex without #)
params.FixedModel   = 'additive'                  // assoc_plink2 model compared as the fixed model
params.GrmPcModel   = 'main'                      // assoc_saige model: SEX + PCs + full GRM (main analysis)
params.GrmOnlyModel = 'detect'                    // assoc_saige model: SEX + full GRM
params.NPcs         = 10

// ---- Tools
// Absolute interpreter, deliberately: `source activate` does not put an env
// first on PATH when another env's bin is exported ahead of it by the login
// profile, so a bare python3 can silently be another env's.
params.python = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'   // pandas, matplotlib, scipy, statsmodels

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

def names(x) {
    return asListOf(x).collect { f -> f.name }.join(',')
}

def toJson(x) {
    return groovy.json.JsonOutput.toJson(x)
}

def readJson(f) {
    return new groovy.json.JsonSlurper().parseText(f.text)
}

// The one sumstats file of a plink2 cohort x model folder.
def plinkSumstats(String cohort) {
    def fs = files("${params.plink2_dir}/${cohort}/01.assoc/${params.FixedModel}/*.glm.logistic")
    if (fs.size() != 1) {
        error "assoc_model_summary | expected one *.glm.logistic in ${params.plink2_dir}/${cohort}/01.assoc/${params.FixedModel}, found ${fs.size()}"
    }
    return file(fs[0].toString(), checkIfExists: true)
}

def modelMap() {
    return "fixed=${params.FixedModel},grm_pc=${params.GrmPcModel},grm_only=${params.GrmOnlyModel}"
}


// ═════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ═════════════════════════════════════════════════════════════════════════════

/* REDUCE_SCAN — one scan to its plotted columns and calibration statistics;
 * lambda_GC and the tail counts are asserted equal to the upstream scan_qc.tsv. */
process REDUCE_SCAN {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort} ${model}"
    publishDir "${params.out_dir}/00.scans", mode: 'copy', pattern: '*.{json,tsv}'

    input:
    tuple val(cohort), val(model), val(upstream_model), path(sumstats), path(scan_qc, stageAs: 'scan_qc.tsv')
    path script
    path lib

    output:
    tuple val(model), val(cohort), path("${cohort}.${model}.npz"), path("${cohort}.${model}.stats.json"),
          path("${cohort}.${model}.top.tsv"), emit: scan

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --sumstats ${sumstats} --scan-qc ${scan_qc} --cohort ${cohort} --model ${model} \\
        --upstream-model ${upstream_model} --out-prefix ${cohort}.${model}
    """
}

/* PC_ASSOC — do the PCs differ between cases and controls? (samples of the mixed-model scans) */
process PC_ASSOC {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag 'pc_assoc'
    publishDir "${params.out_dir}/02.pc_effect", mode: 'copy'

    input:
    tuple val(cohorts), path(pheno, stageAs: 'pheno_*.tsv')
    path script
    path lib

    output:
    path 'pc_case_association.tsv', emit: tsv

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cohorts ${cohorts.join(',')} --pheno ${names(pheno)} --n-pcs ${params.NPcs} \\
        --out pc_case_association.tsv
    """
}

/* FIG_FIXED_RANDOM — question 1: fixed vs random model. */
process FIG_FIXED_RANDOM {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'fixed_vs_random'
    publishDir "${params.out_dir}/01.fixed_vs_random", mode: 'copy', pattern: '{*.png,README.md,*.tsv}'

    input:
    path scans
    tuple val(cohorts), path(plink_peaks, stageAs: 'plink_peaks_*.tsv'), path(saige_peaks, stageAs: 'saige_peaks_*.tsv')
    path plot_style
    path plot_paired
    path plot_calibration
    path scan_tables
    path write_readme
    path lib

    output:
    path 'fixed_vs_random.png', emit: png
    path 'calibration.png'
    path 'README.md'
    path 'slide.md', emit: slide
    path '*.tsv'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${scan_tables} --cohorts ${cohorts.join(',')} --models fixed,grm_pc \\
        --peaks fixed=${names(plink_peaks)} --peaks grm_pc=${names(saige_peaks)} --upstream-models ${modelMap()} \\
        --out-tail tail_counts.tsv --out-maf lambda_by_maf.tsv --out-leads gw_leads.tsv
    ${params.python} ${plot_paired} --cohorts ${cohorts.join(',')} --cohort-labels '${params.CohortLabels}' \\
        --cohort-colors ${params.CohortColors} --models fixed,grm_pc \\
        --peaks fixed=${names(plink_peaks)} --peaks grm_pc=${names(saige_peaks)} --upstream-models ${modelMap()} \\
        --plot-style ${plot_style} --out-png fixed_vs_random.png --out-json fixed_vs_random.json
    ${params.python} ${plot_calibration} --cohorts ${cohorts.join(',')} --cohort-labels '${params.CohortLabels}' \\
        --cohort-colors ${params.CohortColors} --models fixed,grm_pc --tail tail_counts.tsv --maf lambda_by_maf.tsv \\
        --plot-style ${plot_style} --out-png calibration.png
    ${params.python} ${write_readme} --family fixed_vs_random --cohorts ${cohorts.join(',')} \\
        --cohort-labels '${params.CohortLabels}' --tail tail_counts.tsv --maf lambda_by_maf.tsv \\
        --fixed-leads gw_leads.tsv --out-dir .
    """
}

/* FIG_PC_EFFECT — question 2: random model with and without PCs (facts; conclusion pending). */
process FIG_PC_EFFECT {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'pc_effect'
    publishDir "${params.out_dir}/02.pc_effect", mode: 'copy', pattern: '{*.png,README.md,*.tsv,*.json}'

    input:
    path scans
    tuple val(cohorts), path(saige_peaks, stageAs: 'saige_peaks_*.tsv')
    path pc_assoc
    path lead_bymodel
    path plot_style
    path plot_paired
    path scan_tables
    path pc_effect_tables
    path write_readme
    path lib

    output:
    path 'grm_pc_vs_grm.png', emit: png
    path 'README.md'
    path 'slide.md', emit: slide
    path '*.tsv'
    path 'pc_effect.json'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${scan_tables} --cohorts ${cohorts.join(',')} --models grm_only,grm_pc \\
        --peaks grm_only=${names(saige_peaks)} --peaks grm_pc=${names(saige_peaks)} --upstream-models ${modelMap()} \\
        --out-tail tail_counts.tsv --out-maf lambda_by_maf.tsv --out-leads gw_leads.tsv
    ${params.python} ${pc_effect_tables} --cohorts ${cohorts.join(',')} --peaks ${names(saige_peaks)} \\
        --upstream-models ${modelMap()} --lead-bymodel ${lead_bymodel} --pc-assoc ${pc_assoc} \\
        --out-peaks peak_changes.tsv --out-or lead_or_stability.tsv --out-json pc_effect.json
    ${params.python} ${plot_paired} --cohorts ${cohorts.join(',')} --cohort-labels '${params.CohortLabels}' \\
        --cohort-colors ${params.CohortColors} --models grm_only,grm_pc \\
        --peaks grm_only=${names(saige_peaks)} --peaks grm_pc=${names(saige_peaks)} --upstream-models ${modelMap()} \\
        --column-titles 'grm_only=Random model: full GRM, no PCs|grm_pc=Random model: full GRM + ${params.NPcs} PCs (main)' \\
        --plot-style ${plot_style} --out-png grm_pc_vs_grm.png --out-json grm_pc_vs_grm.json
    ${params.python} ${write_readme} --family pc_effect --cohorts ${cohorts.join(',')} \\
        --cohort-labels '${params.CohortLabels}' --maf lambda_by_maf.tsv --pc-effect pc_effect.json \\
        --pc-assoc ${pc_assoc} --peak-changes peak_changes.tsv --lead-or lead_or_stability.tsv --out-dir .
    """
}

/* FIG_MAIN — question 3: the main analysis in the three cohorts. */
process FIG_MAIN {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'main_results'
    publishDir "${params.out_dir}/03.main_results", mode: 'copy', pattern: '{*.png,README.md,main_loci.tsv}'

    input:
    path scans
    tuple val(cohorts), path(saige_peaks, stageAs: 'saige_peaks_*.tsv')
    path crosscohort
    path plot_style
    path plot_paired
    path plot_main_summary
    path scan_tables
    path write_readme
    path lib

    output:
    path 'main_manhattan.png', emit: png
    path 'main_summary.png'
    path 'README.md'
    path 'slide.md', emit: slide
    path 'main_loci.tsv'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${scan_tables} --cohorts ${cohorts.join(',')} --models grm_pc \\
        --peaks grm_pc=${names(saige_peaks)} --upstream-models ${modelMap()} \\
        --out-tail tail_counts.tsv --out-maf lambda_by_maf.tsv --out-leads gw_leads.tsv
    ${params.python} ${plot_paired} --cohorts ${cohorts.join(',')} --cohort-labels '${params.CohortLabels}' \\
        --cohort-colors ${params.CohortColors} --models grm_pc --peaks grm_pc=${names(saige_peaks)} \\
        --upstream-models ${modelMap()} --plot-style ${plot_style} \\
        --out-png main_manhattan.png --out-json main_manhattan.json
    ${params.python} ${plot_main_summary} --cohorts ${cohorts.join(',')} --cohort-labels '${params.CohortLabels}' \\
        --cohort-colors ${params.CohortColors} --model grm_pc --tail tail_counts.tsv --crosscohort ${crosscohort} \\
        --plot-style ${plot_style} --out-png main_summary.png --out-loci main_loci.tsv
    ${params.python} ${write_readme} --family main_results --cohorts ${cohorts.join(',')} \\
        --cohort-labels '${params.CohortLabels}' --tail tail_counts.tsv --loci main_loci.tsv --out-dir .
    """
}

/* INDEX — results/README.md: the three questions and their slide summaries. */
process INDEX {
    executor 'local'
    tag 'index'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path s1, stageAs: 'slide_1.md'
    path s2, stageAs: 'slide_2.md'
    path s3, stageAs: 'slide_3.md'
    path write_readme
    path lib

    output:
    path 'README.md', emit: readme

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${write_readme} --family index --slides ${s1} ${s2} ${s3} --out-dir .
    """
}

/* WRITE_RUN_MANIFEST — provenance, written once the index exists. */
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
    // ── upstream guard: nothing is submitted unless both components fit ───────
    def cohortList = asList(params.Cohorts)
    def sFile = file("${params.saige_dir}/_run_info/run_manifest.json")
    def pFile = file("${params.plink2_dir}/_run_info/run_manifest.json")
    [sFile, pFile].each { f ->
        if (!f.exists()) {
            error "assoc_model_summary | no upstream manifest at ${f}"
        }
    }
    def sm = readJson(sFile)
    def pm = readJson(pFile)
    if (sm.component != 'assoc_saige' || pm.component != 'assoc_plink2') {
        error "assoc_model_summary | --saige_dir must hold an assoc_saige run and --plink2_dir an assoc_plink2 run (found ${sm.component}, ${pm.component})"
    }
    [params.GrmPcModel, params.GrmOnlyModel].each { m ->
        if (!(m in sm.models)) {
            error "assoc_model_summary | assoc_saige model '${m}' not in the upstream run (${sm.models})"
        }
    }
    if (!(params.FixedModel in pm.models)) {
        error "assoc_model_summary | assoc_plink2 model '${params.FixedModel}' not in the upstream run (${pm.models})"
    }
    def pcCovar = "SEX,${(1..params.NPcs).collect { i -> "PC${i}_AVG" }.join(',')}"
    if (sm.model_covariates[params.GrmPcModel] != pcCovar || sm.model_covariates[params.GrmOnlyModel] != 'SEX') {
        error "assoc_model_summary | assoc_saige covariates differ from SEX + ${params.NPcs} PCs / SEX: ${sm.model_covariates}"
    }
    if (pm.covariates != "SEX,PC1_AVG-PC${params.NPcs}_AVG") {
        error "assoc_model_summary | assoc_plink2 covariates '${pm.covariates}' are not SEX + ${params.NPcs} PCs"
    }
    [sm, pm].each { m ->
        def notUp = cohortList.findAll { c -> !(c in m.cohorts) }
        if (notUp) {
            error "assoc_model_summary | cohort(s) ${notUp} not in ${m.component} (${m.cohorts})"
        }
    }
    ['CohortLabels', 'CohortColors'].each { k ->
        if (asList(params[k]).size() != cohortList.size()) {
            error "assoc_model_summary | --${k} needs ${cohortList.size()} comma-separated entries"
        }
    }
    log.info "assoc_model_summary | fixed: ${params.plink2_dir} (${params.FixedModel}; ${pm.covariates})"
    log.info "assoc_model_summary | random: ${params.saige_dir} (${params.GrmPcModel} = ${sm.model_labels[params.GrmPcModel]}; ${params.GrmOnlyModel} = ${sm.model_labels[params.GrmOnlyModel]}; P ${sm.p_value_source})"
    log.info "assoc_model_summary | cohorts ${cohortList.join(', ')} | out ${params.out_dir}"

    // ── scans ─────────────────────────────────────────────────────────────────
    def lib = script_file('common.py')
    def scans = cohortList.collectMany { c ->
        [tuple(c, 'fixed', params.FixedModel, plinkSumstats(c),
               file("${params.plink2_dir}/${c}/02.scan/scan_qc.tsv", checkIfExists: true)),
         tuple(c, 'grm_pc', params.GrmPcModel,
               file("${params.saige_dir}/${c}/01.assoc/${params.GrmPcModel}/sumstats.tsv", checkIfExists: true),
               file("${params.saige_dir}/${c}/02.scan/scan_qc.tsv", checkIfExists: true)),
         tuple(c, 'grm_only', params.GrmOnlyModel,
               file("${params.saige_dir}/${c}/01.assoc/${params.GrmOnlyModel}/sumstats.tsv", checkIfExists: true),
               file("${params.saige_dir}/${c}/02.scan/scan_qc.tsv", checkIfExists: true))]
    }
    REDUCE_SCAN(channel.fromList(scans), script_file('reduce_scan.py'), lib)
    def reduced = REDUCE_SCAN.out.scan

    def plinkPeaks = cohortList.collect { c -> file("${params.plink2_dir}/${c}/03.peaks/model_peaks_annotation.tsv", checkIfExists: true) }
    def saigePeaks = cohortList.collect { c -> file("${params.saige_dir}/${c}/03.peaks/model_peaks_annotation.tsv", checkIfExists: true) }
    def style = shared_file('plot_style.py')

    PC_ASSOC(channel.value(tuple(cohortList, cohortList.collect { c ->
                 file("${params.saige_dir}/${c}/00.prep/pheno_covar.tsv", checkIfExists: true) })),
             script_file('pc_assoc.py'), lib)

    FIG_FIXED_RANDOM(
        reduced.filter { m, _c, _n, _s, _t -> m in ['fixed', 'grm_pc'] }.map { _m, _c, n, s, t -> [n, s, t] }.flatten().collect(),
        channel.value(tuple(cohortList, plinkPeaks, saigePeaks)), style,
        script_file('plot_paired.py'), script_file('plot_calibration.py'), script_file('scan_tables.py'),
        script_file('write_readme.py'), lib)

    FIG_PC_EFFECT(
        reduced.filter { m, _c, _n, _s, _t -> m in ['grm_only', 'grm_pc'] }.map { _m, _c, n, s, t -> [n, s, t] }.flatten().collect(),
        channel.value(tuple(cohortList, saigePeaks)), PC_ASSOC.out.tsv,
        file("${params.saige_dir}/_comparison/tables/lead_bymodel.tsv", checkIfExists: true), style,
        script_file('plot_paired.py'), script_file('scan_tables.py'), script_file('pc_effect_tables.py'),
        script_file('write_readme.py'), lib)

    FIG_MAIN(
        reduced.filter { m, _c, _n, _s, _t -> m == 'grm_pc' }.map { _m, _c, n, s, t -> [n, s, t] }.flatten().collect(),
        channel.value(tuple(cohortList, saigePeaks)),
        file("${params.saige_dir}/_comparison/tables/lead_crosscohort.tsv", checkIfExists: true), style,
        script_file('plot_paired.py'), script_file('plot_main_summary.py'), script_file('scan_tables.py'),
        script_file('write_readme.py'), lib)

    INDEX(FIG_FIXED_RANDOM.out.slide, FIG_PC_EFFECT.out.slide, FIG_MAIN.out.slide, script_file('write_readme.py'), lib)

    def manifest = """{
  "component": "assoc_model_summary",
  "session": ${toJson(workflow.sessionId.toString())},
  "saige_results": ${toJson(params.saige_dir.toString())},
  "plink2_results": ${toJson(params.plink2_dir.toString())},
  "cohorts": ${toJson(cohortList)},
  "cohort_labels": ${toJson(asList(params.CohortLabels))},
  "cohort_colors": ${toJson(asList(params.CohortColors))},
  "models": {"fixed": ${toJson("assoc_plink2:" + params.FixedModel)}, "grm_pc": ${toJson("assoc_saige:" + params.GrmPcModel)}, "grm_only": ${toJson("assoc_saige:" + params.GrmOnlyModel)}},
  "fixed_covariates": ${toJson(pm.covariates)},
  "saige_covariates": ${toJson(sm.model_covariates)},
  "saige_p_value_source": ${toJson(sm.p_value_source)},
  "n_pcs": ${params.NPcs},
  "python": ${toJson(params.python)}
}"""
    WRITE_RUN_MANIFEST(INDEX.out.readme.map { _f -> 'done' }, channel.value(manifest))
}
