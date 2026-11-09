nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// hla_records — one workbook recording the HLA typing (hla.typing.v2) and the HLA
// association (analysis/assoc_hla.v2), in the order and with the definitions,
// terms and numbers of their two decks.
//
// Integrative: it reads finished runs and the decks' data and re-runs no test.
// What it adds:
//   · Typing part (T1–T8): the run and the genomes; DRB3/4/5 against DRB1;
//     completeness, read evidence and tiers per gene; the comparison with ToMMo
//     61KJPN and 1KG JPT; the marker QC record; P groups; every allele and
//     residue column per cohort.
//   · Association part (A1–A12): design and models; thresholds and their exact
//     max-T benchmark; Q1–Q3 (association, conditioning, the GWAS lead); λGC;
//     the DQB1 57 V summary; the DVT comparisons; the reported alleles; the
//     appendix; every single-marker and whole-position test.
//   · Contents (every slide mapped to its sheet) and Sources and checks.
// Every record is re-checked three ways before or after it is written: against
// the decks' data (DECK_CHECK), by an independent recomputation from the
// per-sample and per-run files (RECHECK), and by re-reading the workbook (AUDIT).
//
//   upstream guard (before anything is submitted)
//   TYPING_VERIFY ─► UPSTREAM ───────────────────────────────────────────────┐
//   TYPING_TABLES ─────────────────────────┬─► DECK_CHECK ─┐                │
//   GWAS_LAMBDA × cohort ─► ASSOC_TABLES ──┼─► RECHECK ────┼─► WRITE_WORKBOOK ─► AUDIT
//   DECK_INDEX ────────────────────────────┴───────────────┘                └─► WRITE_RUN_MANIFEST
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
params.project_dir  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir     = "${params.project_dir}/integration/hla_records"
params.script_dir   = "${params.base_dir}/scripts"
params.typing_dir   = "${params.project_dir}/hla.typing.v2"                     // the typing component (read only)
params.assoc_dir    = "${params.project_dir}/analysis/assoc_hla.v2"             // the association component (read only)
params.saige_dir    = "${params.project_dir}/analysis/assoc_saige"              // the genome-wide scan (GWAS lead, λGC)
params.given57v_dir = "${params.project_dir}/tmp/hla_dvtno_given57v"            // supplementary: without DVT given 57 V
params.model_inputs = "${params.project_dir}/wgs.auto.par/results/12_model_inputs" // sex, DVT history
params.cram_summary = "${params.project_dir}/tmp/cram.v6/cram.v6.summary.csv"   // WGS_Platform per genome
params.typing_deck  = "${params.typing_dir}/deck"                               // hla_typing.pptx, build_slides.js, stats.json
params.assoc_deck   = "${params.project_dir}/tmp/ppt/hla_assoc"                 // hla_assoc.pptx, build_slides.js, data

// ---- Workbook
params.Cohorts = 'narrow_mainland,intermediate_mainland,full_mainland'          // nested, as upstream

// ---- Tools (absolute: an activated env is not first on PATH at this site)
params.python = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'     // pandas, scipy, xlsxwriter, openpyxl, python-pptx
params.plink2 = '/home/b/b37974/plink2'                                          // GEC recomputation (A2)

// ---- Helpers
// Scripts are staged as path inputs, never referenced as bare param paths, so
// that a changed script invalidates its task on -resume.
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

def asList(v) {
    return v.toString().split(',').collect { x -> x.trim() }.findAll { x -> x }
}

def toJson(x) {
    return groovy.json.JsonOutput.toJson(x)
}


// ═════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ═════════════════════════════════════════════════════════════════════════════

/* TYPING_VERIFY — the typing component's own verify.sh, run read only (it writes nothing). */
process TYPING_VERIFY {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag 'typing verify'
    publishDir "${params.out_dir}/00.upstream", mode: 'copy'

    input:
    val stamp

    output:
    path 'typing_verify.log', emit: log

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    bash ${params.typing_dir}/verify.sh > typing_verify.log 2>&1
    """
}

/* UPSTREAM — the upstream guard in detail, and the record of every file read. */
process UPSTREAM {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag 'upstream'
    publishDir "${params.out_dir}/00.upstream", mode: 'copy'

    input:
    path typing_log
    val stamp
    path script

    output:
    path 'upstream.json', emit: json

    script:
    def S = "${params.saige_dir}/results"
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing ${params.typing_dir} --assoc ${params.assoc_dir} \\
        --typing-verify-log ${typing_log} --given57v ${params.given57v_dir} \\
        --typing-deck ${params.typing_deck} --assoc-deck ${params.assoc_deck} \\
        --cram-summary ${params.cram_summary} --model-inputs ${params.model_inputs} \\
        --gwas-sumstats ${asList(params.Cohorts).collect { c -> "${S}/${c}/01.assoc/main/sumstats.tsv" }.join(' ')} \\
        --out upstream.json
    """
}

/* TYPING_TABLES — the tables of sheets T1–T8. */
process TYPING_TABLES {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag 'typing tables'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    val stamp
    path script
    path lib

    output:
    path '01.typing', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing ${params.typing_dir} --model-inputs ${params.model_inputs} --out 01.typing
    """
}

/* GWAS_LAMBDA — λGC of the genome-wide scan in and outside the MHC, one cohort. */
process GWAS_LAMBDA {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag "${cohort}"

    input:
    tuple val(cohort), path(sumstats)
    path script

    output:
    path "gwas_lambda.${cohort}.json", emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --sumstats ${sumstats} --cohort ${cohort} --out gwas_lambda.${cohort}.json
    """
}

/* ASSOC_TABLES — the tables of sheets A1–A12 (the GEC recomputation runs plink2). */
process ASSOC_TABLES {
    executor 'slurm'
    queue "${params.queue}"
    time '1h'
    tag 'association tables'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path lambdas
    val stamp
    path script
    path lib

    output:
    path '02.association', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --assoc ${params.assoc_dir} --typing ${params.typing_dir} \\
        --model-inputs ${params.model_inputs} --cram ${params.cram_summary} --saige ${params.saige_dir} \\
        --given57v ${params.given57v_dir} --gwas-lambda ${lambdas} --plink2 ${params.plink2} --out 02.association
    """
}

/* DECK_INDEX — number, section and title of every slide of the two decks. */
process DECK_INDEX {
    executor 'slurm'
    queue "${params.queue}"
    time '15m'
    tag 'decks'
    publishDir "${params.out_dir}/03.decks", mode: 'copy'

    input:
    val stamp
    path script

    output:
    path 'decks.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing-pptx ${params.typing_deck}/hla_typing.pptx \\
        --typing-js ${params.typing_deck}/build_slides.js --assoc-pptx ${params.assoc_deck}/hla_assoc.pptx \\
        --assoc-js ${params.assoc_deck}/build_slides.js --out decks.json
    """
}

/* DECK_CHECK — every quantity shared with a slide against the deck's data. */
process DECK_CHECK {
    executor 'slurm'
    queue "${params.queue}"
    time '15m'
    tag 'deck agreement'
    publishDir "${params.out_dir}/03.decks", mode: 'copy'

    input:
    path typing, stageAs: 'in/01.typing'
    path assoc, stageAs: 'in/02.association'
    val stamp
    path script

    output:
    path 'deck_check.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing in/01.typing --assoc in/02.association \\
        --typing-stats ${params.typing_deck}/stats.json --assoc-stats ${params.assoc_deck}/stats.json \\
        --assoc-results ${params.assoc_deck}/data/results_all.json \\
        --assoc-pgroups ${params.assoc_deck}/data/pos57_pgroups.tsv --out deck_check.json
    """
}

/* RECHECK — every record recomputed by another route from the per-sample and per-run files. */
process RECHECK {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag 'recheck'
    publishDir "${params.out_dir}/04.checks", mode: 'copy'

    input:
    path typing, stageAs: 'in/01.typing'
    path assoc, stageAs: 'in/02.association'
    val stamp
    path script

    output:
    path 'recheck.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing in/01.typing --assoc in/02.association \\
        --typing-results ${params.typing_dir}/results --assoc-results ${params.assoc_dir}/results --out recheck.json
    """
}

/* WRITE_WORKBOOK — the workbook (presentation only; refused unless every check passed). */
process WRITE_WORKBOOK {
    executor 'slurm'
    queue "${params.queue}"
    time '30m'
    tag 'workbook'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path typing, stageAs: 'in/01.typing'
    path assoc, stageAs: 'in/02.association'
    path decks
    path deck_check
    path recheck
    path upstream
    path slide_map
    path script
    path lib

    output:
    path 'hla_records.xlsx', emit: xlsx
    path 'workbook.json', emit: json
    path 'workbook_cells.json.gz', emit: cells

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --typing in/01.typing --assoc in/02.association --decks ${decks} \\
        --slide-map ${slide_map} --checks ${deck_check} ${recheck} --upstream ${upstream} \\
        --out hla_records.xlsx --meta workbook.json --cells workbook_cells.json.gz
    """
}

/* AUDIT — the written workbook re-read: cells, P format, notes, text, XML. */
process AUDIT {
    executor 'slurm'
    queue "${params.queue}"
    time '1h'
    tag 'audit'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path xlsx
    path meta
    path cells
    path typing, stageAs: 'in/01.typing'
    path assoc, stageAs: 'in/02.association'
    path script
    path lib

    output:
    path 'audit.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --xlsx ${xlsx} --meta ${meta} --cells ${cells} --typing in/01.typing \\
        --assoc in/02.association --sample-ids ${params.typing_dir}/results/01.calls/genotypes.tsv --out audit.json
    """
}

/* WRITE_RUN_MANIFEST — provenance, written once the audit has passed. */
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
    def manFile = file("${params.assoc_dir}/results/_run_info/run_manifest.json")
    if (!manFile.exists()) {
        error "hla_records | no association run manifest at ${manFile}"
    }
    def man = new groovy.json.JsonSlurper().parseText(manFile.text)
    if (man.component != 'assoc_hla.v2') {
        error "hla_records | ${manFile} was written by '${man.component}', not assoc_hla.v2"
    }
    def notUp = cohortList.findAll { c -> !(c in man.cohorts) }
    if (notUp) {
        error "hla_records | cohort(s) not in the association run: ${notUp.join(', ')}"
    }
    def logs = files("${params.assoc_dir}/logs/verify.run*.log")
    if (!logs) {
        error "hla_records | no association verify log in ${params.assoc_dir}/logs"
    }
    def vlog = logs.max { f -> (f.name =~ /run(\d+)/)[0][1] as int }
    def vm = (vlog.text =~ /(\d+) passed, (\d+) failed/)
    if (!vm.find() || (vm.group(2) as int) != 0) {
        error "hla_records | ${vlog.name} is not a clean verification"
    }
    if (vlog.lastModified() < manFile.lastModified()) {
        error "hla_records | ${vlog.name} is older than the association run manifest: run its verify.sh first"
    }
    def typSummary = file("${params.typing_dir}/results/06.assoc_inputs/summary.json", checkIfExists: true)
    def deckFiles = ["${params.typing_deck}/hla_typing.pptx", "${params.typing_deck}/build_slides.js",
                     "${params.typing_deck}/stats.json", "${params.assoc_deck}/hla_assoc.pptx",
                     "${params.assoc_deck}/build_slides.js", "${params.assoc_deck}/stats.json",
                     "${params.assoc_deck}/data/results_all.json", "${params.assoc_deck}/data/pos57_pgroups.tsv"]
        .collect { p -> file(p, checkIfExists: true) }
    // a stamp of the upstream state: a new upstream run or a rebuilt deck re-runs every task on -resume
    def stamp = ([manFile, vlog, typSummary] + deckFiles).collect { f -> "${f.name}:${f.lastModified()}" }.join(';')

    def lib = [script_file('deckfmt.py'), script_file('hr_common.py'), script_file('hr_xlsx.py')]
    def byName = { a, b -> a.name <=> b.name }

    // ── STAGE 1 — upstream verification and the record of the files read ───────
    TYPING_VERIFY(stamp)
    UPSTREAM(TYPING_VERIFY.out.log, stamp, script_file('upstream.py'))

    // ── STAGE 2 — the tables of the two parts ─────────────────────────────────
    TYPING_TABLES(stamp, script_file('typing_tables.py'), lib)
    GWAS_LAMBDA(channel.fromList(cohortList.collect { c ->
                    tuple(c, file("${params.saige_dir}/results/${c}/01.assoc/main/sumstats.tsv", checkIfExists: true)) }),
                script_file('gwas_lambda.py'))
    ASSOC_TABLES(GWAS_LAMBDA.out.json.collect(sort: byName), stamp, script_file('assoc_tables.py'), lib)
    DECK_INDEX(stamp, script_file('deck_index.py'))

    // ── STAGE 3 — checks of every record, then the workbook and its audit ──────
    DECK_CHECK(TYPING_TABLES.out.dir, ASSOC_TABLES.out.dir, stamp, script_file('deck_check.py'))
    RECHECK(TYPING_TABLES.out.dir, ASSOC_TABLES.out.dir, stamp, script_file('recheck.py'))
    WRITE_WORKBOOK(TYPING_TABLES.out.dir, ASSOC_TABLES.out.dir, DECK_INDEX.out.json, DECK_CHECK.out.json,
                   RECHECK.out.json, UPSTREAM.out.json, file("${params.base_dir}/assets/slide_map.tsv", checkIfExists: true),
                   script_file('write_workbook.py'), lib)
    AUDIT(WRITE_WORKBOOK.out.xlsx, WRITE_WORKBOOK.out.json, WRITE_WORKBOOK.out.cells, TYPING_TABLES.out.dir,
          ASSOC_TABLES.out.dir, script_file('audit_workbook.py'), lib)

    def manifest = """{
  "component": "hla_records",
  "session": ${toJson(workflow.sessionId.toString())},
  "typing_dir": ${toJson(params.typing_dir.toString())},
  "assoc_dir": ${toJson(params.assoc_dir.toString())},
  "assoc_verify_log": ${toJson(vlog.name)},
  "saige_dir": ${toJson(params.saige_dir.toString())},
  "given57v_dir": ${toJson(params.given57v_dir.toString())},
  "model_inputs": ${toJson(params.model_inputs.toString())},
  "cram_summary": ${toJson(params.cram_summary.toString())},
  "typing_deck": ${toJson(params.typing_deck.toString())},
  "assoc_deck": ${toJson(params.assoc_deck.toString())},
  "cohorts": ${toJson(cohortList)},
  "upstream_stamp": ${toJson(stamp)},
  "python": ${toJson(params.python)},
  "plink2": ${toJson(params.plink2)},
  "nextflow": ${toJson(nextflow.version.toString())}
}"""
    WRITE_RUN_MANIFEST(AUDIT.out.json, manifest)
}
