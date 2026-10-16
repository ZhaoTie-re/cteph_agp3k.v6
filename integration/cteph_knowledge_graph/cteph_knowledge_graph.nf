nextflow.enable.dsl = 2

// ═════════════════════════════════════════════════════════════════════════════
// cteph_knowledge_graph — CTEPH, the conditions related to it, and the genetics
// of both, as a dated snapshot that can be repeated (one run per --FetchDate,
// each compared with the snapshot before it).
//
// External knowledge only: Open Targets (GraphQL), PubMed (E-utilities) and
// Europe PMC open-access full texts; nothing of this project's own results.
//   · Related conditions, each with its evidence: clinical (a curated seed of
//     guideline / study sentences, every quote checked against its PubMed
//     abstract), ontology (EFO / MONDO), and shared genetic signal (GWAS
//     colocalisation with a CTEPH locus).
//   · Genetics of CTEPH and of each condition: GWAS studies and credible sets
//     (lead variant, locus-to-gene genes, ancestry), curated / Mendelian gene
//     evidence, and the genetic studies of CTEPH in the literature with their
//     design, population, findings and methods (extracted by an isolated
//     language model, every item checked against the source text).
//   · Mechanisms: text-mined CTEPH targets by Reactome theme, with the source
//     sentences and grounded summaries; East Asian evidence throughout.
//   · A knowledge base (Markdown + JSON) and a workbook, both organised by
//     question (the literature grouped by the gene themes of
//     config/gene_themes.yaml), a graph (nodes / edges / JSON), and the
//     changes since the previous snapshot.
//
//   FETCH_CTEPH ─┬─────────────────────────────► RESOLVE ─► FETCH_DISEASE × condition ─┐
//   SEED_PUBMED ─┘                                                                       │
//   SEARCH_PUBMED ─► SCREEN ─► FETCH_FULLTEXT ─► EXTRACT ───────────────────────────────┤
//                                                    TIDY_GRAPH ◄────────────────────────┘
//                                                        ├─► SUMMARISE ─┐
//                                                        └─► DIFF ──────┴─► KNOWLEDGE ─► WRITE_TABLE ─► WRITE_RUN_MANIFEST
//
// Every process is local: compute nodes have no internet. Language-model answers
// are cached outside the work directory (llm_cache_dir), so a later snapshot
// calls the model only for new PubMed records or changed inputs.
// Full methodology: docs/METHODS.md. Output dictionary: docs/OUTPUTS.md.
// ═════════════════════════════════════════════════════════════════════════════

// =============================================================================
// SITE CONFIGURATION — override on the command line (`--key value`).
// FetchDate, out_dir and run_info_dir live in nextflow.config.
// =============================================================================
// ---- Paths
params.project_dir   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir      = "${params.project_dir}/integration/cteph_knowledge_graph"
params.script_dir    = "${params.base_dir}/scripts"
params.seed          = "${params.base_dir}/config/related_conditions.yaml"  // curated, PubMed-verified relations
params.gene_themes   = "${params.base_dir}/config/gene_themes.yaml"        // gene themes of the CTEPH literature
params.llm_cache_dir = "${params.base_dir}/cache/llm"                       // model answers, kept across snapshots
params.llm_workdir   = '/tmp/cteph_kg_llm'                                  // the model runs here, away from the project

// ---- Open Targets
params.ot_api        = 'https://api.platform.opentargets.org/api/v4/graphql'
params.DiseaseCore   = 'EFO_0022954'   // chronic thromboembolic pulmonary hypertension
params.ClppMin       = 0.01            // colocalisation with a CTEPH locus: eCAVIAR CLPP ...
params.H4Min         = 0.8             // ... or coloc H4
params.MendelianMin  = 0.5             // smallest curated / Mendelian source score kept
params.L2gMin        = 0.5             // top locus-to-gene score for a gene-condition link
params.EasMin        = 0.5             // East Asian share of a GWAS sample to flag it
params.PageSize      = 100
params.EvidencePage  = 500

// ---- PubMed: the genetic studies of CTEPH
params.PubmedQuery   = '("chronic thromboembolic pulmonary hypertension"[tiab] OR CTEPH[tiab]) AND (gene*[tiab] OR genetic*[tiab] OR polymorphism*[tiab] OR variant*[tiab] OR mutation*[tiab] OR genome-wide[tiab] OR GWAS[tiab] OR HLA[tiab] OR sequencing[tiab] OR genotyp*[tiab] OR allele*[tiab] OR exome[tiab])'

// ---- Language model (isolated `claude -p`)
params.claude        = '/home/b/b37974/.local/bin/claude'
params.LlmModel      = ''      // '' = the CLI's default model
params.LlmEffort     = ''
params.LlmParallel   = 3
params.LlmBatch      = 5       // abstracts per screening call
params.Prescreen     = true    // drop reviews / records without a genetic-study term before the model reads them
params.LlmLimit      = 0       // smoke run: the first N PubMed records only (0 = all)
params.MaxThemes     = 0       // smoke run: the N largest mechanism themes only, no methods / findings summaries
params.SentencesPerPaper = 3
params.ThemeTargets  = 12      // targets whose sentences feed a theme summary
params.MinThemeTargets = 3

// ---- Tools (absolute: an activated env is not first on PATH at this site)
params.python = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'  // pandas, xlsxwriter, openpyxl, yaml

// ---- Helpers
def script_file(String name) {
    return file("${params.script_dir}/${name}", checkIfExists: true)
}

def toJson(x) {
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(x))
}

def llmArgs() {
    return "--cache-dir '${params.llm_cache_dir}' --workdir '${params.llm_workdir}' --model '${params.LlmModel}' " +
        "--effort '${params.LlmEffort}' --parallel ${params.LlmParallel} --claude '${params.claude}'"
}


/* FETCH_CTEPH — the core disease: targets, text-mined and GWAS evidence, target pathways,
 * its GWAS credible sets with their colocalisations. The fetch date is in the task hash. */
process FETCH_CTEPH {
    executor 'local'
    tag "${params.DiseaseCore}"
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    val fetch_date
    path script
    path lib

    output:
    path '00.raw', emit: raw

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode cteph --api '${params.ot_api}' --disease ${params.DiseaseCore} \\
        --fetch-date ${fetch_date} --page-size ${params.PageSize} --evidence-page ${params.EvidencePage} --out-dir 00.raw
    """
}

/* SEED_PUBMED — the PubMed records of every PMID the curated seed cites. */
process SEED_PUBMED {
    executor 'local'
    tag 'seed'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    path seed
    val fetch_date
    path script
    path lib

    output:
    path 'seed_pubmed.json.gz', emit: records

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${fetch_date}
    ${params.python} ${script} --mode seed --seed ${seed} --out seed_pubmed.json.gz
    """
}

/* RESOLVE — the related conditions and their evidence (clinical quotes checked, ontology, colocalisation). */
process RESOLVE {
    executor 'local'
    tag 'conditions'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path raw
    path seed
    path seed_pubmed
    path script
    path lib
    path llm

    output:
    path '01.conditions', emit: dir
    path '01.conditions/fetch_list.tsv', emit: fetch_list

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --seed ${seed} --seed-pubmed ${seed_pubmed} --cteph-gwas ${raw}/gwas.json.gz \\
        --api '${params.ot_api}' --clpp-min ${params.ClppMin} --h4-min ${params.H4Min} --out-dir 01.conditions
    """
}

/* FETCH_DISEASE — one related condition's GWAS and curated / Mendelian gene evidence. */
process FETCH_DISEASE {
    executor 'local'
    maxForks 3
    tag "${disease}"
    publishDir "${params.out_dir}/00.raw/diseases", mode: 'copy'

    input:
    tuple val(disease), val(mendelian)
    val fetch_date
    path script
    path lib

    output:
    path "disease.${disease}.json.gz", emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    # fetch date: ${fetch_date}
    ${params.python} ${script} --mode disease --api '${params.ot_api}' --disease ${disease} --mendelian ${mendelian} \\
        --threshold ${params.MendelianMin} --page-size ${params.PageSize} --out-dir .
    """
}

/* SEARCH_PUBMED — the PubMed records the CTEPH-genetics query finds, up to the fetch date. */
process SEARCH_PUBMED {
    executor 'local'
    tag 'search'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    val fetch_date
    path script
    path lib

    output:
    path 'pubmed_search.json.gz', emit: records

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > query.txt <<'Q'
${params.PubmedQuery}
Q
    ${params.python} ${script} --mode search --query "\$(cat query.txt)" --max-date ${fetch_date} --out pubmed_search.json.gz
    """
}

/* SCREEN — pass A of the extraction (abstracts); lists the included records with a PMCID. */
process SCREEN {
    executor 'local'
    tag 'screen'

    input:
    path records
    path script
    path lib
    path llm

    output:
    path 'fulltext_ids.txt', emit: ids

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --records ${records} ${llmArgs()} --api '${params.ot_api}' --batch ${params.LlmBatch} \\
        --limit ${params.LlmLimit} --prescreen ${params.Prescreen} --fulltext-list fulltext_ids.txt --out-dir screen
    """
}

/* FETCH_FULLTEXT — the Methods sections of the included open-access full texts (Europe PMC). */
process FETCH_FULLTEXT {
    executor 'local'
    tag 'fulltext'
    publishDir "${params.out_dir}/00.pubmed", mode: 'copy'

    input:
    path ids
    path script
    path lib

    output:
    path 'fulltext.json.gz', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --mode fulltext --ids ${ids} --out fulltext.json.gz
    """
}

/* EXTRACT — passes A (cached) and B: the genetic studies of CTEPH, their findings and methods. */
process EXTRACT {
    executor 'local'
    tag 'extract'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path records
    path fulltext
    path script
    path lib
    path llm

    output:
    path '02.literature', emit: dir

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --records ${records} --fulltext ${fulltext} ${llmArgs()} --api '${params.ot_api}' \\
        --batch ${params.LlmBatch} --limit ${params.LlmLimit} --prescreen ${params.Prescreen} --out-dir 02.literature
    """
}

/* TIDY_GRAPH — the long tables and the graph. */
process TIDY_GRAPH {
    executor 'local'
    tag 'tables'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path raw
    path conditions
    path diseases
    path literature
    path script
    path lib

    output:
    path '03.tables', emit: tables

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --cteph-raw ${raw} --conditions ${conditions} --diseases ${diseases} \\
        --extract ${literature} --sentences ${params.SentencesPerPaper} --eas-min ${params.EasMin} \\
        --l2g-min ${params.L2gMin} --out-dir 03.tables
    """
}

/* SUMMARISE — grounded summaries: mechanism themes, methods and findings of the genetic studies. */
process SUMMARISE {
    executor 'local'
    tag 'summaries'
    publishDir "${params.out_dir}/04.summaries", mode: 'copy'

    input:
    path tables
    path script
    path lib
    path llm

    output:
    path 'summaries.tsv', emit: tsv
    path 'summaries.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} ${llmArgs()} --min-targets ${params.MinThemeTargets} \\
        --top-targets ${params.ThemeTargets} --max-themes ${params.MaxThemes} --out summaries.tsv --out-json summaries.json
    """
}

/* DIFF — the changes since the newest earlier snapshot. */
process DIFF {
    executor 'local'
    tag 'changes'
    publishDir "${params.out_dir}/05.changes", mode: 'copy'

    input:
    path tables
    val current
    val earlier
    path script
    path lib

    output:
    path 'changes.tsv', emit: tsv
    path 'changes.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} --current ${current} --snapshots-root '${params.snapshots_root}' \\
        --earlier '${earlier}' --out changes.tsv --out-json changes.json
    """
}

/* KNOWLEDGE — the knowledge base of the snapshot. */
process KNOWLEDGE {
    executor 'local'
    tag 'knowledge'
    publishDir "${params.out_dir}/06.knowledge", mode: 'copy'

    input:
    path tables
    path summaries
    path changes
    path changes_json
    path seed
    path seed_pubmed
    path themes
    path script
    path lib

    output:
    path 'cteph_knowledge_graph.md', emit: md
    path 'cteph_knowledge_graph.json', emit: json

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --seed ${seed} --seed-pubmed ${seed_pubmed} --themes ${themes} \\
        --l2g-min ${params.L2gMin} --out-md cteph_knowledge_graph.md --out-json cteph_knowledge_graph.json
    """
}

/* WRITE_TABLE — the workbook, re-opened and checked before it is published. */
process WRITE_TABLE {
    executor 'local'
    tag 'workbook'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path tables
    path summaries
    path changes
    path changes_json
    path knowledge_json
    path themes
    val workbook
    path script
    path lib

    output:
    path "${workbook}", emit: xlsx
    path 'workbook.json', emit: stats

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    ${params.python} ${script} --tables ${tables} --summaries ${summaries} --changes ${changes} \\
        --changes-json ${changes_json} --knowledge-json ${knowledge_json} --themes ${themes} \\
        --l2g-min ${params.L2gMin} --eas-min ${params.EasMin} --out ${workbook} --out-stats workbook.json
    """
}

/* WRITE_RUN_MANIFEST — provenance, and the pointer to the newest snapshot. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.run_info_dir}", mode: 'copy', pattern: 'run_manifest.json'
    publishDir "${params.results_root}", mode: 'copy', pattern: 'latest.json'

    input:
    path tables
    path workbook_stats
    path changes_json
    val run_json
    val current
    val earlier
    val workbook
    path script
    path lib

    output:
    path 'run_manifest.json'
    path 'latest.json'

    script:
    """
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    cat > run.json <<'JSON'
${run_json}
JSON
    ${params.python} ${script} --run-json run.json --tidy ${tables}/tidy.json --workbook-stats ${workbook_stats} \\
        --changes-json ${changes_json} --snapshot ${current} --earlier '${earlier}' \\
        --snapshots-root '${params.snapshots_root}' --workbook ${workbook} --out run_manifest.json \\
        --out-latest latest.json
    """
}


// ═════════════════════════════════════════════════════════════════════════════
// WORKFLOW
// ═════════════════════════════════════════════════════════════════════════════
workflow {
    def fetchDate = params.FetchDate ? params.FetchDate.toString() : ''
    if (!(fetchDate ==~ /\d{4}-\d{2}-\d{2}/)) {
        error "cteph_knowledge_graph | --FetchDate YYYY-MM-DD is required (the snapshot name), got '${fetchDate}'"
    }
    if (!params.out_dir.toString().endsWith("/snapshots/${fetchDate}")) {
        error "cteph_knowledge_graph | out_dir ${params.out_dir} is not the snapshot of ${fetchDate}"
    }
    def smoke = (params.LlmLimit as int) > 0 || (params.MaxThemes as int) > 0
    if (smoke && !fetchDate.endsWith('-smoke') && !params.out_dir.toString().contains('smoke')) {
        log.warn "cteph_knowledge_graph | smoke settings (LlmLimit ${params.LlmLimit}, MaxThemes ${params.MaxThemes}): this snapshot is partial"
    }
    def root = file(params.snapshots_root)
    def earlier = root.exists()
        ? root.listFiles().findAll { d -> d.isDirectory() && d.name != fetchDate && file("${d}/_run_info/run_manifest.json").exists() }
            .collect { d -> d.name }.sort()
        : []
    def earlierCsv = earlier.join(',')
    def workbook = "cteph_knowledge_graph_${fetchDate}.xlsx"
    log.info "cteph_knowledge_graph | snapshot ${fetchDate} | earlier snapshots: ${earlier ? earlierCsv : 'none'} | out ${params.out_dir}"

    def lib = script_file('ot_common.py')
    def llm = script_file('llm.py')
    def seed = file(params.seed, checkIfExists: true)
    def themes = file(params.gene_themes, checkIfExists: true)
    FETCH_CTEPH(fetchDate, script_file('fetch.py'), lib)
    SEED_PUBMED(seed, fetchDate, script_file('pubmed.py'), lib)
    RESOLVE(FETCH_CTEPH.out.raw, seed, SEED_PUBMED.out.records, script_file('resolve.py'), lib, llm)
    def fetchList = RESOLVE.out.fetch_list.splitCsv(sep: '\t').map { r -> [r[0], r[1]] }
    FETCH_DISEASE(fetchList, fetchDate, script_file('fetch.py'), lib)
    SEARCH_PUBMED(fetchDate, script_file('pubmed.py'), lib)
    SCREEN(SEARCH_PUBMED.out.records, script_file('extract_genetics.py'), lib, llm)
    FETCH_FULLTEXT(SCREEN.out.ids, script_file('pubmed.py'), lib)
    EXTRACT(SEARCH_PUBMED.out.records, FETCH_FULLTEXT.out.json, script_file('extract_genetics.py'), lib, llm)
    TIDY_GRAPH(FETCH_CTEPH.out.raw, RESOLVE.out.dir, FETCH_DISEASE.out.json.collect(), EXTRACT.out.dir,
               script_file('tidy_graph.py'), lib)
    SUMMARISE(TIDY_GRAPH.out.tables, script_file('summarise.py'), lib, llm)
    DIFF(TIDY_GRAPH.out.tables, fetchDate, earlierCsv, script_file('diff_snapshots.py'), lib)
    KNOWLEDGE(TIDY_GRAPH.out.tables, SUMMARISE.out.tsv, DIFF.out.tsv, DIFF.out.json, seed, SEED_PUBMED.out.records,
              themes, script_file('build_knowledge.py'), lib)
    WRITE_TABLE(TIDY_GRAPH.out.tables, SUMMARISE.out.tsv, DIFF.out.tsv, DIFF.out.json, KNOWLEDGE.out.json, themes,
                workbook, script_file('write_table.py'), lib)

    def scripts = ['ot_common.py', 'llm.py', 'fetch.py', 'pubmed.py', 'resolve.py', 'extract_genetics.py',
                   'tidy_graph.py', 'summarise.py', 'diff_snapshots.py', 'build_knowledge.py', 'write_table.py',
                   'manifest.py']
    def run = [
        launched  : workflow.start.toString(),
        session   : workflow.sessionId.toString(),
        command   : workflow.commandLine,
        nextflow  : workflow.nextflow.version.toString(),
        params    : [FetchDate: fetchDate, DiseaseCore: params.DiseaseCore, ClppMin: params.ClppMin, H4Min: params.H4Min,
                     MendelianMin: params.MendelianMin, L2gMin: params.L2gMin, EasMin: params.EasMin,
                     PubmedQuery: params.PubmedQuery, LlmModel: params.LlmModel, LlmEffort: params.LlmEffort,
                     LlmBatch: params.LlmBatch, LlmLimit: params.LlmLimit, MaxThemes: params.MaxThemes,
                     Prescreen: params.Prescreen,
                     SentencesPerPaper: params.SentencesPerPaper, ThemeTargets: params.ThemeTargets,
                     MinThemeTargets: params.MinThemeTargets, seed_md5: seed.text.md5(), gene_themes_md5: themes.text.md5()],
        smoke     : smoke,
        earlier   : earlier,
        script_md5: scripts.collectEntries { s -> [(s): script_file(s).text.md5()] },
    ]
    WRITE_RUN_MANIFEST(TIDY_GRAPH.out.tables, WRITE_TABLE.out.stats, DIFF.out.json, toJson(run), fetchDate, earlierCsv,
                       workbook, script_file('manifest.py'), lib)
}
