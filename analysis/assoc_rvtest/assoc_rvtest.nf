#!/usr/bin/env nextflow
// =============================================================================
// assoc_rvtest — gene-based rare-variant association (RVTESTS), three cohorts.
//
// Reads the callset and the minor-allele-count threshold that `tuning.rv` has
// already decided, runs SKAT-O and CMC over impact strata, tiers the genes by
// Benjamini-Hochberg, and follows up every hit gene with a per-gene deep dive.
//
// WHAT THIS COMPONENT DOES NOT DECIDE
//   The design is fully confounded — cases are 30x (DNBSeq/NovaSeq), controls
//   15x (HiSeqX), with zero platform overlap. `tuning.rv` quantifies that and
//   returns the minAC at which the ADJUSTED apparent group effect reaches a
//   statistical null. Inheriting its threshold is what makes a burden test here
//   interpretable at all, and a hit is still a CANDIDATE, not a confirmed
//   association. See docs/METHODS.md.
//
// PARALLELISM: everything from PREPARE to the association runs PER CHROMOSOME.
//   Genes never span a chromosome and rvtest tests genes independently, so the
//   split is exact rather than an approximation — the concatenated result is
//   bit-identical to a whole-genome run. The ONE thing that cannot be split is
//   the Benjamini-Hochberg correction, which is a property of the whole scan;
//   MERGE_ASSOC therefore sits between the per-chromosome tests and GENE_SCAN,
//   and the FDR is computed once, genome-wide.
//   Measured: the largest chromosome carries 9.0 % of variants, so the critical
//   path drops from ~203 min to ~19 min per cohort.
//
// Tiers, and why they are FDR rather than a fixed P
//   significant  BH FDR < params.FdrSignificant (0.05)
//   suggestive   BH FDR < params.FdrSuggestive  (0.10)
//   A scan tests 164 genes under HIGH and ~8.6 k under MODERATE+HIGH, so a
//   fixed P would mean a different multiple-testing burden per stratum. BH
//   adapts; the price is that the IMPLIED P differs between strata (measured:
//   1.3e-4 vs 4.5e-5), which every figure states explicitly.
//
// Layout
//   assoc_rvtest.nf   parameters + processes + workflow (this file)
//   nextflow.config   SLURM `--rsc` requests only (this site rejects cpus/memory)
//   scripts/          engine-specific scripts
//   ../_shared/scripts/  style, figure documents, figure catalogues
// =============================================================================

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
//  PARAMETERS
// -----------------------------------------------------------------------------
params.project_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir          = "${params.project_dir}/analysis/assoc_rvtest"
params.script_dir        = "${params.base_dir}/scripts"
// Shared with the other association components — style, figure documents and the
// family catalogue. They describe association RESULTS, not the engine that made
// them, so they live once at analysis/_shared.
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.out_dir           = "${params.base_dir}/results"
params.model_inputs      = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"

// ---- Cohorts and chromosomes ------------------------------------------------
params.Cohorts      = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.ValidCohorts = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
// The unit of parallelism: BALANCED BINS of whole chromosomes, not raw
// chromosomes. Any grouping of whole chromosomes is exact, because no gene spans
// one; so the only question is how many jobs to pay for. chr2 alone carries
// 9.05 % of variants and is therefore the wall-clock floor no split can beat.
// Twelve bins reach exactly that floor (max bin 9.05 %) while running 216 RVTEST
// jobs instead of the 396 that one-bin-per-chromosome would need — same speed,
// 45 % fewer jobs competing for the queue. Packed longest-first from the
// full_mainland call set; the ordering is a genome property and is stable.
params.ChromBins    = ['2', '3,22', '6,20', '7,18', '8,17', '5,19',
                       '4,21', '11,14', '9,13', '12,16', '10,15', '1']

// ---- The tuning.rv contract -------------------------------------------------
// tuning.rv owns sample QC, depth-differential variant QC and the minAC choice.
// This component consumes its output rather than re-deriving any of it.
params.TuningDir    = "${params.project_dir}/tuning.rv/results"
params.CallsetTpl   = "${params.TuningDir}/@@COHORT@@/02.callset_filter/filtered"
params.MinACTpl     = "${params.TuningDir}/@@COHORT@@/05.qc_collect/minac_recommendation.tsv"
// null = inherit the recommendation, per cohort. A number overrides it for every
// cohort; do that only to answer a "what if" question, never for the headline.
params.MinAC        = null

// ---- Phenotype / covariates -------------------------------------------------
params.PhenoTpl     = "${params.model_inputs}/@@COHORT@@/phenotype/pheno.tsv"
params.CovarTpl     = "${params.model_inputs}/@@COHORT@@/covariates/bbj_mainland_pc.sex.tsv"
params.PcLabel      = 'bbj_mainland'
params.NPcs         = 10
params.PhenoName    = 'pheno1'
// MUST match the prepared covariate header. rvtest_prepare_tools.py lowercases
// every column, so these are lowercase while the source file is uppercase.
// rvtest EXITS 0 when a --covar-name column is absent, silently fitting an
// uncovaried model, so RVTEST asserts the columns exist before it runs.
params.CovarName    = 'sex,' + (1..params.NPcs).collect { "pc${it}_avg" }.join(',')
params.CovarLabel   = "SEX + ${params.NPcs} ${params.PcLabel} PCs"

// ---- Reference resources ----------------------------------------------------
params.refFlatPath  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz'
params.snpeffDir    = '/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/output/snpEff.v6.index'
params.tommoVcf     = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz'
params.plink2       = '/home/b/b37974/plink2_alpha6/plink2'
params.tool_dir     = file(params.plink2).parent.toString()

// ---- Test matrix ------------------------------------------------------------
params.impactFilters = [
    [ tag: 'low_moderate_high', values: 'LOW MODERATE HIGH' ],
    [ tag: 'moderate_high',     values: 'MODERATE HIGH'     ],
    [ tag: 'high',              values: 'HIGH'              ],
]
// Zeggini dropped: measured against CMC it returned a bit-identical P for 63 % /
// 76 % / 96 % of genes (low_moderate_high / moderate_high / high), because the
// two differ only when a sample carries two rare variants in one gene, which at
// this MAF almost never happens. One kernel test and one burden test is the
// information the pair carried. See docs/METHODS.md.
params.testMethods = [
    [ tag: 'skato', opt: '--kernel skato' ],
    [ tag: 'cmc',   opt: '--burden cmc'   ],
]

// ---- Analysis options -------------------------------------------------------
params.performNorm     = false   // bcftools normalisation during PREPARE
params.allowRefAltSwap = false
params.MinNumVar       = 3       // minimum variants per gene to retain a test
params.FdrSignificant  = 0.05
params.FdrSuggestive   = 0.10
params.LabelGenes      = 12      // hit genes named on a scan figure
// Which scan the cross-cohort figure's per-gene panels show. Every stratum and
// method reaches the TABLE; a single figure cannot hold six scans x three
// cohorts, so it shows one and says which.
params.PeakStratum     = 'moderate_high'
params.PeakMethod      = 'skato'
params.LabelCrossCohort = 18     // genes on the cross-cohort panel, strongest first

// ---- Environment ------------------------------------------------------------
params.conda_env    = 'cteph_geno_pro'
params.rvtestBin    = '/home/b/b37974/rvtests/executable'

// ---- Provenance -------------------------------------------------------------
params.run_info_dir = "${params.out_dir}/_run_info"

// -----------------------------------------------------------------------------
//  HELPERS
// -----------------------------------------------------------------------------
def tpl(String template, String cohort) {
    return template.replace('@@COHORT@@', cohort)
}

// The tier a gene belongs to, derived from its id. gene_scan.py builds gene ids
// as `sig###_<GENE>` / `sug###_<GENE>`, and that prefix is the only tier
// information reachable INSIDE a process: a publishDir closure sees process
// inputs, and adding `tier` as an input would change every task hash. Asserts
// rather than guessing, so a change to the id convention fails loudly here.
def tierOf(String gene_id) {
    if (gene_id.startsWith('sig')) return 'significant'
    if (gene_id.startsWith('sug')) return 'suggestive'
    error "cannot determine the tier of '${gene_id}' — expected the sig/sug prefix gene_scan.py writes"
}

// Resolve the minAC for one cohort at DAG-build time, so the value that actually
// ran is knowable from the manifest rather than only from a log.
def minacFor(String cohort) {
    if (params.MinAC != null) return params.MinAC as Integer
    def f = file(tpl(params.MinACTpl, cohort))
    if (!f.exists()) {
        error "tuning.rv has not produced a minAC recommendation for '${cohort}'.\n" +
              "  expected: ${f}\n" +
              "  run tuning.rv first, or pass --MinAC <n> to override."
    }
    def row = f.readLines().find { it.startsWith('recommended_minac') }
    if (!row) error "no 'recommended_minac' key in ${f}"
    return row.split('\t')[1].trim() as Integer
}

// Which bin a chromosome was packed into. The per-gene fan-out knows a gene's
// CHROMOSOME (from the scan table) but every upstream artefact is keyed by BIN,
// so without this the join silently matches nothing for any chromosome that
// shares a bin.
def chromToBin() {
    def m = [:]
    params.ChromBins.each { b -> b.split(',').each { c -> m[c.trim()] = b.replace(',', '_') } }
    return m
}

// Figure family -> its directory under figures/. The catalogue publishes beside
// the figures it documents, so the two must agree.
params.FigureDirs = ['scan': '01.scan', 'gene': '02.gene']

// =============================================================================
//  PROCESSES
// =============================================================================

/* STEP 0 · SPLIT_CALLSET — one chromosome, minAC applied.
 *   Sample QC and depth-differential variant QC are ALREADY in tuning.rv's
 *   filtered set; re-applying them here would be both wrong and impossible (the
 *   old exclude lists no longer exist). The only filter left is the minAC that
 *   tuning.rv chose. */
process SPLIT_CALLSET {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}:bin${bin}"

    publishDir { "${params.out_dir}/${cohort}/00.callset" }, mode: 'copy', pattern: '*.log'

    input:
    tuple val(cohort), val(chroms), val(bin), val(minac), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), val(bin), path("bin${bin}.bed"), path("bin${bin}.bim"),
          path("bin${bin}.fam"), emit: plink
    path("bin${bin}.split.log")

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    ${params.plink2} --bfile ${bed.baseName} \\
        --chr ${chroms} --mac ${minac} \\
        --threads ${task.ext.threads ?: 1} \\
        --make-bed --out bin${bin} > /dev/null
    mv bin${bin}.log bin${bin}.split.log
    """
}

/* STEP 1 · PREPARE — PLINK -> bgzipped VCF, plus pheno/covar/refFlat in rvtest
 * format. The pheno/covar work repeats per chromosome; it costs seconds and
 * keeps the process self-contained. */
process PREPARE {
    executor 'slurm'
    queue    'gr10478b'
    time     '8h'
    tag      "${cohort}:bin${bin}"

    input:
    tuple val(cohort), val(bin), path(bed), path(bim), path(fam),
          path(pheno), path(covar), path(ref_flat)

    output:
    tuple val(cohort), val(bin), path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf
    tuple val(cohort), val(bin), path('*pheno_rvt.tsv'), path('*covar_rvt.tsv'), emit: pheno_cov
    tuple val(cohort), val(bin), path('*.nochr.txt.gz'),                        emit: refflat
    path('rvtest_prepare.log')

    script:
    def norm_flag = params.performNorm ? '--norm' : ''
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/rvtest_prepare_main.py \\
        --bed-prefix ${bed.baseName} \\
        --pheno-path ${pheno} \\
        --covar-path ${covar} \\
        --refflat-path ${ref_flat} \\
        --threads ${task.ext.threads ?: 1} \\
        --verbose --log-file rvtest_prepare.log ${norm_flag}
    """
}

/* STEP 2 · ANNOTATE — snpEff, per chromosome. */
process ANNOTATE {
    executor 'slurm'
    queue    'gr10478b'
    time     '8h'
    tag      "${cohort}:bin${bin}"

    input:
    tuple val(cohort), val(bin), path(vcf), path(tbi)

    output:
    tuple val(cohort), val(bin), path('*.snpeff.vcf.gz'), path('*.snpeff.vcf.gz.tbi'),
          emit: annotated
    tuple val(cohort), path('*.tsv'), emit: stats
    path('*.log')

    script:
    def swap_flag = params.allowRefAltSwap ? '--allow-ref-alt-swap' : ''
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/snpeff_anno_main.py \\
        --vcf-path ${vcf} \\
        --snpeff-dir ${params.snpeffDir} \\
        --parallel --max-workers ${task.ext.threads ?: 1} \\
        --threads ${task.ext.threads ?: 1} \\
        --keep-cache ${swap_flag}
    """
}

/* STEP 3 · IMPACT_FILTER — one VCF per (chromosome x impact stratum). */
process IMPACT_FILTER {
    executor 'slurm'
    queue    'gr10478b'
    time     '4h'
    tag      "${cohort}:bin${bin}:${stratum}"

    input:
    tuple val(cohort), val(bin), val(stratum), val(values), path(vcf), path(tbi)

    output:
    tuple val(cohort), val(bin), val(stratum), path("*.${stratum}.vcf.gz"),
          path("*.${stratum}.vcf.gz.tbi"), emit: filtered
    path('*.summary.json'), optional: true

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/info_filter_main.py \\
        --input ${vcf} \\
        --info-key impact --values ${values} \\
        --out-prefix bin${bin}.${stratum} \\
        --threads ${task.ext.threads ?: 1} \\
        --check-chr-prefix --no-keep-chr-prefix --verbose
    """
}

/* STEP 4 · RVTEST — the association, per (chromosome x stratum x method).
 *   The --covar-name assertion is not optional: rvtest exits 0 when a named
 *   covariate column is missing and silently fits an uncovaried model, so an
 *   unchecked run can look perfect and be wrong. */
process RVTEST {
    executor 'slurm'
    queue    'gr10478b'
    time     '12h'
    tag      "${cohort}:bin${bin}:${stratum}:${method}"

    input:
    tuple val(cohort), val(bin), val(stratum), path(vcf), path(tbi),
          val(method), val(method_opt), path(pheno), path(covar), path(refflat)

    output:
    tuple val(cohort), val(stratum), val(method), path('*.assoc'), emit: assoc
    path('*.log')

    script:
    """
    export PATH=${params.rvtestBin}:\$PATH
    head -1 ${covar} | tr '\\t' '\\n' > .covar_header
    for c in \$(echo '${params.CovarName}' | tr ',' ' '); do
        grep -qx "\$c" .covar_header || {
            echo "ERROR: covariate column '\$c' is absent from ${covar}." >&2
            echo "rvtest would exit 0 and fit an uncovaried model." >&2
            exit 1
        }
    done
    rvtest \\
        --inVcf ${vcf} \\
        --pheno ${pheno} --pheno-name ${params.PhenoName} \\
        --covar ${covar} --covar-name ${params.CovarName} \\
        --geneFile ${refflat} \\
        --out bin${bin}.${stratum}.${method} \\
        --noweb --numThread ${task.ext.threads ?: 1} \\
        ${method_opt}
    """
}

/* STEP 5 · MERGE_ASSOC — the barrier the chromosome split requires.
 *   Benjamini-Hochberg is a property of the WHOLE scan, so the per-chromosome
 *   tables must be one table before any tier is assigned. Concatenating with a
 *   single header is the only thing this step does. */
process MERGE_ASSOC {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${stratum}:${method}"

    publishDir { "${params.out_dir}/${cohort}/04.assoc/${stratum}/${method}" }, mode: 'copy'

    input:
    tuple val(cohort), val(stratum), val(method), path(parts, stageAs: 'part_*.assoc')

    output:
    tuple val(cohort), val(stratum), val(method), path('merged.assoc'), emit: assoc

    script:
    """
    set -euo pipefail
    first=\$(ls part_*.assoc | head -1)
    head -1 "\$first" > merged.assoc
    for f in part_*.assoc; do tail -n +2 "\$f" >> merged.assoc; done
    echo "[MERGE_ASSOC] ${cohort} ${stratum}/${method}: \$(ls part_*.assoc | wc -l) chromosome parts, \$((\$(wc -l < merged.assoc) - 1)) gene rows"
    """
}

/* STEP 6 · GENE_SCAN — NumVar filter, genome-wide BH, tier assignment. */
process GENE_SCAN {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${stratum}:${method}"

    publishDir { "${params.out_dir}/${cohort}/04.assoc/${stratum}/${method}" }, mode: 'copy'

    input:
    tuple val(cohort), val(stratum), val(method), path(assoc)

    output:
    tuple val(cohort), val(stratum), val(method), path('gene_scan.tsv'),
          path('scan_qc.tsv'), emit: scan
    tuple val(cohort), val(stratum), val(method), path('gene_hits.tsv'), emit: hits

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/rvtest_post_process.py \\
        --input ${assoc} --output filtered.assoc \\
        --num-var-threshold ${params.MinNumVar}
    python3 ${params.script_dir}/gene_scan.py \\
        --assoc filtered.assoc \\
        --cohort ${cohort} --stratum ${stratum} --method ${method} \\
        --fdr-significant ${params.FdrSignificant} \\
        --fdr-suggestive ${params.FdrSuggestive} \\
        --out-dir .
    """
}

/* STEP 7 · PLOT_SCAN — gene-level Manhattan + QQ. */
process PLOT_SCAN {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${stratum}:${method}"

    publishDir { "${params.out_dir}/${cohort}/figures/01.scan" }, mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(stratum), val(method), path(scan), path(qc)

    output:
    path("scan.${stratum}.${method}.png")
    tuple val(cohort), path("scan.${stratum}.${method}.stats.json"), emit: stats

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_gene_scan.py \\
        --scan ${scan} --qc ${qc} \\
        --cohort ${cohort} --stratum ${stratum} --method ${method} \\
        --covar-label '${params.CovarLabel}' \\
        --label-genes ${params.LabelGenes} \\
        --out-png scan.${stratum}.${method}.png
    """
}

/* STEP 8 · GENE_DETAIL — the per-gene fan-out over every tiered gene.
 *   Was `check_gene/`, a hand-driven shell wrapper with a hardcoded cohort. It
 *   is the same question the pipeline should answer for every hit, so it runs
 *   here instead of being remembered. */
process GENE_DETAIL {
    executor 'slurm'
    queue    'gr10478b'
    time     '2h'
    tag      "${cohort}:${stratum}/${method}:${gene_id}"

    publishDir { "${params.out_dir}/${cohort}/05.genes/${tierOf(gene_id)}/${stratum}/${method}/${gene_id}" },
               mode: 'copy', pattern: '*.{summary.txt,sample_details.tsv}'

    input:
    tuple val(cohort), val(gene_id), val(gene), val(stratum), val(method),
          path(assoc), path(vcf), path(tbi), path(bed), path(bim), path(fam),
          path(pheno), path(covar), path(refflat)

    output:
    tuple val(cohort), val(stratum), val(method), val(gene_id), val(gene),
          path('*.summary.txt'), emit: summary
    path('*.sample_details.tsv'), optional: true

    script:
    """
    source activate ${params.conda_env}
    export PATH=${params.rvtestBin}:\$PATH
    python3 ${params.script_dir}/gene_detail.py \\
        --gene ${gene} \\
        --assoc-file ${assoc} \\
        --vcf-file ${vcf} \\
        --plink-prefix ${bed.baseName} \\
        --tommo-vcf ${params.tommoVcf} \\
        --pheno-file ${pheno} --covar-file ${covar} \\
        --covar-name ${params.CovarName} \\
        --refflat-file ${refflat} \\
        --plink2-path ${params.plink2} \\
        --group-name ${stratum}.${method} \\
        --out-dir . --out-log ${gene}.${stratum}.${method}.summary.txt
    """
}


/* STEP 8b · PLOT_GENE_DETAIL — one figure per tiered gene. */
process PLOT_GENE_DETAIL {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${stratum}/${method}:${gene_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/02.gene/${tierOf(gene_id)}" },
               mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(stratum), val(method), val(gene_id), val(gene), path(summary)

    output:
    path("gene.${stratum}.${method}.${gene_id}.png")
    tuple val(cohort), path("gene.${stratum}.${method}.${gene_id}.stats.json"), emit: stats

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_gene_detail.py \
        --summary ${summary} --gene ${gene} --gene-id ${gene_id} \
        --cohort ${cohort} --stratum ${stratum} --method ${method} \
        --out-png gene.${stratum}.${method}.${gene_id}.png
    """
}

/* STEP 8c · CATALOGUE_FIGURES — ONE document per figure family per cohort.
 *   A fan-out over N genes produces N figures; a sidecar each would be N
 *   documents whose prose is identical and whose numbers are a dozen lines
 *   apart. The prose belongs to the family, so it is written once and the
 *   per-gene numbers become rows of one table. */
process CATALOGUE_FIGURES {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      "${cohort}:${family}"

    publishDir { "${params.out_dir}/${cohort}/figures/${params.FigureDirs[family]}" },
               mode: 'copy'

    input:
    tuple val(cohort), val(family), path(stats)

    output:
    path('README.md')

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/figure_catalogue.py \
        --family ${family} --cohort ${cohort} \
        --stats-dir . --figure-dir ${params.FigureDirs[family]} \
        --out README.md
    """
}

/* STEP 10 · CROSS_COHORT_GENES — the union of tiered genes, in every cohort. */
process CROSS_COHORT_GENES {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'all'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    input:
    tuple val(cohorts), val(strata), val(methods),
          path(scans, stageAs: 'scan_*.tsv'), path(qcs, stageAs: 'qc_*.tsv')

    output:
    tuple path('gene_crosscohort.tsv'), path('scan_qc_all.tsv'), emit: tables
    path('gene_hits_all.tsv')

    script:
    def s = [cohorts, strata, methods, scans].transpose()
        .collect { c, st, m, f -> "--scan '${c}=${st}=${m}=${f}'" }.join(' ')
    def q = [cohorts, strata, methods, qcs].transpose()
        .collect { c, st, m, f -> "--qc '${c}=${st}=${m}=${f}'" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/cross_cohort_genes.py \
        ${s} ${q} \
        --cohort-order ${params.Cohorts.join(',')} \
        --out-genes gene_crosscohort.tsv \
        --out-scans scan_qc_all.tsv \
        --out-hits gene_hits_all.tsv
    """
}

/* STEP 11 · COMPARE_COHORTS — the three cohorts side by side. */
process COMPARE_COHORTS {
    executor 'slurm'
    queue    'gr10478b'
    time     '1h'
    tag      'all'

    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy'

    input:
    tuple path(genes), path(scans)

    output:
    path('cohort_compare_genes.*')

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_cohort_compare_genes.py \
        --genes ${genes} --scans ${scans} \
        --cohort-order ${params.Cohorts.join(',')} \
        --peak-stratum ${params.PeakStratum} --peak-method ${params.PeakMethod} \
        --label-genes ${params.LabelCrossCohort} \
        --out-png cohort_compare_genes.png
    """
}

/* STEP 9 · WRITE_RUN_MANIFEST — what actually ran, including the inherited
 * minAC, so the run is reconstructable without reading tuning.rv's logs. */
process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'
    publishDir "${params.run_info_dir}", mode: 'copy'

    input:
    val manifest_json

    output:
    path('run_manifest.json')

    script:
    """
    cat > run_manifest.json <<'MANIFEST'
${manifest_json}
MANIFEST
    """
}

// =============================================================================
//  WORKFLOW
// =============================================================================
workflow {

    // -- validate the cohort list ---------------------------------------------
    def cohortList = (params.Cohorts instanceof List) ? params.Cohorts
                     : params.Cohorts.toString().split(',').collect { it.trim() }.findAll { it }
    cohortList.each { c ->
        if (!params.ValidCohorts.contains(c)) {
            error "Invalid cohort '${c}'. Must be one of: ${params.ValidCohorts.join(', ')}"
        }
    }

    // -- resolve the tuning.rv contract, per cohort ---------------------------
    def minac = cohortList.collectEntries { [(it): minacFor(it)] }
    log.info "[assoc_rvtest] cohorts: ${cohortList.join(', ')}"
    minac.each { c, v ->
        log.info "[assoc_rvtest]   ${c}: minAC = ${v}" +
                 (params.MinAC != null ? '  (OVERRIDDEN on the command line)' : '  (from tuning.rv)')
    }

    ch_cohort = channel.fromList(cohortList)

    ch_callset = ch_cohort.map { c ->
        def p = tpl(params.CallsetTpl, c)
        tuple(c, minac[c],
              file("${p}.bed", checkIfExists: true),
              file("${p}.bim", checkIfExists: true),
              file("${p}.fam", checkIfExists: true))
    }

    ch_pheno_covar = ch_cohort.map { c ->
        tuple(c,
              file(tpl(params.PhenoTpl, c), checkIfExists: true),
              file(tpl(params.CovarTpl, c), checkIfExists: true))
    }
    ch_refflat = channel.value(file(params.refFlatPath, checkIfExists: true))

    // -- [0] one chromosome at a time -----------------------------------------
    SPLIT_CALLSET(
        ch_callset
            .combine(channel.fromList(params.ChromBins))
            .map { c, ac, bed, bim, fam, chroms ->
                   tuple(c, chroms, chroms.replace(',', '_'), ac, bed, bim, fam) })

    // -- [1] prepare -----------------------------------------------------------
    PREPARE(SPLIT_CALLSET.out.plink
        .combine(ch_pheno_covar, by: 0)
        .combine(ch_refflat)
        .map { c, bin, bed, bim, fam, ph, cv, rf -> tuple(c, bin, bed, bim, fam, ph, cv, rf) })

    // pheno/covar and refFlat are bin-independent: every bin writes an identical
    // copy. Take them from ONE NAMED BIN, never from a groupTuple. groupTuple's
    // arrival order is not deterministic, so picking element [0] chose a different
    // bin's copy on each run; the copies live at different work paths, so that
    // changed RVTEST's input hash and re-ran all 216 association tasks on every
    // -resume. Naming the bin makes the choice stable by construction.
    def CANON_BIN = params.ChromBins[0].replace(',', '_')
    ch_prep_pc = PREPARE.out.pheno_cov
        .filter { c, bin, ph, cv -> bin == CANON_BIN }
        .map    { c, bin, ph, cv -> tuple(c, ph, cv) }
    ch_prep_rf = PREPARE.out.refflat
        .filter { c, bin, rf -> bin == CANON_BIN }
        .map    { c, bin, rf -> tuple(c, rf) }

    // -- [2] annotate ----------------------------------------------------------
    ANNOTATE(PREPARE.out.vcf)

    // -- [3] impact strata -----------------------------------------------------
    ch_strata = channel.fromList(params.impactFilters.collect { [it.tag, it.values] })
    IMPACT_FILTER(ANNOTATE.out.annotated
        .combine(ch_strata)
        .map { c, bin, vcf, tbi, tag, vals -> tuple(c, bin, tag, vals, vcf, tbi) })

    // -- [4] the association ---------------------------------------------------
    ch_methods = channel.fromList(params.testMethods.collect { [it.tag, it.opt] })
    RVTEST(IMPACT_FILTER.out.filtered
        .combine(ch_methods)
        .combine(ch_prep_pc, by: 0)
        .combine(ch_prep_rf, by: 0)
        .map { c, bin, stratum, vcf, tbi, m, opt, ph, cv, rf ->
               tuple(c, bin, stratum, vcf, tbi, m, opt, ph, cv, rf) })

    // -- [5] the barrier: one table per scan before any FDR --------------------
    // groupTuple(by: [0, 1, 2]) with three PLAIN STRING keys. A tuple used as a
    // key matches nothing and the channel silently never emits.
    MERGE_ASSOC(RVTEST.out.assoc.groupTuple(by: [0, 1, 2]))

    // -- [6..7] tiers and the scan figure --------------------------------------
    GENE_SCAN(MERGE_ASSOC.out.assoc)
    PLOT_SCAN(GENE_SCAN.out.scan)

    // -- [8] per-gene fan-out --------------------------------------------------
    // Each hit gene needs FOUR things joined to it, and each join is keyed on a
    // different slice of the identity, so the tuple is reshaped between them:
    //   (cohort, stratum, method) -> the merged assoc table it was called from
    //   (cohort, bin, stratum)    -> the impact-filtered VCF holding its variants
    //   (cohort, bin)             -> that bin's PLINK set
    //   (cohort)                  -> pheno / covar / refFlat
    // groupTuple/combine keys are PLAIN STRINGS throughout; a tuple used as a key
    // matches nothing and the channel silently never emits.
    def BIN_OF = chromToBin()
    ch_gene = GENE_SCAN.out.hits
        .splitCsv(header: true, sep: '\t', elem: 3)
        .map { c, stratum, method, r ->
               tuple(c, stratum, method, BIN_OF[r.chrom as String], r.gene_id, r.Gene) }

    GENE_DETAIL(
        ch_gene
            .combine(MERGE_ASSOC.out.assoc, by: [0, 1, 2])
            .map { c, s, m, bin, gid, gene, assoc ->
                   tuple(c, bin, s, m, gid, gene, assoc) }
            .combine(IMPACT_FILTER.out.filtered, by: [0, 1, 2])
            .map { c, bin, s, m, gid, gene, assoc, vcf, tbi ->
                   tuple(c, bin, s, m, gid, gene, assoc, vcf, tbi) }
            .combine(SPLIT_CALLSET.out.plink, by: [0, 1])
            .combine(ch_prep_pc, by: 0)
            .combine(ch_prep_rf, by: 0)
            .map { c, bin, s, m, gid, gene, assoc, vcf, tbi, bed, bim, fam, ph, cv, rf ->
                   tuple(c, gid, gene, s, m, assoc, vcf, tbi, bed, bim, fam, ph, cv, rf) })

    // -- [8b] one figure per tiered gene ---------------------------------------
    PLOT_GENE_DETAIL(GENE_DETAIL.out.summary)

    // -- [8c] one document per figure family per cohort ------------------------
    // groupTuple(by: [0, 1]) with two PLAIN STRING keys. A tuple used as the key
    // matches nothing and the channel silently never emits.
    CATALOGUE_FIGURES(
        PLOT_GENE_DETAIL.out.stats.map { c, s -> tuple(c, 'gene', s) }
            .mix(PLOT_SCAN.out.stats.map { c, s -> tuple(c, 'scan', s) })
            .groupTuple(by: [0, 1]))

    // -- [10..11] across cohorts ------------------------------------------------
    // A barrier, deliberately: the gene list is the UNION over cohorts, so no
    // comparison can start until every cohort has finished its scans.
    ch_all_scans = GENE_SCAN.out.scan
        .toList()
        .map { rows -> tuple(rows.collect { it[0] }, rows.collect { it[1] },
                             rows.collect { it[2] }, rows.collect { it[3] },
                             rows.collect { it[4] }) }
        .filter { it[0] }
    CROSS_COHORT_GENES(ch_all_scans)
    COMPARE_COHORTS(CROSS_COHORT_GENES.out.tables)

    // -- [9] manifest ----------------------------------------------------------
    def manifest = """{
  "component": "assoc_rvtest",
  "cohorts": ${groovy.json.JsonOutput.toJson(cohortList)},
  "minac": ${groovy.json.JsonOutput.toJson(minac)},
  "minac_source": "${params.MinAC != null ? 'command line override' : 'tuning.rv recommendation'}",
  "chromosome_bins": ${groovy.json.JsonOutput.toJson(params.ChromBins)},
  "impact_strata": ${groovy.json.JsonOutput.toJson(params.impactFilters*.tag)},
  "test_methods": ${groovy.json.JsonOutput.toJson(params.testMethods*.tag)},
  "fdr_significant": ${params.FdrSignificant},
  "fdr_suggestive": ${params.FdrSuggestive},
  "min_num_var": ${params.MinNumVar},
  "pheno_name": "${params.PhenoName}",
  "covar_name": "${params.CovarName}",
  "covar_label": "${params.CovarLabel}"
}"""
    WRITE_RUN_MANIFEST(channel.value(manifest))
}
