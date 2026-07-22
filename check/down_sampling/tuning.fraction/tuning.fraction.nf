nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Per-platform down-sampling fractions
//
//   down_sampling.nf keeps a fraction of each sample's reads so that it matches
//   the depth of the cohort's reference platform. This pipeline measures what
//   that fraction has to be — one per SEQUENCING PLATFORM, not one per
//   Target_Depth group, because a group is not one population: its platforms sit
//   at very different depths in the regions we analyse, so a single group-wide
//   fraction would over-down-sample some and under-down-sample others.
//
//       fraction(P) = median CRAM depth of the baseline platform
//                     ------------------------------------------
//                     median CRAM depth of platform P
//
//   Every depth is measured on the CRAMs with `samtools coverage`, restricted to
//   params.regions — that is exactly the read pool `samtools view -s` subsamples.
//   The sample sheet's Observed_Depth is deliberately NOT used: it is genome-wide
//   and computed upstream from the FASTQs, so it does not describe what the
//   down-sampling acts on.
//
//   Samples of refined_core without a CRAM cannot be measured and are excluded;
//   the count is reported per platform so the exclusion is never silent.
//
//   Stages:
//     00  CRAM_REGION_DP      per sample chunk: aligned depth in params.regions,
//                             with/without duplicates, plus read length
//     01  PLATFORM_FRACTIONS  one fraction per platform -> tsv + log + figure
// ─────────────────────────────────────────────────────────────────────────────

// ── Inputs ───────────────────────────────────────────────────────────────────
// Sample sheet: id, Target_Depth group, platform, CRAM path. Used for metadata
// and CRAM paths only — no depth is taken from it.
params.cram_info           = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tmp/cram.v6/cram.v6.summary.csv'
// Analysis cohort whose samples define the comparison set.
params.refined_core_prefix = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/15_random_model_prep/refined_core/cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.random_model'
params.fasta               = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

// ── Baseline ─────────────────────────────────────────────────────────────────
// The platform every other platform is brought down to. HiSeqX 15x is the
// cohort's reference: it carries the overwhelming majority of the samples, so
// matching it is what keeps the cohort homogeneous.
params.baseline_platform   = 'HiSeqX 15x'

// ── Outputs / runtime ────────────────────────────────────────────────────────
params.results_dir         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/tuning.fraction/results'
params.script_dir          = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/tuning.fraction/scripts'
params.conda_env_activate  = 'cteph_geno_pro'
// Samples per CRAM_REGION_DP task, and parallel samples inside one task.
params.cram_chunk_size     = 200
params.cram_threads        = 12

// ── Sample-sheet columns ─────────────────────────────────────────────────────
params.sample_id_col       = 'ID_JHRPv6'
params.target_dp_col       = 'Target_Depth'
params.platform_col        = 'WGS_Platform'
params.cram_path_col       = 'Cram_Path'
params.cram_found_col      = 'Cram_Found'

// ── Target regions (keep in sync with down_sampling.nf) ──────────────────────
params.regions             = [
    'chr16:53703963-54121941',
    'chr4:76306733-76311130',
]
params.regions_label       = 'target regions'

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// The sample sheet has been regenerated as TAB- and as COMMA-separated at
// different times (the .csv name says nothing). Guessing wrong does not throw:
// every field lands in column 0, every lookup returns null, the filters match
// nothing and the run "succeeds" with no tasks. So sniff the header instead of
// hard-coding, and fail loudly if the expected columns are gone.
def sniffSep(path) {
    def header = file(path).withReader { r -> r.readLine() }
    if (header == null) {
        error "Empty sample sheet: ${path}"
    }
    return header.count('\t') >= header.count(',') && header.contains('\t') ? '\t' : ','
}

def readTable(path, required) {
    def sep = sniffSep(path)
    def lines = file(path).readLines().findAll { l -> l.trim() }
    def hdr = lines[0].split(sep == '\t' ? /\t/ : /,/, -1)*.trim()
    def missing = required.findAll { c -> !hdr.contains(c) }
    if (missing) {
        error "Sample sheet ${path} (separator '${sep == '\t' ? '\\t' : ','}') is missing " +
              "column(s) ${missing}. Found: ${hdr}"
    }
    log.info "[tuning.fraction] sample sheet: separator '${sep == '\t' ? '\\t' : ','}', " +
             "${hdr.size()} columns, ${lines.size() - 1} rows"
    return lines.drop(1).collect { l ->
        def f = l.split(sep == '\t' ? /\t/ : /,/, -1)
        [hdr, f].transpose().collectEntries { k, v -> [(k): v?.trim()] }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: aligned read depth of the target regions, straight from the CRAMs,
//   for one chunk of samples. Reports depth with and without duplicates, and the
//   read length, so both are measured rather than assumed.
// ─────────────────────────────────────────────────────────────────────────────
process CRAM_REGION_DP {

    tag "chunk_${chunk.baseName}"

    executor 'slurm'
    queue    'gr10478b'
    time     '6h'

    publishDir "${params.results_dir}/00_cram_region_dp", mode: 'symlink'

    input:
    path chunk
    path script

    output:
    path "${chunk.baseName}.cramdp.tsv", emit: dp

    script:
    def region_list = (params.regions instanceof List ? params.regions : "${params.regions}".split(','))
    def regions_csv = region_list.collect { r -> r.toString().trim() }.findAll { r -> r }.join(',')
    """
    export PATH=/LARGE0/gr10478/b37974/software/samtools/bin:/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${script} \\
        ${chunk} \\
        '${regions_csv}' \\
        ${params.fasta} \\
        ${params.cram_threads} \\
        ${chunk.baseName}.cramdp.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: one down-sampling fraction per platform, relative to the baseline.
//   Emits a machine-readable TSV for down_sampling.nf, plus a log and a figure.
// ─────────────────────────────────────────────────────────────────────────────
process PLATFORM_FRACTIONS {

    tag "${params.baseline_platform}"

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/01_platform_fractions", mode: 'symlink'

    input:
    path cram_tsvs
    path script
    // Taken as path inputs, not referenced through params: the sample sheet is
    // regenerated in place, so its path is stable while its content is not. Only a
    // path input is content-hashed — via params, -resume would serve a stale result.
    path cram_info
    path refined_fam

    output:
    path "platform_fractions.tsv",  emit: fractions
    path "platform_fractions.log",  emit: log
    path "platform_fractions.png",  emit: figure

    script:
    def region_list = (params.regions instanceof List ? params.regions : "${params.regions}".split(','))
    def regions_csv = region_list.collect { r -> r.toString().trim() }.findAll { r -> r }.join(',')
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    python3 ${script} \\
        --cram-dp ${cram_tsvs} \\
        --regions '${regions_csv}' \\
        --cram-info ${cram_info} \\
        --refined-fam ${refined_fam} \\
        --baseline-platform "${params.baseline_platform}" \\
        --sample-id-col "${params.sample_id_col}" \\
        --target-dp-col "${params.target_dp_col}" \\
        --platform-col "${params.platform_col}" \\
        --regions-label "${params.regions_label}" \\
        --out-log platform_fractions.log \\
        --out-tsv platform_fractions.tsv \\
        --out-fig platform_fractions.png
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // Scripts are passed as inputs, not just referenced by path, so that editing
    // one invalidates the task hash and -resume actually re-runs it.
    cramdepth_script = file("${params.script_dir}/cram_region_depth.sh", checkIfExists: true)
    fractions_script = file("${params.script_dir}/platform_fractions.py", checkIfExists: true)

    // refined_core sample set (IID = 2nd fam column). The fam is small, so read it
    // directly. Split on any whitespace: plink1 writes space-, plink2 tab-delimited.
    def refined_iids = file("${params.refined_core_prefix}.fam", checkIfExists: true)
        .readLines()
        .findAll { line -> line.trim() }
        .collect { line -> line.trim().split(/\s+/)[1] } as Set
    log.info "[tuning.fraction] refined_core samples: ${refined_iids.size()}"

    // Only samples with a CRAM can be measured. Report the shortfall per platform
    // here as well, so a missing-CRAM platform is visible before the run finishes.
    def cram_rows = readTable(params.cram_info,
        [params.sample_id_col, params.cram_found_col, params.cram_path_col,
         params.target_dp_col, params.platform_col])
    def in_core = cram_rows.findAll { row -> refined_iids.contains(row[params.sample_id_col]) }
    def with_cram = in_core.findAll { row -> row[params.cram_found_col] == 'True' }
    if (!with_cram) {
        error "No refined_core sample has a CRAM in ${params.cram_info}. Every depth in " +
              "this pipeline is measured from CRAM, so there would be nothing to measure."
    }
    def missing = in_core.size() - with_cram.size()
    log.info "[tuning.fraction] refined_core with CRAM: ${with_cram.size()}" +
             (missing ? " (excluded, no CRAM: ${missing})" : "")
    with_cram.groupBy { row -> row[params.platform_col] }.each { pf, rs ->
        def miss = in_core.count { row -> row[params.platform_col] == pf } - rs.size()
        log.info "[tuning.fraction]   ${pf}: ${rs.size()}" + (miss ? " (+${miss} without CRAM)" : "")
    }
    if (!with_cram.any { row -> row[params.platform_col] == params.baseline_platform }) {
        error "Baseline platform '${params.baseline_platform}' has no refined_core sample " +
              "with a CRAM — every fraction is relative to it, so it must be measurable."
    }

    // ── 00: aligned depth from the CRAMs, in chunks ───────────────────────────
    cram_chunks_ch = channel
        .fromList(with_cram.collect { row ->
            "${row[params.sample_id_col]}\t${row[params.cram_path_col]}\n" })
        .collectFile(name: 'cram_list.txt', sort: true)
        .splitText(by: params.cram_chunk_size, file: true)

    cram_dp_ch = CRAM_REGION_DP(cram_chunks_ch, cramdepth_script).dp

    // ── 01: one fraction per platform ─────────────────────────────────────────
    PLATFORM_FRACTIONS(cram_dp_ch.collect(), fractions_script,
                       file(params.cram_info, checkIfExists: true),
                       file("${params.refined_core_prefix}.fam", checkIfExists: true))
}
