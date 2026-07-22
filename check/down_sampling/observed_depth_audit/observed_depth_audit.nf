nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Observed_Depth audit — the sheet's depth assumed a 150 bp read length
//
//   The sheet labels DNBSeq-T7 samples ~30x, but their CRAMs carry ~18x.
//   Observed_Depth is computed upstream from the FASTQs and copied verbatim into
//   the sheet, so it never saw the CRAM. This pipeline measures the CRAMs
//   independently and tests one explanation:
//
//       Observed_Depth = reads x 150 / G, with the read length hard-coded to 150
//       for every sample instead of each sample's real length. T7 reads are
//       100 bp, so its depth is overstated by 150/100 = 1.5x.
//
//   The test writes both depths as a read density rho times a read length:
//
//       CRAM depth = rho x (real read length)
//       Observed   = rho x (whatever length was assumed)
//
//   Dividing cancels rho, and with it the genome size:
//
//       implied_readlen = readlen x Observed / CRAM depth
//           hard-coded 150   -> implied == 150 for every platform
//           real read length -> implied == measured
//
//   That cancellation is what makes this cheap: the audit needs a read DENSITY,
//   not a genome-wide read total, so a few Mb of probe regions settles it in
//   seconds per sample. A whole-CRAM `samtools stats` pass costs ~1.5 min each.
//
//   Only T7 can discriminate — every other platform is a ~150 bp library, so for
//   them both hypotheses predict the same thing. They are the calibration
//   control, and the log says so rather than presenting them as confirmation.
//
//   Stages:
//     00  CRAM_REGION_STATS    per sample: depth (dup-kept and dup-free) and read
//                              length over the probe regions
//     01  OBSERVED_DEPTH_AUDIT the implied-read-length test -> tsv + log + figure
//         WRITE_REPORT         the whole story as README.md, published alongside
// ─────────────────────────────────────────────────────────────────────────────

// ── Inputs ───────────────────────────────────────────────────────────────────
// The sheet supplies Observed_Depth (the number under audit), the platform label
// and CRAM paths. Every other quantity is measured from the CRAM.
params.cram_info          = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/tmp/cram.v6/cram.v6.summary.csv'
params.fasta              = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'

// ── Probe regions ────────────────────────────────────────────────────────────
// Ten 1 Mb windows on ten different chromosomes, standing in for the genome-wide
// read density that Observed_Depth claims to describe. Mid-arm and away from the
// centromeres, so no window sits in a region of unusual mappability; spreading
// them over ten chromosomes keeps any single locus from driving the answer.
// The effect under test is a 1.5x platform-wide shift, not a local one.
params.regions            = [
    'chr1:100000000-101000000',
    'chr2:100000000-101000000',
    'chr3:100000000-101000000',
    'chr4:100000000-101000000',
    'chr5:100000000-101000000',
    'chr6:100000000-101000000',
    'chr7:100000000-101000000',
    'chr8:100000000-101000000',
    'chr10:100000000-101000000',
    'chr11:100000000-101000000',
]

// ── Sampling ─────────────────────────────────────────────────────────────────
// Per platform; platforms with fewer samples contribute all of them. Five is
// enough: the effect under test is a 1.5x per-platform systematic, read length has
// no within-platform spread at all (fixed-length libraries), and depth varies only
// a few percent. There is no reason to buy precision the claim does not need.
params.samples_per_platform = 5
params.sample_seed          = 42

// ── Outputs / runtime ────────────────────────────────────────────────────────
params.results_dir        = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/observed_depth_audit/results'
params.script_dir         = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/down_sampling/observed_depth_audit/scripts'
params.conda_env_activate = 'cteph_geno_pro'
params.stats_chunk_size   = 5
params.stats_threads      = 5

// ── Sheet columns ────────────────────────────────────────────────────────────
params.sample_id_col      = 'ID_JHRPv6'
params.observed_dp_col    = 'Observed_Depth'
params.platform_col       = 'WGS_Platform'
params.cram_path_col      = 'Cram_Path'
params.cram_found_col     = 'Cram_Found'

// ─────────────────────────────────────────────────────────────────────────────
// Helpers — the sheet has shipped both TAB- and COMMA-separated (the .csv name
// says nothing), and guessing wrong does not throw: every field lands in column
// 0, lookups return null, filters match nothing, and the run "succeeds" with no
// tasks. So sniff the header and fail loudly if the columns are gone.
// ─────────────────────────────────────────────────────────────────────────────
def sniffSep(path) {
    def header = file(path).withReader { r -> r.readLine() }
    if (header == null) error "Empty sample sheet: ${path}"
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
    log.info "[audit] sample sheet: separator '${sep == '\t' ? '\\t' : ','}', " +
             "${hdr.size()} columns, ${lines.size() - 1} rows"
    return lines.drop(1).collect { l ->
        def f = l.split(sep == '\t' ? /\t/ : /,/, -1)
        [hdr, f].transpose().collectEntries { k, v -> [(k): v?.trim()] }
    }
}

def regionsCsv() {
    def rl = (params.regions instanceof List ? params.regions : "${params.regions}".split(','))
    return rl.collect { r -> r.toString().trim() }.findAll { r -> r }.join(',')
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: depth and read length over the probe regions.
// ─────────────────────────────────────────────────────────────────────────────
process CRAM_REGION_STATS {

    tag "chunk_${chunk.baseName}"

    executor 'slurm'
    queue    'gr10478b'
    time     '4h'

    publishDir "${params.results_dir}/00_cram_region_stats", mode: 'symlink'

    input:
    path chunk
    path script

    output:
    path "${chunk.baseName}.regstats.tsv", emit: stats

    script:
    """
    export PATH=/LARGE0/gr10478/b37974/software/samtools/bin:/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    bash ${script} \\
        ${chunk} \\
        '${regionsCsv()}' \\
        ${params.fasta} \\
        ${params.stats_threads} \\
        ${chunk.baseName}.regstats.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the implied-read-length test, plus log and figure.
// ─────────────────────────────────────────────────────────────────────────────
process OBSERVED_DEPTH_AUDIT {

    executor 'slurm'
    queue    'gr10478b'
    time     '2h'

    publishDir "${params.results_dir}/01_audit", mode: 'symlink'

    input:
    path region_tsvs
    path script
    // A path input, not params: the sheet is regenerated in place, so its path is
    // stable while its content is not. Only a path input is content-hashed —
    // through params, -resume would happily serve a stale audit.
    path cram_info

    output:
    path "observed_depth_audit.tsv",         emit: tsv
    path "observed_depth_audit.samples.tsv", emit: samples
    path "observed_depth_audit.log",         emit: log
    path "observed_depth_audit.png",         emit: figure

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    python3 ${script} \\
        --region-stats ${region_tsvs} \\
        --cram-info ${cram_info} \\
        --regions '${regionsCsv()}' \\
        --sample-id-col "${params.sample_id_col}" \\
        --observed-dp-col "${params.observed_dp_col}" \\
        --platform-col "${params.platform_col}" \\
        --out-log observed_depth_audit.log \\
        --out-tsv observed_depth_audit.tsv \\
        --out-samples-tsv observed_depth_audit.samples.tsv \\
        --out-fig observed_depth_audit.png
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Process: the findings as a standalone markdown document.
//   Generated, not hand-written, so the prose cannot drift away from the numbers:
//   every value in it is read back out of the audit's own tsv.
// ─────────────────────────────────────────────────────────────────────────────
process WRITE_REPORT {

    executor 'slurm'
    queue    'gr10478b'
    time     '1h'

    // Published alongside the analysis it describes, so the README's link to the
    // figure is just a sibling filename and the directory needs no second copy.
    publishDir "${params.results_dir}/01_audit", mode: 'symlink'

    input:
    path tsv
    path samples_tsv
    path log_file
    // Taken as an input for its name and for the dependency — the report links to
    // this figure — but never copied: it is already published in this directory.
    path figure
    path script

    output:
    path "README.md", emit: report

    script:
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env_activate}
    set -euo pipefail

    python3 ${script} \\
        --tsv ${tsv} \\
        --samples-tsv ${samples_tsv} \\
        --log ${log_file} \\
        --figure-name ${figure} \\
        --regions '${regionsCsv()}' \\
        --samples-per-platform ${params.samples_per_platform} \\
        --sample-seed ${params.sample_seed} \\
        --out-md README.md
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    stats_script  = file("${params.script_dir}/cram_region_stats.sh", checkIfExists: true)
    audit_script  = file("${params.script_dir}/observed_depth_audit.py", checkIfExists: true)
    report_script = file("${params.script_dir}/write_report.py", checkIfExists: true)

    def rows = readTable(params.cram_info,
        [params.sample_id_col, params.observed_dp_col, params.platform_col,
         params.cram_path_col, params.cram_found_col])

    // Auditable = has a CRAM to measure and a sheet depth to measure it against.
    def usable = rows.findAll { r ->
        r[params.cram_found_col] == 'True' && r[params.observed_dp_col]
    }
    if (!usable) {
        error "No sample has both a CRAM and an ${params.observed_dp_col} in " +
              "${params.cram_info}; there is nothing to audit."
    }

    // Random subset per platform, seeded so the run is reproducible.
    def picked = []
    usable.groupBy { r -> r[params.platform_col] }.each { pf, rs ->
        def shuffled = new ArrayList(rs)
        Collections.shuffle(shuffled, new Random(params.sample_seed))
        def take = Math.min(params.samples_per_platform, shuffled.size())
        picked += shuffled.take(take)
        log.info "[audit]   ${pf}: ${take} of ${rs.size()}"
    }
    log.info "[audit] auditing ${picked.size()} samples over ${regionsCsv().split(',').size()} probe regions"

    // ── 00: measure the CRAMs ────────────────────────────────────────────────
    chunks_ch = channel
        .fromList(picked.collect { r ->
            "${r[params.sample_id_col]}\t${r[params.cram_path_col]}\n" })
        .collectFile(name: 'audit_list.txt', sort: true)
        .splitText(by: params.stats_chunk_size, file: true)

    stats_ch = CRAM_REGION_STATS(chunks_ch, stats_script).stats

    // ── 01: the test ─────────────────────────────────────────────────────────
    audit = OBSERVED_DEPTH_AUDIT(stats_ch.collect(), audit_script,
                                 file(params.cram_info, checkIfExists: true))

    // ── 02: the write-up ─────────────────────────────────────────────────────
    WRITE_REPORT(audit.tsv, audit.samples, audit.log, audit.figure, report_script)
}
