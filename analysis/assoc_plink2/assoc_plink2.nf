
nextflow.enable.dsl = 2

// =============================================================================
// assoc_plink2 — additive-primary association scan with two-tier peak calling
//
// 3 cohorts x 3 genetic models = 9 genome-wide scans. ADDITIVE is the primary
// and the ONLY model that defines peaks. Dominant and recessive are reported as
// scans — calibration, their own genome-wide loci, their own annotation rows and
// their own Manhattan labels — and nothing fans out from them.
//
// Peaks come in two tiers, called once at the suggestive threshold and then
// labelled, because the thresholds are nested:
//     genome_wide   lead P < 5e-8   -> full follow-up (ADDITIVE only)
//     suggestive    lead P < 1e-5   -> annotation table + display only
//
// Peaks are called for ALL THREE models so every scan figure can draw its own,
// but only the ADDITIVE peaks drive the fan-out below.
//
// params.Cohorts is an ORDERED list of sample sets. Where those sets are nested
// they are reported side by side and never ranked, and every peak lead is
// reported in EVERY cohort (CROSS_COHORT) — a nested design always has an
// estimate in the larger sets, and a blank cell is indistinguishable from a
// missing one.
//
// This is a fixed-effects scaffold: it establishes the scale of the signal and
// produces candidate peaks. Residual inflation it reports is fine-scale
// structure that global PCs do not absorb; the follow-up that addresses it is a
// GRM random-effect model (SAIGE / REGENIE), out of scope here. There is no
// power / minimum-detectable-effect component — without replication its only
// load-bearing statement is that a peak sits at the threshold, which is what P
// already says. See docs/METHODS.md.
//
// COST NOTE. RUN_ASSOC is by far the most expensive process. Its task hash
// depends on the interpolated script text and the input file contents, so any
// edit to that process block, to cohortInput(), or to the VALUES of the params
// they resolve (plink2, conda_env, FirthMode, out_dir, covarNameArg,
// covarTag) re-runs every scan. Everything else in this file is cheap to change.
// =============================================================================

// =============================================================================
// SITE CONFIGURATION — everything below is environment- or project-specific.
// Override on the command line (`--key value`) or with `-params-file cfg.yaml`;
// nothing outside this block needs to change to run the component elsewhere.
// =============================================================================

// -----------------------------------------------------------------------------
// Paths
// -----------------------------------------------------------------------------
params.project_dir   = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.base_dir      = "${params.project_dir}/analysis/assoc_plink2"
params.script_dir    = "${params.base_dir}/scripts"
params.out_dir       = "${params.base_dir}/results"
params.model_inputs  = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"

params.plink2        = '/home/b/b37974/plink2'
params.plink19       = '/home/b/b37974/plink'
params.tabix         = '/home/b/b37974/htslib-1.9/tabix'
// Prepended to PATH in the processes that shell out to helper binaries. Defaults to
// the directory holding plink2, which is where these tools normally sit together.
// EXCEPTION: RUN_ASSOC still carries this path as a literal. Its task hash is
// computed from its script block, so editing it re-runs every scan; a project
// starting fresh should switch that one line to ${params.tool_dir} too, and a
// project with an existing run should leave it alone. Nothing else in the file
// hard-codes a path.
params.tool_dir      = file(params.plink2).parent.toString()

params.conda_env     = 'cteph_geno_pro'      // python + plink runtime
params.conda_env_r   = 'r_work'              // single R runtime: susieR, EnsDb, rtracklayer
// Absolute interpreter, deliberately. `source activate r_work` on top of an
// already-active env leaves the FIRST env's Rscript ahead on PATH, which silently
// runs the R scripts against a library that has none of these packages.
params.rscript       = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'

// -----------------------------------------------------------------------------
// Design
// -----------------------------------------------------------------------------
// Sample sets, in the order they should be reported. Where they are nested,
// list them narrowest first — figures use this order top to bottom.
params.Cohorts       = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.Models        = ['additive', 'dominant', 'recessive']
params.PeakModel     = 'additive'            // the ONLY model that defines peaks

// Covariates: SEX + the first NPcs principal components. The covariate file is
// expected to carry PC1_AVG.. columns; --covar-name selects the first NPcs.
// PcLabel names the space those PCs were computed in and is used for the file
// name and for the human-readable CovarLabel printed on every figure.
params.PcLabel       = 'bbj_mainland'
params.NPcs          = 10
params.CovarFile     = "${params.PcLabel}_pc.sex.tsv"

params.FirthMode     = 'no-firth'            // '' to let plink2 fall back to Firth
params.PGenomeWide   = 5e-8
params.PSuggestive   = 1e-5
params.PeakFlank     = 250000                // half-width for distance-based merging
params.MaxCondRounds = 5
// How many suggestive peak leads a scan figure labels, smallest P first. All
// genome-wide leads are always labelled; a full suggestive tier does not fit a
// 7.2 in genomic axis, so the cap is stated in the caption rather than left for
// the reader to infer from a crowded panel.
params.LabelSuggestive = 10
// gwaslab annotation style for the scan figures. 'expand' repels label
// positions symmetrically and sets the text vertical; 'right' is gwaslab's own
// default greedy sweep at 40 degrees. RepelForce is the minimum separation as a
// fraction of the plotted span (gwaslab's default is 0.03).
params.AnnoStyle     = 'auto'
params.RepelForce    = 0.03

// How the covariate set is NAMED on figures and in sidecars. Derived from PcLabel
// and NPcs so it follows them by default; override when the covariates are not
// simply SEX + PCs. Figures used to hard-code this string, which meant a re-used
// pipeline reported covariates it had never fitted.
params.CovarLabel    = "SEX + ${params.NPcs} ${params.PcLabel} PCs"

// Display names for the three LD roles the regional figure draws. The KEYS are
// fixed identifiers for the roles; the VALUES name whichever resource this
// project configured above, and are the only place those names appear.
params.LdPanelLabels = ['cohort'   : 'in-sample cohort LD',
                        'tommo'    : 'ToMMo 54KJPN',
                        '1000g_eas': '1000 Genomes EAS (n = 504)']

// -----------------------------------------------------------------------------
// Fine-mapping / annotation resources
// -----------------------------------------------------------------------------
params.SusieL          = 10
params.SusieCoverage   = 0.95
params.SusieMinAbsCorr = 0.5

params.ref_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension'
params.PopR2Template = "${params.ref_dir}/ToMMo_60KJPN/co-occurrence/tommo-54kjpn-20230828-GRCh38-autosome-chr@@CHROM@@-plink-r2.tsv.gz"
// rsIDs: any tabix-indexed VCF on the same genome build whose ID column carries
// the rsID, looked up by chr:pos:REF:ALT. A dbSNP build works; so does any
// population AF VCF that carries IDs.
params.rsid_vcf      = "${params.ref_dir}/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
// Functional consequence comes from a PRE-COMPUTED snpEff index: one
// tabix-indexed TSV per chromosome, keyed chr:pos:REF:ALT. Nothing here runs
// snpEff itself; point this at whatever index the project already maintains.
params.snpeff_index_dir = '/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/output/snpEff.v6.index'
params.RefPanelBfile = "${params.ref_dir}/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all"
params.recomb_bw     = "${params.ref_dir}/cteph_agp3k/review_analysis/06.regional_plot_rev1/info/recomb1000GAvg.bw"

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
def covarNameArg() {
    return "SEX,PC1_AVG-PC${params.NPcs}_AVG"
}

def covarTag() {
    return "sex_pc1_${params.NPcs}_${params.PcLabel}"
}

def panelLabelArgs() {
    return params.LdPanelLabels.collect { k, v -> "--panel-label '${k}=${v}'" }.join(' ')
}

// plink2 --glm keyword for a model ('' = additive, plink2's default).
def modelFlag(String m) {
    return m == 'additive' ? '' : m
}

def cohortInput(String cohort) {
    def gdir = file("${params.model_inputs}/${cohort}/genotype/fixed_model")
    def bed  = gdir.list().findAll { it.endsWith('.maf_ge_threshold.bed') }
    if (!bed) error "No .maf_ge_threshold.bed under ${gdir} for cohort ${cohort}"
    def pfx  = "${gdir}/" + bed[0].replaceAll(/\.bed$/, '')
    return [ cohort,
             file("${pfx}.bed"), file("${pfx}.bim"), file("${pfx}.fam"),
             file("${params.model_inputs}/${cohort}/phenotype/pheno.tsv"),
             file("${params.model_inputs}/${cohort}/covariates/${params.CovarFile}") ]
}

// =============================================================================
// Processes
// =============================================================================

process RUN_ASSOC {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}" }, mode: 'symlink',
               pattern: '*.glm.logistic*'
    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}" }, mode: 'copy',
               pattern: '*.log'

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), path(pheno), path(covar), val(model)

    output:
    tuple val(cohort), val(model), path("*.glm.logistic*"), emit: glm
    path("*.log"),                                          emit: log

    script:
    def pfx  = bed.baseName
    def out  = "${cohort}.${model}.${covarTag()}"
    def mflg = modelFlag(model)
    """
    export PATH=/home/b/b37974/:\$PATH
    source activate ${params.conda_env}

    ${params.plink2} \
        --bfile ${pfx} \
        --pheno ${pheno} --pheno-name PHENO1 \
        --covar ${covar} --covar-name ${covarNameArg()} \
        --glm ${mflg} omit-ref ${params.FirthMode} hide-covar \
        --ci 0.95 \
        --threads ${task.ext.threads ?: 1} \
        --out ${out}
    """
}

process SCAN_PEAKS {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/02.scan" }, mode: 'copy', pattern: 'scan_qc.tsv'
    publishDir { "${params.out_dir}/${cohort}/03.peaks" }, mode: 'copy',
               pattern: '{peaks.tsv,lead_variants.tsv,model_peaks.tsv,*.sumstat.tsv}'

    input:
    tuple val(cohort), val(models), path(glms), path(fam)

    output:
    tuple val(cohort), path('scan_qc.tsv'),        emit: qc
    tuple val(cohort), path('peaks.tsv'),          emit: peaks
    tuple val(cohort), path('model_peaks.tsv'),    emit: model_peaks
    tuple val(cohort), path('lead_variants.tsv'),  emit: leads
    tuple val(cohort), path('*.sumstat.tsv'),      emit: sumstats, optional: true

    script:
    def pairs = [models, glms].transpose().collect { m, f -> "--glm ${m}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/call_peaks.py \
        ${pairs} \
        --cohort ${cohort} --fam ${fam} \
        --peak-model ${params.PeakModel} \
        --p-genomewide ${params.PGenomeWide} \
        --p-suggestive ${params.PSuggestive} \
        --peak-flank ${params.PeakFlank} \
        --out-dir .
    """
}

process PLOT_SCAN {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/figures/01.scan" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), path(glm), path(model_peaks)

    output:
    path("scan.${model}.*")

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_manhattan_qq.py \
        --glm ${glm} --cohort ${cohort} --model ${model} \
        --model-peaks ${model_peaks} \
        --label-suggestive ${params.LabelSuggestive} \
        --anno-style ${params.AnnoStyle} --repel-force ${params.RepelForce} \
        --covar-label '${params.CovarLabel}' --pc-label '${params.PcLabel}' --n-pcs ${params.NPcs} \
        --alpha ${params.PGenomeWide} --suggestive ${params.PSuggestive} \
        --out-png scan.${model}.png
    """
}

process ANNOTATE_LEADS {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(cohort), path(peaks), path(model_peaks), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path('lead_annotation.tsv'),        emit: annotation
    tuple val(cohort), path('model_peaks_annotation.tsv'), emit: model_peaks

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.script_dir}/annotate_leads.py \
        --peaks ${peaks} --model-peaks ${model_peaks} \
        --bfile ${bed.baseName} --cohort ${cohort} \
        --plink2 ${params.plink2} --tabix ${params.tabix} \
        --rsid-vcf ${params.rsid_vcf} \
        --rscript ${params.rscript} --gene-script ${params.script_dir}/gene_annotate.R \
        --snpeff-index ${params.snpeff_index_dir} \
        --threads ${task.ext.threads ?: 1} \
        --out lead_annotation.tsv --out-model-peaks model_peaks_annotation.tsv
    """
}

process CONDITIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(chrom), val(start), val(end), val(lead_id),
          path(bed), path(bim), path(fam), path(pheno), path(covar)

    output:
    tuple val(cohort), val(peak_id), path("${peak_id}.rounds.tsv"),
          path("${peak_id}.signals.tsv"), path("${peak_id}.round[0-9]*.tsv"), emit: cond

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.script_dir}/conditional_stepwise.py \
        --plink2 ${params.plink2} \
        --bfile ${bed.baseName} \
        --pheno ${pheno} --pheno-name PHENO1 \
        --covar ${covar} --covar-name ${covarNameArg()} \
        --model ${params.PeakModel} --firth-mode ${params.FirthMode} \
        --cohort ${cohort} --locus-id ${peak_id} \
        --chrom ${chrom} --start ${start} --end ${end} --lead-id ${lead_id} \
        --p-threshold ${params.PGenomeWide} --max-rounds ${params.MaxCondRounds} \
        --threads ${task.ext.threads ?: 1} \
        --out-dir .
    """
}

process LD_SOURCES {
    executor 'slurm'
    queue 'gr10478b'
    time '8h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(chrom), val(start), val(end), val(lead_id),
          path(sumstat), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), val(peak_id), val(lead_id), path(sumstat),
          path("${peak_id}.ld_matrix.tsv"), path("${peak_id}.ld_matrix.vars"),
          path("${peak_id}.ld_cohort.tsv"), path("${peak_id}.ld_tommo.tsv"),
          path("${peak_id}.ld_1000g_eas.tsv"), path("${peak_id}.ld_coverage.tsv"),
          path("${peak_id}.exons.tsv"), path("${peak_id}.recomb.tsv"), emit: ld

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.script_dir}/ld_sources.py \
        --plink ${params.plink19} --plink2 ${params.plink2} --tabix ${params.tabix} \
        --bfile ${bed.baseName} --ref-panel-bfile ${params.RefPanelBfile} \
        --pop-r2-template '${params.PopR2Template}' \
        ${panelLabelArgs()} \
        --sumstat ${sumstat} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --chrom ${chrom} --start ${start} --end ${end} \
        --memory-mb 16000 --threads ${task.ext.threads ?: 1} \
        --out-dir .

    ${params.rscript} ${params.script_dir}/region_tracks.R \
        --chrom ${chrom} --start ${start} --end ${end} \
        --recomb-bw ${params.recomb_bw} \
        --out-exons ${peak_id}.exons.tsv \
        --out-recomb ${peak_id}.recomb.tsv
    """
}

process SUSIE {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(sumstat),
          path(ld_matrix), path(ld_vars), val(n_gwas)

    output:
    tuple val(cohort), val(peak_id), val(lead_id),
          path("${peak_id}.pip.tsv"), path("${peak_id}.cs.tsv"),
          path("${peak_id}.susie.json"), emit: susie

    script:
    """
    ${params.rscript} ${params.script_dir}/susie_finemap.R \
        --sumstat ${sumstat} \
        --ld-matrix ${ld_matrix} --ld-vars ${ld_vars} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --n ${n_gwas} --L ${params.SusieL} \
        --coverage ${params.SusieCoverage} --min-abs-corr ${params.SusieMinAbsCorr} \
        --out-pip ${peak_id}.pip.tsv \
        --out-cs ${peak_id}.cs.tsv \
        --out-json ${peak_id}.susie.json
    """
}

process CS_VARIANTS {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(pip), path(cs), path(json), path(sumstat)

    output:
    path("${peak_id}.cs_variants.tsv"), emit: table

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.script_dir}/cs_variants.py \
        --pip ${pip} --cs ${cs} --sumstat ${sumstat} \
        --cohort ${cohort} --peak-id ${peak_id} --lead-id ${lead_id} \
        --tabix ${params.tabix} --rsid-vcf ${params.rsid_vcf} \
        --snpeff-index ${params.snpeff_index_dir} \
        --rscript ${params.rscript} --gene-script ${params.script_dir}/gene_annotate.R \
        --out ${peak_id}.cs_variants.tsv
    """
}

process COLLECT_CS {
    executor 'local'
    tag 'cs'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    input:
    path(tables, stageAs: 'cs_*.tsv')

    output:
    path('cs_variants_all.tsv')

    script:
    """
    source activate ${params.conda_env}
    python3 - <<'PYEOF'
import glob, pandas as pd
fs = sorted(glob.glob('cs_*.tsv'))
parts = [d for d in (pd.read_csv(f, sep='\t') for f in fs) if len(d)]
out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
if len(out):
    out = out.sort_values(['cohort', 'peak_id', 'cs', 'pip_rank'])
out.to_csv('cs_variants_all.tsv', sep='\t', index=False)
print(f'[collect_cs] {len(fs)} peak file(s) -> {len(out)} credible-set variants')
PYEOF
    """
}

process PLOT_REGIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/02.regional" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(sumstat),
          path(ld_cohort), path(ld_tommo), path(ld_eas), path(exons), path(recomb)

    output:
    path("regional.${peak_id}.*")

    script:
    """
    source activate ${params.conda_env}
    mkdir -p ld && cp ${ld_cohort} ${ld_tommo} ${ld_eas} ld/
    python3 ${params.script_dir}/plot_regional.py \
        --sumstat ${sumstat} --ld-dir ld \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --exons ${exons} --recomb ${recomb} --alpha ${params.PGenomeWide} \
        ${panelLabelArgs()} \
        --out-png regional.${peak_id}.png
    """
}

process PLOT_FINEMAP {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/03.finemap" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(pip), path(cs), path(json), path(ld_cohort)

    output:
    path("finemap.${peak_id}.*")

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_finemap.py \
        --pip ${pip} --cs ${cs} --json ${json} --ld-cohort ${ld_cohort} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --out-png finemap.${peak_id}.png
    """
}

process PLOT_CONDITIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/04.conditional" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(rounds), path(signals), path(round_files)

    output:
    path("conditional.${peak_id}.*")

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_conditional.py \
        --cond-dir . --rounds ${rounds} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --covar-label '${params.CovarLabel}' --pc-label '${params.PcLabel}' --n-pcs ${params.NPcs} \
        --alpha ${params.PGenomeWide} \
        --out-png conditional.${peak_id}.png
    """
}

process CROSS_COHORT {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    // All cohorts at once: the union of their peak leads is the variant list, so
    // this cannot be a per-cohort process. The three .bed/.bim/.fam share a file
    // name across cohorts, hence stageAs — the parallel numbering is what keeps
    // geno_N.bed/.bim/.fam a valid plink triple.
    input:
    tuple val(cohorts), path(peaks, stageAs: 'peaks_*.tsv'),
          path(glms,  stageAs: 'glm_*.tsv'),
          path(beds,  stageAs: 'geno_*.bed'),
          path(bims,  stageAs: 'geno_*.bim'),
          path(fams,  stageAs: 'geno_*.fam')

    output:
    path('lead_crosscohort.tsv'), emit: table

    script:
    def pk = [cohorts, peaks].transpose().collect { c, f -> "--peaks ${c}=${f}" }.join(' ')
    def gl = [cohorts, glms].transpose().collect  { c, f -> "--glm ${c}=${f}" }.join(' ')
    def bf = [cohorts, beds].transpose().collect  { c, f -> "--bfile ${c}=${f.baseName}" }.join(' ')
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.script_dir}/cross_cohort.py \
        ${pk} ${gl} ${bf} \
        --cohort-order ${params.Cohorts.join(',')} \
        --plink2 ${params.plink2} --tabix ${params.tabix} \
        --rsid-vcf ${params.rsid_vcf} \
        --rscript ${params.rscript} --gene-script ${params.script_dir}/gene_annotate.R \
        --snpeff-index ${params.snpeff_index_dir} \
        --threads ${task.ext.threads ?: 1} \
        --out lead_crosscohort.tsv
    """
}

process COMPARE_COHORTS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/tables",  mode: 'copy', pattern: '*.tsv'
    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy', pattern: 'cohort_compare.*'

    input:
    tuple val(cohorts), path(qc, stageAs: 'qc_*.tsv'), path(peaks, stageAs: 'peaks_*.tsv'),
          path(anns, stageAs: 'ann_*.tsv'), path(mpk, stageAs: 'mpk_*.tsv'),
          path(crosscohort)

    output:
    path('cohort_compare.*')
    path('scan_qc_all.tsv')
    path('peaks_all.tsv')
    path('lead_annotation_all.tsv')
    path('model_peaks_all.tsv')

    script:
    def q = [cohorts, qc].transpose().collect    { c, f -> "--scan-qc ${c}=${f}" }.join(' ')
    def p = [cohorts, peaks].transpose().collect { c, f -> "--peaks ${c}=${f}" }.join(' ')
    def a = [cohorts, anns].transpose().collect  { c, f -> "--annotation ${c}=${f}" }.join(' ')
    def h = [cohorts, mpk].transpose().collect   { c, f -> "--model-peaks ${c}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_cohort_compare.py \
        ${q} ${p} ${a} ${h} \
        --crosscohort ${crosscohort} \
        --cohort-order ${params.Cohorts.join(',')} \
        --covar-label '${params.CovarLabel}' --n-pcs ${params.NPcs} \
        --model ${params.PeakModel} \
        --alpha ${params.PGenomeWide} --alpha-suggestive ${params.PSuggestive} \
        --out-png cohort_compare.png \
        --out-scan scan_qc_all.tsv \
        --out-peaks peaks_all.tsv \
        --out-annotation lead_annotation_all.tsv \
        --out-model-peaks model_peaks_all.tsv
    """
}

process COMPARE_MANHATTAN {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy'

    // All cohorts at once: the rows share one genomic axis, one y-limit and one
    // data height, none of which can be decided from a single cohort. The three
    // .tsv inputs share a file name across cohorts, hence stageAs.
    input:
    tuple val(cohorts), path(glms, stageAs: 'glm_*.tsv'),
          path(qc,  stageAs: 'qc_*.tsv'),
          path(mpk, stageAs: 'mpk_*.tsv')

    output:
    path('cohort_manhattan.*')

    script:
    def g = [cohorts, glms].transpose().collect { c, f -> "--glm ${c}=${f}" }.join(' ')
    def q = [cohorts, qc].transpose().collect   { c, f -> "--scan-qc ${c}=${f}" }.join(' ')
    def h = [cohorts, mpk].transpose().collect  { c, f -> "--model-peaks ${c}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_cohort_manhattan.py \
        ${g} ${q} ${h} \
        --cohort-order ${params.Cohorts.join(',')} \
        --model ${params.PeakModel} \
        --alpha ${params.PGenomeWide} --suggestive ${params.PSuggestive} \
        --covar-label '${params.CovarLabel}' --pc-label '${params.PcLabel}' --n-pcs ${params.NPcs} \
        --anno-style ${params.AnnoStyle} --repel-force ${params.RepelForce} \
        --out-png cohort_manhattan.png
    """
}

process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'

    publishDir "${params.out_dir}/_run_info", mode: 'copy'

    input:
    val(_trigger)          // presence only — the manifest is built from params

    output:
    path('run_manifest.json')

    script:
    """
    cat > run_manifest.json <<'JSON'
{
  "component": "assoc_plink2",
  "cohorts": ${groovy.json.JsonOutput.toJson(params.Cohorts)},
  "models": ${groovy.json.JsonOutput.toJson(params.Models)},
  "peak_model": "${params.PeakModel}",
  "covariates": "${covarNameArg()}",
  "pc_label": "${params.PcLabel}",
  "n_pcs": ${params.NPcs},
  "firth_mode": "${params.FirthMode}",
  "p_genomewide": ${params.PGenomeWide},
  "p_suggestive": ${params.PSuggestive},
  "peak_flank_bp": ${params.PeakFlank},
  "label_suggestive": ${params.LabelSuggestive},
  "anno_style": "${params.AnnoStyle}",
  "repel_force": ${params.RepelForce},
  "susie": {"L": ${params.SusieL}, "coverage": ${params.SusieCoverage}, "min_abs_corr": ${params.SusieMinAbsCorr}},
  "model_inputs": "${params.model_inputs}",
  "plink2": "${params.plink2}",
  "plink19": "${params.plink19}",
  "rsid_vcf": "${params.rsid_vcf}",
  "snpeff_index_dir": "${params.snpeff_index_dir}",
  "ref_panel_bfile": "${params.RefPanelBfile}",
  "pop_r2_template": "${params.PopR2Template}",
  "recomb_bw": "${params.recomb_bw}",
  "conda_env": "${params.conda_env}",
  "conda_env_r": "${params.conda_env_r}",
  "rscript": "${params.rscript}"
}
JSON
    """
}

// =============================================================================
// Workflow
// =============================================================================
workflow {

    ch_cohort = channel.fromList(params.Cohorts).map { cohortInput(it) }
    ch_geno   = ch_cohort.map { c, bed, bim, fam, _ph, _cv -> tuple(c, bed, bim, fam) }
    ch_fam    = ch_cohort.map { c, _bed, _bim, fam, _ph, _cv -> tuple(c, fam) }

    // ── 9 genome-wide scans (cached; see the header) ───────────────────────
    RUN_ASSOC(ch_cohort.combine(channel.fromList(params.Models)))

    // plink2 writes ONE association file (hide-covar leaves a single TEST), so
    // collapse the emitted list to that file before anything reads it.
    ch_glm = RUN_ASSOC.out.glm.map { c, m, f -> tuple(c, m, f instanceof List ? f[0] : f) }

    // ── per cohort: QC + sample sizes for all three models, peaks from the
    //    additive one, genome-wide hits from all three ─────────────────────
    SCAN_PEAKS(ch_glm.groupTuple(by: 0, size: params.Models.size())
                     .combine(ch_fam, by: 0))

    // ── annotation: additive peak leads (both tiers) AND every model's
    //    genome-wide hits, in one pass over the union ──────────────────────
    ANNOTATE_LEADS(SCAN_PEAKS.out.peaks
        .combine(SCAN_PEAKS.out.model_peaks, by: 0)
        .combine(ch_geno, by: 0))

    // ── one figure per scan: Manhattan + its own peaks + QQ + effect/MAF ────
    // Each scan draws ITS OWN peaks, gene-labelled, so this runs after
    // ANNOTATE_LEADS — the symbols live in model_peaks_annotation.tsv.
    PLOT_SCAN(ch_glm.combine(ANNOTATE_LEADS.out.model_peaks, by: 0))

    // ── fan out over GENOME-WIDE peaks only ────────────────────────────────
    // Suggestive peaks get the annotation table and the landscape figure;
    // they receive no per-peak follow-up, so the filter happens here and no
    // downstream process needs to know about tiers.
    ch_peak = SCAN_PEAKS.out.leads
        .splitCsv(header: true, sep: '\t', elem: 1)
        .filter { _c, r -> r.tier == 'genome_wide' }
        .map { c, r -> tuple(c, r.peak_id, r.chrom, r.start as long, r.end as long, r.lead_id) }

    CONDITIONAL(ch_peak.combine(ch_cohort, by: 0))

    // Per-peak summary statistics arrive as a per-cohort list; key each file by
    // its own peak_id so the joins below cannot pair the wrong peak.
    ch_sumstat = SCAN_PEAKS.out.sumstats
        .flatMap { c, fs -> (fs instanceof List ? fs : [fs])
                            .collect { tuple(c, it.name.replaceAll(/\.sumstat\.tsv$/, ''), it) } }

    LD_SOURCES(ch_peak.combine(ch_sumstat, by: [0, 1]).combine(ch_geno, by: 0))

    // GWAS analysis N per cohort — SuSiE's `n` is the sample size the summary
    // statistics come from, NOT the size of the LD reference.
    ch_n = ch_cohort.map { c, _bed, _bim, fam, _ph, _cv ->
        tuple(c, fam.readLines().count { it.trim() && (it.split(/\s+/)[5] in ['1', '2']) })
    }

    SUSIE(LD_SOURCES.out.ld
        .map { c, pid, lead, ss, mat, vars, _ldc, _ldt, _lde, _cov, _ex, _rec ->
               tuple(c, pid, lead, ss, mat, vars) }
        .combine(ch_n, by: 0))

    ch_sumstat_by_peak = LD_SOURCES.out.ld
        .map { c, pid, _lead, ss, _mat, _vars, _ldc, _ldt, _lde, _cov, _ex, _rec ->
               tuple(c, pid, ss) }
    CS_VARIANTS(SUSIE.out.susie.combine(ch_sumstat_by_peak, by: [0, 1]))
    COLLECT_CS(CS_VARIANTS.out.table.collect(sort: true))

    PLOT_REGIONAL(LD_SOURCES.out.ld
        .map { c, pid, lead, ss, _mat, _vars, ldc, ldt, lde, _cov, ex, rec ->
               tuple(c, pid, lead, ss, ldc, ldt, lde, ex, rec) })

    ch_ldc = LD_SOURCES.out.ld
        .map { c, pid, _lead, _ss, _mat, _vars, ldc, _ldt, _lde, _cov, _ex, _rec ->
               tuple(c, pid, ldc) }
    PLOT_FINEMAP(SUSIE.out.susie.combine(ch_ldc, by: [0, 1]))

    PLOT_CONDITIONAL(CONDITIONAL.out.cond
        .combine(ch_peak.map { c, pid, _ch, _s, _e, lead -> tuple(c, pid, lead) }, by: [0, 1])
        .map { c, pid, rounds, signals, rfiles, lead ->
               tuple(c, pid, lead, rounds, signals, rfiles) })

    // ── every peak lead of every cohort, reported in EVERY cohort ──────────
    // A barrier, deliberately: the variant list is the UNION over cohorts, so no
    // cohort can be processed until all of them have called their peaks.
    ch_cross = SCAN_PEAKS.out.peaks
        .combine(ch_glm.filter { _c, m, _f -> m == params.PeakModel }
                       .map { c, _m, f -> tuple(c, f) }, by: 0)
        .combine(ch_geno, by: 0)
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] }, o.collect { it[2] },
                  o.collect { it[3] }, o.collect { it[4] }, o.collect { it[5] })
        }
    CROSS_COHORT(ch_cross)

    // ── one cross-cohort comparison figure, additive only ──────────────────
    ch_all = SCAN_PEAKS.out.qc
        .combine(SCAN_PEAKS.out.peaks, by: 0)
        .combine(ANNOTATE_LEADS.out.annotation, by: 0)
        .combine(ANNOTATE_LEADS.out.model_peaks, by: 0)
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] }, o.collect { it[2] },
                  o.collect { it[3] }, o.collect { it[4] })
        }
    COMPARE_COHORTS(ch_all.combine(CROSS_COHORT.out.table))

    // ── the three additive scans stacked on ONE genomic axis ───────────────
    // A barrier for the same reason as CROSS_COHORT: the shared offset map, the
    // shared y-limit and the union of genome-wide loci are all properties of the
    // three cohorts together, so none of them exists until all three have run.
    ch_man = ch_glm.filter { _c, m, _f -> m == params.PeakModel }
        .map { c, _m, f -> tuple(c, f) }
        .combine(SCAN_PEAKS.out.qc, by: 0)
        .combine(ANNOTATE_LEADS.out.model_peaks, by: 0)
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] },
                  o.collect { it[2] }, o.collect { it[3] })
        }
    COMPARE_MANHATTAN(ch_man)

    WRITE_RUN_MANIFEST(channel.of('go'))
}
