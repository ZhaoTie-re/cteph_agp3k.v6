
nextflow.enable.dsl = 2

// =============================================================================
// assoc_saige — mixed-model association scan with two-tier peak calling
//
// The random-effect counterpart of assoc_plink2. Same design, same peak rule,
// same figures, same tables — a full genetic relationship matrix replaces the
// fixed set of principal components as the way population structure and
// relatedness are absorbed, so the RANDOM-MODEL genotypes are used and relateds
// are retained rather than removed.
//
// N sample sets x 2 covariate models = the scans. The two models are:
//     main     the reported analysis      SEX + params.NPcs PCs, full GRM
//     detect   a calibration probe        SEX only, no PCs,   full GRM
// `detect` exists to answer one question — how much of the structure the GRM
// absorbs on its own, and how much the PCs still remove on top of it. It is
// reported (its own scan figure, its own rows in model_peaks.tsv, the
// model-comparison figure) and NOTHING fans out from it, exactly the treatment
// assoc_plink2 gives its non-primary genetic models.
//
// Peaks come in two tiers, called once at the suggestive threshold and then
// labelled, because the thresholds are nested:
//     genome_wide   lead P < PGenomeWide   -> full per-peak follow-up
//     suggestive    lead P < PSuggestive   -> the same, judged against PSuggestive
// Both tiers fan out; outputs are separated under <tier>/ so they never mix.
// Only the PEAK MODEL defines peaks at all — that restriction is unchanged.
//
// SAIGE emits BOTH p-values in a single pass — p.value (saddlepoint-corrected
// where Is.SPA) and p.value.NA (normal approximation). params.PValue therefore
// selects between them in TO_SUMSTATS, downstream of every expensive stage:
// switching it re-runs the adapter and the figures and NEVER re-runs SAIGE.
// Which one was used is recorded in results/_run_info/run_manifest.json.
//
// COST NOTE. FIT_NULL and RUN_ASSOC are by far the most expensive processes; a
// full GRM on a few thousand samples is hours, and there is one per sample set
// per model. Their task hashes are computed from their script blocks, so an edit
// there — including adding a comment, since Nextflow hashes that text and `//`
// is not a shell comment — re-runs the scan. Everything presentational is
// consumed downstream of TO_SUMSTATS and is free to change.
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
params.base_dir      = "${params.project_dir}/analysis/assoc_saige"
// Engine-specific scripts (this component only).
params.script_dir    = "${params.base_dir}/scripts"
// Scripts shared with the other association components — peak calling, annotation,
// LD, fine-mapping and every figure. They are about association RESULTS, not about
// which engine produced them, so they live once at analysis/_shared.
params.shared_script_dir = "${params.project_dir}/analysis/_shared/scripts"
params.out_dir       = "${params.base_dir}/results"
params.model_inputs  = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"
params.info_dir      = "${params.project_dir}/info"

params.plink2        = '/home/b/b37974/plink2'
params.plink19       = '/home/b/b37974/plink'
params.tabix         = '/home/b/b37974/htslib-1.9/tabix'
// Prepended to PATH in the processes that shell out to helper binaries.
params.tool_dir      = file(params.plink2).parent.toString()

params.conda_env     = 'cteph_geno_pro'      // python + plink runtime
params.conda_env_r   = 'r_work'              // susieR, EnsDb, rtracklayer
// Absolute interpreters, deliberately. `source activate <env>` on top of an
// already-active env leaves the FIRST env's Rscript ahead on PATH, which silently
// runs against a library that has none of the required packages. This is not
// hypothetical: SAIGE's own step scripts start `#!/usr/bin/env Rscript`.
params.rscript       = '/home/b/b37974/anaconda3/envs/r_work/bin/Rscript'
params.saige_rscript = '/home/b/b37974/anaconda3/envs/saige/bin/Rscript'
params.saige_step1   = '/home/b/b37974/anaconda3/envs/saige/bin/step1_fitNULLGLMM.R'
params.saige_step2   = '/home/b/b37974/anaconda3/envs/saige/bin/step2_SPAtests.R'

// -----------------------------------------------------------------------------
// Design
// -----------------------------------------------------------------------------
// Sample sets, in the order they should be reported. Where they are nested,
// list them narrowest first — figures use this order top to bottom.
params.Cohorts       = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']

// The covariate models. KEYS name the model everywhere downstream; VALUES are the
// SAIGE --covarColList for it. `detect` deliberately drops the PCs so the GRM is
// the only thing absorbing structure.
params.NPcs          = 10
params.PcLabel       = 'bbj_mainland'
params.Models        = ['main', 'detect']
params.PeakModel     = 'main'                // the ONLY model that defines peaks
params.ModelCovars   = ['main'  : (['SEX'] + (1..params.NPcs).collect { "PC${it}_AVG" }).join(','),
                        'detect': 'SEX']
params.CovarFile     = "${params.PcLabel}_pc.sex.tsv"
params.PhenoName     = 'PHENO1'
params.SampleIdCol   = 'IID'

// How each model's covariate set is NAMED on figures and in sidecars.
params.ModelLabels   = ['main'  : "SEX + ${params.NPcs} ${params.PcLabel} PCs, full GRM",
                        'detect': 'SEX only, full GRM']

// Which SAIGE p-value the whole downstream uses. Applied in TO_SUMSTATS, so
// changing it never re-runs SAIGE. Recorded in run_manifest.json.
//
// DEFAULT: no_spa — the normal-approximation score test (SAIGE's p.value.NA).
// The saddlepoint correction exists to control the tail when case/control
// imbalance is severe; this study's imbalance was judged not to warrant it. The
// choice is consequential rather than cosmetic: it decides whether the strongest
// locus clears 5e-8, so it is recorded in the manifest and stated in METHODS.
params.PValue        = 'no_spa'              // 'spa' | 'no_spa'

params.PGenomeWide   = 5e-8
params.PSuggestive   = 1e-5
params.PeakFlank     = 250000                // half-width for distance-based merging
params.MaxCondRounds = 5
// How many suggestive peaks the scan figure NAMES *and* the fan-out FOLLOWS UP,
// smallest P first. Genome-wide peaks are always both. One parameter, not two,
// because the two must agree: at 1e-5 the suggestive tier holds tens of peaks
// per cohort against a comparable number expected by chance, so the tier is a
// description of the scan's shape rather than a list of findings, and only its
// strongest members earn a regional / fine-map / conditional figure.
//
// THE INVARIANT this buys: a suggestive locus named on a scan figure always has
// downstream figures, and one that is not named never does. Splitting this into
// separate label and follow-up counts breaks that silently, which is why it is
// passed to call_peaks.py and plot_manhattan_qq.py from the same place.
params.MaxSuggestive = 10
params.AnnoStyle     = 'auto'
params.RepelForce    = 0.03

// -----------------------------------------------------------------------------
// GRM construction
// -----------------------------------------------------------------------------
// Markers for the full GRM: LD-pruned, with the known high-LD regions excluded.
// The GRM describes relatedness and structure, so correlated and inversion-prone
// regions would let a handful of loci dominate it.
params.HighLdFile    = "${params.info_dir}/high-LD-regions-hg38-GRCh38_modified.txt"
params.PruneWindow   = 500                   // kb
params.PruneStep     = 50
params.PruneR2       = 0.2
params.VarianceRatioMarkers = 200
params.Autosomes     = (1..22).collect { it as String }

// Reproducibility. Nothing here samples without a seed: the GRM is built from
// EVERY pruned marker, plink2's pruning is deterministic, and SAIGE's step 1
// calls set.seed() itself before drawing the markers it estimates the variance
// ratio from. RandomSeed RECORDS that value for provenance — SAIGE exposes no
// option to override it, so changing this does not change the fit.
params.RandomSeed    = 1

// Relatedness read-out: report pairs at or above this GRM RELATIONSHIP value.
// 0.0884 is the third-degree boundary on the GRM scale (twice the 0.0442 KING
// kinship boundary, because a GRM relationship is 2 x kinship).
params.MinRelationship = 0.0884
// The relatedness-pruned sample set to compare against, under model_inputs —
// what a fixed-effects design must drop and this one keeps. Blank disables it.
params.CompareGenotypeDir = 'fixed_model'

// -----------------------------------------------------------------------------
// Fine-mapping / annotation resources  (identical to assoc_plink2)
// -----------------------------------------------------------------------------
params.SusieL          = 10
params.SusieCoverage   = 0.95
params.SusieMinAbsCorr = 0.5

params.ref_dir       = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension'
params.PopR2Template = "${params.ref_dir}/ToMMo_60KJPN/co-occurrence/tommo-54kjpn-20230828-GRCh38-autosome-chr@@CHROM@@-plink-r2.tsv.gz"
params.rsid_vcf      = "${params.ref_dir}/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
params.snpeff_index_dir = '/LARGE1/gr10478/platform/JHRPv6/workspace/pipeline/output/snpEff.v6.index'
params.RefPanelBfile = "${params.ref_dir}/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all"
params.recomb_bw     = "${params.ref_dir}/cteph_agp3k/review_analysis/06.regional_plot_rev1/info/recomb1000GAvg.bw"
params.LdPanelLabels = ['cohort'   : 'in-sample cohort LD',
                        'tommo'    : 'ToMMo 54KJPN',
                        '1000g_eas': '1000 Genomes EAS (n = 504)']

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
// Figure family -> its directory under figures/. The catalogue publishes beside
// the figures it documents, so the two must agree.
params.FigureDirs = ['regional': '02.regional', 'finemap': '03.finemap',
                     'conditional': '04.conditional']

def panelLabelArgs() {
    return params.LdPanelLabels.collect { k, v -> "--panel-label '${k}=${v}'" }.join(' ')
}

// The tier a peak belongs to, derived from its id. call_peaks.py builds peak ids
// as `gw###_<chrom>_<pos>` / `sg###_<chrom>_<pos>`, and that prefix is the only
// tier information reachable INSIDE a process: a publishDir closure sees process
// inputs, and adding `tier` as an input would change every task hash and re-run
// the entire existing fan-out purely to name a directory. Asserts rather than
// guessing, so a change to the id convention fails loudly here.
def tierOf(String peak_id) {
    if (peak_id.startsWith('gw')) return 'genome_wide'
    if (peak_id.startsWith('sg')) return 'suggestive'
    error "cannot determine the tier of peak '${peak_id}' — expected the gw/sg prefix call_peaks.py writes"
}

// Conditional analysis is judged against the threshold that DEFINED the peak.
// Judging a suggestive peak against the genome-wide line would make round 0 fail
// for every one of them, and the table would report zero signals everywhere —
// meaningless rather than empty.
def peakThreshold(String peak_id) {
    return tierOf(peak_id) == 'genome_wide' ? params.PGenomeWide : params.PSuggestive
}

// How many PCs a model actually fits, and the space they came from. `detect`
// fits none, so its figures must not name a PC space it never used.
def nPcsFor(String model) {
    return covarsFor(model).split(',').findAll { it ==~ /(?i)PC\d+_AVG/ }.size()
}

def pcLabelFor(String model) {
    return nPcsFor(model) > 0 ? params.PcLabel : ''
}

def covarsFor(String model) {
    def c = params.ModelCovars[model]
    if (!c) error "No covariate list configured for model '${model}' (params.ModelCovars)"
    return c
}

// The random-model genotypes: relateds RETAINED, because the GRM is what absorbs
// them. This is the difference from assoc_plink2, which uses fixed_model.
def compareFam(String cohort) {
    if (!params.CompareGenotypeDir) return null
    def d = file("${params.model_inputs}/${cohort}/genotype/${params.CompareGenotypeDir}")
    if (!d.exists()) return null
    def f = d.list().findAll { it.endsWith('.fam') }
    return f ? file("${d}/${f[0]}") : null
}

def cohortInput(String cohort) {
    def gdir = file("${params.model_inputs}/${cohort}/genotype/random_model")
    def bed  = gdir.list().findAll { it.endsWith('.bed') }
    if (!bed) error "No random_model .bed under ${gdir} for cohort ${cohort}"
    def pfx   = bed[0] - '.bed'
    def pheno = file("${params.model_inputs}/${cohort}/phenotype/pheno.tsv")
    def covar = file("${params.model_inputs}/${cohort}/covariates/${params.CovarFile}")
    if (!pheno.exists()) error "Missing phenotype for ${cohort}: ${pheno}"
    if (!covar.exists()) error "Missing covariates for ${cohort}: ${covar}"
    return tuple(cohort,
                 file("${gdir}/${pfx}.bed"), file("${gdir}/${pfx}.bim"), file("${gdir}/${pfx}.fam"),
                 pheno, covar)
}

// =============================================================================
// Stage 1 — inputs the GRM and the scan are built from (per cohort)
// =============================================================================

process PREP_PHENO_COV {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), path(pheno), path(covar)

    output:
    tuple val(cohort), path('pheno_covar.tsv'), emit: table

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/merge_pheno_cov.py \
        --pheno ${pheno} --cov ${covar} \
        --pheno-col ${params.PhenoName} \
        --id-col ${params.SampleIdCol} \
        --keep-fam ${fam} \
        --out pheno_covar.tsv
    """
}

process PRUNE_MARKERS {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy', pattern: '*.prune.in'

    input:
    tuple val(cohort), path(bed), path(bim), path(fam)

    output:
    tuple val(cohort), path('grm_markers.bed'), path('grm_markers.bim'),
          path('grm_markers.fam'), path('grm.prune.in'), emit: markers

    script:
    def hi = params.HighLdFile
    """
    export PATH=${params.tool_dir}:\$PATH
    set -euo pipefail
    awk 'BEGIN{OFS="\\t"} !/^#/ && NF>=3 {c=\$1; sub(/^chr/,"",c); print c, \$2, \$3, "highld"}' \
        ${hi} > high_ld.bed
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \
        --exclude bed1 high_ld.bed \
        --indep-pairwise ${params.PruneWindow} ${params.PruneStep} ${params.PruneR2} \
        --out grm --threads ${task.ext.threads ?: 1}
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \
        --extract grm.prune.in --make-bed --out grm_markers \
        --threads ${task.ext.threads ?: 1}
    """
}

process SPLIT_BY_CHROM {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:chr${chrom}"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), val(chrom)

    output:
    tuple val(cohort), val(chrom), path("chr${chrom}.bed"), path("chr${chrom}.bim"),
          path("chr${chrom}.fam"), emit: geno

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \
        --chr ${chrom} --make-bed --out chr${chrom} \
        --threads ${task.ext.threads ?: 1}
    """
}

// =============================================================================
// Stage 2 — SAIGE.  EXPENSIVE.  Do not edit these two script blocks casually.
// =============================================================================

process FIT_NULL {
    executor 'slurm'
    queue 'gr10478b'
    time '72h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), path(gbed), path(gbim), path(gfam), path(pheno_covar)

    output:
    tuple val(cohort), val(model), path('null.rda'), path('null.varianceRatio.txt'), emit: nullModel
    path('null.fit.log')

    script:
    def covars = covarsFor(model)
    """
    set -euo pipefail
    ${params.saige_rscript} ${params.saige_step1} \
        --plinkFile=${gbed.baseName} \
        --phenoFile=${pheno_covar} \
        --phenoCol=${params.PhenoName} \
        --covarColList=${covars} \
        --sampleIDColinphenoFile=${params.SampleIdCol} \
        --traitType=binary \
        --outputPrefix=null \
        --nThreads=${task.ext.threads ?: 1} \
        --isDiagofKinSetAsOne=True \
        --numRandomMarkerforVarianceRatio=${params.VarianceRatioMarkers} \
        --skipVarianceRatioEstimation=FALSE \
        --useSparseGRMtoFitNULL=FALSE \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE \
        --isCovariateOffset=FALSE 2>&1 | tee null.fit.log
    """
}

process RUN_ASSOC {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "${cohort}:${model}:chr${chrom}"

    input:
    tuple val(cohort), val(model), val(chrom), path(bed), path(bim), path(fam),
          path(rda), path(vr)

    output:
    tuple val(cohort), val(model), val(chrom), path("chr${chrom}.assoc.txt"), emit: assoc

    script:
    """
    set -euo pipefail
    ${params.saige_rscript} ${params.saige_step2} \
        --bedFile=${bed} \
        --bimFile=${bim} \
        --famFile=${fam} \
        --AlleleOrder=alt-first \
        --chrom=${chrom} \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${vr} \
        --SAIGEOutputFile=chr${chrom}.assoc.txt \
        --is_Firth_beta=FALSE \
        --LOCO=TRUE \
        --is_output_moreDetails=TRUE
    """
}

// =============================================================================
// Stage 3 — merge, then adapt to the canonical schema.  CHEAP: everything
// presentational, including the SPA choice, is decided from here down.
// =============================================================================

process MERGE_ASSOC {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), val(chroms), path(parts, stageAs: 'part_*.txt')

    output:
    tuple val(cohort), val(model), path('saige.assoc.tsv'), emit: merged

    script:
    def ordered = [chroms, parts].transpose()
        .sort { a, b -> (a[0] as Integer) <=> (b[0] as Integer) }
        .collect { it[1].name }.join(' ')
    """
    set -euo pipefail
    first=\$(echo ${ordered} | awk '{print \$1}')
    head -1 \$first > saige.assoc.tsv
    for f in ${ordered}; do tail -n +2 \$f >> saige.assoc.tsv; done
    """
}

process TO_SUMSTATS {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/01.assoc/${model}" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), path(merged), path(bim)

    output:
    tuple val(cohort), val(model), path('sumstats.tsv'), emit: sumstats
    path('saige_extra.tsv')

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/saige_to_sumstats.py \
        --saige ${merged} --bim ${bim} \
        --p-value ${params.PValue} \
        --cohort ${cohort} --model ${model} \
        --out-sumstats sumstats.tsv --out-extra saige_extra.tsv
    """
}

process RELATEDNESS {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(gbed), path(gbim), path(gfam), path(fam), path(cmpfam)

    output:
    tuple val(cohort), path('relatedness_summary.tsv'), path('relatedness_pairs.tsv'),
          path('relatedness_hist.tsv'), path('relatedness_diag.tsv'), emit: rel

    script:
    def cmp = cmpfam.name == 'NO_FILE' ? '' : "--compare-fam ${cmpfam}"
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    ${params.plink2} --bfile ${gbed.baseName} --make-rel triangle bin4 \
        --out grm --threads ${task.ext.threads ?: 1}
    python3 ${params.script_dir}/relatedness.py \
        --rel-bin grm.rel.bin --rel-id grm.rel.id \
        --cohort ${cohort} --fam ${fam} ${cmp} \
        --min-relationship ${params.MinRelationship} \
        --out-pairs relatedness_pairs.tsv \
        --out-summary relatedness_summary.tsv \
        --out-hist relatedness_hist.tsv \
        --out-diag relatedness_diag.tsv
    """
}

process COMPARE_RELATEDNESS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy', pattern: 'relatedness.*'
    publishDir "${params.out_dir}/_comparison/tables",  mode: 'copy', pattern: '*_all.tsv'

    input:
    tuple val(cohorts), path(summ, stageAs: 'sum_*.tsv'),
          path(prs, stageAs: 'pairs_*.tsv'), path(hst, stageAs: 'hist_*.tsv'),
          path(dia, stageAs: 'diag_*.tsv')

    output:
    path('relatedness.*')
    path('relatedness_all.tsv')

    script:
    def s = [cohorts, summ].transpose().collect { c, f -> "--summary ${c}=${f}" }.join(' ')
    def r = [cohorts, prs].transpose().collect  { c, f -> "--pairs ${c}=${f}" }.join(' ')
    def h = [cohorts, hst].transpose().collect  { c, f -> "--hist ${c}=${f}" }.join(' ')
    def d = [cohorts, dia].transpose().collect  { c, f -> "--diag ${c}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_relatedness.py \
        ${s} ${r} ${h} ${d} \
        --cohort-order ${params.Cohorts.join(',')} \
        --min-relationship ${params.MinRelationship} \
        --out-png relatedness.png --out-table relatedness_all.tsv
    """
}

// =============================================================================
// Stage 4 — shared downstream.  Identical to assoc_plink2 from here on, because
// it reads the canonical schema TO_SUMSTATS produces.
// =============================================================================

process SCAN_PEAKS {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}"

    publishDir { "${params.out_dir}/${cohort}/02.scan" }, mode: 'copy', pattern: 'scan_qc.tsv'
    publishDir { "${params.out_dir}/${cohort}/03.peaks" }, mode: 'copy',
               pattern: '{peaks.tsv,lead_variants.tsv,model_peaks.tsv,*.sumstat.tsv}'

    // Every model's adapter output is called sumstats.tsv, so they must be staged
    // under distinct names; the numbering follows the input list order, which is
    // the same order `models` is in.
    input:
    tuple val(cohort), val(models), path(sumstats, stageAs: 'ss_*.tsv'), path(fam)

    output:
    tuple val(cohort), path('scan_qc.tsv'),        emit: qc
    tuple val(cohort), path('peaks.tsv'),          emit: peaks
    tuple val(cohort), path('model_peaks.tsv'),    emit: model_peaks
    tuple val(cohort), path('lead_variants.tsv'),  emit: leads
    tuple val(cohort), path('*.sumstat.tsv'),      emit: sumstats, optional: true

    script:
    def pairs = [models, sumstats].transpose().collect { m, f -> "--glm ${m}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/call_peaks.py \
        ${pairs} \
        --cohort ${cohort} --fam ${fam} \
        --peak-model ${params.PeakModel} \
        --p-genomewide ${params.PGenomeWide} \
        --p-suggestive ${params.PSuggestive} \
        --max-suggestive-followup ${params.MaxSuggestive} \
        --peak-flank ${params.PeakFlank} \
        --out-dir .
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
    python3 ${params.shared_script_dir}/annotate_leads.py \
        --peaks ${peaks} --model-peaks ${model_peaks} \
        --bfile ${bed.baseName} --cohort ${cohort} \
        --p-genomewide ${params.PGenomeWide} --p-suggestive ${params.PSuggestive} \
        --plink2 ${params.plink2} --tabix ${params.tabix} \
        --rsid-vcf ${params.rsid_vcf} \
        --rscript ${params.rscript} --gene-script ${params.shared_script_dir}/gene_annotate.R \
        --snpeff-index ${params.snpeff_index_dir} \
        --threads ${task.ext.threads ?: 1} \
        --out lead_annotation.tsv --out-model-peaks model_peaks_annotation.tsv
    """
}

process PLOT_SCAN {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${cohort}:${model}"

    publishDir { "${params.out_dir}/${cohort}/figures/01.scan" }, mode: 'copy'

    input:
    tuple val(cohort), val(model), path(sumstat), path(model_peaks)

    output:
    path("scan.${model}.*")

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/plot_manhattan_qq.py \
        --glm ${sumstat} --cohort ${cohort} --model ${model} \
        --model-peaks ${model_peaks} \
        --label-suggestive ${params.MaxSuggestive} \
        --anno-style ${params.AnnoStyle} --repel-force ${params.RepelForce} \
        --covar-label '${params.ModelLabels[model]}' \
        --pc-label '${pcLabelFor(model)}' --n-pcs ${nPcsFor(model)} \
        --alpha ${params.PGenomeWide} --suggestive ${params.PSuggestive} \
        --peak-flank ${params.PeakFlank} \
        --out-png scan.${model}.png
    """
}

process CONDITIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${tierOf(peak_id)}/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(chrom), val(start), val(end), val(lead_id),
          path(bed), path(bim), path(fam), path(rda), path(vr)

    output:
    tuple val(cohort), val(peak_id), path("${peak_id}.rounds.tsv"),
          path("${peak_id}.signals.tsv"), path("${peak_id}.round[0-9]*.tsv"), emit: cond

    script:
    """
    source activate ${params.conda_env}
    export PATH=${params.tool_dir}:\$PATH
    python3 ${params.script_dir}/saige_conditional.py \
        --saige-rscript ${params.saige_rscript} --saige-step2 ${params.saige_step2} \
        --plink2 ${params.plink2} --threads ${task.ext.threads ?: 1} \
        --bed ${bed} --bim ${bim} --fam ${fam} \
        --rda ${rda} --variance-ratio ${vr} \
        --cohort ${cohort} --locus-id ${peak_id} \
        --chrom ${chrom} --start ${start} --end ${end} --lead-id ${lead_id} \
        --p-value ${params.PValue} \
        --p-threshold ${peakThreshold(peak_id)} --max-rounds ${params.MaxCondRounds} \
        --out-dir .
    """
}

process LD_SOURCES {
    executor 'slurm'
    queue 'gr10478b'
    time '8h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${tierOf(peak_id)}/${peak_id}" }, mode: 'copy'

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
    python3 ${params.shared_script_dir}/ld_sources.py \
        --plink ${params.plink19} --plink2 ${params.plink2} --tabix ${params.tabix} \
        --bfile ${bed.baseName} --ref-panel-bfile ${params.RefPanelBfile} \
        --pop-r2-template '${params.PopR2Template}' \
        ${panelLabelArgs()} \
        --sumstat ${sumstat} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --chrom ${chrom} --start ${start} --end ${end} \
        --memory-mb 16000 --threads ${task.ext.threads ?: 1} \
        --out-dir .

    ${params.rscript} ${params.shared_script_dir}/region_tracks.R \
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

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${tierOf(peak_id)}/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(sumstat),
          path(ld_matrix), path(ld_vars), val(n_gwas)

    output:
    tuple val(cohort), val(peak_id), val(lead_id),
          path("${peak_id}.pip.tsv"), path("${peak_id}.cs.tsv"),
          path("${peak_id}.susie.json"), emit: susie

    script:
    """
    ${params.rscript} ${params.shared_script_dir}/susie_finemap.R \
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

    publishDir { "${params.out_dir}/${cohort}/03.peaks/${tierOf(peak_id)}/${peak_id}" }, mode: 'copy'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(pip), path(cs), path(json), path(sumstat)

    output:
    path("${peak_id}.cs_variants.tsv"), emit: table

    script:
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/cs_variants.py \
        --pip ${pip} --cs ${cs} --sumstat ${sumstat} \
        --cohort ${cohort} --peak-id ${peak_id} --lead-id ${lead_id} \
        --tabix ${params.tabix} --rsid-vcf ${params.rsid_vcf} \
        --snpeff-index ${params.snpeff_index_dir} \
        --rscript ${params.rscript} --gene-script ${params.shared_script_dir}/gene_annotate.R \
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

    publishDir { "${params.out_dir}/${cohort}/figures/02.regional/${tierOf(peak_id)}" }, mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(sumstat),
          path(ld_cohort), path(ld_tommo), path(ld_eas), path(exons), path(recomb),
          path(lead_annotation)

    output:
    path("regional.${peak_id}.png"), emit: png
    tuple val(cohort), path("regional.${peak_id}.stats.json"), emit: stats

    script:
    """
    source activate ${params.conda_env}
    mkdir -p ld && cp ${ld_cohort} ${ld_tommo} ${ld_eas} ld/
    python3 ${params.shared_script_dir}/plot_regional.py \
        --sumstat ${sumstat} --ld-dir ld \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --exons ${exons} --recomb ${recomb} --alpha ${params.PGenomeWide} \
        --alpha-suggestive ${params.PSuggestive} \
        ${panelLabelArgs()} \
        --lead-annotation ${lead_annotation} \
        --out-png regional.${peak_id}.png
    """
}

process PLOT_FINEMAP {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/03.finemap/${tierOf(peak_id)}" }, mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(pip), path(cs), path(json), path(ld_cohort),
          path(lead_annotation)

    output:
    path("finemap.${peak_id}.png"), emit: png
    tuple val(cohort), path("finemap.${peak_id}.stats.json"), emit: stats

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/plot_finemap.py \
        --pip ${pip} --cs ${cs} --json ${json} --ld-cohort ${ld_cohort} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --alpha ${params.PGenomeWide} \
        --alpha-suggestive ${params.PSuggestive} \
        --lead-annotation ${lead_annotation} \
        --out-png finemap.${peak_id}.png
    """
}

process PLOT_CONDITIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag "${cohort}:${peak_id}"

    publishDir { "${params.out_dir}/${cohort}/figures/04.conditional/${tierOf(peak_id)}" }, mode: 'copy', pattern: '*.png'

    input:
    tuple val(cohort), val(peak_id), val(lead_id), path(rounds), path(signals), path(round_files),
          path(lead_annotation)

    output:
    path("conditional.${peak_id}.png"), emit: png
    tuple val(cohort), path("conditional.${peak_id}.stats.json"), emit: stats

    script:
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/plot_conditional.py \
        --cond-dir . --rounds ${rounds} \
        --cohort ${cohort} --locus-id ${peak_id} --lead-id ${lead_id} \
        --covar-label '${params.ModelLabels[params.PeakModel]}' \
        --pc-label '${pcLabelFor(params.PeakModel)}' --n-pcs ${nPcsFor(params.PeakModel)} \
        --alpha ${params.PGenomeWide} \
        --alpha-suggestive ${params.PSuggestive} \
        --tier ${tierOf(peak_id)} \
        --lead-annotation ${lead_annotation} \
        --out-png conditional.${peak_id}.png
    """
}

// One document per figure family per cohort, replacing the per-locus sidecars.
// A fan-out over N peaks produced 3N figures; a sidecar each meant 3N documents
// whose prose was identical and whose numbers were six lines apart. The prose is
// a property of the family, so it is written once and the numbers become rows.
process CATALOGUE_FIGURES {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${cohort}:${family}"

    publishDir { "${params.out_dir}/${cohort}/figures/${params.FigureDirs[family]}" }, mode: 'copy'

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
        --alpha ${params.PGenomeWide} --alpha-suggestive ${params.PSuggestive} \
        --pc-label '${params.PcLabel}' --n-pcs ${params.NPcs} \
        ${panelLabelArgs()} \
        --out README.md
    """
}

process CROSS_COHORT {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/tables", mode: 'copy'

    // All cohorts at once: the union of their peak leads is the variant list, so
    // this cannot be a per-cohort process. The .bed/.bim/.fam share a file name
    // across cohorts, hence stageAs — the parallel numbering is what keeps
    // geno_N.bed/.bim/.fam a valid plink triple.
    input:
    tuple val(cohorts), path(peaks, stageAs: 'peaks_*.tsv'),
          path(sumstats, stageAs: 'ss_*.tsv'),
          path(beds,  stageAs: 'geno_*.bed'),
          path(bims,  stageAs: 'geno_*.bim'),
          path(fams,  stageAs: 'geno_*.fam')
    // A SEPARATE input, not part of the tuple: `combine` flattens a collected
    // list into the tuple element-wise, so the paths arrived as N elements and
    // the declaration no longer matched.
    path(signals, stageAs: 'sig_*.tsv')

    output:
    path('lead_crosscohort.tsv'), emit: table

    script:
    def pk = [cohorts, peaks].transpose().collect    { c, f -> "--peaks ${c}=${f}" }.join(' ')
    def gl = [cohorts, sumstats].transpose().collect { c, f -> "--glm ${c}=${f}" }.join(' ')
    def bf = [cohorts, beds].transpose().collect     { c, f -> "--bfile ${c}=${f.baseName}" }.join(' ')
    """
    export PATH=${params.tool_dir}:\$PATH
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/cross_cohort.py \
        ${pk} ${gl} ${bf} \
        --cohort-order ${params.Cohorts.join(',')} \
        --p-genomewide ${params.PGenomeWide} --p-suggestive ${params.PSuggestive} \
        --peak-flank ${params.PeakFlank} \
        --signals ${signals} \
        --plink2 ${params.plink2} --tabix ${params.tabix} \
        --rsid-vcf ${params.rsid_vcf} \
        --rscript ${params.rscript} --gene-script ${params.shared_script_dir}/gene_annotate.R \
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
    python3 ${params.shared_script_dir}/plot_cohort_compare.py \
        ${q} ${p} ${a} ${h} \
        --crosscohort ${crosscohort} \
        --cohort-order ${params.Cohorts.join(',')} \
        --covar-label '${params.ModelLabels[params.PeakModel]}' --n-pcs ${nPcsFor(params.PeakModel)} \
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

    input:
    tuple val(cohorts), path(sumstats, stageAs: 'ss_*.tsv'),
          path(qc,  stageAs: 'qc_*.tsv'),
          path(mpk, stageAs: 'mpk_*.tsv')

    output:
    path('cohort_manhattan.*')

    script:
    def g = [cohorts, sumstats].transpose().collect { c, f -> "--glm ${c}=${f}" }.join(' ')
    def q = [cohorts, qc].transpose().collect       { c, f -> "--scan-qc ${c}=${f}" }.join(' ')
    def h = [cohorts, mpk].transpose().collect      { c, f -> "--model-peaks ${c}=${f}" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.shared_script_dir}/plot_cohort_manhattan.py \
        ${g} ${q} ${h} \
        --cohort-order ${params.Cohorts.join(',')} \
        --model ${params.PeakModel} \
        --alpha ${params.PGenomeWide} --suggestive ${params.PSuggestive} \
        --anno-style ${params.AnnoStyle} --repel-force ${params.RepelForce} \
        --covar-label '${params.ModelLabels[params.PeakModel]}' \
        --pc-label '${pcLabelFor(params.PeakModel)}' --n-pcs ${nPcsFor(params.PeakModel)} \
        --out-png cohort_manhattan.png
    """
}

process COMPARE_MODELS {
    executor 'slurm'
    queue 'gr10478b'
    time '2h'
    tag 'all'

    publishDir "${params.out_dir}/_comparison/figures", mode: 'copy', pattern: 'model_compare.*'
    publishDir "${params.out_dir}/_comparison/tables",  mode: 'copy', pattern: '*.tsv'

    // Every cohort x model at once: the figure's whole content is the comparison
    // between the reported model and the calibration probe, in each cohort.
    input:
    tuple val(keys), path(sumstats, stageAs: 'ss_*.tsv'),
          path(qc,  stageAs: 'qc_*.tsv'),
          path(crosscohort)

    output:
    path('model_compare.*')
    path('lead_bymodel.tsv')

    script:
    def s = [keys, sumstats].transpose().collect { k, f -> "--sumstat ${k}=${f}" }.join(' ')
    def q = [keys, qc].transpose().collect       { k, f -> "--scan-qc ${k}=${f}" }.join(' ')
    def lab = params.ModelLabels.collect { k, v -> "--model-label '${k}=${v}'" }.join(' ')
    """
    source activate ${params.conda_env}
    python3 ${params.script_dir}/plot_model_compare.py \
        ${s} ${q} ${lab} \
        --crosscohort ${crosscohort} \
        --cohort-order ${params.Cohorts.join(',')} \
        --model-order ${params.Models.join(',')} \
        --peak-model ${params.PeakModel} \
        --alpha ${params.PGenomeWide} --suggestive ${params.PSuggestive} \
        --out-png model_compare.png --out-table lead_bymodel.tsv
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
    def pcol = params.PValue == 'spa' ? 'p.value' : 'p.value.NA'
    """
    cat > run_manifest.json <<'JSON'
{
  "component": "assoc_saige",
  "engine": "SAIGE step1 fitNULLGLMM + step2 SPAtests, full GRM",
  "genotypes": "random_model (relateds retained; the GRM absorbs them)",
  "cohorts": ${groovy.json.JsonOutput.toJson(params.Cohorts)},
  "models": ${groovy.json.JsonOutput.toJson(params.Models)},
  "peak_model": "${params.PeakModel}",
  "model_covariates": ${groovy.json.JsonOutput.toJson(params.ModelCovars)},
  "model_labels": ${groovy.json.JsonOutput.toJson(params.ModelLabels)},
  "pc_label": "${params.PcLabel}",
  "n_pcs": ${params.NPcs},
  "p_value_source": "${params.PValue}",
  "p_value_column": "${pcol}",
  "p_value_note": "SAIGE emits both in one pass; this is chosen in TO_SUMSTATS, downstream of the scan, so switching it never re-runs SAIGE.",
  "firth": "is_Firth_beta=FALSE",
  "loco": true,
  "sparse_grm": false,
  "p_genomewide": ${params.PGenomeWide},
  "p_suggestive": ${params.PSuggestive},
  "peak_flank_bp": ${params.PeakFlank},
  "max_suggestive": ${params.MaxSuggestive},
  "anno_style": "${params.AnnoStyle}",
  "repel_force": ${params.RepelForce},
  "grm_markers": {"prune_window_kb": ${params.PruneWindow}, "prune_step": ${params.PruneStep}, "prune_r2": ${params.PruneR2}, "high_ld_excluded": "${params.HighLdFile}"},
  "variance_ratio_markers": ${params.VarianceRatioMarkers},
  "random_seed": ${params.RandomSeed},
  "random_seed_note": "set by step1_fitNULLGLMM.R itself before the variance-ratio markers are drawn. The GRM uses every pruned marker and is not sampled, and plink2 pruning is deterministic, so the fit is reproducible.",
  "min_grm_relationship_reported": ${params.MinRelationship},
  "susie": {"L": ${params.SusieL}, "coverage": ${params.SusieCoverage}, "min_abs_corr": ${params.SusieMinAbsCorr}},
  "model_inputs": "${params.model_inputs}",
  "saige_step1": "${params.saige_step1}",
  "saige_step2": "${params.saige_step2}",
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
    ch_bim    = ch_cohort.map { c, _bed, bim, _fam, _ph, _cv -> tuple(c, bim) }

    // ── inputs the GRM and the scan are built from ─────────────────────────
    PREP_PHENO_COV(ch_cohort)
    PRUNE_MARKERS(ch_geno)
    SPLIT_BY_CHROM(ch_geno.combine(channel.fromList(params.Autosomes)))

    // ── what the GRM is accounting for ─────────────────────────────────────
    // On exactly the markers the GRM is built from, so it describes the matrix
    // the null model fits. Independent of FIT_NULL, so it never delays it.
    RELATEDNESS(PRUNE_MARKERS.out.markers
        .map { c, bed, bim, fam, _prune -> tuple(c, bed, bim, fam) }
        .combine(ch_fam, by: 0)
        .map { c, bed, bim, gfam, fam ->
               def cmp = compareFam(c)
               tuple(c, bed, bim, gfam, fam, cmp ?: file('NO_FILE')) })

    ch_rel = RELATEDNESS.out.rel
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] },
                  o.collect { it[2] }, o.collect { it[3] }, o.collect { it[4] })
        }
        .filter { it[0] }
    COMPARE_RELATEDNESS(ch_rel)

    // ── SAIGE step 1: one null model per cohort x covariate model ──────────
    ch_null_in = PRUNE_MARKERS.out.markers
        .map { c, bed, bim, fam, _prune -> tuple(c, bed, bim, fam) }
        .combine(PREP_PHENO_COV.out.table, by: 0)
        .combine(channel.fromList(params.Models))
        .map { c, bed, bim, fam, pc, m -> tuple(c, m, bed, bim, fam, pc) }
    FIT_NULL(ch_null_in)

    // ── SAIGE step 2: one job per cohort x model x chromosome ──────────────
    // The per-chromosome genotypes are shared by BOTH models of a cohort, which
    // is why SPLIT_BY_CHROM is keyed on the cohort alone.
    ch_assoc_in = SPLIT_BY_CHROM.out.geno
        .combine(FIT_NULL.out.nullModel, by: 0)
        .map { c, chrom, bed, bim, fam, m, rda, vr -> tuple(c, m, chrom, bed, bim, fam, rda, vr) }
    RUN_ASSOC(ch_assoc_in)

    // groupTuple by the TWO key columns, not by a composite tuple key: a tuple
    // used as the `by` element never emits, and the failure is silent — the
    // channel simply stays empty and every process downstream of it is skipped.
    ch_merged = RUN_ASSOC.out.assoc.groupTuple(by: [0, 1], size: params.Autosomes.size())
    MERGE_ASSOC(ch_merged)

    // ── the adapter: SAIGE -> the canonical schema, and the P choice ───────
    TO_SUMSTATS(MERGE_ASSOC.out.merged.combine(ch_bim, by: 0)
                    .map { c, m, merged, bim -> tuple(c, m, merged, bim) })
    ch_sumstats = TO_SUMSTATS.out.sumstats                       // (cohort, model, file)
    ch_peak_ss  = ch_sumstats.filter { _c, m, _f -> m == params.PeakModel }
                             .map { c, _m, f -> tuple(c, f) }

    // ── everything below is the assoc_plink2 downstream, unchanged ─────────
    SCAN_PEAKS(ch_sumstats.map { c, m, f -> tuple(c, m, f) }
                   .groupTuple(by: 0, size: params.Models.size())
                   .combine(ch_fam, by: 0))

    ANNOTATE_LEADS(SCAN_PEAKS.out.peaks
        .combine(SCAN_PEAKS.out.model_peaks, by: 0)
        .combine(ch_geno, by: 0))

    PLOT_SCAN(ch_sumstats.combine(ANNOTATE_LEADS.out.model_peaks, by: 0))

    // ── fan out over GENOME-WIDE peaks of the peak model only ──────────────
    ch_peak = SCAN_PEAKS.out.leads
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { c, r -> tuple(c, r.peak_id, r.chrom, r.start as long, r.end as long, r.lead_id) }

    // EVERY peak is computed; only the marked ones are drawn. The split exists
    // because conditional analysis is what decides whether two nearby leads are
    // one locus or two (cross_cohort.py), and that judgement must be available
    // for every peak — otherwise it would depend on whether the peak happened to
    // be followed up, which is circular. `figure` is call_peaks.py's flag.
    ch_drawn = SCAN_PEAKS.out.leads
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { c, r ->
            // call_peaks.py is invoked by absolute path, so editing it does NOT
            // invalidate SCAN_PEAKS. A stale cached lead_variants.tsv would leave
            // `figure` null, every plot channel empty, and no error — fail here.
            if (!r.containsKey('figure')) {
                error "lead_variants.tsv has no `figure` column — SCAN_PEAKS is cached from an " +
                      "older call_peaks.py. Remove its work dirs to force it."
            }
            tuple(c, r.peak_id, (r.figure as Integer))
        }
        .filter { _c, _p, f -> f == 1 }
        .map { c, pid, _f -> tuple(c, pid) }

    // Conditional analysis re-uses the fitted null model and its variance ratio,
    // so no GRM is ever re-fitted for a conditioning round.
    ch_cond_null = FIT_NULL.out.nullModel
        .filter { _c, m, _rda, _vr -> m == params.PeakModel }
        .map { c, _m, rda, vr -> tuple(c, rda, vr) }
    ch_chrom_geno = SPLIT_BY_CHROM.out.geno
        .map { c, chrom, bed, bim, fam -> tuple(c, chrom as String, bed, bim, fam) }
    CONDITIONAL(ch_peak
        .map { c, pid, chrom, s, e, lead -> tuple(c, (chrom as String).replaceAll('^chr', ''), pid, chrom, s, e, lead) }
        .combine(ch_chrom_geno, by: [0, 1])
        .map { c, _k, pid, chrom, s, e, lead, bed, bim, fam ->
               tuple(c, pid, chrom, s, e, lead, bed, bim, fam) }
        .combine(ch_cond_null, by: 0)
        .map { c, pid, chrom, s, e, lead, bed, bim, fam, rda, vr ->
               tuple(c, pid, chrom, s, e, lead, bed, bim, fam, rda, vr) })

    ch_sumstat = SCAN_PEAKS.out.sumstats
        .flatMap { c, fs -> (fs instanceof List ? fs : [fs])
                            .collect { tuple(c, it.name.replaceAll(/\.sumstat\.tsv$/, ''), it) } }

    LD_SOURCES(ch_peak.combine(ch_sumstat, by: [0, 1]).combine(ch_geno, by: 0))

    // SuSiE's `n` is the size of the sample the summary statistics come from.
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

    // The lead annotation is per COHORT while these channels are per (cohort, peak),
    // so `by: 0` is the join. combine APPENDS, which is why `path(lead_annotation)`
    // is last in each input tuple. The locus figures now depend on ANNOTATE_LEADS,
    // which already runs and takes about a minute.
    ch_ann = ANNOTATE_LEADS.out.annotation

    PLOT_REGIONAL(LD_SOURCES.out.ld
        .map { c, pid, lead, ss, _mat, _vars, ldc, ldt, lde, _cov, ex, rec ->
               tuple(c, pid, lead, ss, ldc, ldt, lde, ex, rec) }
        .combine(ch_drawn, by: [0, 1])
        .combine(ch_ann, by: 0))

    ch_ldc = LD_SOURCES.out.ld
        .map { c, pid, _lead, _ss, _mat, _vars, ldc, _ldt, _lde, _cov, _ex, _rec ->
               tuple(c, pid, ldc) }
    PLOT_FINEMAP(SUSIE.out.susie.combine(ch_ldc, by: [0, 1]).combine(ch_drawn, by: [0, 1])
        .combine(ch_ann, by: 0))

    PLOT_CONDITIONAL(CONDITIONAL.out.cond
        .combine(ch_peak.map { c, pid, _ch, _s, _e, lead -> tuple(c, pid, lead) }, by: [0, 1])
        .combine(ch_drawn, by: [0, 1])
        .map { c, pid, rounds, signals, rfiles, lead ->
               tuple(c, pid, lead, rounds, signals, rfiles) }
        .combine(ch_ann, by: 0))

    // ── one document per figure family per cohort ─────────────────────────
    // groupTuple(by: [0, 1]) with two PLAIN STRING keys. A tuple used as the key
    // matches nothing and the channel silently never emits, which is how this
    // pattern fails: no error, just an empty downstream.
    CATALOGUE_FIGURES(
        PLOT_REGIONAL.out.stats.map { c, s -> tuple(c, 'regional', s) }
            .mix(PLOT_FINEMAP.out.stats.map { c, s -> tuple(c, 'finemap', s) })
            .mix(PLOT_CONDITIONAL.out.stats.map { c, s -> tuple(c, 'conditional', s) })
            .groupTuple(by: [0, 1]))

    // ── cross-cohort, on the peak model ───────────────────────────────────
    // Every peak's independent-signal list, from EVERY cohort. This is what makes
    // the locus definition data-driven: two leads within a peak half-width are
    // one locus unless some cohort's conditional analysis lists both as distinct
    // signals. A barrier by necessity — no locus can be named until every peak
    // has been conditioned.
    ch_signals = CONDITIONAL.out.cond.map { _c, _pid, _r, sig, _rf -> sig }.collect()

    ch_cross = SCAN_PEAKS.out.peaks
        .combine(ch_peak_ss, by: 0)
        .combine(ch_geno, by: 0)
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] }, o.collect { it[2] },
                  o.collect { it[3] }, o.collect { it[4] }, o.collect { it[5] })
        }
        .filter { it[0] }
    CROSS_COHORT(ch_cross, ch_signals)

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
        .filter { it[0] }
    COMPARE_COHORTS(ch_all.combine(CROSS_COHORT.out.table))

    ch_man = ch_peak_ss
        .combine(SCAN_PEAKS.out.qc, by: 0)
        .combine(ANNOTATE_LEADS.out.model_peaks, by: 0)
        .toList()
        .map { rows ->
            def o = params.Cohorts.collect { c -> rows.find { it[0] == c } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] },
                  o.collect { it[2] }, o.collect { it[3] })
        }
        .filter { it[0] }
    COMPARE_MANHATTAN(ch_man)

    // ── the reported model against its calibration probe ───────────────────
    // Keyed on "cohort/model" so one channel carries every cell of the grid.
    // The variants are cohort_compare's, read from CROSS_COHORT's `panel` column,
    // so the two figures are about the same loci. CROSS_COHORT already runs
    // before this stage (COMPARE_COHORTS consumes the same table), so this adds
    // no new dependency edge.
    ch_models = ch_sumstats
        .combine(SCAN_PEAKS.out.qc, by: 0)
        .map { c, m, ss, qc -> tuple("${c}/${m}", ss, qc) }
        .toList()
        .map { rows ->
            def order = []
            params.Cohorts.each { c -> params.Models.each { m -> order << "${c}/${m}" } }
            def o = order.collect { k -> rows.find { it[0] == k } }.findAll { it }
            tuple(o.collect { it[0] }, o.collect { it[1] }, o.collect { it[2] })
        }
        .filter { it[0] }
        .combine(CROSS_COHORT.out.table)
    COMPARE_MODELS(ch_models)

    WRITE_RUN_MANIFEST(channel.of('go'))
}
