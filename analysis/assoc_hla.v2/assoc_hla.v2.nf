nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════════════════════
// assoc_hla.v2 — HLA allele and amino-acid residue association on the
// hla.typing.v2 matrices: four case definitions -- all cases, cases with DVT
// and cases without DVT, each against all controls, and a case-only contrast
// of cases with DVT against cases without -- three nested cohorts, two marker
// classes, a fixed-effect and a random-effect logistic model, both adjusted
// for sex and 10 BBJ-mainland PCs.
//
// Five design points, in the order they matter:
//
//   1. THE MARKERS ARE hla.typing.v2's, UNCHANGED. Each cohort's allele and
//      residue matrices (06.assoc_inputs; 143/144/142 allele and 672/681/679
//      residue columns) are read as exported -- md5 checked against the
//      export's summary.json -- and written as hard calls (blank = missing)
//      to one PLINK set per case set x cohort x class, holding exactly that
//      case set's samples. Both models read that set. No marker is added,
//      dropped or re-filtered, so every family is fixed before any P exists.
//   2. THE TWO MODELS ARE THE GENOME-WIDE SCAN'S TWO MODELS.
//        fixed  : plink2 --glm logistic, no-firth, on the fixed_model sample
//                 set (relatives removed), complete case per marker.
//        random : SAIGE step 2, LOCO, missing calls mean-imputed; primary P
//                 the unadjusted score P (p.value.NA), as in the genome-wide
//                 SAIGE scan, with the saddle-point P kept beside it.
//      Both carry SEX + PC1-PC10 (bbj_mainland).
//   3. THE NULL MODEL FOLLOWS THE CASE SET. For 'all' it is the genome-wide
//      scan's own (analysis/assoc_saige, 'main'), REUSED. For a DVT set the
//      samples or the outcome change, so FIT_NULL fits one here with the genome-wide scan's
//      settings and GRM markers (its published grm.prune.in, on the same
//      random_model genotypes). CHECK_NULL proves every null model -- reused
//      or fitted -- has exactly the case set's samples, labels and covariates
//      before any random-model test runs.
//   4. THRESHOLDS ARE FIXED BEFORE ANY P, per (case set, cohort, class), at a
//      family-wise 0.05 over the class's tests (a two-residue mirror pair is
//      one test): LD-adjusted 0.05 / Me with Me from GEC v0.2 (PRIMARY), the
//      exact max-T threshold from the covariate-adjusted genotype correlation
//      (VALIDATION), and conventional Bonferroni 0.05 / tests (SENSITIVITY).
//      Both models use them; no BH (LD makes BH calls depend on how many
//      correlated markers move together). Not corrected across case sets,
//      cohorts or models. lambda_GC is reported, never applied.
//   5. CONSISTENCY, NOT REPLICATION. The random model is primary; a marker is
//      reported as significant in n/3 cohorts (same direction), per threshold
//      and per class (1/3 = "Single cohort"). Cohorts are nested; models share
//      samples. The two DVT sets share their controls, so their side-by-side
//      estimates are not a test of difference; 'dvt_case_only' (PHENO1
//      recoded within cases: 1 = DVT yes, 0 = DVT no) is that test.
//
// Per cohort:   CASE_SETS + NULL_MODEL_INFO('all') → VALIDATE_INPUTS
// Per case set x cohort:
//
//   PREP ─┬────────────────────────────────→ BUILD_MARKERS (x2 classes) ─┬→ EFFECTIVE_TESTS ─┐
//         │                                                               ├→ ASSOC_FIXED ─────┤
//         └→ [FIT_NULL, DVT sets only] → NULL_MODEL_INFO → CHECK_NULL ───┴→ ASSOC_RANDOM ────┴→ COLLECT
//
// ...then per case set: COMPARE; across case sets: CASE_SET_COMPARE; and,
// once everything below is done, FIGURES (one task per figure family).
//
// Conditional analyses, on a pooled PLINK set per case set x cohort of
// alleles + residues + the genome-wide random-model scan's MHC lead(s):
//
//   GWAS_LEADS → EXTRACT_LEAD (per cohort) → BUILD_POOLED ─→ STEPWISE (per case set x model x class,
//                                                             the three cohorts jointly)
//                                                          └→ CONDITIONAL (per case set x cohort x model)
//                                                             → CONDITIONAL_COMPARE
//
// independent signals within each class (entry at the GEC threshold in any
// cohort; lead = the smallest largest-P over the cohorts), HLA markers given
// the GWAS lead, and the GWAS lead given HLA markers and the signal leads --
// with the same engines as the scans.
//
// Full methodology: docs/METHODS.md. Output dictionary: docs/OUTPUTS.md.
// Every check that gates the run: docs/VERIFICATION.md.
// ═══════════════════════════════════════════════════════════════════════════

// ── Design ─────────────────────────────────────────────────────────────────
params.Cohorts      = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.ValidCohorts = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
params.Classes      = ['allele', 'residue']
params.CaseSets     = ['all', 'dvt_yes', 'dvt_no', 'dvt_case_only']

// ── Paths ──────────────────────────────────────────────────────────────────
params.project_dir  = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
params.script_dir   = "${projectDir}/scripts"
params.out_dir      = "${projectDir}/results"

// hla.typing.v2 step-06 export: the matrices, their marker table and md5s.
params.typing_dir   = "${params.project_dir}/hla.typing.v2/results/06.assoc_inputs"
params.MatrixTpl    = "${params.typing_dir}/@@COHORT@@/@@CLASS@@_dosage.tsv"
// ...and its residue matching record: which reference allele each 2-field
// allele was read from, and that allele's residues (carrier alleles of a signal).
params.AlleleMap    = "${params.project_dir}/hla.typing.v2/results/05.exports/allele_match_map.tsv"
params.ResidueRef   = "${params.project_dir}/hla.typing.v2/results/05.exports/residue_reference"

// Phenotype, covariates and genotypes: the genome-wide scans' model inputs.
params.model_inputs = "${params.project_dir}/wgs.auto.par/results/12_model_inputs"
params.PhenoTpl     = "${params.model_inputs}/@@COHORT@@/phenotype/pheno.tsv"
params.CovarTpl     = "${params.model_inputs}/@@COHORT@@/covariates/bbj_mainland_pc.sex.tsv"
params.FixedFamDir  = "${params.model_inputs}/@@COHORT@@/genotype/fixed_model"
params.FixedFamGlob = '*.maf_ge_threshold.fam'
params.RandomGenoDir  = "${params.model_inputs}/@@COHORT@@/genotype/random_model"
params.RandomGenoGlob = '*.random_model.bed'

// DVT status of the cases: the v6 sample sheet.
params.v6_info      = "${params.project_dir}/info/cteph_agp3k.v6.20260507.xlsx"
params.InfoIdCol    = 'ID_JHRPv6'
params.DvtCol       = 'DVT_Presence'

// The genome-wide SAIGE scan: its null model ('main', reused for 'all') and the
// LD-pruned GRM marker list its null models were fitted on (reused for the DVT sets).
params.NullTpl      = "${params.project_dir}/analysis/assoc_saige/results/@@COHORT@@/01.assoc/main/null"
params.PruneInTpl   = "${params.project_dir}/analysis/assoc_saige/results/@@COHORT@@/00.prep/grm.prune.in"

// The genome-wide random-model scan's MHC signals to relate to the HLA markers:
// every peak overlapping the MHC window whose lead has P < GwasLeadP in every
// cohort (the three-cohort rule), read from the scan's own tables.
params.GwasPeaksTpl    = "${params.project_dir}/analysis/assoc_saige/results/@@COHORT@@/03.peaks/peaks.tsv"
params.GwasSumstatsTpl = "${params.project_dir}/analysis/assoc_saige/results/@@COHORT@@/01.assoc/main/sumstats.tsv"
params.GwasLeadP    = 1e-5                     // the scan's suggestive line
params.MhcStart     = 28510120                 // GRCh38 extended MHC (chr6)
params.MhcEnd       = 33480577

// Gene spans (GRCh38) for marker placement and gene order only.
params.refgene      = '/LARGE0/gr10478/bin/annovar/humandb/hg38_refGene.txt'

// ── Tools (absolute: compute_env/bin precedes activated envs on PATH) ──────
params.python        = '/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python3'
params.plink2        = '/home/b/b37974/plink2'
params.saige_rscript = '/home/b/b37974/anaconda3/envs/saige/bin/Rscript'
params.saige_step1   = '/home/b/b37974/anaconda3/envs/saige/bin/step1_fitNULLGLMM.R'
params.saige_step2   = '/home/b/b37974/anaconda3/envs/saige/bin/step2_SPAtests.R'
params.java          = '/usr/bin/java'                       // Java 8; GEC v0.2
params.gec_jar       = "${projectDir}/tools/gec-0.2/gec.jar"  // archived; see tools/gec-0.2/SOURCE.md
params.GecSha256     = '9d6e60013fb6e6512a7da33ae99faae8daa363fc1d0491b479aaf0b1a109e456'

// ── Model ──────────────────────────────────────────────────────────────────
params.PhenoName    = 'PHENO1'                 // 0 = control, 1 = case
params.NPcs         = 10
params.PcLabel      = 'bbj_mainland'
params.CovarName    = 'SEX,PC1_AVG-PC10_AVG'   // plink2 range form
params.FirthMode    = 'no-firth'               // as in the genome-wide fixed scan
params.ImputeMethod = 'mean'                   // SAIGE; the blanks are undetermined copies
params.MaxMissing   = 0.15                     // SAIGE default; every marker is >= 95% called
params.MinMAC       = 0.5                      // SAIGE default; every marker has MAF >= 0.5%
params.VarianceRatioMarkers = 200              // FIT_NULL: as the genome-wide null models

// ── Conditional analyses ───────────────────────────────────────────────────
params.CondMaxSignals = 5                      // stepwise: at most this many signals per class
params.CondNTag     = 10                       // lead | HLA: top markers by |r| and by marginal P

// ── Significance ───────────────────────────────────────────────────────────
params.Alpha        = 0.05                     // family-wise error rate of all three thresholds
params.MaxTDraws    = 400000                   // max-T: multivariate-normal draws
params.MaxTSeed     = 20261023
params.MaxTBootDraws = 10000                   // benchmark: parametric-bootstrap replicates (fixed-model samples)
params.MaxTBootBatch = 500                     //   phenotype columns per plink2 run

// ── Figure families, in reading order ──────────────────────────────────────
// 01-05 main figures, 06-08 supplementary; one README.md per family.
params.Families = [['landscape', '01.landscape'], ['signals', '02.signals'], ['conditioning', '03.conditioning'],
                   ['gwas_link', '04.gwas_link'], ['dvt', '05.dvt'], ['thresholds', '06.thresholds'],
                   ['models', '07.models'], ['qq', '08.qq']]

// ═══════════════════════════════════════════════════════════════════════════
// helpers
// ═══════════════════════════════════════════════════════════════════════════

def tpl(String t, String cohort)  { t.replace('@@COHORT@@', cohort) }
def script_file(String name)      { file("${params.script_dir}/${name}", checkIfExists: true) }

// The expanded covariate list, as SAIGE stores it and the checks compare.
def covarCols() {
    return ['SEX'] + (1..params.NPcs).collect { i -> "PC${i}_AVG".toString() }
}

def oneFile(String dir, String glob) {
    def hits = files("${dir}/${glob}")
    if (hits.size() != 1) error "expected one ${glob} under ${dir}, found ${hits.size()}"
    return hits[0]
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 0 — inputs
// ═══════════════════════════════════════════════════════════════════════════

/* GENE_COORDS — GRCh38 span of each typed gene, from RefSeq. Orders the genes
 * and gives each marker a unique chromosome-6 position; no P depends on it. */
process GENE_COORDS {
    executor 'local'
    tag 'coords'

    publishDir "${params.out_dir}/_reference", mode: 'copy'

    input:
    path refgene
    val genes
    path script

    output:
    path 'hla_gene_coords.tsv', emit: coords

    script:
    """
    ${params.python} ${script} --refgene ${refgene} --genes ${genes} --out hla_gene_coords.tsv
    """
}

/* CASE_SETS — DVT status (yes / no / unknown) of every sample, from the sheet. */
process CASE_SETS {
    executor 'local'
    tag 'dvt'

    publishDir "${params.out_dir}/_reference", mode: 'copy'

    input:
    path info
    path script

    output:
    path 'case_sets.tsv', emit: tsv

    script:
    """
    ${params.python} ${script} --info ${info} --id-col ${params.InfoIdCol} \
        --dvt-col ${params.DvtCol} --out case_sets.tsv
    """
}

/* NULL_MODEL_INFO — what a SAIGE null model was fitted on (reused or fitted). */
process NULL_MODEL_INFO {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${case_set}:${cohort}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), path(rda), path(vr)
    path script

    output:
    tuple val(case_set), val(cohort), path(rda), path(vr),
          path('null_model.samples.tsv'), path('null_model.json'), emit: info

    script:
    """
    ${params.saige_rscript} ${script} --rda ${rda} --variance-ratio ${vr} --out-prefix null_model
    """
}

/* VALIDATE_INPUTS — per cohort, before any test: matrices = the typing export
 * (md5), rows = phenotype samples in order, covariates complete, the reused
 * null model = these samples, labels and covariates, DVT status covers every case. */
process VALIDATE_INPUTS {
    executor 'local'
    tag "${cohort}"

    publishDir { "${params.out_dir}/all/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(cohort), path(allele), path(residue), path(pheno), path(covar), path(fam),
          path(null_samples), path(null_json)
    path summary
    path markers
    path case_sets
    path script

    output:
    tuple val(cohort), path('validation.json'), emit: ok

    script:
    """
    ${params.python} ${script} --cohort ${cohort} --summary ${summary} --markers ${markers} \
        --allele ${allele} --residue ${residue} \
        --pheno ${pheno} --pheno-col ${params.PhenoName} \
        --covar ${covar} --covar-cols ${covarCols().join(',')} \
        --fixed-fam ${fam} --null-samples ${null_samples} --null-json ${null_json} \
        --case-sets ${case_sets} --out validation.json
    """
}

/* PREP — the phenotype and covariates of both models for one case set. */
process PREP {
    executor 'local'
    tag "${case_set}:${cohort}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/00.prep" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), path(pheno), path(covar), path(fam), path(validation)
    path case_sets
    path script

    output:
    tuple val(case_set), val(cohort), path('pheno_covar.fixed.tsv'), path('pheno_covar.random.tsv'),
          path('fixed.keep'), path('pheno_coding.json'), emit: files

    script:
    """
    ${params.python} ${script} --cohort ${cohort} --case-set ${case_set} --case-sets ${case_sets} \
        --pheno ${pheno} --pheno-col ${params.PhenoName} \
        --covar ${covar} --covar-cols ${covarCols().join(',')} --fixed-fam ${fam} \
        --out-fixed pheno_covar.fixed.tsv --out-random pheno_covar.random.tsv \
        --out-keep fixed.keep --out-json pheno_coding.json
    """
}

/* FIT_NULL — DVT sets only: the SAIGE null model on the set's samples, with the
 * genome-wide scan's settings and GRM markers (its grm.prune.in on the same
 * random_model genotypes; high-LD regions already excluded there). */
process FIT_NULL {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${case_set}:${cohort}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/00.prep/null_model" }, mode: 'copy',
               pattern: '{null.rda,null.varianceRatio.txt,null.fit.log}'

    input:
    tuple val(case_set), val(cohort), path(pheno_random), path(bed), path(bim), path(fam), path(prune_in)

    output:
    tuple val(case_set), val(cohort), path('null.rda'), path('null.varianceRatio.txt'), emit: nullModel
    path 'null.fit.log'

    script:
    """
    set -euo pipefail
    tail -n +2 ${pheno_random} | cut -f1,2 > samples.keep
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} \
        --extract ${prune_in} --keep samples.keep --make-bed --out grm_markers \
        --threads ${task.ext.threads ?: 1}
    sed '1s/^#FID/FID/' ${pheno_random} > pheno_covar.null.tsv
    ${params.saige_rscript} ${params.saige_step1} \
        --plinkFile=grm_markers \
        --phenoFile=pheno_covar.null.tsv \
        --phenoCol=${params.PhenoName} \
        --covarColList=${covarCols().join(',')} \
        --sampleIDColinphenoFile=IID \
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

/* CHECK_NULL — every null model (reused or fitted) against its case set. */
process CHECK_NULL {
    executor 'local'
    tag "${case_set}:${cohort}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/00.prep" }, mode: 'copy',
               pattern: 'null_check.json'

    input:
    tuple val(case_set), val(cohort), path(rda), path(vr), path(null_samples), path(null_json),
          path(pheno_random)
    path script

    output:
    tuple val(case_set), val(cohort), path(rda), path(vr), emit: ok
    path 'null_check.json'

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --pheno ${pheno_random} \
        --pheno-col ${params.PhenoName} --covar-cols ${covarCols().join(',')} \
        --null-samples ${null_samples} --null-json ${null_json} --out null_check.json
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 1 — markers
// ═══════════════════════════════════════════════════════════════════════════

/* BUILD_MARKERS — the case set's rows of the matrix -> VCF of hard calls
 * (REF A = absent, ALT P = present) -> one PLINK set both engines read; plus
 * per-marker counts per sample set. */
process BUILD_MARKERS {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${case_set}:${cohort}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/01.markers/${cls}" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), val(cls), path(matrix), path(pheno_random), path(keep)
    path coords
    path script

    output:
    tuple val(case_set), val(cohort), val(cls), path("${cls}.bed"), path("${cls}.bim"), path("${cls}.fam"), emit: plink
    tuple val(case_set), val(cohort), val(cls), path('marker_info.tsv'), emit: info
    path "${cls}.log"

    script:
    """
    set -euo pipefail
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --class ${cls} \
        --matrix ${matrix} --coords ${coords} --pheno ${pheno_random} --pheno-col ${params.PhenoName} \
        --fixed-keep ${keep} --out-vcf markers.vcf --out-info marker_info.tsv
    ${params.plink2} --vcf markers.vcf --double-id --make-bed \
        --threads ${task.ext.threads ?: 1} --out ${cls}
    rm -f markers.vcf
    """
}

/* EFFECTIVE_TESTS — the class's tests (mirror pairs merged) and its three
 * family-wise thresholds, from the genotypes alone: GEC v0.2 on the case set's
 * unrelated (fixed-model) samples, max-T on its random-model samples. */
process EFFECTIVE_TESTS {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${case_set}:${cohort}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/01.markers/${cls}" }, mode: 'copy',
               pattern: '{effective_tests,tests,gec}.*'

    input:
    tuple val(case_set), val(cohort), val(cls), path(bed), path(bim), path(fam), path(pheno_random), path(keep)
    path gec_jar
    path script
    path vocab

    output:
    tuple val(case_set), val(cohort), val(cls), path("effective_tests.${cls}.tsv"), path("tests.${cls}.tsv"),
          emit: tables
    path "gec.${cls}.*"

    script:
    def t = task.ext.threads ?: 1
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${t} OPENBLAS_NUM_THREADS=${t} MKL_NUM_THREADS=${t}
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --class ${cls} \
        --bfile ${bed.baseName} --pheno-random ${pheno_random} --keep ${keep} \
        --covar-cols ${covarCols().join(',')} --plink2 ${params.plink2} \
        --java ${params.java} --gec-jar ${gec_jar} --gec-sha256 ${params.GecSha256} \
        --alpha ${params.Alpha} --draws ${params.MaxTDraws} --seed ${params.MaxTSeed} --threads ${t} \
        --out-thresholds effective_tests.${cls}.tsv --out-tests tests.${cls}.tsv
    mv gec.log gec.${cls}.log
    mv gec.block.txt gec.${cls}.block.txt
    """
}

/* MAXT_BOOT — BENCHMARK, no call depends on it: the exact family-wise threshold
 * by parametric bootstrap on the case set's fixed-model samples (null
 * phenotypes drawn from the fitted null logistic model), for two statistics on
 * the same replicates: plink2's Wald test with ASSOC_FIXED's settings, and the
 * score test without SPA (the random model's statistic without the GRM); and
 * the family-wise error rate the GEC and Bonferroni thresholds actually give. */
process MAXT_BOOT {
    executor 'slurm'
    queue 'eb'          // the general queue (1 h limit): about 20 min per residue set, 5 min per allele set
    time '1h'
    tag "${case_set}:${cohort}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/01.markers/${cls}" }, mode: 'copy',
               pattern: 'maxt_boot.*'

    input:
    tuple val(case_set), val(cohort), val(cls), path(bed), path(bim), path(fam), path(pheno_fixed), path(keep),
          path(thresholds), path(tests)
    path script
    path vocab

    output:
    tuple val(case_set), val(cohort), val(cls), path("maxt_boot.${cls}.tsv"), emit: res
    path "maxt_boot.${cls}.minp.tsv.gz"
    path "maxt_boot.${cls}.observed.tsv"

    script:
    def t = task.ext.threads ?: 1
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${t} OPENBLAS_NUM_THREADS=${t} MKL_NUM_THREADS=${t}
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --class ${cls} \
        --bfile ${bed.baseName} --pheno-fixed ${pheno_fixed} --keep ${keep} \
        --tests ${tests} --thresholds ${thresholds} --pheno-col ${params.PhenoName} \
        --covar-cols ${covarCols().join(',')} --covar-name ${params.CovarName} --firth-mode ${params.FirthMode} \
        --plink2 ${params.plink2} --alpha ${params.Alpha} \
        --draws ${params.MaxTBootDraws} --batch ${params.MaxTBootBatch} --seed ${params.MaxTSeed} --threads ${t} \
        --out maxt_boot.${cls}.tsv --out-minp maxt_boot.${cls}.minp.tsv.gz \
        --out-observed maxt_boot.${cls}.observed.tsv
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 2 — association
// ═══════════════════════════════════════════════════════════════════════════

/* ASSOC_FIXED — plink2 logistic regression, the genome-wide fixed scan's
 * settings. --1: PHENO1 is 0/1 and plink2 otherwise reads 0 as missing. */
process ASSOC_FIXED {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${case_set}:${cohort}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/02.assoc/fixed/${cls}" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), val(cls), path(bed), path(bim), path(fam), path(pheno_fixed), path(keep)

    output:
    tuple val(case_set), val(cohort), val(cls), val('fixed'), path("fixed.${cls}.PHENO1.glm.logistic"),
          path("fixed.${cls}.log"), emit: res

    script:
    """
    ${params.plink2} --bfile ${bed.baseName} --keep ${keep} \
        --pheno ${pheno_fixed} --pheno-name ${params.PhenoName} --1 \
        --covar ${pheno_fixed} --covar-name ${params.CovarName} \
        --glm omit-ref ${params.FirthMode} hide-covar cols=+a1countcc,+totallelecc \
        --ci 0.95 \
        --threads ${task.ext.threads ?: 1} \
        --out fixed.${cls}
    """
}

/* ASSOC_RANDOM — SAIGE step 2 on the case set's checked null model (LOCO: the
 * chromosome-6-left-out predictor). */
process ASSOC_RANDOM {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${case_set}:${cohort}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/02.assoc/random/${cls}" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), val(cls), path(bed), path(bim), path(fam), path(rda), path(vr)

    output:
    tuple val(case_set), val(cohort), val(cls), val('random'), path("random.${cls}.saige.txt"),
          path("random.${cls}.log"), emit: res

    script:
    """
    set -o pipefail
    ${params.saige_rscript} ${params.saige_step2} \
        --bedFile=${bed} --bimFile=${bim} --famFile=${fam} \
        --AlleleOrder=alt-first \
        --chrom=6 \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${vr} \
        --SAIGEOutputFile=random.${cls}.saige.txt \
        --is_Firth_beta=FALSE \
        --LOCO=TRUE \
        --is_output_moreDetails=TRUE \
        --minMAC=${params.MinMAC} \
        --maxMissing=${params.MaxMissing} \
        --impute_method=${params.ImputeMethod} 2>&1 | tee random.${cls}.log
    """
}

/* COLLECT — one schema, the calls at the three thresholds, a one-row scan
 * summary; stops the run if the engine analysed other samples or markers than
 * intended. */
process COLLECT {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${case_set}:${cohort}:${cls}:${model}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/03.tables" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), val(cls), val(model), path(engine_out), path(engine_log),
          path(info), path(thresholds), path(tests), path(coding)
    path script

    output:
    tuple val(case_set), path("assoc.${model}.${cls}.tsv"), emit: assoc
    tuple val(case_set), path("scan_qc.${model}.${cls}.tsv"), emit: qc

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --class ${cls} --model ${model} \
        --engine-out ${engine_out} --engine-log ${engine_log} \
        --marker-info ${info} --coding ${coding} --n-covariates ${covarCols().size()} \
        --thresholds ${thresholds} --tests ${tests} \
        --out-assoc assoc.${model}.${cls}.tsv --out-qc scan_qc.${model}.${cls}.tsv
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 3 — across cohorts and models, then across case sets
// ═══════════════════════════════════════════════════════════════════════════

process COMPARE {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${case_set}"

    publishDir { "${params.out_dir}/${case_set}/_comparison" }, mode: 'copy'

    input:
    tuple val(case_set), path(assoc, stageAs: 'assoc_??/*'), path(qc, stageAs: 'qc_??/*')
    path script

    output:
    tuple val(case_set), path('assoc_all.tsv'), path('scan_qc_all.tsv'), path('model_concordance.tsv'),
          path('marker_summary.tsv'), path('hits.tsv'), emit: tables

    script:
    """
    ${params.python} ${script} --assoc ${assoc} --qc ${qc} --out-dir .
    """
}

/* CASE_SET_COMPARE — all cases, DVT yes and DVT no side by side. */
process CASE_SET_COMPARE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'case_sets'

    publishDir "${params.out_dir}/_case_sets", mode: 'copy'

    input:
    path assoc, stageAs: 'assoc_?/*'
    path qc,    stageAs: 'qc_?/*'
    path script
    path vocab

    output:
    path '*.tsv'

    script:
    """
    ${params.python} ${script} --assoc ${assoc} --qc ${qc} --out-dir .
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 4 — conditional analyses
// ═══════════════════════════════════════════════════════════════════════════

/* GWAS_LEADS — the MHC leads of the genome-wide random-model scan that pass the
 * three-cohort rule; every MHC peak is listed with its reason. */
process GWAS_LEADS {
    executor 'local'
    tag 'gwas_leads'

    publishDir "${params.out_dir}/_conditional", mode: 'copy'

    input:
    val cohorts
    path peaks, stageAs: 'peaks_?.tsv'
    path sumstats, stageAs: 'sumstats_?.tsv'
    path script

    output:
    path 'gwas_leads.tsv', emit: table
    path 'gwas_leads.ids', emit: ids

    script:
    """
    ${params.python} ${script} --cohorts ${cohorts} --peaks ${peaks} --sumstats ${sumstats} \
        --start ${params.MhcStart} --end ${params.MhcEnd} --p ${params.GwasLeadP} \
        --out gwas_leads.tsv --out-ids gwas_leads.ids
    """
}

/* EXTRACT_LEAD — the lead variant(s) from the cohort's random_model genotypes
 * (the genotypes the genome-wide scan tested), as a VCF named by IID. */
process EXTRACT_LEAD {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${cohort}"

    input:
    tuple val(cohort), path(bed), path(bim), path(fam), path(_prune)
    path ids

    output:
    tuple val(cohort), path('lead.vcf'), emit: vcf

    script:
    """
    ${params.plink2} --bed ${bed} --bim ${bim} --fam ${fam} --extract ${ids} \
        --export vcf id-paste=iid --threads ${task.ext.threads ?: 1} --out lead
    n=\$(grep -vc '^#' lead.vcf || true)
    [ "\$n" -eq "\$(wc -l < ${ids})" ] || { echo "lead.vcf has \$n variant(s), expected \$(wc -l < ${ids})" >&2; exit 1; }
    """
}

/* BUILD_POOLED — alleles + residues + lead(s) in one PLINK set (sentinel first;
 * see build_pooled.py), and the LD of every HLA marker with each lead. */
process BUILD_POOLED {
    executor 'slurm'
    queue 'gr10478b'
    time '30m'
    tag "${case_set}:${cohort}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/04.conditional/pooled" }, mode: 'copy'

    input:
    tuple val(case_set), val(cohort), path(allele), path(residue), path(pheno_random), path(lead_vcf)
    path coords
    path script

    output:
    tuple val(case_set), val(cohort), path('pooled.bed'), path('pooled.bim'), path('pooled.fam'),
          path('pooled_markers.tsv'), emit: pooled
    path 'pooled.log'

    script:
    """
    set -euo pipefail
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --allele ${allele} \
        --residue ${residue} --lead-vcf ${lead_vcf} --coords ${coords} --pheno ${pheno_random} \
        --pheno-col ${params.PhenoName} --out-vcf pooled.vcf --out-markers pooled_markers.tsv
    ${params.plink2} --vcf pooled.vcf --double-id --make-bed \
        --threads ${task.ext.threads ?: 1} --out pooled
    rm -f pooled.vcf
    """
}

/* STEPWISE — independent signals of one class for one (case set, model),
 * selected jointly over the three cohorts (stepwise.py), with the scan's own
 * engine on each cohort's pooled set. Every per-cohort input is a list in the
 * order of `cohorts`, staged as c1/, c2/, c3/. */
process STEPWISE {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${case_set}:${model}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/_conditional/${model}" }, mode: 'copy',
               pattern: '{steps,stepwise,signals,leads}.*'

    input:
    tuple val(case_set), val(model), val(cls), val(cohorts),
          path(bed, stageAs: 'c?/*'), path(bim, stageAs: 'c?/*'), path(fam, stageAs: 'c?/*'),
          path(markers, stageAs: 'c?/*'), path(pheno_fixed, stageAs: 'c?/*'), path(keep, stageAs: 'c?/*'),
          path(coding, stageAs: 'c?/*'), path(rda, stageAs: 'c?/*'), path(vr, stageAs: 'c?/*'),
          path(thresholds, stageAs: 'c?/*')
    path allele_map
    path residue_ref
    path script
    path engine
    path vocab

    output:
    tuple val(case_set), val(model), val(cls), path("steps.${cls}.tsv"), path("signals.${cls}.tsv"),
          path("leads.${cls}.txt"), emit: res
    path "stepwise.${cls}.tsv", emit: long

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --model ${model} --class ${cls} \
        --cohorts ${cohorts.join(',')} --bed ${bed} --markers ${markers} --coding ${coding} \
        --pheno-fixed ${pheno_fixed} --keep ${keep} --rda ${rda} --vr ${vr} --thresholds ${thresholds} \
        --allele-map ${allele_map} --residue-ref ${residue_ref} --covar-name ${params.CovarName} \
        --plink2 ${params.plink2} --saige-rscript ${params.saige_rscript} --saige-step2 ${params.saige_step2} \
        --threads ${task.ext.threads ?: 1} --max-signals ${params.CondMaxSignals} --out-dir .
    """
}

/* CONDITIONAL — the GWAS lead(s) against the HLA markers for one (case set,
 * cohort, model): HLA | lead, lead | HLA (LD tags, top P, the signal leads of
 * both classes from STEPWISE), with the scan's own engine. */
process CONDITIONAL {
    executor 'slurm'
    queue 'gr10478b'
    time '4h'
    tag "${case_set}:${cohort}:${model}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/04.conditional/${model}" }, mode: 'copy',
               pattern: '{hla_given_lead.tsv,lead_given_hla.tsv,conditional.json}'

    input:
    tuple val(case_set), val(cohort), val(model), path(bed), path(bim), path(fam), path(markers),
          path(pheno_fixed), path(keep), path(coding), path(rda), path(vr), path(signal_leads), path(position_leads)
    path script
    path engine

    output:
    tuple val(case_set), val(cohort), val(model), path('hla_given_lead.tsv'), path('lead_given_hla.tsv'), emit: res
    path 'conditional.json'

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --model ${model} \
        --bfile ${bed.baseName} --markers ${markers} --coding ${coding} \
        --pheno-fixed ${pheno_fixed} --keep ${keep} --covar-name ${params.CovarName} \
        --rda ${rda} --vr ${vr} --signal-leads ${signal_leads} --position-leads ${position_leads} \
        --plink2 ${params.plink2} --saige-rscript ${params.saige_rscript} --saige-step2 ${params.saige_step2} \
        --threads ${task.ext.threads ?: 1} --n-tag ${params.CondNTag} --out-dir .
    """
}

/* HLA_GIVEN_SIGNAL — every HLA marker (and the GWAS lead) given each signal
 * lead of the case set and model, and given each whole position the lead is
 * conditioned on, with the scan's own engine on the pooled set
 * (hla_given_signal.py). It completes HLA | signal: STEPWISE steps only within
 * a class that has a signal, so e.g. the alleles given the residue lead come
 * from here. Same inputs as CONDITIONAL. */
process HLA_GIVEN_SIGNAL {
    executor 'slurm'
    queue 'eb'          // the general queue (1 h limit): a few engine runs per task
    time '1h'
    tag "${case_set}:${cohort}:${model}"

    publishDir { "${params.out_dir}/${case_set}/${cohort}/04.conditional/${model}" }, mode: 'copy',
               pattern: 'hla_given_signal.tsv'

    input:
    tuple val(case_set), val(cohort), val(model), path(bed), path(bim), path(fam), path(markers),
          path(pheno_fixed), path(keep), path(coding), path(rda), path(vr), path(signal_leads), path(position_leads)
    path script
    path conditional
    path engine

    output:
    tuple val(case_set), val(cohort), val(model), path('hla_given_signal.tsv'), emit: res

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --cohort ${cohort} --model ${model} \
        --bfile ${bed.baseName} --markers ${markers} --coding ${coding} \
        --pheno-fixed ${pheno_fixed} --keep ${keep} --covar-name ${params.CovarName} \
        --rda ${rda} --vr ${vr} --signal-leads ${signal_leads} --position-leads ${position_leads} \
        --plink2 ${params.plink2} --saige-rscript ${params.saige_rscript} --saige-step2 ${params.saige_step2} \
        --threads ${task.ext.threads ?: 1} --out-dir .
    """
}

/* OMNIBUS_GIVEN — the omnibus tests of every residue position given ONE marker
 * at a time: each GWAS MHC lead and each residue signal lead of the case set
 * (either model), as one more covariate (omnibus_given.py, with omnibus.py's
 * own test). The marginal scan is recomputed and must equal OMNIBUS's. */
process OMNIBUS_GIVEN {
    executor 'slurm'
    queue 'eb'          // the general queue (1 h limit): a few minutes per task
    time '1h'
    tag "${case_set}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/_omnibus" }, mode: 'copy', pattern: 'omnibus_given.*.tsv'

    input:
    tuple val(case_set), val(cls), val(cohorts),
          path(bed, stageAs: 'c?/*'), path(bim, stageAs: 'c?/*'), path(fam, stageAs: 'c?/*'),
          path(pheno_fixed, stageAs: 'c?/*'), path(keep, stageAs: 'c?/*'),
          path(pbed, stageAs: 'p?/*'), path(pbim, stageAs: 'p?/*'), path(pfam, stageAs: 'p?/*'),
          path(signal_leads, stageAs: 'sl?/*'), path(omnibus)
    path leads
    path script
    path omnibus_py
    path vocab

    output:
    tuple val(case_set), val(cls), path("omnibus_given.${cls}.tsv"), emit: res

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --class ${cls} --cohorts ${cohorts.join(',')} \
        --bed ${bed} --pheno-fixed ${pheno_fixed} --keep ${keep} --pooled ${pbed} \
        --leads ${leads} --signal-leads ${signal_leads} --omnibus ${omnibus} --pheno-col ${params.PhenoName} \
        --covar-cols ${covarCols().join(',')} --plink2 ${params.plink2} --alpha ${params.Alpha} \
        --threads ${task.ext.threads ?: 1}
    """
}

process CONDITIONAL_COMPARE {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag 'conditional'

    publishDir "${params.out_dir}/_conditional", mode: 'copy'

    input:
    path steps,    stageAs: 'steps_??/*'
    path signals,  stageAs: 'signals_??/*'
    path hla_given, stageAs: 'hla_given_??/*'
    path lead_given, stageAs: 'lead_given_??/*'
    path markers,  stageAs: 'markers_??/*'
    path gwas_leads
    path script
    path vocab

    output:
    tuple path('signals_all.tsv'), path('steps_all.tsv'), path('lead_given_hla_all.tsv'),
          path('hla_given_lead_all.tsv'), emit: tables
    path '*.tsv'

    script:
    """
    ${params.python} ${script} --steps ${steps} --signals ${signals} --hla-given ${hla_given} \
        --lead-given ${lead_given} --markers ${markers} --gwas-leads ${gwas_leads} --out-dir .
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// STAGE 5 — omnibus tests of residue positions (a main analysis; fixed model, Bonferroni)
// ═══════════════════════════════════════════════════════════════════════════

/* OMNIBUS — MAIN ANALYSIS: each residue position tested as a whole
 * (likelihood-ratio test of all its residues, fixed model on the case set's
 * unrelated samples), Bonferroni over the positions; independent positions
 * selected stepwise over the three cohorts, conditioning on every residue of
 * the entered positions (omnibus.py). Residue scans only: gene-level omnibus
 * tests over alleles are not run (METHODS §8). Per-cohort inputs are lists in
 * the order of `cohorts`, staged as c1/, c2/, c3/. */
process OMNIBUS {
    executor 'slurm'
    queue 'eb'          // the general queue (1 h limit): under a minute per task
    time '1h'
    tag "${case_set}:${cls}"

    publishDir { "${params.out_dir}/${case_set}/_omnibus" }, mode: 'copy'

    input:
    tuple val(case_set), val(cls), val(cohorts),
          path(bed, stageAs: 'c?/*'), path(bim, stageAs: 'c?/*'), path(fam, stageAs: 'c?/*'),
          path(pheno_fixed, stageAs: 'c?/*'), path(keep, stageAs: 'c?/*')
    path script
    path vocab

    output:
    tuple val(case_set), val(cls), path("omnibus.${cls}.tsv"), path("omnibus_signals.${cls}.tsv"),
          path("omnibus_steps.${cls}.tsv"), emit: res
    path "omnibus_stepwise.${cls}.tsv"

    script:
    """
    ${params.python} ${script} --case-set ${case_set} --class ${cls} --cohorts ${cohorts.join(',')} \
        --bed ${bed} --pheno-fixed ${pheno_fixed} --keep ${keep} --pheno-col ${params.PhenoName} \
        --covar-cols ${covarCols().join(',')} --plink2 ${params.plink2} --alpha ${params.Alpha} \
        --max-signals ${params.CondMaxSignals} --threads ${task.ext.threads ?: 1}
    """
}

/* OMNIBUS_SUMMARY — the omnibus tests of every case set and class stacked and
 * set beside the single-marker scans: the best single marker of each unit in
 * the fixed and the random model, and which unit holds each random-model
 * signal lead (omnibus_summary.py). */
process OMNIBUS_SUMMARY {
    executor 'local'
    tag 'omnibus'

    publishDir "${params.out_dir}/_omnibus", mode: 'copy'

    input:
    path omnibus, stageAs: 'om_??/*'
    path signals, stageAs: 'sg_??/*'
    path steps, stageAs: 'st_??/*'
    path assoc, stageAs: 'as_??/*'
    path single_signals
    val cohorts
    path script
    path vocab

    output:
    tuple path('omnibus_all.tsv'), path('omnibus_units.tsv'), path('omnibus_signals_all.tsv'),
          path('omnibus_steps_all.tsv'), emit: tables

    script:
    """
    ${params.python} ${script} --omnibus ${omnibus} --signals ${signals} --steps ${steps} \
        --assoc ${assoc} --single-signals ${single_signals} --cohorts ${cohorts}
    """
}

/* FIGURES — one figure family per task (figures.py), then its one README
 * (catalogue.py), from every case set's tables, every family's thresholds and
 * the conditional analyses. Families 01-05 are the main figures, 06-08
 * supplementary; all in results/figures/0N.<family>/. */
process FIGURES {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    tag "${family}"

    publishDir { "${params.out_dir}/figures/${dir}" }, mode: 'copy', pattern: '{*.png,README.md}'

    input:
    tuple val(family), val(dir)
    path assoc,      stageAs: 'assoc_?/*'
    path summary,    stageAs: 'summary_?/*'
    path thresholds, stageAs: 'thr_??/*'
    path stepwise,   stageAs: 'sw_??/*'
    tuple path(signals), path(steps), path(lead_given), path(hla_given)
    path gwas_leads
    path coords
    path figures_script
    path catalogue_script
    path vocab

    output:
    path '*.png', optional: true
    path '*.stats.json'
    path 'README.md', emit: readme

    script:
    def data = "--assoc ${assoc} --summary ${summary} --thresholds ${thresholds} --stepwise ${stepwise} " +
               "--signals ${signals} --steps ${steps} --lead-given ${lead_given} --hla-given ${hla_given} " +
               "--gwas-leads ${gwas_leads} --coords ${coords}"
    """
    export MPLCONFIGDIR=\$PWD/.mpl
    ${params.python} ${figures_script} --family ${family} ${data} \
        --covar-label 'sex + PC1-PC${params.NPcs} (${params.PcLabel})' --out-dir .
    ${params.python} ${catalogue_script} --family ${family} ${data} \
        --stats *.stats.json --figure-dir ${dir} --out README.md
    """
}

process WRITE_RUN_MANIFEST {
    executor 'local'
    tag 'manifest'

    publishDir "${params.out_dir}/_run_info", mode: 'copy', pattern: 'run_manifest.json'
    publishDir "${params.out_dir}",           mode: 'copy', pattern: 'README.md'

    input:
    val manifest
    val index
    val done

    output:
    path 'run_manifest.json'
    path 'README.md'

    script:
    """
    cat > run_manifest.json <<'JSON'
${manifest}
JSON
    cat > README.md <<'MD'
${index}
MD
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// Workflow
//
// Every `by:` below is a plain index or index list, never a composite tuple
// key: a tuple used as a `by:` element matches nothing and the channel simply
// never emits. verify.sh checks the per-process task counts in trace.txt
// against this DAG's cardinality, which is how an empty channel is caught.
// ═══════════════════════════════════════════════════════════════════════════

workflow {

    // ── validation, before anything is submitted ───────────────────────────
    def cohortList = (params.Cohorts instanceof List)
        ? params.Cohorts
        : params.Cohorts.toString().split(',').collect { x -> x.trim() }.findAll { x -> x }
    def bad = cohortList.findAll { ch -> !(ch in params.ValidCohorts) }
    if (bad) error "unknown cohort(s): ${bad.join(', ')}. Valid: ${params.ValidCohorts.join(', ')}"
    def badSets = params.CaseSets.findAll { s -> !(s in ['all', 'dvt_yes', 'dvt_no', 'dvt_case_only']) }
    if (badSets) error "unknown case set(s): ${badSets.join(', ')}. Valid: all, dvt_yes, dvt_no, dvt_case_only"
    if (!('all' in params.CaseSets)) error "CaseSets must include 'all': VALIDATE_INPUTS checks the reused null model"
    if (params.ImputeMethod != 'mean') {
        error "ImputeMethod is '${params.ImputeMethod}'. Only 'mean' is supported: COLLECT asserts it from\n" +
              "  the SAIGE log, and SAIGE's default best_guess would set every undetermined copy to\n" +
              "  the commonest genotype, i.e. to 'absent' for every rare marker."
    }

    def summaryFile = file("${params.typing_dir}/summary.json", checkIfExists: true)
    def markersFile = file("${params.typing_dir}/markers.tsv", checkIfExists: true)
    def typed = new groovy.json.JsonSlurper().parseText(summaryFile.text)
    def genes = typed.genes.join(',')

    log.info "assoc_hla.v2 | case sets: ${params.CaseSets.join(', ')}"
    log.info "assoc_hla.v2 | cohorts: ${cohortList.join(', ')} | classes: ${params.Classes.join(', ')}"
    log.info "assoc_hla.v2 | covariates: ${params.CovarName} (${params.PcLabel})"
    log.info "assoc_hla.v2 | fixed: plink2 --glm ${params.FirthMode} | random: SAIGE step 2, impute ${params.ImputeMethod}"

    ch_cohort = channel.fromList(cohortList)
    ch_cohort_set = ch_cohort.combine(channel.fromList(params.CaseSets))       // (cohort, case_set)

    // ── STAGE 0 ────────────────────────────────────────────────────────────
    GENE_COORDS(file(params.refgene, checkIfExists: true), genes, script_file('gene_coords.py'))
    CASE_SETS(file(params.v6_info, checkIfExists: true), script_file('case_sets.py'))

    ch_null_reused = ch_cohort.map { c ->
        def p = tpl(params.NullTpl, c)
        tuple('all', c, file("${p}.rda", checkIfExists: true), file("${p}.varianceRatio.txt", checkIfExists: true))
    }

    ch_inputs = ch_cohort.map { c ->
        tuple(c, file(tpl(params.PhenoTpl, c), checkIfExists: true),
                 file(tpl(params.CovarTpl, c), checkIfExists: true),
                 oneFile(tpl(params.FixedFamDir, c), params.FixedFamGlob))
    }
    ch_matrix = ch_cohort.flatMap { c ->
        params.Classes.collect { cls ->
            tuple(c, cls, file(tpl(params.MatrixTpl, c).replace('@@CLASS@@', cls), checkIfExists: true))
        }
    }

    // NULL_MODEL_INFO runs twice, through two thin sub-workflows: once on the
    // reused null models (VALIDATE_INPUTS needs it before PREP) and once on the
    // null models FIT_NULL fits after PREP -- one invocation cannot serve both
    // without a cycle in the DAG.
    ch_validate_base = ch_matrix.filter { _c, cls, _m -> cls == 'allele' }.map { c, _cls, m -> tuple(c, m) }
        .join(ch_matrix.filter { _c, cls, _m -> cls == 'residue' }.map { c, _cls, m -> tuple(c, m) }, by: 0)
        .join(ch_inputs, by: 0)

    ch_random_geno = ch_cohort.map { c ->
        def bed = oneFile(tpl(params.RandomGenoDir, c), params.RandomGenoGlob)
        def stem = bed.toString().replaceAll(/\.bed$/, '')
        tuple(c, bed, file("${stem}.bim", checkIfExists: true), file("${stem}.fam", checkIfExists: true),
                 file(tpl(params.PruneInTpl, c), checkIfExists: true))
    }

    NULL_INFO_ALL(ch_null_reused, script_file('null_model_info.R'))

    VALIDATE_INPUTS(
        ch_validate_base.join(NULL_INFO_ALL.out.map { _s, c, _r, _v, ns, nj -> tuple(c, ns, nj) }, by: 0),
        summaryFile, markersFile, CASE_SETS.out.tsv, script_file('validate_inputs.py'))

    ch_prep_in = ch_cohort_set
        .combine(ch_inputs.join(VALIDATE_INPUTS.out.ok, by: 0), by: 0)
        .map { c, s, ph, cv, fam, val -> tuple(s, c, ph, cv, fam, val) }
    PREP(ch_prep_in, CASE_SETS.out.tsv, script_file('prep_pheno_cov.py'))

    ch_random_pheno = PREP.out.files.map { s, c, _fx, rd, _keep, _js -> tuple(s, c, rd) }

    // FIT_NULL for the DVT sets only.
    ch_fit_in = ch_random_pheno.filter { s, _c, _rd -> s != 'all' }
        .map { s, c, rd -> tuple(c, s, rd) }
        .combine(ch_random_geno, by: 0)
        .map { c, s, rd, bed, bim, fam, prune -> tuple(s, c, rd, bed, bim, fam, prune) }
    FIT_NULL(ch_fit_in)
    NULL_INFO_FIT(FIT_NULL.out.nullModel, script_file('null_model_info.R'))

    ch_check_in = NULL_INFO_ALL.out.mix(NULL_INFO_FIT.out)
        .combine(ch_random_pheno, by: [0, 1])
    CHECK_NULL(ch_check_in, script_file('check_null.py'))

    // ── STAGE 1 ────────────────────────────────────────────────────────────
    ch_build_in = PREP.out.files.map { s, c, _fx, rd, keep, _js -> tuple(c, s, rd, keep) }
        .combine(ch_matrix, by: 0)
        .map { c, s, rd, keep, cls, m -> tuple(s, c, cls, m, rd, keep) }
    BUILD_MARKERS(ch_build_in, GENE_COORDS.out.coords, script_file('build_markers.py'))
    EFFECTIVE_TESTS(
        BUILD_MARKERS.out.plink.combine(
            PREP.out.files.map { s, c, _fx, rd, keep, _js -> tuple(s, c, rd, keep) }, by: [0, 1]),
        file(params.gec_jar, checkIfExists: true), script_file('effective_tests.py'), script_file('vocab.py'))
    // BENCHMARK (no call depends on it): exact max-T by parametric bootstrap on the fixed-model samples
    MAXT_BOOT(
        BUILD_MARKERS.out.plink
            .combine(PREP.out.files.map { s, c, fx, _rd, keep, _js -> tuple(s, c, fx, keep) }, by: [0, 1])
            .combine(EFFECTIVE_TESTS.out.tables, by: [0, 1, 2]),
        script_file('maxt_boot.py'), script_file('vocab.py'))

    // ── STAGE 2 ────────────────────────────────────────────────────────────
    ASSOC_FIXED(BUILD_MARKERS.out.plink.combine(
        PREP.out.files.map { s, c, fx, _rd, keep, _js -> tuple(s, c, fx, keep) }, by: [0, 1]))
    ASSOC_RANDOM(BUILD_MARKERS.out.plink.combine(CHECK_NULL.out.ok, by: [0, 1]))

    ch_collect_in = ASSOC_FIXED.out.res.mix(ASSOC_RANDOM.out.res)
        .combine(BUILD_MARKERS.out.info, by: [0, 1, 2])
        .combine(EFFECTIVE_TESTS.out.tables, by: [0, 1, 2])
        .combine(PREP.out.files.map { s, c, _fx, _rd, _keep, js -> tuple(s, c, js) }, by: [0, 1])
    COLLECT(ch_collect_in, script_file('collect.py'))

    // ── STAGE 3 ────────────────────────────────────────────────────────────
    ch_compare_in = COLLECT.out.assoc.groupTuple(by: 0)
        .join(COLLECT.out.qc.groupTuple(by: 0), by: 0)
    COMPARE(ch_compare_in, script_file('compare.py'))

    CASE_SET_COMPARE(
        COMPARE.out.tables.map { _s, a, _q, _c, _m, _h -> a }.collect(),
        COMPARE.out.tables.map { _s, _a, q, _c, _m, _h -> q }.collect(),
        script_file('case_set_compare.py'), script_file('vocab.py'))

    // ── STAGE 4 — conditional analyses ─────────────────────────────────────
    GWAS_LEADS(cohortList.join(','),
               cohortList.collect { c -> file(tpl(params.GwasPeaksTpl, c), checkIfExists: true) },
               cohortList.collect { c -> file(tpl(params.GwasSumstatsTpl, c), checkIfExists: true) },
               script_file('gwas_leads.py'))
    EXTRACT_LEAD(ch_random_geno, GWAS_LEADS.out.ids)

    ch_pooled_in = ch_random_pheno
        .map { s, c, rd -> tuple(c, s, rd) }
        .combine(ch_matrix.filter { _c, cls, _m -> cls == 'allele' }.map { c, _cls, m -> tuple(c, m) }, by: 0)
        .combine(ch_matrix.filter { _c, cls, _m -> cls == 'residue' }.map { c, _cls, m -> tuple(c, m) }, by: 0)
        .combine(EXTRACT_LEAD.out.vcf, by: 0)
        .map { c, s, rd, al, re, vcf -> tuple(s, c, al, re, rd, vcf) }
    BUILD_POOLED(ch_pooled_in, GENE_COORDS.out.coords, script_file('build_pooled.py'))

    // (case set, cohort, pooled set, phenotype files, null model)
    ch_cond_base = BUILD_POOLED.out.pooled
        .combine(PREP.out.files.map { s, c, fx, _rd, keep, js -> tuple(s, c, fx, keep, js) }, by: [0, 1])
        .combine(CHECK_NULL.out.ok, by: [0, 1])

    // STEPWISE: one task per (case set, model, class) over the three cohorts.
    // groupTuple's order is arbitrary, so every list is put in cohortList order
    // (stable inputs keep -resume's cache hits).
    ch_step_in = ch_cond_base
        .combine(EFFECTIVE_TESTS.out.tables.map { s, c, cls, th, _t -> tuple(s, c, cls, th) }, by: [0, 1])
        .combine(channel.of('fixed', 'random'))
        .map { s, c, bed, bim, fam, mk, fx, keep, js, rda, vr, cls, th, model ->
            tuple(s, model, cls, c, bed, bim, fam, mk, fx, keep, js, rda, vr, th) }
        .groupTuple(by: [0, 1, 2], size: cohortList.size())
        .map { row ->
            def idx = cohortList.collect { c -> row[3].indexOf(c) }
            [row[0], row[1], row[2]] + (3..<row.size()).collect { j -> idx.collect { i -> row[j][i] } } }
    STEPWISE(ch_step_in, file(params.AlleleMap, checkIfExists: true), file(params.ResidueRef, checkIfExists: true),
             script_file('stepwise.py'), script_file('engine.py'), script_file('vocab.py'))

    // ── STAGE 5 — omnibus tests of residue positions (a main analysis) ─────
    // run before CONDITIONAL, which conditions the GWAS lead on the omnibus position leads
    // one task per case set over the three cohorts, lists in cohortList order; residues only
    ch_omni_in = BUILD_MARKERS.out.plink
        .filter { _s, _c, cls, _bed, _bim, _fam -> cls == 'residue' }
        .combine(PREP.out.files.map { s, c, fx, _rd, keep, _js -> tuple(s, c, fx, keep) }, by: [0, 1])
        .map { s, c, cls, bed, bim, fam, fx, keep -> tuple(s, cls, c, bed, bim, fam, fx, keep) }
        .groupTuple(by: [0, 1], size: cohortList.size())
        .map { row ->
            def idx = cohortList.collect { c -> row[2].indexOf(c) }
            [row[0], row[1]] + (2..<row.size()).collect { j -> idx.collect { i -> row[j][i] } } }
    OMNIBUS(ch_omni_in, script_file('omnibus.py'), script_file('vocab.py'))

    // CONDITIONAL: per (case set, cohort, model), with that case set and
    // model's signal leads of both classes, and the case set's omnibus position leads.
    ch_leads = STEPWISE.out.res.map { s, m, _cls, _st, _sg, ld -> tuple(s, m, ld) }
        .groupTuple(by: [0, 1], size: params.Classes.size())
        .map { s, m, lds -> tuple(s, m, lds.sort { f -> f.name }) }
    ch_cond_in = ch_cond_base
        .combine(channel.of('fixed', 'random'))
        .map { s, c, bed, bim, fam, mk, fx, keep, js, rda, vr, model ->
            tuple(s, model, c, bed, bim, fam, mk, fx, keep, js, rda, vr) }
        .combine(ch_leads, by: [0, 1])
        .map { s, model, c, bed, bim, fam, mk, fx, keep, js, rda, vr, lds ->
            tuple(s, c, model, bed, bim, fam, mk, fx, keep, js, rda, vr, lds) }
        .combine(OMNIBUS.out.res.map { s, _cls, _om, sg, _st -> tuple(s, sg) }, by: 0)
    CONDITIONAL(ch_cond_in, script_file('conditional.py'), script_file('engine.py'))
    HLA_GIVEN_SIGNAL(ch_cond_in, script_file('hla_given_signal.py'), script_file('conditional.py'),
                     script_file('engine.py'))

    // OMNIBUS_GIVEN: per case set, the residue positions given each GWAS lead and each residue signal lead
    // (both models' leads.residue.txt, fixed then random); pooled sets in cohortList order.
    ch_pooled_set = BUILD_POOLED.out.pooled.map { s, c, bed, bim, fam, _mk -> tuple(s, c, bed, bim, fam) }
        .groupTuple(by: 0, size: cohortList.size())
        .map { row ->
            def idx = cohortList.collect { c -> row[1].indexOf(c) }
            [row[0]] + (1..<row.size()).collect { j -> idx.collect { i -> row[j][i] } } }
    ch_res_leads = STEPWISE.out.res.filter { _s, _m, cls, _st, _sg, _ld -> cls == 'residue' }
        .map { s, m, _cls, _st, _sg, ld -> tuple(s, m, ld) }
        .groupTuple(by: 0, size: 2)
        .map { s, ms, lds -> tuple(s, ['fixed', 'random'].collect { x -> lds[ms.indexOf(x)] }) }
    ch_omni_given_in = ch_omni_in
        .combine(ch_pooled_set, by: 0)
        .map { s, cls, cs, bed, bim, fam, fx, keep, cs2, pbed, pbim, pfam ->
            if (cs != cs2) error "OMNIBUS_GIVEN: cohort order differs for ${s}"
            tuple(s, cls, cs, bed, bim, fam, fx, keep, pbed, pbim, pfam) }
        .combine(ch_res_leads, by: 0)
        .combine(OMNIBUS.out.res.map { s, cls, om, _sg, _st -> tuple(s, cls, om) }, by: [0, 1])
    OMNIBUS_GIVEN(ch_omni_given_in, GWAS_LEADS.out.ids, script_file('omnibus_given.py'), script_file('omnibus.py'),
                  script_file('vocab.py'))

    CONDITIONAL_COMPARE(
        STEPWISE.out.res.map { _s, _m, _cls, st, _sg, _ld -> st }.collect(),
        STEPWISE.out.res.map { _s, _m, _cls, _st, sg, _ld -> sg }.collect(),
        CONDITIONAL.out.res.map { _s, _c, _m, hg, _lg -> hg }.collect(),
        CONDITIONAL.out.res.map { _s, _c, _m, _hg, lg -> lg }.collect(),
        BUILD_POOLED.out.pooled.map { _s, _c, _b, _bi, _f, mk -> mk }.collect(),
        GWAS_LEADS.out.table, script_file('conditional_compare.py'), script_file('vocab.py'))

    // ── STAGE 5 — omnibus tests of residue positions: the summary (OMNIBUS runs before CONDITIONAL) ──
    OMNIBUS_SUMMARY(
        OMNIBUS.out.res.map { _s, _cls, om, _sg, _st -> om }.collect(),
        OMNIBUS.out.res.map { _s, _cls, _om, sg, _st -> sg }.collect(),
        OMNIBUS.out.res.map { _s, _cls, _om, _sg, st -> st }.collect(),
        COMPARE.out.tables.map { _s, a, _q, _c, _m, _h -> a }.collect(),
        CONDITIONAL_COMPARE.out.tables.map { sg, _st, _lg, _hg -> sg },
        cohortList.join(','), script_file('omnibus_summary.py'), script_file('vocab.py'))

    // ── figures: one task per family, over every case set ──────────────────
    // collect(sort: true): the staged order, and so the task hash, is stable.
    FIGURES(
        channel.fromList(params.Families).map { f -> tuple(f[0], f[1]) },
        COMPARE.out.tables.map { _s, a, _q, _c, _m, _h -> a }.collect(sort: true),
        COMPARE.out.tables.map { _s, _a, _q, _c, m, _h -> m }.collect(sort: true),
        EFFECTIVE_TESTS.out.tables.map { _s, _c, _cls, th, _t -> th }.collect(sort: true),
        STEPWISE.out.long.collect(sort: true),
        CONDITIONAL_COMPARE.out.tables,
        GWAS_LEADS.out.table, GENE_COORDS.out.coords,
        script_file('figures.py'), script_file('catalogue.py'), script_file('vocab.py'))

    // ── provenance ─────────────────────────────────────────────────────────
    def manifest = """{
  "component": "assoc_hla.v2",
  "case_sets": ${groovy.json.JsonOutput.toJson(params.CaseSets)},
  "case_set_definition": "all = every case; dvt_yes / dvt_no = cases whose ${params.DvtCol} in ${params.v6_info} is Yes / No (Unknown or blank: in 'all' only); these three against all of the cohort's controls. dvt_case_only = the dvt_yes and dvt_no cases only, no controls, PHENO1 recoded 1 = DVT yes / 0 = DVT no: the test of a difference between the two",
  "cohorts": ${groovy.json.JsonOutput.toJson(cohortList)},
  "cohort_nesting": "narrow_mainland is a strict subset of intermediate_mainland is a strict subset of full_mainland; they differ mainly in controls. Agreement across cohorts is a sample-selection sensitivity analysis, not replication.",
  "markers": "hla.typing.v2 results/06.assoc_inputs/<cohort>/{allele,residue}_dosage.tsv, read unchanged (md5 checked against summary.json); 2-field alleles and amino-acid residues at IMGT positions; 0/1/2 copies, blank = undetermined",
  "markers_tested": ${groovy.json.JsonOutput.toJson(cohortList.collectEntries { c -> [(c): [allele: typed.cohorts[c].alleles_tested, residue: typed.cohorts[c].residues_tested]] })},
  "coding": "VCF hard calls, REF A = copy does not carry the marker, ALT P = copy carries it; effects are per copy carrying the marker",
  "positions": "chromosome 6, gene start (RefSeq hg38_refGene, union of chr6 transcripts) + rank of the marker within its gene; labels, not physical residue positions",
  "phenotype": "${params.PhenoName}, 0 = control / 1 = case (plink2 --1)",
  "covariates": "${params.CovarName} (${params.PcLabel})",
  "fixed_model": "plink2 --glm omit-ref ${params.FirthMode} hide-covar --ci 0.95 on the 12_model_inputs fixed_model sample set (relatives removed) restricted to the case set; complete case per marker",
  "random_model": "SAIGE step 2 (--LOCO=TRUE --chrom=6 --AlleleOrder=alt-first --is_Firth_beta=FALSE --impute_method=${params.ImputeMethod} --maxMissing=${params.MaxMissing} --minMAC=${params.MinMAC}) on the case set's samples",
  "null_models": "all: analysis/assoc_saige 'main' null model, reused. dvt_yes / dvt_no / dvt_case_only: fitted here with SAIGE step 1 (binary, full GRM, LOCO, isDiagofKinSetAsOne, ${params.VarianceRatioMarkers} variance-ratio markers, SEX + 10 PCs) on the random_model genotypes restricted to analysis/assoc_saige 00.prep/grm.prune.in and the case set's samples. Every null model checked against its case set (CHECK_NULL).",
  "random_model_p": "p = p.value.NA (no SPA), as in the genome-wide SAIGE scan; p_spa = p.value kept as a secondary column; se = 1/sqrt(var), consistent with p.value.NA; OR = exp(BETA)",
  "tests": "per (case set, cohort, class): markers polymorphic in the case set's random-model samples; two residues at one IMGT position whose dosages are complements in every sample (a mirror pair) count as one test",
  "significance": "Per (case set, cohort, class), family-wise ${params.Alpha} over the tests, the same for both models, fixed before any P: PRIMARY LD-adjusted p < ${params.Alpha} / Me, Me = effective number of independent tests from GEC v0.2 (Li et al. 2012, standalone release archived in tools/gec-0.2, sha256 checked; --maf 0, on the case set's fixed-model = unrelated samples); VALIDATION exact max-T threshold (the ${params.Alpha} quantile of the minimum P under a multivariate normal with the correlation of the covariate-adjusted, mean-imputed dosages of the random-model samples; ${params.MaxTDraws} draws, seed ${params.MaxTSeed}); SENSITIVITY conventional Bonferroni p < ${params.Alpha} / tests. No Benjamini-Hochberg. Not corrected across case sets, cohorts or models. lambda_GC reported, never applied.",
  "benchmark": "MAXT_BOOT, no call depends on it: per (case set, cohort, class), ${params.MaxTBootDraws} null phenotypes drawn by parametric bootstrap from the null logistic model (${params.PhenoName} ~ 1 + SEX + PC1-10) fitted on the case set's fixed-model = unrelated samples; on every replicate, the tests (one per mirror pair) rerun as plink2 Wald tests with the fixed model's settings and as score tests without SPA (blanks mean-imputed; the random model's statistic without the GRM); exact threshold = the ${params.Alpha} quantile of the smallest P; the family-wise error rate of the GEC and Bonferroni thresholds = the share of replicates whose smallest P is below them; seed ${params.MaxTSeed}",
  "omnibus": "MAIN ANALYSIS, OMNIBUS: per case set, each residue position tested as a whole in the fixed model (the case set's unrelated samples, ${params.PhenoName} ~ 1 + SEX + PC1-10); gene-level omnibus tests over alleles are not run: categories = the position's residues plus 'unlisted' (copies the matrix does not list); a category < 1 % joins 'rare', a 'rare' pool < 0.5 % joins the reference (the most frequent category); likelihood-ratio test with df = the rank the unit adds; Bonferroni ${params.Alpha} / units tested per cohort; a separated or unfinished fit is never significant; independent units selected stepwise over the three cohorts with the rule of STEPWISE, conditioning on every category of the entered units; summary beside the single-marker scans in _omnibus/",
  "consistency": "random model (primary): per threshold and class, the number of the three nested cohorts in which a marker is significant with the same direction of effect, written n/3 (1/3 = Single cohort); a sample-selection sensitivity analysis, not replication. The fixed model is shown for effect sizes and model concordance",
  "tools": {"plink2": "${params.plink2}", "saige_step1": "${params.saige_step1}", "saige_step2": "${params.saige_step2}", "python": "${params.python}"},
  "conditional": "with the scan's own engine on a pooled PLINK set per (case set, cohort) (alleles + residues + the GWAS MHC lead[s] + a SENTINEL first record that absorbs a SAIGE --condition defect on the first marker). (1) STEPWISE per (case set, model, class), the three cohorts jointly: candidates = markers significant at the cohort's GEC threshold in at least one cohort, conditional on the signals entered; lead = the candidate with the smallest largest-P over the cohorts; at most ${params.CondMaxSignals} signals; members = markers significant before a lead entered and in no cohort after; carrier alleles of a residue lead from hla.typing.v2 allele_match_map + residue_reference. (2) CONDITIONAL per (case set, cohort, model): every HLA marker given each GWAS lead; each GWAS lead given single HLA markers (top ${params.CondNTag} by |r|, top ${params.CondNTag} by marginal P, the signal leads of both classes) and given the signal leads jointly; and given WHOLE POSITIONS (every residue column but the most frequent) for the positions of the residue signal leads and the omnibus position leads. (3) HLA_GIVEN_SIGNAL per (case set, cohort, model): every HLA marker and the GWAS lead given each signal lead and each of those whole positions (completes HLA | signal: STEPWISE steps only within a class that has a signal). (4) OMNIBUS_GIVEN per case set: every residue position given one marker at a time, each GWAS lead and each residue signal lead of either model, as one more covariate in OMNIBUS's own test (missing calls mean-imputed); the recomputed marginal scan must equal OMNIBUS's",
  "gwas_leads": "analysis/assoc_saige 'main' (SAIGE, SEX + 10 bbj_mainland PCs, p.value.NA): peaks overlapping chr6:${params.MhcStart}-${params.MhcEnd} whose lead has P < ${params.GwasLeadP} in every cohort; genotypes from 12_model_inputs random_model. Missing lead calls are mean-imputed here (best_guess in the genome-wide scan), so the lead's marginal P here differs slightly from the scan's",
  "independence": "designed and run without reading any file, script, null model or result of analysis/assoc_hla",
  "figures": "results/figures/0N.<family>/, one README.md per family. Main: 01.landscape (bars.<case_set>.png: one bar per test spanning its P over the three cohorts, coloured by n/3 at GEC; manhattan.<case_set>.png: the cohorts as rows), 02.signals, 03.conditioning (conditioning.png, carriers.png), 04.gwas_link, 05.dvt. Supplementary: 06.thresholds, 07.models, 08.qq. Random model throughout except 07.models"
}"""
    def index = """# assoc_hla.v2 results

| layer | where |
|---|---|
| per case set (`all`, `dvt_yes`, `dvt_no`, `dvt_case_only`) | `<case_set>/` — same layout in each (below) |
| all results, one schema | `<case_set>/_comparison/assoc_all.tsv` |
| markers significant in any cohort, with n/3 per threshold | `<case_set>/_comparison/hits.tsv` (all markers: `marker_summary.tsv`) |
| tests and thresholds (GEC, max-T, Bonferroni) | `<case_set>/<cohort>/01.markers/<class>/effective_tests.<class>.tsv`, `tests.<class>.tsv`, `gec.<class>.log` |
| benchmark: exact max-T by parametric bootstrap, fixed-model samples | `<case_set>/<cohort>/01.markers/<class>/maxt_boot.<class>.tsv` (+ `.minp.tsv.gz` per replicate, `.observed.tsv`) |
| scan summaries | `<case_set>/_comparison/scan_qc_all.tsv`; fixed vs random: `model_concordance.tsv` |
| figures | `figures/` — main `01.landscape/`, `02.signals/`, `03.conditioning/`, `04.gwas_link/`, `05.dvt/`; supplementary `06.thresholds/`, `07.models/`, `08.qq/`; one `README.md` each |
| per-cohort stages | `<case_set>/<cohort>/00.prep`, `01.markers/<class>`, `02.assoc/<model>/<class>`, `03.tables`, `04.conditional` |
| DVT yes against DVT no | `_case_sets/` — `dvt_contrast.tsv`, `dvt_hits.tsv`, `scan_qc_case_sets.tsv`; the case-only test itself is `dvt_case_only/`; figures `figures/02.signals/`, `05.dvt/` |
| conditional analyses | `_conditional/` — `signals_all.tsv` (independent signals), `consistency.tsv`, `summary.tsv`, `steps_all.tsv`, `hla_given_lead_all.tsv`, `lead_given_hla_all.tsv`, `ld_lead_hla.tsv`, `gwas_leads.tsv`; figures `figures/03.conditioning/`, `04.gwas_link/`; per walk: `<case_set>/_conditional/<model>/`; per GWAS-lead run: `<case_set>/<cohort>/04.conditional/<model>/` (`hla_given_lead.tsv`, `lead_given_hla.tsv`; every HLA marker given each signal lead and whole position: `hla_given_signal.tsv`) |
| omnibus tests of residue positions (a main analysis; fixed model, Bonferroni) | `_omnibus/` — `omnibus_units.tsv` (each unit beside its best single markers), `omnibus_signals_all.tsv`, `omnibus_all.tsv`, `omnibus_steps_all.tsv`; per case set: `<case_set>/_omnibus/` (positions given each GWAS lead and residue signal lead: `omnibus_given.residue.tsv`) |
| reference tables | `_reference/hla_gene_coords.tsv`, `_reference/case_sets.tsv` |
| run record | `_run_info/run_manifest.json`, `trace.txt`, `report.html`, `timeline.html`, `dag.html` |

Case sets: ${params.CaseSets.join(', ')} (the first three against all controls; dvt_case_only = DVT yes vs DVT no within cases). Cohorts: ${cohortList.join(', ')}. Classes: ${params.Classes.join(', ')}. Models: random (SAIGE, primary), fixed (plink2). Thresholds: LD-adjusted GEC (primary), exact max-T (validation), conventional Bonferroni (sensitivity); consistency = significant in n/3 cohorts. Column dictionaries: `../docs/OUTPUTS.md`; methods: `../docs/METHODS.md`."""
    WRITE_RUN_MANIFEST(channel.value(manifest), channel.value(index),
                       FIGURES.out.readme.collect().map { rs -> rs.size() })
}

// The same process, used twice: once for the reused null models, once for the
// fitted ones (a process can be invoked only once per workflow).
workflow NULL_INFO_ALL {
    take:
    ch
    script
    main:
    NULL_MODEL_INFO(ch, script)
    emit:
    NULL_MODEL_INFO.out.info
}

workflow NULL_INFO_FIT {
    take:
    ch
    script
    main:
    NULL_MODEL_INFO(ch, script)
    emit:
    NULL_MODEL_INFO.out.info
}
